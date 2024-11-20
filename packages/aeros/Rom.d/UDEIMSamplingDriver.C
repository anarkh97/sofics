#include "UDEIMSamplingDriver.h"

#include "ElementSamplingDriver.h"

#include "VecNodeDof6Conversion.h"
#include "NodalRestrictionMapping.h"
#include "ConnectivityUtils.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"
#include "BasisOps.h"

#include "RenumberingUtils.h"
#include "MeshDesc.h"

#include "DistrBasisFile.h"
#include <Timers.d/StaticTimers.h>

#include <Driver.d/Domain.h>
#include <Driver.d/SysState.h>
#include <Driver.d/GeoSource.h>
#include <Math.d/DBSparseMatrix.h>
#include <Utils.d/dofset.h>
#include <Utils.d/DistHelper.h>

#ifdef USE_EIGEN3
#include <Eigen/Dense>
#endif
#include <cmath>
#include <utility>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <fstream>

namespace Rom {

UDEIMSamplingDriver::UDEIMSamplingDriver(Domain *domain) :
  SingleDomainDynamic(domain),
  nodeDofMap(NULL),
  converter(NULL)
{}

void
UDEIMSamplingDriver::solve() {

  SingleDomainDynamic::preProcess();
  if(domain->solInfo().newmarkBeta == 0) {
    domain->assembleNodalInertiaTensors(melArray);
  }

  converter = new VecNodeDof6Conversion(*domain->getCDSA());
  nodeDofMap = new VecNodeDof6Map(*domain->getCDSA());

  const int podSizeMax = domain->solInfo().maxSizePodRom; 
  bool normalized = true;

  domain->createKelArray(kelArrayCopy);                    // initialize linear stiffness matrix

    VecBasis unassembledForceBuf;                          //container for unassembled force snapshots/svd
    VecBasis assembledForceBuf;                            //container for assembled force svd
    std::vector<int> umaskIndices;                         //unassembled indices container
    std::vector<int> amaskIndices;                         //assembled indices container
    std::set<int> selectedElemRank;                        //selected element rank container for unpacked element numbers
    std::vector<double> singularVals;                      //singular Values for error estimator
    std::vector<std::pair<int,int> > elemRankDOFContainer; //container for each selected element's slected DOF

    readInBasis(podBasis_, BasisId::STATE, BasisId::POD, podSizeMax, normalized);                              //get POD projection basis

    if(domain->solInfo().computeForceSnap){
      writeUnassembledForceSnap(unassembledForceBuf);                                        //get force snapshots
    } else {
      readUnassembledForceSnap(unassembledForceBuf,singularVals);                                              //read basis from file
      assembleBasisVectors(assembledForceBuf, unassembledForceBuf);
      computeInterpIndices(unassembledForceBuf,umaskIndices);                                                  //get UDEIM indices
      computeAssembledIndices(umaskIndices,amaskIndices,selectedElemRank,elemRankDOFContainer);                //map unassembled indices to assembled indices
      computeAndWriteUDEIMBasis(unassembledForceBuf,assembledForceBuf,umaskIndices,amaskIndices,singularVals); //compute and write to file the UDEIM POD basis
      writeSampledMesh(amaskIndices,selectedElemRank,elemRankDOFContainer);                                    //write UDEIM sampled mesh to file
      if(domain->solInfo().statePodRomFile.size() ==1 && domain->solInfo().velocPodRomFile.size() ==1)
        computeErrorBound(umaskIndices);
    }
}

void
UDEIMSamplingDriver::computeErrorBound(std::vector<int> &umaskIndices){

#ifdef USE_EIGEN3
  std::cout << "Computing Error Bound" << std::endl;
  VecBasis aForceSnapshots;
  readInBasis(aForceSnapshots, BasisId::STATE, BasisId::SNAPSHOTS,0);

  std::ifstream uForceFile;
  FileNameInfo fileInfo;
  std::string fileName = BasisFileId(fileInfo, BasisId::VELOCITY, BasisId::SNAPSHOTS);
  uForceFile.open(fileName.c_str());

  double squareNorm     = 0.0;
  double squareNormDiff = 0.0;

  if(uForceFile.is_open()){
    std::cout << "reading unassembled snapshots from " << fileName.c_str() << std::endl;
    int lengthVec; int numSnaps;
    if(uForceFile.good()){
      uForceFile >> lengthVec; std::cout << "vector length = " << lengthVec << std::endl; 
      uForceFile >> numSnaps;  std::cout << "number of snapshots = " << numSnaps << std::endl;
    }else{
      throw std::runtime_error("... Bad File ...");
    }

    for(int column = 0; column != numSnaps; column++){
      Vector dummy1(podBasis_.numVec());
      Vector uforceSnapshot(unassembledVecInfo());
      Eigen::Matrix<double,Eigen::Dynamic,1>  dummy2(podBasis_.numVec());
      Eigen::Matrix<double,Eigen::Dynamic,1>  buffer(umaskIndices.size());
 
      podBasis_.reduce(aForceSnapshots[column],dummy1);
      squareNorm += dummy1.squareNorm();
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > dummy1Map(dummy1.data(),dummy1.size());

      double timeStamp;
      uForceFile >> timeStamp; 
      std::cout << "\r timestamp = " << timeStamp;
      for(int i = 0; i != lengthVec; i++){
        if(uForceFile.good()){
          uForceFile >> uforceSnapshot[i];
        }else{
          throw std::runtime_error("... Bad File ...");
        }
      }

      for(int row = 0; row != buffer.rows(); row++){
        buffer(row) = uforceSnapshot[umaskIndices[row]];
      }

      dummy2 = compressedDBTranspose*buffer;
      dummy1Map -= dummy2;
      squareNormDiff += dummy1Map.squaredNorm();
    }

    squareNormDiff = sqrt(squareNormDiff);
    squareNorm     = sqrt(squareNorm);
 
    std::cout << "\n||V^T(F-U(P^TU)^-1(P^T)F)||_F/||V^TF||_F = " << squareNormDiff/squareNorm << std::endl;

  }else{
    throw std::runtime_error("... Bad File ...");
  }

#endif
}

void
UDEIMSamplingDriver::readInBasis(VecBasis &podBasis, BasisId::Type type, BasisId::Level level, int podSizeMax, bool normalized)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level);
 if(normalized) fileName.append(".massorthonormalized");
 BasisInputStream<6> in(fileName, *converter);
 if (podSizeMax != 0) {
   std::cout << "reading in " << podSizeMax << " vectors from " << fileName.c_str() << std::endl;
   readVectors(in, podBasis, podSizeMax);
 } else {
   std::cout << "reading in all vectors from " << fileName.c_str() << std::endl;
   readVectors(in, podBasis);
 }
}

template <typename Scalar>
void
UDEIMSamplingDriver::writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type, BasisId::Level level)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level);
 BasisOutputStream<6> outputNormalized(fileName, *converter, false);
 filePrint(stderr, " ... Writing basis to file %s ...\n", fileName.c_str());
 for (int iVec = 0; iVec < OutputBasis.vectorCount(); ++iVec) {
   if(singularValue.size() > 0)
     outputNormalized << std::make_pair(double(singularValue[iVec]), OutputBasis[iVec]);
   else
     outputNormalized << OutputBasis[iVec];
 }
}

void
UDEIMSamplingDriver::writeUnassembledForceSnap(VecBasis &unassembledForceBasis) 
{
#ifdef USE_EIGEN3
  //First read in state snapshots 
  VecBasis displac;
  std::vector<double> timeStamps;
  std::vector<int> snapshotCounts;
  readAndProjectSnapshots(BasisId::STATE, converter->vectorSize(), podBasis_, converter,
                          snapshotCounts, timeStamps, displac);

  //read in velocity snapshots if provided
  VecBasis *veloc_;
  if(!domain->solInfo().velocPodRomFile.empty()) {
   std::vector<double> velTimeStamps;
   std::vector<int> velSnapshotCounts;
   veloc_ = new VecBasis;
   readAndProjectSnapshots(BasisId::VELOCITY, converter->vectorSize(), podBasis_, converter,
                           velSnapshotCounts, velTimeStamps, *veloc_);
   if(velSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent velocity snapshots\n";
  } else { veloc_ = NULL;}

  //read in acceleration snapshots if provided
  VecBasis *accel_;
  if(!domain->solInfo().accelPodRomFile.empty()) {
   std::vector<double> accTimeStamps;
   std::vector<int> accSnapshotCounts;
   accel_ = new VecBasis;
   readAndProjectSnapshots(BasisId::ACCELERATION, converter->vectorSize(), podBasis_, converter,
                           accSnapshotCounts, accTimeStamps, *accel_);
   if(accSnapshotCounts != snapshotCounts) std::cerr << " *** WARNING: inconsistent acceleration snapshots\n";
  } else { accel_ = NULL; } 

  //Now build forcevectors
  std::vector<double> SVs;
  buildForceArray(unassembledForceBasis,displac,veloc_,accel_,timeStamps,snapshotCounts);
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > unassembledMap(unassembledForceBasis.data(),unassembledForceBasis.size(),unassembledForceBasis.numVectors());

  //output unassembled force snapshots for use in error bound computation
  std::ofstream snpFile;
  FileNameInfo fileInfo;
  std::string snpFileName = BasisFileId(fileInfo, BasisId::FORCE, BasisId::SNAPSHOTS);
  snpFileName.erase(snpFileName.end()-3,snpFileName.end());
  snpFileName += "snap";
  snpFile.open(snpFileName.c_str());
 
  snpFile << std::setprecision(16);
  snpFile << unassembledMap.rows() << " " << unassembledMap.cols() << std::endl;
  for(int column = 0; column != unassembledMap.cols(); column++){
    snpFile << "\t" << timeStamps[column] << std::endl;
    snpFile << unassembledMap.col(column) << std::endl;
  }
  snpFile.close();

  OrthoForceSnap(unassembledForceBasis,SVs); //orthogonalize unassembled vectors
  //output orthogonalized force snapshots for use in UDEIM basis
  std::ofstream svdFile;
  std::string svdFileName = BasisFileId(fileInfo, BasisId::FORCE, BasisId::SNAPSHOTS);
  svdFile.open(svdFileName.c_str());

  svdFile << uDOFaDOFmap.size() << " " << unassembledMap.rows() << " " << unassembledMap.cols() << std::endl;
  for(std::map<int, std::pair<int,int> >::const_iterator it = uDOFaDOFmap.begin(); it != uDOFaDOFmap.end(); it++)
    svdFile << it->first << " " << it->second.first << " " << it->second.second << std::endl;
 
  svdFile << std::setprecision(16); 
  for(int column = 0; column != unassembledMap.cols(); column++){
    svdFile << SVs[column] << std::endl;
    svdFile << unassembledMap.col(column) << std::endl;
  }
  svdFile.close();
#endif
}
 
void
UDEIMSamplingDriver::assembleBasisVectors(VecBasis &assembledForceBasis, VecBasis &unassembledForceBasis) 
{ 
  assembledForceBasis.dimensionIs(unassembledForceBasis.vectorCount(), SingleDomainDynamic::solVecInfo());
  //assemble the unassembledSVD

  for(int i = 0; i != assembledForceBasis.vectorCount(); i++)  
    for(int j = 0; j != assembledForceBasis.size(); j++) 
      assembledForceBasis[i][j] = 0;

  for (std::map<int, std::pair<int,int> >::const_iterator it = uDOFaDOFmap.begin(); it != uDOFaDOFmap.end(); it++){
    int iele = it->second.first;
    int idof = it->second.second;
    int aDofNum = domain->getCDSA()->getRCN((*domain->getAllDOFs())[iele][idof]); 
    for(int vec = 0; vec != unassembledForceBasis.numVectors(); vec++)
      assembledForceBasis[vec][aDofNum] += unassembledForceBasis[vec][it->first];
  }
}

void
UDEIMSamplingDriver::readUnassembledForceSnap(VecBasis &unassembledForceBasis, std::vector<double> &SVs){
  std::ifstream svdFile;
  FileNameInfo fileInfo;
  std::string svdFileName = BasisFileId(fileInfo, BasisId::FORCE, BasisId::SNAPSHOTS);
  svdFile.open(svdFileName.c_str());

  int mapSize;
  int numRows;
  int numCols = domain->solInfo().forcePodSize;
  int maxCols;

  std::cout << "reading in unassembled force snapshots from " << svdFileName.c_str() << std::endl;

  if(svdFile.is_open()){
 
   if(svdFile.good())
     svdFile >> mapSize;
   else
       throw std::runtime_error("... Bad File ...");

   if(svdFile.good())
     svdFile >> numRows; 
    else
       throw std::runtime_error("... Bad File ...");

   if(svdFile.good())
     svdFile >> maxCols;
    else
       throw std::runtime_error("... Bad File ...");
   
   std::cout << "Map size = " << mapSize << " Unassembled Length = " << numRows << " Maximum Basis size = " << maxCols << std::endl;
   unassembledForceBasis.dimensionIs(numCols,numRows);
    
   int first, second1, second2;
   for(int i = 0; i != mapSize; i++){
     if(svdFile.good())
       svdFile >> first; 
     else
       throw std::runtime_error("... Bad File ...");

     if(svdFile.good())
       svdFile >> second1; 
     else
       throw std::runtime_error("... Bad File ...");

     if(svdFile.good())
       svdFile >> second2;
     else
       throw std::runtime_error("... Bad File ...");

     uDOFaDOFmap.insert(std::make_pair(first,std::make_pair(second1,second2)));
   }
   
   std::cout << "Singular Values:" << std::endl;
   double element, singularValue;
   for(int col = 0; col != numCols; col++){
     if(svdFile.good())
       svdFile >> singularValue;
     else
       throw std::runtime_error("... Bad File ...");

     SVs.push_back(singularValue);
     std::cout << singularValue << " ";
     for(int row = 0; row != numRows; row++){
       svdFile >> element;
       unassembledForceBasis[col][row] = element;
     }
   }
   std::cout << std::endl;
   svdFile >> singularValue;
   SVs.push_back(singularValue);
 } else {
   throw std::runtime_error("... File not Open ...");
 }
}

void
UDEIMSamplingDriver::computeInterpIndices(VecBasis &forceBasis, std::vector<int> &maskIndices) {
  //member function for determining the sampling indicies for UDEIM
  //result is (U*(P^T*U)^-1) where U is the left singular vectors of force snapshots and P
  //is the column selection matrix
  
#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > forceMatrix(forceBasis.data(),forceBasis.vectorSize(),forceBasis.vectorCount());

  int maxCoeffSlot;

  {
   Eigen::Matrix<double,Eigen::Dynamic,1> firstCol(forceMatrix.rows());
   firstCol = forceMatrix.col(0);
   //take absolute value of first column
   for(int row = 0; row != firstCol.rows(); row++)
     if(firstCol(row) < 0.0)
       firstCol(row) = -1.0*firstCol(row);
   //find maximum component
   firstCol.maxCoeff(&maxCoeffSlot);
   maskIndices.push_back(maxCoeffSlot);
  }

  //start loop to compute mask indicies 
  for(int i = 1; i < domain->solInfo().forcePodSize; ++i){ //loop starts at 1 i.e. the 2nd column
    filePrint(stderr,"\r %4.2f%% complete", double(i)/double(domain->solInfo().forcePodSize)*100.);

    //allocate space for P^T*U and P^T*u_i
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Umasked(maskIndices.size(),maskIndices.size());
    Eigen::Matrix<double,Eigen::Dynamic,1>              u_i_masked(maskIndices.size());

    for(int j = 0; j < maskIndices.size(); ++j) {//select proper rows and columns of force basis
      Umasked.row(j) = forceMatrix.block(maskIndices[j],0,1,maskIndices.size()); //(P^T*U) is square
      u_i_masked(j) = forceMatrix(maskIndices[j],i);//mask next column over
    }

    Eigen::FullPivLU< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > luOfUmasked(Umasked); //invert row reduced basis

    Eigen::Matrix<double,Eigen::Dynamic,1> residual(forceBasis.vectorSize());

    if(luOfUmasked.isInvertible())
      residual = forceMatrix.col(i) - forceMatrix.leftCols(maskIndices.size())*luOfUmasked.inverse()*u_i_masked;
    else //check for invertability as a safegaurd against misuse of snapshots
      throw std::runtime_error("... Matrix Not Invertible ...");

    //take absolute value of residual vector
    for(int row = 0; row != residual.rows(); row++)
      if(residual(row) < 0.0)
        residual(row) = -1.0*residual(row);

    residual.maxCoeff(&maxCoeffSlot); //find indice of maximum component
    maskIndices.push_back(maxCoeffSlot);
  }
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);

  Eigen::Map< Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> > indSol(maskIndices.data(),maskIndices.size(),1);
  std::cout << "unassembled selected indices:" << indSol.transpose() << std::endl;
#endif
}

void
UDEIMSamplingDriver::getFullElemIndices(int selectedElem,std::set<int> &uAuxilary, std::set<int> &aAuxilary){
//this function will add all dofs for a given element. 

  int NumElemDof =  kelArrayCopy[selectedElem].dim();

  //loop over all dofs for given element. Determing which of those DOFS correspond to the chosen node
  for(int ElemDof = 0; ElemDof != NumElemDof; ElemDof++){
    int candidateAssemInd = domain->getCDSA()->getRCN((*domain->getAllDOFs())[selectedElem][ElemDof]);
    if(candidateAssemInd >= 0){ //make sure indice is in domain
    aAuxilary.insert(candidateAssemInd);//if so add assembled and unassembled indices to auxilary sets 
    for(std::map<int,std::pair<int,int> >::const_iterator it = uDOFaDOFmap.begin(); it != uDOFaDOFmap.end(); it++)//find unassembled rank for elem/dof pair
      if(it->second.first == selectedElem && it->second.second == ElemDof)
        uAuxilary.insert(it->first);
    }
  }
}

void
UDEIMSamplingDriver::getFullNodeIndices(int selectedElem,int assembledInd, std::set<int> &uAuxilary, std::set<int> &aAuxilary){
//this functiona dds all dofs for a selected node
  int selectedNode = nodeDofMap->nodeDof(assembledInd).nodeRank;
 
  int NumElemDof =  kelArrayCopy[selectedElem].dim();

  //loop over all dofs for given element. Determing which of those DOFS correspond to the chosen node
  for(int ElemDof = 0; ElemDof != NumElemDof; ElemDof++){
    int candidateAssemInd = domain->getCDSA()->getRCN((*domain->getAllDOFs())[selectedElem][ElemDof]);//see which assembled indice this belongs to
    if(candidateAssemInd >= 0){ //make sure indice is in domain
      int candidateNode = nodeDofMap->nodeDof(candidateAssemInd).nodeRank;//get corresponding node
      if(candidateNode == selectedNode){//see if current elem dof is in given node
        aAuxilary.insert(candidateAssemInd);//if so add assembled and unassembled indices to auxilary sets 
        for(std::map<int,std::pair<int,int> >::const_iterator it = uDOFaDOFmap.begin(); it != uDOFaDOFmap.end(); it++)//find unassembled rank for elem/dof pair
 	  if(it->second.first == selectedElem && it->second.second == ElemDof)
            uAuxilary.insert(it->first);
      }
    }
  }
}

void 
UDEIMSamplingDriver::computeAssembledIndices(std::vector<int> &umaskIndices, std::vector<int> &amaskIndices, std::set<int> &selectedElemRank, std::vector<std::pair<int,int> > &elemRankDOFContainer)
{
 //convert selected unassembled indices to assembled indices for compatibility with online DEIM code
 //with mixed element types, each element could have a different number of dofs

 std::set<int> uAuxilaryIndices;
 std::set<int> aAuxilaryIndices;

 for(int i = 0; i != umaskIndices.size(); i++){
                      //unassembled to assembled DOF map
   int selectedElem = uDOFaDOFmap.at(umaskIndices[i]).first;
   int selectedDOF  = uDOFaDOFmap.at(umaskIndices[i]).second;

   int assembledInd = domain->getCDSA()->getRCN((*domain->getAllDOFs())[selectedElem][selectedDOF]);
   amaskIndices.push_back(assembledInd);   //assembled indices container

   if(domain->solInfo().selectFullNode)  
     getFullNodeIndices(selectedElem,assembledInd,uAuxilaryIndices,aAuxilaryIndices);

   if(domain->solInfo().selectFullElem)
     getFullElemIndices(selectedElem,uAuxilaryIndices,aAuxilaryIndices);

   selectedElemRank.insert(selectedElem);  //unpacked element container
   elemRankDOFContainer.push_back(std::make_pair(selectedElem,selectedDOF)); //map from selected element to selected DOF
 }

  if(domain->solInfo().selectFullNode || domain->solInfo().selectFullElem){
    for(int i = 0; i != amaskIndices.size(); i++)
      aAuxilaryIndices.erase(aAuxilaryIndices.erase(amaskIndices[i]));//remove redundant auxilary indices from set
    for(int i = 0; i != umaskIndices.size(); i++)
      uAuxilaryIndices.erase(uAuxilaryIndices.erase(umaskIndices[i]));

   for(std::set<int>::const_iterator it = uAuxilaryIndices.begin(); it != uAuxilaryIndices.end(); it++){
      int selectedAuxElem = uDOFaDOFmap.at(*it).first;
      int selectedAuxDOF  = uDOFaDOFmap.at(*it).second;

      elemRankDOFContainer.push_back(std::make_pair(selectedAuxElem,selectedAuxDOF)); //map from selected element to selected DOF
   }

   std::copy(aAuxilaryIndices.begin(),aAuxilaryIndices.end(),std::back_inserter(amaskIndices));//append auxilarry indices to end to selected indices list
   std::copy(uAuxilaryIndices.begin(),uAuxilaryIndices.end(),std::back_inserter(umaskIndices));
  }


#ifdef USE_EIGEN3
  Eigen::Map< Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> > indSol(amaskIndices.data(),amaskIndices.size(),1);
  std::cout << "assembled selected indices:" << indSol.transpose() << std::endl;
#endif

}

void
UDEIMSamplingDriver::computeAndWriteUDEIMBasis(VecBasis &unassembledForceBuf,VecBasis &assembledForceBuf, std::vector<int> &umaskIndices, std::vector<int> &amaskIndices, std::vector<double> &singularVals)
{
  //member function for computing and writing the UDEIM basis P*(P^T*U)^-T*U^T*V
  //where V is the mass-orthogonal POD basis, U is the left Singular vectors of the force snapshots
  //and P column selection matrix derived form the sampled indicies computed above
  int maxDeimBasisSize = domain->solInfo().maxDeimBasisSize;//set max basis size in case we want to solve least squares problem
  if(maxDeimBasisSize == 0)
    maxDeimBasisSize = unassembledForceBuf.size();//if not, we solve a square system

#ifdef USE_EIGEN3
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > podMap(podBasis_.data(),podBasis_.size(), podBasis_.vectorCount());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > uforceMap(unassembledForceBuf.data(), unassembledForceBuf.size(), unassembledForceBuf.vectorCount());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > aforceMap(assembledForceBuf.data(),assembledForceBuf.size(),assembledForceBuf.vectorCount());

  compressedDBTranspose(podBasis_.numVectors(),umaskIndices.size());
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> rowReducedUFM(umaskIndices.size(),maxDeimBasisSize);

  //initialize (P^T*U_unassembled) were U_unassembled is the svd of the unassembled force snapshots
  for(int row = 0; row != umaskIndices.size(); row++)
    rowReducedUFM.row(row) = uforceMap.block(umaskIndices[row],0,1,maxDeimBasisSize);

  //compute W^T = V^T*U_assembled*(P^T*U_unassembled)^-1
  Eigen::JacobiSVD< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > SVDOfUmasked(rowReducedUFM,Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix<double,Eigen::Dynamic,1> invSVs(SVDOfUmasked.nonzeroSingularValues());

  for(int i = 0; i != invSVs.rows(); i++) invSVs(i) = 1.0/SVDOfUmasked.singularValues()(i);

  std::cout << "condition Number of (P^T*U) = " << SVDOfUmasked.singularValues()(0)/SVDOfUmasked.singularValues()(SVDOfUmasked.nonzeroSingularValues()-1) << std::endl;
  std::cout << "||(P^T*U)^-1)||_2 = " << invSVs(SVDOfUmasked.nonzeroSingularValues()-1) << std::endl;

  compressedDBTranspose = podMap.transpose()*aforceMap.leftCols(maxDeimBasisSize)*SVDOfUmasked.matrixV()*invSVs.asDiagonal()*SVDOfUmasked.matrixU().transpose();
  //we are computing the transpose of the basis

  std::cout << "compressed Basis" << std::endl;
  std::cout << compressedDBTranspose.transpose() << std::endl;

  //expand rows W^T*P^T
  std::cout << " num of rows = " << compressedDBTranspose.rows() << " num of cols = " << compressedDBTranspose.cols() << std::endl;
#endif
}

void
UDEIMSamplingDriver::writeSampledMesh(std::vector<int> &maskIndices, std::set<int> &selectedElemRank, std::vector<std::pair<int,int> > &elemRankDOFContainer) {

  //get nodes & dofs belonging to sampled indices (global numbering)
  std::set<int> selectedNodeSet; //To ensure unicity, since several locations (i.e. vector indices) can correspond to one node
  std::vector<std::pair<int,int> > compressedNodeKey; //need separate vector container for repeated nodes with different selected dofs
  for (std::vector<int>::const_iterator it = maskIndices.begin(); it != maskIndices.end(); ++it) {
    selectedNodeSet.insert(nodeDofMap->nodeDof(*it).nodeRank);
    int selectedNode = nodeDofMap->nodeDof(*it).nodeRank;
    int selectedNodeDof = std::log10(nodeDofMap->nodeDof(*it).dofId)/std::log10(2.0); //this function call is returning 2^dof for some reason, call log to normalize
    compressedNodeKey.push_back(std::make_pair(selectedNode,selectedNodeDof));
  }

  //print nodes to screen
  filePrint(stderr,"selected Nodes:");
  for(std::set<int>::iterator it = selectedNodeSet.begin(); it != selectedNodeSet.end(); it++)
    filePrint(stderr," %d", *it);
  filePrint(stderr,"\n");

  {//more standard output information
   DofSetArray *cdsa = domain->getCDSA();
   for(std::vector<std::pair<int,int> >::iterator it = compressedNodeKey.begin(); it != compressedNodeKey.end(); it++)
     std::cout << "node: " << it->first+1 << " dof1 = " << cdsa->firstdof(it->first)+1 << " slot = " << it->second << std::endl;
  }

  // Determine mapping between elements and nodes
  std::unique_ptr<Connectivity> elemToNode(new Connectivity(geoSource->getElemSet()->asSet()));

  //fill element container with all packed elements 
  std::vector<int> packedToInput(elementCount());
  {
   Elemset &inputElemSet = *(geoSource->getElemSet());
   for (int iElem = 0, iElemEnd = inputElemSet.size(); iElem != iElemEnd; ++iElem) {
     Element *elem = inputElemSet[iElem];
     if (elem) {
       const int iPackElem = domain->glToPackElem(iElem);
       if(iPackElem >= 0) {
         assert(iPackElem < packedToInput.size());
         packedToInput[iPackElem] = iElem;
       }
     }
   }
  }

  //construct element weight solution vector
  std::vector<double> solution(elementCount(),0.0);

  //fill weight map with selected elements
  std::vector<int> sampleElemIds;
  sampleElemIds.reserve(selectedElemRank.size());
  std::map<int, double> weights;
  for(std::set<int>::const_iterator it = selectedElemRank.begin(), it_end = selectedElemRank.end(); it != it_end; ++it){
    const int elemRank = packedToInput[*it];
    weights.insert(std::make_pair(elemRank, 1.0));
    solution[*it] = 1.0;
    sampleElemIds.push_back(elemRank);
  }
 
  //output Full Weights for compatibility with full mesh hyperreduction (i.e. old method)
  {
   // 1. write the ATTRIBUTES data
   outputFullWeights(solution, packedToInput);

   // 2. append the SNSLOT and UDBASIS data
   const std::string fileName = domain->solInfo().reducedMeshFile;
   const std::ios_base::openmode mode = std::ios_base::out | std::ios_base::app;
   std::ofstream weightOut(fileName.c_str(), mode);
   weightOut.precision(std::numeric_limits<double>::digits10+1);

   weightOut << "*\n";
   weightOut << "SNSLOT\n";
/*   int counter = 0;
   for(std::vector<std::pair<int,int> >::iterator it = compressedNodeKey.begin(); it != compressedNodeKey.end(); it++){
     int selElem = elemRankDOFContainer[counter].first;
     int selDOF  = elemRankDOFContainer[counter].second;
     //output assembled indices in form of a node plus a node dof and an element with and element dof (to keep neighbor elements for adding to the same dof)
     weightOut << it->first + 1 << " " << it->second << " " << packedToInput[selElem] + 1 << " " << selDOF << std::endl; 
     counter++;
   }*/
   for(int counter = 0; counter != elemRankDOFContainer.size(); counter++){
     int selElem = elemRankDOFContainer[counter].first;
     int selDOF  = elemRankDOFContainer[counter].second;
     //output basis column, element, element DOF
     weightOut << counter << " " << packedToInput[selElem] + 1 << " " << selDOF << std::endl;
   }

   weightOut << "*\nUDBASIS\n"  ;
#ifdef USE_EIGEN3
   weightOut << compressedDBTranspose.rows() << " " <<  compressedDBTranspose.cols() << std::endl;
   for(int row = 0; row < compressedDBTranspose.rows(); ++row){
     weightOut << compressedDBTranspose.row(row).transpose() << std::endl;
   } 
#endif
  }

  // compute the reduced forces (constant only). 
  Vector constForceFull(SingleDomainDynamic::solVecInfo());
  getConstForce(constForceFull);
  Vector constForceRed(podBasis_.vectorCount());
  reduce(podBasis_, constForceFull,  constForceRed);
  bool reduce_f = (constForceFull.norm() != 0);

  // compute the reduced initial conditions
  Vector d0Full(SingleDomainDynamic::solVecInfo()),
         v0Full(SingleDomainDynamic::solVecInfo());
  Vector tmp(SingleDomainDynamic::solVecInfo());
  SysState<Vector> inState(d0Full, v0Full, tmp, tmp);
  SingleDomainDynamic::getInitState(inState);
  Vector d0Red(podBasis_.vectorCount()),
         v0Red(podBasis_.vectorCount());
  bool reduce_idis = (d0Full.norm() != 0),
       reduce_ivel = (v0Full.norm() != 0);
  if(domain->solInfo().useMassNormalizedBasis || domain->solInfo().newmarkBeta == 0) {
    AllOps<double> allOps;
    if(reduce_idis || reduce_ivel) {
      //std::cerr << "building mass matrix\n";
      allOps.M = domain->constructDBSparseMatrix<double>();
      domain->makeSparseOps(allOps, 0, 0, 0);
    }
    if(reduce_idis) {
      allOps.M->mult(d0Full, tmp);
      reduce(podBasis_, tmp, d0Red);
    }
    if(reduce_ivel) {
      allOps.M->mult(v0Full, tmp);
      reduce(podBasis_, tmp, v0Red);
    }
  }
  else {
    if(reduce_idis) reduce(podBasis_, d0Full, d0Red);
    if(reduce_ivel) reduce(podBasis_, v0Full, v0Red);
  }
 
  // construct and print renumbered mesh
  const FileNameInfo fileInfo;
  const MeshRenumbering meshRenumbering(sampleElemIds.begin(), sampleElemIds.end(), *elemToNode, true);
  const MeshDesc reducedMesh(domain, geoSource, meshRenumbering, weights);
  outputMeshFile(fileInfo, reducedMesh, podBasis_.vectorCount());

  // initialize reduced mesh output file stream
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), std::ios_base::app);
  meshOut.precision(std::numeric_limits<double>::digits10+1);
  if(domain->solInfo().reduceFollower) meshOut << "EXTFOL\n";

  // reduced stiffness
  meshOut << "*\nREDSTIFF\n"  ;

  for(int column=0; column<podBasis_.vectorCount(); ++column){
    Vector columnOfKtimesV(solVecInfo());
    Vector columnOfRedK(podBasis_.vectorCount());
    columnOfKtimesV.zero();
    domain->getKtimesU(podBasis_[column], bcx, columnOfKtimesV, 1.0, kelArrayCopy);
    reduce(podBasis_, columnOfKtimesV, columnOfRedK);
    for(int i=0; i<podBasis_.vectorCount(); ++i){
      meshOut << columnOfRedK[i] << std::endl;}
  }

#ifdef USE_EIGEN3
  meshOut << "*\nUDBASIS\n"  ;
  meshOut << compressedDBTranspose.rows() << " " << compressedDBTranspose.cols() << std::endl;
  for(int row = 0; row < compressedDBTranspose.rows(); ++row){
    meshOut << compressedDBTranspose.row(row).transpose() << std::endl;
  }
#endif

  // output the reduced forces
  if(reduce_f) {
    meshOut << "*\nFORCES\nMODAL\n";
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " "  << constForceRed[i] << std::endl; 
  }


  // output the reduced initial conditions
  if(reduce_idis) {
    meshOut << "*\nIDISPLACEMENTS\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << d0Red[i] << std::endl;
  }
  if(reduce_ivel) {
    meshOut << "*\nIVELOCITIES\nMODAL\n";
    meshOut.precision(std::numeric_limits<double>::digits10+1);
    for(int i=0; i<podBasis_.vectorCount(); ++i)
      meshOut << i+1 << " " << v0Red[i] << std::endl;
  }

#ifdef USE_EIGEN3
  // build and output compressed basis
  podBasis_.makeSparseBasis(meshRenumbering.reducedNodeIds(), domain->getCDSA());
  {
    std::string PODfilename = BasisFileId(fileInfo, BasisId::STATE, BasisId::POD);
    std::string DEIMfilename = BasisFileId(fileInfo, BasisId::FORCE, BasisId::ROB);
    PODfilename.append(".reduced");
    DEIMfilename.erase(DEIMfilename.end()-6,DEIMfilename.end());
    DEIMfilename.append(".reduced.udeim");
    if(domain->solInfo().newmarkBeta == 0 || domain->solInfo().useMassNormalizedBasis) PODfilename.append(".massorthonormalized");
    filePrint(stderr," ... Writing compressed POD/UDEIM basis to file %s and %s ...\n", PODfilename.c_str(), DEIMfilename.c_str());
    DofSetArray reduced_dsa(reducedMesh.nodes().size(), const_cast<Elemset&>(reducedMesh.elements()));
    ConstrainedDSA reduced_cdsa(reduced_dsa, reducedMesh.dirichletBConds().size(), const_cast<BCond*>(&reducedMesh.dirichletBConds()[0]));
    VecNodeDof6Conversion converter(reduced_cdsa);
    BasisOutputStream<6> PODoutput(PODfilename, converter, false);
    BasisOutputStream<6> DEIMoutput(DEIMfilename, converter, false);

    for (int iVec = 0; iVec < podBasis_.vectorCount(); ++iVec) {
      PODoutput << podBasis_.compressedBasis().col(iVec);
    }

    std::map<int,int> nodeRenum(meshRenumbering.nodeRenumbering());
    std::map<int,int> elemRenum(meshRenumbering.elemRenumbering());

    meshOut << "*\nSNSLOT\n";
    int counter = 0;
    for(std::vector<std::pair<int,int> >::iterator it = compressedNodeKey.begin(); it != compressedNodeKey.end(); it++){
     int selElem = elemRankDOFContainer[counter].first;
     int selDOF  = elemRankDOFContainer[counter].second;
     //output assembled indices in form of a node plus a node dof and an element with and element dof (to keep neighbor elements for adding to the same dof)
     meshOut << nodeRenum[it->first] + 1 << " " << it->second << " " << elemRenum[packedToInput[selElem]]  + 1 << " " << selDOF <<std::endl;    
     counter++;
    }

  }
#endif 

}

void
UDEIMSamplingDriver::buildForceArray(VecBasis &unassembledForceBasis,const VecBasis &displac,const VecBasis *veloc,
                                     const VecBasis *accel,std::vector<double> timeStamps_,std::vector<int> snapshotCounts_)
{//this member function is for converting state snapshots to force snapshots in the absence of precollected force snapshots from model I
  //most of the code is copied from assembleTrainingData in ElementSamplingDriver.C
  std::vector<double>::iterator timeStampIt = timeStamps_.begin();

  //initialize assembled and unassembled snapshot containers
  unassembledForceBasis.dimensionIs(displac.vectorCount(), unassembledVecInfo());
  
  //temporary working arrays
  Vector unassembledTarget(unassembledVecInfo(),0.0);

  int iSnap = 0;
  for(int i = 0; i < snapshotCounts_.size(); i++) {//loop over snapshot sets

    for(int jSnap = 0; jSnap != snapshotCounts_[i]; ++iSnap, ++jSnap){//loop over snapshots in set
      filePrint(stderr,"\r %4.2f%% complete, time = %f", double(iSnap)/double(std::accumulate(snapshotCounts_.begin(),snapshotCounts_.end(),0))*100.,*timeStampIt);

      geomState->explicitUpdate(domain->getNodes(), displac[iSnap]);
      if(veloc){ geomState->setVelocity((*veloc)[iSnap]);} //just set the velocity at the nodes
      if(accel){ geomState->setAcceleration((*accel)[iSnap]);} //just set the acceleration at the nodes

      Vector dsp(solVecInfo());
      dsp = displac[iSnap];

      //special function call to get unassembled snapshots and mapping from unassembled to assembled DOFS
      SingleDomainDynamic::getUnassembledNonLinearInternalForce(dsp,unassembledTarget,uDOFaDOFmap,kelArrayCopy,*timeStampIt,jSnap);
   
      unassembledForceBasis[iSnap] = unassembledTarget;
      timeStampIt++;

      }

    }
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);
}

void UDEIMSamplingDriver::OrthoForceSnap(VecBasis &forceBasis,std::vector<double> &SVs)
{
#ifdef USE_EIGEN3
  std::cout << "  ... Orthogonalizing Snapshots ..." << std::endl;
  SVs.resize(forceBasis.numVectors());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > SingularValueMap(SVs.data(),forceBasis.numVectors());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ForceMap(forceBasis.data(), forceBasis.size(), forceBasis.numVectors());
  Eigen::JacobiSVD<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ForceSVD(ForceMap, Eigen::ComputeThinU);

  ForceMap         = ForceSVD.matrixU();
  SingularValueMap = ForceSVD.singularValues();
#endif
}

int UDEIMSamplingDriver::elementCount() const {
  return domain->numElements();
}

//member functions for fascilitating the reading of snapshot data
int 
UDEIMSamplingDriver::snapSize(BasisId::Type type, std::vector<int> &snapshotCounts)
{
  // compute the number of snapshots that will be used for a given skipping strategy
  FileNameInfo fileInfo;
  snapshotCounts.clear();
  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  for(int i = 0; i < FileNameInfo::size(type, BasisId::SNAPSHOTS); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    DistrBasisInputFile in(fileName);
    const double N = in.stateCount() - skipOffSet;
    const int singleSnapshotCount = (N > 0) ? 1+(N-1)/skipFactor : 0;
    snapshotCounts.push_back(singleSnapshotCount);
  }
  // return the total
  return std::accumulate(snapshotCounts.begin(), snapshotCounts.end(), 0);
}

void 
UDEIMSamplingDriver::readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
                                              const VecNodeDof6Conversion *vecDofConversion,
                                              std::vector<int> &snapshotCounts, std::vector<double> &timeStamps, VecBasis &config)
{
  const int snapshotCount = snapSize(type, snapshotCounts);
  filePrint(stderr, " ... Reading in %d %s Snapshots ...\n", snapshotCount, toString(type).c_str());

  config.dimensionIs(snapshotCount, vectorSize);
  timeStamps.clear();
  timeStamps.reserve(snapshotCount);

  const int skipFactor = std::max(domain->solInfo().skipPodRom, 1); // skipFactor must be >= 1
  const int skipOffSet = std::max(domain->solInfo().skipOffSet, 0); // skipOffSet must be >= 0
  const int podVectorCount = podBasis.vectorCount();
  Vector snapshot(vectorSize);
  Vector podComponents(podVectorCount);
  const FileNameInfo fileInfo;

  int offset = 0;
  for(int i = 0; i < FileNameInfo::size(type, BasisId::SNAPSHOTS); i++) {
    std::string fileName = BasisFileId(fileInfo, type, BasisId::SNAPSHOTS, i);
    filePrint(stderr, " ... Processing File: %s ...\n", fileName.c_str());
    BasisInputStream<6> in(fileName, *vecDofConversion);

    int count = 0;
    int skipCounter = skipFactor - skipOffSet;
    while(count < snapshotCounts[i]) {
      std::pair<double, double *> data;
      data.second = snapshot.data();
      in >> data;
      assert(in);
      if(skipCounter == skipFactor) {
    //    expand(podBasis, reduce(podBasis, snapshot, podComponents), config[offset+count]);
        config[offset+count] = snapshot;
        timeStamps.push_back(data.first);
        skipCounter = 1;
        ++count;
        filePrint(stderr, "\r ... timeStamp = %7.2e, %4.2f%% complete ...", data.first, double(count)/snapshotCounts[i]*100);
      }
      else {
        ++skipCounter;
      }
    }

    filePrint(stderr,"\n");
    offset += snapshotCounts[i];
  }

  assert(timeStamps.size() == snapshotCount);
}

int UDEIMSamplingDriver::unassembledVecInfo(){

  //compute number of unassembled indices on an element by element base
  //to minimize the amount of storage required
  //rather than just returning maxNumElementDOFs*numberOfElements
  int unassembledLength = 0;

  for(int iele = 0; iele != elementCount(); iele++){
    for(int idof = 0; idof != kelArrayCopy[iele].dim(); idof++){
      int uDofNum = domain->getCDSA()->getRCN((*domain->getAllDOFs())[iele][idof]);
      if(uDofNum >=0)
        unassembledLength += 1;
    }
  }
 
//  std::cout<<"unassembledLength = "<<unassembledLength<<std::endl;
  return unassembledLength;

}

} /* end namespace Rom */

Rom::DriverInterface *udeimSamplingDriverNew(Domain *domain) {
  return new Rom::UDEIMSamplingDriver(domain);
}
