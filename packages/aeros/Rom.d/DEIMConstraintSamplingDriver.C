#include "DEIMConstraintSamplingDriver.h"

#ifdef USE_EIGEN3
#include "ElementSamplingDriver.h"

#include "SvdOrthogonalization.h"
#include "VecNodeDof6Conversion.h"
#include "NodalRestrictionMapping.h"
#include "ConnectivityUtils.h"
#include "BasisFileStream.h"
#include "FileNameInfo.h"
#include "SimpleBuffer.h"
#include "BasisOps.h"

#include "NonnegativeMatrixFactorization.h"
#include <Element.d/MpcElement.d/MpcElement.h>

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

#include <Eigen/Dense>
#include <cmath>
#include <utility>
#include <algorithm>
#include <numeric>

namespace Rom {

DEIMConstraintSamplingDriver::DEIMConstraintSamplingDriver(Domain *domain) :
  SingleDomainDynamic(domain),
  converter6(NULL),
  converter1(NULL)
{}

void
DEIMConstraintSamplingDriver::solve() {

  SingleDomainDynamic::preProcess();
  if(domain->solInfo().newmarkBeta == 0) {
    domain->assembleNodalInertiaTensors(melArray);
  }

  converter6 = new VecNodeDof6Conversion(*domain->getCDSA());
  converter1 = new VecNodeDof1Conversion(domain->getNumCTC());

  std::cout << "size of constraint vec = " << domain->getNumCTC() << std::endl;

  bool normalized = true;

  VecBasis constraintBasis;

  if(domain->solInfo().computeConstraintSnap){         // construct and write constraint snapshots to file
//    readInBasis(podBasis, BasisId::STATE, BasisId::POD, true, 0, podSizeMax); // read in POD basis for snapsot projection
    writeProjConstraintSnap(); // constraint snapshots are written to file specified under OUTPUT: constvct
  } 

  if(domain->solInfo().orthogConstraintSnap) {  // orthogonalize constraint snapshots
    if(domain->solInfo().computeConstraintSnap){
      readInBasis(constraintBasis, BasisId::STATE, BasisId::SNAPSHOTS, false, 1); // read in contact snapshots under second entry of trnvct command
    } else {
      readInBasis(constraintBasis, BasisId::STATE, BasisId::SNAPSHOTS, false, 0); // read in contact snapshots under first entry of trnvct command
    }
    if(domain->solInfo().filterSnapshotRows){
     readInBasis(dualBasis, BasisId::DUALSTATE, BasisId::POD, false);
     filterRows(constraintBasis);
    }
    std::vector<double> SVs;
    OrthoConstraintSnap(constraintBasis,SVs);
    writeBasisToFile(constraintBasis, SVs, BasisId::STATE, BasisId::ROB,false); // write constraint basis to file specified under robdataf
  } 

  if ((domain->solInfo().computeConstraintSnap && domain->solInfo().orthogConstraintSnap) || 
      (!domain->solInfo().computeConstraintSnap && !domain->solInfo().orthogConstraintSnap)) { // selected interpolation indices and write DEIM basis
    readInBasis(dualBasis, BasisId::DUALSTATE, BasisId::POD, false); // Lagrange multiplier ROB is read from dualbasis specified under DEIM
    std::vector<int> maskIndices;
    std::vector<double> timeStamps;
    computeInterpIndices(constraintBasis, maskIndices);
    std::cout << " Constraint Basis size = " << constraintBasis.numVectors() << " x " << constraintBasis.size() << std::endl;
    std::cout << " dual Basis size       = " << dualBasis.numVectors() << " x " << dualBasis.size() << std::endl;
    computeAndWriteDEIMBasis(constraintBasis, maskIndices); // deim basis is written to file specified under dualbasis with ".deim" appended
    writeSampledMesh(maskIndices); // only output part of reduced mesh that contains contact specification
    if(domain->solInfo().statePodRomFile.size() > 0)
      computeErrorBound();
  }

}

void
DEIMConstraintSamplingDriver::computeErrorBound(){

  std::cout << " ... Computing Error Bound ... " << std::endl;
  VecBasis constraintSnapshots;
  std::vector<double> timeStamps;
  readInBasis(constraintSnapshots, BasisId::CONSTRAINT, BasisId::SNAPSHOTS,false);

  double squareNorm     = 0.0;
  double squareNormDiff = 0.0;

  for(int column = 0; column != constraintSnapshots.numVectors(); column++){
    Vector dummy1(deimBasis.numVec());
    Vector dummy2(dualBasis.numVec());
    deimBasis.reduce(constraintSnapshots[column],dummy1);
    dualBasis.reduce(constraintSnapshots[column],dummy2);

    squareNorm = dummy2.squareNorm();
    dummy2 = dummy2 - dummy1;
    squareNormDiff += dummy2.squareNorm();  
  }
  
  squareNormDiff = sqrt(squareNormDiff);
  squareNorm     = sqrt(squareNorm);

  std::cout << " ... ||V^T(F-U(P^TU)^-1(P^T)F)||_F/||V^TF||_F = " << squareNormDiff/squareNorm << " ... " <<  std::endl;

}

void
DEIMConstraintSamplingDriver::readInBasis(VecBasis &basis, BasisId::Type type, BasisId::Level level, bool vectorQuant, int i, int podSizeMax, bool normalized)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level, i);
 if(normalized) fileName.append(".massorthonormalized");

 if(vectorQuant){ // read in vector valued data
   BasisInputStream<6> in(fileName, *converter6);
   if (podSizeMax != 0) {
     std::cout << " ... reading in " << podSizeMax << " vectors from " << fileName.c_str() << " ... " << std::endl;
     readVectors(in, basis, podSizeMax);
   } else {
     std::cout << " ... reading in all vectors from " << fileName.c_str() << " ... " << std::endl;
     readVectors(in, basis);
   }
 } else {       // read in scalar valued data
     BasisInputStream<1> in(fileName, *converter1);
     if (podSizeMax != 0) {
       std::cout << " ... reading in " << podSizeMax << " vectors from " << fileName.c_str() << " ... " << std::endl;
       readVectors(in, basis, podSizeMax);
       if(type == BasisId::CONSTRAINT){ // read in next singular value to compute error bound
         Vector dummyV(basis.vectorInfo());
         double dummyD = 0;
         std::pair<double,Vector> foo;
         foo = std::make_pair(dummyD,dummyV);
         in >> foo;
         MPOSingularValue = foo.first;
       }
     } else {
       std::cout << " ... reading in all vectors from " << fileName.c_str() << " ... " << std::endl;
       readVectors(in, basis);
     }
 }

}

template <typename Scalar>
void
DEIMConstraintSamplingDriver::writeBasisToFile(const VecBasis &OutputBasis, std::vector<Scalar> singularValue, BasisId::Type type, BasisId::Level level, bool vectorQuant)
{
 FileNameInfo fileInfo;
 std::string fileName = BasisFileId(fileInfo, type, level);
 if(vectorQuant){  // write vector valued data
   BasisOutputStream<6> outputNormalized(fileName, *converter6, false);
   filePrint(stderr, " ... Writing basis to file %s ...\n", fileName.c_str());
   if(singularValue.size() > 0){
     for (int iVec = 0; iVec < singularValue.size(); ++iVec) {
       outputNormalized << std::make_pair(double(singularValue[iVec]), OutputBasis[iVec]);
     }
   } else {
     for (int iVec = 0; iVec < OutputBasis.numVec(); ++iVec) {
       outputNormalized << OutputBasis[iVec];
     }
   }
 } else {         // write scalar valued data
   BasisOutputStream<1> outputNormalized(fileName, *converter1, false);
   filePrint(stderr, " ... Writing basis to file %s ...\n", fileName.c_str());
   if(singularValue.size() > 0){
     for (int iVec = 0; iVec < singularValue.size(); ++iVec) {
       outputNormalized << std::make_pair(double(singularValue[iVec]), OutputBasis[iVec]);
     }
   } else {
     for (int iVec = 0; iVec < OutputBasis.numVec(); ++iVec) {
       outputNormalized << OutputBasis[iVec];
     }
   }
 }
}

void
DEIMConstraintSamplingDriver::writeProjConstraintSnap() 
{
  //First read in state snapshots 
  VecBasis displac;
  readInBasis(displac, BasisId::STATE, BasisId::SNAPSHOTS, true, 0);
//  readAndProjectSnapshots(BasisId::STATE, converter6->vectorSize(), podBasis, converter6,
//                          snapshotCounts, timeStamps, displac);

  //Now build constraint vectors
  buildConstraintArray(displac);
}

void
DEIMConstraintSamplingDriver::computeInterpIndices(VecBasis &constraintBasis, std::vector<int> &maskIndices) {
  //member function for determining the sampling indicies for DEIM
  //result is (U*(P^T*U)^-1) where U is the left singular vectors of constraint snapshots and P
  //is the column selection matrix

  int constraintPodSizeMax = domain->solInfo().constraintPodSize;

  readInBasis(constraintBasis, BasisId::CONSTRAINT, BasisId::POD, false, 0, constraintPodSizeMax); // read in Constraint ROB specified under cpodrb 
  Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > constraintMatrix(constraintBasis.data(),constraintBasis.vectorSize(),constraintBasis.vectorCount());
  
  std::set<int> auxilaryIndices;//container for all dofs of selected nodes. Pseudo-GNAT type hyper-reduction
 
  int maxCoeffSlot;

  {
   Eigen::Matrix<double,Eigen::Dynamic,1> firstCol(constraintMatrix.rows());
   firstCol = constraintMatrix.col(0);
   //take absolute value of first column
   for(int row = 0; row != firstCol.rows(); row++)
       firstCol(row) = firstCol(row)*firstCol(row);
   //find maximum component
   firstCol.maxCoeff(&maxCoeffSlot);
   maskIndices.push_back(maxCoeffSlot);
  }

  //start loop to compute mask indicies 
  for(int i = 1; i < domain->solInfo().constraintPodSize; ++i){ //loop starts at 1 i.e. the 2nd column
    filePrint(stderr,"\r %4.2f%% complete", double(i)/double(domain->solInfo().constraintPodSize)*100.);

    //allocate space for P^T*U and P^T*u_i
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Umasked(maskIndices.size(),maskIndices.size());
    Eigen::Matrix<double,Eigen::Dynamic,1>              u_i_masked(maskIndices.size());

    for(int j = 0; j < maskIndices.size(); ++j) {//select proper rows and columns of constraint basis
      Umasked.row(j) = constraintMatrix.block(maskIndices[j],0,1,maskIndices.size()); //(P^T*U) is square
      u_i_masked(j) = constraintMatrix(maskIndices[j],maskIndices.size());//mask next column over
    }

    Eigen::FullPivLU< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > luOfUmasked(Umasked); //invert row reduced basis

    Eigen::Matrix<double,Eigen::Dynamic,1> residual(constraintBasis.vectorSize());

    if(luOfUmasked.isInvertible())
      residual = constraintMatrix.col(maskIndices.size()) - constraintMatrix.leftCols(maskIndices.size())*luOfUmasked.inverse()*u_i_masked;
    else
      throw std::runtime_error("... Matrix Not Invertible ...");

    //take absolute value of residual vector
    for(int row = 0; row != residual.rows(); row++)
      residual(row) = residual(row)*residual(row);

    residual.maxCoeff(&maxCoeffSlot); //find indice of maximum component

    maskIndices.push_back(maxCoeffSlot);
     
  }
  filePrint(stderr,"\r %4.2f%% complete\n", 100.);

  Eigen::Map< Eigen::Matrix<int,Eigen::Dynamic,Eigen::Dynamic> > indSol(maskIndices.data(),maskIndices.size(),1);
  std::cout << "selected indices of Constraint Function: \n" << indSol.transpose() << std::endl;
}


void
DEIMConstraintSamplingDriver::computeAndWriteDEIMBasis(VecBasis &constraintBasis, std::vector<int> &maskIndices)
{
  //member function for computing and writing the DEIM basis P*(P^T*U)^-T*U^T*V
  //where V is the mass-orthogonal POD basis, U is the left Singular vectors of the constraint snapshots
  //and P column selection matrix derived form the sampled indicies computed above
  int maxDeimBasisSize = domain->solInfo().maxDeimBasisSize;//if this is less than the number of indices, then we solve a least squares problem
  if(maxDeimBasisSize == 0)
    maxDeimBasisSize = maskIndices.size();

  // need indices to be ordered so that multiplication is done correctly
  std::set<int> indiceSet;
  for(std::vector<int>::iterator it = maskIndices.begin(); it !=  maskIndices.end(); ++it )
    indiceSet.insert(*it);

  deimBasis.dimensionIs(dualBasis.vectorCount(),dualBasis.vectorInfo());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dualMap(dualBasis.data(),dualBasis.size(), dualBasis.vectorCount());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > constraintMap(constraintBasis.data(),constraintBasis.size(), constraintBasis.vectorCount());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > deimMap(deimBasis.data(),deimBasis.size(), deimBasis.vectorCount());
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> compressedDBTranspose(dualBasis.numVectors(),maskIndices.size());
  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> rowReducedFM(maskIndices.size(),maxDeimBasisSize);

  //initialize (P^T*U)
  int row = 0; 
  for(std::set<int>::iterator it = indiceSet.begin(); it !=  indiceSet.end(); ++it ){
    rowReducedFM.row(row) = constraintMap.block(*it,0,1,maxDeimBasisSize); 
    row++;
  }

  //compute W^T = V^T*U*(P^T*U)^-1
  Eigen::JacobiSVD< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > SVDOfUmasked(rowReducedFM,Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix<double,Eigen::Dynamic,1> invSVs(SVDOfUmasked.nonzeroSingularValues());

  for(int i = 0; i != invSVs.rows(); i++) invSVs(i) = 1.0/SVDOfUmasked.singularValues()(i);

  std::cout << " ... condition Number of (P^T*U) = " << SVDOfUmasked.singularValues()(0)/SVDOfUmasked.singularValues()(SVDOfUmasked.nonzeroSingularValues()-1) << " ... " << std::endl;
  std::cout << " ... ||(P^T*U)^-1)||_2 = " << invSVs(SVDOfUmasked.nonzeroSingularValues()-1) << " ... " << std::endl;
  std::cout << " ... Sigma_m+1 = " << MPOSingularValue << " ... " << std::endl;
  std::cout << " ... E(f) = " << invSVs(SVDOfUmasked.nonzeroSingularValues()-1)*MPOSingularValue << " ... " << std::endl;

  compressedDBTranspose = dualMap.transpose()*constraintMap.leftCols(maxDeimBasisSize)*SVDOfUmasked.matrixV()*invSVs.asDiagonal()*SVDOfUmasked.matrixU().transpose();
  //we are computing the transpose of the basis

  std::cout << "compressed Basis" << std::endl;
  std::cout << compressedDBTranspose.transpose() << std::endl;

  //initialize deim basis to all zeros
  for(int i = 0; i != deimMap.rows(); i++)
    for(int j = 0; j != deimMap.cols(); j++)
      deimMap(i,j) = 0.0;

  //also open filestream for compressed deim basis
  FileNameInfo fileInfo;
  std::string fileName = BasisFileId(fileInfo, BasisId::DUALSTATE, BasisId::ROB); // filename is set under OUTPUT command robdataf
  fileName += ".compressed";
  VecNodeDof1Conversion converterReduced(compressedDBTranspose.cols());
  BasisOutputStream<1> output(fileName, converterReduced, false);
  std::cout << "output file vector size = " << output.vectorSize() << std::endl;

  std::cout << " ... writing compressed basis file to " << fileName.c_str() << " ... " << std::endl;
  //fist output the compressed basis
  for(int col = 0; col != compressedDBTranspose.rows(); ++col){
    Eigen::Matrix<double,Eigen::Dynamic,1> buffer(compressedDBTranspose.cols());
    buffer = compressedDBTranspose.row(col); // eigen is stored column major, use buffer to stream data correctly
    output << std::make_pair(double(col)+1.0,buffer.data());
  }

  //then output the expanded basis
  //expand rows W^T*P^T
  std::cout << " num of rows = " << compressedDBTranspose.rows() << " num of cols = " << compressedDBTranspose.cols() << std::endl;
  row = 0;
  for(std::set<int>::iterator it = indiceSet.begin(); it !=  indiceSet.end(); ++it ){
    deimMap.row(*it) = compressedDBTranspose.col(row);//use .col member to get rows of transposed basis
    row++;
  }

  std::vector<double> dummySVs; 
  writeBasisToFile(deimBasis, dummySVs, BasisId::DUALSTATE, BasisId::ROB,false);
}

void
DEIMConstraintSamplingDriver::writeSampledMesh(std::vector<int> &maskIndices) {
// keep only constraint functions that correspond to selected indices

  FileNameInfo fileInfo;
  const std::ios_base::openmode mode = std::ios_base::out;
  std::ofstream meshOut(getMeshFilename(fileInfo).c_str(), mode);
  meshOut << "TOPOLOGY" << std::endl;
  meshOut << "*constraints needed for DEIM" << std::endl; 

  Elemset &constraintEset = geoSource->getPackedEsetConstraintElementIeq();

  // need indices to be ordered so that multiplication is done correctly
  std::set<int> indiceSet;
  for(std::vector<int>::iterator it = maskIndices.begin(); it !=  maskIndices.end(); ++it )
    indiceSet.insert(*it);

  int counter = 1; 
  for(std::set<int>::iterator it = indiceSet.begin(); it !=  indiceSet.end(); ++it ){ // loop over the selected indices
    int numLMPC = 0;
    for(int i=0; i<geoSource->getNumConstraintElementsIeq(); ++i) { //loop over constraint functions to find which element the current indice belongs to
      Element *ele = constraintEset[i];
      int n = ele->getNumMPCs(); // number of constraints specified by this element
      int nn = ele->numNodes();  // number of nodes in this constraint element
      int *nvec = ele->allNodes();  // vector of node number -1 
      int tnum = ele->getElementType(); // element type
      int glnum = ele->getGlNum();
      for(int j = 0; j < n; ++j) {
        if(*it == numLMPC){
          meshOut << domain->numElements() + counter << " "; // output new element number
          meshOut << tnum << " ";                            // then output element type
          for( int w = 0;  w < nn - 1; ++w){
             meshOut << nvec[w] + 1 << " ";                  // lastly the element nodes
          }
          meshOut << std::endl;
          counter++;
        }     
        numLMPC++; 
      }
    }
  }
  meshOut << "ATTRIBUTES" << std::endl;
  meshOut << domain->numElements() + 1 << " " << domain->numElements() + counter - 1 << " 2" << std::endl;

}

void
DEIMConstraintSamplingDriver::buildConstraintArray(const VecBasis &displac)
{//this member function is for converting state snapshots to constraint function snapshots in the absence of precollected snapshots from model I

  // set up constraint functions
  if(geoSource->getNumConstraintElementsIeq() > 0) {

    // if lagrange multiplier snapshots are included, use them to filter out stateshots
    VecBasis *lagrangeMP = NULL;
    if(domain->solInfo().statePodRomFile.size() > 2){
      lagrangeMP = new VecBasis;
      std::cout << " ... Selective Snapshot Procedure ...\n";
      readInBasis(*lagrangeMP, BasisId::STATE, BasisId::SNAPSHOTS, false, 2);
    }


    // set up file to dump constraint snapshots into
    FileNameInfo fileInfo;
    std::string fileName = BasisFileId(fileInfo, BasisId::CONSTRAINT, BasisId::SNAPSHOTS); // filename is set under OUTPUT command constvct
    BasisOutputStream<1> output(fileName, *converter1, false);

    if(lagrangeMP){
      if(lagrangeMP->numVectors() != displac.numVectors()){
        std::cout << " ... Lagrange snapshots do not correspond to state snapshots ... " << std::endl;
        std::cout << "lagrangeMP->numVectors() = " << lagrangeMP->numVectors() << std::endl;
        std::cout << "displac.numVectors() = " << displac.numVectors() << std::endl;
     }
    }

    std::cout << " ... writing constraint snapshots to " <<  fileName.c_str() << " ... " << std::endl;

    // allocate space for constraint snapshot
    Eigen::Matrix<double,Eigen::Dynamic,1> ConstraintSnap(output.vectorSize());
  
    int iSnap = 0;

    for(int jSnap = 0; jSnap != displac.numVectors(); ++iSnap, ++jSnap){ // loop over snapshots
      bool notZeroColumn = true;
      if(lagrangeMP) { // if we have read in the lagrange multipliers, see if the current column has any active constraints
        if((*lagrangeMP)[jSnap].norm() < 1e-16){
          notZeroColumn = false;
        }
      }

      if(notZeroColumn){
        filePrint(stderr,"\r %4.2f%% complete", double(iSnap)/double(displac.numVectors())*100.);
        geomState->explicitUpdate(domain->getNodes(), displac[iSnap]);

        Elemset &eset = geoSource->getPackedEsetConstraintElementIeq();
        ResizeArray<LMPCons *> lmpc(0);
        int numLMPC = 0;
        for(int i=0; i<geoSource->getNumConstraintElementsIeq(); ++i) { //loop over constraint functions
          Element *ele = eset[i];
          static_cast<MpcElement*>(ele)->update(geomState, *geomState, domain->getNodes(), double(jSnap));
          int n = ele->getNumMPCs();
          LMPCons **l = ele->getMPCs();
          for(int j = 0; j < n; ++j) {
            lmpc[numLMPC++] = l[j];
          }
          delete [] l;
        }
        ConstraintSnap.setZero();
        for(int i=0; i<numLMPC; ++i) { // get right hand side which is Constraint evaluated at current state
          ConstraintSnap[i] = -1.0*lmpc[i]->rhs.r_value;
        }
        output << std::make_pair(1.0,ConstraintSnap.data()); // dump constraint to snapshot file
      }
    }
    filePrint(stderr,"\r %4.2f%% complete\n", 100.);

  } else {

    std::cout << "No constraints specified\n" << std::endl; 

  }

}

void DEIMConstraintSamplingDriver::filterRows(VecBasis &constraintBasis) {
  // zero out rows of constraint snapshots that correspond to zeros of the dual basis
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ConstraintMap(constraintBasis.data(), constraintBasis.size(), constraintBasis.numVectors());
  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > dualMap(dualBasis.data(), dualBasis.size(), dualBasis.numVectors());

  std::cout << " ... Filtering Constraint Snapshots ... " << std::endl;
  for (int row = 0; row < dualMap.rows(); ++row){
    if(dualMap.row(row).norm() <= 1e-16) // row is all zeros
      ConstraintMap.row(row).setZero();
  }

}

void DEIMConstraintSamplingDriver::OrthoConstraintSnap(VecBasis &constraintBasis,std::vector<double> &SVs)
{
  if(domain->solInfo().positiveElements){ // use nonnegative matrix factorization to create the constraint basis 
    std::cout << " ... Constructing Positive Constraint Subspace ... " << std::endl;
    // set up nonnegative matrix solver 
    NonnegativeMatrixFactorization solver(std::min(constraintBasis.numVectors(), constraintBasis.size()), domain->solInfo().use_nmf);
    solver.maxIterIs(domain->solInfo().nmfMaxIter);
    solver.toleranceIs(domain->solInfo().nmfTol);
    solver.numRandInitIs(domain->solInfo().nmfRandInit);
    int maxBasisDimension = domain->solInfo().maxSizePodRom + (domain->solInfo().nmfDelROBDim)*(domain->solInfo().nmfNumROBDim-1);
    solver.matrixSizeIs(constraintBasis.size(), constraintBasis.numVectors());
    solver.robSizeIs(constraintBasis.size(), maxBasisDimension);
    // read in vectors to solver
    int colCounter = 0;
    for (int iCol = 0; iCol < constraintBasis.numVectors(); ++iCol){
      double *buffer = solver.matrixCol(colCounter);
      std::copy(constraintBasis[iCol].data(),constraintBasis[iCol].data()+constraintBasis.size(),buffer);
      colCounter++;
    }
  
    // loop over increasing basis size
    int finalDim; 
    for (int iBasis=0; iBasis < domain->solInfo().nmfNumROBDim; ++iBasis) {
      int orthoBasisDim = domain->solInfo().maxSizePodRom + iBasis*domain->solInfo().nmfDelROBDim;
      filePrint(stderr, " ... Computation of a positive basis of size %d ...\n", orthoBasisDim);
      solver.basisDimensionIs(orthoBasisDim);
      if (iBasis==0)
        solver.solve(0);
      else
        solver.solve(orthoBasisDim-domain->solInfo().nmfDelROBDim);
      finalDim = orthoBasisDim;
    }
   
    // copy back into container
    for (int iVec = 0; iVec < finalDim; ++iVec) {
     constraintBasis[iVec] = solver.robCol(iVec);
     SVs.push_back(1.0); // dummy singular values 
    }

  } else { // user standard pod
#ifdef USE_EIGEN3
    std::cout << " ... Orthogonalizing Constraint Snapshots ..." << std::endl;
    Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ConstraintMap(constraintBasis.data(), constraintBasis.size(), constraintBasis.numVectors());
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> filteredConstraintMat; 
    if(domain->solInfo().filterSnapshotRows){ // only do svd on non-zero rows of constraint snapshot matrix
      std::cout << " ... Masked SVD ... " << std::endl;
      int nnzrows = 0;
      for(int row = 0; row < ConstraintMap.rows(); ++row){ // count the nonzero rows
        if(ConstraintMap.row(row).norm() >  1e-16) //non-zero
          nnzrows++;
      }

      // allocate space for dense section
      SVs.resize(std::min(constraintBasis.numVectors(),nnzrows));
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > SingularValueMap(SVs.data(),SVs.size());
      filteredConstraintMat.resize(nnzrows,ConstraintMap.cols());

      int counter = 0;
      for(int row = 0; row < ConstraintMap.rows(); ++row){ // get non-zero rows
        if(ConstraintMap.row(row).norm() > 1e-16){
          filteredConstraintMat.row(counter) = ConstraintMap.row(row);
          counter++;
        }
      }

      Eigen::JacobiSVD<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ConstraintSVD(filteredConstraintMat, Eigen::ComputeThinU);

      // expand the orthogonal basis
      counter = 0;
      for(int row = 0; row < ConstraintMap.rows(); ++row){
        if(ConstraintMap.row(row).norm() > 1e-16){
          ConstraintMap.block(row,0,1,SVs.size()) = ConstraintSVD.matrixU().row(counter);
          counter++;
        }
      }

      SingularValueMap = ConstraintSVD.singularValues();
      std::cout << " ... Orthogonal basis is size " << ConstraintSVD.matrixU().cols() << " by " << ConstraintSVD.matrixU().rows() << " ... " << std::endl;
    } else { // do svd on all rows
      SVs.resize(std::min(constraintBasis.numVectors(),constraintBasis.size()));
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,1> > SingularValueMap(SVs.data(),SVs.size());
      Eigen::JacobiSVD<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > ConstraintSVD(ConstraintMap, Eigen::ComputeThinU);
      ConstraintMap.block(0,0,ConstraintSVD.matrixU().rows(),ConstraintSVD.matrixU().cols()) = ConstraintSVD.matrixU();
      SingularValueMap = ConstraintSVD.singularValues();
      std::cout << " ... Orthogonal basis is size " << ConstraintSVD.matrixU().cols() << " by " << ConstraintSVD.matrixU().rows() << " ... " << std::endl;
    }
#endif
  }
}

//member functions for fascilitating the reading of snapshot data
int 
DEIMConstraintSamplingDriver::snapSize(BasisId::Type type, std::vector<int> &snapshotCounts)
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
DEIMConstraintSamplingDriver::readAndProjectSnapshots(BasisId::Type type, const int vectorSize, VecBasis &podBasis,
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
        expand(podBasis, reduce(podBasis, snapshot, podComponents), config[offset+count]);
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

} /* end namespace Rom */
#endif

Rom::DriverInterface *deimConstraintSamplingDriverNew(Domain *domain) {
#ifdef USE_EIGEN3
  return new Rom::DEIMConstraintSamplingDriver(domain);
#else
  std::cerr << " *** ERROR: DEIMConstraintSamplingDriver requires AERO-S configured with the Eigen library. Exiting...\n";
  exit(-1);
#endif
}
