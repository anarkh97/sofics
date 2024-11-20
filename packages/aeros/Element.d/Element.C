#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <numeric>

#include <Element.d/Element.h>
#include <Math.d/matrix.h>
#include <Math.d/Vector.h>
#include <Math.d/FullSquareMatrix.h>
#include <Utils.d/pstress.h>
extern std::map<int, double> weightList;
extern std::map<int, double> fieldWeightList;

double Element::weight() const
{
	auto it = weightList.find(getElementType());
	double weight = 1.0;
	double trueWeight = 1.0;
	if(it != weightList.end())
		weight = it->second;

	auto it2 = fieldWeightList.find((int)getCategory());
	if(it2 != fieldWeightList.end()) {
		double ratio =
			it2->second/std::accumulate(fieldWeightList.begin(), fieldWeightList.end(),
				0, [](double x, const std::pair<int, double>& y) {
					return x + y.second;
				});
		weight *= ratio;
	}
	return weight;
//	std::map<int, double>::iterator it1 = weightList.find(etype);
//	if(it1 != weightList.end()) {
//		ele->setWeight(it1->second);
//		ele->setTrueWeight(it1->second);
//	}
//
//	// adjust weight using FWEI if defined
//	if(!fieldWeightList.empty()) {
//		std::map<int,double>::iterator it2 = fieldWeightList.find((int)ele->getCategory());
//		if(it2 != fieldWeightList.end()) {
//			double weight = ele->weight();
//			double trueWeight = ele->trueWeight();
//			double ratio = it2->second/std::accumulate(fieldWeightList.begin(), fieldWeightList.end(), 0, weight_add());
//			ele->setWeight(weight*ratio);
//			ele->setTrueWeight(trueWeight*ratio);
//		}
//	}
}

double Element::trueWeight() const
{
	return weight();
}

void
Element::setCompositeData(int, int, double*, double*, double*)
{
  fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for non composite element type %d\n", getElementType());
}

double *
Element::setCompositeData2(int, int, double*, double*, CoordSet&, double)
{
  fprintf(stderr," *** WARNING: Attempting to define composite attributes\n"
                 "              for non composite element type %d\n", getElementType());
  return 0;
}

void
Element::getCFrame(CoordSet& cs, double cFrame[3][3]) const
{
  cFrame[0][0] = cFrame[1][1] = cFrame[2][2] = 1.;
  cFrame[0][1] = cFrame[0][2] = cFrame[1][0] = cFrame[1][2] = cFrame[2][0] = cFrame[2][1] = 0.;
}

void
Element::getVonMisesInt(CoordSet&, Vector&, double&, double&, int,
			double&, double&, double* dT)
{
  assert(0);
}

double
Element::weight(CoordSet& cs, double *gravityAcceleration)
{
  if(prop == NULL || gravityAcceleration == NULL) return 0.0;

  double mass = getMass(cs);
  double gravAccNorm = sqrt(gravityAcceleration[0]*gravityAcceleration[0] +
                            gravityAcceleration[1]*gravityAcceleration[1] +
                            gravityAcceleration[2]*gravityAcceleration[2]);
  return mass*gravAccNorm;
}

double
Element::getWeightThicknessSensitivity(CoordSet& cs, double *gravityAcceleration)
{
  if(prop == NULL || gravityAcceleration == NULL) return 0.0;

  double massSensitivity = getMassThicknessSensitivity(cs);
  double gravAccNorm = sqrt(gravityAcceleration[0]*gravityAcceleration[0] +
                            gravityAcceleration[1]*gravityAcceleration[1] +
                            gravityAcceleration[2]*gravityAcceleration[2]);
  return massSensitivity*gravAccNorm;
}

void 
Element::getWeightNodalCoordinateSensitivity(Vector& dwdx, CoordSet& cs, double *gravityAcceleration) 
{ 
  dwdx.zero();
  fprintf(stderr," *** WARNING: getWeightNodalCoordinateSensitivity is not implemented for element type %d\n", getElementType());
}

void
Element::getVonMises(Vector &stress, Vector &weight,
		     CoordSet &, Vector &, int, int, double*,
		     double, double, int)
{
  if(!isConstraintElement() && !isSpring())
    fprintf(stderr," *** WARNING: getVonMises not implemented for element type %d\n", getElementType());
  stress.zero();
  weight.zero();
}

extern "C" {
  void _FORTRAN(vmelmvc)(DComplex*, const int &, const int &, const int &, const int &, const int &);
  void _FORTRAN(strainvmc)(DComplex*, const int &, const int &, const int &, const int &);
}

void
Element::getVonMises(ComplexVector &stress, Vector &weight,
                     CoordSet &cs, ComplexVector &disp, int strInd, int surface,
                     double *ndTemps, double ylayer, double zlayer, int avgnum)
{
 // split displacement into real & imaginary vectors to compute independently the real
 // and imaginary components of the normal & shear stresses/strains
 // these can then be combined into a complex vector & used to compute the von mises stress or strain
 Vector realDisp(disp.size());
 for(int i=0; i<disp.size(); ++i) realDisp[i] = disp[i].real();
 Vector imagDisp(disp.size());
 for(int i=0; i<disp.size(); ++i) imagDisp[i] = disp[i].imag();

 if( ((strInd >= 0) && (strInd <=5)) || ((strInd >= 7) && (strInd <=12)) ) { // not von mises
   Vector realStress(stress.size(), 0.0);
   getVonMises(realStress, weight, cs, realDisp, strInd, surface, ndTemps, ylayer, zlayer, avgnum);
   Vector imagStress(stress.size(), 0.0);
   getVonMises(imagStress, weight, cs, imagDisp, strInd, surface, ndTemps, ylayer, zlayer, avgnum);
   for(int i=0; i<stress.size(); ++i) stress[i] = DComplex(realStress[i], imagStress[i]);
 }
 else if((strInd==6) || (strInd==13)) { // von mises
   if((ylayer!=0) || (zlayer!=0)) {
     fprintf(stderr," *** WARNING: complex getVonMises not implemented for ylayer %f zlayer %f \n",ylayer,zlayer);
     stress.zero();
     weight.zero();
     return;
   }
   int maxgus = numNodes();
   int si = (strInd == 6) ? 0 : 1;
   FullM realStress(maxgus,9); realStress.zero();
   getAllStress(realStress, weight, cs, realDisp, si, surface, ndTemps);
   FullM imagStress(maxgus,9); imagStress.zero();
   getAllStress(imagStress, weight, cs, imagDisp, si, surface, ndTemps);
   GenFullM<DComplex> complexStress(maxgus,7);  complexStress.zero();
   for(int i=0; i<maxgus; ++i)
     for(int j=0; j<6; ++j) complexStress[i][j] = DComplex(realStress[i][j], imagStress[i][j]);

   if(strInd == 6) _FORTRAN(vmelmvc)(complexStress.data(),maxgus,7,1,1,maxgus);
   else _FORTRAN(strainvmc)(complexStress.data(),maxgus,7,1,maxgus);

   for(int i=0; i<maxgus; ++i) stress[i] = complexStress[i][6];
 }
 else {
   if(!isConstraintElement() && !isSpring())
     fprintf(stderr," *** WARNING: complex getVonMises not implemented for strInd %d\n",strInd);
   stress.zero();
   weight.zero();
 }
}

void
Element::getAllStress(FullM &stress, Vector &weight, CoordSet &cs,
                      Vector &elDisp, int strInd, int surface,
                      double *ndTemps)
{
  if(isConstraintElement() || isSpring()) {
    stress.zero();
    weight.zero();
  }
  else {
    // Get components of symmetric stress or strain tensors one at a time by calling Element::getVonMises.
    Vector s(numNodes());
    for(int i=0; i<6; ++i) {
      int index = (strInd == 0) ? i : i+7;
      getVonMises(s, weight, cs, elDisp, index, surface, ndTemps);
      for(int j=0; j<numNodes(); ++j) stress[j][i] = s[j];
    }

    // Get element principal stress or strain for each node
    double svec[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    double pvec[3] = {0.0,0.0,0.0};
    for(int i=0; i<numNodes(); ++i) {
      for(int j=0; j<6; ++j) {
        svec[j] = stress[i][j];
      }
      // Convert engineering to tensor strains
      if(strInd != 0) {
        svec[3] /= 2;
        svec[4] /= 2;
        svec[5] /= 2;
      }
      pstress(svec,pvec);
      for(int j=0; j<3; ++j) {
        stress[i][j+6] = pvec[j];
      }
    }
  }
}

void
Element::getAllStress(FullMC &stress, Vector &weight,
                      CoordSet &cs, ComplexVector &disp, int strInd, int surface,
                      double *ndTemps)
{
 // split displacement into real & imaginary vectors to compute independently the real
 // and imaginary components of the normal & shear stresses/strains
 // these can then be combined into a complex matrix
 Vector realDisp(disp.size());
 for(int i=0; i<disp.size(); ++i) realDisp[i] = disp[i].real();
 Vector imagDisp(disp.size());
 for(int i=0; i<disp.size(); ++i) imagDisp[i] = disp[i].imag();

 FullM realStress(numNodes(),9); realStress.zero();
 getAllStress(realStress, weight, cs, realDisp, strInd, surface, ndTemps);
 FullM imagStress(numNodes(),9); imagStress.zero();
 getAllStress(imagStress, weight, cs, imagDisp, strInd, surface, ndTemps);
 for(int i=0; i<numNodes(); ++i) {
   for(int j=0; j<6; ++j) stress[i][j] = DComplex(realStress[i][j], imagStress[i][j]);

   // Get Element Principals for each node without averaging
   complex<double> svec[6];
   for (int j=0; j<6; ++j) svec[j] = stress[i][j];
   // Convert Engineering to Tensor Strains
   if(strInd != 0) {
     svec[3] /= 2;
     svec[4] /= 2;
     svec[5] /= 2;
   }
   complex<double> pvec[3];
   pstress(svec,pvec);
   for (int j=0; j<3; ++j) {
     stress[i][j+6] = pvec[j];
   }
 }
}

PrioInfo Element::examine(int sub, MultiFront *mf)
{
  fprintf(stderr," *** ERROR: Element type %d cannot be decomposed since examine function is not implemented \n", getElementType());
  PrioInfo p;
  p.isReady = false;
  return p;
}

void
Element::getGravityForce(CoordSet&, double *, Vector &force, int, GeomState *)
{
  if(!isConstraintElement() && !isSpring() && getCategory() != Element::Thermal)
    fprintf(stderr," *** WARNING: Gravity force not implemented for element (%d), type %d\n", getGlNum()+1, getElementType());
  force.zero();
}

void
Element::getGravityForceThicknessSensitivity(CoordSet&, double *, Vector &forceSen, int, GeomState *)
{
  forceSen.zero();
}

void
Element::getGravityForceNodalCoordinateSensitivity(CoordSet& cs, double *gravityAcceleration,
                                                   GenFullM<double> &dGfdx, int gravflg, GeomState *geomState)
{
  if(!isConstraintElement() && !isSpring() && getCategory() != Element::Thermal)
    fprintf(stderr," *** WARNING: Gravity force sensitivity not implemented for element (%6d), type %3d\n", getGlNum()+1, getElementType());
  dGfdx.zero();
}

void
Element::getThermalForce(CoordSet&, Vector &, Vector &force, int glflag,
                         GeomState *)
{
  force.zero();
}

void
Element::getThermalForceAdj(CoordSet& cs, Vector &ndH, Vector &b, int glflag)
{
  // assuming that the Element::getThermalForce is a linear function, i.e. force = A*ndT,
  // compute the action of the adjoint, i.e. ndH = A.transpose()*b

  // default implementation:
  Vector xi(numNodes());
  Vector Ai(numDofs());
  for(int i=0; i<numNodes(); ++i) {
    xi.zero();
    xi[i] = 1;
    getThermalForce(cs, xi, Ai, glflag);
    ndH[i] = Ai*b;
  }
}

void
Element::computeHeatFluxes(Vector &heatflux, CoordSet&, Vector &, int)
{
  if(!isConstraintElement() && !isSpring())
    fprintf(stderr," *** WARNING: Heat Fluxes not implemented for element type %d\n", getElementType());
  heatflux.zero();
}

void
Element::computeSloshDisp(Vector &fluidDispSlosh, CoordSet&, Vector &, int)
{
  if (getElementType() != 302) {
    fprintf(stderr," *** WARNING: Fluid Displacements not implemented for element type %d\n", getElementType());
    fluidDispSlosh.zero();
  }
}

void
Element::computeSloshDispAll(Vector &fluidDispSlosh, CoordSet&, Vector &)
{
  if (getElementType() != 302) {
    fprintf(stderr," *** WARNING: Fluid Displacements not implemented for element type %d\n", getElementType());
    fluidDispSlosh.zero();
  }
}

void
Element::getIntrnForce(Vector &elForce, CoordSet&, double *, int,double *)
{
  //if(!isConstraintElement() && !isSpring())
  //  fprintf(stderr," *** WARNING: Internal force not implemented for element type %d\n", getElementType());
  elForce.zero();
}

int
Element::getFace(int iFace, int *fn)
{
  fprintf(stderr," *** WARNING: getFace not implemented for element type %d\n", getElementType());
  return 0;
}

void
Element::computePressureForce(CoordSet&, Vector& elPressureForce,
                              GeomState *, int cflg, double time)
{
  if(!isConstraintElement() && !isSpring())
    fprintf(stderr," *** WARNING: Pressure force not implemented for element type %d\n", getElementType());
  elPressureForce.zero();
}

int
Element::getTopNumber() const
{
 return -1;
}

void
Element::addFaces(PolygonSet *)
{
  fprintf (stderr, "\n\n ERROR: addFaces not implemented for this element !!!"
                   " \n   exiting... \n\n");
  exit (1);
}

void
Element::setMaterial(NLMaterial *)
{
  fprintf(stderr, " *** WARNING: Trying to use Material on unsupported element!\n");
}

Corotator *
Element::getCorotator(CoordSet &, double *, int, int)
{
  fprintf(stderr, " *** WARNING: Corotator not implemented for element %d\n", glNum+1);
  return 0;
}

void
Element::getCG(CoordSet &cset, double &xcg, double &ycg, double &zcg)
{
  int *myNodes;
  int preAllocNodes[125];
  int nnd = numNodes()-numInternalNodes();
  if(numNodes() > 125)
    myNodes = new int[numNodes()];
  else
    myNodes = preAllocNodes;

  nodes(myNodes);

  xcg = ycg = zcg = 0.0;

  int i;
  for(i = 0; i < nnd; ++i) {
    Node *n = cset[myNodes[i]];
    xcg += n->x;
    ycg += n->y;
    zcg += n->z;
  }
  xcg /= nnd;
  ycg /= nnd;
  zcg /= nnd;
  if(numNodes() > 125)
    delete[] myNodes;
}

// END DEC

FullSquareMatrix Element::stiffness(const CoordSet&, double* kel, int cmflg) const
{
  FullSquareMatrix result(1, kel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

FullSquareMatrix Element::massMatrix(const CoordSet&, double* mel, int cmflg) const
{
  FullSquareMatrix result(1, mel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

FullSquareMatrixC
Element::stiffness(const CoordSet&, complex<double> *d) const
{ 
  return FullSquareMatrixC();
}

FullSquareMatrixC
Element::massMatrix(const CoordSet&, complex<double> *d) const
{
  return FullSquareMatrixC();
}

FullSquareMatrix
Element::dampingMatrix(CoordSet& cs, double *m, int cmflg)
{
  FullSquareMatrix ret(1,m);
  ret.setSize(numDofs());
  ret.zero();
  return ret;
}

FullSquareMatrix
Element::imStiffness(CoordSet& cs, double *m, int cmflg)
{
  fprintf(stderr, " *** WARNING: Element imStiffness not implmented for this element\n");
  FullSquareMatrix ret(4,m);
  ret.zero();
  return ret;
}


void Element::aRubberStiffnessDerivs(CoordSet& cs, complex<double> *d, int n,
                                     double omega) {
  fprintf(stderr, " *** WARNING: Element aRubberStiffnessDerivs not implemented for this element\n");
}

void Element::lumpMatrix(FullSquareMatrix& m, std::vector<double>& factors) const
{
  double MM = 0.0; // total mass of element
  double MD = 0.0; // mass of the diagonal
  const int dim = m.dim();
  for(int i = 0; i < dim; ++i) {
    MD += m[i][i];
    for(int j = 0; j < i; ++j) { 
      MM += m[i][j]; 
    }
  }
  MM = MD + 2.0*MM;
  if (MD > 0) {
    const double factor = MM/MD;
    for(int i = 0; i < dim; ++i) {
      m[i][i] *= factor;
      factors.push_back(m[i][i]/MM); // PJSA store this for getGravityForce
      for(int j = 0; j < i; ++j) {
        m[i][j] = m[j][i] = 0.0;
      }
    }
  }
  return;
}

// mratio = 1.0 for consistent
// mratio = 0.0 for lumped
FullSquareMatrix Element::massMatrix(const CoordSet& cs, double* m, double mratio) const
{
  if(mratio == 0.0 && getMassType() == 1) { // in this case, get the consistent mass matrix and lump it using diagonal scaling
    FullSquareMatrix result = massMatrix(cs, m, 1);
    std::vector<double> factors;
    lumpMatrix(result, factors);
    return result;
  }
  else return massMatrix(cs, m, int(mratio));
}

FullSquareMatrixC Element::complexStiffness(CoordSet&, DComplex* kel, int cmflg)
{
  FullSquareMatrixC result(1, kel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

FullSquareMatrixC Element::complexDampingMatrix(CoordSet&, DComplex* cel, int cmflg)
{
  FullSquareMatrixC result(1, cel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

FullSquareMatrixC Element::complexMassMatrix(CoordSet&, DComplex* mel, double mratio)
{
  FullSquareMatrixC result(1, mel);
  result.setSize(numDofs());
  result.zero();
  return result;
}

#include <Element.d/Helm.d/HelmElement.h>

bool Element::isFluidElement() { return dynamic_cast<HelmElement *>(this); }

double
Element::computeStabilityTimeStep(FullSquareMatrix &K, FullSquareMatrix &M, CoordSet &cs, GeomState *gs,
                                  double stable_tol, int stable_maxit)
{
  if(prop && (prop->rho != 0 || getMass(cs) != 0)) {

      using std::sqrt;
      using std::abs;
      double eigmax;
      double relTol    = stable_tol; // stable_tol default is 1.0e-3
      double preeigmax = 0.0;

      int numdofs = K.dim();
      int maxIte  = stable_maxit; // stable_maxit default is 100

      Vector v(numdofs);
      Vector z(numdofs);

// Starts from an arbitrary array.
      int i,j;
      for (i=0; i<numdofs; ++i)
        v[i] = (double) (i+1) / (double) numdofs;

// Power iteration loop

      for (i=0; i<maxIte; ++i) {
        z.zero();
        K.multiply(v,z,1);

        for (j=0; j< numdofs; ++j)
          z[j] /= M[j][j];

// Normalize

        double zmax = z[0];
        for (j=1; j< numdofs; ++j)
          if (abs(z[j])>zmax) zmax = abs(z[j]);

        eigmax = zmax;

        v = (1.0/zmax)*z;

        if ( abs(eigmax - preeigmax) < relTol*abs(preeigmax) ) break;

        preeigmax = eigmax;
      }

      // compute stability maximum time step
      double sdt = 2.0 / sqrt(eigmax);

      return sdt;
  }
  else { // phantom or other element without mass
      return std::numeric_limits<double>::infinity();
  }
}

void
Element::getStiffnessThicknessSensitivity(CoordSet&, FullSquareMatrix &dStiffdThick, int)
{
  dStiffdThick.zero();
}

void
Element::getStiffnessNodalCoordinateSensitivity(FullSquareMatrix *&dStiffdx, CoordSet &cs)
{
  for(int i=0; i<numNodes()*3; ++i) dStiffdx[i].zero();
  fprintf(stderr," *** WARNING: getStiffnessNodalCoordinateSensitivity is not implemented for element type %d\n", getElementType());
}

void
Element::getVonMisesThicknessSensitivity(Vector &dStdThick, Vector &weight, CoordSet&, Vector&,
                                         int, int, double *, int, double, double)
{
  weight = 1;
  dStdThick.zero();
}

void
Element::getVonMisesThicknessSensitivity(ComplexVector &dStdThick, ComplexVector &weight, CoordSet&,
                                         ComplexVector&, int, int, double *, int, double, double)
{
  weight = DComplex(1,0);
  dStdThick.zero();
}

void
Element::getVonMisesDisplacementSensitivity(GenFullM<double> &dStdDisp, Vector &weight, GenFullM<double> *, CoordSet&,
                                            Vector&, int, int, double *, int, double, double)
{
  dStdDisp.zero();
  weight.zero();
  fprintf(stderr," *** WARNING: getVonMisesDisplacementSensitivity is not implemented for element type %d\n", getElementType());
}

void
Element::getVonMisesDisplacementSensitivity(GenFullM<DComplex> &dStdDisp, ComplexVector &weight, GenFullM<DComplex> *,
                                            CoordSet&, ComplexVector&, int, int, double *,
                                            int, double, double)
{
  dStdDisp.zero();
  weight.zero();
  fprintf(stderr," *** WARNING: getVonMisesDisplacementSensitivity is not implemented for element type %d\n", getElementType());
}

void
Element::getVonMisesNodalCoordinateSensitivity(GenFullM<double> &dStdx, Vector &weight, CoordSet&, Vector&,
                                               int, int, double *, int, double, double)
{
  dStdx.zero();
  weight.zero();
  fprintf(stderr," *** WARNING: getVonMisesNodalCoordinateSensitivity is not implemented for element type %d\n", getElementType());
}
