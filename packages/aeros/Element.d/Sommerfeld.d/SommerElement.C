#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <Utils.d/dbg_alloca.h>
#include <Element.d/Sommerfeld.d/SommerElement.h>
#include <Element.d/Helm.d/HelmElement.h>
#include <Driver.d/Domain.h>
#include <Corotational.d/GeomState.h>

bool SommerElement::first = true; //HB
#include <Utils.d/DistHelper.h>

void
SommerElement::renum(const int *table) {
	int i;
	int *nn = getNodes();

	for (i = 0; i < numNodes(); ++i)
		nn[i] = table[nn[i]];
}

void
SommerElement::renum(EleRenumMap &table) {
	int i;
	int *nn = getNodes();

	for (i = 0; i < numNodes(); ++i)
		nn[i] = table[nn[i]];
}


SommerElement *SommerElement::clone() {
	fprintf(stderr, " *** WARNING: SommerElement::clone not implemented\n");
	return 0;
}

int *
SommerElement::nodes(int *p) const {
	if (p == 0) p = new int[numNodes()];
	auto nn = getNodes();
	int i;
	for (i = 0; i < numNodes(); i++) {
		p[i] = nn[i];
	}
	return p;
}

void SommerElement::flipNormal() {
	int *nds = (int *) dbg_alloca(sizeof(int) * numNodes());
	nodes(nds);
	int *nn = getNodes();
	int i;
	int j = numNodes();
	for (i = 0; i < numNodes(); i++) {
		nn[i] = nds[--j];
	}
}


void
SommerElement::findBothEle(const Connectivity *nodeToElem, int *eleTouch,
                           int *eleCount, int myNum, const Elemset *eset, int *ie) {
	ie[0] = ie[1] = -1;
	int i, iele;
	int nNodes = numNodes();
	int *nn = getNodes();
	for (i = 0; i < nNodes; i++) {
		for (iele = 0; iele < nodeToElem->num(nn[i]); iele++) {
			int eleNum = (*nodeToElem)[nn[i]][iele];
			if (eset) {
				if ((*eset)[eleNum]->isSommerElement()) continue;
				if ((*eset)[eleNum]->isFsiElement()) continue;
			}

			if (eleTouch[eleNum] != myNum) {
				eleTouch[eleNum] = myNum;
				eleCount[eleNum] = 1;
			} else {
				eleCount[eleNum]++;
				if (eleCount[eleNum] == nNodes) {
					if (ie[0] == -1) ie[0] = eleNum;
					else {
						ie[1] = eleNum;
						return;
					}
				}
			}
		}
	}
	fprintf(stderr, "\n");
}


int
SommerElement::findAndSetBothEle(const CoordSet &cs, Elemset &eset,
                                 const Connectivity *nodeToElem, int *eleTouch, int *eleCount, int myNum) {
	int ie[2];
	findBothEle(nodeToElem, eleTouch, eleCount, myNum, &eset, ie);

	if (ie[0] == -1) {
		fprintf(stderr, "SommerElement::findAndSetBothEle could not find the"
				" corresponding element 1.\n");
		return 0;
	}
	el = eset[ie[0]];
	if (ie[1] == -1) {
		fprintf(stderr, "SommerElement::findAndSetBothEle did not find the"
				" corresponding element 2.\n");
	}
	if (ie[1] != -1) el2 = eset[ie[1]];

	HelmElement *he = dynamic_cast<HelmElement *>(el);
	int isFluid = 1;
	if (he == 0) isFluid = 0;
//  fprintf(stderr,"SommerElement::findAndSetBothEle he %d\n",isFluid);

	if (he) soundSpeed = el->getProperty()->soundSpeed;
	else if (ie[1] != -1) soundSpeed = el2->getProperty()->soundSpeed;
	else {
		fprintf(stderr, "SommerElement::findAndSetBothEle could not set ss.\n");
		return 0;
	}


	int *nodes = (int *) alloca(sizeof(int) * el->numNodes());
	el->nodes(nodes);
	double x, y, z;
	x = y = z = 0.0;
	int iEle;
	for (iEle = 0; iEle < el->numNodes(); iEle++) {
		Node nd = cs.getNode(nodes[iEle]);
		x += nd.x;
		y += nd.y;
		z += nd.z;
	}
	x = x / el->numNodes();
	y = y / el->numNodes();
	z = z / el->numNodes();

	double normal[3];
	getNormal(cs, normal);

	int *sNodes = getNodes();
	Node nd = cs.getNode(sNodes[0]);
	double d = -1.0 * (normal[0] * nd.x + normal[1] * nd.y + normal[2] * nd.z);

	if (normal[0] * x + normal[1] * y + normal[2] * z + d > 0.0) {
		if (isFluid) return 1;
		else return -1;
	} else {
		if (isFluid) return -1;
		else return 1;
	}
}


int
SommerElement::findEle(const Connectivity *nodeToElem, int *eleTouch,
                       int *eleCount, int myNum, Elemset *eset, int it) {
	int i, iele;
	int nNodes = numNodes();
	int *nn = getNodes();
	for (i = 0; i < nNodes; i++) {
		for (iele = 0; iele < nodeToElem->num(nn[i]); iele++) {
			int eleNum = (*nodeToElem)[nn[i]][iele];
      if(iele > 0 && eleNum == (*nodeToElem)[nn[i]][iele-1]) continue;
			if (eset) {
				if ((*eset)[eleNum]->isSommerElement()) continue;
				if ((*eset)[eleNum]->isFsiElement()) continue;
// RT
				HelmElement *he = dynamic_cast<HelmElement *>((*eset)[eleNum]);
				if (it == 1) {
					if (he == 0) continue;
					if (!he->isFluid()) continue;
				}
// RT
				if (it == 2) {
					if (he != 0) continue;
				}
			}

			if (eleTouch[eleNum] != myNum) {
				eleTouch[eleNum] = myNum;
				eleCount[eleNum] = 1;
			} else {
				eleCount[eleNum]++;
				if (eleCount[eleNum] == nNodes) {
					iEle = eleNum;
					return eleNum;
				}
			}
		}
	}
	return -1;
}

int
SommerElement::findAndSetEle(const CoordSet &cs, Elemset &eset,
                             const Connectivity *nodeToElem, int *eleTouch, int *eleCount, int myNum, int it) {
	iEle = findEle(nodeToElem, eleTouch, eleCount, myNum, &eset, it);
	if (iEle == -1) {
		fprintf(stderr, "SommerElement::findAndSetEle could not find the"
				" corresponding element.\n");
		return 0;
	}
	el = eset[iEle];
	soundSpeed = el->getProperty()->soundSpeed;

	int *nodes = (int *) dbg_alloca(sizeof(int) * el->numNodes());
	el->nodes(nodes);
	double x, y, z;
	x = y = z = 0.0;
	for (iEle = 0; iEle < el->numNodes(); iEle++) {
		Node nd = cs.getNode(nodes[iEle]);
		x += nd.x;
		y += nd.y;
		z += nd.z;
	}
	x = x / el->numNodes();
	y = y / el->numNodes();
	z = z / el->numNodes();

	double normal[3];
	getNormal(cs, normal);

	int *sNodes = getNodes();
	Node nd = cs.getNode(sNodes[0]);
	double d = -1.0 * (normal[0] * nd.x + normal[1] * nd.y + normal[2] * nd.z);

	if (normal[0] * x + normal[1] * y + normal[2] * z + d > 0.0) {
		return 1;
	} else {
		return -1;
	}
}

void
SommerElement::center(CoordSet &cs, double *c) {
	//returns the cordinates of an element in c - assumption: 3D problem
	int size = (this)->numNodes();
	int *nn = getNodes();
	int i;
	c[0] = 0;
	c[1] = 0;
	c[2] = 0;
	for (i = 0; i < size; i++) {
		Node nd = cs.getNode(nn[i]);
		c[0] += nd.x;
		c[1] += nd.y;
		c[2] += nd.z;
	}
	c[0] = c[0] / size;
	c[1] = c[1] / size;
	c[2] = c[2] / size;
}

void
SommerElement::checkAndSetNormal(CoordSet &cs) {
	int *nodes = (int *) dbg_alloca(sizeof(int) * el->numNodes());
	el->nodes(nodes);
	double x, y, z;
	x = y = z = 0.0;
	int i;
	for (i = 0; i < el->numNodes(); i++) {
		Node nd = cs.getNode(nodes[i]);
		x += nd.x;
		y += nd.y;
		z += nd.z;
	}
	x = x / el->numNodes();
	y = y / el->numNodes();
	z = z / el->numNodes();
	double normal[3];
	getNormal(cs, normal);

	int *sNodes = getNodes();
	Node nd = cs.getNode(sNodes[0]);
	double d = -1.0 * (normal[0] * nd.x + normal[1] * nd.y + normal[2] * nd.z);

	if (normal[0] * x + normal[1] * y + normal[2] * z + d > 0.0) {
		flipNormal();
	}
}

FullSquareMatrix SommerElement::sommerMatrix(CoordSet &cs) const {
	int nd = numDofs();
	return sommerMatrix(cs, new double[nd * nd]);
}

FullSquareMatrix SommerElement::turkelMatrix(CoordSet &cs) const {
	int nd = numDofs();
	return turkelMatrix(cs, new double[nd * nd]);
}

FullSquareMatrix SommerElement::refinedSommerMatrix(CoordSet &cs) {
	int nd = numDofs();
	return sommerMatrix(cs, new double[nd * nd]);
}

/*
// not used anymore - developped for the ABC time domain - JFD
FullSquareMatrix SommerElement::surfStiffMatrix(CoordSet &cs)
{
 int nd = numDofs();
 return surfStiffMatrix(cs,new double [nd*nd]);
}
*/

FullSquareMatrix SommerElement::HSommerMatrix(const CoordSet &cs) const {
	int nd = numDofs();
	return HSommerMatrix(cs, new double[nd * nd]);
}

/*
// not used anymore - developped for ABC time domain - JFD
FullSquareMatrix SommerElement::HKSommerMatrix(CoordSet &cs)
{
 int nd = numDofs();
 return HKSommerMatrix(cs,new double [nd*nd]);
}
*/

FullSquareMatrix SommerElement::interfMatrixConsistent(CoordSet &cs) {
	int nd = numDofs();
	return interfMatrixConsistent(cs, new double[nd * nd]);
}

FullSquareMatrix SommerElement::interfMatrixLumped(CoordSet &cs) {
	int nd = numDofs();
	return interfMatrixLumped(cs, new double[nd * nd]);
}

FullRectMatrix SommerElement::transferMatrix(CoordSet &cs, double *d) {

	fprintf(stderr, "SommerElement::transferMatrix should not be called.\n");
	FullRectMatrix M(numDofs(), numDofs(), d);
	return M;
}

FullSquareMatrix SommerElement::turkelMatrix(CoordSet &cs, double *d) const {

	fprintf(stderr, "SommerElement::turkelMatrix should not be called.\n");
	FullSquareMatrix M(numDofs(), d);
	return M;
}

FullSquareMatrix SommerElement::refinedSommerMatrix(CoordSet &cs, double *d) {

	//fprintf(stderr,"SommerElement::refinedSommerMatrix should not be called.\n");
	//FullSquareMatrix M(numDofs(),d);
	//return M; //JTG
	return sommerMatrix(cs, d);
}

/*
// not used anymore - developped for ABC time domain - JFD
FullSquareMatrix SommerElement::surfStiffMatrix(CoordSet &cs, double *d) {

 fprintf(stderr,"SommerElement::surfStiffMatrix should not be called.\n");
 FullSquareMatrix M(numDofs(),d);
 return M;
}
*/

FullSquareMatrix SommerElement::HSommerMatrix(const CoordSet &cs, double *d) const {

	fprintf(stderr, "SommerElement::HSommerMatrix should not be called.\n");
	fprintf(stderr, "Note that if the ABC is a flat surface, ATDARB==0.0 is equivalent to ATDARB=1.0\n");
	fprintf(stderr, "as ATDARB=1.0 involves the curvature that is 0 on a flat surface.\n");
	FullSquareMatrix M(numDofs(), d);
	return M;
}

/*
// not used anymore - developped for ABC time domain - JFD
FullSquareMatrix SommerElement::HKSommerMatrix(CoordSet &cs, double *d) {

 fprintf(stderr,"SommerElement::HKSommerMatrix should not be called.\n");
 FullSquareMatrix M(numDofs(),d);
 return M;
}
*/

FullSquareMatrix SommerElement::interfMatrixConsistent(CoordSet &cs, double *d) {

	sommerMatrix(cs, d);
// double kappa = el->getProperty()->kappaHelm;
// for (int i=0; i<numDofs()*numDofs(); i++) d[i] *= kappa;
	FullSquareMatrix M(numDofs(), d);

	return M;
}


FullSquareMatrix SommerElement::interfMatrixLumped(CoordSet &cs, double *d) {

	sommerMatrix(cs, d);
	FullSquareMatrix M(numDofs(), d);
// double kappa = el->getProperty()->kappaHelm;
	for (int i = 0; i < numDofs(); i++) {
		double dd = 0.0;
		for (int j = 0; j < numDofs(); j++) {
			dd += M[i][j];
			M[i][j] = 0.0;
		}
//   M[i][i] = kappa*dd;
	}
	return M;
}

void SommerElement::neumVector(CoordSet &cs, Vector &cv, int pflag, GeomState *, double) {

	if (numDofs() != numNodes()) {
		fprintf(stderr, "SommerElement::neumVector called for a nonscalar element.\n"
				"Terminating.\n");
		exit(-1);
	}

	int i;
	double *x = (double *) dbg_alloca(sizeof(double) * numNodes());
	double *y = (double *) dbg_alloca(sizeof(double) * numNodes());
	double *z = (double *) dbg_alloca(sizeof(double) * numNodes());

	for (i = 0; i < numNodes(); i++) {
		Node nd = cs.getNode(getNode(i));
		x[i] = nd.x;
		y[i] = nd.y;
		z[i] = nd.z;
	}
	double *dmatrix = (double *) dbg_alloca(sizeof(double) * numDofs() * numDofs());
	FullSquareMatrix matrix = sommerMatrix(cs, dmatrix);
	for (i = 0; i < numDofs(); i++) {
		int j;
		cv[i] = 0.0;
		for (j = 0; j < numDofs(); j++) {
			cv[i] -= matrix[i][j];// sommerMatrix is negative definite
		}
	}
}

void SommerElement::neumVectorJacobian(CoordSet &cs, FullSquareMatrix &kel, int pflag, GeomState *, double) {

	kel.zero();
}

void SommerElement::neumVector(CoordSet &cs, ComplexVector &cv, double k,
                               double dx, double dy, double dz, int pflag) {

	if (numDofs() != numNodes()) {
		fprintf(stderr, "SommerElement::neumVector called for a nonscalar element.\n"
				"Terminating.\n");
		exit(-1);
	}

	int i;
	double *x = (double *) dbg_alloca(sizeof(double) * numNodes());
	double *y = (double *) dbg_alloca(sizeof(double) * numNodes());
	double *z = (double *) dbg_alloca(sizeof(double) * numNodes());

	for (i = 0; i < numNodes(); i++) {
		Node nd = cs.getNode(getNode(i));
		x[i] = nd.x;
		y[i] = nd.y;
		z[i] = nd.z;
	}

	double normal[3];
	getNormal(cs, normal);

	double *dmatrix = (double *) dbg_alloca(sizeof(double) * numDofs() * numDofs());
	FullSquareMatrix matrix = sommerMatrix(cs, dmatrix);

	ComplexD *dudn = (ComplexD *) dbg_alloca(sizeof(ComplexD) * numDofs());
	//HB:to compute the source term due to the incident wave in coupled pb (working with total field)
	bool isCoupled = (dom) ? dom->solInfo().isCoupled : false;
	ComplexD tmp = (isCoupled) ? ComplexD(0.0, k * (normal[0] * dx + normal[1] * dy + normal[2] * dz) - k)
	                           : ComplexD(0.0, k * (normal[0] * dx + normal[1] * dy + normal[2] * dz));
	for (i = 0; i < numNodes(); i++) {
		dudn[i] = tmp * exp(ComplexD(0.0, k * (dx * x[i] + dy * y[i] + dz * z[i])));
	}

	for (i = 0; i < numDofs(); i++) {
		int j;
		cv[i] = ComplexD(0.0, 0.0);
		for (j = 0; j < numDofs(); j++) {
			// RT: sommerMatrix is negative definite
			cv[i] -= matrix[i][j] * dudn[j];
		}
	}

}

void SommerElement::neumVectorDeriv(CoordSet &cs, ComplexVector &cv, double k,
                                    double dx, double dy, double dz, int n,
                                    int pflag) {
	if (numDofs() != numNodes()) {
		fprintf(stderr, "SommerElement::neumVectorDeriv called for a nonscalar element.\n"
				"Terminating.\n");
		exit(-1);
	}
	if (pflag != 0) fprintf(stderr, "Point source not implemented in SommerElement::neumVectorDeriv .\n");

	int i;
	double *x = (double *) dbg_alloca(sizeof(double) * numNodes());
	double *y = (double *) dbg_alloca(sizeof(double) * numNodes());
	double *z = (double *) dbg_alloca(sizeof(double) * numNodes());

	for (i = 0; i < numNodes(); i++) {
		Node nd = cs.getNode(getNode(i));
		x[i] = nd.x;
		y[i] = nd.y;
		z[i] = nd.z;
	}

	double normal[3];
	getNormal(cs, normal);

	double *dmatrix = (double *) dbg_alloca(sizeof(double) * numDofs() * numDofs());
	FullSquareMatrix matrix = sommerMatrix(cs, dmatrix);

	ComplexD *dudn = (ComplexD *) dbg_alloca(sizeof(ComplexD) * numDofs());

	// f = bk e^(ak)
	// f^(1) = b( e^(ak) + ak e^(ak) ) = b(1+ak)e^(ak)
	// f^(2) = b(a e^(ak) + (1+ak)a e^(ak)) = ba(2+ak)e^(ak)
	// f^(3) = ba(a + (2+ak)a) e^(ak) = ba^2(3+ak)e^(ak)

	bool isCoupled = (dom) ? dom->solInfo().isCoupled : false;
	ComplexD b = (isCoupled) ? ComplexD(0.0, normal[0] * dx + normal[1] * dy + normal[2] * dz - 1.0)
	                         : ComplexD(0.0, normal[0] * dx + normal[1] * dy + normal[2] * dz);

	for (i = 0; i < numNodes(); i++) {
		ComplexD a = ComplexD(0.0, dx * x[i] + dy * y[i] + dz * z[i]);
		dudn[i] = b * pow(a, n - 1) * (double(n) + a * k) * exp(a * k);
	}

	for (i = 0; i < numDofs(); i++) {
		int j;
		cv[i] = ComplexD(0.0, 0.0);
		for (j = 0; j < numDofs(); j++) {
			// RT: sommerMatrix is negative definite
			cv[i] -= matrix[i][j] * dudn[j];
		}
	}
}


/* 
void SommerElement::neumVector(CoordSet& cs, ComplexVector& cv, double k,
                              double dx, double dy, double dz, int order) {
 
 if (numDofs()!=numNodes()) {
   fprintf(stderr,"SommerElement::neumVector called for a nonscalar element.n"
                  "Terminating.\n");
   exit(-1);
 }

 int i;
 double *x = (double*)dbg_alloca(sizeof(double)*numNodes());
 double *y = (double*)dbg_alloca(sizeof(double)*numNodes());
 double *z = (double*)dbg_alloca(sizeof(double)*numNodes());

 for(i=0; i<numNodes(); i++) {
   Node nd = cs.getNode(getNode(i));
   x[i] = nd.x;
   y[i] = nd.y;
   z[i] = nd.z;
 }

 double normal[3];
 getNormal(cs,normal);
 
 double *dmatrix = (double*)dbg_alloca(sizeof(double)*numDofs()*numDofs());
 FullSquareMatrix matrix = sommerMatrix(cs,dmatrix);

 double fact = 1.0;
 for(i=2;i<=order;i++) fact *= double(i);

 ComplexD *dudn = (ComplexD*)dbg_alloca(sizeof(ComplexD)*numDofs());
 ComplexD ind(0.0,normal[0]*dx+normal[1]*dy+normal[2]*dz);
 for(i=0; i<numNodes(); i++) {
   ComplexD ixd(0.0,dx*x[i]+dy*y[i]+dz*z[i]);
   if (order==0)
     dudn[i] = ind*k * exp( k*ixd );
   else if (order==1)
     dudn[i] = ind * ( k*ixd + ComplexD(order,0.0) ) * exp(k*ixd);
   else
     dudn[i] = ind * ( k*pow(ixd,order) + double(order)*pow(ixd,order-1) ) *
               exp(k*ixd)/fact;
 }

 for(i=0;i<numDofs();i++) {
   int j;
   cv[i] = ComplexD(0.0,0.0);
   for(j=0;j<numDofs();j++) {
// RT: sommerMatrix is negative definite
     cv[i] -= matrix[i][j] * dudn[j];
   }
 }
}
*/


void SommerElement::sommerVector(CoordSet &cs, ComplexVector &cv,
                                 ComplexVector &cf) {

	fprintf(stderr, "SommerElement::sommerVector not implemented.\n");
}


void SommerElement::btVector(CoordSet &cs, ComplexVector &cv,
                             ComplexVector &cf) {

	fprintf(stderr, "SommerElement::btVector not implemented.\n");
}


void SommerElement::ffpNeum(int ndir, ComplexD *ffp,
                            CoordSet &cs, ComplexD *u,
                            double k, double (*dir)[3], double *idir) {

	int nNodes = numNodes();
	int nDofs = numDofs();
	if (nDofs != nNodes) {
		fprintf(stderr, "SommerElement::ffpDir called for a nonscalar element.\n"
				"Terminating.\n");
		exit(-1);
	}

	int i;
	double *x = (double *) dbg_alloca(sizeof(double) * nNodes);
	double *y = (double *) dbg_alloca(sizeof(double) * nNodes);
	double *z = (double *) dbg_alloca(sizeof(double) * nNodes);

	for (i = 0; i < nNodes; i++) {
		Node nd = cs.getNode(getNode(i));
		x[i] = nd.x;
		y[i] = nd.y;
		z[i] = nd.z;
	}

	double normal[3];
	getNormal(cs, normal);

	double *dmatrix = (double *) dbg_alloca(sizeof(double) * nDofs * nDofs);
	FullSquareMatrix matrix = sommerMatrix(cs, dmatrix);

	ComplexD tmp = ffpCoef(k);

	double &dx = idir[0], &dy = idir[1], &dz = idir[2];
	ComplexD *dudn = (ComplexD *) dbg_alloca(nDofs * sizeof(ComplexD));
	ComplexD tmp2(0.0, k * (normal[0] * dx + normal[1] * dy + normal[2] * dz));
	for (i = 0; i < nNodes; i++)
		dudn[i] = exp(ComplexD(0.0, k * (dx * x[i] + dy * y[i] + dz * z[i]))) * tmp2;


	ComplexD *solPart = (ComplexD *) dbg_alloca(nDofs * sizeof(ComplexD));
	ComplexD *expPart = (ComplexD *) dbg_alloca(nDofs * sizeof(ComplexD));

	for (int iDir = 0; iDir < ndir; iDir++) {
		double &dirx = dir[iDir][0], &diry = dir[iDir][1], &dirz = dir[iDir][2];
		for (i = 0; i < numNodes(); ++i) {
			solPart[i] = ComplexD(0.0, k *
			                           (dirx * normal[0] + diry * normal[1] + dirz * normal[2])) * u[i] + dudn[i];
			expPart[i] = exp(-ComplexD(0.0, k * (dirx * x[i] + diry * y[i] + dirz * z[i])));
		}
		int j;
		for (i = 0; i < nDofs; ++i)
			for (j = 0; j < nDofs; ++j)
				// RT: sommerMatrix returns negative coefficients
				ffp[iDir] -= tmp * solPart[i] * matrix[i][j] * expPart[j];
	}
}


void SommerElement::ffpDir(int ndir, ComplexD *ffp,
                           CoordSet &cs,
                           ComplexD *u, ComplexD *dudn,
                           double k, double (*dir)[3], double *idir) {

	fprintf(stderr, "Old SommerElement::ffpDir not implemented.\n");
}


void SommerElement::ffpDir(int ndir, ComplexD *ffp,
                           CoordSet &cs, ComplexD *u,
                           double k, double (*dir)[3], double *idir) {

	int nNodes = numNodes();
	int nDofs = numDofs();
	if (nDofs != nNodes) {
		fprintf(stderr, "SommerElement::ffpDir called for a nonscalar element.\n"
				"Terminating.\n");
		exit(-1);
	}

	int i;
	double *x = (double *) dbg_alloca(sizeof(double) * nNodes);
	double *y = (double *) dbg_alloca(sizeof(double) * nNodes);
	double *z = (double *) dbg_alloca(sizeof(double) * nNodes);

	for (i = 0; i < nNodes; i++) {
		Node nd = cs.getNode(getNode(i));
		x[i] = nd.x;
		y[i] = nd.y;
		z[i] = nd.z;
	}

	double normal[3];
	getNormal(cs, normal);


	double *dmatrix = (double *) dbg_alloca(sizeof(double) * nDofs * nDofs);
	FullSquareMatrix matrix = sommerMatrix(cs, dmatrix);

	ComplexD tmp = ffpCoef(k);

	ComplexD *solPart = (ComplexD *) dbg_alloca(nDofs * sizeof(ComplexD));
	ComplexD *expPart = (ComplexD *) dbg_alloca(nDofs * sizeof(ComplexD));

	for (int iDir = 0; iDir < ndir; iDir++) {
		double &dirx = dir[iDir][0], &diry = dir[iDir][1], &dirz = dir[iDir][2];
		for (i = 0; i < nDofs; ++i) {
			solPart[i] = ComplexD(0.0, k *
			                           (dirx * normal[0] + diry * normal[1] + dirz * normal[2])) * u[i];
			expPart[i] = exp(-ComplexD(0.0, k * (dirx * x[i] + diry * y[i] + dirz * z[i])));
		}
		int j;
		for (i = 0; i < nDofs; ++i)
			for (j = 0; j < nDofs; ++j)
				// RT: sommerMatrix returns negative coefficients
				ffp[iDir] -= tmp * solPart[i] * matrix[i][j] * expPart[j];
	}

}

void SommerElement::BT2(CoordSet &cs, double *e, double *f, double *g,
                        double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d) {
	fprintf(stderr, "SommerElement::BT2 not implemented.\n");
}

void SommerElement::BT2n(CoordSet &cs, double *e, double *f, double *g,
                         double (*tau1)[3], double (*tau2)[3], double k, ComplexD *d, int n) {
	fprintf(stderr, "SommerElement::BT2n not implemented.\n");
}

void SommerElement::sphereBT2(CoordSet &cs, double r, double k, ComplexD *d) {
	fprintf(stderr, "SommerElement::sphereBT2 not implemented.\n");
}

void SommerElement::ellipsoidBT2(CoordSet &cs, double a, double b, double k, ComplexD *d) {
	fprintf(stderr, "SommerElement::ellipsoidBT2 not implemented.\n");
}

void SommerElement::sommerMatrixEllipsoid(CoordSet &cs, double kappa, double H[3], double K[3], ComplexD *d) {
	fprintf(stderr, "SommerElement::sommerMatrixEllipsoid not implemented.\n");
}

double SommerElement::getSize(CoordSet &cs) {
	fprintf(stderr, "SommerElement::getSize not implemented.\n");
	return 0.0;
}

FullSquareMatrixC
SommerElement::turkelMatrix(CoordSet &, double, int) {
	fprintf(stderr, " *** WARNING: Attempting to use turkelMatrix\n");
	return FullSquareMatrixC(1);
}


FullSquareMatrixC
SommerElement::turkelMatrix(CoordSet &, double, int, DComplex *) {
	fprintf(stderr, " *** WARNING: Attempting to use turkelMatrix\n");
	return FullSquareMatrixC(1);
}


GenStackFSFullMatrix<double>
SommerElement::wetInterfaceMatrix(CoordSet &, double *) {
	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::wetInterfaceMatrix\n");
	return GenStackFSFullMatrix<double>(0, 0, 0);
}


void
SommerElement::wetInterfaceLMPC(CoordSet &, LMPCons *, int) {
	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::wetInterfaceLMPC\n");
}


void
SommerElement::ffp(CoordSet &cs, int numFFP, double *dirFFP,
                   complex<double> *sol, complex<double> *ffpv, bool direction) {

	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::ffp\n");
}

int SommerElement::numWetDofs() const {
	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::numWetDofs\n");
	return 0;
}

int SommerElement::numSolidDofs() const {
	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::numSolidDofs\n");
	return 0;
}

int *SommerElement::wetDofs(DofSetArray &, int *p) const {
	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::wetDofs\n");
	return 0;
}

int *SommerElement::solidDofs(DofSetArray &, int *p) const {
	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::solidDofs\n");
	return 0;
}

void SommerElement::wetInterfaceVector(CoordSet &cs, ComplexVector &cv,
                                       double kappa, double dx, double dy, double dz,
                                       int dero, int pflag) {
	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::wetInterfaceVector\n");
}


void SommerElement::wetInterfaceVector(CoordSet &cs, ComplexVector &cv,
                                       complex<double> (*diri)[3], complex<double> *coefi) {
	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::wetInterfaceVector\n");
}


void SommerElement::wetInterfaceVectorDeriv(CoordSet &cs, ComplexVector &cv,
                                            complex<double> (*diri)[3], complex<double> *coefi,
                                            complex<double> *kappa, int dero) {
	fprintf(stderr, " *** WARNING: Attempting to use SommerElement::wetInterfaceVectorVector\n");
}


int
SommerElement::dim() const {
	fprintf(stderr, " *** WARNING: SommerElement::dim not implemented\n");
	return 0;
}

ComplexD
SommerElement::ffpCoef(double k) const {
	fprintf(stderr, " *** WARNING: SommerElement::ffpCoef not implemented\n");
	return 0.0;
}

void
SommerElement::getNormal(const CoordSet &, double[3]) const {
	fprintf(stderr, " *** WARNING: SommerElement::getNormal not implemented\n");
}


int SommerElement::ludcmp(FullSquareMatrix A, int n, int *indx) const// LU decomposition
//form "Numerical recipes in C", cambridge
{
	int i, imax, j, k;
	double big, dum, sum, temp;
	double *vv = new double[n];
	for (i = 0; i < n; i++) vv[i] = 0.0;
	//*d=1.0;
	for (i = 0; i < n; i++) {
		big = 0.0;
		for (j = 0; j < n; j++)
			if ((temp = fabs(A[i][j])) > big) big = temp;
		if (big == 0.0) fprintf(stderr, "Singular Matrix in ludcmp-SommerElement.C\n");
		vv[i] = 1.0 / big;
	}
	for (j = 0; j < n; j++) {
		for (i = 0; i < j; i++) {
			sum = A[i][j];
			for (k = 0; k < i; k++) sum -= A[i][k] * A[k][j];
			A[i][j] = sum;
		}
		big = 0.0;
		for (i = j; i < n; i++) {
			sum = A[i][j];
			for (k = 0; k < j; k++)
				sum -= A[i][k] * A[k][j];
			A[i][j] = sum;
			if ((dum = vv[i] * fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if (j != imax) {
			for (k = 0; k < n; k++) {
				dum = A[imax][k];
				A[imax][k] = A[j][k];
				A[j][k] = dum;
			}
			//*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if (A[j][j] == 0.0) A[j][j] = 1.0e-10;
		if (j != (n - 1)) {
			dum = 1.0 / (A[j][j]);
			for (i = j + 1; i < n; i++) A[i][j] *= dum;
		}
	}
	delete[] vv;
	vv = NULL;
	return 0;
}

void SommerElement::lubksb(FullSquareMatrix A, int n, int *indx, double *b) const//LU factorisation
//form "Numerical recipes in C", cambridge
{
	int i, ii = 0, ip, j;
	double sum;
	for (i = 0; i < n; i++) {
		ip = indx[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii)
			for (j = ii - 1; j <= i - 1; j++)
				sum -= A[i][j] * b[j];
		else if (sum) ii = i + 1;
		b[i] = sum;
	}
	for (i = n - 1; i >= 0; i--) {
		sum = b[i];
		for (j = i + 1; j < n; j++) sum -= A[i][j] * b[j];
		b[i] = sum / A[i][i];
	}
}

void SommerElement::invert(FullSquareMatrix A, FullSquareMatrix B) const {
	int n = A.dim();
	//copy A in Acopy
	double *valA = new double[n * n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			valA[i + n * j] = A[i][j];
	FullSquareMatrix Acopy(n, valA);
	int *ind = new int[n];
	ludcmp(Acopy, n, ind);
	double *vect = new double[n];
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			vect[j] = 0.0;
		vect[i] = 1.0;
		lubksb(Acopy, n, ind, vect);
		for (int j = 0; j < n; j++)
			B[i][j] = vect[j];
	}
	delete[] valA;
	valA = NULL;
	delete[] ind;
	ind = NULL;
	delete[] vect;
	vect = NULL;
}

//Functions that depend on the order of the absorbing boundary condition
//functions dofs, markDofs and numDofs are defined in each element.
#include<Driver.d/Domain.h>

void
SommerElement::markDofs(DofSetArray &dsa) const {
	fprintf(stderr, "SommerElement::markDofs not implemented.\n");
}

FullSquareMatrix
SommerElement::massMatrix(const CoordSet &cs, double *d, int cmflg) const {
	//This function calls subfunctions and assembles the matrix using the submatrixes produced.
	double numD = dom->solInfo().ATDARBFlag;
	if (numD == 1.5) {
		std::cerr << "SommerElement::massMatrix, ATDARBFlag>1 does not work" << std::endl;
/*
    FullSquareMatrix ret(2*numNodes(),d);
    ret.zero();

    bool morecompl = false;//this flag is for mass, damping and stiffness matrixes
    if (morecompl) {//add bottom right block matrix = coef*tB
      double ss = el->getProperty()->ss;
      int nNo = numNodes();
      int i, j;
      //double coef = (dom->solInfo().newmarkBeta)*(dom->solInfo().dt)/(dom->solInfo().newmarkGamma);
      double coefa = 1.0;
      double coefb = 3/(16*coefa);
      FullSquareMatrix temp = surfStiffMatrix(cs);
      for (i = 0 ; i<nNo; i++)
        for (j = 0 ; j<nNo; j++)
          ret[i+nNo][j+nNo] = -coefb*ss/2*temp[j][i];
      temp.zero();
      temp = HKSommerMatrix(cs);
      for (i = 0 ; i<nNo; i++)
        for (j = 0 ; j<nNo; j++)
          ret[i+nNo][j+nNo] += coefb*ss/2*temp[j][i];// two - cancel
    }
    return ret;
*/
	} else {
		FullSquareMatrix ret(numNodes(), d);
		ret.zero();
		return ret;
	}
	return FullSquareMatrix(0);
}

FullSquareMatrix SommerElement::stiffness(const CoordSet &cs, double *Ks, int flg) const {
	double numD = dom->solInfo().ATDARBFlag;
	int nNo = numNodes();
	if (numD == 1.5) {
		std::cerr << "SommerElement::stiffness, ATDARBFlag>1 does not work" << std::endl;
/*
    double ss = el->getProperty()->ss;
    int i, j;
    FullSquareMatrix ret(2*nNo,Ks);
    ret.zero();
    FullSquareMatrix temp = HSommerMatrix(cs);
    for (i = 0 ; i<nNo; i++)
      for (j = 0 ; j<nNo; j++)
        ret[i][j] = temp[i][j];

    bool Complx = false;
    if (!Complx) {
      temp.zero();
      temp = surfStiffMatrix(cs);
      temp -= HKSommerMatrix(cs);
      bool SurfHKposdef =false;//make the HK matrix positive definite by moving the negative eignevalues
      if (SurfHKposdef) {
        double tempthr = 0.01;//threshold
        double *tempVa = (double*)alloca(sizeof(double)*nNo);
        double *tempd = (double*)alloca(sizeof(double)*nNo*nNo);
        double *temptemp = (double*)alloca(sizeof(double)*nNo*nNo);
        for (i = 0 ; i < nNo*nNo ; i++) {tempd[i] = 0.0; temptemp[i] = 0.0;}
        FullSquareMatrix tempD(nNo,tempd);
        FullSquareMatrix tempTemp(nNo,temptemp);
        temp.eigenV((double*)tempVa);//gets eigenvalues and eigenvectors (in temp)
        //set the diagonal matrix
        for (i = 0 ; i < nNo ; i++) {
          if (tempVa[i]<=0) {
            tempD[i][i] = tempthr;
            if (fabs(tempVa[i])>tempthr) fprintf(stderr," SommerElement.C: large change in eigenvalues\n");
          }
          else
            tempD[i][i] = tempVa[i];
        }
        //multiply the matrixes
        tempD.multiply(temp,tempTemp);
        invert(temp,tempD);
        tempD.multiply(tempTemp,temp);
    }
      temp *= ss/2;
      for (i = 0 ; i<nNo; i++)
        for (j = 0 ; j<nNo; j++) {
          ret[i+nNo][j] = temp[j][i];
          ret[i][j+nNo] = temp[i][j];
        }
    }

    bool morecompl1=false;//flag for the stiffness matrix only
    if (morecompl1) {//add the bottom right
      double coef = 1.0;
      for (i = 0 ; i < nNo ; i++)
        for (j = 0 ; j < nNo ; j++)
          ret[i+nNo][j+nNo] = -coef*ret[i+nNo][j];
    }
    bool morecompl=false;//flag for the mass, damping and stiffness matrixes
    if (morecompl) {//add the bottom right and top left block matrix
      //double coef = (dom->solInfo().newmarkBeta)*(dom->solInfo().dt)/(dom->solInfo().newmarkGamma);
      double coefa = 1.0;
      double coefb = 3/(16*coefa);
      for (i = 0 ; i < nNo ; i++)
        for (j = 0 ; j < nNo ; j++) {
          ret[i][j] = -coefb*ret[i][j+nNo];
          ret[i+nNo][j+nNo] = -coefa*ret[i+nNo][j];
        }
    }
    return ret;
*/
	} else if (numD == 1.0) {
		FullSquareMatrix ret(nNo, Ks);
		ret = HSommerMatrix(cs, Ks);
		return ret;
	} else {
		FullSquareMatrix ret(nNo, Ks);
		ret.zero();
		return ret;
	}
	return FullSquareMatrix(0);
}

FullSquareMatrix SommerElement::dampingMatrix(CoordSet &cs, double *Cs, int flg) {
	double numD = dom->solInfo().ATDARBFlag;
	double ss = el->getProperty()->ss;
	int nNo = numNodes();
	if (numD == 1.5) {
		std::cerr << "SommerElement::dampingMatrix ATDARBFlag>1 has not been implemented" << std::endl;
/*
    // this function is commented because it diverges ... 
    int i, j;
    FullSquareMatrix ret(2*nNo,Cs);
    ret.zero();
    FullSquareMatrix temp = refinedSommerMatrix(cs);
    temp /= -ss;
    for (i = 0 ; i<nNo; i++)
      for (j = 0 ; j<nNo; j++)
        ret[i][j] = temp[i][j];
    temp.zero();
    //temp = turkelMatrix(cs);//this is an option for comparison
    temp = surfStiffMatrix(cs);
    temp -= HKSommerMatrix(cs);
    bool SurfHKposdef =false;
    bool Complx = false;
    double tempcoef = (Complx) ? 1.0 : -1.0;
    if (SurfHKposdef) {
      double tempthr = 0.01;//threshold
      double *tempVa = (double*)alloca(sizeof(double)*nNo);
      double *tempd = (double*)alloca(sizeof(double)*nNo*nNo);
      double *temptemp = (double*)alloca(sizeof(double)*nNo*nNo);
      for (i = 0 ; i < nNo*nNo ; i++) {tempd[i] = 0.0; temptemp[i] = 0.0;}
      FullSquareMatrix tempD(nNo,tempd);
      FullSquareMatrix tempTemp(nNo,temptemp);
      temp.eigenV((double*)tempVa);//gets eigenvalues and eigenvectors (in temp)
      //set the diagonal matrix
      for (i = 0 ; i < nNo ; i++) {
        if (tempVa[i]<=0) {
          tempD[i][i] = tempthr;
          if (fabs(tempVa[i])>tempthr) fprintf(stderr," SommerElement.C: large change in eigenvalues\n");
        }
        else
          tempD[i][i] = tempVa[i];
      }
      //multiply the matrixes
      tempD.multiply(temp,tempTemp);
      invert(temp,tempD);
      tempD.multiply(tempTemp,temp);
    }
    temp *= ss/2;
    for (i = 0 ; i<nNo; i++)
      for (j = 0 ; j<nNo; j++)
        ret[i+nNo][j+nNo] = tempcoef*temp[j][i];//transpose

    bool morecompl = false;//flag for mass, damping and stiffness matrix
    if (morecompl) {//add top right and bottom left block matrix
      //double coef = (dom->solInfo().newmarkBeta)*(dom->solInfo().dt)/(dom->solInfo().newmarkGamma);
      double coefa = 1.0;
      double coefb = 3/(16*coefa);
      for (i = 0 ; i < nNo ; i++)
        for (j = 0 ; j < nNo ; j++) {
          ret[i][j+nNo] = -coefb*ret[j+nNo][i+nNo];//transpose
          ret[i+nNo][j] = -coefb*ret[i+nNo][j+nNo];
        }
    }
    return ret;
*/
	} else {
		FullSquareMatrix ret(nNo, Cs);
		/*ret = refinedSommerMatrix(cs,Cs); //JTG*/
		ret = sommerMatrix(cs,
		                   Cs); // PJSA 10/27/09 replaced call to refinedSommerMatrix with sommerMatrix after observing incorrect solution.
		// refinedSommerMatrix functions in TriangleSommerBC.C and QuadSommerBC.C need to be tested and debugged before using
		ret /= -ss;
		return ret;
	}
	return FullSquareMatrix(0);
}

/*

FullSquareMatrix SommerElement::imStiffness(CoordSet &cs, double *Ks, int flg)
{
  //this function is a test for stabilizing the high order absorbing boundary condition
//fprintf(stderr,"\nimStiffnessM\n");
  double numD = dom->solInfo().ATDARBFlag;
  double ss = el->getProperty()->ss;
  int nNo = numNodes();
  int i, j;
  if (numD==1.5) {
    FullSquareMatrix ret(2*nNo,Ks);
    ret.zero();
    FullSquareMatrix temp = surfStiffMatrix(cs);
    temp -= HKSommerMatrix(cs);
    bool SurfHKposdef =false;
    if (SurfHKposdef) {
      double tempthr = 0.01;//threshold
      double *tempVa = (double*)alloca(sizeof(double)*nNo);
      double *tempd = (double*)alloca(sizeof(double)*nNo*nNo);
      double *temptemp = (double*)alloca(sizeof(double)*nNo*nNo);
      for (i = 0 ; i < nNo*nNo ; i++) {tempd[i] = 0.0; temptemp[i] = 0.0;}
      FullSquareMatrix tempD(nNo,tempd);//temp^-1*tempD*temp > 0
      FullSquareMatrix tempTemp(nNo,temptemp);
      temp.eigenV((double*)tempVa);//gets eigenvalues and eigenvectors (in temp)
      //set the diagonal matrix
      for (i = 0 ; i < nNo ; i++) {
        if (tempVa[i]<=0) {
          tempD[i][i] = tempthr;
          if (fabs(tempVa[i])>tempthr) fprintf(stderr," SommerElement.C: large change in eigenvalues\n");
        }
        else
          tempD[i][i] = tempVa[i];
      }
      //multiply the matrixes
      tempD.multiply(temp,tempTemp);
      invert(temp,tempD);
      tempD.multiply(tempTemp,temp);
    }
    temp *= ss/2;
    for (i = 0 ; i<nNo; i++)
      for (j = 0 ; j<nNo; j++) {
        ret[i+nNo][j] = temp[j][i];//transpose
        ret[i][j+nNo] = temp[i][j];
        //ret[i+nNo][j+nNo] = 0.0;//transpose
      }
    return ret;
  }
  else {

    FullSquareMatrix ret(nNo,Ks);
    ret.zero();
    fprintf(stderr," SommerElement.C: not implemented for ATDARB<1.5\n");
    return ret;
  }
}
*/

