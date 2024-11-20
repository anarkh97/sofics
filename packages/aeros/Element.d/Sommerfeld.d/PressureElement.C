#include <Utils.d/dofset.h>
#include <Corotational.d/GeomState.h>
#include <Mortar.d/MortarDefines.h>
#include <iostream>

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::PressureElement(int *_nn, PressureBCond *_pbc)
		: nNodes(FaceElementType::NumberOfNodes), pbc(_pbc) {
	nn = new int[nNodes];
	for (int i = 0; i < nNodes; ++i)
		nn[i] = _nn[i];
	addTerms(DofSet::XYZdisp);
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
void
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::addTerms(DofSet nodalDofs) {
	terms.clear();
	nterms = 0;
	for (int i = 0; i < nNodes; ++i) {
		int j = 0;
		for (int k = 0;; ++k)
			if (nodalDofs.contains(1 << k)) {
				BCond t;
				t.setData(nn[i], k, 0.0);
				terms.push_back(t);
				nterms++;
				j++;
				if (j == nodalDofs.count()) break;
			}
	}
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::~PressureElement() {
	delete pbc;
	delete[] nn;
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
int
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::numNodes() const {
	return nNodes;
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
void
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::renum(const int *table) {
	for (int i = 0; i < numNodes(); ++i)
		if (nn[i] > -1)
			nn[i] = table[nn[i]];
	for (int i = 0; i < nterms; ++i)
		terms[i].nnum = table[terms[i].nnum];
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
void
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::renum(EleRenumMap &table) {
	for (int i = 0; i < numNodes(); ++i)
		if (nn[i] > -1)
			nn[i] = table[nn[i]];
	for (int i = 0; i < nterms; ++i)
		terms[i].nnum = table[terms[i].nnum];
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
int *
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::nodes(int *p) const {
	if (p == 0) p = new int[numNodes()];
	for (int i = 0; i < numNodes(); ++i) p[i] = nn[i];
	return p;
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
int
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::numDofs() const {
	return nterms;
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
int *
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::dofs(DofSetArray &dsa, int *p) const  {
	if (p == 0) p = new int[numDofs()];
	for (int i = 0; i < nterms; i++)
		dsa.number(terms[i].nnum, 1 << terms[i].dofnum, p + i);
	return p;
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
void
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::markDofs(DofSetArray &dsa) const {
	for (int i = 0; i < nterms; i++)
		dsa.mark(terms[i].nnum, 1 << terms[i].dofnum);
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
void
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>::neumVector(CoordSet &c0, Vector &F,
                                                                                             int, GeomState *c1,
                                                                                             double t) {
	// construct face element & local coordsets
	int nodes[FaceElementType::NumberOfNodes];
	typedef std::vector<NodeTemplate<double> *> CoordSetType;
	CoordSetType face_c0(FaceElementType::NumberOfNodes); // reference configuration
	CoordSetType face_c1(FaceElementType::NumberOfNodes); // current configuration
	for (int inode = 0; inode < FaceElementType::NumberOfNodes; ++inode) {
		nodes[inode] = inode;
		face_c0[inode] = new NodeTemplate<double>(c0[nn[inode]]->x, c0[nn[inode]]->y, c0[nn[inode]]->z);
		face_c1[inode] = (c1) ? new NodeTemplate<double>((*c1)[nn[inode]].x, (*c1)[nn[inode]].y, (*c1)[nn[inode]].z)
		                      : face_c0[inode];
	}
	FaceElementType *FaceEl = new FaceElementType(nodes);

	// quadrature rule
	int deg = (pbc->conwep && pbc->conwepswitch) ? VariableDegree : ConstantDegree;  // quadrature rule degree
	QuadratureRule c(deg);

	// local variables to be computed at each integration point
	double xi[2];                             // abscissa
	double weight;                            // weight
	double N[FaceElementType::NumberOfNodes]; // shape function values
	double normal[3];
	double Xi[3], Normal[3];                  // reference coordinates and unit normal
	double pressure = pbc->val;

	F.zero();

	// loop over the integration points
	for (int ip = 0; ip < c.getN(); ++ip) {

		// get the integration point abscissa and weight
		c.getAbscissaAndWeight(ip, xi, weight);

		// compute shape functions
		FaceEl->GetShapeFctVal(N, xi);

		// compute the surface normal
		FaceEl->GetIsoParamMappingNormalJacobianProduct(normal, xi, face_c1);

		if (pbc->conwep && pbc->conwepswitch) {
			FaceEl->GetUnitNormal(Normal, xi, face_c0);
			FaceEl->LocalToGlobalCoord(Xi, xi, face_c0);
			pressure = pbc->val + BlastLoading::ComputeGaussPointPressure(Xi, Normal, t, *(pbc->conwep));
		}

		for (int k = 0; k < FaceElementType::NumberOfNodes; ++k) {
			double tmp = N[k] * weight * pressure;
			F[3 * k + 0] += tmp * normal[0];
			F[3 * k + 1] += tmp * normal[1];
			F[3 * k + 2] += tmp * normal[2];
		}
	}

	delete FaceEl;
	for (int inode = 0; inode < FaceElementType::NumberOfNodes; ++inode) {
		delete face_c0[inode];
		if (c1) delete face_c1[inode];
	}
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
void
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>::neumVectorJacobian(CoordSet &c0,
                                                                                                     FullSquareMatrix &Ktan,
                                                                                                     int, GeomState *c1,
                                                                                                     double t) {
	// construct face element & local coordsets
	int nodes[FaceElementType::NumberOfNodes];
	typedef std::vector<NodeTemplate<double> *> CoordSetType;
	CoordSetType face_c0(FaceElementType::NumberOfNodes); // reference configuration
	CoordSetType face_c1(FaceElementType::NumberOfNodes); // current configuration
	for (int inode = 0; inode < FaceElementType::NumberOfNodes; ++inode) {
		nodes[inode] = inode;
		face_c0[inode] = new NodeTemplate<double>(c0[nn[inode]]->x, c0[nn[inode]]->y, c0[nn[inode]]->z);
		face_c1[inode] = (c1) ? new NodeTemplate<double>((*c1)[nn[inode]].x, (*c1)[nn[inode]].y, (*c1)[nn[inode]].z)
		                      : face_c0[inode];
	}
	FaceElementType *FaceEl = new FaceElementType(nodes);

	// quadrature rule
	int deg = (pbc->conwep && pbc->conwepswitch) ? VariableDegree : ConstantDegree;  // quadrature rule degree
	QuadratureRule c(deg);

	// local variables to be computed at each integration point
	double xi[2];                             // abscissa
	double weight;                            // weight
	double N[FaceElementType::NumberOfNodes]; // shape function values
	double dJNormal[3 * FaceElementType::NumberOfNodes][3]; // derivatives of normal
	double Xi[3], Normal[3];                  // reference coordinates and unit normal
	double pressure = pbc->val;

	Ktan.zero();

	// loop over the integration points
	for (int ip = 0; ip < c.getN(); ++ip) {

		// get the integration point abscissa and weight
		c.getAbscissaAndWeight(ip, xi, weight);

		// compute shape functions
		FaceEl->GetShapeFctVal(N, xi);

		// compute the partial derivatives of the surface normal w.r.t. the nodal coordinates
		FaceEl->GetdJNormal(dJNormal, xi, face_c1);

		if (pbc->conwep && pbc->conwepswitch) {
			FaceEl->GetUnitNormal(Normal, xi, face_c0);
			FaceEl->LocalToGlobalCoord(Xi, xi, face_c0);
			pressure = pbc->val + BlastLoading::ComputeGaussPointPressure(Xi, Normal, t, *(pbc->conwep));
		}

		for (int k = 0; k < FaceElementType::NumberOfNodes; ++k) {
			double tmp = N[k] * weight * pressure;
			for (int l = 0; l < 3 * FaceElementType::NumberOfNodes; ++l) {
				Ktan[3 * k + 0][l] += tmp * dJNormal[l][0];
				Ktan[3 * k + 1][l] += tmp * dJNormal[l][1];
				Ktan[3 * k + 2][l] += tmp * dJNormal[l][2];
			}
		}
	}

	delete FaceEl;
	for (int inode = 0; inode < FaceElementType::NumberOfNodes; ++inode) {
		delete face_c0[inode];
		if (c1) delete face_c1[inode];
	}
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
FullSquareMatrix
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::sommerMatrix(CoordSet &cs, double *d) const {
	FullSquareMatrix sommerM(numDofs(), d);
	sommerM.zero();

	return sommerM;
}

template<typename FaceElementType, typename QuadratureRule, int ConstantDegree, int VariableDegree, int eType>
int
PressureElement<FaceElementType, QuadratureRule, ConstantDegree, VariableDegree, eType>
::findAndSetEle(const CoordSet &cs, Elemset &eset, const Connectivity *nodeToElem, int *eleTouch, int *eleCount,
                int myNum,
                int it) {
	// overriding SommerElement::findAndSetEle because the normal should not be reversed
	this->iEle = findEle(nodeToElem, eleTouch, eleCount, myNum, &eset, it);
	if (iEle == -1) {
		std::cerr << "PressureElement::findAndSetEle could not find the corresponding element.\n";
		return 0;
	}
	el = eset[iEle];
	return -1;
}
