#ifndef KRYLOVPROJECTOR_H_
#define KRYLOVPROJECTOR_H_

template <class Scalar> class GenVector;
typedef GenVector<double> Vector;

template<class Scalar, class AnyVector>
class KrylovProjector 
{
   int maxVecStorage;  // maximum # of vectors to store
   int numKrylov;      // number of krylov spaces
   int numSystem;      // number of solved systems
   int lastIndex;
   int *kindex;
   int numDofs;

   Scalar *wKw;         // dot product
   Scalar *wMat;        // direction
   Scalar *kWMat;       // projection
   Scalar *y;

public:
   // Constructor
   KrylovProjector(int numDof, int maxDim);

   void project(const AnyVector *resid, AnyVector *precResid);
   int  initialization(AnyVector *rhs, AnyVector *x0);

   int  addDirection(const AnyVector &w, const AnyVector &Kw, Scalar dotProd);
   void newKrylov();
   void ortho(AnyVector *direction, int iter);
   void zeroKrylov();
   int  size() { return numKrylov; }
   void newSystem();

   // Functions added to access Krylov space
   double  getDotProduct(int i) { return wKw[i]; }
   double* getDirection(int i)  { return (wMat  + i*numDofs); }
   double* getProjection(int i) { return (kWMat + i*numDofs); }

};

#ifdef _TEMPLATE_FIX_
#include <Solvers.d/KProject.C>
#endif

#endif
