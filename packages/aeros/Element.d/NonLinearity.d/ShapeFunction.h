#ifndef _SHAPEFUNCTION_H_
#define _SHAPEFUNCTION_H_
#include <Utils.d/NodeSpaceArray.h>
#include <Math.d/Vector.h>

class Node;

class ShapeFunction
{
   protected:
     int numdofs;

   public:
     ShapeFunction(int n) : numdofs(n) {}
     virtual Tensor *getGradUInstance();
     virtual Tensor *getDgradUDqkInstance();
     virtual Tensor *getValInstance() = 0;
     virtual Tensor *getNodesCoordinatesInstance();
     virtual Tensor *getDisplacementsInstance();
     virtual Tensor *getLocalDerivativesInstance();
     virtual void getLocalDerivatives(Tensor *localDerivatives, double xi[3]) = 0;
     virtual void getValues(Tensor *val, double xi[3]) = 0;
     virtual void getGlobalGrads(Tensor *gradU, Tensor *dgradUdqk, double *jac, 
                                 Node *nodes, double xi[3], Vector &disp);
     virtual void getNodesCoordinates(Node *nodes, Tensor *nodescoordinates);
     virtual void getDisplacements(Vector &disp, Tensor *displacements);
     virtual void getGlobalGrads(Tensor *gradU, Tensor *dgradUdqk, double *jac,
                                 Tensor *nodescoordinates, double xi[3], Tensor *displacements,
                                 Tensor *localderivatives);
     virtual void getGradU(Tensor *gradU, Node *nodes, double xi[3], Vector &disp);
     virtual void getGradU(Tensor *gradU, double *jac, Node *nodes, double xi[3], Vector &disp);
     virtual void getJacobianDeterminant(double *jac, Node *nodes, double xi[3]);
     virtual double interpolateScalar(double *_q, double _xi[3]);
};

template <class TensorTypes>
class GenShapeFunction
{
   public:
     GenShapeFunction() {}
     virtual void getGlobalGrads(typename TensorTypes::GradUTensor &gradU,
                                 typename TensorTypes::GradUDerivTensor &dgradUdqk, double *jac,
                                 Node *nodes, double xi[3], Vector &disp) = 0;
     virtual void getGradU(typename TensorTypes::GradUTensor &gradU,
                           Node *nodes, double xi[3], Vector &disp) = 0;
     virtual double interpolateScalar(double *_q, double _xi[3]) = 0;
};

#endif
