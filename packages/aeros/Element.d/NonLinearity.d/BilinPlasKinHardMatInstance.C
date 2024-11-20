#include <Element.d/NonLinearity.d/BilinPlasKinHardMat.C>

#define ELASPLASKINHARDMAT_INSTANTIATION_HELPER(e) \
template \
ElasPlasKinHardMat<e>::ElasPlasKinHardMat(StructProp*);\
\
template \
void ElasPlasKinHardMat<e>::getStress(Tensor*, Tensor&, double*, double);\
\
template \
void ElasPlasKinHardMat<e>::getElasticity(Tensor*) const;\
\
template \
void ElasPlasKinHardMat<e>::getTangentMaterial(Tensor*, Tensor&, double*, double);\
\
template \
void ElasPlasKinHardMat<e>::getStressAndTangentMaterial(Tensor*, Tensor*, Tensor&, double*, double);\
\
template \
void ElasPlasKinHardMat<e>::updateStates(Tensor&, Tensor&, double*, double);\
\
template \
void ElasPlasKinHardMat<e>::integrate(Tensor*, Tensor*, Tensor&, Tensor&, double*, double*, double,\
                                      Tensor*, double) const;\
\
template \
void ElasPlasKinHardMat<e>::integrate(Tensor*, Tensor&, Tensor&, double*, double*, double,\
                                      Tensor*, double) const;\
\
template \
void ElasPlasKinHardMat<e>::initStates(double*);\
\
template \
StrainEvaluator* ElasPlasKinHardMat<e>::getStrainEvaluator() const;\
\
template \
bool ElasPlasKinHardMat<e>::getBackStress(double*, Tensor*);\
\
template \
bool ElasPlasKinHardMat<e>::getPlasticStrain(double*, Tensor*);\
\
template \
double ElasPlasKinHardMat<e>::getStrainEnergyDensity(Tensor&, double*, double);\
\
template \
double ElasPlasKinHardMat<e>::getDissipatedEnergy(double*);

ELASPLASKINHARDMAT_INSTANTIATION_HELPER(0);
ELASPLASKINHARDMAT_INSTANTIATION_HELPER(1);
ELASPLASKINHARDMAT_INSTANTIATION_HELPER(2);

