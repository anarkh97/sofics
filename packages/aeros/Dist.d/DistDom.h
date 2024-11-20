#ifndef _DIST_DOM_
#define _DIST_DOM_

#include <Driver.d/Domain.h>
#include <Driver.d/SubDomain.h>
#include <Utils.d/MyComplex.h>
#include <Driver.d/DecDomain.h>

template<class Scalar>
class GenDistrDomain : virtual public GenDecDomain<Scalar> 
{
  private:
    DistrInfo masterInfo;         // used to create master dist vecs
    int **masterFlag;             // master flag for distributed printout
    int *numFlags;
    int *nodeOffsets;             // offsets for distributed writing
    int *elemNodeOffsets;
    int *elemOffsets;
    int *numRes;		  // number of solutions
    FSCommPattern<Scalar> *nodePat;
    int x;
    DistVec<Scalar> *masterStress;
    Connectivity *clusToCpu;
  public:
    GenDistrDomain(Domain *d);
    virtual ~GenDistrDomain();

    void clean();

    void postProcessing(GenDistrVector<Scalar> &u, GenDistrVector<Scalar> &f, double eigV = 0.0,
                        GenDistrVector<Scalar> *aeroF = 0, int x = 0, GenMDDynamMat<Scalar> *dynOps = 0,
                        SysState<GenDistrVector<Scalar> > *distState = 0, int ndflag = 0); 
    void postProcessing(DistrGeomState *u, GenDistrVector<Scalar> &, Corotator ***, double x = 0,
                        SysState<GenDistrVector<Scalar> > *distState = 0, GenDistrVector<Scalar> *aeroF = 0,
                        DistrGeomState *refState = 0, GenDistrVector<Scalar> *reactions = 0,
                        GenMDDynamMat<Scalar> *dynOps = 0, GenDistrVector<Scalar> *resF = 0);
    virtual void forceContinuity(GenDistrVector<Scalar> &u);

    void setsizeSfemStress(int fileNumber);
    int getsizeSfemStress() { return this->sizeSfemStress; }
    Scalar * getSfemStress(int fileNumber);
    void updateSfemStress(Scalar* str, int fileNumber);

  private:
    void initialize();
    void initPostPro();
    void makeNodePat();
    void makeMasterInfo();
    void createMasterFlag();
    void createOutputOffsets();
    template<int dim>
    void getPrimal(DistSVec<Scalar, dim> &disps, DistSVec<Scalar, dim> &masterDisps,
                   double time, int x, int fileNumber, int ndof, int startdof);//DofSet::max_known_nonL_dof
    void getAeroForceScalar(DistSVec<Scalar, 6> &aerof, DistSVec<Scalar, 6> &masterAeroF,
                            double time, int x, int fileNumber, int dof);
    void getStressStrain(GenDistrVector<Scalar> &u, double time, int x, int fileNumber, int Findex, int printFlag=0);
    void getStressStrain(GenDistrVector<Scalar> &sol, int fileNumber, int stressIndex, double time, int printFlag=0) 
            { getStressStrain(sol, time, 0, fileNumber, stressIndex, printFlag); }
    void getElementStressStrain(GenDistrVector<Scalar> &, double, int, int, int, int printFlag=0);
    void getPrincipalStress(GenDistrVector<Scalar> &, double, int, int, int);
    void getElementPrincipalStress(GenDistrVector<Scalar> &, double, int, int, int);
    void getElementForce(GenDistrVector<Scalar> &, double, int, int, int);
    void getElementAttr(int, int, double);
    void getStressStrain(DistrGeomState *gs, Corotator ***allCorot, double time,
                         int x, int fileNumber, int Findex, DistrGeomState *refState);
    void getElementStressStrain(DistrGeomState *gs, Corotator ***allCorot, double time,
                                int iter, int fileNumber, int Findex, DistrGeomState *refState);
    void getPrincipalStress(DistrGeomState *gs, Corotator ***allCorot, double time,
                            int x, int fileNumber, int strIndex, DistrGeomState *refState);
    void getElementPrincipalStress(DistrGeomState *gs, Corotator ***allCorot, double time,
                                   int x, int fileNumber, int strIndex, DistrGeomState *refState);
    void getElementForce(DistrGeomState *gs, Corotator ***allCorot, double time, int x,
                         int fileNumber, int Findex);
    void unify(DistSVec<Scalar, 11> &vec);
    void getDeletedElements(int iOut);
};

#ifdef _TEMPLATE_FIX_
  #include <Dist.d/DistDom.C>
#endif

#endif
