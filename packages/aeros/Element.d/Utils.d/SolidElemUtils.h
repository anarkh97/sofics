#ifndef _SOLIDELEMUTILS_H_
#define _SOLIDELEMUTILS_H_

template <class Scalar> class GenFullSquareMatrix;
typedef GenFullSquareMatrix<double> FullSquareMatrix;

// For stiffness & mass matrices 3D solid elements 
void addBtCBtoK3DSolid(FullSquareMatrix &K, double (*DShape)[3], double C[6][6], double alpha, int nnodes, int* ls);
void addNtDNtoM3DSolid(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls, double (*D)[3] = 0);
void addDBtCDBtodKdx3DSolid(FullSquareMatrix *&dKdx, double (*DShape)[3], double DDShape[4][3][12], double C[6][6], double alpha, int nnodes, int* ls);

void addBtBtoK3DHelm(FullSquareMatrix &K, double (*DShape)[3], double alpha, int nnodes, int* ls);
void addNtNtoM3DHelm(FullSquareMatrix &M, double* Shape, double alpha, int nnodes, int* ls);

int checkJacobian(double *J, int *jSign, int elId, const char* mssg= 0, double atol = 0.0, bool stop=true, FILE* file=stderr);

// For anisotropic constitutive matrix and thermal expansion coefficients
void rotateConstitutiveMatrix(double *_Cin, double *T33, double Cout[6][6]);
void rotateVector(double *_w, double *cFrame, double *_alpha);

// For stresses & strains evaluation in case of ansitropic constitutive matrix
void computeEngStrain3DSolid(double Strain[6], double (*DShape)[3], double* U, int nnodes, int* ls=0);
void computeStress3DSolid(double Stress[6], double Strain[6], double C[6][6]);
void computeStressAndEngStrain3DSolid(double Stress[6], double Strain[6], double C[6][6], double (*DShape)[3], double* U, int nnodes, int* ls=0);
double computeVonMisesStress(double Stress[6]);
double computeVonMisesStrain(double Strain[6]);
 
#endif
