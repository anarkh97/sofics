#ifndef SCBLACS_H_
#define SCBLACS_H_


#ifdef SCARRAYS_DEV
#include "linkfc.h"
#else
#include <Utils.d/linkfc.h>
#endif


#ifdef __cplusplus
extern "C" {
    void _FORTRAN(blacs_pinfo)(int *mypnum, int *nprocs);
    void _FORTRAN(blacs_get)(int *context, int *what, int *val);
    void _FORTRAN(blacs_gridinit)(int *context, char *order, int *nprow, int *npcol);
    void _FORTRAN(blacs_gridinfo)(int *context, int *nprow, int *npcol, int *myrow, int *mycol);
    void _FORTRAN(blacs_gridexit)(int *context);
    void _FORTRAN(blacs_barrier)(int *context, char *scope);
    void _FORTRAN(blacs_exit)(int *);
    void _FORTRAN(igsum2d)(int *context, char *scope, char *top, int *m, int *n, int *A, int *lda, int *rdest, int *cdest);
    void _FORTRAN(dgsum2d)(int *context, char *scope, char *top, int *m, int *n, double *A, int *lda, int *rdest, int *cdest);
    void _FORTRAN(igamx2d)(int *context, char *scope, char *top, int *m, int *n, int *A, int *lda, int *ra, int *ca,
                           int *rcflag, int *rdest, int *cdest);
    void _FORTRAN(igamn2d)(int *context, char *scope, char *top, int *m, int *n, int *A, int *lda, int *ra, int *ca,
                           int *rcflag, int *rdest, int *cdest);
    void _FORTRAN(dgebs2d)(int *context, char *scope, char *top, int *m, int *n, double *A, int *lda);
    void _FORTRAN(dgebr2d)(int *context, char *scope, char *top, int *m, int *n, double *A, int *lda, int *rsrc, int *csrc);
    void _FORTRAN(igebs2d)(int *context, char *scope, char *top, int *m, int *n, int *A, int *lda);
    void _FORTRAN(igebr2d)(int *context, char *scope, char *top, int *m, int *n, int *A, int *lda, int *rsrc, int *csrc);
    int _FORTRAN(blacs_pnum)(int *context, int *iprow, int *jpcol);
}
#endif

#endif /* SCBLACS_H_ */
