************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    12-27-1994 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   LSTATS ... gather statistics about Cholesky factorization  *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine prints statistics on sparse Cholesky
*     factorization.
*
*   ----------
*   Arguments:
*   ----------
*
*     NSUPER    (input) integer
*               Number of supernodes.  (N = XSUPER(NSUPER+1)-1.)
*
*     XSUPER    (input) integer array, dimension NSUPER+1
*               The supernodal partition.
*
*     XLINDX    (input) integer array, dimension NSUPER+1
*               Pointers to column structure of the Cholesky factor.
*
*     XLNZ      (input) integer array, dimension N+1
*               Pointers to nonzero entries in the Cholesky factor.
*
*     TMPSIZ    (input) integer
*               Size of work space in TEMP required by BLKLDL.
*
*     RWSIZE    (input) integer
*               Size of work space in RWORK required by BLKLDL.
*
*     OUTUNT    (input) integer
*               Output unit.
*
************************************************************************
*
      SUBROUTINE  LSTATS  ( NSUPER, XSUPER, XLINDX, XLNZ  , TMPSIZ,
     &                      RWSIZE, OUTUNT                          )
*       
************************************************************************
*                           
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             NSUPER, OUTUNT, RWSIZE, TMPSIZ
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XSUPER(*), XLINDX(*), XLNZ(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             J     , JLEN  , JSIZE , JSUPER, MAXSUP,
     &                      N     , NCOLS , NOFNZ , NOFSUB, SUPSZE
        DOUBLE PRECISION    FCTOPS, SLVOPS
*
*       --------------
*       Parameters ...
*       --------------
        REAL                ZERO
        PARAMETER       (   ZERO = 0.0D0 )
*
************************************************************************
*
        N = XSUPER(NSUPER+1) - 1
*
        WRITE (OUTUNT,*)  ' '
*       --------------------------------------------------------
*       Determine the number of nonzero entries in the Cholesky
*       factor and the number of row indices in representing the
*       supernodal structure.
*       --------------------------------------------------------
        NOFNZ  = XLNZ(N+1) - 1
        NOFSUB = XLINDX(NSUPER+1) - 1
        WRITE (OUTUNT,1) 
     &      'NUMBER OF SUPERNODES                  = ', NSUPER
        WRITE (OUTUNT,1) 
     &      'NUMBER OF NONZEROS IN L               = ', NOFNZ
        WRITE (OUTUNT,1) 
     &      'NUMBER OF SUBSCRIPTS IN L             = ', NOFSUB
*
*       -------------------------------------------------------
*       Determine the largest supernode in the Cholesky factor.
*       -------------------------------------------------------
        MAXSUP = 0
        SUPSZE = 0
        DO  JSUPER = 1, NSUPER
*           ---------------------------------------------------
*           NCOLS is the number of columns in supernode JSUPER.
*           ---------------------------------------------------
            NCOLS = XSUPER(JSUPER+1) - XSUPER(JSUPER)
            IF  ( NCOLS .GT. MAXSUP )  MAXSUP = NCOLS
*
*           ---------------------------------------------------
*           JSIZE is the number of nonzero entries in supernode
*           JSUPER.
*           ---------------------------------------------------
            JLEN  = XLINDX(JSUPER+1) - XLINDX(JSUPER)
            JSIZE = ((2*JLEN - NCOLS + 1)*NCOLS)/2
            IF  ( JSIZE .GT. SUPSZE )  SUPSZE = JSIZE
        END DO
        WRITE (OUTUNT,1) 
     &      'LARGEST SUPERNODE BY COLUMNS          = ', MAXSUP
        WRITE (OUTUNT,1) 
     &      'LARGEST SUPERNODE BY NONZEROS         = ', SUPSZE
*
        WRITE (OUTUNT,1) 
     &      'SIZE OF WORK SPACE IN TEMP            = ', TMPSIZ
        WRITE (OUTUNT,1) 
     &      'SIZE OF WORK SPACE IN RWORK           = ', RWSIZE
*
*       ---------------------------
*       Determine operation counts.
*       ---------------------------
        SLVOPS = ZERO
        FCTOPS = ZERO
        DO  J = 1, N
            JLEN   = XLNZ(J+1) - XLNZ(J)
            SLVOPS = SLVOPS + 2*JLEN - 1
            FCTOPS = FCTOPS + JLEN**2 - 1
        END DO
        SLVOPS = 2*SLVOPS
        WRITE (OUTUNT,2) 
     &      'FACTORIZATION OPERATION COUNT         = ', FCTOPS
        WRITE (OUTUNT,2) 
     &      'TRIANGULAR SOLN OPERATION COUNT       = ', SLVOPS
*
    1   FORMAT ( A40, I10 )
    2   FORMAT ( A40, 1PD20.10 )
*
        RETURN
*
*       --------------
*       End of LSTATS.
*       --------------
      END
