************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    01-14-1999 (Michel Lesoinne and Kendall H. Pierson)
*
************************************************************************
************************************************************************
*****   ZZEROMAT ... input numerical values  (complex)             *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*   This subroutine zeroes the matrix (LNZ) and creates the INVSUPER
*   array.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*   
*     NSUPER    (input) integer
*               Number of supernodes.
*
*     XSUPER    (input) integer array, dimension NSUPER+1
*               The supernodal partition.
*
*     XLNZ      (input) integer array, dimension N+1
*               Pointers to nonzero entries in the Cholesky factor.
*   
*     LNZ       (output) COMPLEX*16 array, dimension XLNZ(N+1)-1
*               Numerical values of nonzero entries in the Cholesky
*               factor, stored by columns.
*
*     INVSUPER  (output) integer array dimension N
*               Keep track of dof to super node map
*               
************************************************************************
*
      SUBROUTINE  ZZEROMAT ( N, NSUPER, XSUPER, XLNZ, LNZ, INVSUPER )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             N     , NSUPER
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XSUPER(*)
        INTEGER             XLNZ(*)
        COMPLEX*16          LNZ(*)
        INTEGER             INVSUPER(*)
*
************************************************************************
*
*       --------------
*       Parameters ...
*       --------------
        COMPLEX*16         ZERO
        PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             II, JSUPER
*
************************************************************************
*
*       ------------------------------
*       Initialize the data structure.
*       ------------------------------
*        DO  II = 1, XLNZ(N+1)-1
*            LNZ(II) = ZERO
*        END DO
*
*
*       ----------------------
*       Fill INVSUPER array
*       ----------------------
        DO  JSUPER = 1, NSUPER
          DO II=XSUPER(JSUPER), XSUPER(JSUPER+1) - 1
              INVSUPER(II) = JSUPER
          END DO
        END DO
*
*       --------------
*       End of ZEROMAT.
*       --------------
      END
