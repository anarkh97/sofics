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
*****   FNTSIZ ... compute work storage size for BLKLDL            *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine determines the size of the work space required
*     by BLKLDL.
*
*   ----------
*   Arguments:
*   ----------
*
*     NSUPER    (input) integer
*               Number of supernodes.
*
*     XSUPER    (input) integer array, dimension NSUPER+1
*               The supernodal partition.
*
*     SNODE     (input) integer array, dimension N
*               Membership of columns in the supernode partition.
*
*     XLINDX    (input) integer array, dimension NSUPER+1
*               Pointers to column structure of the Cholesky factor.
*
*     LINDX     (input) integer array, dimension XSUPER(NSUPER+1)-1
*               Row indices of nonzero entries in the Cholesky factor,
*               stored by columns, in a compressed representation.
*
*     TMPSIZ    (output) integer
*               Size of work space in TEMP required by BLKLDL.
*
*     RWSIZE    (output) integer
*               Size of work space in RWORK required by BLKLDL.
*
************************************************************************
*
      SUBROUTINE  FNTSIZ  ( NSUPER, XSUPER, SNODE , XLINDX, LINDX ,
     &                      TMPSIZ, RWSIZE                          )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             NSUPER, RWSIZE, TMPSIZ
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XSUPER(*)
        INTEGER             SNODE (*)
        INTEGER             XLINDX(*), LINDX(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             BOUND , CLEN  , CURSUP, I     , IBEGIN,
     &                      IEND  , KSUP  , LENGTH, NCOLS , NXTSUP, 
     &                      TSIZE , WIDTH
*
*       ----------------------
*       Intrinsic Functions ...
*       ----------------------
        INTRINSIC           MAX
*
************************************************************************
*
*       ----------------------------------------------------------
*       Returns size of temporary array used by BLKLDL subroutine.
*       Note that the value returned is an estimate, though it is
*       usually tight.
*       ----------------------------------------------------------
*
*       ----------------------------------------------------------
*       Compute size of temporary storage vector needed by BLKLDL.
*       ----------------------------------------------------------
        TMPSIZ = 0
        RWSIZE = 0
        DO  KSUP = NSUPER, 1, -1
            NCOLS  = XSUPER(KSUP+1) - XSUPER(KSUP)
            IBEGIN = XLINDX(KSUP) + NCOLS
            IEND   = XLINDX(KSUP+1) - 1
            LENGTH = IEND - IBEGIN + 1
            BOUND  = LENGTH*LENGTH
            RWSIZE = MAX(RWSIZE,NCOLS*(LENGTH+NCOLS))
            IF  ( LENGTH .GT. 0 .AND. BOUND .GT. TMPSIZ )  THEN
                CURSUP = SNODE(LINDX(IBEGIN))
                CLEN   = XLINDX(CURSUP+1) - XLINDX(CURSUP)
                WIDTH  = 0
                DO  I = IBEGIN, IEND
                    NXTSUP = SNODE(LINDX(I))
                    IF  ( NXTSUP .EQ. CURSUP )  THEN
                        WIDTH = WIDTH + 1
                        IF  ( I .EQ. IEND )  THEN
                            IF  ( CLEN .GT. LENGTH )  THEN
                                TSIZE  = LENGTH*WIDTH
                                TMPSIZ = MAX(TSIZE,TMPSIZ)
                            END IF
                        END IF
                    ELSE
                        IF  ( CLEN .GT. LENGTH )  THEN
                            TSIZE  = LENGTH*WIDTH
                            TMPSIZ = MAX(TSIZE,TMPSIZ)
                        END IF
                        LENGTH = LENGTH - WIDTH
                        BOUND  = LENGTH*LENGTH
                        IF  ( BOUND .LE. TMPSIZ )  GO TO 100
                        WIDTH  = 1
                        CURSUP = NXTSUP
                        CLEN   = XLINDX(CURSUP+1) - XLINDX(CURSUP)
                    END IF
                END DO
            END IF
  100       CONTINUE
        END DO
        RETURN
*
*       --------------
*       End of FNTSIZ.
*       --------------
      END
