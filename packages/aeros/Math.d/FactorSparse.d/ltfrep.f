************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    09-23-98 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   LTFREP ... generate the full matrix                        *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine generates a full representation of a symmetric
*     matrix from a representation of the lower triangular part of
*     the matrix.  This uses more space than necessary, but it saves
*     time in some of the steps, particularly in inputting numerical
*     values into the data structure.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     COLPTR    (input/output) integer array, dimension N+1
*               Pointers to the column structure of the input matrix.
*
*     ROWIDX    (input/output) integer array, dimension COLPTR(N+1)-1
*               Row indices of nonzero entries in the input matrix,
*               stored by columns.
*
*     VALUES    (input/output) double precision array, dimension
*                   COLPTR(N+1)-1
*               Numerical values of nonzero entries in the input
*               matrix, stored by columns.
*
*     TMPPTR    (temporary) integer array, dimension N+1
*               Copy of COLPTR.
*
*     TMPIDX    (temporary) integer array, dimension COLPTR(N+1)-1
*               Copy of ROWIDX.
*
*     TMPVAL    (temporary) double precision array, dimension
*                   COLPTR(N+1)-1
*               Copy of VALUES.
*
*     TEMP      (temporary) integer array, dimension N
*               Nonzero counts for the columns.
*
************************************************************************
*
      SUBROUTINE  LTFREP  ( N     , COLPTR, ROWIDX, VALUES, TMPPTR,
     &                      TMPIDX, TMPVAL, TEMP                    )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             N
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             COLPTR(*), ROWIDX(*)
        DOUBLE PRECISION    VALUES(*)
        INTEGER             TMPPTR(*), TMPIDX(*)
        DOUBLE PRECISION    TMPVAL(*)
        INTEGER             TEMP(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             COLBEG, COLEND, I     , II    , IX    ,
     &                      J
        DOUBLE PRECISION    V
*
************************************************************************
*
*       ---------------
*       Initialization.
*       ---------------
        DO  J = 1, N
            TEMP(J) = 0
        END DO
        DO  J = 1, N+1
            TMPPTR(J) = COLPTR(J)
        END DO
        DO  J = 1, COLPTR(N+1) - 1
            TMPIDX(J) = ROWIDX(J)
            TMPVAL(J) = VALUES(J)
        END DO
*       -----------------------------------
*       Count number of nonzero entries
*       in each column of the input matrix.
*       -----------------------------------
        COLBEG = TMPPTR(1)
        DO  J = 1, N
            COLEND = TMPPTR(J+1)
            DO  II = COLBEG, COLEND-1
                I       = TMPIDX(II)
                TEMP(J) = TEMP(J) + 1
                IF  ( I .NE. J )  TEMP(I) = TEMP(I) + 1
            END DO
            COLBEG = COLEND
        END DO
*       ---------------------------
*       Create new column pointers.
*       ---------------------------
        COLPTR(1) = 1
        DO  J = 1, N
            COLPTR(J+1) = COLPTR(J) + TEMP(J)
            TEMP(J)     = COLPTR(J)
        END DO
*       -----------------------------------------------------------
*       Expand data structure to include the upper triangular part.
*       No attempt is made to sort the row indices in each column.
*       -----------------------------------------------------------
        COLBEG = TMPPTR(1)
        DO  J = 1, N
            COLEND = TMPPTR(J+1)
            DO  II = COLBEG, COLEND-1
*
                I = TMPIDX(II)
                V = TMPVAL(II)
*
                IX         = TEMP(J)
                ROWIDX(IX) = I
                VALUES(IX) = V
                TEMP(J)    = IX + 1
*
                IF  ( I .NE. J )  THEN
                    IX         = TEMP(I)
                    ROWIDX(IX) = J
                    VALUES(IX) = V
                    TEMP(I)    = IX + 1
                END IF
*
            END DO
            COLBEG = COLEND
        END DO
        RETURN
*
*       --------------
*       End of LTFREP.
*       --------------
      END
