************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    09-18-87 (Barry W. Peyton)
*   Modified:   09-23-98 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   GETADJ ... generate the full adjacency structure           *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine extracts the adjacency structure from the full
*     representation of the input matrix, which is assumed to be
*     stored by columns.
*
*   ----------
*   Arguments:
*   ----------
*
*     N         (input) integer
*               Number of equations.
*
*     COLPTR    (input) integer array, dimension N+1
*               Pointers to the column structure of the input matrix.
*
*     ROWIDX    (input) integer array, dimension COLPTR(N+1)-1
*               Row indices of nonzero entries in the input matrix,
*               stored by columns.
*
*     XADJ      (output) integer array, dimension N+1
*               Pointers to the adjacency structure.
*
*     ADJNCY    (output) integer array, dimension XADJ(N+1)-1
*               The adjacency structure.
*
************************************************************************
*
      SUBROUTINE  GETADJ  ( N     , COLPTR, ROWIDX, XADJ  , ADJNCY  )
*
************************************************************************
*
*       --------------------
*       Scaler Arguments ...
*       --------------------
        INTEGER             N
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             COLPTR(*), ROWIDX(*)
        INTEGER             XADJ(*)  , ADJNCY(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             COLBEG, COLEND, CURENT, J     , K     ,
     &                      ROWIND
*
************************************************************************
*
        COLBEG  = COLPTR(1)
        CURENT  = 1
        XADJ(1) = 1
        DO  J = 1, N
            COLEND = COLPTR(J+1)
            DO  K = COLBEG, COLEND-1
                ROWIND = ROWIDX(K)
                IF  ( ROWIND .NE. J )  THEN
                    ADJNCY(CURENT) = ROWIND
                    CURENT         = CURENT + 1
                END IF
            END DO
            XADJ(J+1) = CURENT
            COLBEG    = COLEND
        END DO
        RETURN
*
*       --------------
*       End of GETADJ.
*       --------------
      END
