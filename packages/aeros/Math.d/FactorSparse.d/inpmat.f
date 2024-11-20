************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    Esmond G. Ng and Barry W. Peyton
*   Created:    09-17-1987 (Barry W. Peyton)
*   Modified:   09-23-1998 (Esmond G. Ng)
*
*   Computational Mathematics and Statistics Section
*   Oak Ridge National Laboratory
*
************************************************************************
************************************************************************
*****   INPMAT ... read in Harwell-Boeing matrix                   *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*     This subroutine reads in the structure and numerical values of
*     a real symmetric matrix stored in the Harwell-Boeing format.  It
*     assumes that the header information has been read in in the main
*     program.  There is the option to skip the numerical values in
*     the input file, in which case the subroutine will generate the
*     numerical values so that the matrix is diagonally dominant.
*
*   ----------
*   Arguments:
*   ----------
*
*     MATUNT    (input) integer
*               Unit from which the matrix is read.
*
*     VALKEY    (input) integer
*               A flag indicating if numerical values from the input
*               file should be read.
*
*               VALKEY = 0 : read numerical values from input file.
*               VALKEY = 1 : generate numerical values.
*
*     PTRFMT    (input) character string
*               Format for reading the column pointers.
*
*     INDFMT    (input) character string
*               Format for reading the row indices.
*
*     VALFMT    (input) character string
*               Format for reading the numerical values.
*
*     N         (output) integer
*               Number of equations.
*
*     NNZERO    (output) integer
*               Number of nonzero entries.
*
*     COLPTR    (output) integer array, dimension N+1
*               Pointers to the column structure of the input matrix.
*
*     ROWIDX    (output) integer array, dimension NNZERO
*               Row indices of nonzero entries in the input matrix,
*               stored by columns.
*
*     VALUES    (output) double precision array, dimension NNZERO
*               Numerical values of nonzero entries in the input
*               matrix, stored by columns.
*
*     TEMP      (temporary) integer array, dimension N
*               TEMP(I) is the number of off-diagonal nonzero entries
*               in column/row I of the input matrix.
*
*     IFLAG     (output) integer
*               Error flag.
*
*                   IFLAG = 0:  no error.
*                   IFLAG = 1:  end of file while reading pointers.
*                   IFLAG = 2:  end of file while reading row indices.
*                   IFLAG = 3:  end of file while reading numerical
*                               values.
*
************************************************************************
*
      SUBROUTINE  INPMAT  ( MATUNT, VALKEY, PTRFMT, INDFMT, VALFMT,
     &                      N     , NNZERO, COLPTR, ROWIDX, VALUES,
     &                      TEMP  , IFLAG                           )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        CHARACTER           PTRFMT*16, INDFMT*16, VALFMT*16
        INTEGER             IFLAG , MATUNT, VALKEY
        INTEGER             N     , NNZERO
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             COLPTR(*), ROWIDX(*)
        DOUBLE PRECISION    VALUES(*)
        INTEGER             TEMP(*)
*
************************************************************************
*
*       --------------
*       Parameters ...
*       --------------
        DOUBLE PRECISION    ONE   , ZERO
        PARAMETER       (   ONE  = 1.0D0,
     &                      ZERO = 0.0D0 )
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             COLBEG, COLEND, I     , J     , JJ    ,
     &                      K
*
************************************************************************
*
        IFLAG = 0
*
*       ------------------------------
*       Read pointers and row indices.
*       ------------------------------
        READ (MATUNT,PTRFMT,END=100)  (COLPTR(K),K=1,N+1)
        READ (MATUNT,INDFMT,END=200)  (ROWIDX(K),K=1,NNZERO)
        IF  ( VALKEY .EQ. 0 )  THEN
*           ----------------------
*           Read numerical values.
*           ----------------------
            READ (MATUNT,VALFMT,END=300)  (VALUES(K),K=1,NNZERO)
        ELSE
*           --------------------------
*           Generate numerical values.
*           --------------------------
            DO  K = 1, N
                TEMP(K) = ZERO
            END DO
*           -----------------------------------------
*           First count the number of nonzero entries
*           in each row/column.
*           -----------------------------------------
            COLBEG = COLPTR(1)
            DO  J = 1, N
                COLEND = COLPTR(J+1)
                DO  JJ = COLBEG, COLEND-1
                    I = ROWIDX(JJ)
                    IF  ( I .NE. J )  THEN
                        TEMP(I) = TEMP(I) + ONE
                        TEMP(J) = TEMP(J) + ONE
                    END IF
                END DO
                COLBEG = COLEND
            END DO
*           ------------------------------------------
*           Assign values so that the resulting matrix
*           is diagonally dominant.
*           ------------------------------------------
            COLBEG = COLPTR(1)
            DO  J = 1, N
                COLEND = COLPTR(J+1)
                DO  JJ = COLBEG, COLEND-1
                    I = ROWIDX(JJ)
                    IF  ( I .EQ. J )  THEN
                        VALUES(JJ) = 2*(TEMP(J) + 1)
                    ELSE
                        VALUES(JJ) = - ONE
                    END IF
                END DO
                COLBEG = COLEND
            END DO
        END IF
*       --------------
*       Normal return.
*       --------------
        RETURN
*
*       ----------------
*       Abnormal return.
*       ----------------
  100   CONTINUE
        IFLAG = 1
        RETURN
*
  200   CONTINUE
        IFLAG = 2
        RETURN
*
  300   CONTINUE
        IFLAG = 3
        RETURN
*
*       --------------
*       End of INPMAT.
*       --------------
      END
