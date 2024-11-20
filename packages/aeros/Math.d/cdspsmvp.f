C
C-------------------------------------------------------------------C
C                                                                   C
C     PERFORMS A COLUMN ORIENTED SPARSE SYMMETRIC_MATRIX-VECTOR     C
C     PRODUCT                                                       C
C                                                                   C
C     INPUT  :                                                      C
C             UNONZ, XUNONZ, ROWU, V                                C
C                                                                   C
C     OUTPUT :                                                      C
C             RES                                                   C
C                                                                   C
C-------------------------------------------------------------------C
C
      subroutine cdspsmvp(ncolu, unonz, xunonz, rowu, v, RES)

      integer xunonz(1), rowu(1)
      integer ncolu
      real*8 unonz(1)
      complex*16  v(1), res(1)
C
C.... LOCAL DECLARATIONS
C
      integer jcol, irow, i, istart,  istop
C
C.... FIRST PASS
C
      call cdspmvp(unonz, xunonz, rowu, ncolu, ncolu, v, res)
C
C.... SECOND PASS
C
      do 100 jcol = 2, ncolu
C
C.... SKIP EMPTY COLUMN
C
      if(xunonz(jcol).eq.0) go to 100
C
C.... LOCATE COLUMN JCOL
C
      call micfincol(jcol,xunonz,ncolu,istart,istop)
C
C.... PERFORM INNER PRODUCT
C
            if(rowu(istop).eq.jcol) istop = istop - 1
            do 200 i  = istart,istop
            irow      = rowu(i)
            res(jcol) = res(jcol) + (unonz(i)*v(irow))
C            res(jcol) = res(jcol) +
C     *        dcmplx(unonz(i)*dreal(v(irow)),unonz(i)*dimag(v(irow)))
200         continue
100   continue
C

      return
      end
C
      subroutine cdspmvp(nonz, xnonz, row, nrow, ncol, v, res)
C-------------------------------------------------------------------C
C                                                                   C 
C     PERFORMS A COLUMN ORIENTED SPARSE MATRIX-VECTOR PRODUCT       C
C                                                                   C
C     INPUT  :                                                      C
C             NONZ,XNONZ,ROW,V                                      C
C                                                                   C
C     OUTPUT :                                                      C
C             RES                                                   C
C                                                                   C
C-------------------------------------------------------------------C
C
      integer xnonz(1),row(1),nrow,ncol
      real*8 nonz(1)
      complex*16  v(1),res(1)
C
C.... LOCAL DECLARATIONS
C
      integer jcol,irow,i,istart,istop
      complex*16  vm
C
C.... INITIALIZE RES
C
      do 10 irow = 1, nrow
      res(irow)  = 0.0d0
10    continue
C
C.... PERFORM A COLUMN ORIENTED PRODUCT
C
      do 100 jcol = 1, ncol
      vm          = v(jcol)
C
C.... EXPLOIT SPARSITY IN V
C
      if(vm.eq.0.d0) go to 100
C
C.... SKIP EMPTY COLUMN 
C
      if(xnonz(jcol).eq.0) go to 100
C
C.... LOCATE COLUMN JCOL
C
      call micfincol(jcol,xnonz,ncol,istart,istop)
C
C.... OUTER PRODUCT
C
             do 200 i = istart,istop
             irow     = row(i)
             res(irow)= res(irow) + (nonz(i)*vm)
C             res(irow) = res(irow)+
C     * dcmplx(nonz(i)*dreal(vm),nonz(i)*dimag(vm))

200          continue
100   continue
C

      return
      end
