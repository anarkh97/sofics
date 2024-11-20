C
C-------------------------------------------------------------------C
C                                                                   C
C AUTHOR: C. FARHAT                                                 C
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
      subroutine micspsmvp(ncolu, unonz, xunonz, rowu, v, RES)

      integer xunonz(*), rowu(*)
      integer ncolu
      real*8  unonz(*), v(*), res(*)
C
C.... LOCAL DECLARATIONS
C
      integer jcol, irow, i, istart,  istop
C
C.... FIRST PASS
C
      call micspmvp(unonz, xunonz, rowu, ncolu, ncolu, v, res)
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
C pjsa 9-7-05: following line is invalid if istop <= 0
C     if(rowu(istop).eq.jcol) istop = istop - 1
      if(istop.gt.0 .and. rowu(istop).eq.jcol) istop = istop - 1
      do 200 i  = istart,istop
        irow      = rowu(i)
        res(jcol) = res(jcol) + unonz(i)*v(irow)
200   continue
100   continue

c      ncolu = 2*ncolu
      return
      end
C
      subroutine micspmvp(nonz, xnonz, row, nrow, ncol, v, res)
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
      integer xnonz(*),row(*),nrow,ncol
      real*8  nonz(*),v(*),res(*)
C
C.... LOCAL DECLARATIONS
C
      integer jcol,irow,i,istart,istop
      real*8  vm
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
        res(irow)= res(irow) + vm*nonz(i)
200   continue
100   continue
C
      return
      end
C
      subroutine micfincol(jcol, xnonz, ncol, istart, istop)
C-------------------------------------------------------------------C
C                                                                   C
C AUTHOR: C. FARHAT                                                 C
C                                                                   C
C     LOCATES COLUMN JCOL IN THE SPARSE DATA STRUCTURE DEFINED      C
C     PARTIALLY BY XNONZ. EMPTY COLUMN IS ASSUMED TO BE FLAGGED     C
C     AS XNONZ = 0.                                                 C 
C                                                                   C
C     INPUT :                                                       C
C             JCOL, XNONZ, NCOL                                     C
C                                                                   C
C     OUTPUT:                                                       C
C             ISTART, ISTOP                                         C
C                                                                   C
C-------------------------------------------------------------------C
C
      integer xnonz(1)
      integer jcol,ncol,istart,istop
C
C-----LOCAL DECLARATIONS
C
      integer jc1,j,istop1
C       
      istart = xnonz(jcol)
      jc1    = jcol + 1
C
      do 100 j = jc1,ncol
        istop1   = xnonz(j)
        if(istop1.ne.0) go to 200
100   continue
      istop1 = xnonz(ncol + 1)
200   istop  = istop1 - 1
C
      return
      end

