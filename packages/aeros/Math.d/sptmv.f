      subroutine sptmv(nonz,xnonz,row,ncol,v,res) 
c---------------------------------------------------------------
c
c     Performs a sparse (transpose matrix)-vector product
c
c     input  :
c             nonz,xnonz,row,v
c
c     output :
c             res
c
c---------------------------------------------------------------
c
      integer xnonz(*),row(*),ncol
      real*8  nonz(*),v(*),res(*)
c
c-----LOCAL DECLARATIONS
c
      integer jcol,irow,i,istart,istop
      real*8  s
c
c
c
      do 10 i    = 1, ncol
      res(i)     = 0.0d0
10    continue
c
      do 100 jcol = 1, ncol
c
c     skip empty column 
c
      if(xnonz(jcol).eq.0) go to 100
c
c     locate column jcol
c
c      call locate(jcol,xnonz,ncol,istart,istop)
      call micfincol(jcol,xnonz,ncol,istart,istop)
c
c     product
c
      s         = 0.0d0
c
             do 200 i = istart,istop
             irow     = row(i)
             s        = s + nonz(i)*v(irow)
200          continue
c
      res(jcol) = s
100   continue                     
c
c
c
      return
      end

