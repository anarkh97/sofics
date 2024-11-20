C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compchk( elm , dbb , dmm , dbm , dmb )
C=====================================================================C
C                                                                     C
C                                                                     C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      integer   elm
      real*8    dbb(3,3) , dmm(3,3) , dbm(3,3) , dmb(3,3)
C
C.....Local Variables
C
      integer   i , j , nrot , poseig , negeig , nileig
      real*8    zero , d(6,6) , eigU(6,6) , eigD(6)
      real*8    eigL1(6) , eigL2(6)
      logical   output , warning , shortwarning
C
C     ----
C     DATA
C     ----
C
      data zero         /0.000000D+00/
C
C.....INITIALIZE THE LOGICAL [OUTPUT] TO ".TRUE." FOR OUTPUT
C
      data output       /.false./
C
C.....INITIALIZE THE LOGICAL [WARNING] TO ".TRUE." FOR LONG WARNING
C
      data warning      /.false./
C
C.....INITIALIZE THE LOGICAL [SHORTWARNING] TO ".TRUE." FOR SHORT WARNING

C
      data shortwarning /.true./
C
C     -----
C     LOGIC
C     -----
C
C.....CLEAR LOCAL STORAGES
C
      do 1001 j=1,6
         do 1002 i=1,6
            d(i,j) = zero
 1002    continue
         do 1003 i=1,6
            eigU(i,j) = zero
 1003    continue
 1001 continue
C
      do 1004 i=1,6
         eigD(i) = zero
 1004 continue
C
      do 1005 i=1,6
         eigL1(i) = zero
 1005 continue
C
      do 1006 i=1,6
         eigL2(i) = zero
 1006 continue
C
      nrot = 0
C
C.....INITIALIZE THE 6 BY 6 CONSTITUTIVE MATRIX
C
      do 2001 j=1,3
         do 2002 i=1,3
            d(i,j) = dbb(i,j)
 2002    continue
 2001 continue
C
      do 2003 j=1,3
         do 2004 i=1,3
            d(i+3,j+3) = dmm(i,j)
 2004    continue
 2003 continue
C
      do 2005 j=1,3
         do 2006 i=1,3
            d(i,j+3) = dbm(i,j)
 2006    continue
 2005 continue
C
      do 2007 j=1,3
         do 2008 i=1,3
            d(i+3,j) = dmb(i,j)
 2008    continue
 2007 continue
C
C.....COMPUTE THE EIGENPAIRS WITH A FULL JACOBI SOLVER
C.....WARNING: THE LOCAL STORAGE [D] WILL BE OVERWRITTEN
C
      call compjac( d    , 6     , 6     , eigD ,
     $              eigU , eigL1 , eigL2 , nrot )
C
C.....ANALYZE THE RESULTS
C
      poseig = 0
      negeig = 0
      nileig = 0
C
      do 3001 i=1,6
         if ( eigD(i).gt.zero ) then
            poseig = poseig + 1
         endif
         if ( eigD(i).lt.zero ) then
            negeig = negeig + 1
         endif
         if ( eigD(i).eq.zero ) then
            nileig = nileig + 1
         endif
 3001 continue
C
C.....OUTPUT TO THE USER
C
      if ( output ) then
C
      write(*,*)
      write(*,*) "+--------------------------------------------+"
      write(*,*) "| EIGENVALUES  OF  THE  CONSTITUTIVE  MATRIX |"
      write(*,*) "+--------------------------------------------+"
      write(*,1) elm
      write(*,*) "+--------------------------------------------+"
C
      do 3002 i=1,6
         write(*,2) i,eigD(i)
         write(*,*) "+--------------------------------------------+"
 3002 continue
C
      endif
C
C.....LONG WARNING MESSAGE
C
      if ( warning ) then
C
      if ( poseig.ne.6 ) then
C
      write(*,*)
      write(*,*) "+--------------------------------------------+"
      write(*,*) "| EIGENVALUES  OF  THE  CONSTITUTIVE  MATRIX |"
      write(*,*) "+--------------------------------------------+"
      write(*,1) elm
      write(*,*) "+--------------------------------------------+"
      write(*,3) poseig
      write(*,*) "+--------------------------------------------+"
      write(*,4) nileig
      write(*,*) "+--------------------------------------------+"
      write(*,5) negeig
      write(*,*) "+--------------------------------------------+"
      write(*,*)
C
      endif
C
      endif
C
C.....SHORT WARNING MESSAGE
C
      if ( shortwarning ) then
C
      if ( poseig.ne.6 ) then
C
      write(*,6) elm,poseig,nileig,negeig
C
      endif
C
      endif
C
C     ------
C     RETURN
C     ------
C
      return
C
C     ------
C     FORMAT
C     ------
C
    1 format("|  COMPOSITE  SHELL  ELEMENT  NUMBER  ",I5,"  |")
    2 format("|   ",I5,"     Eigenvalue ",E16.9,"    |")
    3 format("|   Number of Positive Eigenvalues: ",I5,"    |")
    4 format("|   Number of Zero     Eigenvalues: ",I5,"    |")
    5 format("|   Number of Negative Eigenvalues: ",I5,"    |")
    6 format(" WARNING: Constitutive Matrix of Element ",I5,
     $       " Has + / 0 / - Eigenvalues: ",I2," /",I2," /",I2)
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPCHK"
C========================C
