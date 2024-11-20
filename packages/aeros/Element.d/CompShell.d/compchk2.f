C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compchk2( elm , ke )
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
      real*8    ke(18,18)
C
C.....Local Variables
C
      integer   i , j , nrot , poseig , negeig , nileig
      real*8    zero , eigU(18,18) , eigD(18) , ke2(18,18)
      real*8    eigL1(18) , eigL2(18) , eps , average , zeroeig
      logical   output , warning , shortwarning
C
C     ----
C     DATA
C     ----
C
      data zero         /0.000000D+00/
      data eps          /1.000000D-10/
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
      do 1001 j=1,18
         do 1002 i=1,18
            eigU(i,j) = zero
 1002    continue
 1001 continue
C
      do 1003 i=1,18
         eigD(i) = zero
 1003 continue
C
      do 1004 i=1,18
         eigL1(i) = zero
 1004 continue
C
      do 1005 i=1,18
         eigL2(i) = zero
 1005 continue
C
      do 1006 j=1,18
         do 1007 i=1,18
            ke2(i,j) = zero
 1007    continue
 1006 continue
C
      average = zero
      zeroeig = zero
C
      nrot    = 0
C
C.....INITIALIZE THE LOCAL STIFFNESS MATRIX
C.....(THE JACOBI ROUTINE OVERWRITES THE INPUT MATRIX)
C
      do 2001 j=1,18
         do 2002 i=1,18
            ke2(i,j) = ke(i,j)
 2002    continue
 2001 continue
C
C.....COMPUTE THE EIGENPAIRS WITH A FULL JACOBI SOLVER
C
      call compjac( ke2  , 18    , 18    , eigD ,
     $              eigU , eigL1 , eigL2 , nrot )
C
C.....ANALYZE THE RESULTS
C
      poseig = 0
      negeig = 0
      nileig = 0
C
      do 3001 i=1,18
         average = average + eigD(i)
 3001 continue
C
      average = average/18.00D+00
C
      zeroeig = eps*average
C
      do 3002 i=1,18
         if ( eigD(i).gt.zero ) then
            if ( eigD(i).gt.zeroeig ) then
               poseig = poseig + 1
            else
               nileig = nileig + 1
            endif
         endif
         if ( eigD(i).lt.zero ) then
            if ( abs(eigD(i)).gt.zeroeig ) then
               negeig = negeig + 1
            else
               nileig = nileig + 1
            endif
         endif
         if ( eigD(i).eq.zero ) then
            nileig = nileig + 1
         endif
 3002 continue
C
C.....OUTPUT TO USER
C
      if ( output ) then
C
      write(*,*)
      write(*,*) "+--------------------------------------------+"
      write(*,*) "| EIGENVALUES  OF  THE  ELEMENTAL  STIFFNESS |"
      write(*,*) "+--------------------------------------------+"
      write(*,1) elm
      write(*,*) "+--------------------------------------------+"
C
      do 4001 i=1,18
         write(*,2) i,eigD(i)
         write(*,*) "+--------------------------------------------+"
 4001 continue
C
      endif
C
C.....LONG WARNING MESSAGE
C
      if ( warning ) then
C
      if ( negeig.ne.0 ) then
C
      write(*,*)
      write(*,*) "+--------------------------------------------+"
      write(*,*) "| EIGENVALUES  OF  THE  ELEMENTAL  STIFFNESS |"
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
      if ( negeig.ne.0 ) then
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
    6 format(" WARNING: Stiffness Matrix of Composite Shell ",I5,
     $       " Has + / 0 / - Eigenvalues: ",I2," /",I2," /",I2)
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPCHK2"
C=========================C
