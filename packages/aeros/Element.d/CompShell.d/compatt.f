C=====================================================================C
C        1         2         3         4         5         6         7C
C234567890123456789012345678901234567890123456789012345678901234567890C
C=====================================================================C
      subroutine compatt( ncmpat    , ncmpfr    , nen1     , numel    ,
     $                    tycmp     , tycfr     , iatrb               )
C=====================================================================C
C THIS SUBROUTINE WILL SET UP POINTERS IN THE ATTRIBUTE ARRAY FOR     C
C QUICK ACCESS OF THE ATTRIBUTE AND FRAME NUMBERS OF THE COMPOSITE    C
C SHELL ELEMENTS.                                                     C
C=====================================================================C
C     AUTHOR :   Francois M. Hemez                                    C
C     DATE   :   June 7, 1995                                         C
C     VERSION:   1.0                                                  C
C=====================================================================C
C                                                                     C
C     -----------------                                               C
C     V A R I A B L E S                                               C
C     -----------------                                               C
C                                                                     C
C     NCMPAT = Number of Input Attributes (Read in the Input File)    C
C                                                                     C
C     NCMPFR = Number of Input Frames (Read in the Input File)        C
C                                                                     C
C     NEN1   = Number of Rows of the Attribute Array (NEN1 = 8)       C
C                                                                     C
C     NUMEL  = Number of Finite Elements in the Mesh                  C
C                                                                     C
C     TYCMP  = Pointers for Storage of Composite Attributes           C
C                                                                     C
C     TYCFR  = Pointers for Storage of Composite Frames               C
C                                                                     C
C     IATRB  = Attribute Array                                        C
C                                                                     C
C=====================================================================C
C
C     -----------
C     DECLARATION
C     -----------
C
C.....Global Variables
C
      integer   ncmpat , ncmpfr , nen1 , numel
      integer   tycmp(3,ncmpat) , tycfr(2,ncmpfr)
      integer   iatrb(nen1,numel)
C
C.....Local Variables
C
      integer   ife , iatt , ilaw , iframe
      integer   lawnumber , typelaw , framenumber , address
C
C     -----
C     LOGIC
C     -----
C
C.....SANITY CHECK ON THE SIZE OF ARRAY [IATRB]
C
      if ( nen1.lt.8 ) go to 100
C
C.....CLEAR THE STORAGE IN THE ATTRIBUTE ARRAY
C
      do 1000 ife=1,numel
         do 1001 iatt=6,8
            iatrb(iatt,ife) = 0
 1001    continue
 1000 continue
C
C.....SET UP POINTERS FOR THE COMPOSITE CONSTITUTIVE LAW
C
      do 2000 ilaw=1,ncmpat
C
C.....GET THE INPUT CONSTITUTIVE LAW NUMBER
C
      lawnumber = tycmp(1,ilaw)
C
C.....GET THE TYPE OF CONSTITUTIVE LAW
C
      typelaw = tycmp(2,ilaw)
C
C.....GET THE ADDRESS IN STORAGE (FIRST COLUMN IN COMPOSITE STORAGE)
C
      address = tycmp(3,ilaw)
C
C.....LOOP ON FINITE ELEMENTS
C
      do 2100 ife=1,numel
C
C.....SET UP POINTERS IF THE ELEMENT IS RELEVANT
C
      if ( iatrb(4,ife).eq.lawnumber ) then
C
C.....STORE THE TYPE OF CONSTITUTIVE LAW
C
      iatrb(6,ife) = typelaw
C
C.....TREATMENT FOR A TYPE "-1" CONSTITUTIVE LAW:
C.....ERROR MESSAGE BECAUSE THE ELEMENT IS NOT A COMPOSITE SHELL
C
      if ( typelaw.eq.-1 ) then
         go to 200
      endif
C
C.....TREATMENT FOR A TYPE "0" CONSTITUTIVE LAW:
C.....DO NOTHING (ISOTROPIC MATERIAL PROPERTIES ASSUMED)
C
      if ( typelaw.eq.0 ) then
      endif
C
C.....TREATMENT FOR A TYPE "1", "2", OR "3" CONSTITUTIVE LAW:
C.....STORE THE CORRESPONDING COLUMN NUMBER
C
      if ( (typelaw.eq.1).or.(typelaw.eq.2).or.(typelaw.eq.3) ) then
         iatrb(7,ife) = address
      endif
C
      endif
C
 2100 continue
C
 2000 continue
C
C.....SET UP POINTERS FOR THE COMPOSITE FRAMES
C
      do 3000 iframe=1,ncmpfr
C
C.....GET THE INPUT FRAME NUMBER
C
      framenumber = tycfr(1,iframe)
C
C.....GET THE ADDRESS (COLUMN NUMBER IN FRAME STORAGE)
C
      address = tycfr(2,iframe)
C
C.....LOOP ON FINITE ELEMENTS
C
      do 3100 ife=1,numel
C
C.....SET UP POINTERS IF THE ELEMENT IS RELEVANT:
C.....STORE THE COLUMN NUMBER OF THE FRAME STORAGE
C
      if ( iatrb(5,ife).eq.framenumber ) then
         iatrb(8,ife) = address
      endif
C
 3100 continue
C
 3000 continue
C
C     ------
C     RETURN
C     ------
C
      return
C
C     ---------------
C     ERROR TREATMENT
C     ---------------
C
C.....ERROR-MESSAGE IF DIMENSION [NEN1] IS SMALLER THAN EIGHT
C
  100 continue
      write(*,*)
      write(*,*) "*** FATAL ERROR in Routine COMPATT ***"
      write(*,*) "*** Dimension NEN1 is Smaller than ***"
      write(*,*) "*** Eight: Can not Store Pointers. ***"
      write(*,*) "*** Check the Memory Allocation... ***"
      stop
C
C.....ERROR-MESSAGE IF THE ELEMENT IS NOT A COMPOSITE SHELL
C
  200 continue
      write(*,*)
      write(*,*) "*** FATAL ERROR in Routine COMPATT  ***"
      write(*,*) "*** The Element Attributes Are Not  ***"
      write(*,*) "*** Consistent w/ a Composite Shell ***"
      write(*,*) "*** Check the Input File...         ***"
      stop
C
C     ---
C     END
C     ---
C
      end
C=end of routine "COMPATT"
C========================C
