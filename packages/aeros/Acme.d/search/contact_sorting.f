C $Id$
      subroutine contact_indexx(n,arrin,indx,ndim)

c  create an index array so that ARRIN(INDX(j)) is in ascending order
c  copied for the numerical recipes book 2nd edition.
      
      implicit real*8 (a-h,o-z)

      dimension arrin(ndim),indx(ndim)
      integer nstack
      parameter (m=7,nstack=100)
      integer istack(nstack)

      do j=1,n
        indx(j) = j
      end do
      jstack = 0
      l = 1
      ir = n
    1 if( ir-l.lt.m )then
        do j = l+1,ir
          indxt=indx(j)
          a=arrin(indxt)
          do i=j-1,l,-1
            if (arrin(indx(i)).le.a) goto 2
            indx(i+1)=indx(i)
          enddo
          i=l-1
    2     indx(i+1)=indxt
        enddo
        if (jstack.eq.0) goto 6
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if (arrin(indx(l)).gt.arrin(indx(ir))) then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if (arrin(indx(l+1)).gt.arrin(indx(ir))) then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if (arrin(indx(l)).gt.arrin(indx(l+1))) then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arrin(indxt)
    3   continue
        i=i+1
        if (arrin(indx(i)).lt.a) goto 3
    4   continue
        j=j-1
        if (arrin(indx(j)).gt.a) goto 4
        if (j.lt.i) goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
    5   indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if (jstack.gt.NSTACK) then
          print*,'NSTACK too small in indexx'
          goto 6
        endif
        if (ir-i+1.ge.j-l) then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
    6 continue

      return
      end




      subroutine contact_indexx_float(n,arrin,indx,ndim)

c  create an index array so that ARRIN(INDX(j)) is in ascending order
c  copied for the numerical recipes book 2nd edition.
      
      implicit real*4 (a-h,o-z)

      dimension arrin(ndim)
      integer indx(ndim)
      integer nstack
      parameter (m=7,nstack=100)
      integer istack(nstack)

      do j=1, n
        indx(j) = j
      end do
      jstack = 0
      l = 1
      ir = n
    1 if( ir-l.lt.m )then
        do j = l+1,ir
          indxt=indx(j)
          a=arrin(indxt)
          do i=j-1,l,-1
            if (arrin(indx(i)).le.a) goto 2
            indx(i+1)=indx(i)
          enddo
          i=l-1
    2     indx(i+1)=indxt
        enddo
        if (jstack.eq.0) goto 6
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if (arrin(indx(l)).gt.arrin(indx(ir))) then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if (arrin(indx(l+1)).gt.arrin(indx(ir))) then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if (arrin(indx(l)).gt.arrin(indx(l+1))) then
          itemp=indx(l)
          indx(l)=indx(l+1)
          indx(l+1)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l+1)
        a=arrin(indxt)
    3   continue
        i=i+1
        if (arrin(indx(i)).lt.a) goto 3
    4   continue
        j=j-1
        if (arrin(indx(j)).gt.a) goto 4
        if (j.lt.i) goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
    5   indx(l+1)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if (jstack.gt.NSTACK) then
          print*,'NSTACK too small in indexx'
          goto 6
        endif
        if (ir-i+1.ge.j-l) then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
    6 continue

      return
      end



      subroutine contact_rank(n,indx,irank,ndim)

c create a rank array from an index array
c copied from the numberical recipes book

      dimension indx(ndim),irank(ndim)
      
      do i=1,n
        irank(indx(i)) = i
      end do
      
      return
      end


      subroutine contact_get_bound(NE,X,IND,NP,XMIN,XMAX,NDIM,ILO,IUP,
     *  ND2,SCR)

      implicit real*8 (a-h,o-z)
      parameter( neblk = 1 )

      DIMENSION X(NDIM),IND(NDIM),XMIN(NE),XMAX(NE),ILO(NE),IUP(NE),
     &  SCR(ND2)
C
C  find the elements in a sorted array x whose values fall in the
C  interval between xmin and xmax. No elements having
C  values equal to xmin or xmax are included. Since the array is
C  sorted, the elements can be specified by the upper and
C  lower element numbers in the range.
C
C     |    X(ILO)   . . . . . . . . . . . . .   X(IUP)    |
C    XMIN                                               XMAX       X>
C
C  It is assumed that the array x has been sorted in increasing order,
C  but the elements have not been moved.
C  The sorted list is determined by the array indx,
C  which positions the original unsorted x array elements
C  in the sorted list. Thus, the 5th element in the sorted list is
C    X(IND(5))
C
C  INPUT
C  X      -  array in unsorted order
C  IND    -  index array giving the element order in the sorted list
C  NP     -  the number of particles in the list
C  XMIN   -  the lower limit of the interval
C  XMAX   -  the upper limit of the interval
C  NDIM   -  the dimension of the arrays
C
C  OUTPUT
C  ILO    -  the first element in the sorted list .gt. xmin
C  IUP    -  the last element in the sorted list .lt. xmax
C
C
      do j = 1, ne, neblk
        n = min(neblk,ne-j+1)
c
c  search to find the first element .ge. xmin
        call contact_srchge(n,x,ind,xmin(j),1,np,ndim,ilo(j))
c
c  search to find the first element .gt. xmax
        call contact_srchgt(n,x,ind,xmax(j),ilo(j),np,ndim,iup(j))
c
c  the previous element is the last one .lt. xmax
        do jj = 1, n
          iup(j+jj-1)=iup(j+jj-1) - 1
        end do
c
      end do
c
      return
      end



      subroutine contact_srchge(ne,X,IND,XV,IMIN,IMAX,NDIM,I)

      implicit real*8 (a-h,o-z)

      DIMENSION X(NDIM),IND(NDIM), XV(*),I(32),IL(32),
     *          IU(32),IT(32),XTST(32),INDX1(32),
     *          INDX2(32)
C
C  perform a binary search to find the element number i
c  of a sorted array for which all elements at i or above are
c  greater or equal to than some value xv,
c  while all elements below i are less than the value.
C
C       X(I-2)     X(I-1)    X(I)    X(I+1)   X(I+2)
C                      XV                             X>
C
C  Assumed that the array x has been sorted in increasing order,
c  but the elements have not been moved.
c  The sorted list is determined by the array indx,
c  which positions the original unsorted x array elements
c  in the sorted list. Thus, the 5th element in the sorted list is
C    X(IND(5))
C
C  INPUT
C  X      -  ARRAY IN UNSORTED ORDER
C  IND    -  INDEX ARRAY GIVING THE ELEMENT ORDER IN THE SORTED LIST
C  XV     -  X VALUE TO TEST AGAINST
C  IMIN   -  THE LOWEST NUMBERED POSITION IN THE SORTED LIST TO TEST
C  IMAX   -  THE HIGHEST NUMBERED POSITION IN THE SORTED LIST TO TEST
C  NDIM   -  THE DIMENSION OF THE ARRAYS
C
C  OUTPUT
C  I      -  THE FIRST POSITION IN THE SORTED LIST .GT. XV
C
C
C  EXECUTE THE BINARY SEARCH
C
      if (imax.eq.0.or.ndim.eq.0) then
         do j = 1, ne
            i(j) = 0
         enddo
         return
      endif
      DO 25 J = 1, NE
        IL(J) = IMIN
        IU(J) = IMAX
        INDX1(J) = J
 25   CONTINUE
        ILOOP = NE
 1000 CONTINUE

      DO 50 JJ = 1, ILOOP
        INDX2(JJ) = INDX1(JJ)
        J = INDX1(JJ)      
        IT(J) =  (IU(J) + IL(J)) / 2 
 50   CONTINUE
      DO 35 J = 1, NE
       XTST(J) = X( IND(IT(J)) )
 35   CONTINUE

      IF ( ILOOP .GT. 64) THEN
C  DIR$ IVDEP
      ILP = 0
      DO 60 JJ = 1, ILOOP
        J = INDX2(JJ)
        IF( XTST(J) .LT. XV(J) )THEN
          IL(J) =IT(J) + 1
        ELSE
          IU(J) =IT(J) - 1
        ENDIF
        IF( IL(J) .LE. IU(J)) THEN
          ILP = ILP + 1
          INDX1(ILP) = J 
        ENDIF
 60   CONTINUE
      ELSE
C  DIR$ SHORT LOOP
C  DIR$ IVDEP
      ILP = 0
      DO 51 JJ = 1, ILOOP
        J = INDX2(JJ)      
        IF( XTST(J) .LT. XV(J) )THEN
          IL(J) =IT(J) + 1
        ELSE
          IU(J) =IT(J) - 1
        ENDIF
        IF( IL(J) .LE. IU(J)) THEN
          ILP = ILP + 1
          INDX1(ILP) = J 
        ENDIF
 51   CONTINUE
      ENDIF

         ILOOP = ILP

      IF(ILOOP .NE. 0 )GO TO 1000
C  range had narrowed to 1 location. However, the point last tested
C  could be above, below, or on the search point. Check for proper case
      DO 200 J = 1, NE
        IF( XTST(J) .LT. XV(J) )THEN
          I(J) = IT(J) + 1
        ELSE
          I(J) = IT(J) 
        ENDIF
 200  CONTINUE

      RETURN
      END
      


      subroutine contact_srchgt(ne,X,IND,XV,IMIN,IMAX,NDIM,I)

      implicit real*8 (a-h,o-z)

      DIMENSION X(NDIM),IND(NDIM), XV(*),I(*),IL(32),
     *          IU(32),IT(32),XTST(32),INDX1(32),
     *          INDX2(32),IMIN(*)
C
C  perform a binary search to find the element number i
c  of a sorted array for which all elements at i or above are
c  greater than some value xv,
c  while all elements below i are less than or equal to the value.
C
C       X(I-2)     X(I-1)    X(I)    X(I+1)   X(I+2)
C                      XV                             X>
C
C  Assumed that the array x has been sorted in increasing order,
c  but the elements have not been moved.
c  The sorted list is determined by the array indx,
c  which positions the original unsorted x array elements
c  in the sorted list. Thus, the 5th element in the sorted list is
C    X(IND(5))
C
C  INPUT
C  X      -  ARRAY IN UNSORTED ORDER
C  IND    -  INDEX ARRAY GIVING THE ELEMENT ORDER IN THE SORTED LIST
C  XV     -  X VALUE TO TEST AGAINST
C  IMIN   -  THE LOWEST NUMBERED POSITION IN THE SORTED LIST TO TEST
C  IMAX   -  THE HIGHEST NUMBERED POSITION IN THE SORTED LIST TO TEST
C  NDIM   -  THE DIMENSION OF THE ARRAYS
C
C  OUTPUT
C  I      -  THE FIRST POSITION IN THE SORTED LIST .GE. XV
C
C
C  EXECUTE THE BINARY SEARCH
C
      if (imax.eq.0.or.ndim.eq.0) then
         do j = 1, ne
            i(j) = 0
         enddo
         return
      endif
      DO 25 J = 1, NE
        IL(J) = IMIN(J)
        IU(J) = IMAX
        INDX1(J) = J
 25   CONTINUE
        ILOOP = NE
 1000 CONTINUE

      DO 50 JJ = 1, ILOOP
        INDX2(JJ) = INDX1(JJ)
        J = INDX1(JJ)      
        IT(J) =  (IU(J) + IL(J)) / 2 
 50   CONTINUE

      DO 35 J = 1, NE
       XTST(J) = X( IND(IT(J)) )
 35   CONTINUE

      IF ( ILOOP .GT. 64) THEN
C  DIR$ IVDEP
      ILP = 0
      DO 60 JJ = 1, ILOOP
        J = INDX2(JJ)
        IF( XTST(J) .LE. XV(J) )THEN
          IL(J) =IT(J) + 1
        ELSE
          IU(J) =IT(J) - 1
        ENDIF
        IF( IL(J) .LE. IU(J)) THEN
          ILP = ILP + 1
          INDX1(ILP) = J 
        ENDIF
 60   CONTINUE
      ELSE
C  DIR$ SHORT LOOP
C  DIR$ IVDEP
      ILP = 0
      DO 51 JJ = 1, ILOOP
        J = INDX2(JJ)      
        IF( XTST(J) .LE. XV(J) )THEN
          IL(J) =IT(J) + 1
        ELSE
          IU(J) =IT(J) - 1
        ENDIF
        IF( IL(J) .LE. IU(J)) THEN
          ILP = ILP + 1
          INDX1(ILP) = J 
        ENDIF
 51   CONTINUE
      ENDIF

         ILOOP = ILP

      IF(ILOOP .NE. 0 )GO TO 1000
C  RANGE HAD NARROWED TO 1 LOCATION. HOWEVER, THE POINT LAST TESTED
C  COULD BE ABOVE, BELOW, OR ON THE SEARCH POINT. CHECK FOR PROPER CASE
      DO 200 J = 1, NE
        IF( XTST(J) .LE. XV(J) )THEN
          I(J) = IT(J) + 1
        ELSE
          I(J) = IT(J) 
        ENDIF
 200  CONTINUE

      RETURN
      END

      subroutine contact_get_bound_t(NE,X,IND,NP,XMIN,XMAX,NDIM,ILO,IUP,
     *  ND2,SCR,indexoffset)

      implicit real*8 (a-h,o-z)
      parameter( neblk = 1 )

      DIMENSION X(NDIM),IND(NDIM),XMIN(NE),XMAX(NE),ILO(NE),IUP(NE),
     &  SCR(ND2)
      integer indexoffset
C
C  find the elements in a sorted array x whose values fall in the
C  interval between xmin and xmax. No elements having
C  values equal to xmin or xmax are included. Since the array is
C  sorted, the elements can be specified by the upper and
C  lower element numbers in the range.
C
C     |    X(ILO)   . . . . . . . . . . . . .   X(IUP)    |
C    XMIN                                               XMAX       X>
C
C  It is assumed that the array x has been sorted in increasing order,
C  but the elements have not been moved.
C  The sorted list is determined by the array indx,
C  which positions the original unsorted x array elements
C  in the sorted list. Thus, the 5th element in the sorted list is
C    X(IND(5))
C
C  INPUT
C  X      -  array in unsorted order
C  IND    -  index array giving the element order in the sorted list
C  NP     -  the number of particles in the list
C  XMIN   -  the lower limit of the interval
C  XMAX   -  the upper limit of the interval
C  NDIM   -  the dimension of the arrays
C
C  OUTPUT
C  ILO    -  the first element in the sorted list .gt. xmin
C  IUP    -  the last element in the sorted list .lt. xmax
C
C
      il = indexoffset+1
      do j = 1, ne, neblk
        n = min(neblk,ne-j+1)
c
c  search to find the first element .ge. xmin
        call contact_srchge(n,x,ind,xmin(j),il,np,ndim,ilo(j))
c
c  search to find the first element .gt. xmax
        call contact_srchgt(n,x,ind,xmax(j),ilo(j),np,ndim,iup(j))
c
c  the previous element is the last one .lt. xmax
        do jj = 1, n
          iup(j+jj-1)=iup(j+jj-1) - 1
        end do
c
      end do
c
      return
      end



      subroutine contact_make_list(NUMPTS,IND,IRNK2,IUP,ILO,
     *  IE,LIST,NLIST,NBLKSZ,NSPC)

C - vector make list (3d) 
C given a list of particles (ie  their index and rank) f1ind
C the list of particles within the bounds set by xmin and xmax
C for the ie'th particle in the vector block

C     numpts integer   number of points to be searched
C     ind    integer   order index
c     irnk2  integer   rank 
c     iup    integer   scracth (nblksz long)
c     ilo    integer   scracth (nblksz long)
C     ie     integer   particle number
C     list   integer   list of found particles
c     nlist  integer   number of particles found
C     nblksz integer   block size of iup and ilo blocks
c     nspc   integer   number of spacial coord. (number of dimensions)
C
C --------------------------------------

      DIMENSION IND(NUMPTS,NSPC),IRNK2(NUMPTS,NSPC,*),
     *  IUP(NBLKSZ,NSPC),ILO(NBLKSZ,NSPC),LIST(NUMPTS)


C loop over all master surfaces and build a list of slave nodes that can
c interact with each master surface
      J = IE
      NLIST = 0
      IF( NSPC .EQ. 1)THEN
c============================== o n e   -  d ======================
        NUM1 = IUP(J,1) - ILO(J,1) + 1
        ILOW = ILO(J,1)
        IUPR = IUP(J,1)      
        DO 101 I1 = ILOW, IUPR
          NLIST = NLIST +1
          LIST(NLIST) = IND(I1,1)
 101    CONTINUE

      ELSE IF( NSPC .EQ. 2 )THEN
c============================== t w o   -  d ======================
        NUM1 = IUP(J,1) - ILO(J,1) + 1
        NUM2 = IUP(J,2) - ILO(J,2) + 1

C do we have a list ? 
        IF( NUM2.LE.0 .OR. NUM1.LE.0 ) RETURN

C select which list is the smallest
        IF( NUM1 .LE. NUM2 )THEN
          IXYZ = 1
          IY = 2
          NUM = NUM1
        ELSE
          IXYZ = 2
          IY = 1
          NUM = NUM2
        ENDIF

        ILOW = ILO(J,IXYZ)
        IUPR = IUP(J,IXYZ)      

C first test
        IF( NUM .GT. 64 ) THEN
          DO 201 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND. 
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              NLIST = NLIST +1
              LIST(NLIST) = IND(I1,IXYZ)
            ENDIF
 201      CONTINUE
        ELSE
C  DIR$ SHORT LOOP
          DO 202 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND. 
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              NLIST = NLIST +1
              LIST(NLIST) = IND(I1,IXYZ)
            ENDIF
 202      CONTINUE
        ENDIF

      ELSE IF( NSPC .EQ. 3 )THEN
c============================== t h r e e   -  d ======================

        NUM1 = IUP(J,1) - ILO(J,1) + 1
        NUM2 = IUP(J,2) - ILO(J,2) + 1
        NUM3 = IUP(J,3) - ILO(J,3) + 1

C do we have a list ? 
        IF( NUM3.LE.0 .OR. NUM2.LE.0 .OR. NUM1.LE.0 ) RETURN

C select which list is the smallest
        IF( NUM1 .LE. NUM2 .AND. NUM1 .LE. NUM3 )THEN
          IXYZ = 1
          IY = 2
          IZ = 3
          NUM = NUM1
        ELSEIF( NUM2 .LE. NUM1 .AND. NUM2 .LE. NUM3 )THEN
          IXYZ = 2
          IY = 1
          IZ = 3
          NUM = NUM2
        ELSE
          IXYZ = 3
          IY = 1
          IZ = 2
          NUM = NUM3
        ENDIF
        ILOW = ILO(J,IXYZ)
        IUPR = IUP(J,IXYZ)      
        if (ilow.eq.0) then
          nlist = 0
          return
        endif
C first test
        DO 311 I1 = ILOW, IUPR
          IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND. 
     *      IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) .AND. 
     *      IRNK2(I1,IXYZ,2) .GE. ILO(J,IZ) .AND.
     *      IRNK2(I1,IXYZ,2) .LE. IUP(J,IZ) )THEN
            NLIST = NLIST + 1
            LIST(NLIST) = IND(I1,IXYZ)
          ENDIF
 311    CONTINUE
      ENDIF

      RETURN
      END
