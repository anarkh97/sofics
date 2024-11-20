      subroutine ffulpvm(neq,pdiag,a,dimblo,nspdim,pivot,listel,nsp,s1,
     &                   s2,v,permut,flag)
      integer neq,dimblo,nspdim,pivot(neq),flag
      real*8 a(*),nsp(neq,nspdim)
      real*8 s1(dimblo,dimblo),s2(dimblo,dimblo),v(dimblo)
      integer pdiag(neq),listel(nspdim),permut(dimblo)
      integer i,j,jj,k,i0
c-----------------------------------------------------------------------
C*************************************************************
C WARNING: This is a modified version of ffulpvm
C          Modification is related to the flag argument
C*************************************************************
C
c  inputs :  neq  =  dimension of the matrix
c           pdiag =  location of diagonal terms in the profile
c             a   =  on input, skyline profile of the matrix, symmetric
c                    LU factorization already performed for neq-dimblo
c                    first rows
c          dimblo =  dimension of the lower diagonal block not yet
c                    factorized
c          nspdim =  nullspace dimension
c          flag   =  flag to prevent nullspace recovery (flag = 11)
c
c  outputs :  listel  =  list of eliminated equations
c               nsp   =  nullspace vectors (need to perform backward
c                        substitution if flag different from 11)
c
c  workspace : v(real),s1(real),s2(real),permut(integer)
c-----------------------------------------------------------------------
c
c  extraction of the lower diagonal block from the profile
c
C      write(6,*) "  neq, dimblo, nspdim : ",neq,dimblo,nspdim
C      write(6,*) "  diagonal entries of Schur complement:"
C      write(6,*) (a(pdiag(i)),i=neq-dimblo+1,neq)
      call vclr(dimblo*dimblo,s1)
      do j=1,dimblo
        jj=neq-dimblo+j
        if(jj.eq.1) then
          i0=max(j-pdiag(jj)+1,1)
        else
          i0=max(j-pdiag(jj)+pdiag(jj-1)+1,1)
        end if
        do i=i0,j
          s1(i,j)=a(pdiag(jj)-j+i)
        end do
      end do
c
c  complement by symmetry
c
      do j=1,dimblo
        do i=j+1,dimblo
          s1(i,j)=s1(j,i)
        end do
      end do
      do i=1,dimblo
        do j=1,dimblo
          s2(i,j)=s1(i,j)
        end do
      end do
c     do i = 1,dimblo
c      write(6,*)'  schur complement: diagonal entry: k =',i,s1(i,i)
c     end do
c
c  test of factorization of the Schur complement to construct the
c  list of eliminated equations
c
      call ldltst(dimblo,s1,nspdim,listel,permut,v)
c
c  modification of the profile
c
c
c  recovery of the upper part of rigid body modes from the profile
c
      if(flag.ne.11) then
         call vclr(neq*nspdim,nsp)
         do j=1,nspdim
           jj=neq-dimblo+listel(j)
           if(jj.eq.1) then
             do i=jj-pdiag(jj)+1,jj
               nsp(i,j)=-a(pdiag(jj)+i-jj)
             end do
           else
             do i=jj-pdiag(jj)+pdiag(jj-1)+1,jj
               nsp(i,j)=-a(pdiag(jj)+i-jj)
             end do
           end if
           do i=jj+1,neq
             if((pdiag(i)-pdiag(i-1)).gt.(i-jj)) then
               nsp(i,j)=-a(pdiag(i)+jj-i)
             end if
           end do
         end do
      end if
c  elimination of columns associated with rigid body modes
      do j=1,nspdim
        jj=neq-dimblo+listel(j)
        if(jj.eq.1) then
          do k=1,pdiag(jj)
            a(k)=0.d0
          end do
        else
          do k=pdiag(jj-1)+1,pdiag(jj)
            a(k)=0.d0
          end do
        end if
      end do
c
c  factorization of the modified Schur complement and copy in the
c  profile of the matrix
c
      do i=1,nspdim
        do j=1,dimblo
          s2(listel(i),j)=0.d0
          s2(j,listel(i))=0.d0
        end do
        s2(listel(i),listel(i))=1.d0
      end do
      call lus(dimblo,s2,v)
c
c  set diagonal entries associated with rigid body modes to 0
c  before copying in the profile of the matrix
c
      do i=1,nspdim
        s2(listel(i),listel(i))=0.d0
      end do      
      do j=1,dimblo
        jj=neq-dimblo+j
        if(jj.eq.1) then
          i0=max(j-pdiag(jj)+1,1)
        else
          i0=max(j-pdiag(jj)+pdiag(jj-1)+1,1)
        end if
        do i=i0,j
          a(pdiag(jj)-j+i)=s2(i,j)
        end do
      end do
c
c  reset diagonal entries associated with rigid body modes to 1
c
      do i=1,nspdim
        s2(listel(i),listel(i))=1.d0
      end do    
c
c  list of eliminated equations in global numbering
c      
      do i=1,nspdim
        listel(i)=listel(i)+neq-dimblo
        pivot(listel(i)) = 0
      end do
cc    write(6,*)'  list of eliminated equations:',(listel(i),i=1,nspdim)
c
c  set to identity of the submatrix of rigid body modes
c
      if(flag.ne.11) then
        do i=1,nspdim
          do j=1,nspdim
            nsp(listel(i),j)=0.d0
          end do
          nsp(listel(i),i)=1.d0
        end do
      endif
c
c  computation of rigid body modes
c
c  complete forward substitution
c
      if(flag.ne.11) then
        do j=1,nspdim
          call desf(dimblo,s2,nsp(neq-dimblo+1,j))
        end do
      endif
      return
      end
c**********************************************************************c
      subroutine desf(dim,a,x) 
      implicit none
      integer dim
      real*8 a(dim,dim),x(dim)
      integer i,j
c
c  forward substitution
c  A is a dense matrix that has been factorized by the
c  symmetric LU method
c
c  forward substitution
c
      x(1)=x(1)/a(1,1)
      do i=2,dim
        do j=1,i-1
          x(i)=x(i)-a(i,j)*x(j)
        end do
        x(i)=x(i)/a(i,i)
      end do
      do i=1,dim
        x(i)=x(i)*a(i,i)
      end do   
      return
      end
c***********************************************************************
      subroutine ldltst(n,a,nspdim,listel,permut,v)
      integer n,nspdim
      real*8 a(n,n),v(n)
      integer listel(nspdim),permut(n)
      integer i,j,k,km,kk
      real*8 piv,dd
c
c  L D Lt factorization with full pivoting in order to find the
c  nspdim equations associated with the smallest pivots
c
c  initialization of permutation
      do i=1,n
        permut(i)=i
      end do
Cc  Crout factorization 
      do k=1,n
Cc  maximum pivot
        piv=a(k,k)
        km=k
        do i=k+1,n
          if(a(i,i).gt.piv) then
            km=i
            piv=a(i,i)
          end if
        end do
c  permutation of rows and columns
c  permut(i) is the initial position of rows and columns number i 
        kk=permut(km)
        permut(km)=permut(k)
        permut(k)=kk
        do i=1,k-1
          v(i)=a(k,i)
          a(k,i)=a(km,i)
          a(km,i)=v(i)
        end do
        do i=k+1,n
          v(i)=a(i,k)
        end do
        do i=k+1,km-1
          a(i,k)=a(km,i)
          a(km,i)=v(i)
        end do
        do i=km+1,n
          a(i,k)=a(i,km)
          a(i,km)=v(i)
        end do
        dd=a(k,k)
        a(k,k)=a(km,km)
        a(km,km)=dd
        piv=a(k,k)
c  computation of the Schur complement
        dd=1.d0/a(k,k)
        do i=k+1,n
          v(i)=a(i,k)*dd
        end do
        do j=k+1,n
          do i=j,n
            a(i,j)=a(i,j)-a(i,k)*v(j)
          end do
        end do
c  recomposition of the column
        do i=k+1,n
          a(i,k)=v(i)
        end do
c
c  return when nspdim steps remain
c
        if(k.eq.(n-nspdim)) then
          do i=1,nspdim
            listel(i)=permut(n-nspdim+i)
          end do
C          write(6,*)'     LDLTST - Assimilated Zero Pivots:'
          do i=1,nspdim
C          write(6,*)'     k = ',i,' pivot = ',a(k+i,k+i)
          end do
          return
        end if
      end do
      return
      end
**********************************************************************c
      subroutine lus(dim,a,v) 
      integer dim
      real*8 a(dim,dim),v(dim)
      integer i,j,k
c
c  symmetric LU factorization 
c
c  L D Lt factorization with storage of U = D Lt or Ut = L D
c
      do i=2,dim
        do j=1,i-1
          do k=1,j-1
            a(i,j)=a(i,j)-(v(k)*a(j,k))
          end do
          IF(a(j,j).eq.0) WRITE(6,*)'LUS - ZERO j = ',j
         v(j)=a(i,j)/a(j,j)
         a(j,i)=a(i,j)
        end do
        do j=1,i-1
          a(i,i)=a(i,i)-(a(i,j)*v(j))
        end do
      end do
      return
      end
cc**********************************************************************c
      subroutine vclr(n,x)
      integer n
      real*8 x(n)
      integer i
      do i=1,n
        x(i)=0.d0
      end do
      return
      end
