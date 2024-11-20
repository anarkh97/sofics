       subroutine transform(xl,yl,zl,xg,yg,zg,str)
C
C**********************************************************************C
C
C Purpose: to form a transformation matrix from local coordinates to
C          global coordinates
C
C input variables:
C      xl = x local unit vector
C      yl = y local unit vector
C      zl = z local unit vector
C      xg = x global unit vector
C      yg = y global unit vector
C      zg = z global unit vector
C      str = stress/strain 6x1 vector
C            sigmaxx, sigmayy, sigmazz, sigma12, sigma23, sigma13
C
C local variables:
C      l1 = direction cosine between xl and xg
C      l2 = direction cosine between xl and yg
C      l3 = direction cosine between xl and zg
C      m1 = direction cosine between yl and xg
C      m2 = direction cosine between yl and yg
C      m3 = direction cosine between yl and zg
C      n1 = direction cosine between zl and xg
C      n2 = direction cosine between zl and yg
C      n3 = direction cosine between zl and zg
C      t  = transformation matrix from local to global
C
C**********************************************************************C
C
C Declarations
C
        double precision xl(3),yl(3),zl(3)
        double precision xg(3),yg(3),zg(3)
        double precision str(6)
C
C Local Declarations
C
        double precision l1,l2,l3
        double precision m1,m2,m3
        double precision n1,n2,n3
        double precision t(6,6), s(6)
C
C Copy stress/strain values to a temporary array
C
        s(1) = str(1)
        s(2) = str(2)
        s(3) = str(3)
        s(4) = str(4)
        s(5) = str(5)
        s(6) = str(6)
C
C Compute direction cosines
C       
        l1 = xg(1)*xl(1) + xg(2)*xl(2) + xg(3)*xl(3)
        l2 = yg(1)*xl(1) + yg(2)*xl(2) + yg(3)*xl(3)
        l3 = zg(1)*xl(1) + zg(2)*xl(2) + zg(3)*xl(3)
        
        m1 = xg(1)*yl(1) + xg(2)*yl(2) + xg(3)*yl(3)
        m2 = yg(1)*yl(1) + yg(2)*yl(2) + yg(3)*yl(3)
        m3 = zg(1)*yl(1) + zg(2)*yl(2) + zg(3)*yl(3)
        
        n1 = xg(1)*zl(1) + xg(2)*zl(2) + xg(3)*zl(3)
        n2 = yg(1)*zl(1) + yg(2)*zl(2) + yg(3)*zl(3)
        n3 = zg(1)*zl(1) + zg(2)*zl(2) + zg(3)*zl(3)
C
C Construct the 6x6 transformation matrix
C       
        t(1,1) = l1*l1
        t(1,2) = m1*m1
        t(1,3) = n1*n1
        t(1,4) = 2.0*l1*m1
        t(1,5) = 2.0*m1*n1
        t(1,6) = 2.0*n1*l1
C       
        t(2,1) = l2*l2
        t(2,2) = m2*m2
        t(2,3) = n2*n2
        t(2,4) = 2.0*l2*m2
        t(2,5) = 2.0*m2*n2
        t(2,6) = 2.0*n2*l2
C       
        t(3,1) = l3*l3
        t(3,2) = m3*m3
        t(3,3) = n3*n3
        t(3,4) = 2.0*l3*m3
        t(3,5) = 2.0*m3*n3
        t(3,6) = 2.0*n3*l3
C       
        t(4,1) = l1*l2
        t(4,2) = m1*m2
        t(4,3) = n1*n2
        t(4,4) = l1*m2 + l2*m1
        t(4,5) = m1*n2 + m2*n1
        t(4,6) = n1*l2 + n2*l1
C       
        t(5,1) = l2*l3
        t(5,2) = m2*m3
        t(5,3) = n2*n3
        t(5,4) = l2*m3 + l3*m2
        t(5,5) = m2*n3 + m3*n2
        t(5,6) = n2*l3 + n3*l2
C       
        t(6,1) = l3*l1
        t(6,2) = m3*m1
        t(6,3) = n3*n1
        t(6,4) = l3*m1 + l1*m3
        t(6,5) = m3*n1 + m1*n3
        t(6,6) = n3*l1 + n1*l3
C
C Perform the multiplication {str'} = T{str}
C       
       str(1) = t(1,1)*s(1) + t(1,2)*s(2) + t(1,3)*s(3) +
     &          t(1,4)*s(4) + t(1,5)*s(5) + t(1,6)*s(6)
     
       str(2) = t(2,1)*s(1) + t(2,2)*s(2) + t(2,3)*s(3) +
     &          t(2,4)*s(4) + t(2,5)*s(5) + t(2,6)*s(6)
     
       str(3) = t(3,1)*s(1) + t(3,2)*s(2) + t(3,3)*s(3) +
     &          t(3,4)*s(4) + t(3,5)*s(5) + t(3,6)*s(6)
     
       str(4) = t(4,1)*s(1) + t(4,2)*s(2) + t(4,3)*s(3) +
     &          t(4,4)*s(4) + t(4,5)*s(5) + t(4,6)*s(6)
     
       str(5) = t(5,1)*s(1) + t(5,2)*s(2) + t(5,3)*s(3) +
     &          t(5,4)*s(4) + t(5,5)*s(5) + t(5,6)*s(6)
     
       str(6) = t(6,1)*s(1) + t(6,2)*s(2) + t(6,3)*s(3) +
     &          t(6,4)*s(4) + t(6,5)*s(5) + t(6,6)*s(6)
     
        
        return
        end
