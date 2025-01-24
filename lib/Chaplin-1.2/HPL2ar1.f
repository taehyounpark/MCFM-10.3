c ------------------------------------------------
      double complex function HPL2ar1(n1,n2,x)
      implicit none
      integer n1,n2,j,bcflag
      double complex x,ris,myi,zp,llzp
      double precision pi, zeta2,ll2,xre

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      bcflag = 0

      j=1+(n2+1)+(n1+1)*3
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      select case(j)
      case(1)                   !-1-1
         
         zp = 1d0-x
         
         ris = -((zp*ll2)/2d0) + (ll2**2)/2d0 + zp**
     &        5*(5d0/384d0 - (ll2)/160d0) + zp**3*(1d0/16d0 - (ll2)/2
     &        4d0) + zp**6*(137d0/23040d0 - (ll2)/384d0) + zp**4*(11d
     &        0/384d0 - (ll2)/64d0) + zp**2*(1d0/8d0 - (ll2)/8d0)
         
      case(2)                   !-10
         
         zp = 1d0-x
         
         ris = -((pi**2)/12d0) + (zp**2)/4d0 + (zp**
     &        3)/6d0 + (5d0*zp**4)/48d0 + (zp**5)/15d0 + (2d0*zp**6)/
     &        45d0
         
      case(3)                   !-11
         
         zp = 1d0-x
         llzp = log(zp)
         
         ris = (pi**2)/12d0 - (ll2**2)/2d0 + zp**5*(
     &        -(1d0/800d0) + (llzp)/160d0) + zp**3*(-(1d0/72d0) + (ll
     &        zp)/24d0) + zp*(-(1d0/2d0) + (llzp)/2d0) + zp**6*(-(1d0
     &        /2304d0) + (llzp)/384d0) + zp**4*(-(1d0/256d0) + (llzp)
     &        /64d0) + zp**2*(-(1d0/16d0) + (llzp)/8d0)
         
      case(4)                   !0-1
         
         zp = 1d0-x
         
         ris = (pi**2)/12d0 - zp*ll2 + zp**2*(1d0/4d
     &        0 - (ll2)/2d0) + zp**3*(5d0/24d0 - (ll2)/3d0) + zp**4*(
     &        1d0/6d0 - (ll2)/4d0) + zp**5*(131d0/960d0 - (ll2)/5d0) 
     &        + zp**6*(661d0/5760d0 - (ll2)/6d0)
         
      case(5)                   !00
         
         zp = 1d0-x
         
         ris = (zp**2)/2d0 + (zp**3)/2d0 + (11d0*zp*
     &        *4)/24d0 + (5d0*zp**5)/12d0 + (137d0*zp**6)/360d0
         
      case(6)                   !01
         
         zp = 1d0-x
         llzp = log(zp)
         
         ris = (pi**2)/6d0 + zp*(-1 + llzp) + zp**2*
     &        (-(1d0/4d0) + (llzp)/2d0) + zp**3*(-(1d0/9d0) + (llzp)/
     &        3d0) + zp**4*(-(1d0/16d0) + (llzp)/4d0) + zp**5*(-(1d0/
     &        25d0) + (llzp)/5d0) + zp**6*(-(1d0/36d0) + (llzp)/6d0)
         
      case(7)                   !1-1
         
         zp = 1d0-x
         llzp = log(zp)
         
         ris = -((pi**2)/12d0) + (zp)/2d0 + (zp**2)/
     &        16d0 + (zp**3)/72d0 + (zp**4)/256d0 + (zp**5)/800d0 + (
     &        zp**6)/2304d0 + (ll2**2)/2d0 - ll2*llzp
         
      case(8)                   !10
         
         zp = 1d0-x
         
         ris = -((pi**2)/6d0) + zp + (zp**2)/4d0 + (
     &        zp**3)/9d0 + (zp**4)/16d0 + (zp**5)/25d0 + (zp**6)/36d0
         
      case(9)                   !11
         
         zp = 1d0-x
         llzp = log(zp)
         
         ris = (llzp**2)/2d0
c     End of expansions around x = +1
      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)
         if (n2.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n2.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n2.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
         endif
      endif        
      
      HPL2ar1=ris
      return
      end function
      
