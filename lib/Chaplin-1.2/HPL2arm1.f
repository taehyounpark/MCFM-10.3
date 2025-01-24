c ------------------------------------------------
      double complex function HPL2arm1(n1,n2,x)
      implicit none
      integer n1,n2,j,bcflag,s,szp
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

         zp = x+1d0
         llzp = log(zp)
         
         ris = (llzp**2)/2d0
         
      case(2)                   !-10
         
         zp = x+1d0
         llzp = log(zp)
         szp = s(zp)
         
         ris = (pi**2)/6d0 - zp - (zp**2)/4d0 - (zp*
     &        *3)/9d0 - (zp**4)/16d0 - (zp**5)/25d0 - (zp**6)/36d0 - 
     &        (zp**7)/49d0 - (zp**8)/64d0 - (zp**9)/81d0 + myi*pi*szp
     &        *llzp
         
      case(3)                   !-11
         
         zp = x+1d0
         llzp = log(zp)
         
         ris = -((pi**2)/12d0) + (zp)/2d0 + (zp**2)/
     &        16d0 + (zp**3)/72d0 + (zp**4)/256d0 + (zp**5)/800d0 + (
     &        zp**6)/2304d0 + (zp**7)/6272d0 + (zp**8)/16384d0 + (zp*
     &        *9)/41472d0 + (ll2**2)/2d0 - ll2*llzp
         
      case(4)                   !0-1
         
         zp = x+1d0
         llzp = log(zp)
         
         ris = -((pi**2)/6d0) + zp*(1 - llzp) + zp**
     &        2*(1d0/4d0 - (llzp)/2d0) + zp**3*(1d0/9d0 - (llzp)/3d0)
     &        + zp**4*(1d0/16d0 - (llzp)/4d0) + zp**5*(1d0/25d0 - (l
     &        lzp)/5d0) + zp**6*(1d0/36d0 - (llzp)/6d0) + zp**7*(1d0/
     &        49d0 - (llzp)/7d0) + zp**8*(1d0/64d0 - (llzp)/8d0) + zp
     &        **9*(1d0/81d0 - (llzp)/9d0)
         
      case(5)                   !00
         
         zp = x+1d0
         szp = s(zp)
         
         ris = -((pi**2)/2d0) - myi*pi*szp*zp + (1d0
     &        /2d0 - (myi*pi*szp)/2d0)*zp**2 + (1d0/2d0 - (myi*pi*szp
     &        )/3d0)*zp**3 + (11d0/24d0 - (myi*pi*szp)/4d0)*zp**4 + (
     &        5d0/12d0 - (myi*pi*szp)/5d0)*zp**5 + (137d0/360d0 - (my
     &        i*pi*szp)/6d0)*zp**6 + (7d0/20d0 - (myi*pi*szp)/7d0)*zp
     &        **7 + (363d0/1120d0 - (myi*pi*szp)/8d0)*zp**8 + (761d0/
     &        2520d0 - (myi*pi*szp)/9d0)*zp**9
         
      case(6)                   !01
         
         zp = x+1d0
         
         ris = -((pi**2)/12d0) + zp*ll2 + zp**2*(-(1
     &        d0/4d0) + (ll2)/2d0) + zp**3*(-(5d0/24d0) + (ll2)/3d0) 
     &        + zp**4*(-(1d0/6d0) + (ll2)/4d0) + zp**5*(-(131d0/960d0
     &        ) + (ll2)/5d0) + zp**6*(-(661d0/5760d0) + (ll2)/6d0) + 
     &        zp**7*(-(1327d0/13440d0) + (ll2)/7d0) + zp**8*(-(1163d0
     &        /13440d0) + (ll2)/8d0) + zp**9*(-(148969d0/1935360d0) +
     &        (ll2)/9d0)
         
      case(7)                   !1-1
         
         zp = x+1d0
         llzp = log(zp)
         
         ris = (pi**2)/12d0 - (ll2**2)/2d0 + zp**5*(
     &        -(1d0/800d0) + (llzp)/160d0) + zp**8*(-(1d0/16384d0) + 
     &        (llzp)/2048d0) + zp**3*(-(1d0/72d0) + (llzp)/24d0) + zp
     &        *(-(1d0/2d0) + (llzp)/2d0) + zp**6*(-(1d0/2304d0) + (ll
     &        zp)/384d0) + zp**9*(-(1d0/41472d0) + (llzp)/4608d0) + z
     &        p**4*(-(1d0/256d0) + (llzp)/64d0) + zp**7*(-(1d0/6272d0
     &        ) + (llzp)/896d0) + zp**2*(-(1d0/16d0) + (llzp)/8d0)
         
      case(8)                   !10
         
         zp = x+1d0
         szp = s(zp)
         
         ris = (pi**2)/12d0 + (myi*pi*szp*zp)/2d0 + 
     &        (-(1d0/4d0) + (myi*pi*szp)/8d0)*zp**2 + (-(1d0/6d0) + (
     &        myi*pi*szp)/24d0)*zp**3 + (-(5d0/48d0) + (myi*pi*szp)/6
     &        4d0)*zp**4 + (-(1d0/15d0) + (myi*pi*szp)/160d0)*zp**5 +
     &        (-(2d0/45d0) + (myi*pi*szp)/384d0)*zp**6 + (-(13d0/420
     &        d0) + (myi*pi*szp)/896d0)*zp**7 + (-(151d0/6720d0) + (m
     &        yi*pi*szp)/2048d0)*zp**8 + (-(16d0/945d0) + (myi*pi*szp
     &        )/4608d0)*zp**9 - myi*pi*szp*ll2
         
      case(9)                   !11
         
         zp = x+1d0
         
         ris = -((zp*ll2)/2d0) + (ll2**2)/2d0 + zp**
     &        5*(5d0/384d0 - (ll2)/160d0) + zp**8*(363d0/286720d0 - (
     &        ll2)/2048d0) + zp**3*(1d0/16d0 - (ll2)/24d0) + zp**6*(1
     &        37d0/23040d0 - (ll2)/384d0) + zp**9*(761d0/1290240d0 - 
     &        (ll2)/4608d0) + zp**4*(11d0/384d0 - (ll2)/64d0) + zp**7
     &        *(7d0/2560d0 - (ll2)/896d0) + zp**2*(1d0/8d0 - (ll2)/8d
     &        0)
c     End of expansions around x = -1
      end select
c     --- set the imaginary part back to zero if it has been modified to
c     --- get the branch cuts right (and should be zero).
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
      
      HPL2arm1=ris
      return
      end function
