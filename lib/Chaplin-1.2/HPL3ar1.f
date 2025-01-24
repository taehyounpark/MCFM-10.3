c ------------------------------------------------
      double complex function HPL3ar1(n1,n2,n3,x)
      implicit none
      integer n1,n2,n3,j,bcflag
      double complex x,ris,myi,zp,llzp
      double precision pi, zeta2, zeta3,ll2,xre

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      bcflag = 0

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9
      ris = dcmplx(0d0,0d0)

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  

c This was file contains the Taylor 
c expansions around x = +1
c The expansion parameter is zp = 1-x
      select case(j)
      
         case(1)            !-1-1-1

            zp = 1d0-x

            ris = -((zp*ll2**2)/4d0) + (ll2**3)/6d0 + z
     &p**4*(-(1d0/64d0) + (11d0*ll2)/384d0 - (ll2**2)/128d0) 
     &+ zp**2*((ll2)/8d0 - (ll2**2)/16d0) + zp**5*(-(7d0/768d
     &0) + (5d0*ll2)/384d0 - (ll2**2)/320d0) + zp**3*(-(1d0/4
     &8d0) + (ll2)/16d0 - (ll2**2)/48d0) + zp**6*(-(5d0/1024d
     &0) + (137d0*ll2)/23040d0 - (ll2**2)/768d0)

         case(2)            !-1-10

            zp = 1d0-x

            ris = (pi**2*zp)/24d0 + (pi**2*zp**2)/96d0 
     &+ (-(1d0/24d0) + (pi**2)/288d0)*zp**3 + ((pi**2)/768d0 
     &- 7d0/192d0)*zp**4 + ((pi**2)/1920d0 - 1d0/40d0)*zp**5 
     &+ (-(23d0/1440d0) + (pi**2)/4608d0)*zp**6 - (pi**2*ll2)
     &/12d0 + (zeta3)/8d0

         case(3)            !-1-11

            zp = 1d0-x
            llzp = log(zp)

            ris = -((ll2**3)/6d0) + zp*(-((pi**2)/24d0)
     & + (ll2**2)/4d0) + zp**3*(-((pi**2)/288d0) + 7d0/96d0 +
     & (ll2**2)/48d0 - (llzp)/16d0) + zp**6*(2213d0/460800d0 
     &- (pi**2)/4608d0 + (ll2**2)/768d0 - (137d0*llzp)/23040d
     &0) + zp**4*(131d0/4608d0 - (pi**2)/768d0 + (ll2**2)/128
     &d0 - (11d0*llzp)/384d0) + zp**5*(-((pi**2)/1920d0) + 53
     &d0/4608d0 + (ll2**2)/320d0 - (5d0*llzp)/384d0) + zp**2*
     &(3d0/16d0 - (pi**2)/96d0 + (ll2**2)/16d0 - (llzp)/8d0) 
     &+ (zeta3)/8d0

         case(4)            !-10-1

            zp = 1d0-x

            ris = -((pi**2*zp)/24d0) + (pi**2*ll2)/12d0
     & + zp**5*(-((pi**2)/1920d0) - 1d0/30d0 + (ll2)/15d0) + 
     &zp**6*(-((pi**2)/4608d0) - 97d0/3840d0 + (2d0*ll2)/45d0
     &) + zp**2*(-((pi**2)/96d0) + (ll2)/4d0) + zp**4*(-(1d0/
     &24d0) - (pi**2)/768d0 + (5d0*ll2)/48d0) + zp**3*(-(1d0/
     &24d0) - (pi**2)/288d0 + (ll2)/6d0) - (zeta3)/4d0

         case(5)            !-100

            zp = 1d0-x

            ris = -((zp**3)/12d0) - (3d0*zp**4)/32d0 - 
     &(zp**5)/12d0 - (5d0*zp**6)/72d0 + (3d0*zeta3)/4d0

         case(6)            !-101

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*zp)/12d0) + (pi**2*ll2)/6d0 
     &+ zp**5*(79d0/1800d0 - (pi**2)/960d0 - (llzp)/15d0) + z
     &p**6*(-((pi**2)/2304d0) + 169d0/7200d0 - (2d0*llzp)/45d
     &0) + zp**2*(-((pi**2)/48d0) + 3d0/8d0 - (llzp)/4d0) + z
     &p**4*(25d0/288d0 - (pi**2)/384d0 - (5d0*llzp)/48d0) + z
     &p**3*(-((pi**2)/144d0) + 13d0/72d0 - (llzp)/6d0) - (5d0
     &*zeta3)/8d0

         case(7)            !-11-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**2*ll2)/12d0 - (ll2**3)/6d0 + zp*
     &*5*((pi**2)/1920d0 - 41d0/4608d0 - (ll2)/800d0 - (ll2**
     &2)/320d0 + (ll2*llzp)/160d0) + zp**3*((pi**2)/288d0 - 5
     &d0/96d0 - (ll2)/72d0 - (ll2**2)/48d0 + (ll2*llzp)/24d0)
     & + zp*((pi**2)/24d0 - (ll2)/2d0 - (ll2**2)/4d0 + (ll2*l
     &lzp)/2d0) + zp**6*((pi**2)/4608d0 - 5269d0/1382400d0 - 
     &(ll2)/2304d0 - (ll2**2)/768d0 + (ll2*llzp)/384d0) + zp*
     &*4*(-(49d0/2304d0) + (pi**2)/768d0 - (ll2)/256d0 - (ll2
     &**2)/128d0 + (ll2*llzp)/64d0) + zp**2*(-(1d0/8d0) + (pi
     &**2)/96d0 - (ll2)/16d0 - (ll2**2)/16d0 + (ll2*llzp)/8d0
     &) - (zeta3)/4d0

         case(8)            !-110

            zp = 1d0-x

            ris = (pi**2*zp)/12d0 + ((pi**2)/48d0 - 1d0
     &/4d0)*zp**2 + ((pi**2)/144d0 - 1d0/8d0)*zp**3 + ((pi**2
     &)/384d0 - 35d0/576d0)*zp**4 + (-(11d0/360d0) + (pi**2)/
     &960d0)*zp**5 + ((pi**2)/2304d0 - 347d0/21600d0)*zp**6 +
     & (pi**2*ll2)/12d0 - zeta3

         case(9)            !-111

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*ll2)/12d0) + (ll2**3)/6d0 + 
     &zp**4*(-(1d0/1024d0) + (llzp)/256d0 - (llzp**2)/128d0) 
     &+ zp**2*(-(1d0/32d0) + (llzp)/16d0 - (llzp**2)/16d0) + 
     &zp**5*(-(1d0/4000d0) + (llzp)/800d0 - (llzp**2)/320d0) 
     &+ zp**3*(-(1d0/216d0) + (llzp)/72d0 - (llzp**2)/48d0) +
     & zp*(-(1d0/2d0) + (llzp)/2d0 - (llzp**2)/4d0) + zp**6*(
     &-(1d0/13824d0) + (llzp)/2304d0 - (llzp**2)/768d0) + (7d
     &0*zeta3)/8d0

         case(10)            !0-1-1

            zp = 1d0-x

            ris = -((zp*ll2**2)/2d0) + zp**5*(-(83d0/19
     &20d0) + (131d0*ll2)/960d0 - (ll2**2)/10d0) + zp**6*(-(1
     &1d0/288d0) + (661d0*ll2)/5760d0 - (ll2**2)/12d0) + zp**
     &2*((ll2)/4d0 - (ll2**2)/4d0) + zp**3*(-(1d0/24d0) + (5d
     &0*ll2)/24d0 - (ll2**2)/6d0) + zp**4*(-(3d0/64d0) + (ll2
     &)/6d0 - (ll2**2)/8d0) + (zeta3)/8d0

         case(11)            !0-10

            zp = 1d0-x

            ris = (pi**2*zp)/12d0 + (pi**2*zp**2)/24d0 
     &+ (-(1d0/12d0) + (pi**2)/36d0)*zp**3 + ((pi**2)/48d0 - 
     &5d0/48d0)*zp**4 + (-(5d0/48d0) + (pi**2)/60d0)*zp**5 + 
     &(-(47d0/480d0) + (pi**2)/72d0)*zp**6 - (3d0*zeta3)/2d0

         case(12)            !0-11

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*ll2)/4d0) + zp*(-((pi**2)/12
     &d0) + (ll2**2)/2d0) + zp**2*(-((pi**2)/24d0) + 3d0/8d0 
     &+ (ll2**2)/4d0 - (llzp)/4d0) + zp**3*(-((pi**2)/36d0) +
     & 37d0/144d0 + (ll2**2)/6d0 - (5d0*llzp)/24d0) + zp**6*(
     &13369d0/115200d0 - (pi**2)/72d0 + (ll2**2)/12d0 - (661d
     &0*llzp)/5760d0) + zp**4*(-((pi**2)/48d0) + 107d0/576d0 
     &+ (ll2**2)/8d0 - (llzp)/6d0) + zp**5*(-((pi**2)/60d0) +
     & 8257d0/57600d0 + (ll2**2)/10d0 - (131d0*llzp)/960d0) +
     & (13d0*zeta3)/8d0

         case(13)            !00-1

            zp = 1d0-x

            ris = -((pi**2*zp)/12d0) + zp**4*(-((pi**2)
     &/48d0) - 11d0/96d0 + (11d0*ll2)/24d0) + zp**2*(-((pi**2
     &)/24d0) + (ll2)/2d0) + zp**3*(-(1d0/12d0) - (pi**2)/36d
     &0 + (ll2)/2d0) + zp**6*(-((pi**2)/72d0) - 731d0/5760d0 
     &+ (137d0*ll2)/360d0) + zp**5*(-((pi**2)/60d0) - 1d0/8d0
     & + (5d0*ll2)/12d0) + (3d0*zeta3)/4d0

         case(14)            !000

            zp = 1d0-x

            ris = -((zp**3)/6d0) - (zp**4)/4d0 - (7d0*z
     &p**5)/24d0 - (5d0*zp**6)/16d0

         case(15)            !001

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*zp)/6d0) + zp**4*(-((pi**2)/
     &24d0) + 131d0/288d0 - (11d0*llzp)/24d0) + zp**2*(-((pi*
     &*2)/12d0) + 3d0/4d0 - (llzp)/2d0) + zp**3*(-((pi**2)/18
     &d0) + 7d0/12d0 - (llzp)/2d0) + zp**6*(-((pi**2)/36d0) +
     & 2213d0/7200d0 - (137d0*llzp)/360d0) + zp**5*(-((pi**2)
     &/30d0) + 53d0/144d0 - (5d0*llzp)/12d0) + zeta3

         case(16)            !01-1

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**2*ll2)/4d0 + zp*((pi**2)/12d0 - 
     &ll2 - (ll2**2)/2d0 + ll2*llzp) + zp**2*((pi**2)/24d0 - 
     &1d0/4d0 - (ll2)/4d0 - (ll2**2)/4d0 + (ll2*llzp)/2d0) + 
     &zp**3*((pi**2)/36d0 - 3d0/16d0 - (ll2)/9d0 - (ll2**2)/6
     &d0 + (ll2*llzp)/3d0) + zp**4*((pi**2)/48d0 - 83d0/576d0
     & - (ll2)/16d0 - (ll2**2)/8d0 + (ll2*llzp)/4d0) + zp**5*
     &(-(1337d0/11520d0) + (pi**2)/60d0 - (ll2)/25d0 - (ll2**
     &2)/10d0 + (ll2*llzp)/5d0) + zp**6*(-(33497d0/345600d0) 
     &+ (pi**2)/72d0 - (ll2)/36d0 - (ll2**2)/12d0 + (ll2*llzp
     &)/6d0) - zeta3

         case(17)            !010

            zp = 1d0-x

            ris = (pi**2*zp)/6d0 + ((pi**2)/12d0 - 1d0/
     &2d0)*zp**2 + ((pi**2)/18d0 - 5d0/12d0)*zp**3 + ((pi**2)
     &/24d0 - 49d0/144d0)*zp**4 + ((pi**2)/30d0 - 41d0/144d0)
     &*zp**5 + ((pi**2)/36d0 - 5269d0/21600d0)*zp**6 - 2*zeta
     &3

         case(18)            !011

            zp = 1d0-x
            llzp = log(zp)

            ris = zp**5*(-(1d0/125d0) + (llzp)/25d0 - (
     &llzp**2)/10d0) + zp**6*(-(1d0/216d0) + (llzp)/36d0 - (l
     &lzp**2)/12d0) + zp*(-1 + llzp - (llzp**2)/2d0) + zp**2*
     &(-(1d0/8d0) + (llzp)/4d0 - (llzp**2)/4d0) + zp**3*(-(1d
     &0/27d0) + (llzp)/9d0 - (llzp**2)/6d0) + zp**4*(-(1d0/64
     &d0) + (llzp)/16d0 - (llzp**2)/8d0) + zeta3

         case(19)            !1-1-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((pi**2*ll2)/12d0) + (zp*ll2)/2d0 + 
     &(ll2**3)/3d0 + zp**2*(-(1d0/16d0) + (ll2)/16d0) + zp**6
     &*(-(137d0/138240d0) + (ll2)/2304d0) + zp**4*(-(11d0/153
     &6d0) + (ll2)/256d0) + zp**3*(-(1d0/48d0) + (ll2)/72d0) 
     &+ zp**5*(-(1d0/384d0) + (ll2)/800d0) - (ll2**2*llzp)/2d
     &0 + (zeta3)/8d0

         case(20)            !1-10

            zp = 1d0-x
            llzp = log(zp)

            ris = -((zp**2)/8d0) - (zp**3)/18d0 - (5d0*
     &zp**4)/192d0 - (zp**5)/75d0 - (zp**6)/135d0 - (pi**2*ll
     &2)/4d0 + (pi**2*llzp)/12d0 + (13d0*zeta3)/8d0

         case(21)            !1-11

            zp = 1d0-x
            llzp = log(zp)

            ris = (pi**2*ll2)/6d0 - (ll2**3)/3d0 - (pi*
     &*2*llzp)/12d0 + (ll2**2*llzp)/2d0 + zp**2*(1d0/16d0 - (
     &llzp)/16d0) + zp**6*(1d0/6912d0 - (llzp)/2304d0) + zp**
     &4*(1d0/512d0 - (llzp)/256d0) + zp*(1 - (llzp)/2d0) + zp
     &**3*(1d0/108d0 - (llzp)/72d0) + zp**5*(1d0/2000d0 - (ll
     &zp)/800d0) - (7d0*zeta3)/4d0

         case(22)            !10-1

            zp = 1d0-x
            llzp = log(zp)

            ris = zp*ll2 + zp**4*(-(1d0/24d0) + (ll2)/1
     &6d0) + zp**5*(-(131d0/4800d0) + (ll2)/25d0) + zp**6*(-(
     &661d0/34560d0) + (ll2)/36d0) + zp**2*(-(1d0/8d0) + (ll2
     &)/4d0) + zp**3*(-(5d0/72d0) + (ll2)/9d0) - (pi**2*llzp)
     &/12d0 - (5d0*zeta3)/8d0

         case(23)            !100

            zp = 1d0-x

            ris = -((zp**2)/4d0) - (zp**3)/6d0 - (11d0*
     &zp**4)/96d0 - (zp**5)/12d0 - (137d0*zp**6)/2160d0 + zet
     &a3

         case(24)            !101

            zp = 1d0-x
            llzp = log(zp)

            ris = zp*(2 - llzp) - (pi**2*llzp)/6d0 + zp
     &**4*(1d0/32d0 - (llzp)/16d0) + zp**5*(2d0/125d0 - (llzp
     &)/25d0) + zp**6*(1d0/108d0 - (llzp)/36d0) + zp**2*(1d0/
     &4d0 - (llzp)/4d0) + zp**3*(2d0/27d0 - (llzp)/9d0) - 2*z
     &eta3

         case(25)            !11-1

            zp = 1d0-x
            llzp = log(zp)

            ris = -((zp)/2d0) - (zp**2)/32d0 - (zp**3)/
     &216d0 - (zp**4)/1024d0 - (zp**5)/4000d0 - (zp**6)/13824
     &d0 - (pi**2*ll2)/12d0 + (ll2**3)/6d0 + (pi**2*llzp)/12d
     &0 - (ll2**2*llzp)/2d0 + (ll2*llzp**2)/2d0 + (7d0*zeta3)
     &/8d0

         case(26)            !110

            zp = 1d0-x
            llzp = log(zp)

            ris = -zp - (zp**2)/8d0 - (zp**3)/27d0 - (z
     &p**4)/64d0 - (zp**5)/125d0 - (zp**6)/216d0 + (pi**2*llz
     &p)/6d0 + zeta3

         case(27)            !111

            zp = 1d0-x
            llzp = log(zp)

            ris = -((llzp**3)/6d0)
c End of expansions around x = +1

      end select
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)
         if (n3.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n3.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1.and.n2.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n3.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1.and.n2.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
         endif
      endif    

      HPL3ar1=ris
      return
      end function
