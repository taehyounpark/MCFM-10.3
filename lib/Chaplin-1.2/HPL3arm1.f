c ------------------------------------------------
      double complex function HPL3arm1(n1,n2,n3,x)
      implicit none
      integer n1,n2,n3,j,bcflag,s,szp
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

      select case (j)

c This was file contains the Taylor 
c expansions around x = -1
c The expansion parameter is zp = x+1

         case(1)            !-1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (llzp**3)/6d0

         case(2)            !-1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -zp - (zp**2)/8d0 - (zp**3)/27d0 - (z
     &p**4)/64d0 - (zp**5)/125d0 - (zp**6)/216d0 - (zp**7)/34
     &3d0 - (zp**8)/512d0 - (zp**9)/729d0 + (pi**2*llzp)/6d0 
     &+ (myi*pi*szp*llzp**2)/2d0 + zeta3

         case(3)            !-1-11

            zp = x+1d0
            llzp = log(zp)

            ris = (zp)/2d0 + (zp**2)/32d0 + (zp**3)/216
     &d0 + (zp**4)/1024d0 + (zp**5)/4000d0 + (zp**6)/13824d0 
     &+ (zp**7)/43904d0 + (zp**8)/131072d0 + (zp**9)/373248d0
     & + (pi**2*ll2)/12d0 - (ll2**3)/6d0 - (pi**2*llzp)/12d0 
     &+ (ll2**2*llzp)/2d0 - (ll2*llzp**2)/2d0 - (7d0*zeta3)/8
     &d0

         case(4)            !-10-1

            zp = x+1d0
            llzp = log(zp)

            ris = zp*(2 - llzp) - (pi**2*llzp)/6d0 + zp
     &**4*(1d0/32d0 - (llzp)/16d0) + zp**5*(2d0/125d0 - (llzp
     &)/25d0) + zp**6*(1d0/108d0 - (llzp)/36d0) + zp**7*(2d0/
     &343d0 - (llzp)/49d0) + zp**2*(1d0/4d0 - (llzp)/4d0) + z
     &p**8*(1d0/256d0 - (llzp)/64d0) + zp**9*(2d0/729d0 - (ll
     &zp)/81d0) + zp**3*(2d0/27d0 - (llzp)/9d0) - 2*zeta3

         case(5)            !-100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (myi*pi**3*szp)/6d0 - myi*pi*szp*zp +
     & (1d0/4d0 - (myi*pi*szp)/4d0)*zp**2 + (1d0/6d0 - (myi*p
     &i*szp)/9d0)*zp**3 + (11d0/96d0 - (myi*pi*szp)/16d0)*zp*
     &*4 + (1d0/12d0 - (myi*pi*szp)/25d0)*zp**5 + (137d0/2160
     &d0 - (myi*pi*szp)/36d0)*zp**6 + (1d0/20d0 - (myi*pi*szp
     &)/49d0)*zp**7 + (363d0/8960d0 - (myi*pi*szp)/64d0)*zp**
     &8 + (761d0/22680d0 - (myi*pi*szp)/81d0)*zp**9 - (pi**2*
     &llzp)/2d0 - zeta3

         case(6)            !-101

            zp = x+1d0
            llzp = log(zp)

            ris = zp*ll2 + zp**4*(-(1d0/24d0) + (ll2)/1
     &6d0) + zp**5*(-(131d0/4800d0) + (ll2)/25d0) + zp**6*(-(
     &661d0/34560d0) + (ll2)/36d0) + zp**7*(-(1327d0/94080d0)
     & + (ll2)/49d0) + zp**2*(-(1d0/8d0) + (ll2)/4d0) + zp**8
     &*(-(1163d0/107520d0) + (ll2)/64d0) + zp**9*(-(148969d0/
     &17418240d0) + (ll2)/81d0) + zp**3*(-(5d0/72d0) + (ll2)/
     &9d0) - (pi**2*llzp)/12d0 - (5d0*zeta3)/8d0

         case(7)            !-11-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*ll2)/6d0) + (ll2**3)/3d0 + (
     &pi**2*llzp)/12d0 - (ll2**2*llzp)/2d0 + zp**8*(-(1d0/655
     &36d0) + (llzp)/16384d0) + zp**2*(-(1d0/16d0) + (llzp)/1
     &6d0) + zp**6*(-(1d0/6912d0) + (llzp)/2304d0) + zp**4*(-
     &(1d0/512d0) + (llzp)/256d0) + zp*(-1 + (llzp)/2d0) + zp
     &**9*(-(1d0/186624d0) + (llzp)/41472d0) + zp**7*(-(1d0/2
     &1952d0) + (llzp)/6272d0) + zp**3*(-(1d0/108d0) + (llzp)
     &/72d0) + zp**5*(-(1d0/2000d0) + (llzp)/800d0) + (7d0*ze
     &ta3)/4d0

         case(8)            !-110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((myi*pi**3*szp)/12d0) + (myi*pi*szp
     &*zp)/2d0 + (-(1d0/8d0) + (myi*pi*szp)/16d0)*zp**2 + (-(
     &1d0/18d0) + (myi*pi*szp)/72d0)*zp**3 + (-(5d0/192d0) + 
     &(myi*pi*szp)/256d0)*zp**4 + (-(1d0/75d0) + (myi*pi*szp)
     &/800d0)*zp**5 + (-(1d0/135d0) + (myi*pi*szp)/2304d0)*zp
     &**6 + (-(13d0/2940d0) + (myi*pi*szp)/6272d0)*zp**7 + (-
     &(151d0/53760d0) + (myi*pi*szp)/16384d0)*zp**8 + (-(16d0
     &/8505d0) + (myi*pi*szp)/41472d0)*zp**9 - (pi**2*ll2)/4d
     &0 + (myi*pi*szp*ll2**2)/2d0 + (pi**2*llzp)/12d0 - myi*p
     &i*szp*ll2*llzp + (13d0*zeta3)/8d0

         case(9)            !-111

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**2*ll2)/12d0 - (zp*ll2)/2d0 - (ll
     &2**3)/3d0 + zp**8*(363d0/2293760d0 - (ll2)/16384d0) + z
     &p**2*(1d0/16d0 - (ll2)/16d0) + zp**6*(137d0/138240d0 - 
     &(ll2)/2304d0) + zp**4*(11d0/1536d0 - (ll2)/256d0) + zp*
     &*9*(761d0/11612160d0 - (ll2)/41472d0) + zp**7*(1d0/2560
     &d0 - (ll2)/6272d0) + zp**3*(1d0/48d0 - (ll2)/72d0) + zp
     &**5*(1d0/384d0 - (ll2)/800d0) + (ll2**2*llzp)/2d0 - (ze
     &ta3)/8d0

         case(10)            !0-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = zp**5*(-(1d0/125d0) + (llzp)/25d0 - (
     &llzp**2)/10d0) + zp**6*(-(1d0/216d0) + (llzp)/36d0 - (l
     &lzp**2)/12d0) + zp**7*(-(1d0/343d0) + (llzp)/49d0 - (ll
     &zp**2)/14d0) + zp**8*(-(1d0/512d0) + (llzp)/64d0 - (llz
     &p**2)/16d0) + zp**9*(-(1d0/729d0) + (llzp)/81d0 - (llzp
     &**2)/18d0) + zp*(-1 + llzp - (llzp**2)/2d0) + zp**2*(-(
     &1d0/8d0) + (llzp)/4d0 - (llzp**2)/4d0) + zp**3*(-(1d0/2
     &7d0) + (llzp)/9d0 - (llzp**2)/6d0) + zp**4*(-(1d0/64d0)
     & + (llzp)/16d0 - (llzp**2)/8d0) + zeta3

         case(11)            !0-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((myi*pi**3*szp)/6d0) + zp*(-((pi**2
     &)/6d0) + myi*pi*szp - myi*pi*szp*llzp) + zp**2*(-((pi**
     &2)/12d0) + 1d0/2d0 + (myi*pi*szp)/4d0 - (myi*pi*szp*llz
     &p)/2d0) + zp**3*(-((pi**2)/18d0) + 5d0/12d0 + (myi*pi*s
     &zp)/9d0 - (myi*pi*szp*llzp)/3d0) + zp**4*(-((pi**2)/24d
     &0) + 49d0/144d0 + (myi*pi*szp)/16d0 - (myi*pi*szp*llzp)
     &/4d0) + zp**5*(-((pi**2)/30d0) + 41d0/144d0 + (myi*pi*s
     &zp)/25d0 - (myi*pi*szp*llzp)/5d0) + zp**6*(-((pi**2)/36
     &d0) + 5269d0/21600d0 + (myi*pi*szp)/36d0 - (myi*pi*szp*
     &llzp)/6d0) + zp**7*(-((pi**2)/42d0) + 767d0/3600d0 + (m
     &yi*pi*szp)/49d0 - (myi*pi*szp*llzp)/7d0) + zp**8*(26668
     &1d0/1411200d0 - (pi**2)/48d0 + (myi*pi*szp)/64d0 - (myi
     &*pi*szp*llzp)/8d0) + zp**9*(-((pi**2)/54d0) + 1077749d0
     &/6350400d0 + (myi*pi*szp)/81d0 - (myi*pi*szp*llzp)/9d0)
     & + 2*zeta3

         case(12)            !0-11

            zp = x+1d0
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
     &)/6d0) + zp**7*(-(5587d0/67200d0) + (pi**2)/84d0 - (ll2
     &)/49d0 - (ll2**2)/14d0 + (ll2*llzp)/7d0) + zp**8*(-(136
     &919d0/1881600d0) + (pi**2)/96d0 - (ll2)/64d0 - (ll2**2)
     &/16d0 + (ll2*llzp)/8d0) + zp**9*((pi**2)/108d0 - 350549
     &39d0/541900800d0 - (ll2)/81d0 - (ll2**2)/18d0 + (ll2*ll
     &zp)/9d0) - zeta3

         case(13)            !00-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**2*zp)/6d0 + zp**4*((pi**2)/24d0 
     &- 131d0/288d0 + (11d0*llzp)/24d0) + zp**2*((pi**2)/12d0
     & - 3d0/4d0 + (llzp)/2d0) + zp**3*((pi**2)/18d0 - 7d0/12
     &d0 + (llzp)/2d0) + zp**6*((pi**2)/36d0 - 2213d0/7200d0 
     &+ (137d0*llzp)/360d0) + zp**8*((pi**2)/48d0 - 647707d0/
     &2822400d0 + (363d0*llzp)/1120d0) + zp**5*((pi**2)/30d0 
     &- 53d0/144d0 + (5d0*llzp)/12d0) + zp**9*((pi**2)/54d0 -
     & 1290829d0/6350400d0 + (761d0*llzp)/2520d0) + zp**7*((p
     &i**2)/42d0 - 947d0/3600d0 + (7d0*llzp)/20d0) - zeta3

         case(14)            !000

            zp = x+1d0
            szp = s(zp)

            ris = -((myi*pi**3*szp)/6d0) + (pi**2*zp)/2
     &d0 + ((pi**2)/4d0 + (myi*pi*szp)/2d0)*zp**2 + (-(1d0/6d
     &0) + (pi**2)/6d0 + (myi*pi*szp)/2d0)*zp**3 + (-(1d0/4d0
     &) + (pi**2)/8d0 + (myi*pi*11d0*szp)/24d0)*zp**4 + ((pi*
     &*2)/10d0 - 7d0/24d0 + (myi*pi*5d0*szp)/12d0)*zp**5 + ((
     &pi**2)/12d0 - 5d0/16d0 + (myi*pi*137d0*szp)/360d0)*zp**
     &6 + ((pi**2)/14d0 - 29d0/90d0 + (myi*pi*7d0*szp)/20d0)*
     &zp**7 + ((pi**2)/16d0 - 469d0/1440d0 + (myi*pi*363d0*sz
     &p)/1120d0)*zp**8 + ((pi**2)/18d0 - 29531d0/90720d0 + (m
     &yi*pi*761d0*szp)/2520d0)*zp**9

         case(15)            !001

            zp = x+1d0

            ris = (pi**2*zp)/12d0 + zp**4*((pi**2)/48d0
     & + 11d0/96d0 - (11d0*ll2)/24d0) + zp**2*((pi**2)/24d0 -
     & (ll2)/2d0) + zp**3*(1d0/12d0 + (pi**2)/36d0 - (ll2)/2d
     &0) + zp**6*((pi**2)/72d0 + 731d0/5760d0 - (137d0*ll2)/3
     &60d0) + zp**8*(3931d0/32256d0 + (pi**2)/96d0 - (363d0*l
     &l2)/1120d0) + zp**5*((pi**2)/60d0 + 1d0/8d0 - (5d0*ll2)
     &/12d0) + zp**9*((pi**2)/108d0 + 42799d0/362880d0 - (761
     &d0*ll2)/2520d0) + zp**7*(721d0/5760d0 + (pi**2)/84d0 - 
     &(7d0*ll2)/20d0) - (3d0*zeta3)/4d0

         case(16)            !01-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*ll2)/4d0) + zp*(-((pi**2)/12
     &d0) + (ll2**2)/2d0) + zp**8*(314543d0/3763200d0 - (pi**
     &2)/96d0 + (ll2**2)/16d0 - (1163d0*llzp)/13440d0) + zp**
     &7*(-((pi**2)/84d0) + 953d0/9800d0 + (ll2**2)/14d0 - (13
     &27d0*llzp)/13440d0) + zp**9*(-((pi**2)/108d0) + 3572057
     &71d0/4877107200d0 + (ll2**2)/18d0 - (148969d0*llzp)/193
     &5360d0) + zp**2*(-((pi**2)/24d0) + 3d0/8d0 + (ll2**2)/4
     &d0 - (llzp)/4d0) + zp**3*(-((pi**2)/36d0) + 37d0/144d0 
     &+ (ll2**2)/6d0 - (5d0*llzp)/24d0) + zp**6*(13369d0/1152
     &00d0 - (pi**2)/72d0 + (ll2**2)/12d0 - (661d0*llzp)/5760
     &d0) + zp**4*(-((pi**2)/48d0) + 107d0/576d0 + (ll2**2)/8
     &d0 - (llzp)/6d0) + zp**5*(-((pi**2)/60d0) + 8257d0/5760
     &0d0 + (ll2**2)/10d0 - (131d0*llzp)/960d0) + (13d0*zeta3
     &)/8d0

         case(17)            !010

            zp = x+1d0
            szp = s(zp)

            ris = -((myi*pi**3*szp)/12d0) + zp*(-((pi**
     &2)/12d0) + myi*pi*szp*ll2) + zp**2*(-((pi**2)/24d0) - (
     &myi*pi*szp)/4d0 + (myi*pi*szp*ll2)/2d0) + zp**3*(1d0/12
     &d0 - (pi**2)/36d0 - (myi*pi*5d0*szp)/24d0 + (myi*pi*szp
     &*ll2)/3d0) + zp**4*(-((pi**2)/48d0) + 5d0/48d0 - (myi*p
     &i*szp)/6d0 + (myi*pi*szp*ll2)/4d0) + zp**5*(5d0/48d0 - 
     &(pi**2)/60d0 - (myi*pi*131d0*szp)/960d0 + (myi*pi*szp*l
     &l2)/5d0) + zp**6*(47d0/480d0 - (pi**2)/72d0 - (myi*pi*6
     &61d0*szp)/5760d0 + (myi*pi*szp*ll2)/6d0) + zp**7*(13d0/
     &144d0 - (pi**2)/84d0 - (myi*pi*1327d0*szp)/13440d0 + (m
     &yi*pi*szp*ll2)/7d0) + zp**8*(3341d0/40320d0 - (pi**2)/9
     &6d0 - (myi*pi*1163d0*szp)/13440d0 + (myi*pi*szp*ll2)/8d
     &0) + zp**9*(13817d0/181440d0 - (pi**2)/108d0 - (myi*pi*
     &148969d0*szp)/1935360d0 + (myi*pi*szp*ll2)/9d0) + (3d0*
     &zeta3)/2d0

         case(18)            !011

            zp = x+1d0

            ris = -((zp*ll2**2)/2d0) + zp**5*(-(83d0/19
     &20d0) + (131d0*ll2)/960d0 - (ll2**2)/10d0) + zp**6*(-(1
     &1d0/288d0) + (661d0*ll2)/5760d0 - (ll2**2)/12d0) + zp**
     &7*(-(5417d0/161280d0) + (1327d0*ll2)/13440d0 - (ll2**2)
     &/14d0) + zp**8*(-(137d0/4608d0) + (1163d0*ll2)/13440d0 
     &- (ll2**2)/16d0) + zp**9*(-(617027d0/23224320d0) + (148
     &969d0*ll2)/1935360d0 - (ll2**2)/18d0) + zp**2*((ll2)/4d
     &0 - (ll2**2)/4d0) + zp**3*(-(1d0/24d0) + (5d0*ll2)/24d0
     & - (ll2**2)/6d0) + zp**4*(-(3d0/64d0) + (ll2)/6d0 - (ll
     &2**2)/8d0) + (zeta3)/8d0

         case(19)            !1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**2*ll2)/12d0 - (ll2**3)/6d0 + zp*
     &*4*(1d0/1024d0 - (llzp)/256d0 + (llzp**2)/128d0) + zp**
     &2*(1d0/32d0 - (llzp)/16d0 + (llzp**2)/16d0) + zp**7*(1d
     &0/43904d0 - (llzp)/6272d0 + (llzp**2)/1792d0) + zp**5*(
     &1d0/4000d0 - (llzp)/800d0 + (llzp**2)/320d0) + zp**8*(1
     &d0/131072d0 - (llzp)/16384d0 + (llzp**2)/4096d0) + zp**
     &3*(1d0/216d0 - (llzp)/72d0 + (llzp**2)/48d0) + zp*(1d0/
     &2d0 - (llzp)/2d0 + (llzp**2)/4d0) + zp**6*(1d0/13824d0 
     &- (llzp)/2304d0 + (llzp**2)/768d0) + zp**9*(1d0/373248d
     &0 - (llzp)/41472d0 + (llzp**2)/9216d0) - (7d0*zeta3)/8d
     &0

         case(20)            !1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (myi*pi**3*szp)/12d0 + (pi**2*ll2)/12
     &d0 - (myi*pi*szp*ll2**2)/2d0 + zp**5*(-(11d0/360d0) + (
     &pi**2)/960d0 - (myi*pi*szp)/800d0 + (myi*pi*szp*llzp)/1
     &60d0) + zp**8*((pi**2)/12288d0 - 9701d0/1881600d0 - (my
     &i*pi*szp)/16384d0 + (myi*pi*szp*llzp)/2048d0) + zp**3*(
     &(pi**2)/144d0 - 1d0/8d0 - (myi*pi*szp)/72d0 + (myi*pi*s
     &zp*llzp)/24d0) + zp*((pi**2)/12d0 - (myi*pi*szp)/2d0 + 
     &(myi*pi*szp*llzp)/2d0) + zp**6*((pi**2)/2304d0 - 347d0/
     &21600d0 - (myi*pi*szp)/2304d0 + (myi*pi*szp*llzp)/384d0
     &) + zp**9*((pi**2)/27648d0 - 209d0/66150d0 - (myi*pi*sz
     &p)/41472d0 + (myi*pi*szp*llzp)/4608d0) + zp**4*((pi**2)
     &/384d0 - 35d0/576d0 - (myi*pi*szp)/256d0 + (myi*pi*szp*
     &llzp)/64d0) + zp**7*(-(149d0/16800d0) + (pi**2)/5376d0 
     &- (myi*pi*szp)/6272d0 + (myi*pi*szp*llzp)/896d0) + zp**
     &2*((pi**2)/48d0 - 1d0/4d0 - (myi*pi*szp)/16d0 + (myi*pi
     &*szp*llzp)/8d0) - zeta3

         case(21)            !1-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*ll2)/12d0) + (ll2**3)/6d0 + 
     &zp**5*(-((pi**2)/1920d0) + 41d0/4608d0 + (ll2)/800d0 + 
     &(ll2**2)/320d0 - (ll2*llzp)/160d0) + zp**8*(-((pi**2)/2
     &4576d0) + 266681d0/361267200d0 + (ll2)/16384d0 + (ll2**
     &2)/4096d0 - (ll2*llzp)/2048d0) + zp**3*(-((pi**2)/288d0
     &) + 5d0/96d0 + (ll2)/72d0 + (ll2**2)/48d0 - (ll2*llzp)/
     &24d0) + zp*(-((pi**2)/24d0) + (ll2)/2d0 + (ll2**2)/4d0 
     &- (ll2*llzp)/2d0) + zp**6*(-((pi**2)/4608d0) + 5269d0/1
     &382400d0 + (ll2)/2304d0 + (ll2**2)/768d0 - (ll2*llzp)/3
     &84d0) + zp**9*(1077749d0/3251404800d0 - (pi**2)/55296d0
     & + (ll2)/41472d0 + (ll2**2)/9216d0 - (ll2*llzp)/4608d0)
     & + zp**4*(49d0/2304d0 - (pi**2)/768d0 + (ll2)/256d0 + (
     &ll2**2)/128d0 - (ll2*llzp)/64d0) + zp**7*(-((pi**2)/107
     &52d0) + 767d0/460800d0 + (ll2)/6272d0 + (ll2**2)/1792d0
     & - (ll2*llzp)/896d0) + zp**2*(1d0/8d0 - (pi**2)/96d0 + 
     &(ll2)/16d0 + (ll2**2)/16d0 - (ll2*llzp)/8d0) + (zeta3)/
     &4d0

         case(22)            !10-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*zp)/12d0) + (pi**2*ll2)/6d0 
     &+ zp**5*(79d0/1800d0 - (pi**2)/960d0 - (llzp)/15d0) + z
     &p**7*(521d0/39200d0 - (pi**2)/5376d0 - (13d0*llzp)/420d
     &0) + zp**6*(-((pi**2)/2304d0) + 169d0/7200d0 - (2d0*llz
     &p)/45d0) + zp**2*(-((pi**2)/48d0) + 3d0/8d0 - (llzp)/4d
     &0) + zp**4*(25d0/288d0 - (pi**2)/384d0 - (5d0*llzp)/48d
     &0) + zp**8*(-((pi**2)/12288d0) + 7493d0/940800d0 - (151
     &d0*llzp)/6720d0) + zp**3*(-((pi**2)/144d0) + 13d0/72d0 
     &- (llzp)/6d0) + zp**9*(-((pi**2)/27648d0) + 3001d0/5953
     &50d0 - (16d0*llzp)/945d0) - (5d0*zeta3)/8d0

         case(23)            !100

            zp = x+1d0
            szp = s(zp)

            ris = (myi*pi**3*szp)/12d0 - (pi**2*zp)/4d0
     & + (-((pi**2)/16d0) - (myi*pi*szp)/4d0)*zp**2 + (1d0/12
     &d0 - (pi**2)/48d0 - (myi*pi*szp)/6d0)*zp**3 + (-((pi**2
     &)/128d0) + 3d0/32d0 - (myi*pi*5d0*szp)/48d0)*zp**4 + (1
     &d0/12d0 - (pi**2)/320d0 - (myi*pi*szp)/15d0)*zp**5 + (5
     &d0/72d0 - (pi**2)/768d0 - (myi*pi*2d0*szp)/45d0)*zp**6 
     &+ (-((pi**2)/1792d0) + 41d0/720d0 - (myi*pi*13d0*szp)/4
     &20d0)*zp**7 + (-((pi**2)/4096d0) + 539d0/11520d0 - (myi
     &*pi*151d0*szp)/6720d0)*zp**8 + (22d0/567d0 - (pi**2)/92
     &16d0 - (myi*pi*16d0*szp)/945d0)*zp**9 + (pi**2*ll2)/2d0
     & - (3d0*zeta3)/4d0

         case(24)            !101

            zp = x+1d0

            ris = -((pi**2*zp)/24d0) + (pi**2*ll2)/12d0
     & + zp**5*(-((pi**2)/1920d0) - 1d0/30d0 + (ll2)/15d0) + 
     &zp**7*(-((pi**2)/10752d0) - 767d0/40320d0 + (13d0*ll2)/
     &420d0) + zp**6*(-((pi**2)/4608d0) - 97d0/3840d0 + (2d0*
     &ll2)/45d0) + zp**2*(-((pi**2)/96d0) + (ll2)/4d0) + zp**
     &4*(-(1d0/24d0) - (pi**2)/768d0 + (5d0*ll2)/48d0) + zp**
     &8*(-((pi**2)/24576d0) - 935d0/64512d0 + (151d0*ll2)/672
     &0d0) + zp**3*(-(1d0/24d0) - (pi**2)/288d0 + (ll2)/6d0) 
     &+ zp**9*(-(2041d0/181440d0) - (pi**2)/55296d0 + (16d0*l
     &l2)/945d0) - (zeta3)/4d0

         case(25)            !11-1

            zp = x+1d0
            llzp = log(zp)

            ris = (ll2**3)/6d0 + zp*((pi**2)/24d0 - (ll
     &2**2)/4d0) + zp**3*((pi**2)/288d0 - 7d0/96d0 - (ll2**2)
     &/48d0 + (llzp)/16d0) + zp**6*(-(2213d0/460800d0) + (pi*
     &*2)/4608d0 - (ll2**2)/768d0 + (137d0*llzp)/23040d0) + z
     &p**8*((pi**2)/24576d0 - 647707d0/722534400d0 - (ll2**2)
     &/4096d0 + (363d0*llzp)/286720d0) + zp**4*(-(131d0/4608d
     &0) + (pi**2)/768d0 - (ll2**2)/128d0 + (11d0*llzp)/384d0
     &) + zp**5*((pi**2)/1920d0 - 53d0/4608d0 - (ll2**2)/320d
     &0 + (5d0*llzp)/384d0) + zp**9*(-(1290829d0/3251404800d0
     &) + (pi**2)/55296d0 - (ll2**2)/9216d0 + (761d0*llzp)/12
     &90240d0) + zp**7*((pi**2)/10752d0 - 947d0/460800d0 - (l
     &l2**2)/1792d0 + (7d0*llzp)/2560d0) + zp**2*(-(3d0/16d0)
     & + (pi**2)/96d0 - (ll2**2)/16d0 + (llzp)/8d0) - (zeta3)
     &/8d0

         case(26)            !110

            zp = x+1d0
            szp = s(zp)

            ris = -((pi**2*ll2)/12d0) + (myi*pi*szp*ll2
     &**2)/2d0 + zp**5*((pi**2)/1920d0 - 1d0/40d0 + (myi*pi*5
     &d0*szp)/384d0 - (myi*pi*szp*ll2)/160d0) + zp**8*(-(1019
     &d0/161280d0) + (pi**2)/24576d0 + (myi*pi*363d0*szp)/286
     &720d0 - (myi*pi*szp*ll2)/2048d0) + zp**3*(-(1d0/24d0) +
     & (pi**2)/288d0 + (myi*pi*szp)/16d0 - (myi*pi*szp*ll2)/2
     &4d0) + zp*((pi**2)/24d0 - (myi*pi*szp*ll2)/2d0) + zp**6
     &*(-(23d0/1440d0) + (pi**2)/4608d0 + (myi*pi*137d0*szp)/
     &23040d0 - (myi*pi*szp*ll2)/384d0) + zp**9*((pi**2)/5529
     &6d0 - 23d0/5670d0 + (myi*pi*761d0*szp)/1290240d0 - (myi
     &*pi*szp*ll2)/4608d0) + zp**4*((pi**2)/768d0 - 7d0/192d0
     & + (myi*pi*11d0*szp)/384d0 - (myi*pi*szp*ll2)/64d0) + z
     &p**7*(-(101d0/10080d0) + (pi**2)/10752d0 + (myi*pi*7d0*
     &szp)/2560d0 - (myi*pi*szp*ll2)/896d0) + zp**2*((pi**2)/
     &96d0 + (myi*pi*szp)/8d0 - (myi*pi*szp*ll2)/8d0) + (zeta
     &3)/8d0

         case(27)            !111

            zp = x+1d0

            ris = (zp*ll2**2)/4d0 - (ll2**3)/6d0 + zp**
     &4*(1d0/64d0 - (11d0*ll2)/384d0 + (ll2**2)/128d0) + zp**
     &2*(-((ll2)/8d0) + (ll2**2)/16d0) + zp**7*(29d0/11520d0 
     &- (7d0*ll2)/2560d0 + (ll2**2)/1792d0) + zp**5*(7d0/768d
     &0 - (5d0*ll2)/384d0 + (ll2**2)/320d0) + zp**8*(469d0/36
     &8640d0 - (363d0*ll2)/286720d0 + (ll2**2)/4096d0) + zp**
     &3*(1d0/48d0 - (ll2)/16d0 + (ll2**2)/48d0) + zp**6*(5d0/
     &1024d0 - (137d0*ll2)/23040d0 + (ll2**2)/768d0) + zp**9*
     &(29531d0/46448640d0 - (761d0*ll2)/1290240d0 + (ll2**2)/
     &9216d0)
c End of expansions around x = -1

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

      HPL3arm1=ris
      return
      end function
