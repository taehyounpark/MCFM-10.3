      double complex function HPL4arm1(n1,n2,n3,n4,x)
      implicit none
      integer n1,n2,n3,n4,j,bcflag,s,szp
      double complex x,ris,myi,zp,llzp,cli4,cli4pt5
      double precision pi, zeta2, zeta3,ll2,xre

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      myi = dcmplx(0d0,1d0)
      ll2 = dlog(2d0)
      cli4pt5 = cli4(dcmplx(0.5d0,0d0))
      bcflag = 0

      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
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

         case(1)            !-1-1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (llzp**4)/24d0

         case(2)            !-1-1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/90d0 - zp - (zp**2)/16d0 - (z
     &p**3)/81d0 - (zp**4)/256d0 - (zp**5)/625d0 - (zp**6)/12
     &96d0 - (zp**7)/2401d0 - (zp**8)/4096d0 - (zp**9)/6561d0
     & + (pi**2*llzp**2)/12d0 + (myi*pi*szp*llzp**3)/6d0 + ll
     &zp*zeta3

         case(3)            !-1-1-11

            zp = x+1d0
            llzp = log(zp)

            ris = (zp)/2d0 + (zp**2)/64d0 + (zp**3)/648
     &d0 + (zp**4)/4096d0 + (zp**5)/20000d0 + (zp**6)/82944d0
     & + (zp**7)/307328d0 + (zp**8)/1048576d0 + (zp**9)/33592
     &32d0 + (pi**2*ll2*llzp)/12d0 - (ll2**3*llzp)/6d0 - (pi*
     &*2*llzp**2)/24d0 + (ll2**2*llzp**2)/4d0 - (ll2*llzp**3)
     &/6d0 - cli4pt5 - (7d0*llzp*zeta3)/8d0

         case(4)            !-1-10-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/30d0) + zp*(3 - llzp) - (pi
     &**2*llzp**2)/12d0 + zp**5*(3d0/625d0 - (llzp)/125d0) + 
     &zp**6*(1d0/432d0 - (llzp)/216d0) + zp**3*(1d0/27d0 - (l
     &lzp)/27d0) + zp**7*(3d0/2401d0 - (llzp)/343d0) + zp**8*
     &(3d0/4096d0 - (llzp)/512d0) + zp**4*(3d0/256d0 - (llzp)
     &/64d0) + zp**9*(1d0/2187d0 - (llzp)/729d0) + zp**2*(3d0
     &/16d0 - (llzp)/8d0) - 2*llzp*zeta3

         case(5)            !-1-100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4)/360d0) - myi*pi*szp*zp + (1
     &d0/8d0 - (myi*pi*szp)/8d0)*zp**2 + (1d0/18d0 - (myi*pi*
     &szp)/27d0)*zp**3 + (11d0/384d0 - (myi*pi*szp)/64d0)*zp*
     &*4 + (1d0/60d0 - (myi*pi*szp)/125d0)*zp**5 + (137d0/129
     &60d0 - (myi*pi*szp)/216d0)*zp**6 + (1d0/140d0 - (myi*pi
     &*szp)/343d0)*zp**7 + (363d0/71680d0 - (myi*pi*szp)/512d
     &0)*zp**8 + (761d0/204120d0 - (myi*pi*szp)/729d0)*zp**9 
     &+ (myi*pi**3*szp*llzp)/6d0 - (pi**2*llzp**2)/4d0 + myi*
     &pi*szp*zeta3 - llzp*zeta3

         case(6)            !-1-101

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/288d0 + zp*ll2 + (pi**2*ll2**
     &2)/24d0 - (ll2**4)/24d0 + zp**5*(-(131d0/24000d0) + (ll
     &2)/125d0) + zp**6*(-(661d0/207360d0) + (ll2)/216d0) + z
     &p**3*(-(5d0/216d0) + (ll2)/27d0) + zp**7*(-(1327d0/6585
     &60d0) + (ll2)/343d0) + zp**8*(-(1163d0/860160d0) + (ll2
     &)/512d0) + zp**4*(-(1d0/96d0) + (ll2)/64d0) + zp**9*(-(
     &148969d0/156764160d0) + (ll2)/729d0) + zp**2*(-(1d0/16d
     &0) + (ll2)/8d0) - (pi**2*llzp**2)/24d0 - cli4pt5 - (7d0
     &*ll2*zeta3)/8d0 - (5d0*llzp*zeta3)/8d0

         case(7)            !-1-11-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**2*ll2*llzp)/6d0) + (ll2**3*llz
     &p)/3d0 + (pi**2*llzp**2)/24d0 - (ll2**2*llzp**2)/4d0 + 
     &zp**4*(-(3d0/4096d0) + (llzp)/1024d0) + zp**8*(-(3d0/10
     &48576d0) + (llzp)/131072d0) + zp**6*(-(1d0/27648d0) + (
     &llzp)/13824d0) + zp**3*(-(1d0/216d0) + (llzp)/216d0) + 
     &zp*(-(3d0/2d0) + (llzp)/2d0) + zp**2*(-(3d0/64d0) + (ll
     &zp)/32d0) + zp**9*(-(1d0/1119744d0) + (llzp)/373248d0) 
     &+ zp**5*(-(3d0/20000d0) + (llzp)/4000d0) + zp**7*(-(3d0
     &/307328d0) + (llzp)/43904d0) + 3*cli4pt5 + (7d0*llzp*ze
     &ta3)/4d0

         case(8)            !-1-110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*11d0)/720d0) + (myi*pi*szp*z
     &p)/2d0 + (-(1d0/16d0) + (myi*pi*szp)/32d0)*zp**2 + (-(1
     &d0/54d0) + (myi*pi*szp)/216d0)*zp**3 + (-(5d0/768d0) + 
     &(myi*pi*szp)/1024d0)*zp**4 + (-(1d0/375d0) + (myi*pi*sz
     &p)/4000d0)*zp**5 + (-(1d0/810d0) + (myi*pi*szp)/13824d0
     &)*zp**6 + (-(13d0/20580d0) + (myi*pi*szp)/43904d0)*zp**
     &7 + (-(151d0/430080d0) + (myi*pi*szp)/131072d0)*zp**8 +
     & (-(16d0/76545d0) + (myi*pi*szp)/373248d0)*zp**9 + (myi
     &*pi**3*szp*ll2)/12d0 - (myi*pi*szp*ll2**3)/6d0 + (ll2**
     &4)/8d0 - (myi*pi**3*szp*llzp)/12d0 - (pi**2*ll2*llzp)/4
     &d0 + (myi*pi*szp*ll2**2*llzp)/2d0 + (pi**2*llzp**2)/24d
     &0 - (myi*pi*szp*ll2*llzp**2)/2d0 + 3*cli4pt5 - (myi*pi*
     &7d0*szp*zeta3)/8d0 + (13d0*llzp*zeta3)/8d0

         case(9)            !-1-111

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/720d0) - (zp*ll2)/2d0 - (pi
     &**2*ll2**2)/12d0 + (ll2**4)/8d0 + zp**4*(11d0/6144d0 - 
     &(ll2)/1024d0) + zp**8*(363d0/18350080d0 - (ll2)/131072d
     &0) + zp**6*(137d0/829440d0 - (ll2)/13824d0) + zp**3*(1d
     &0/144d0 - (ll2)/216d0) + zp**2*(1d0/32d0 - (ll2)/32d0) 
     &+ zp**9*(761d0/104509440d0 - (ll2)/373248d0) + zp**5*(1
     &d0/1920d0 - (ll2)/4000d0) + zp**7*(1d0/17920d0 - (ll2)/
     &43904d0) + (pi**2*ll2*llzp)/12d0 - (ll2**3*llzp)/3d0 + 
     &(ll2**2*llzp**2)/4d0 + ll2*zeta3 - (llzp*zeta3)/8d0

         case(10)            !-10-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/30d0 + zp**8*(-(3d0/4096d0) +
     & (llzp)/256d0 - (llzp**2)/128d0) + zp**9*(-(1d0/2187d0)
     & + (2d0*llzp)/729d0 - (llzp**2)/162d0) + zp**3*(-(1d0/2
     &7d0) + (2d0*llzp)/27d0 - (llzp**2)/18d0) + zp*(-3 + 2*l
     &lzp - (llzp**2)/2d0) + zp**4*(-(3d0/256d0) + (llzp)/32d
     &0 - (llzp**2)/32d0) + zp**5*(-(3d0/625d0) + (2d0*llzp)/
     &125d0 - (llzp**2)/50d0) + zp**6*(-(1d0/432d0) + (llzp)/
     &108d0 - (llzp**2)/72d0) + zp**2*(-(3d0/16d0) + (llzp)/4
     &d0 - (llzp**2)/8d0) + zp**7*(-(3d0/2401d0) + (2d0*llzp)
     &/343d0 - (llzp**2)/98d0) + llzp*zeta3

         case(11)            !-10-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*7d0)/360d0 - (myi*pi**3*szp*ll
     &zp)/6d0 + zp*(-((pi**2)/6d0) + 2*myi*pi*szp - myi*pi*sz
     &p*llzp) + zp**4*(49d0/576d0 - (pi**2)/96d0 + (myi*pi*sz
     &p)/32d0 - (myi*pi*szp*llzp)/16d0) + zp**5*(-((pi**2)/15
     &0d0) + 41d0/720d0 + (myi*pi*2d0*szp)/125d0 - (myi*pi*sz
     &p*llzp)/25d0) + zp**6*(-((pi**2)/216d0) + 5269d0/129600
     &d0 + (myi*pi*szp)/108d0 - (myi*pi*szp*llzp)/36d0) + zp*
     &*7*(-((pi**2)/294d0) + 767d0/25200d0 + (myi*pi*2d0*szp)
     &/343d0 - (myi*pi*szp*llzp)/49d0) + zp**2*(-((pi**2)/24d
     &0) + 1d0/4d0 + (myi*pi*szp)/4d0 - (myi*pi*szp*llzp)/4d0
     &) + zp**8*(266681d0/11289600d0 - (pi**2)/384d0 + (myi*p
     &i*szp)/256d0 - (myi*pi*szp*llzp)/64d0) + zp**9*(-((pi**
     &2)/486d0) + 1077749d0/57153600d0 + (myi*pi*2d0*szp)/729
     &d0 - (myi*pi*szp*llzp)/81d0) + zp**3*(-((pi**2)/54d0) +
     & 5d0/36d0 + (myi*pi*2d0*szp)/27d0 - (myi*pi*szp*llzp)/9
     &d0) - 2*myi*pi*szp*zeta3 + 2*llzp*zeta3

         case(12)            !-10-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/1440d0) + (pi**2*ll2**2)/24
     &d0 - (ll2**4)/24d0 + (pi**2*ll2*llzp)/4d0 + zp*((pi**2)
     &/12d0 - 2*ll2 - (ll2**2)/2d0 + ll2*llzp) + zp**4*((pi**
     &2)/192d0 - 83d0/2304d0 - (ll2)/32d0 - (ll2**2)/32d0 + (
     &ll2*llzp)/16d0) + zp**5*((pi**2)/300d0 - 1337d0/57600d0
     & - (2d0*ll2)/125d0 - (ll2**2)/50d0 + (ll2*llzp)/25d0) +
     & zp**6*(-(33497d0/2073600d0) + (pi**2)/432d0 - (ll2)/10
     &8d0 - (ll2**2)/72d0 + (ll2*llzp)/36d0) + zp**7*(-(5587d
     &0/470400d0) + (pi**2)/588d0 - (2d0*ll2)/343d0 - (ll2**2
     &)/98d0 + (ll2*llzp)/49d0) + zp**2*((pi**2)/48d0 - 1d0/8
     &d0 - (ll2)/4d0 - (ll2**2)/8d0 + (ll2*llzp)/4d0) + zp**8
     &*(-(136919d0/15052800d0) + (pi**2)/768d0 - (ll2)/256d0 
     &- (ll2**2)/128d0 + (ll2*llzp)/64d0) + zp**9*(-(35054939
     &d0/4877107200d0) + (pi**2)/972d0 - (2d0*ll2)/729d0 - (l
     &l2**2)/162d0 + (ll2*llzp)/81d0) + zp**3*((pi**2)/108d0 
     &- 1d0/16d0 - (2d0*ll2)/27d0 - (ll2**2)/18d0 + (ll2*llzp
     &)/9d0) - cli4pt5 + (7d0*ll2*zeta3)/4d0 - llzp*zeta3

         case(13)            !-100-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/72d0) + (pi**2*zp)/6d0 + zp
     &**5*(-(13d0/144d0) + (pi**2)/150d0 + (llzp)/12d0) + zp*
     &*7*((pi**2)/294d0 - 161d0/3600d0 + (llzp)/20d0) + zp**6
     &*((pi**2)/216d0 - 8009d0/129600d0 + (137d0*llzp)/2160d0
     &) + zp**2*((pi**2)/24d0 - 1d0/2d0 + (llzp)/4d0) + zp**3
     &*(-(1d0/4d0) + (pi**2)/54d0 + (llzp)/6d0) + zp**9*((pi*
     &*2)/486d0 - 167101d0/6350400d0 + (761d0*llzp)/22680d0) 
     &+ zp**8*((pi**2)/384d0 - 190513d0/5644800d0 + (363d0*ll
     &zp)/8960d0) + zp**4*(-(41d0/288d0) + (pi**2)/96d0 + (11
     &d0*llzp)/96d0) - llzp*zeta3

         case(14)            !-1000

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*13d0)/180d0) + (pi**2*zp)/2d
     &0 + ((pi**2)/8d0 + (myi*pi*szp)/4d0)*zp**2 + (-(1d0/18d
     &0) + (pi**2)/18d0 + (myi*pi*szp)/6d0)*zp**3 + (-(1d0/16
     &d0) + (pi**2)/32d0 + (myi*pi*11d0*szp)/96d0)*zp**4 + ((
     &pi**2)/50d0 - 7d0/120d0 + (myi*pi*szp)/12d0)*zp**5 + ((
     &pi**2)/72d0 - 5d0/96d0 + (myi*pi*137d0*szp)/2160d0)*zp*
     &*6 + (-(29d0/630d0) + (pi**2)/98d0 + (myi*pi*szp)/20d0)
     &*zp**7 + ((pi**2)/128d0 - 469d0/11520d0 + (myi*pi*363d0
     &*szp)/8960d0)*zp**8 + ((pi**2)/162d0 - 29531d0/816480d0
     & + (myi*pi*761d0*szp)/22680d0)*zp**9 - (myi*pi**3*szp*l
     &lzp)/6d0 - myi*pi*szp*zeta3

         case(15)            !-1001

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*11d0)/360d0) + (pi**2*zp)/12
     &d0 - (pi**2*ll2**2)/12d0 + (ll2**4)/12d0 + zp**5*((pi**
     &2)/300d0 + 1d0/40d0 - (ll2)/12d0) + zp**7*(103d0/5760d0
     & + (pi**2)/588d0 - (ll2)/20d0) + zp**6*((pi**2)/432d0 +
     & 731d0/34560d0 - (137d0*ll2)/2160d0) + zp**2*((pi**2)/4
     &8d0 - (ll2)/4d0) + zp**3*((pi**2)/108d0 + 1d0/36d0 - (l
     &l2)/6d0) + zp**9*(42799d0/3265920d0 + (pi**2)/972d0 - (
     &761d0*ll2)/22680d0) + zp**8*(3931d0/258048d0 + (pi**2)/
     &768d0 - (363d0*ll2)/8960d0) + zp**4*((pi**2)/192d0 + 11
     &d0/384d0 - (11d0*ll2)/96d0) + 2*cli4pt5 + (7d0*ll2*zeta
     &3)/4d0 - (3d0*llzp*zeta3)/4d0

         case(16)            !-101-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/160d0) - (pi**2*ll2**2)/8d0
     & + (ll2**4)/8d0 + zp*(-((pi**2)/12d0) + (ll2**2)/2d0) -
     & (pi**2*ll2*llzp)/4d0 + zp**8*(7401d0/627200d0 - (pi**2
     &)/768d0 + (ll2**2)/128d0 - (1163d0*llzp)/107520d0) + zp
     &**9*(398917091d0/43893964800d0 - (pi**2)/972d0 + (ll2**
     &2)/162d0 - (148969d0*llzp)/17418240d0) + zp**4*(-((pi**
     &2)/192d0) + 131d0/2304d0 + (ll2**2)/32d0 - (llzp)/24d0)
     & + zp**5*(-((pi**2)/300d0) + 9829d0/288000d0 + (ll2**2)
     &/50d0 - (131d0*llzp)/4800d0) + zp**6*(-((pi**2)/432d0) 
     &+ 46717d0/2073600d0 + (ll2**2)/72d0 - (661d0*llzp)/3456
     &0d0) + zp**3*(-((pi**2)/108d0) + 47d0/432d0 + (ll2**2)/
     &18d0 - (5d0*llzp)/72d0) + zp**2*(-((pi**2)/48d0) + 1d0/
     &4d0 + (ll2**2)/8d0 - (llzp)/8d0) + zp**7*(52379d0/32928
     &00d0 - (pi**2)/588d0 + (ll2**2)/98d0 - (1327d0*llzp)/94
     &080d0) + 3*cli4pt5 + (13d0*llzp*zeta3)/8d0

         case(17)            !-1010

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*17d0)/1440d0 + zp*(-((pi**2)/1
     &2d0) + myi*pi*szp*ll2) + zp**4*(-((pi**2)/192d0) + 5d0/
     &192d0 - (myi*pi*szp)/24d0 + (myi*pi*szp*ll2)/16d0) + zp
     &**5*(-((pi**2)/300d0) + 1d0/48d0 - (myi*pi*131d0*szp)/4
     &800d0 + (myi*pi*szp*ll2)/25d0) + zp**6*(-((pi**2)/432d0
     &) + 47d0/2880d0 - (myi*pi*661d0*szp)/34560d0 + (myi*pi*
     &szp*ll2)/36d0) + zp**7*(13d0/1008d0 - (pi**2)/588d0 - (
     &myi*pi*1327d0*szp)/94080d0 + (myi*pi*szp*ll2)/49d0) + z
     &p**2*(-((pi**2)/48d0) - (myi*pi*szp)/8d0 + (myi*pi*szp*
     &ll2)/4d0) + zp**8*(3341d0/322560d0 - (pi**2)/768d0 - (m
     &yi*pi*1163d0*szp)/107520d0 + (myi*pi*szp*ll2)/64d0) + z
     &p**9*(13817d0/1632960d0 - (pi**2)/972d0 - (myi*pi*14896
     &9d0*szp)/17418240d0 + (myi*pi*szp*ll2)/81d0) + zp**3*(-
     &((pi**2)/108d0) + 1d0/36d0 - (myi*pi*5d0*szp)/72d0 + (m
     &yi*pi*szp*ll2)/9d0) - (myi*pi**3*szp*llzp)/12d0 - (myi*
     &pi*5d0*szp*zeta3)/8d0 + (3d0*llzp*zeta3)/2d0

         case(18)            !-1011

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/480d0 - (zp*ll2**2)/2d0 + zp*
     &*8*(-(137d0/36864d0) + (1163d0*ll2)/107520d0 - (ll2**2)
     &/128d0) + zp**9*(-(617027d0/209018880d0) + (148969d0*ll
     &2)/17418240d0 - (ll2**2)/162d0) + zp**3*(-(1d0/72d0) + 
     &(5d0*ll2)/72d0 - (ll2**2)/18d0) + zp**4*(-(3d0/256d0) +
     & (ll2)/24d0 - (ll2**2)/32d0) + zp**5*(-(83d0/9600d0) + 
     &(131d0*ll2)/4800d0 - (ll2**2)/50d0) + zp**6*(-(11d0/172
     &8d0) + (661d0*ll2)/34560d0 - (ll2**2)/72d0) + zp**2*((l
     &l2)/8d0 - (ll2**2)/8d0) + zp**7*(-(5417d0/1128960d0) + 
     &(1327d0*ll2)/94080d0 - (ll2**2)/98d0) + (llzp*zeta3)/8d
     &0

         case(19)            !-11-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**2*ll2*llzp)/12d0 - (ll2**3*llzp)
     &/6d0 + zp**7*(3d0/307328d0 - (llzp)/21952d0 + (llzp**2)
     &/12544d0) + zp**3*(1d0/216d0 - (llzp)/108d0 + (llzp**2)
     &/144d0) + zp**5*(3d0/20000d0 - (llzp)/2000d0 + (llzp**2
     &)/1600d0) + zp**8*(3d0/1048576d0 - (llzp)/65536d0 + (ll
     &zp**2)/32768d0) + zp**2*(3d0/64d0 - (llzp)/16d0 + (llzp
     &**2)/32d0) + zp**6*(1d0/27648d0 - (llzp)/6912d0 + (llzp
     &**2)/4608d0) + zp*(3d0/2d0 - llzp + (llzp**2)/4d0) + zp
     &**4*(3d0/4096d0 - (llzp)/512d0 + (llzp**2)/512d0) + zp*
     &*9*(1d0/1119744d0 - (llzp)/186624d0 + (llzp**2)/82944d0
     &) - 3*cli4pt5 - (7d0*llzp*zeta3)/8d0

         case(20)            !-11-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/96d0 - (myi*pi**3*szp*ll2)/6d
     &0 - (pi**2*ll2**2)/24d0 + (myi*pi*szp*ll2**3)/3d0 - (ll
     &2**4)/8d0 + (myi*pi**3*szp*llzp)/12d0 + (pi**2*ll2*llzp
     &)/12d0 - (myi*pi*szp*ll2**2*llzp)/2d0 + zp**8*(-(9701d0
     &/15052800d0) + (pi**2)/98304d0 - (myi*pi*szp)/65536d0 +
     & (myi*pi*szp*llzp)/16384d0) + zp**2*(-(1d0/8d0) + (pi**
     &2)/96d0 - (myi*pi*szp)/16d0 + (myi*pi*szp*llzp)/16d0) +
     & zp**6*((pi**2)/13824d0 - 347d0/129600d0 - (myi*pi*szp)
     &/6912d0 + (myi*pi*szp*llzp)/2304d0) + zp**4*((pi**2)/15
     &36d0 - 35d0/2304d0 - (myi*pi*szp)/512d0 + (myi*pi*szp*l
     &lzp)/256d0) + zp*((pi**2)/12d0 - myi*pi*szp + (myi*pi*s
     &zp*llzp)/2d0) + zp**9*((pi**2)/248832d0 - 209d0/595350d
     &0 - (myi*pi*szp)/186624d0 + (myi*pi*szp*llzp)/41472d0) 
     &+ zp**7*(-(149d0/117600d0) + (pi**2)/37632d0 - (myi*pi*
     &szp)/21952d0 + (myi*pi*szp*llzp)/6272d0) + zp**3*(-(1d0
     &/24d0) + (pi**2)/432d0 - (myi*pi*szp)/108d0 + (myi*pi*s
     &zp*llzp)/72d0) + zp**5*(-(11d0/1800d0) + (pi**2)/4800d0
     & - (myi*pi*szp)/2000d0 + (myi*pi*szp*llzp)/800d0) - 3*c
     &li4pt5 + (myi*pi*7d0*szp*zeta3)/4d0 - llzp*zeta3

         case(21)            !-11-11

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/160d0 + (pi**2*ll2**2)/8d0 - 
     &(ll2**4)/8d0 - (pi**2*ll2*llzp)/12d0 + (ll2**3*llzp)/6d
     &0 + zp**8*(-((pi**2)/196608d0) + 266681d0/2890137600d0 
     &+ (ll2)/65536d0 + (ll2**2)/32768d0 - (ll2*llzp)/16384d0
     &) + zp**2*(1d0/16d0 - (pi**2)/192d0 + (ll2)/16d0 + (ll2
     &**2)/32d0 - (ll2*llzp)/16d0) + zp**6*(-((pi**2)/27648d0
     &) + 5269d0/8294400d0 + (ll2)/6912d0 + (ll2**2)/4608d0 -
     & (ll2*llzp)/2304d0) + zp**4*(-((pi**2)/3072d0) + 49d0/9
     &216d0 + (ll2)/512d0 + (ll2**2)/512d0 - (ll2*llzp)/256d0
     &) + zp*(-((pi**2)/24d0) + ll2 + (ll2**2)/4d0 - (ll2*llz
     &p)/2d0) + zp**9*(1077749d0/29262643200d0 - (pi**2)/4976
     &64d0 + (ll2)/186624d0 + (ll2**2)/82944d0 - (ll2*llzp)/4
     &1472d0) + zp**7*(-((pi**2)/75264d0) + 767d0/3225600d0 +
     & (ll2)/21952d0 + (ll2**2)/12544d0 - (ll2*llzp)/6272d0) 
     &+ zp**3*(5d0/288d0 - (pi**2)/864d0 + (ll2)/108d0 + (ll2
     &**2)/144d0 - (ll2*llzp)/72d0) + zp**5*(41d0/23040d0 - (
     &pi**2)/9600d0 + (ll2)/2000d0 + (ll2**2)/1600d0 - (ll2*l
     &lzp)/800d0) - 2*ll2*zeta3 + (llzp*zeta3)/4d0

         case(22)            !-110-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*29d0)/1440d0 - (pi**2*zp)/12d0
     & + (pi**2*ll2**2)/24d0 - (ll2**4)/8d0 + (pi**2*ll2*llzp
     &)/6d0 + zp**6*(-((pi**2)/13824d0) + 667d0/129600d0 - (l
     &lzp)/135d0) + zp**3*(17d0/216d0 - (pi**2)/432d0 - (llzp
     &)/18d0) + zp**7*(-((pi**2)/37632d0) + 2083d0/823200d0 -
     & (13d0*llzp)/2940d0) + zp**8*(6757d0/5017600d0 - (pi**2
     &)/98304d0 - (151d0*llzp)/53760d0) + zp**4*(-((pi**2)/15
     &36d0) + 65d0/2304d0 - (5d0*llzp)/192d0) + zp**5*(-((pi*
     &*2)/4800d0) + 103d0/9000d0 - (llzp)/75d0) + zp**9*(-((p
     &i**2)/248832d0) + 4121d0/5358150d0 - (16d0*llzp)/8505d0
     &) + zp**2*(1d0/4d0 - (pi**2)/96d0 - (llzp)/8d0) - 3*cli
     &4pt5 - (5d0*llzp*zeta3)/8d0

         case(23)            !-1100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*17d0)/360d0 - (pi**2*zp)/4d0 +
     & (-((pi**2)/32d0) - (myi*pi*szp)/8d0)*zp**2 + (-((pi**2
     &)/144d0) + 1d0/36d0 - (myi*pi*szp)/18d0)*zp**3 + (3d0/1
     &28d0 - (pi**2)/512d0 - (myi*pi*5d0*szp)/192d0)*zp**4 + 
     &(-((pi**2)/1600d0) + 1d0/60d0 - (myi*pi*szp)/75d0)*zp**
     &5 + (-((pi**2)/4608d0) + 5d0/432d0 - (myi*pi*szp)/135d0
     &)*zp**6 + (-((pi**2)/12544d0) + 41d0/5040d0 - (myi*pi*1
     &3d0*szp)/2940d0)*zp**7 + (-((pi**2)/32768d0) + 539d0/92
     &160d0 - (myi*pi*151d0*szp)/53760d0)*zp**8 + (22d0/5103d
     &0 - (pi**2)/82944d0 - (myi*pi*16d0*szp)/8505d0)*zp**9 -
     & (myi*pi**3*szp*ll2)/4d0 - (pi**2*ll2**2)/6d0 - (ll2**4
     &)/12d0 + (myi*pi**3*szp*llzp)/12d0 + (pi**2*ll2*llzp)/2
     &d0 - 2*cli4pt5 + (myi*pi*13d0*szp*zeta3)/8d0 - (3d0*llz
     &p*zeta3)/4d0

         case(24)            !-1101

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/24d0 - (pi**2*zp)/24d0 + (pi*
     &*2*ll2**2)/8d0 - (ll2**4)/6d0 + zp**6*(-((pi**2)/27648d
     &0) - 97d0/23040d0 + (ll2)/135d0) + zp**3*(-(1d0/72d0) -
     & (pi**2)/864d0 + (ll2)/18d0) + zp**7*(-((pi**2)/75264d0
     &) - 767d0/282240d0 + (13d0*ll2)/2940d0) + zp**8*(-((pi*
     &*2)/196608d0) - 935d0/516096d0 + (151d0*ll2)/53760d0) +
     & zp**4*(-((pi**2)/3072d0) - 1d0/96d0 + (5d0*ll2)/192d0)
     & + zp**5*(-(1d0/150d0) - (pi**2)/9600d0 + (ll2)/75d0) +
     & zp**9*(-(2041d0/1632960d0) - (pi**2)/497664d0 + (16d0*
     &ll2)/8505d0) + zp**2*(-((pi**2)/192d0) + (ll2)/8d0) + (
     &pi**2*ll2*llzp)/12d0 - 4*cli4pt5 - (21d0*ll2*zeta3)/8d0
     & - (llzp*zeta3)/4d0

         case(25)            !-111-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/288d0) + (pi**2*ll2**2)/24d
     &0 - (ll2**4)/8d0 + zp*((pi**2)/24d0 - (ll2**2)/4d0) + (
     &ll2**3*llzp)/6d0 + zp**6*((pi**2)/27648d0 - 8009d0/8294
     &400d0 - (ll2**2)/4608d0 + (137d0*llzp)/138240d0) + zp**
     &4*((pi**2)/3072d0 - 41d0/4608d0 - (ll2**2)/512d0 + (11d
     &0*llzp)/1536d0) + zp**2*((pi**2)/192d0 - 1d0/8d0 - (ll2
     &**2)/32d0 + (llzp)/16d0) + zp**7*(-(161d0/460800d0) + (
     &pi**2)/75264d0 - (ll2**2)/12544d0 + (llzp)/2560d0) + zp
     &**8*(-(190513d0/1445068800d0) + (pi**2)/196608d0 - (ll2
     &**2)/32768d0 + (363d0*llzp)/2293760d0) + zp**5*(-(13d0/
     &4608d0) + (pi**2)/9600d0 - (ll2**2)/1600d0 + (llzp)/384
     &d0) + zp**3*(-(1d0/32d0) + (pi**2)/864d0 - (ll2**2)/144
     &d0 + (llzp)/48d0) + zp**9*(-(167101d0/3251404800d0) + (
     &pi**2)/497664d0 - (ll2**2)/82944d0 + (761d0*llzp)/11612
     &160d0) - (llzp*zeta3)/8d0

         case(26)            !-1110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*7d0)/360d0) + (myi*pi**3*szp
     &*ll2)/12d0 + (pi**2*ll2**2)/12d0 - (myi*pi*szp*ll2**3)/
     &3d0 + (ll2**4)/12d0 + zp**8*(-(1019d0/1290240d0) + (pi*
     &*2)/196608d0 + (myi*pi*363d0*szp)/2293760d0 - (myi*pi*s
     &zp*ll2)/16384d0) + zp**2*((pi**2)/192d0 + (myi*pi*szp)/
     &16d0 - (myi*pi*szp*ll2)/16d0) + zp**6*((pi**2)/27648d0 
     &- 23d0/8640d0 + (myi*pi*137d0*szp)/138240d0 - (myi*pi*s
     &zp*ll2)/2304d0) + zp**4*((pi**2)/3072d0 - 7d0/768d0 + (
     &myi*pi*11d0*szp)/1536d0 - (myi*pi*szp*ll2)/256d0) + zp*
     &((pi**2)/24d0 - (myi*pi*szp*ll2)/2d0) + zp**9*((pi**2)/
     &497664d0 - 23d0/51030d0 + (myi*pi*761d0*szp)/11612160d0
     & - (myi*pi*szp*ll2)/41472d0) + zp**7*(-(101d0/70560d0) 
     &+ (pi**2)/75264d0 + (myi*pi*szp)/2560d0 - (myi*pi*szp*l
     &l2)/6272d0) + zp**3*(-(1d0/72d0) + (pi**2)/864d0 + (myi
     &*pi*szp)/48d0 - (myi*pi*szp*ll2)/72d0) + zp**5*(-(1d0/2
     &00d0) + (pi**2)/9600d0 + (myi*pi*szp)/384d0 - (myi*pi*s
     &zp*ll2)/800d0) - (pi**2*ll2*llzp)/12d0 + (myi*pi*szp*ll
     &2**2*llzp)/2d0 + 2*cli4pt5 - (myi*pi*szp*zeta3)/8d0 + (
     &llzp*zeta3)/8d0

         case(27)            !-1111

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/90d0) - (pi**2*ll2**2)/12d0
     & + (zp*ll2**2)/4d0 + (ll2**4)/6d0 + zp**7*(29d0/80640d0
     & - (ll2)/2560d0 + (ll2**2)/12544d0) + zp**3*(1d0/144d0 
     &- (ll2)/48d0 + (ll2**2)/144d0) + zp**5*(7d0/3840d0 - (l
     &l2)/384d0 + (ll2**2)/1600d0) + zp**8*(469d0/2949120d0 -
     & (363d0*ll2)/2293760d0 + (ll2**2)/32768d0) + zp**2*(-((
     &ll2)/16d0) + (ll2**2)/32d0) + zp**6*(5d0/6144d0 - (137d
     &0*ll2)/138240d0 + (ll2**2)/4608d0) + zp**4*(1d0/256d0 -
     & (11d0*ll2)/1536d0 + (ll2**2)/512d0) + zp**9*(29531d0/4
     &18037760d0 - (761d0*ll2)/11612160d0 + (ll2**2)/82944d0)
     & - (ll2**3*llzp)/6d0 + cli4pt5 + ll2*zeta3

         case(28)            !0-1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/90d0) + zp**2*(1d0/16d0 - (
     &llzp)/8d0 + (llzp**2)/8d0 - (llzp**3)/12d0) + zp**3*(1d
     &0/81d0 - (llzp)/27d0 + (llzp**2)/18d0 - (llzp**3)/18d0)
     & + zp**4*(1d0/256d0 - (llzp)/64d0 + (llzp**2)/32d0 - (l
     &lzp**3)/24d0) + zp**5*(1d0/625d0 - (llzp)/125d0 + (llzp
     &**2)/50d0 - (llzp**3)/30d0) + zp**6*(1d0/1296d0 - (llzp
     &)/216d0 + (llzp**2)/72d0 - (llzp**3)/36d0) + zp**7*(1d0
     &/2401d0 - (llzp)/343d0 + (llzp**2)/98d0 - (llzp**3)/42d
     &0) + zp**8*(1d0/4096d0 - (llzp)/512d0 + (llzp**2)/128d0
     & - (llzp**3)/48d0) + zp**9*(1d0/6561d0 - (llzp)/729d0 +
     & (llzp**2)/162d0 - (llzp**3)/54d0) + zp*(1 - llzp + (ll
     &zp**2)/2d0 - (llzp**3)/6d0)

         case(29)            !0-1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4)/72d0) + zp*((pi**2)/6d0 - m
     &yi*pi*szp - (pi**2*llzp)/6d0 + myi*pi*szp*llzp - (myi*p
     &i*szp*llzp**2)/2d0 - zeta3) + myi*pi*szp*zeta3 + zp**2*
     &((pi**2)/24d0 + 1d0/2d0 - (myi*pi*szp)/8d0 - (pi**2*llz
     &p)/12d0 + (myi*pi*szp*llzp)/4d0 - (myi*pi*szp*llzp**2)/
     &4d0 - (zeta3)/2d0) + zp**3*((pi**2)/54d0 + 3d0/8d0 - (m
     &yi*pi*szp)/27d0 - (pi**2*llzp)/18d0 + (myi*pi*szp*llzp)
     &/9d0 - (myi*pi*szp*llzp**2)/6d0 - (zeta3)/3d0) + zp**4*
     &(251d0/864d0 + (pi**2)/96d0 - (myi*pi*szp)/64d0 - (pi**
     &2*llzp)/24d0 + (myi*pi*szp*llzp)/16d0 - (myi*pi*szp*llz
     &p**2)/8d0 - (zeta3)/4d0) + zp**5*((pi**2)/150d0 + 407d0
     &/1728d0 - (myi*pi*szp)/125d0 - (pi**2*llzp)/30d0 + (myi
     &*pi*szp*llzp)/25d0 - (myi*pi*szp*llzp**2)/10d0 - (zeta3
     &)/5d0) + zp**6*((pi**2)/216d0 + 256103d0/1296000d0 - (m
     &yi*pi*szp)/216d0 - (pi**2*llzp)/36d0 + (myi*pi*szp*llzp
     &)/36d0 - (myi*pi*szp*llzp**2)/12d0 - (zeta3)/6d0) + zp*
     &*7*((pi**2)/294d0 + 4081d0/24000d0 - (myi*pi*szp)/343d0
     & - (pi**2*llzp)/42d0 + (myi*pi*szp*llzp)/49d0 - (myi*pi
     &*szp*llzp**2)/14d0 - (zeta3)/7d0) + zp**8*((pi**2)/384d
     &0 + 9822481d0/65856000d0 - (myi*pi*szp)/512d0 - (pi**2*
     &llzp)/48d0 + (myi*pi*szp*llzp)/64d0 - (myi*pi*szp*llzp*
     &*2)/16d0 - (zeta3)/8d0) + zp**9*((pi**2)/486d0 + 787084
     &73d0/592704000d0 - (myi*pi*szp)/729d0 - (pi**2*llzp)/54
     &d0 + (myi*pi*szp*llzp)/81d0 - (myi*pi*szp*llzp**2)/18d0
     & - (zeta3)/9d0)

         case(30)            !0-1-11

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/80d0 - (pi**2*ll2**2)/12d0 - 
     &(ll2**4)/24d0 - cli4pt5 - (7d0*ll2*zeta3)/8d0 + zp**2*(
     &-((pi**2)/48d0) - 1d0/4d0 - (pi**2*ll2)/24d0 + (ll2)/8d
     &0 + (ll2**2)/8d0 + (ll2**3)/12d0 + (pi**2*llzp)/24d0 - 
     &(ll2*llzp)/4d0 - (ll2**2*llzp)/4d0 + (ll2*llzp**2)/4d0 
     &+ (7d0*zeta3)/16d0) + zp**3*(-((pi**2)/108d0) - 17d0/96
     &d0 + (ll2)/27d0 - (pi**2*ll2)/36d0 + (ll2**2)/18d0 + (l
     &l2**3)/18d0 + (pi**2*llzp)/36d0 - (ll2*llzp)/9d0 - (ll2
     &**2*llzp)/6d0 + (ll2*llzp**2)/6d0 + (7d0*zeta3)/24d0) +
     & zp**4*(-((pi**2)/192d0) - 463d0/3456d0 - (pi**2*ll2)/4
     &8d0 + (ll2)/64d0 + (ll2**2)/32d0 + (ll2**3)/24d0 + (pi*
     &*2*llzp)/48d0 - (ll2*llzp)/16d0 - (ll2**2*llzp)/8d0 + (
     &ll2*llzp**2)/8d0 + (7d0*zeta3)/32d0) + zp**5*(-(14843d0
     &/138240d0) - (pi**2)/300d0 + (ll2)/125d0 - (pi**2*ll2)/
     &60d0 + (ll2**2)/50d0 + (ll2**3)/30d0 + (pi**2*llzp)/60d
     &0 - (ll2*llzp)/25d0 - (ll2**2*llzp)/10d0 + (ll2*llzp**2
     &)/10d0 + (7d0*zeta3)/40d0) + zp**6*(-(1856239d0/2073600
     &0d0) - (pi**2)/432d0 + (ll2)/216d0 - (pi**2*ll2)/72d0 +
     & (ll2**2)/72d0 + (ll2**3)/36d0 + (pi**2*llzp)/72d0 - (l
     &l2*llzp)/36d0 - (ll2**2*llzp)/12d0 + (ll2*llzp**2)/12d0
     & + (7d0*zeta3)/48d0) + zp**8*(-((pi**2)/768d0) - 636802
     &727d0/9483264000d0 + (ll2)/512d0 - (pi**2*ll2)/96d0 + (
     &ll2**2)/128d0 + (ll2**3)/48d0 + (pi**2*llzp)/96d0 - (ll
     &2*llzp)/64d0 - (ll2**2*llzp)/16d0 + (ll2*llzp**2)/16d0 
     &+ (7d0*zeta3)/64d0) + zp**9*(-(81511906681d0/1365590016
     &000d0) - (pi**2)/972d0 - (pi**2*ll2)/108d0 + (ll2)/729d
     &0 + (ll2**2)/162d0 + (ll2**3)/54d0 + (pi**2*llzp)/108d0
     & - (ll2*llzp)/81d0 - (ll2**2*llzp)/18d0 + (ll2*llzp**2)
     &/18d0 + (7d0*zeta3)/72d0) + zp**7*(-(1856489d0/24192000
     &d0) - (pi**2)/588d0 + (ll2)/343d0 - (pi**2*ll2)/84d0 + 
     &(ll2**2)/98d0 + (ll2**3)/42d0 + (pi**2*llzp)/84d0 - (ll
     &2*llzp)/49d0 - (ll2**2*llzp)/14d0 + (ll2*llzp**2)/14d0 
     &+ (zeta3)/8d0) + zp*(-((pi**2)/12d0) + ll2 - (pi**2*ll2
     &)/12d0 + (ll2**2)/2d0 + (ll2**3)/6d0 + (pi**2*llzp)/12d
     &0 - ll2*llzp - (ll2**2*llzp)/2d0 + (ll2*llzp**2)/2d0 + 
     &(7d0*zeta3)/8d0)

         case(31)            !0-10-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/120d0 + zp**2*(-((pi**2)/24d0
     &) - 5d0/4d0 + (pi**2*llzp)/12d0 + (llzp)/2d0 + zeta3) +
     & zp*(-((pi**2)/6d0) + (pi**2*llzp)/6d0 + 2*zeta3) + zp*
     &*4*(-(1151d0/1728d0) - (pi**2)/96d0 + (pi**2*llzp)/24d0
     & + (49d0*llzp)/144d0 + (zeta3)/2d0) + zp**6*(-((pi**2)/
     &216d0) - 17653d0/40500d0 + (pi**2*llzp)/36d0 + (5269d0*
     &llzp)/21600d0 + (zeta3)/3d0) + zp**3*(-((pi**2)/54d0) -
     & 8d0/9d0 + (pi**2*llzp)/18d0 + (5d0*llzp)/12d0 + (2d0*z
     &eta3)/3d0) + zp**8*(-((pi**2)/384d0) - 127203607d0/3951
     &36000d0 + (266681d0*llzp)/1411200d0 + (pi**2*llzp)/48d0
     & + (zeta3)/4d0) + zp**5*(-((pi**2)/150d0) - 2281d0/4320
     &d0 + (pi**2*llzp)/30d0 + (41d0*llzp)/144d0 + (2d0*zeta3
     &)/5d0) + zp**7*(-((pi**2)/294d0) - 93371d0/252000d0 + (
     &pi**2*llzp)/42d0 + (767d0*llzp)/3600d0 + (2d0*zeta3)/7d
     &0) + zp**9*(-((pi**2)/486d0) - 2276013631d0/8001504000d
     &0 + (pi**2*llzp)/54d0 + (1077749d0*llzp)/6350400d0 + (2
     &d0*zeta3)/9d0)

         case(32)            !0-100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/20d0 + 2*myi*pi*szp*zeta3 + z
     &p*(-((pi**2)/2d0) - (myi*pi**3*szp)/6d0 + (pi**2*llzp)/
     &2d0 + zeta3) + zp**2*(-((pi**2)/8d0) - (myi*pi**3*szp)/
     &12d0 + (myi*pi*szp)/2d0 + (pi**2*llzp)/4d0 + (zeta3)/2d
     &0) + zp**3*(-(1d0/12d0) - (pi**2)/18d0 - (myi*pi**3*szp
     &)/18d0 + (myi*pi*5d0*szp)/12d0 + (pi**2*llzp)/6d0 + (ze
     &ta3)/3d0) + zp**4*(-((pi**2)/32d0) - 5d0/48d0 - (myi*pi
     &**3*szp)/24d0 + (myi*pi*49d0*szp)/144d0 + (pi**2*llzp)/
     &8d0 + (zeta3)/4d0) + zp**5*(-(17d0/160d0) - (pi**2)/50d
     &0 - (myi*pi**3*szp)/30d0 + (myi*pi*41d0*szp)/144d0 + (p
     &i**2*llzp)/10d0 + (zeta3)/5d0) + zp**6*(-(59d0/576d0) -
     & (pi**2)/72d0 - (myi*pi**3*szp)/36d0 + (myi*pi*5269d0*s
     &zp)/21600d0 + (pi**2*llzp)/12d0 + (zeta3)/6d0) + zp**7*
     &(-(2929d0/30240d0) - (pi**2)/98d0 - (myi*pi**3*szp)/42d
     &0 + (myi*pi*767d0*szp)/3600d0 + (pi**2*llzp)/14d0 + (ze
     &ta3)/7d0) + zp**8*(-((pi**2)/128d0) - 629d0/6912d0 + (m
     &yi*pi*266681d0*szp)/1411200d0 - (myi*pi**3*szp)/48d0 + 
     &(pi**2*llzp)/16d0 + (zeta3)/8d0) + zp**9*(-((pi**2)/162
     &d0) - 185921d0/2177280d0 - (myi*pi**3*szp)/54d0 + (myi*
     &pi*1077749d0*szp)/6350400d0 + (pi**2*llzp)/18d0 + (zeta
     &3)/9d0)

         case(33)            !0-101

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*71d0)/1440d0 + (pi**2*ll2**2)/
     &6d0 - (ll2**4)/6d0 - 4*cli4pt5 - (7d0*ll2*zeta3)/2d0 + 
     &zp**2*(-((pi**2)/48d0) - (ll2)/2d0 + (pi**2*llzp)/24d0 
     &+ (5d0*zeta3)/16d0) + zp**3*(-((pi**2)/108d0) + 1d0/24d
     &0 - (5d0*ll2)/12d0 + (pi**2*llzp)/36d0 + (5d0*zeta3)/24
     &d0) + zp**4*(-((pi**2)/192d0) + 7d0/144d0 - (49d0*ll2)/
     &144d0 + (pi**2*llzp)/48d0 + (5d0*zeta3)/32d0) + zp**6*(
     &-((pi**2)/432d0) + 3793d0/86400d0 - (5269d0*ll2)/21600d
     &0 + (pi**2*llzp)/72d0 + (5d0*zeta3)/48d0) + zp**7*(4882
     &1d0/1209600d0 - (pi**2)/588d0 - (767d0*ll2)/3600d0 + (p
     &i**2*llzp)/84d0 + (5d0*zeta3)/56d0) + zp**8*(2511659d0/
     &67737600d0 - (pi**2)/768d0 - (266681d0*ll2)/1411200d0 +
     & (pi**2*llzp)/96d0 + (5d0*zeta3)/64d0) + zp**9*(1041298
     &1d0/304819200d0 - (pi**2)/972d0 - (1077749d0*ll2)/63504
     &00d0 + (pi**2*llzp)/108d0 + (5d0*zeta3)/72d0) + zp**5*(
     &-((pi**2)/300d0) + 17d0/360d0 - (41d0*ll2)/144d0 + (pi*
     &*2*llzp)/60d0 + (zeta3)/8d0) + zp*(-((pi**2)/12d0) + (p
     &i**2*llzp)/12d0 + (5d0*zeta3)/8d0)

         case(34)            !0-11-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*7d0)/288d0) + (pi**2*ll2**2)
     &/8d0 + (ll2**4)/8d0 + 3*cli4pt5 + zp**7*((pi**2)/588d0 
     &+ 14001083d0/84672000d0 + (pi**2*ll2)/42d0 - (ll2**2)/9
     &8d0 - (ll2**3)/21d0 - (5587d0*llzp)/67200d0 - (pi**2*ll
     &zp)/84d0 + (ll2**2*llzp)/14d0 - (zeta3)/4d0) + zp**3*((
     &pi**2)/108d0 + 5d0/12d0 + (pi**2*ll2)/18d0 - (ll2**2)/1
     &8d0 - (ll2**3)/9d0 - (pi**2*llzp)/36d0 - (3d0*llzp)/16d
     &0 + (ll2**2*llzp)/6d0 - (7d0*zeta3)/12d0) + zp**4*((pi*
     &*2)/192d0 + 2101d0/6912d0 + (pi**2*ll2)/24d0 - (ll2**2)
     &/32d0 - (ll2**3)/12d0 - (pi**2*llzp)/48d0 - (83d0*llzp)
     &/576d0 + (ll2**2*llzp)/8d0 - (7d0*zeta3)/16d0) + zp**5*
     &((pi**2)/300d0 + 82237d0/345600d0 + (pi**2*ll2)/30d0 - 
     &(ll2**2)/50d0 - (ll2**3)/15d0 - (1337d0*llzp)/11520d0 -
     & (pi**2*llzp)/60d0 + (ll2**2*llzp)/10d0 - (7d0*zeta3)/2
     &0d0) + zp**6*((pi**2)/432d0 + 505931d0/2592000d0 + (pi*
     &*2*ll2)/36d0 - (ll2**2)/72d0 - (ll2**3)/18d0 - (33497d0
     &*llzp)/345600d0 - (pi**2*llzp)/72d0 + (ll2**2*llzp)/12d
     &0 - (7d0*zeta3)/24d0) + zp**8*(169983053d0/1185408000d0
     & + (pi**2)/768d0 + (pi**2*ll2)/48d0 - (ll2**2)/128d0 - 
     &(ll2**3)/24d0 - (136919d0*llzp)/1881600d0 - (pi**2*llzp
     &)/96d0 + (ll2**2*llzp)/16d0 - (7d0*zeta3)/32d0) + zp**9
     &*(86419598141d0/682795008000d0 + (pi**2)/972d0 + (pi**2
     &*ll2)/54d0 - (ll2**2)/162d0 - (ll2**3)/27d0 - (pi**2*ll
     &zp)/108d0 - (35054939d0*llzp)/541900800d0 + (ll2**2*llz
     &p)/18d0 - (7d0*zeta3)/36d0) + zp*((pi**2)/12d0 + (pi**2
     &*ll2)/6d0 - (ll2**2)/2d0 - (ll2**3)/3d0 - (pi**2*llzp)/
     &12d0 + (ll2**2*llzp)/2d0 - (7d0*zeta3)/4d0) + zp**2*((p
     &i**2)/48d0 + 5d0/8d0 + (pi**2*ll2)/12d0 - (ll2**2)/8d0 
     &- (ll2**3)/6d0 - (pi**2*llzp)/24d0 - (llzp)/4d0 + (ll2*
     &*2*llzp)/4d0 - (7d0*zeta3)/8d0)

         case(35)            !0-110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*11d0)/480d0) + (myi*pi**3*sz
     &p*ll2)/4d0 - (pi**2*ll2**2)/6d0 + (ll2**4)/6d0 + 4*cli4
     &pt5 - myi*pi*szp*zeta3 + zp**2*((pi**2)/48d0 + (myi*pi*
     &*3*szp)/24d0 - (myi*pi*szp)/4d0 + (pi**2*ll2)/8d0 - (my
     &i*pi*szp*ll2)/4d0 - (myi*pi*szp*ll2**2)/4d0 - (pi**2*ll
     &zp)/24d0 + (myi*pi*szp*ll2*llzp)/2d0 - (13d0*zeta3)/16d
     &0) + zp**3*((pi**2)/108d0 + 1d0/24d0 + (myi*pi**3*szp)/
     &36d0 - (myi*pi*3d0*szp)/16d0 + (pi**2*ll2)/12d0 - (myi*
     &pi*szp*ll2)/9d0 - (myi*pi*szp*ll2**2)/6d0 - (pi**2*llzp
     &)/36d0 + (myi*pi*szp*ll2*llzp)/3d0 - (13d0*zeta3)/24d0)
     & + zp**4*((pi**2)/192d0 + 13d0/288d0 + (myi*pi**3*szp)/
     &48d0 - (myi*pi*83d0*szp)/576d0 + (pi**2*ll2)/16d0 - (my
     &i*pi*szp*ll2)/16d0 - (myi*pi*szp*ll2**2)/8d0 - (pi**2*l
     &lzp)/48d0 + (myi*pi*szp*ll2*llzp)/4d0 - (13d0*zeta3)/32
     &d0) + zp**5*(119d0/2880d0 + (pi**2)/300d0 - (myi*pi*133
     &7d0*szp)/11520d0 + (myi*pi**3*szp)/60d0 + (pi**2*ll2)/2
     &0d0 - (myi*pi*szp*ll2)/25d0 - (myi*pi*szp*ll2**2)/10d0 
     &- (pi**2*llzp)/60d0 + (myi*pi*szp*ll2*llzp)/5d0 - (13d0
     &*zeta3)/40d0) + zp**6*((pi**2)/432d0 + 3167d0/86400d0 -
     & (myi*pi*33497d0*szp)/345600d0 + (myi*pi**3*szp)/72d0 +
     & (pi**2*ll2)/24d0 - (myi*pi*szp*ll2)/36d0 - (myi*pi*szp
     &*ll2**2)/12d0 - (pi**2*llzp)/72d0 + (myi*pi*szp*ll2*llz
     &p)/6d0 - (13d0*zeta3)/48d0) + zp**7*(1403d0/43200d0 + (
     &pi**2)/588d0 - (myi*pi*5587d0*szp)/67200d0 + (myi*pi**3
     &*szp)/84d0 + (pi**2*ll2)/28d0 - (myi*pi*szp*ll2)/49d0 -
     & (myi*pi*szp*ll2**2)/14d0 - (pi**2*llzp)/84d0 + (myi*pi
     &*szp*ll2*llzp)/7d0 - (13d0*zeta3)/56d0) + zp**8*(490589
     &d0/16934400d0 + (pi**2)/768d0 - (myi*pi*136919d0*szp)/1
     &881600d0 + (myi*pi**3*szp)/96d0 + (pi**2*ll2)/32d0 - (m
     &yi*pi*szp*ll2)/64d0 - (myi*pi*szp*ll2**2)/16d0 - (pi**2
     &*llzp)/96d0 + (myi*pi*szp*ll2*llzp)/8d0 - (13d0*zeta3)/
     &64d0) + zp**9*(3972277d0/152409600d0 + (pi**2)/972d0 + 
     &(myi*pi**3*szp)/108d0 - (myi*pi*35054939d0*szp)/5419008
     &00d0 + (pi**2*ll2)/36d0 - (myi*pi*szp*ll2)/81d0 - (myi*
     &pi*szp*ll2**2)/18d0 - (pi**2*llzp)/108d0 + (myi*pi*szp*
     &ll2*llzp)/9d0 - (13d0*zeta3)/72d0) + zp*((pi**2)/12d0 +
     & (myi*pi**3*szp)/12d0 + (pi**2*ll2)/4d0 - myi*pi*szp*ll
     &2 - (myi*pi*szp*ll2**2)/2d0 - (pi**2*llzp)/12d0 + myi*p
     &i*szp*ll2*llzp - (13d0*zeta3)/8d0)

         case(36)            !0-111

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*7d0)/288d0) - (pi**2*5d0*ll2
     &**2)/24d0 + (ll2**4)/12d0 + 2*cli4pt5 + (21d0*ll2*zeta3
     &)/8d0 + zp**2*(-((pi**2*ll2)/24d0) + (ll2)/4d0 + (ll2**
     &2)/8d0 + (ll2**3)/6d0 - (ll2**2*llzp)/4d0 + (zeta3)/16d
     &0) + zp**3*(-(1d0/48d0) - (pi**2*ll2)/36d0 + (3d0*ll2)/
     &16d0 + (ll2**2)/18d0 + (ll2**3)/9d0 - (ll2**2*llzp)/6d0
     & + (zeta3)/24d0) + zp**4*(-(1d0/48d0) - (pi**2*ll2)/48d
     &0 + (83d0*ll2)/576d0 + (ll2**2)/32d0 + (ll2**3)/12d0 - 
     &(ll2**2*llzp)/8d0 + (zeta3)/32d0) + zp**5*(-(139d0/7680
     &d0) + (1337d0*ll2)/11520d0 - (pi**2*ll2)/60d0 + (ll2**2
     &)/50d0 + (ll2**3)/15d0 - (ll2**2*llzp)/10d0 + (zeta3)/4
     &0d0) + zp**6*(-(143d0/9216d0) + (33497d0*ll2)/345600d0 
     &- (pi**2*ll2)/72d0 + (ll2**2)/72d0 + (ll2**3)/18d0 - (l
     &l2**2*llzp)/12d0 + (zeta3)/48d0) + zp**7*(-(13007d0/967
     &680d0) + (5587d0*ll2)/67200d0 - (pi**2*ll2)/84d0 + (ll2
     &**2)/98d0 + (ll2**3)/21d0 - (ll2**2*llzp)/14d0 + (zeta3
     &)/56d0) + zp**8*(-(13061d0/1105920d0) + (136919d0*ll2)/
     &1881600d0 - (pi**2*ll2)/96d0 + (ll2**2)/128d0 + (ll2**3
     &)/24d0 - (ll2**2*llzp)/16d0 + (zeta3)/64d0) + zp**9*(-(
     &5861129d0/557383680d0) - (pi**2*ll2)/108d0 + (35054939d
     &0*ll2)/541900800d0 + (ll2**2)/162d0 + (ll2**3)/27d0 - (
     &ll2**2*llzp)/18d0 + (zeta3)/72d0) + zp*(-((pi**2*ll2)/1
     &2d0) + (ll2**2)/2d0 + (ll2**3)/3d0 - (ll2**2*llzp)/2d0 
     &+ (zeta3)/8d0)

         case(37)            !00-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/360d0 - zp*zeta3 + zp**2*(7d0
     &/8d0 - (3d0*llzp)/4d0 + (llzp**2)/4d0 - (zeta3)/2d0) + 
     &zp**3*(41d0/72d0 - (7d0*llzp)/12d0 + (llzp**2)/4d0 - (z
     &eta3)/3d0) + zp**4*(1397d0/3456d0 - (131d0*llzp)/288d0 
     &+ (11d0*llzp**2)/48d0 - (zeta3)/4d0) + zp**5*(2671d0/86
     &40d0 - (53d0*llzp)/144d0 + (5d0*llzp**2)/24d0 - (zeta3)
     &/5d0) + zp**6*(322493d0/1296000d0 - (2213d0*llzp)/7200d
     &0 + (137d0*llzp**2)/720d0 - (zeta3)/6d0) + zp**7*(10464
     &1d0/504000d0 - (947d0*llzp)/3600d0 + (7d0*llzp**2)/40d0
     & - (zeta3)/7d0) + zp**8*(140539517d0/790272000d0 - (647
     &707d0*llzp)/2822400d0 + (363d0*llzp**2)/2240d0 - (zeta3
     &)/8d0) + zp**9*(2486560891d0/16003008000d0 - (1290829d0
     &*llzp)/6350400d0 + (761d0*llzp**2)/5040d0 - (zeta3)/9d0
     &)

         case(38)            !00-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/30d0 + zp*((myi*pi**3*szp)/6d
     &0 - 2*zeta3) + zp**2*((pi**2)/12d0 + (myi*pi**3*szp)/12
     &d0 - (myi*pi*3d0*szp)/4d0 + (myi*pi*szp*llzp)/2d0 - zet
     &a3) - myi*pi*szp*zeta3 + zp**4*((pi**2*11d0)/144d0 - 11
     &d0/48d0 + (myi*pi**3*szp)/24d0 - (myi*pi*131d0*szp)/288
     &d0 + (myi*pi*11d0*szp*llzp)/24d0 - (zeta3)/2d0) + zp**6
     &*((pi**2*137d0)/2160d0 - 37d0/144d0 + (myi*pi**3*szp)/3
     &6d0 - (myi*pi*2213d0*szp)/7200d0 + (myi*pi*137d0*szp*ll
     &zp)/360d0 - (zeta3)/3d0) + zp**3*((pi**2)/12d0 - 1d0/6d
     &0 + (myi*pi**3*szp)/18d0 - (myi*pi*7d0*szp)/12d0 + (myi
     &*pi*szp*llzp)/2d0 - (2d0*zeta3)/3d0) + zp**8*((pi**2*12
     &1d0)/2240d0 - 43171d0/172800d0 + (myi*pi**3*szp)/48d0 -
     & (myi*pi*647707d0*szp)/2822400d0 + (myi*pi*363d0*szp*ll
     &zp)/1120d0 - (zeta3)/4d0) + zp**5*(-(181d0/720d0) + (pi
     &**2*5d0)/72d0 + (myi*pi**3*szp)/30d0 - (myi*pi*53d0*szp
     &)/144d0 + (myi*pi*5d0*szp*llzp)/12d0 - (2d0*zeta3)/5d0)
     & + zp**7*(-(38569d0/151200d0) + (pi**2*7d0)/120d0 + (my
     &i*pi**3*szp)/42d0 - (myi*pi*947d0*szp)/3600d0 + (myi*pi
     &*7d0*szp*llzp)/20d0 - (2d0*zeta3)/7d0) + zp**9*((pi**2*
     &761d0)/15120d0 - 9261559d0/38102400d0 + (myi*pi**3*szp)
     &/54d0 - (myi*pi*1290829d0*szp)/6350400d0 + (myi*pi*761d
     &0*szp*llzp)/2520d0 - (2d0*zeta3)/9d0)

         case(39)            !00-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*19d0)/1440d0) + (7d0*ll2*zet
     &a3)/4d0 + zp*(-((pi**2*ll2)/4d0) + zeta3) + zp**2*(-((p
     &i**2)/24d0) + (3d0*ll2)/4d0 - (pi**2*ll2)/8d0 + (ll2**2
     &)/4d0 - (ll2*llzp)/2d0 + (zeta3)/2d0) + zp**3*(1d0/12d0
     & - (pi**2)/24d0 - (pi**2*ll2)/12d0 + (7d0*ll2)/12d0 + (
     &ll2**2)/4d0 - (ll2*llzp)/2d0 + (zeta3)/3d0) + zp**4*(-(
     &(pi**2*11d0)/288d0) + 7d0/64d0 - (pi**2*ll2)/16d0 + (13
     &1d0*ll2)/288d0 + (11d0*ll2**2)/48d0 - (11d0*ll2*llzp)/2
     &4d0 + (zeta3)/4d0) + zp**5*(-((pi**2*5d0)/144d0) + 67d0
     &/576d0 - (pi**2*ll2)/20d0 + (53d0*ll2)/144d0 + (5d0*ll2
     &**2)/24d0 - (5d0*ll2*llzp)/12d0 + (zeta3)/5d0) + zp**6*
     &(-((pi**2*137d0)/4320d0) + 893d0/7680d0 - (pi**2*ll2)/2
     &4d0 + (2213d0*ll2)/7200d0 + (137d0*ll2**2)/720d0 - (137
     &d0*ll2*llzp)/360d0 + (zeta3)/6d0) + zp**7*(274607d0/241
     &9200d0 - (pi**2*7d0)/240d0 - (pi**2*ll2)/28d0 + (947d0*
     &ll2)/3600d0 + (7d0*ll2**2)/40d0 - (7d0*ll2*llzp)/20d0 +
     & (zeta3)/7d0) + zp**8*(2123381d0/19353600d0 - (pi**2*12
     &1d0)/4480d0 - (pi**2*ll2)/32d0 + (647707d0*ll2)/2822400
     &d0 + (363d0*ll2**2)/2240d0 - (363d0*ll2*llzp)/1120d0 + 
     &(zeta3)/8d0) + zp**9*(-((pi**2*761d0)/30240d0) + 804796
     &9d0/76204800d0 - (pi**2*ll2)/36d0 + (1290829d0*ll2)/635
     &0400d0 + (761d0*ll2**2)/5040d0 - (761d0*ll2*llzp)/2520d
     &0 + (zeta3)/9d0)

         case(40)            !000-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/90d0) + zp*zeta3 + zp**2*(-
     &((pi**2)/12d0) + (zeta3)/2d0) + zp**3*(-((pi**2)/12d0) 
     &+ 11d0/36d0 - (llzp)/6d0 + (zeta3)/3d0) + zp**4*(-((pi*
     &*2*11d0)/144d0) + 19d0/48d0 - (llzp)/4d0 + (zeta3)/4d0)
     & + zp**5*(599d0/1440d0 - (pi**2*5d0)/72d0 - (7d0*llzp)/
     &24d0 + (zeta3)/5d0) + zp**6*(-((pi**2*137d0)/2160d0) + 
     &79d0/192d0 - (5d0*llzp)/16d0 + (zeta3)/6d0) + zp**7*(-(
     &(pi**2*7d0)/120d0) + 3343d0/8400d0 - (29d0*llzp)/90d0 +
     & (zeta3)/7d0) + zp**8*(-((pi**2*121d0)/2240d0) + 21977d
     &0/57600d0 - (469d0*llzp)/1440d0 + (zeta3)/8d0) + zp**9*
     &(-((pi**2*761d0)/15120d0) + 83359739d0/228614400d0 - (2
     &9531d0*llzp)/90720d0 + (zeta3)/9d0)

         case(41)            !0000

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4)/24d0 + (myi*pi**3*szp*zp)/6d0
     & + (-((pi**2)/4d0) + (myi*pi**3*szp)/12d0)*zp**2 + (-((
     &pi**2)/4d0) + (myi*pi**3*szp)/18d0 - (myi*pi*szp)/6d0)*
     &zp**3 + (1d0/24d0 - (pi**2*11d0)/48d0 + (myi*pi**3*szp)
     &/24d0 - (myi*pi*szp)/4d0)*zp**4 + (1d0/12d0 - (pi**2*5d
     &0)/24d0 + (myi*pi**3*szp)/30d0 - (myi*pi*7d0*szp)/24d0)
     &*zp**5 + (17d0/144d0 - (pi**2*137d0)/720d0 + (myi*pi**3
     &*szp)/36d0 - (myi*pi*5d0*szp)/16d0)*zp**6 + (-((pi**2*7
     &d0)/40d0) + 7d0/48d0 + (myi*pi**3*szp)/42d0 - (myi*pi*2
     &9d0*szp)/90d0)*zp**7 + (-((pi**2*363d0)/2240d0) + 967d0
     &/5760d0 - (myi*pi*469d0*szp)/1440d0 + (myi*pi**3*szp)/4
     &8d0)*zp**8 + (-((pi**2*761d0)/5040d0) + 89d0/480d0 + (m
     &yi*pi**3*szp)/54d0 - (myi*pi*29531d0*szp)/90720d0)*zp**
     &9

         case(42)            !0001

            zp = x+1d0

            ris = -((pi**4*7d0)/720d0) + (3d0*zp*zeta3)
     &/4d0 + zp**9*(-(112391d0/1451520d0) - (pi**2*761d0)/302
     &40d0 + (29531d0*ll2)/90720d0 + (zeta3)/12d0) + zp**4*(-
     &((pi**2*11d0)/288d0) - 1d0/48d0 + (ll2)/4d0 + (3d0*zeta
     &3)/16d0) + zp**5*(-(19d0/480d0) - (pi**2*5d0)/144d0 + (
     &7d0*ll2)/24d0 + (3d0*zeta3)/20d0) + zp**7*(-(2591d0/403
     &20d0) - (pi**2*7d0)/240d0 + (29d0*ll2)/90d0 + (3d0*zeta
     &3)/28d0) + zp**8*(-(23d0/320d0) - (pi**2*121d0)/4480d0 
     &+ (469d0*ll2)/1440d0 + (3d0*zeta3)/32d0) + zp**3*(-((pi
     &**2)/24d0) + (ll2)/6d0 + (zeta3)/4d0) + zp**6*(-((pi**2
     &*137d0)/4320d0) - 31d0/576d0 + (5d0*ll2)/16d0 + (zeta3)
     &/8d0) + zp**2*(-((pi**2)/24d0) + (3d0*zeta3)/8d0)

         case(43)            !001-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/180d0) - (pi**2*ll2**2)/12d
     &0 + (ll2**4)/12d0 + 2*cli4pt5 + zp**2*((pi**2)/24d0 + (
     &pi**2*ll2)/8d0 - (ll2**2)/4d0 - (13d0*zeta3)/16d0) + zp
     &**3*((pi**2)/24d0 - 11d0/72d0 + (pi**2*ll2)/12d0 - (ll2
     &**2)/4d0 + (llzp)/12d0 - (13d0*zeta3)/24d0) + zp**4*(-(
     &215d0/1152d0) + (pi**2*11d0)/288d0 + (pi**2*ll2)/16d0 -
     & (11d0*ll2**2)/48d0 + (11d0*llzp)/96d0 - (13d0*zeta3)/3
     &2d0) + zp**5*((pi**2*5d0)/144d0 - 181d0/960d0 + (pi**2*
     &ll2)/20d0 - (5d0*ll2**2)/24d0 + (llzp)/8d0 - (13d0*zeta
     &3)/40d0) + zp**6*(-(2321d0/12800d0) + (pi**2*137d0)/432
     &0d0 + (pi**2*ll2)/24d0 - (137d0*ll2**2)/720d0 + (731d0*
     &llzp)/5760d0 - (13d0*zeta3)/48d0) + zp**7*((pi**2*7d0)/
     &240d0 - 138503d0/806400d0 + (pi**2*ll2)/28d0 - (7d0*ll2
     &**2)/40d0 + (721d0*llzp)/5760d0 - (13d0*zeta3)/56d0) + 
     &zp**8*(-(182923d0/1128960d0) + (pi**2*121d0)/4480d0 + (
     &pi**2*ll2)/32d0 - (363d0*ll2**2)/2240d0 + (3931d0*llzp)
     &/32256d0 - (13d0*zeta3)/64d0) + zp**9*((pi**2*761d0)/30
     &240d0 - 139798291d0/914457600d0 + (pi**2*ll2)/36d0 - (7
     &61d0*ll2**2)/5040d0 + (42799d0*llzp)/362880d0 - (13d0*z
     &eta3)/72d0) + zp*((pi**2*ll2)/4d0 - (13d0*zeta3)/8d0)

         case(44)            !0010

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4*7d0)/240d0 - (myi*pi*3d0*szp*z
     &eta3)/4d0 + zp**3*((pi**2)/24d0 + (myi*pi*szp)/12d0 + (
     &myi*pi**3*szp)/36d0 - (myi*pi*szp*ll2)/2d0 - (zeta3)/2d
     &0) + zp**5*((pi**2*5d0)/144d0 - 3d0/80d0 + (myi*pi**3*s
     &zp)/60d0 + (myi*pi*szp)/8d0 - (myi*pi*5d0*szp*ll2)/12d0
     & - (3d0*zeta3)/10d0) + zp**7*(-(187d0/3360d0) + (pi**2*
     &7d0)/240d0 + (myi*pi*721d0*szp)/5760d0 + (myi*pi**3*szp
     &)/84d0 - (myi*pi*7d0*szp*ll2)/20d0 - (3d0*zeta3)/14d0) 
     &+ zp**8*((pi**2*121d0)/4480d0 - 691d0/11520d0 + (myi*pi
     &*3931d0*szp)/32256d0 + (myi*pi**3*szp)/96d0 - (myi*pi*3
     &63d0*szp*ll2)/1120d0 - (3d0*zeta3)/16d0) + zp*((myi*pi*
     &*3*szp)/12d0 - (3d0*zeta3)/2d0) + zp**6*((pi**2*137d0)/
     &4320d0 - 7d0/144d0 + (myi*pi**3*szp)/72d0 + (myi*pi*731
     &d0*szp)/5760d0 - (myi*pi*137d0*szp*ll2)/360d0 - (zeta3)
     &/4d0) + zp**2*((pi**2)/24d0 + (myi*pi**3*szp)/24d0 - (m
     &yi*pi*szp*ll2)/2d0 - (3d0*zeta3)/4d0) + zp**9*(-(2521d0
     &/40320d0) + (pi**2*761d0)/30240d0 + (myi*pi**3*szp)/108
     &d0 + (myi*pi*42799d0*szp)/362880d0 - (myi*pi*761d0*szp*
     &ll2)/2520d0 - (zeta3)/6d0) + zp**4*((pi**2*11d0)/288d0 
     &- 1d0/48d0 + (myi*pi**3*szp)/48d0 + (myi*pi*11d0*szp)/9
     &6d0 - (myi*pi*11d0*szp*ll2)/24d0 - (3d0*zeta3)/8d0)

         case(45)            !0011

            zp = x+1d0

            ris = -((pi**4)/48d0) - (pi**2*ll2**2)/12d0
     & + (ll2**4)/12d0 + 2*cli4pt5 - (zp*zeta3)/8d0 + (7d0*ll
     &2*zeta3)/4d0 + zp**2*((ll2**2)/4d0 - (zeta3)/16d0) + zp
     &**3*(-((ll2)/12d0) + (ll2**2)/4d0 - (zeta3)/24d0) + zp*
     &*4*(1d0/96d0 - (11d0*ll2)/96d0 + (11d0*ll2**2)/48d0 - (
     &zeta3)/32d0) + zp**5*(17d0/960d0 - (ll2)/8d0 + (5d0*ll2
     &**2)/24d0 - (zeta3)/40d0) + zp**6*(253d0/11520d0 - (731
     &d0*ll2)/5760d0 + (137d0*ll2**2)/720d0 - (zeta3)/48d0) +
     & zp**7*(979d0/40320d0 - (721d0*ll2)/5760d0 + (7d0*ll2**
     &2)/40d0 - (zeta3)/56d0) + zp**8*(10943d0/430080d0 - (39
     &31d0*ll2)/32256d0 + (363d0*ll2**2)/2240d0 - (zeta3)/64d
     &0) + zp**9*(4703d0/181440d0 - (42799d0*ll2)/362880d0 + 
     &(761d0*ll2**2)/5040d0 - (zeta3)/72d0)

         case(46)            !01-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*11d0)/720d0 - (ll2**4)/8d0 - 3
     &*cli4pt5 + zp**2*(-(7d0/16d0) - (pi**2*ll2)/24d0 + (ll2
     &**3)/12d0 + (3d0*llzp)/8d0 - (llzp**2)/8d0 + (7d0*zeta3
     &)/16d0) + zp**3*(-(227d0/864d0) - (pi**2*ll2)/36d0 + (l
     &l2**3)/18d0 + (37d0*llzp)/144d0 - (5d0*llzp**2)/48d0 + 
     &(7d0*zeta3)/24d0) + zp**4*(-(1247d0/6912d0) - (pi**2*ll
     &2)/48d0 + (ll2**3)/24d0 + (107d0*llzp)/576d0 - (llzp**2
     &)/12d0 + (7d0*zeta3)/32d0) + zp**5*(-(470159d0/3456000d
     &0) - (pi**2*ll2)/60d0 + (ll2**3)/30d0 + (8257d0*llzp)/5
     &7600d0 - (131d0*llzp**2)/1920d0 + (7d0*zeta3)/40d0) + z
     &p**6*(-(2257309d0/20736000d0) - (pi**2*ll2)/72d0 + (ll2
     &**3)/36d0 + (13369d0*llzp)/115200d0 - (661d0*llzp**2)/1
     &1520d0 + (7d0*zeta3)/48d0) + zp**8*(-(183970943d0/23708
     &16000d0) - (pi**2*ll2)/96d0 + (ll2**3)/48d0 + (314543d0
     &*llzp)/3763200d0 - (1163d0*llzp**2)/26880d0 + (7d0*zeta
     &3)/64d0) + zp**9*(-(833624776009d0/12290310144000d0) - 
     &(pi**2*ll2)/108d0 + (ll2**3)/54d0 + (357205771d0*llzp)/
     &4877107200d0 - (148969d0*llzp**2)/3870720d0 + (7d0*zeta
     &3)/72d0) + zp**7*(-(107435801d0/1185408000d0) - (pi**2*
     &ll2)/84d0 + (ll2**3)/42d0 + (953d0*llzp)/9800d0 - (1327
     &d0*llzp**2)/26880d0 + (zeta3)/8d0) + zp*(-((pi**2*ll2)/
     &12d0) + (ll2**3)/6d0 + (7d0*zeta3)/8d0)

         case(47)            !01-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*13d0)/1440d0 - (myi*pi**3*szp*
     &ll2)/4d0 + (pi**2*ll2**2)/6d0 - (ll2**4)/6d0 - 4*cli4pt
     &5 + (myi*pi*13d0*szp*zeta3)/8d0 + zp*(-((myi*pi**3*szp)
     &/12d0) - (pi**2*ll2)/12d0 + (myi*pi*szp*ll2**2)/2d0 + z
     &eta3) + zp**2*(-((pi**2)/24d0) - (myi*pi**3*szp)/24d0 +
     & (myi*pi*3d0*szp)/8d0 - (pi**2*ll2)/24d0 + (myi*pi*szp*
     &ll2**2)/4d0 - (myi*pi*szp*llzp)/4d0 + (zeta3)/2d0) + zp
     &**3*(1d0/12d0 - (pi**2*5d0)/144d0 - (myi*pi**3*szp)/36d
     &0 + (myi*pi*37d0*szp)/144d0 - (pi**2*ll2)/36d0 + (myi*p
     &i*szp*ll2**2)/6d0 - (myi*pi*5d0*szp*llzp)/24d0 + (zeta3
     &)/3d0) + zp**4*(-((pi**2)/36d0) + 3d0/32d0 - (myi*pi**3
     &*szp)/48d0 + (myi*pi*107d0*szp)/576d0 - (pi**2*ll2)/48d
     &0 + (myi*pi*szp*ll2**2)/8d0 - (myi*pi*szp*llzp)/6d0 + (
     &zeta3)/4d0) + zp**5*(251d0/2880d0 - (pi**2*131d0)/5760d
     &0 - (myi*pi**3*szp)/60d0 + (myi*pi*8257d0*szp)/57600d0 
     &- (pi**2*ll2)/60d0 + (myi*pi*szp*ll2**2)/10d0 - (myi*pi
     &*131d0*szp*llzp)/960d0 + (zeta3)/5d0) + zp**6*(1343d0/1
     &7280d0 - (pi**2*661d0)/34560d0 + (myi*pi*13369d0*szp)/1
     &15200d0 - (myi*pi**3*szp)/72d0 - (pi**2*ll2)/72d0 + (my
     &i*pi*szp*ll2**2)/12d0 - (myi*pi*661d0*szp*llzp)/5760d0 
     &+ (zeta3)/6d0) + zp**7*(2977d0/43200d0 - (pi**2*1327d0)
     &/80640d0 - (myi*pi**3*szp)/84d0 + (myi*pi*953d0*szp)/98
     &00d0 - (pi**2*ll2)/84d0 + (myi*pi*szp*ll2**2)/14d0 - (m
     &yi*pi*1327d0*szp*llzp)/13440d0 + (zeta3)/7d0) + zp**8*(
     &29711d0/483840d0 - (pi**2*1163d0)/80640d0 + (myi*pi*314
     &543d0*szp)/3763200d0 - (myi*pi**3*szp)/96d0 - (pi**2*ll
     &2)/96d0 + (myi*pi*szp*ll2**2)/16d0 - (myi*pi*1163d0*szp
     &*llzp)/13440d0 + (zeta3)/8d0) + zp**9*(-((pi**2*148969d
     &0)/11612160d0) + 8406389d0/152409600d0 - (myi*pi**3*szp
     &)/108d0 + (myi*pi*357205771d0*szp)/4877107200d0 - (pi**
     &2*ll2)/108d0 + (myi*pi*szp*ll2**2)/18d0 - (myi*pi*14896
     &9d0*szp*llzp)/1935360d0 + (zeta3)/9d0)

         case(48)            !01-11

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*7d0)/720d0 + (pi**2*ll2**2)/4d
     &0 - (21d0*ll2*zeta3)/8d0 + zp**3*(-(1d0/24d0) + (pi**2*
     &5d0)/288d0 + (pi**2*ll2)/36d0 - (37d0*ll2)/144d0 - (5d0
     &*ll2**2)/48d0 - (ll2**3)/18d0 + (5d0*ll2*llzp)/24d0 - (
     &zeta3)/12d0) + zp**4*(-(17d0/384d0) + (pi**2)/72d0 + (p
     &i**2*ll2)/48d0 - (107d0*ll2)/576d0 - (ll2**2)/12d0 - (l
     &l2**3)/24d0 + (ll2*llzp)/6d0 - (zeta3)/16d0) + zp**5*((
     &pi**2*131d0)/11520d0 - 457d0/11520d0 + (pi**2*ll2)/60d0
     & - (8257d0*ll2)/57600d0 - (131d0*ll2**2)/1920d0 - (ll2*
     &*3)/30d0 + (131d0*ll2*llzp)/960d0 - (zeta3)/20d0) + zp*
     &*6*((pi**2*661d0)/69120d0 - 955d0/27648d0 - (13369d0*ll
     &2)/115200d0 + (pi**2*ll2)/72d0 - (661d0*ll2**2)/11520d0
     & - (ll2**3)/36d0 + (661d0*ll2*llzp)/5760d0 - (zeta3)/24
     &d0) + zp**7*((pi**2*1327d0)/161280d0 - 291769d0/9676800
     &d0 + (pi**2*ll2)/84d0 - (953d0*ll2)/9800d0 - (1327d0*ll
     &2**2)/26880d0 - (ll2**3)/42d0 + (1327d0*ll2*llzp)/13440
     &d0 - (zeta3)/28d0) + zp**8*((pi**2*1163d0)/161280d0 - 2
     &9407d0/1105920d0 - (314543d0*ll2)/3763200d0 + (pi**2*ll
     &2)/96d0 - (1163d0*ll2**2)/26880d0 - (ll2**3)/48d0 + (11
     &63d0*ll2*llzp)/13440d0 - (zeta3)/32d0) + zp**9*((pi**2*
     &148969d0)/23224320d0 - 231350923d0/9754214400d0 + (pi**
     &2*ll2)/108d0 - (357205771d0*ll2)/4877107200d0 - (148969
     &d0*ll2**2)/3870720d0 - (ll2**3)/54d0 + (148969d0*ll2*ll
     &zp)/1935360d0 - (zeta3)/36d0) + zp*((pi**2*ll2)/12d0 - 
     &(ll2**3)/6d0 - (zeta3)/4d0) + zp**2*((pi**2)/48d0 + (pi
     &**2*ll2)/24d0 - (3d0*ll2)/8d0 - (ll2**2)/8d0 - (ll2**3)
     &/12d0 + (ll2*llzp)/4d0 - (zeta3)/8d0)

         case(49)            !010-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/480d0 + zp**2*((pi**2)/24d0 -
     & (pi**2*ll2)/12d0 + (5d0*zeta3)/16d0) + zp**3*((pi**2*5
     &d0)/144d0 - 11d0/72d0 - (pi**2*ll2)/18d0 + (llzp)/12d0 
     &+ (5d0*zeta3)/24d0) + zp**4*((pi**2)/36d0 - 95d0/576d0 
     &- (pi**2*ll2)/24d0 + (5d0*llzp)/48d0 + (5d0*zeta3)/32d0
     &) + zp**6*((pi**2*661d0)/34560d0 - 941d0/7200d0 - (pi**
     &2*ll2)/36d0 + (47d0*llzp)/480d0 + (5d0*zeta3)/48d0) + z
     &p**7*(-(4d0/35d0) + (pi**2*1327d0)/80640d0 - (pi**2*ll2
     &)/42d0 + (13d0*llzp)/144d0 + (5d0*zeta3)/56d0) + zp**8*
     &(-(1137251d0/11289600d0) + (pi**2*1163d0)/80640d0 - (pi
     &**2*ll2)/48d0 + (3341d0*llzp)/40320d0 + (5d0*zeta3)/64d
     &0) + zp**9*((pi**2*148969d0)/11612160d0 - 20502379d0/22
     &8614400d0 - (pi**2*ll2)/54d0 + (13817d0*llzp)/181440d0 
     &+ (5d0*zeta3)/72d0) + zp**5*(-(43d0/288d0) + (pi**2*131
     &d0)/5760d0 - (pi**2*ll2)/30d0 + (5d0*llzp)/48d0 + (zeta
     &3)/8d0) + zp*(-((pi**2*ll2)/6d0) + (5d0*zeta3)/8d0)

         case(50)            !0100

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4)/80d0 + (myi*pi*3d0*szp*zeta3)
     &/2d0 + zp**9*((pi**2*148969d0)/3870720d0 - 37d0/768d0 +
     & (myi*pi*13817d0*szp)/181440d0 - (myi*pi**3*szp)/108d0 
     &- (pi**2*ll2)/18d0 + (zeta3)/12d0) + zp**4*((pi**2)/12d
     &0 - 1d0/48d0 - (myi*pi**3*szp)/48d0 + (myi*pi*5d0*szp)/
     &48d0 - (pi**2*ll2)/8d0 + (3d0*zeta3)/16d0) + zp**5*((pi
     &**2*131d0)/1920d0 - 17d0/480d0 + (myi*pi*5d0*szp)/48d0 
     &- (myi*pi**3*szp)/60d0 - (pi**2*ll2)/10d0 + (3d0*zeta3)
     &/20d0) + zp**7*((pi**2*1327d0)/26880d0 - 95d0/2016d0 + 
     &(myi*pi*13d0*szp)/144d0 - (myi*pi**3*szp)/84d0 - (pi**2
     &*ll2)/14d0 + (3d0*zeta3)/28d0) + zp**8*((pi**2*1163d0)/
     &26880d0 - 557d0/11520d0 + (myi*pi*3341d0*szp)/40320d0 -
     & (myi*pi**3*szp)/96d0 - (pi**2*ll2)/16d0 + (3d0*zeta3)/
     &32d0) + zp**3*((pi**2*5d0)/48d0 + (myi*pi*szp)/12d0 - (
     &myi*pi**3*szp)/36d0 - (pi**2*ll2)/6d0 + (zeta3)/4d0) + 
     &zp*(-((myi*pi**3*szp)/12d0) - (pi**2*ll2)/2d0 + (3d0*ze
     &ta3)/4d0) + zp**6*(-(25d0/576d0) + (pi**2*661d0)/11520d
     &0 + (myi*pi*47d0*szp)/480d0 - (myi*pi**3*szp)/72d0 - (p
     &i**2*ll2)/12d0 + (zeta3)/8d0) + zp**2*((pi**2)/8d0 - (m
     &yi*pi**3*szp)/24d0 - (pi**2*ll2)/4d0 + (3d0*zeta3)/8d0)

         case(51)            !0101

            zp = x+1d0

            ris = (pi**4*13d0)/288d0 + (pi**2*ll2**2)/6
     &d0 - (ll2**4)/6d0 - 4*cli4pt5 - (7d0*ll2*zeta3)/2d0 + z
     &p**3*((pi**2*5d0)/288d0 - (ll2)/12d0 - (pi**2*ll2)/36d0
     & + (zeta3)/12d0) + zp**4*((pi**2)/72d0 + 1d0/96d0 - (pi
     &**2*ll2)/48d0 - (5d0*ll2)/48d0 + (zeta3)/16d0) + zp**5*
     &((pi**2*131d0)/11520d0 + 1d0/60d0 - (5d0*ll2)/48d0 - (p
     &i**2*ll2)/60d0 + (zeta3)/20d0) + zp**6*((pi**2*661d0)/6
     &9120d0 + 7d0/360d0 - (47d0*ll2)/480d0 - (pi**2*ll2)/72d
     &0 + (zeta3)/24d0) + zp**7*((pi**2*1327d0)/161280d0 + 10
     &9d0/5376d0 - (13d0*ll2)/144d0 - (pi**2*ll2)/84d0 + (zet
     &a3)/28d0) + zp**8*((pi**2*1163d0)/161280d0 + 12979d0/64
     &5120d0 - (3341d0*ll2)/40320d0 - (pi**2*ll2)/96d0 + (zet
     &a3)/32d0) + zp**9*((pi**2*148969d0)/23224320d0 + 56591d
     &0/2903040d0 - (13817d0*ll2)/181440d0 - (pi**2*ll2)/108d
     &0 + (zeta3)/36d0) + zp*(-((pi**2*ll2)/12d0) + (zeta3)/4
     &d0) + zp**2*((pi**2)/48d0 - (pi**2*ll2)/24d0 + (zeta3)/
     &8d0)

         case(52)            !011-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/80d0 - (pi**2*ll2**2)/24d0 - 
     &(ll2**4)/12d0 - 2*cli4pt5 + zp**2*(-((pi**2)/48d0) + (l
     &l2**2)/8d0 - (ll2**3)/12d0 + (zeta3)/16d0) + zp**3*(11d
     &0/144d0 - (pi**2*5d0)/288d0 + (5d0*ll2**2)/48d0 - (ll2*
     &*3)/18d0 - (llzp)/24d0 + (zeta3)/24d0) + zp**4*(-((pi**
     &2)/72d0) + 59d0/768d0 + (ll2**2)/12d0 - (ll2**3)/24d0 -
     & (3d0*llzp)/64d0 + (zeta3)/32d0) + zp**5*(-((pi**2*131d
     &0)/11520d0) + 7651d0/115200d0 + (131d0*ll2**2)/1920d0 -
     & (ll2**3)/30d0 - (83d0*llzp)/1920d0 + (zeta3)/40d0) + z
     &p**6*(65d0/1152d0 - (pi**2*661d0)/69120d0 + (661d0*ll2*
     &*2)/11520d0 - (ll2**3)/36d0 - (11d0*llzp)/288d0 + (zeta
     &3)/48d0) + zp**7*(-((pi**2*1327d0)/161280d0) + 1092631d
     &0/22579200d0 + (1327d0*ll2**2)/26880d0 - (ll2**3)/42d0 
     &- (5417d0*llzp)/161280d0 + (zeta3)/56d0) + zp**8*(-((pi
     &**2*1163d0)/161280d0) + 7763d0/184320d0 + (1163d0*ll2**
     &2)/26880d0 - (ll2**3)/48d0 - (137d0*llzp)/4608d0 + (zet
     &a3)/64d0) + zp**9*(-((pi**2*148969d0)/23224320d0) + 217
     &6291643d0/58525286400d0 + (148969d0*ll2**2)/3870720d0 -
     & (ll2**3)/54d0 - (617027d0*llzp)/23224320d0 + (zeta3)/7
     &2d0) + zp*(-((ll2**3)/6d0) + (zeta3)/8d0)

         case(53)            !0110

            zp = x+1d0
            szp = s(zp)

            ris = -((pi**4)/288d0) + (myi*pi*szp*zeta3)
     &/8d0 + zp**2*(-((pi**2)/48d0) + (pi**2*ll2)/24d0 + (myi
     &*pi*szp*ll2)/4d0 - (myi*pi*szp*ll2**2)/4d0 - (zeta3)/16
     &d0) + zp**3*(-((pi**2*5d0)/288d0) - (myi*pi*szp)/24d0 +
     & (pi**2*ll2)/36d0 + (myi*pi*5d0*szp*ll2)/24d0 - (myi*pi
     &*szp*ll2**2)/6d0 - (zeta3)/24d0) + zp**4*(-((pi**2)/72d
     &0) + 1d0/96d0 - (myi*pi*3d0*szp)/64d0 + (pi**2*ll2)/48d
     &0 + (myi*pi*szp*ll2)/6d0 - (myi*pi*szp*ll2**2)/8d0 - (z
     &eta3)/32d0) + zp**5*(-((pi**2*131d0)/11520d0) + 1d0/64d
     &0 - (myi*pi*83d0*szp)/1920d0 + (pi**2*ll2)/60d0 + (myi*
     &pi*131d0*szp*ll2)/960d0 - (myi*pi*szp*ll2**2)/10d0 - (z
     &eta3)/40d0) + zp**6*(11d0/640d0 - (pi**2*661d0)/69120d0
     & - (myi*pi*11d0*szp)/288d0 + (pi**2*ll2)/72d0 + (myi*pi
     &*661d0*szp*ll2)/5760d0 - (myi*pi*szp*ll2**2)/12d0 - (ze
     &ta3)/48d0) + zp**7*(-((pi**2*1327d0)/161280d0) + 49d0/2
     &880d0 - (myi*pi*5417d0*szp)/161280d0 + (pi**2*ll2)/84d0
     & + (myi*pi*1327d0*szp*ll2)/13440d0 - (myi*pi*szp*ll2**2
     &)/14d0 - (zeta3)/56d0) + zp**8*(-((pi**2*1163d0)/161280
     &d0) + 2603d0/161280d0 - (myi*pi*137d0*szp)/4608d0 + (pi
     &**2*ll2)/96d0 + (myi*pi*1163d0*szp*ll2)/13440d0 - (myi*
     &pi*szp*ll2**2)/16d0 - (zeta3)/64d0) + zp**9*(-((pi**2*1
     &48969d0)/23224320d0) + 809d0/53760d0 - (myi*pi*617027d0
     &*szp)/23224320d0 + (pi**2*ll2)/108d0 + (myi*pi*148969d0
     &*szp*ll2)/1935360d0 - (myi*pi*szp*ll2**2)/18d0 - (zeta3
     &)/72d0) + zp*((pi**2*ll2)/12d0 - (myi*pi*szp*ll2**2)/2d
     &0 - (zeta3)/8d0)

         case(54)            !0111

            zp = x+1d0

            ris = -((pi**4)/90d0) - (pi**2*ll2**2)/24d0
     & + (zp*ll2**3)/6d0 + (ll2**4)/24d0 + zp**2*(-((ll2**2)/
     &8d0) + (ll2**3)/12d0) + zp**3*((ll2)/24d0 - (5d0*ll2**2
     &)/48d0 + (ll2**3)/18d0) + zp**4*(-(1d0/192d0) + (3d0*ll
     &2)/64d0 - (ll2**2)/12d0 + (ll2**3)/24d0) + zp**5*(-(7d0
     &/960d0) + (83d0*ll2)/1920d0 - (131d0*ll2**2)/1920d0 + (
     &ll2**3)/30d0) + zp**6*(-(35d0/4608d0) + (11d0*ll2)/288d
     &0 - (661d0*ll2**2)/11520d0 + (ll2**3)/36d0) + zp**7*(-(
     &155d0/21504d0) + (5417d0*ll2)/161280d0 - (1327d0*ll2**2
     &)/26880d0 + (ll2**3)/42d0) + zp**8*(-(2441d0/368640d0) 
     &+ (137d0*ll2)/4608d0 - (1163d0*ll2**2)/26880d0 + (ll2**
     &3)/48d0) + zp**9*(-(19997d0/3317760d0) + (617027d0*ll2)
     &/23224320d0 - (148969d0*ll2**2)/3870720d0 + (ll2**3)/54
     &d0) + cli4pt5 + (7d0*ll2*zeta3)/8d0

         case(55)            !1-1-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = zp**8*(-(1d0/1048576d0) + (llzp)/1310
     &72d0 - (llzp**2)/32768d0 + (llzp**3)/12288d0) + zp*(-(1
     &d0/2d0) + (llzp)/2d0 - (llzp**2)/4d0 + (llzp**3)/12d0) 
     &+ zp**3*(-(1d0/648d0) + (llzp)/216d0 - (llzp**2)/144d0 
     &+ (llzp**3)/144d0) + zp**6*(-(1d0/82944d0) + (llzp)/138
     &24d0 - (llzp**2)/4608d0 + (llzp**3)/2304d0) + zp**9*(-(
     &1d0/3359232d0) + (llzp)/373248d0 - (llzp**2)/82944d0 + 
     &(llzp**3)/27648d0) + zp**4*(-(1d0/4096d0) + (llzp)/1024
     &d0 - (llzp**2)/512d0 + (llzp**3)/384d0) + zp**2*(-(1d0/
     &64d0) + (llzp)/32d0 - (llzp**2)/32d0 + (llzp**3)/48d0) 
     &+ zp**7*(-(1d0/307328d0) + (llzp)/43904d0 - (llzp**2)/1
     &2544d0 + (llzp**3)/5376d0) + zp**5*(-(1d0/20000d0) + (l
     &lzp)/4000d0 - (llzp**2)/1600d0 + (llzp**3)/960d0) + cli
     &4pt5

         case(56)            !1-1-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4)/720d0 + (myi*pi**3*szp*ll2)/1
     &2d0 - (myi*pi*szp*ll2**3)/6d0 + (ll2**4)/24d0 + cli4pt5
     & - (myi*pi*7d0*szp*zeta3)/8d0 - (ll2*zeta3)/8d0 + zp**5
     &*(-(329d0/17280d0) - (pi**2)/4800d0 + (myi*pi*szp)/4000
     &d0 + (pi**2*llzp)/960d0 - (myi*pi*szp*llzp)/800d0 + (my
     &i*pi*szp*llzp**2)/320d0 + (zeta3)/160d0) + zp**8*(-(183
     &91283d0/9483264000d0) - (pi**2)/98304d0 + (myi*pi*szp)/
     &131072d0 + (pi**2*llzp)/12288d0 - (myi*pi*szp*llzp)/163
     &84d0 + (myi*pi*szp*llzp**2)/4096d0 + (zeta3)/2048d0) + 
     &zp**3*(-((pi**2)/432d0) - 5d0/48d0 + (myi*pi*szp)/216d0
     & + (pi**2*llzp)/144d0 - (myi*pi*szp*llzp)/72d0 + (myi*p
     &i*szp*llzp**2)/48d0 + (zeta3)/24d0) + zp*(-((pi**2)/12d
     &0) + (myi*pi*szp)/2d0 + (pi**2*llzp)/12d0 - (myi*pi*szp
     &*llzp)/2d0 + (myi*pi*szp*llzp**2)/4d0 + (zeta3)/2d0) + 
     &zp**6*(-((pi**2)/13824d0) - 44581d0/5184000d0 + (myi*pi
     &*szp)/13824d0 + (pi**2*llzp)/2304d0 - (myi*pi*szp*llzp)
     &/2304d0 + (myi*pi*szp*llzp**2)/768d0 + (zeta3)/384d0) +
     & zp**9*(-(20706533d0/21337344000d0) - (pi**2)/248832d0 
     &+ (myi*pi*szp)/373248d0 + (pi**2*llzp)/27648d0 - (myi*p
     &i*szp*llzp)/41472d0 + (myi*pi*szp*llzp**2)/9216d0 + (ze
     &ta3)/4608d0) + zp**4*(-((pi**2)/1536d0) - 151d0/3456d0 
     &+ (myi*pi*szp)/1024d0 + (pi**2*llzp)/384d0 - (myi*pi*sz
     &p*llzp)/256d0 + (myi*pi*szp*llzp**2)/128d0 + (zeta3)/64
     &d0) + zp**7*(-((pi**2)/37632d0) - 48581d0/12096000d0 + 
     &(myi*pi*szp)/43904d0 + (pi**2*llzp)/5376d0 - (myi*pi*sz
     &p*llzp)/6272d0 + (myi*pi*szp*llzp**2)/1792d0 + (zeta3)/
     &896d0) + zp**2*(-(1d0/4d0) - (pi**2)/96d0 + (myi*pi*szp
     &)/32d0 + (pi**2*llzp)/48d0 - (myi*pi*szp*llzp)/16d0 + (
     &myi*pi*szp*llzp**2)/16d0 + (zeta3)/8d0)

         case(57)            !1-1-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/288d0) - (pi**2*ll2**2)/24d
     &0 + (ll2**4)/24d0 + (7d0*ll2*zeta3)/8d0 + zp**7*(4081d0
     &/3072000d0 + (pi**2)/75264d0 + (pi**2*ll2)/10752d0 - (l
     &l2)/43904d0 - (ll2**2)/12544d0 - (ll2**3)/5376d0 - (pi*
     &*2*llzp)/10752d0 + (ll2*llzp)/6272d0 + (ll2**2*llzp)/17
     &92d0 - (ll2*llzp**2)/1792d0 - (zeta3)/1024d0) + zp**5*(
     &407d0/55296d0 + (pi**2)/9600d0 + (pi**2*ll2)/1920d0 - (
     &ll2)/4000d0 - (ll2**2)/1600d0 - (ll2**3)/960d0 - (pi**2
     &*llzp)/1920d0 + (ll2*llzp)/800d0 + (ll2**2*llzp)/320d0 
     &- (ll2*llzp**2)/320d0 - (7d0*zeta3)/1280d0) + zp**8*((p
     &i**2)/196608d0 + 9822481d0/16859136000d0 - (ll2)/131072
     &d0 + (pi**2*ll2)/24576d0 - (ll2**2)/32768d0 - (ll2**3)/
     &12288d0 - (pi**2*llzp)/24576d0 + (ll2*llzp)/16384d0 + (
     &ll2**2*llzp)/4096d0 - (ll2*llzp**2)/4096d0 - (7d0*zeta3
     &)/16384d0) + zp*((pi**2)/24d0 + (pi**2*ll2)/24d0 - (ll2
     &)/2d0 - (ll2**2)/4d0 - (ll2**3)/12d0 - (pi**2*llzp)/24d
     &0 + (ll2*llzp)/2d0 + (ll2**2*llzp)/4d0 - (ll2*llzp**2)/
     &4d0 - (7d0*zeta3)/16d0) + zp**3*(3d0/64d0 + (pi**2)/864
     &d0 - (ll2)/216d0 + (pi**2*ll2)/288d0 - (ll2**2)/144d0 -
     & (ll2**3)/144d0 - (pi**2*llzp)/288d0 + (ll2*llzp)/72d0 
     &+ (ll2**2*llzp)/48d0 - (ll2*llzp**2)/48d0 - (7d0*zeta3)
     &/192d0) + zp**6*((pi**2)/27648d0 + 256103d0/82944000d0 
     &- (ll2)/13824d0 + (pi**2*ll2)/4608d0 - (ll2**2)/4608d0 
     &- (ll2**3)/2304d0 - (pi**2*llzp)/4608d0 + (ll2*llzp)/23
     &04d0 + (ll2**2*llzp)/768d0 - (ll2*llzp**2)/768d0 - (7d0
     &*zeta3)/3072d0) + zp**9*((pi**2)/497664d0 + 78708473d0/
     &303464448000d0 - (ll2)/373248d0 + (pi**2*ll2)/55296d0 -
     & (ll2**2)/82944d0 - (ll2**3)/27648d0 - (pi**2*llzp)/552
     &96d0 + (ll2*llzp)/41472d0 + (ll2**2*llzp)/9216d0 - (ll2
     &*llzp**2)/9216d0 - (7d0*zeta3)/36864d0) + zp**4*(251d0/
     &13824d0 + (pi**2)/3072d0 - (ll2)/1024d0 + (pi**2*ll2)/7
     &68d0 - (ll2**2)/512d0 - (ll2**3)/384d0 - (pi**2*llzp)/7
     &68d0 + (ll2*llzp)/256d0 + (ll2**2*llzp)/128d0 - (ll2*ll
     &zp**2)/128d0 - (7d0*zeta3)/512d0) + zp**2*((pi**2)/192d
     &0 + 1d0/8d0 - (ll2)/32d0 + (pi**2*ll2)/96d0 - (ll2**2)/
     &32d0 - (ll2**3)/48d0 - (pi**2*llzp)/96d0 + (ll2*llzp)/1
     &6d0 + (ll2**2*llzp)/16d0 - (ll2*llzp**2)/16d0 - (7d0*ze
     &ta3)/64d0)

         case(58)            !1-10-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*19d0)/1440d0) + (pi**2*ll2**
     &2)/24d0 + (ll2**4)/24d0 + cli4pt5 + zp*((pi**2)/12d0 - 
     &(pi**2*llzp)/12d0 - zeta3) + (ll2*zeta3)/4d0 + zp**8*(1
     &0723549d0/2370816000d0 + (pi**2)/98304d0 - (pi**2*llzp)
     &/12288d0 - (9701d0*llzp)/1881600d0 - (zeta3)/1024d0) + 
     &zp**3*((pi**2)/432d0 + 1d0/4d0 - (pi**2*llzp)/144d0 - (
     &llzp)/8d0 - (zeta3)/12d0) + zp**6*((pi**2)/13824d0 + 51
     &521d0/2592000d0 - (pi**2*llzp)/2304d0 - (347d0*llzp)/21
     &600d0 - (zeta3)/192d0) + zp**9*(24451813d0/10668672000d
     &0 + (pi**2)/248832d0 - (pi**2*llzp)/27648d0 - (209d0*ll
     &zp)/66150d0 - (zeta3)/2304d0) + zp**4*((pi**2)/1536d0 +
     & 709d0/6912d0 - (pi**2*llzp)/384d0 - (35d0*llzp)/576d0 
     &- (zeta3)/32d0) + zp**7*((pi**2)/37632d0 + 393707d0/423
     &36000d0 - (149d0*llzp)/16800d0 - (pi**2*llzp)/5376d0 - 
     &(zeta3)/448d0) + zp**2*(5d0/8d0 + (pi**2)/96d0 - (pi**2
     &*llzp)/48d0 - (llzp)/4d0 - (zeta3)/4d0) + zp**5*(1909d0
     &/43200d0 + (pi**2)/4800d0 - (11d0*llzp)/360d0 - (pi**2*
     &llzp)/960d0 - (zeta3)/80d0)

         case(59)            !1-100

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*41d0)/1440d0) + (myi*pi**3*s
     &zp*ll2)/12d0 + (pi**2*ll2**2)/4d0 - myi*pi*szp*zeta3 - 
     &(3d0*ll2*zeta3)/4d0 + zp**5*((pi**2)/1600d0 + 5d0/192d0
     & - (myi*pi*11d0*szp)/360d0 + (myi*pi**3*szp)/960d0 - (p
     &i**2*llzp)/320d0 - (zeta3)/160d0) + zp**8*((pi**2)/3276
     &8d0 + 4669d0/552960d0 + (myi*pi**3*szp)/12288d0 - (myi*
     &pi*9701d0*szp)/1881600d0 - (pi**2*llzp)/4096d0 - (zeta3
     &)/2048d0) + zp**3*((pi**2)/144d0 + 1d0/24d0 + (myi*pi**
     &3*szp)/144d0 - (myi*pi*szp)/8d0 - (pi**2*llzp)/48d0 - (
     &zeta3)/24d0) + zp*((pi**2)/4d0 + (myi*pi**3*szp)/12d0 -
     & (pi**2*llzp)/4d0 - (zeta3)/2d0) + zp**6*(41d0/2304d0 +
     & (pi**2)/4608d0 + (myi*pi**3*szp)/2304d0 - (myi*pi*347d
     &0*szp)/21600d0 - (pi**2*llzp)/768d0 - (zeta3)/384d0) + 
     &zp**9*(10457d0/1741824d0 + (pi**2)/82944d0 + (myi*pi**3
     &*szp)/27648d0 - (myi*pi*209d0*szp)/66150d0 - (pi**2*llz
     &p)/9216d0 - (zeta3)/4608d0) + zp**4*((pi**2)/512d0 + 7d
     &0/192d0 + (myi*pi**3*szp)/384d0 - (myi*pi*35d0*szp)/576
     &d0 - (pi**2*llzp)/128d0 - (zeta3)/64d0) + zp**7*((pi**2
     &)/12544d0 + 2941d0/241920d0 - (myi*pi*149d0*szp)/16800d
     &0 + (myi*pi**3*szp)/5376d0 - (pi**2*llzp)/1792d0 - (zet
     &a3)/896d0) + zp**2*((pi**2)/32d0 + (myi*pi**3*szp)/48d0
     & - (myi*pi*szp)/4d0 - (pi**2*llzp)/16d0 - (zeta3)/8d0)

         case(60)            !1-101

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*11d0)/240d0) - (pi**2*ll2**2
     &)/8d0 + (ll2**4)/6d0 + 4*cli4pt5 + (13d0*ll2*zeta3)/4d0
     & + zp**5*(-(31d0/2880d0) + (pi**2)/9600d0 + (11d0*ll2)/
     &360d0 - (pi**2*llzp)/1920d0 - (zeta3)/256d0) + zp**8*((
     &pi**2)/196608d0 - 744197d0/270950400d0 + (9701d0*ll2)/1
     &881600d0 - (pi**2*llzp)/24576d0 - (5d0*zeta3)/16384d0) 
     &+ zp*((pi**2)/24d0 - (pi**2*llzp)/24d0 - (5d0*zeta3)/16
     &d0) + zp**3*(-(1d0/48d0) + (pi**2)/864d0 + (ll2)/8d0 - 
     &(pi**2*llzp)/288d0 - (5d0*zeta3)/192d0) + zp**6*((pi**2
     &)/27648d0 - 73d0/10800d0 + (347d0*ll2)/21600d0 - (pi**2
     &*llzp)/4608d0 - (5d0*zeta3)/3072d0) + zp**9*((pi**2)/49
     &7664d0 - 555271d0/304819200d0 + (209d0*ll2)/66150d0 - (
     &pi**2*llzp)/55296d0 - (5d0*zeta3)/36864d0) + zp**4*(-(1
     &9d0/1152d0) + (pi**2)/3072d0 + (35d0*ll2)/576d0 - (pi**
     &2*llzp)/768d0 - (5d0*zeta3)/512d0) + zp**2*((pi**2)/192
     &d0 + (ll2)/4d0 - (pi**2*llzp)/96d0 - (5d0*zeta3)/64d0) 
     &+ zp**7*(-(10313d0/2419200d0) + (pi**2)/75264d0 + (149d
     &0*ll2)/16800d0 - (pi**2*llzp)/10752d0 - (5d0*zeta3)/716
     &8d0)

         case(61)            !1-11-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/1440d0 - (pi**2*ll2**2)/24d0 
     &+ (ll2**4)/24d0 + (ll2*zeta3)/4d0 + zp**7*(-((pi**2)/75
     &264d0) - 93371d0/32256000d0 - (pi**2*ll2)/5376d0 + (ll2
     &**2)/12544d0 + (ll2**3)/2688d0 + (pi**2*llzp)/10752d0 +
     & (767d0*llzp)/460800d0 - (ll2**2*llzp)/1792d0 + (zeta3)
     &/512d0) + zp**6*(-(17653d0/2592000d0) - (pi**2)/27648d0
     & - (pi**2*ll2)/2304d0 + (ll2**2)/4608d0 + (ll2**3)/1152
     &d0 + (pi**2*llzp)/4608d0 + (5269d0*llzp)/1382400d0 - (l
     &l2**2*llzp)/768d0 + (7d0*zeta3)/1536d0) + zp**9*(-(2276
     &013631d0/4096770048000d0) - (pi**2)/497664d0 - (pi**2*l
     &l2)/27648d0 + (ll2**2)/82944d0 + (ll2**3)/13824d0 + (10
     &77749d0*llzp)/3251404800d0 + (pi**2*llzp)/55296d0 - (ll
     &2**2*llzp)/9216d0 + (7d0*zeta3)/18432d0) + zp**4*(-(115
     &1d0/27648d0) - (pi**2)/3072d0 - (pi**2*ll2)/384d0 + (ll
     &2**2)/512d0 + (ll2**3)/192d0 + (49d0*llzp)/2304d0 + (pi
     &**2*llzp)/768d0 - (ll2**2*llzp)/128d0 + (7d0*zeta3)/256
     &d0) + zp**2*(-((pi**2)/192d0) - 5d0/16d0 - (pi**2*ll2)/
     &48d0 + (ll2**2)/32d0 + (ll2**3)/24d0 + (llzp)/8d0 + (pi
     &**2*llzp)/96d0 - (ll2**2*llzp)/16d0 + (7d0*zeta3)/32d0)
     & + zp**5*(-(2281d0/138240d0) - (pi**2)/9600d0 - (pi**2*
     &ll2)/960d0 + (ll2**2)/1600d0 + (ll2**3)/480d0 + (pi**2*
     &llzp)/1920d0 + (41d0*llzp)/4608d0 - (ll2**2*llzp)/320d0
     & + (7d0*zeta3)/640d0) + zp**8*(-(127203607d0/1011548160
     &00d0) - (pi**2)/196608d0 - (pi**2*ll2)/12288d0 + (ll2**
     &2)/32768d0 + (ll2**3)/6144d0 + (pi**2*llzp)/24576d0 + (
     &266681d0*llzp)/361267200d0 - (ll2**2*llzp)/4096d0 + (7d
     &0*zeta3)/8192d0) + zp*(-((pi**2)/24d0) - (pi**2*ll2)/12
     &d0 + (ll2**2)/4d0 + (ll2**3)/6d0 + (pi**2*llzp)/24d0 - 
     &(ll2**2*llzp)/4d0 + (7d0*zeta3)/8d0) + zp**3*(-((pi**2)
     &/864d0) - 1d0/9d0 - (pi**2*ll2)/144d0 + (ll2**2)/144d0 
     &+ (ll2**3)/72d0 + (pi**2*llzp)/288d0 + (5d0*llzp)/96d0 
     &- (ll2**2*llzp)/48d0 + (7d0*zeta3)/96d0)

         case(62)            !1-110

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4)/360d0) - (myi*pi**3*szp*ll2
     &)/12d0 - (pi**2*ll2**2)/24d0 + (myi*pi*szp*ll2**3)/6d0 
     &+ (myi*pi*szp*zeta3)/4d0 + ll2*zeta3 + zp**5*(-(49d0/57
     &60d0) - (pi**2)/9600d0 - (myi*pi**3*szp)/1920d0 + (myi*
     &pi*41d0*szp)/4608d0 - (pi**2*ll2)/640d0 + (myi*pi*szp*l
     &l2)/800d0 + (myi*pi*szp*ll2**2)/320d0 + (pi**2*llzp)/19
     &20d0 - (myi*pi*szp*ll2*llzp)/160d0 + (13d0*zeta3)/1280d
     &0) + zp**8*(-((pi**2)/196608d0) - 374123d0/270950400d0 
     &- (myi*pi**3*szp)/24576d0 + (myi*pi*266681d0*szp)/36126
     &7200d0 - (pi**2*ll2)/8192d0 + (myi*pi*szp*ll2)/16384d0 
     &+ (myi*pi*szp*ll2**2)/4096d0 + (pi**2*llzp)/24576d0 - (
     &myi*pi*szp*ll2*llzp)/2048d0 + (13d0*zeta3)/16384d0) + z
     &p*(-((pi**2)/24d0) - (myi*pi**3*szp)/24d0 - (pi**2*ll2)
     &/8d0 + (myi*pi*szp*ll2)/2d0 + (myi*pi*szp*ll2**2)/4d0 +
     & (pi**2*llzp)/24d0 - (myi*pi*szp*ll2*llzp)/2d0 + (13d0*
     &zeta3)/16d0) + zp**3*(-(1d0/48d0) - (pi**2)/864d0 - (my
     &i*pi**3*szp)/288d0 + (myi*pi*5d0*szp)/96d0 - (pi**2*ll2
     &)/96d0 + (myi*pi*szp*ll2)/72d0 + (myi*pi*szp*ll2**2)/48
     &d0 + (pi**2*llzp)/288d0 - (myi*pi*szp*ll2*llzp)/24d0 + 
     &(13d0*zeta3)/192d0) + zp**6*(-((pi**2)/27648d0) - 1609d
     &0/345600d0 - (myi*pi**3*szp)/4608d0 + (myi*pi*5269d0*sz
     &p)/1382400d0 - (pi**2*ll2)/1536d0 + (myi*pi*szp*ll2)/23
     &04d0 + (myi*pi*szp*ll2**2)/768d0 + (pi**2*llzp)/4608d0 
     &- (myi*pi*szp*ll2*llzp)/384d0 + (13d0*zeta3)/3072d0) + 
     &zp**9*(-((pi**2)/497664d0) - 469253d0/609638400d0 + (my
     &i*pi*1077749d0*szp)/3251404800d0 - (myi*pi**3*szp)/5529
     &6d0 - (pi**2*ll2)/18432d0 + (myi*pi*szp*ll2)/41472d0 + 
     &(myi*pi*szp*ll2**2)/9216d0 + (pi**2*llzp)/55296d0 - (my
     &i*pi*szp*ll2*llzp)/4608d0 + (13d0*zeta3)/36864d0) + zp*
     &*4*(-(17d0/1152d0) - (pi**2)/3072d0 + (myi*pi*49d0*szp)
     &/2304d0 - (myi*pi**3*szp)/768d0 - (pi**2*ll2)/256d0 + (
     &myi*pi*szp*ll2)/256d0 + (myi*pi*szp*ll2**2)/128d0 + (pi
     &**2*llzp)/768d0 - (myi*pi*szp*ll2*llzp)/64d0 + (13d0*ze
     &ta3)/512d0) + zp**2*(-((pi**2)/192d0) + (myi*pi*szp)/8d
     &0 - (myi*pi**3*szp)/96d0 - (pi**2*ll2)/32d0 + (myi*pi*s
     &zp*ll2)/16d0 + (myi*pi*szp*ll2**2)/16d0 + (pi**2*llzp)/
     &96d0 - (myi*pi*szp*ll2*llzp)/8d0 + (13d0*zeta3)/64d0) +
     & zp**7*(-(6107d0/2419200d0) - (pi**2)/75264d0 - (myi*pi
     &**3*szp)/10752d0 + (myi*pi*767d0*szp)/460800d0 - (pi**2
     &*ll2)/3584d0 + (myi*pi*szp*ll2)/6272d0 + (myi*pi*szp*ll
     &2**2)/1792d0 + (pi**2*llzp)/10752d0 - (myi*pi*szp*ll2*l
     &lzp)/896d0 + (13d0*zeta3)/7168d0)

         case(63)            !1-111

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/30d0 + (pi**2*ll2**2)/6d0 - (
     &ll2**4)/6d0 - 3*cli4pt5 - (23d0*ll2*zeta3)/8d0 + zp**5*
     &(17d0/5120d0 + (pi**2*ll2)/1920d0 - (41d0*ll2)/4608d0 -
     & (ll2**2)/1600d0 - (ll2**3)/480d0 + (ll2**2*llzp)/320d0
     & - (zeta3)/1280d0) + zp**8*(629d0/1769472d0 + (pi**2*ll
     &2)/24576d0 - (266681d0*ll2)/361267200d0 - (ll2**2)/3276
     &8d0 - (ll2**3)/6144d0 + (ll2**2*llzp)/4096d0 - (zeta3)/
     &16384d0) + zp*((pi**2*ll2)/24d0 - (ll2**2)/4d0 - (ll2**
     &3)/6d0 + (ll2**2*llzp)/4d0 - (zeta3)/16d0) + zp**3*(1d0
     &/96d0 + (pi**2*ll2)/288d0 - (5d0*ll2)/96d0 - (ll2**2)/1
     &44d0 - (ll2**3)/72d0 + (ll2**2*llzp)/48d0 - (zeta3)/192
     &d0) + zp**6*(59d0/36864d0 + (pi**2*ll2)/4608d0 - (5269d
     &0*ll2)/1382400d0 - (ll2**2)/4608d0 - (ll2**3)/1152d0 + 
     &(ll2**2*llzp)/768d0 - (zeta3)/3072d0) + zp**9*(185921d0
     &/1114767360d0 - (1077749d0*ll2)/3251404800d0 + (pi**2*l
     &l2)/55296d0 - (ll2**2)/82944d0 - (ll2**3)/13824d0 + (ll
     &2**2*llzp)/9216d0 - (zeta3)/36864d0) + zp**4*(5d0/768d0
     & - (49d0*ll2)/2304d0 + (pi**2*ll2)/768d0 - (ll2**2)/512
     &d0 - (ll2**3)/192d0 + (ll2**2*llzp)/128d0 - (zeta3)/512
     &d0) + zp**2*(-((ll2)/8d0) + (pi**2*ll2)/96d0 - (ll2**2)
     &/32d0 - (ll2**3)/24d0 + (ll2**2*llzp)/16d0 - (zeta3)/64
     &d0) + zp**7*(2929d0/3870720d0 + (pi**2*ll2)/10752d0 - (
     &767d0*ll2)/460800d0 - (ll2**2)/12544d0 - (ll2**3)/2688d
     &0 + (ll2**2*llzp)/1792d0 - (zeta3)/7168d0)

         case(64)            !10-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/288d0) - (pi**2*ll2**2)/24d
     &0 + (ll2**4)/24d0 + cli4pt5 + (zp*zeta3)/2d0 - (ll2*zet
     &a3)/8d0 + zp**5*(-(12017d0/432000d0) + (79d0*llzp)/1800
     &d0 - (llzp**2)/30d0 + (zeta3)/160d0) + zp**8*(-(2783246
     &3d0/9483264000d0) + (7493d0*llzp)/940800d0 - (151d0*llz
     &p**2)/13440d0 + (zeta3)/2048d0) + zp**3*(-(71d0/432d0) 
     &+ (13d0*llzp)/72d0 - (llzp**2)/12d0 + (zeta3)/24d0) + z
     &p**6*(-(64861d0/5184000d0) + (169d0*llzp)/7200d0 - (llz
     &p**2)/45d0 + (zeta3)/384d0) + zp**9*(-(293914637d0/1920
     &36096000d0) + (3001d0*llzp)/595350d0 - (8d0*llzp**2)/94
     &5d0 + (zeta3)/4608d0) + zp**4*(-(113d0/1728d0) + (25d0*
     &llzp)/288d0 - (5d0*llzp**2)/96d0 + (zeta3)/64d0) + zp**
     &7*(-(3505829d0/592704000d0) + (521d0*llzp)/39200d0 - (1
     &3d0*llzp**2)/840d0 + (zeta3)/896d0) + zp**2*(-(7d0/16d0
     &) + (3d0*llzp)/8d0 - (llzp**2)/8d0 + (zeta3)/8d0)

         case(65)            !10-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = -((pi**4*17d0)/480d0) + (myi*pi**3*sz
     &p*ll2)/6d0 - (pi**2*ll2**2)/6d0 + (ll2**4)/6d0 + 4*cli4
     &pt5 - (myi*pi*5d0*szp*zeta3)/8d0 + (3d0*ll2*zeta3)/2d0 
     &+ zp*(-((myi*pi**3*szp)/12d0) + zeta3) + zp**8*(-((pi**
     &2*151d0)/40320d0) + 42371d0/1382400d0 - (myi*pi**3*szp)
     &/12288d0 + (myi*pi*7493d0*szp)/940800d0 - (myi*pi*151d0
     &*szp*llzp)/6720d0 + (zeta3)/1024d0) + zp**3*(1d0/12d0 -
     & (pi**2)/36d0 - (myi*pi**3*szp)/144d0 + (myi*pi*13d0*sz
     &p)/72d0 - (myi*pi*szp*llzp)/6d0 + (zeta3)/12d0) + zp**6
     &*(-((pi**2)/135d0) + 179d0/3456d0 - (myi*pi**3*szp)/230
     &4d0 + (myi*pi*169d0*szp)/7200d0 - (myi*pi*2d0*szp*llzp)
     &/45d0 + (zeta3)/192d0) + zp**9*(735253d0/30481920d0 - (
     &pi**2*8d0)/2835d0 - (myi*pi**3*szp)/27648d0 + (myi*pi*3
     &001d0*szp)/595350d0 - (myi*pi*16d0*szp*llzp)/945d0 + (z
     &eta3)/2304d0) + zp**4*(1d0/12d0 - (pi**2*5d0)/288d0 + (
     &myi*pi*25d0*szp)/288d0 - (myi*pi**3*szp)/384d0 - (myi*p
     &i*5d0*szp*llzp)/48d0 + (zeta3)/32d0) + zp**7*(-((pi**2*
     &13d0)/2520d0) + 23963d0/604800d0 + (myi*pi*521d0*szp)/3
     &9200d0 - (myi*pi**3*szp)/5376d0 - (myi*pi*13d0*szp*llzp
     &)/420d0 + (zeta3)/448d0) + zp**2*(-((pi**2)/24d0) - (my
     &i*pi**3*szp)/48d0 + (myi*pi*3d0*szp)/8d0 - (myi*pi*szp*
     &llzp)/4d0 + (zeta3)/4d0) + zp**5*(-((pi**2)/90d0) + 97d
     &0/1440d0 + (myi*pi*79d0*szp)/1800d0 - (myi*pi**3*szp)/9
     &60d0 - (myi*pi*szp*llzp)/15d0 + (zeta3)/80d0)

         case(66)            !10-11

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4*7d0)/180d0 - (pi**2*ll2**2)/12
     &d0 - (ll2**4)/6d0 - 4*cli4pt5 - (13d0*ll2*zeta3)/8d0 + 
     &zp**5*((pi**2)/180d0 - 173d0/5760d0 + (pi**2*ll2)/640d0
     & - (79d0*ll2)/1800d0 - (ll2**2)/30d0 + (ll2*llzp)/15d0 
     &- (zeta3)/160d0) + zp**8*(-(479389d0/38707200d0) + (pi*
     &*2*151d0)/80640d0 + (pi**2*ll2)/8192d0 - (7493d0*ll2)/9
     &40800d0 - (151d0*ll2**2)/13440d0 + (151d0*ll2*llzp)/672
     &0d0 - (zeta3)/2048d0) + zp**3*(-(1d0/24d0) + (pi**2)/72
     &d0 - (13d0*ll2)/72d0 + (pi**2*ll2)/96d0 - (ll2**2)/12d0
     & + (ll2*llzp)/6d0 - (zeta3)/24d0) + zp*((pi**2*ll2)/8d0
     & - (zeta3)/2d0) + zp**6*((pi**2)/270d0 - 3067d0/138240d
     &0 + (pi**2*ll2)/1536d0 - (169d0*ll2)/7200d0 - (ll2**2)/
     &45d0 + (2d0*ll2*llzp)/45d0 - (zeta3)/384d0) + zp**9*(-(
     &1164053d0/121927680d0) + (pi**2*4d0)/2835d0 + (pi**2*ll
     &2)/18432d0 - (3001d0*ll2)/595350d0 - (8d0*ll2**2)/945d0
     & + (16d0*ll2*llzp)/945d0 - (zeta3)/4608d0) + zp**4*(-(5
     &d0/128d0) + (pi**2*5d0)/576d0 + (pi**2*ll2)/256d0 - (25
     &d0*ll2)/288d0 - (5d0*ll2**2)/96d0 + (5d0*ll2*llzp)/48d0
     & - (zeta3)/64d0) + zp**7*(-(39751d0/2419200d0) + (pi**2
     &*13d0)/5040d0 + (pi**2*ll2)/3584d0 - (521d0*ll2)/39200d
     &0 - (13d0*ll2**2)/840d0 + (13d0*ll2*llzp)/420d0 - (zeta
     &3)/896d0) + zp**2*((pi**2)/48d0 + (pi**2*ll2)/32d0 - (3
     &d0*ll2)/8d0 - (ll2**2)/8d0 + (ll2*llzp)/4d0 - (zeta3)/8
     &d0)

         case(67)            !100-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/60d0 + (pi**2*ll2**2)/12d0 - 
     &(ll2**4)/12d0 - 2*cli4pt5 - (zp*zeta3)/2d0 - (3d0*ll2*z
     &eta3)/4d0 + zp**5*(-(317d0/2880d0) + (pi**2)/90d0 + (ll
     &zp)/12d0 - (zeta3)/160d0) + zp**8*((pi**2*151d0)/40320d
     &0 - 41419d0/921600d0 + (539d0*llzp)/11520d0 - (zeta3)/2
     &048d0) + zp**3*((pi**2)/36d0 - 11d0/72d0 + (llzp)/12d0 
     &- (zeta3)/24d0) + zp**6*((pi**2)/135d0 - 187d0/2304d0 +
     & (5d0*llzp)/72d0 - (zeta3)/384d0) + zp**9*(-(6297983d0/
     &182891520d0) + (pi**2*8d0)/2835d0 + (22d0*llzp)/567d0 -
     & (zeta3)/4608d0) + zp**4*(-(55d0/384d0) + (pi**2*5d0)/2
     &88d0 + (3d0*llzp)/32d0 - (zeta3)/64d0) + zp**7*((pi**2*
     &13d0)/2520d0 - 3451d0/57600d0 + (41d0*llzp)/720d0 - (ze
     &ta3)/896d0) + zp**2*((pi**2)/24d0 - (zeta3)/8d0)

         case(68)            !1000

            zp = x+1d0
            szp = s(zp)

            ris = -((pi**4*23d0)/720d0) - (myi*pi**3*sz
     &p*zp)/12d0 + ((pi**2)/8d0 - (myi*pi**3*szp)/48d0)*zp**2
     & + ((pi**2)/12d0 + (myi*pi*szp)/12d0 - (myi*pi**3*szp)/
     &144d0)*zp**3 + (-(1d0/48d0) + (pi**2*5d0)/96d0 - (myi*p
     &i**3*szp)/384d0 + (myi*pi*3d0*szp)/32d0)*zp**4 + (-(1d0
     &/30d0) + (pi**2)/30d0 + (myi*pi*szp)/12d0 - (myi*pi**3*
     &szp)/960d0)*zp**5 + (-(11d0/288d0) + (pi**2)/45d0 - (my
     &i*pi**3*szp)/2304d0 + (myi*pi*5d0*szp)/72d0)*zp**6 + (-
     &(13d0/336d0) + (pi**2*13d0)/840d0 - (myi*pi**3*szp)/537
     &6d0 + (myi*pi*41d0*szp)/720d0)*zp**7 + ((pi**2*151d0)/1
     &3440d0 - 427d0/11520d0 - (myi*pi**3*szp)/12288d0 + (myi
     &*pi*539d0*szp)/11520d0)*zp**8 + (-(14d0/405d0) + (pi**2
     &*8d0)/945d0 - (myi*pi**3*szp)/27648d0 + (myi*pi*22d0*sz
     &p)/567d0)*zp**9 + (myi*pi**3*szp*ll2)/6d0 - (myi*pi*3d0
     &*szp*zeta3)/4d0

         case(69)            !1001

            zp = x+1d0

            ris = -((pi**4)/288d0) - (3d0*zp*zeta3)/8d0
     & + (3d0*ll2*zeta3)/4d0 + zp**3*((pi**2)/72d0 - (ll2)/12
     &d0 - (zeta3)/32d0) + zp**4*((pi**2*5d0)/576d0 + 1d0/96d
     &0 - (3d0*ll2)/32d0 - (3d0*zeta3)/256d0) + zp**2*((pi**2
     &)/48d0 - (3d0*zeta3)/32d0) + zp**7*(47d0/2880d0 + (pi**
     &2*13d0)/5040d0 - (41d0*ll2)/720d0 - (3d0*zeta3)/3584d0)
     & + zp**6*((pi**2)/270d0 + 13d0/768d0 - (5d0*ll2)/72d0 -
     & (zeta3)/512d0) + zp**9*(481d0/35840d0 + (pi**2*4d0)/28
     &35d0 - (22d0*ll2)/567d0 - (zeta3)/6144d0) + zp**5*((pi*
     &*2)/180d0 + 1d0/64d0 - (ll2)/12d0 - (3d0*zeta3)/640d0) 
     &+ zp**8*((pi**2*151d0)/80640d0 + 1379d0/92160d0 - (539d
     &0*ll2)/11520d0 - (3d0*zeta3)/8192d0)

         case(70)            !101-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4*5d0)/144d0) + (pi**2*ll2**2)
     &/12d0 + (ll2**4)/6d0 + 4*cli4pt5 + ll2*zeta3 + zp**5*(-
     &((pi**2)/180d0) + 1367d0/28800d0 - (pi**2*ll2)/640d0 + 
     &(ll2**2)/30d0 - (llzp)/30d0 + (13d0*zeta3)/1280d0) + zp
     &**8*(306053d0/18063360d0 - (pi**2*151d0)/80640d0 - (pi*
     &*2*ll2)/8192d0 + (151d0*ll2**2)/13440d0 - (935d0*llzp)/
     &64512d0 + (13d0*zeta3)/16384d0) + zp*(-((pi**2*ll2)/8d0
     &) + (13d0*zeta3)/16d0) + zp**3*(11d0/144d0 - (pi**2)/72
     &d0 - (pi**2*ll2)/96d0 + (ll2**2)/12d0 - (llzp)/24d0 + (
     &13d0*zeta3)/192d0) + zp**6*(-((pi**2)/270d0) + 7639d0/2
     &30400d0 - (pi**2*ll2)/1536d0 + (ll2**2)/45d0 - (97d0*ll
     &zp)/3840d0 + (13d0*zeta3)/3072d0) + zp**9*(23078341d0/1
     &828915200d0 - (pi**2*4d0)/2835d0 - (pi**2*ll2)/18432d0 
     &+ (8d0*ll2**2)/945d0 - (2041d0*llzp)/181440d0 + (13d0*z
     &eta3)/36864d0) + zp**4*(19d0/288d0 - (pi**2*5d0)/576d0 
     &- (pi**2*ll2)/256d0 + (5d0*ll2**2)/96d0 - (llzp)/24d0 +
     & (13d0*zeta3)/512d0) + zp**2*(-((pi**2)/48d0) - (pi**2*
     &ll2)/32d0 + (ll2**2)/8d0 + (13d0*zeta3)/64d0) + zp**7*(
     &3671d0/156800d0 - (pi**2*13d0)/5040d0 - (pi**2*ll2)/358
     &4d0 + (13d0*ll2**2)/840d0 - (767d0*llzp)/40320d0 + (13d
     &0*zeta3)/7168d0)

         case(71)            !1010

            zp = x+1d0
            szp = s(zp)

            ris = -((pi**4*11d0)/288d0) + (myi*pi**3*sz
     &p*ll2)/12d0 - (pi**2*ll2**2)/6d0 + (ll2**4)/6d0 + 4*cli
     &4pt5 - (myi*pi*szp*zeta3)/4d0 + 2*ll2*zeta3 + zp**3*(-(
     &(pi**2)/72d0) - (myi*pi*szp)/24d0 - (myi*pi**3*szp)/288
     &d0 + (myi*pi*szp*ll2)/6d0 + (zeta3)/16d0) + zp**6*(17d0
     &/1152d0 - (pi**2)/270d0 - (myi*pi**3*szp)/4608d0 - (myi
     &*pi*97d0*szp)/3840d0 + (myi*pi*2d0*szp*ll2)/45d0 + (zet
     &a3)/256d0) + zp**9*(14081d0/1451520d0 - (pi**2*4d0)/283
     &5d0 - (myi*pi*2041d0*szp)/181440d0 - (myi*pi**3*szp)/55
     &296d0 + (myi*pi*16d0*szp*ll2)/945d0 + (zeta3)/3072d0) +
     & zp**4*(-((pi**2*5d0)/576d0) + 1d0/96d0 - (myi*pi*szp)/
     &24d0 - (myi*pi**3*szp)/768d0 + (myi*pi*5d0*szp*ll2)/48d
     &0 + (3d0*zeta3)/128d0) + zp**2*(-((pi**2)/48d0) - (myi*
     &pi**3*szp)/96d0 + (myi*pi*szp*ll2)/4d0 + (3d0*zeta3)/16
     &d0) + zp**7*(179d0/13440d0 - (pi**2*13d0)/5040d0 - (myi
     &*pi**3*szp)/10752d0 - (myi*pi*767d0*szp)/40320d0 + (myi
     &*pi*13d0*szp*ll2)/420d0 + (3d0*zeta3)/1792d0) + zp**5*(
     &-((pi**2)/180d0) + 7d0/480d0 - (myi*pi**3*szp)/1920d0 -
     & (myi*pi*szp)/30d0 + (myi*pi*szp*ll2)/15d0 + (3d0*zeta3
     &)/320d0) + zp**8*(-((pi**2*151d0)/80640d0) + 1057d0/921
     &60d0 - (myi*pi**3*szp)/24576d0 - (myi*pi*935d0*szp)/645
     &12d0 + (myi*pi*151d0*szp*ll2)/6720d0 + (3d0*zeta3)/4096
     &d0) + zp*(-((myi*pi**3*szp)/24d0) + (3d0*zeta3)/4d0)

         case(72)            !1011

            zp = x+1d0

            ris = (pi**4)/30d0 + (pi**2*ll2**2)/8d0 - (
     &ll2**4)/8d0 - 3*cli4pt5 + (zp*zeta3)/16d0 - (11d0*ll2*z
     &eta3)/4d0 + zp**5*(-(13d0/1920d0) + (ll2)/30d0 - (ll2**
     &2)/30d0 + (zeta3)/1280d0) + zp**8*(-(2321d0/516096d0) +
     & (935d0*ll2)/64512d0 - (151d0*ll2**2)/13440d0 + (zeta3)
     &/16384d0) + zp**3*((ll2)/24d0 - (ll2**2)/12d0 + (zeta3)
     &/192d0) + zp**6*(-(37d0/5760d0) + (97d0*ll2)/3840d0 - (
     &ll2**2)/45d0 + (zeta3)/3072d0) + zp**9*(-(157d0/43008d0
     &) + (2041d0*ll2)/181440d0 - (8d0*ll2**2)/945d0 + (zeta3
     &)/36864d0) + zp**4*(-(1d0/192d0) + (ll2)/24d0 - (5d0*ll
     &2**2)/96d0 + (zeta3)/512d0) + zp**2*(-((ll2**2)/8d0) + 
     &(zeta3)/64d0) + zp**7*(-(221d0/40320d0) + (767d0*ll2)/4
     &0320d0 - (13d0*ll2**2)/840d0 + (zeta3)/7168d0)

         case(73)            !11-1-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/720d0 + (ll2**4)/24d0 - (ll2*
     &zeta3)/8d0 + zp**7*(104641d0/64512000d0 + (pi**2*ll2)/1
     &0752d0 - (ll2**3)/5376d0 - (947d0*llzp)/460800d0 + (7d0
     &*llzp**2)/5120d0 - (zeta3)/1024d0) + zp**5*(2671d0/2764
     &80d0 + (pi**2*ll2)/1920d0 - (ll2**3)/960d0 - (53d0*llzp
     &)/4608d0 + (5d0*llzp**2)/768d0 - (7d0*zeta3)/1280d0) + 
     &zp**8*(140539517d0/202309632000d0 + (pi**2*ll2)/24576d0
     & - (ll2**3)/12288d0 - (647707d0*llzp)/722534400d0 + (36
     &3d0*llzp**2)/573440d0 - (7d0*zeta3)/16384d0) + zp*((pi*
     &*2*ll2)/24d0 - (ll2**3)/12d0 - (7d0*zeta3)/16d0) + zp**
     &3*(41d0/576d0 + (pi**2*ll2)/288d0 - (ll2**3)/144d0 - (7
     &d0*llzp)/96d0 + (llzp**2)/32d0 - (7d0*zeta3)/192d0) + z
     &p**6*(322493d0/82944000d0 + (pi**2*ll2)/4608d0 - (ll2**
     &3)/2304d0 - (2213d0*llzp)/460800d0 + (137d0*llzp**2)/46
     &080d0 - (7d0*zeta3)/3072d0) + zp**9*(2486560891d0/81935
     &40096000d0 + (pi**2*ll2)/55296d0 - (ll2**3)/27648d0 - (
     &1290829d0*llzp)/3251404800d0 + (761d0*llzp**2)/2580480d
     &0 - (7d0*zeta3)/36864d0) + zp**4*(1397d0/55296d0 + (pi*
     &*2*ll2)/768d0 - (ll2**3)/384d0 - (131d0*llzp)/4608d0 + 
     &(11d0*llzp**2)/768d0 - (7d0*zeta3)/512d0) + zp**2*(7d0/
     &32d0 + (pi**2*ll2)/96d0 - (ll2**3)/48d0 - (3d0*llzp)/16
     &d0 + (llzp**2)/16d0 - (7d0*zeta3)/64d0)

         case(74)            !11-10

            zp = x+1d0
            llzp = log(zp)
            szp = s(zp)

            ris = (pi**4*7d0)/288d0 + (pi**2*ll2**2)/24
     &d0 + (myi*pi*szp*ll2**3)/6d0 - (ll2**4)/12d0 - 2*cli4pt
     &5 - (myi*pi*szp*zeta3)/8d0 - (13d0*ll2*zeta3)/8d0 + zp*
     &*5*(-(107d0/5760d0) + (pi**2*5d0)/2304d0 + (myi*pi**3*s
     &zp)/1920d0 - (myi*pi*53d0*szp)/4608d0 + (pi**2*ll2)/192
     &0d0 - (myi*pi*szp*ll2**2)/320d0 + (myi*pi*5d0*szp*llzp)
     &/384d0 - (zeta3)/160d0) + zp**8*(-(115543d0/38707200d0)
     & + (pi**2*121d0)/573440d0 + (myi*pi**3*szp)/24576d0 - (
     &myi*pi*647707d0*szp)/722534400d0 + (pi**2*ll2)/24576d0 
     &- (myi*pi*szp*ll2**2)/4096d0 + (myi*pi*363d0*szp*llzp)/
     &286720d0 - (zeta3)/2048d0) + zp**3*(-(1d0/24d0) + (pi**
     &2)/96d0 + (myi*pi**3*szp)/288d0 - (myi*pi*7d0*szp)/96d0
     & + (pi**2*ll2)/288d0 - (myi*pi*szp*ll2**2)/48d0 + (myi*
     &pi*szp*llzp)/16d0 - (zeta3)/24d0) + zp*((myi*pi**3*szp)
     &/24d0 + (pi**2*ll2)/24d0 - (myi*pi*szp*ll2**2)/4d0 - (z
     &eta3)/2d0) + zp**6*((pi**2*137d0)/138240d0 - 79d0/7680d
     &0 - (myi*pi*2213d0*szp)/460800d0 + (myi*pi**3*szp)/4608
     &d0 + (pi**2*ll2)/4608d0 - (myi*pi*szp*ll2**2)/768d0 + (
     &myi*pi*137d0*szp*llzp)/23040d0 - (zeta3)/384d0) + zp**9
     &*((pi**2*761d0)/7741440d0 - 983419d0/609638400d0 - (myi
     &*pi*1290829d0*szp)/3251404800d0 + (myi*pi**3*szp)/55296
     &d0 + (pi**2*ll2)/55296d0 - (myi*pi*szp*ll2**2)/9216d0 +
     & (myi*pi*761d0*szp*llzp)/1290240d0 - (zeta3)/4608d0) + 
     &zp**4*((pi**2*11d0)/2304d0 - 1d0/32d0 - (myi*pi*131d0*s
     &zp)/4608d0 + (myi*pi**3*szp)/768d0 + (pi**2*ll2)/768d0 
     &- (myi*pi*szp*ll2**2)/128d0 + (myi*pi*11d0*szp*llzp)/38
     &4d0 - (zeta3)/64d0) + zp**7*(-(13441d0/2419200d0) + (pi
     &**2*7d0)/15360d0 + (myi*pi**3*szp)/10752d0 - (myi*pi*94
     &7d0*szp)/460800d0 + (pi**2*ll2)/10752d0 - (myi*pi*szp*l
     &l2**2)/1792d0 + (myi*pi*7d0*szp*llzp)/2560d0 - (zeta3)/
     &896d0) + zp**2*((pi**2)/48d0 - (myi*pi*3d0*szp)/16d0 + 
     &(myi*pi**3*szp)/96d0 + (pi**2*ll2)/96d0 - (myi*pi*szp*l
     &l2**2)/16d0 + (myi*pi*szp*llzp)/8d0 - (zeta3)/8d0)

         case(75)            !11-11

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/30d0) - (pi**2*ll2**2)/8d0 
     &+ (ll2**4)/12d0 + 3*cli4pt5 + (11d0*ll2*zeta3)/4d0 + zp
     &**6*(-((pi**2*137d0)/276480d0) + 37d0/9216d0 + (2213d0*
     &ll2)/460800d0 - (pi**2*ll2)/4608d0 + (137d0*ll2**2)/460
     &80d0 + (ll2**3)/2304d0 - (137d0*ll2*llzp)/23040d0 + (ze
     &ta3)/1536d0) + zp**9*(-((pi**2*761d0)/15482880d0) + 926
     &1559d0/19508428800d0 + (1290829d0*ll2)/3251404800d0 - (
     &pi**2*ll2)/55296d0 + (761d0*ll2**2)/2580480d0 + (ll2**3
     &)/27648d0 - (761d0*ll2*llzp)/1290240d0 + (zeta3)/18432d
     &0) + zp**4*(-((pi**2*11d0)/4608d0) + 11d0/768d0 + (131d
     &0*ll2)/4608d0 - (pi**2*ll2)/768d0 + (11d0*ll2**2)/768d0
     & + (ll2**3)/384d0 - (11d0*ll2*llzp)/384d0 + (zeta3)/256
     &d0) + zp**2*(-((pi**2)/96d0) + (3d0*ll2)/16d0 - (pi**2*
     &ll2)/96d0 + (ll2**2)/16d0 + (ll2**3)/48d0 - (ll2*llzp)/
     &8d0 + (zeta3)/32d0) + zp**7*(38569d0/19353600d0 - (pi**
     &2*7d0)/30720d0 - (pi**2*ll2)/10752d0 + (947d0*ll2)/4608
     &00d0 + (7d0*ll2**2)/5120d0 + (ll2**3)/5376d0 - (7d0*ll2
     &*llzp)/2560d0 + (zeta3)/3584d0) + zp**5*(181d0/23040d0 
     &- (pi**2*5d0)/4608d0 - (pi**2*ll2)/1920d0 + (53d0*ll2)/
     &4608d0 + (5d0*ll2**2)/768d0 + (ll2**3)/960d0 - (5d0*ll2
     &*llzp)/384d0 + (zeta3)/640d0) + zp**8*(-((pi**2*121d0)/
     &1146880d0) + 43171d0/44236800d0 - (pi**2*ll2)/24576d0 +
     & (647707d0*ll2)/722534400d0 + (363d0*ll2**2)/573440d0 +
     & (ll2**3)/12288d0 - (363d0*ll2*llzp)/286720d0 + (zeta3)
     &/8192d0) + zp*(-((pi**2*ll2)/24d0) + (ll2**3)/12d0 + (z
     &eta3)/8d0) + zp**3*(-((pi**2)/192d0) + 1d0/48d0 - (pi**
     &2*ll2)/288d0 + (7d0*ll2)/96d0 + (ll2**2)/32d0 + (ll2**3
     &)/144d0 - (ll2*llzp)/16d0 + (zeta3)/96d0)

         case(76)            !110-1

            zp = x+1d0
            llzp = log(zp)

            ris = -((pi**4)/480d0) - (pi**2*ll2**2)/12d
     &0 + (5d0*ll2*zeta3)/8d0 + zp**5*(-((pi**2*5d0)/2304d0) 
     &+ 77d0/2400d0 + (pi**2*ll2)/960d0 - (llzp)/40d0 - (zeta
     &3)/256d0) + zp**8*(232819d0/45158400d0 - (pi**2*121d0)/
     &573440d0 + (pi**2*ll2)/12288d0 - (1019d0*llzp)/161280d0
     & - (5d0*zeta3)/16384d0) + zp*((pi**2*ll2)/12d0 - (5d0*z
     &eta3)/16d0) + zp**3*(11d0/144d0 - (pi**2)/96d0 + (pi**2
     &*ll2)/144d0 - (llzp)/24d0 - (5d0*zeta3)/192d0) + zp**6*
     &(-((pi**2*137d0)/138240d0) + 169d0/9600d0 + (pi**2*ll2)
     &/2304d0 - (23d0*llzp)/1440d0 - (5d0*zeta3)/3072d0) + zp
     &**9*(40487d0/14288400d0 - (pi**2*761d0)/7741440d0 + (pi
     &**2*ll2)/27648d0 - (23d0*llzp)/5670d0 - (5d0*zeta3)/368
     &64d0) + zp**4*(-((pi**2*11d0)/2304d0) + 127d0/2304d0 + 
     &(pi**2*ll2)/384d0 - (7d0*llzp)/192d0 - (5d0*zeta3)/512d
     &0) + zp**2*(-((pi**2)/48d0) + (pi**2*ll2)/48d0 - (5d0*z
     &eta3)/64d0) + zp**7*(13423d0/1411200d0 - (pi**2*7d0)/15
     &360d0 + (pi**2*ll2)/5376d0 - (101d0*llzp)/10080d0 - (5d
     &0*zeta3)/7168d0)

         case(77)            !1100

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4)/48d0 - (myi*pi**3*szp*ll2)/12
     &d0 - (pi**2*ll2**2)/6d0 - (ll2**4)/12d0 - 2*cli4pt5 + (
     &myi*pi*szp*zeta3)/8d0 - ll2*zeta3 + zp**3*(-((pi**2)/32
     &d0) - (myi*pi*szp)/24d0 + (myi*pi**3*szp)/288d0 + (pi**
     &2*ll2)/48d0 - (zeta3)/32d0) + zp**4*(-((pi**2*11d0)/768
     &d0) + 1d0/96d0 + (myi*pi**3*szp)/768d0 - (myi*pi*7d0*sz
     &p)/192d0 + (pi**2*ll2)/128d0 - (3d0*zeta3)/256d0) + zp*
     &*2*(-((pi**2)/16d0) + (myi*pi**3*szp)/96d0 + (pi**2*ll2
     &)/16d0 - (3d0*zeta3)/32d0) + zp**7*(167d0/16128d0 - (pi
     &**2*7d0)/5120d0 - (myi*pi*101d0*szp)/10080d0 + (myi*pi*
     &*3*szp)/10752d0 + (pi**2*ll2)/1792d0 - (3d0*zeta3)/3584
     &d0) + zp**6*(29d0/2304d0 - (pi**2*137d0)/46080d0 - (myi
     &*pi*23d0*szp)/1440d0 + (myi*pi**3*szp)/4608d0 + (pi**2*
     &ll2)/768d0 - (zeta3)/512d0) + zp**9*(2569d0/414720d0 - 
     &(pi**2*761d0)/2580480d0 + (myi*pi**3*szp)/55296d0 - (my
     &i*pi*23d0*szp)/5670d0 + (pi**2*ll2)/9216d0 - (zeta3)/61
     &44d0) + zp**5*(-((pi**2*5d0)/768d0) + 13d0/960d0 + (myi
     &*pi**3*szp)/1920d0 - (myi*pi*szp)/40d0 + (pi**2*ll2)/32
     &0d0 - (3d0*zeta3)/640d0) + zp**8*(-((pi**2*363d0)/57344
     &0d0) + 497d0/61440d0 - (myi*pi*1019d0*szp)/161280d0 + (
     &myi*pi**3*szp)/24576d0 + (pi**2*ll2)/4096d0 - (3d0*zeta
     &3)/8192d0) + zp*((myi*pi**3*szp)/24d0 + (pi**2*ll2)/4d0
     & - (3d0*zeta3)/8d0)

         case(78)            !1101

            zp = x+1d0

            ris = -((pi**4)/30d0) - (pi**2*ll2**2)/6d0 
     &+ (ll2**4)/8d0 + 3*cli4pt5 + (23d0*ll2*zeta3)/8d0 + zp*
     &*6*(-((pi**2*137d0)/276480d0) - 31d0/5760d0 + (23d0*ll2
     &)/1440d0 + (pi**2*ll2)/4608d0 - (zeta3)/1536d0) + zp**9
     &*(-(43d0/20160d0) - (pi**2*761d0)/15482880d0 + (pi**2*l
     &l2)/55296d0 + (23d0*ll2)/5670d0 - (zeta3)/18432d0) + zp
     &**4*(-(1d0/192d0) - (pi**2*11d0)/4608d0 + (pi**2*ll2)/7
     &68d0 + (7d0*ll2)/192d0 - (zeta3)/256d0) + zp**2*(-((pi*
     &*2)/96d0) + (pi**2*ll2)/96d0 - (zeta3)/32d0) + zp**7*(-
     &(221d0/53760d0) - (pi**2*7d0)/30720d0 + (101d0*ll2)/100
     &80d0 + (pi**2*ll2)/10752d0 - (zeta3)/3584d0) + zp**5*(-
     &(1d0/160d0) - (pi**2*5d0)/4608d0 + (pi**2*ll2)/1920d0 +
     & (ll2)/40d0 - (zeta3)/640d0) + zp**8*(-((pi**2*121d0)/1
     &146880d0) - 7709d0/2580480d0 + (1019d0*ll2)/161280d0 + 
     &(pi**2*ll2)/24576d0 - (zeta3)/8192d0) + zp*((pi**2*ll2)
     &/24d0 - (zeta3)/8d0) + zp**3*(-((pi**2)/192d0) + (ll2)/
     &24d0 + (pi**2*ll2)/288d0 - (zeta3)/96d0)

         case(79)            !111-1

            zp = x+1d0
            llzp = log(zp)

            ris = (pi**4)/90d0 + (pi**2*ll2**2)/24d0 - 
     &(ll2**4)/12d0 - cli4pt5 - (7d0*ll2*zeta3)/8d0 + zp**5*(
     &-(599d0/46080d0) + (pi**2*5d0)/4608d0 - (5d0*ll2**2)/76
     &8d0 + (ll2**3)/960d0 + (7d0*llzp)/768d0 - (zeta3)/1280d
     &0) + zp**8*((pi**2*121d0)/1146880d0 - 21977d0/14745600d
     &0 - (363d0*ll2**2)/573440d0 + (ll2**3)/12288d0 + (469d0
     &*llzp)/368640d0 - (zeta3)/16384d0) + zp*((ll2**3)/12d0 
     &- (zeta3)/16d0) + zp**3*((pi**2)/192d0 - 11d0/288d0 - (
     &ll2**2)/32d0 + (ll2**3)/144d0 + (llzp)/48d0 - (zeta3)/1
     &92d0) + zp**6*((pi**2*137d0)/276480d0 - 79d0/12288d0 - 
     &(137d0*ll2**2)/46080d0 + (ll2**3)/2304d0 + (5d0*llzp)/1
     &024d0 - (zeta3)/3072d0) + zp**9*((pi**2*761d0)/15482880
     &d0 - 83359739d0/117050572800d0 - (761d0*ll2**2)/2580480
     &d0 + (ll2**3)/27648d0 + (29531d0*llzp)/46448640d0 - (ze
     &ta3)/36864d0) + zp**4*((pi**2*11d0)/4608d0 - 19d0/768d0
     & - (11d0*ll2**2)/768d0 + (ll2**3)/384d0 + (llzp)/64d0 -
     & (zeta3)/512d0) + zp**2*((pi**2)/96d0 - (ll2**2)/16d0 +
     & (ll2**3)/48d0 - (zeta3)/64d0) + zp**7*(-(3343d0/107520
     &0d0) + (pi**2*7d0)/30720d0 - (7d0*ll2**2)/5120d0 + (ll2
     &**3)/5376d0 + (29d0*llzp)/11520d0 - (zeta3)/7168d0)

         case(80)            !1110

            zp = x+1d0
            szp = s(zp)

            ris = (pi**4)/90d0 + (pi**2*ll2**2)/12d0 - 
     &(myi*pi*szp*ll2**3)/6d0 - (ll2**4)/24d0 - cli4pt5 - ll2
     &*zeta3 + zp**5*(-(11d0/1920d0) + (pi**2*5d0)/4608d0 + (
     &myi*pi*7d0*szp)/768d0 - (pi**2*ll2)/1920d0 - (myi*pi*5d
     &0*szp*ll2)/384d0 + (myi*pi*szp*ll2**2)/320d0 + (zeta3)/
     &1280d0) + zp**8*((pi**2*121d0)/1146880d0 - 563d0/286720
     &d0 + (myi*pi*469d0*szp)/368640d0 - (pi**2*ll2)/24576d0 
     &- (myi*pi*363d0*szp*ll2)/286720d0 + (myi*pi*szp*ll2**2)
     &/4096d0 + (zeta3)/16384d0) + zp*(-((pi**2*ll2)/24d0) + 
     &(myi*pi*szp*ll2**2)/4d0 + (zeta3)/16d0) + zp**3*((pi**2
     &)/192d0 + (myi*pi*szp)/48d0 - (pi**2*ll2)/288d0 - (myi*
     &pi*szp*ll2)/16d0 + (myi*pi*szp*ll2**2)/48d0 + (zeta3)/1
     &92d0) + zp**6*(-(103d0/23040d0) + (pi**2*137d0)/276480d
     &0 + (myi*pi*5d0*szp)/1024d0 - (pi**2*ll2)/4608d0 - (myi
     &*pi*137d0*szp*ll2)/23040d0 + (myi*pi*szp*ll2**2)/768d0 
     &+ (zeta3)/3072d0) + zp**9*(-(203d0/165888d0) + (pi**2*7
     &61d0)/15482880d0 + (myi*pi*29531d0*szp)/46448640d0 - (p
     &i**2*ll2)/55296d0 - (myi*pi*761d0*szp*ll2)/1290240d0 + 
     &(myi*pi*szp*ll2**2)/9216d0 + (zeta3)/36864d0) + zp**4*(
     &-(1d0/192d0) + (pi**2*11d0)/4608d0 + (myi*pi*szp)/64d0 
     &- (pi**2*ll2)/768d0 - (myi*pi*11d0*szp*ll2)/384d0 + (my
     &i*pi*szp*ll2**2)/128d0 + (zeta3)/512d0) + zp**2*((pi**2
     &)/96d0 - (pi**2*ll2)/96d0 - (myi*pi*szp*ll2)/8d0 + (myi
     &*pi*szp*ll2**2)/16d0 + (zeta3)/64d0) + zp**7*(-(493d0/1
     &61280d0) + (pi**2*7d0)/30720d0 + (myi*pi*29d0*szp)/1152
     &0d0 - (pi**2*ll2)/10752d0 - (myi*pi*7d0*szp*ll2)/2560d0
     & + (myi*pi*szp*ll2**2)/1792d0 + (zeta3)/7168d0)

         case(81)            !1111

            zp = x+1d0

            ris = -((zp*ll2**3)/12d0) + (ll2**4)/24d0 +
     & zp**8*(967d0/1474560d0 - (469d0*ll2)/368640d0 + (363d0
     &*ll2**2)/573440d0 - (ll2**3)/12288d0) + zp**3*(-((ll2)/
     &48d0) + (ll2**2)/32d0 - (ll2**3)/144d0) + zp**6*(17d0/9
     &216d0 - (5d0*ll2)/1024d0 + (137d0*ll2**2)/46080d0 - (ll
     &2**3)/2304d0) + zp**9*(89d0/245760d0 - (29531d0*ll2)/46
     &448640d0 + (761d0*ll2**2)/2580480d0 - (ll2**3)/27648d0)
     & + zp**4*(1d0/384d0 - (ll2)/64d0 + (11d0*ll2**2)/768d0 
     &- (ll2**3)/384d0) + zp**2*((ll2**2)/16d0 - (ll2**3)/48d
     &0) + zp**7*(7d0/6144d0 - (29d0*ll2)/11520d0 + (7d0*ll2*
     &*2)/5120d0 - (ll2**3)/5376d0) + zp**5*(1d0/384d0 - (7d0
     &*ll2)/768d0 + (5d0*ll2**2)/768d0 - (ll2**3)/960d0)
c End of expansions around x = -1

      end select

c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero).
      if (bcflag.eq.1) then
         xre = dreal(x)

         if (n4.eq.0.and.xre.gt.0d0) then
            if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c
         else if (n4.eq.1.and.xre.lt.1d0) then
            if (n1.ne.-1.and.n2.ne.-1.and.n3.ne.-1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.gt.-1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
c            
         else if (n4.eq.-1.and.xre.gt.-1d0) then
            if (n1.ne.1.and.n2.ne.1.and.n3.ne.1) then
               ris = dcmplx(dreal(ris),0d0)
            else if (xre.lt.1d0) then
               ris = dcmplx(dreal(ris),0d0)
            endif
            
         endif
      endif

      HPL4arm1=ris
      return
      end function
