      double complex function HPL4atm1(n1,n2,n3,n4)
      implicit none
      integer n1,n2,n3,n4,j
      double complex ris,myi,cli4pt5,cli4
      double precision pi, zeta2, zeta3,zeta4,ll2

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)

      ll2 = dlog(2d0)
      cli4pt5 = cli4(dcmplx(0.5d0,0d0))

      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      ris = dcmplx(0d0,0d0)

      if (j.ge.28) then
         select case (j)
         case(28)               !0-1-1-1
            ris = -pi**4/90d0
         case(29)               !0-1-10
            ris = -pi**4/72d0 + myi*pi*zeta3
         case(30)               !0-1-11
            ris = pi**4/80d0 - (pi**2*ll2**2)/12d0 
     &           - ll2**4/24d0 - cli4pt5 
     &           - (7*ll2*zeta3)/8d0
         case(31)               !0-10-1
            ris = pi**4/120d0
         case(32)               !0-100
            ris = pi**4/20d0 + 2d0*myi*pi*zeta3
         case(33)               !0-101
            ris = (71*pi**4)/1440d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 4*cli4pt5 
     &           - (7*ll2*zeta3)/2d0
         case(34)               !0-11-1
            ris = (-7*pi**4)/288d0 + (pi**2*ll2**2)/8d0 
     &           + ll2**4/8d0 + 3*cli4pt5
         case(35)               !0-110
            ris = (-71*pi**4)/1440d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + myi*pi*((pi**2*ll2)/4d0 
     &           - zeta3) + (7*ll2*zeta3)/2d0 
     &           - 2*((-19*pi**4)/1440d0 
     &           + (7*ll2*zeta3)/4d0)
         case(36)               !0-111
            ris = (-7*pi**4)/288d0 - (5*pi**2*ll2**2)/24d0 
     &           + ll2**4/12d0 + 2*cli4pt5 
     &           + (21*ll2*zeta3)/8d0
         case(37)               !00-1-1
            ris = pi**4/360d0
         case(38)               !00-10
            ris = pi**4/30d0 - myi*pi*zeta3
         case(39)               !00-11
            ris = (-19*pi**4)/1440d0 + (7*ll2*zeta3)/4d0
         case(40)               !000-1
            ris = -pi**4/90d0
         case(41)               !0000
            ris = pi**4/24d0
         case(42)               !0001
            ris = (-7*pi**4)/720d0
         case(43)               !001-1
            ris = -pi**4/180d0 - (pi**2*ll2**2)/12d0 
     &           + ll2**4/12d0 + 2*cli4pt5
         case(44)               !0010
            ris = (7*pi**4)/240d0 - myi*3d0/4d0*pi*zeta3
         case(45)               !0011
            ris = -pi**4/48d0 - (pi**2*ll2**2)/12d0 
     &           + ll2**4/12d0 + 2*cli4pt5 
     &           + (7*ll2*zeta3)/4d0
         case(46)               !01-1-1
            ris = (11*pi**4)/720d0 - ll2**4/8d0 - 3*cli4pt5
         case(47)               !01-10
            ris = -pi**4/480d0 - 2*(-pi**4/180d0 
     &           - (pi**2*ll2**2)/12d0 + ll2**4/12d0 
     &           + 2*cli4pt5) 
     &           + myi*pi*(-(pi**2*ll2)/4d0 + (13*zeta3)/8d0)
         case(48)               !01-11
            ris = (7*pi**4)/720d0 + (pi**2*ll2**2)/4d0 
     &           - (21*ll2*zeta3)/8d0
         case(49)               !010-1
            ris = pi**4/480d0
         case(50)               !0100
            ris = pi**4/80d0 + myi*3d0/2d0*pi*zeta3
         case(51)               !0101
            ris = (13*pi**4)/288d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 4*cli4pt5 
     &           - (7*ll2*zeta3)/2d0
         case(52)               !011-1
            ris = pi**4/80d0 - (pi**2*ll2**2)/24d0 
     &           - ll2**4/12d0 - 2*cli4pt5
         case(53)               !0110
            ris = (-13*pi**4)/288d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + myi/8d0*pi*zeta3 
     &           + (7*ll2*zeta3)/2d0 - 2*(-pi**4/48d0 
     &           - (pi**2*ll2**2)
     &           /12d0 + ll2**4/12d0 + 2*cli4pt5 
     &           + (7*ll2*zeta3)/4d0)
         case(54)               !0111
            ris = -pi**4/90d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + cli4pt5 
     &           + (7*ll2*zeta3)/8d0
         case(55)               !1-1-1-1
            ris = cli4pt5
         case(56)               !1-1-10
            ris = pi**4/720d0 + ll2**4/24d0 + cli4pt5 
     &           - myi*pi*(-(pi**2*ll2)/12d0 + ll2**3/6d0 
     &           + (7*zeta3)/8d0) - (ll2*zeta3)/8d0
         case(57)               !1-1-11
            ris = -pi**4/288d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + (7*ll2*zeta3)/8d0
         case(58)               !1-10-1
            ris = (-19*pi**4)/1440d0 + (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + cli4pt5 
     &           + (ll2*zeta3)/4d0
         case(59)               !1-100
            ris = (19*pi**4)/1440d0 - (pi**2*(pi**2/12d0 
     &           - ll2**2/2d0))/2d0 
     &           - myi*pi*((pi**2*ll2)/6d0 
     &           - (5*zeta3)/8d0) - (3*ll2*zeta3)/4d0 
     &           - myi*pi*(-(pi**2*ll2)/4d0 + (13*zeta3)/8d0)
         case(60)               !1-101
            ris = (-11*pi**4)/240d0 - (pi**2*ll2**2)/8d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + (13*ll2*zeta3)/4d0
         case(61)               !1-11-1
            ris = pi**4/1440d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + (ll2*zeta3)/4d0
         case(62)               !1-110
            ris = -pi**4/360d0 - (pi**2*ll2**2)/24d0 
     &           - myi*pi*((pi**2*ll2)/12d0 - ll2**3/6d0 
     &           - zeta3/4d0) 
     &           + ll2*zeta3
         case(63)               !1-111
            ris = pi**4/30d0 + (pi**2*ll2**2)/6d0 
     &           - ll2**4/6d0 - 3*cli4pt5 
     &           - (23*ll2*zeta3)/8d0
         case(64)               !10-1-1
            ris = -pi**4/288d0 - (pi**2*ll2**2)/24d0 
     &           + ll2**4/24d0 + cli4pt5 
     &           - (ll2*zeta3)/8d0
         case(65)               !10-10
            ris = -pi**4/480d0 + myi*pi*((pi**2*ll2)/6d0 
     &           - (5*zeta3)/8d0) - 2*(pi**4/60d0 
     &           + (pi**2*ll2**2)/12d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - (3*ll2*zeta3)/4d0)
         case(66)               !10-11
            ris = (7*pi**4)/180d0 - (pi**2*ll2**2)/12d0 
     &           - ll2**4/6d0 - 4*cli4pt5 
     &           - (13*ll2*zeta3)/8d0
         case(67)               !100-1
            ris = pi**4/60d0 + (pi**2*ll2**2)/12d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - (3*ll2*zeta3)/4d0
         case(68)               !1000
            ris = (-23*pi**4)/720d0 + myi/6d0*pi**3*ll2 
     &           - myi*3d0/4d0*pi*zeta3
         case(69)               !1001
            ris = -pi**4/288d0 + (3*ll2*zeta3)/4d0
         case(70)               !101-1
            ris = (-5*pi**4)/144d0 + (pi**2*ll2**2)/12d0 
     &           + ll2**4/6d0 + 4*cli4pt5 + ll2*zeta3
         case(71)               !1010
            ris = (-13*pi**4)/288d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/6d0 + 4*cli4pt5 
     &           + myi*pi*((pi**2*ll2)
     &           /12d0 - zeta3/4d0) + (7*ll2*zeta3)/2d0 
     &           - 2*(-pi**4/288d0 
     &           + (3*ll2*zeta3)/4d0)
         case(72)               !1011
            ris = pi**4/30d0 + (pi**2*ll2**2)/8d0 
     &           - ll2**4/8d0 - 3*cli4pt5 
     &           - (11*ll2*zeta3)/4d0
         case(73)               !11-1-1
            ris = pi**4/720d0 + ll2**4/24d0 
     &           - (ll2*zeta3)/8d0
         case(74)               !11-10
            ris = (7*pi**4)/288d0 + (pi**2*ll2**2)/24d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - myi*pi*(-ll2**3/6d0 
     &           + zeta3/8d0) - (13*ll2*zeta3)/8d0
         case(75)               !11-11
            ris = -pi**4/30d0 - (pi**2*ll2**2)/8d0 
     &           + ll2**4/12d0 + 3*cli4pt5 
     &           + (11*ll2*zeta3)/4d0
         case(76)               !110-1
            ris = -pi**4/480d0 - (pi**2*ll2**2)/12d0 
     &           + (5*ll2*zeta3)/8d0
         case(77)               !1100
            ris = pi**4/48d0 - (pi**2*ll2**2)/6d0 
     &           - ll2**4/12d0 - 2*cli4pt5 
     &           - myi*pi*((pi**2*ll2)
     &           /12d0 - zeta3/4d0) - myi/8d0*pi*zeta3 
     &           - ll2*zeta3
         case(78)               !1101
            ris = -pi**4/30d0 - (pi**2*ll2**2)/6d0 
     &           + ll2**4/8d0 + 3*cli4pt5 
     &           + (23*ll2*zeta3)/8d0
         case(79)               !111-1
            ris = pi**4/90d0 + (pi**2*ll2**2)/24d0 
     &           - ll2**4/12d0 - cli4pt5 
     &           - (7*ll2*zeta3)/8d0
         case(80)               !1110
            ris = pi**4/90d0 + (pi**2*ll2**2)/12d0 
     &           - myi/6d0*pi*ll2**3 - ll2**4/24d0 
     &           - cli4pt5 
     &           - ll2*zeta3
         case(81)               !1111
            ris = ll2**4/24d0
         end select
      else
         print*, ""
         print*, "****************"
         print*, "ERROR in HPL4: "
         print*, "HPL4(",n1,",",n2,",",n3,",",n4
     &        ,",-1) is divergent!"
         print*, "Aborting..."
         print*,"****************"
         stop
      endif
      HPL4atm1=ris
      return
      end function
