      double complex function HPL4else(n1, n2, n3, n4, x)
      implicit  none 
      double precision pi, zeta2, zeta3, zeta4, xre
      double complex x, ris, myi
c      double complex ccli2,cli3
      double complex cli4
      double complex basis1,basis2,basis3,basis4,basis5,basis6,basis7
      double complex basis8,basis9,basis10,basis11,basis12,basis13
      double complex basis14,basis15,basis16,basis17,basis18
      double complex basis3_1,basis3_2,basis3_3,basis3_4,basis3_5
      double complex basis3_6,basis3_7,basis3_8
      double complex basis2_1,basis2_2,basis2_3
      integer n1,n2,n3,n4,j,bcflag,bcflag_save
      double complex ll2,cli4pt5,ll1px,ll1mx,llx,b16,b17,b18

c     #####################################################
c     basis1(x) = cli4(x) 
c     basis2(x) = cli4(-x)
c     basis3(x) = cli4(1-x)
c     basis4(x) = cli4(1/(1+x)) 
c     basis5(x) = cli4(x/(x-1))
c     basis6(x) = cli4(x/(x+1)) 
c     basis7(x) = cli4((1+x)/2) 
c     basis8(x) = cli4((1-x)/2)
c     basis9(x) = cli4((1-x)/(1+x))
c     basis10(x) = cli4((x-1)/(x+1))
c     basis11(x) = cli4(2x/(1+x))
c     basis12(x) = cli4(2x/(x-1)) 
c     basis13(x) = cli4(1-x^2) = cli4_sbc 
c     basis14(x) = cli4(x^2/(x^2-1)) 
c     basis15(x) = cli4(4x/(1+x)^2) = cli4_sbc_2  
c     basis16(x) = ch2m2(x) 
c     basis17(x) = ch21m1(x) 
c     basis18(x) = ch21m1(-x) 
c     #####################################################
c     basis3_1(z) = cli3(z) 
c     basis3_2(z) = cli3(-z)
c     basis3_3(z) = cli3(1-z)
c     basis3_4(z) = cli3(1/(1+z)) 
c     basis3_5(z) = cli3((1+z)/2) 
c     basis3_6(z) = cli3((1-z)/2) 
c     basis3_7(z) = cli3((1-z)/(1+z)) 
c     basis3_8(z) = cli3(2z/(z-1))
c     basis2_1(z) = ccli2(z)
c     basis2_2(z) = ccli2(-z)) 
c     basis2_3(z) = ccli2((1-z)/2)
c     #####################################################
c     #####################################################
    
      pi=3.1415926535897932385d0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)

      ll2 = dlog(2d0)
      cli4pt5 = cli4(dcmplx(0.5d0,0d0))

      j=1+(n4+1)*1+(n3+1)*3+(n2+1)*9+(n1+1)*27
      
      ris=dcmplx(0d0,0d0)
      bcflag = 0

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      ll1px = log(1d0+x)
      ll1mx = log(1d0-x)
      llx = log(x)

      select case (j)
            
      case(1)                !-1-1-1-1
            
         ris = ll1px**4/24d0
            
      case(2)                   !-1-1-10
         
         ris = -pi**4/90d0 + basis4(x) - (pi**2*ll1px**2)
     &/12d0 + ll1px**4/24d0 + ll1px*zeta3
         
      case(3)                   !-1-1-11
         
         ris = -cli4pt5 + basis7(x) 
     &+ (pi**2*ll2
     &*ll1px)/12d0 - (ll2**3*ll1px)/6d0 - (pi**2
     &*ll1px**2)/24d0 + (ll2**2*ll1px**2)/4d0 
     &- (ll2*ll1px**3)/6d0 - (7*ll1px*zeta3)/8d0

      case(4)                   !-1-10-1

         ris = pi**4/30d0 - 3*basis4(x) - basis3_4(x)
     &*ll1px + (pi**2*ll1px**2)/12d0 + ll1px**4/24d0 
     &- 2*ll1px*zeta3

      case(5)                   !-1-100

         ris = pi**4/90d0 - basis2(x) - basis4(x) 
     &- basis6(x) - basis3_4(x)*llx - (pi**2*llx
     &*ll1px)/6d0 + (pi**2*ll1px**2)/12d0 - (llx**2
     &*ll1px**2)/4d0 + (llx*ll1px**3)/3d0 
     &- ll1px**4/12d0 + llx*zeta3 - ll1px*zeta3

      case(6)                   !-1-101

         ris = pi**4/480d0 + basis3(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 - basis13(x)/4d0 
     &+ basis3_4(x)*ll1mx + (pi**2*ll1mx*ll1px)
     &/12d0 + (pi**2*ll1px**2)/12d0 - (ll1mx*ll1px**3)
     &/6d0 - (7d0*ll1mx*zeta3)/8d0 - (5d0*ll1px*zeta3)/8d0

      case(7)                   !-1-11-1

         ris = 3d0*cli4pt5 
     &- 3d0*basis7(x) 
     &+ basis3_5(x)*ll1px 
     &- (pi**2*ll2*ll1px)/6d0 
     &+ (ll2**3*ll1px)/3d0 + (pi**2*ll1px**2)/24d0 
     &- (ll2**2*ll1px**2)/4d0 + (7*ll1px*zeta3)/4d0

      case(8)                   !-1-110

         ris = pi**4/90d0 + 3*cli4pt5 + basis2(x)/2d0 
     &- basis1(x)/2d0 - basis15(x)/4d0 - basis4(x) 
     &+ basis11(x) - 3*basis7(x) 
     &+ basis3_5(x)*llx + (pi**2*ll2*llx)/12d0 
     &- (ll2**3*llx)/6d0 - (pi**2*ll2*ll1px)/4d0 
     &+ (ll2**3*ll1px)/2d0 - (pi**2*llx*ll1px)/12d0 
     &+ (ll2**2*llx*ll1px)/2d0 + (5*pi**2*ll1px**2)
     &/24d0 - (3*ll2**2*ll1px**2)/4d0 - (ll2*llx
     &*ll1px**2)/2d0 + (ll2*ll1px**3)/2d0 + (llx
     &*ll1px**3)/6d0 - ll1px**4/6d0 - (7*llx*zeta3)/8d0 
     &+ (13*ll1px*zeta3)/8d0

      case(9)                   !-1-111

         ris = (-7*pi**4)/720d0 - basis8(x) 
     &- basis10(x) + basis7(x) 
     &- basis3_5(x)*ll1mx - (pi**2*ll2*ll1mx)/6d0 
     &+ (ll2**3*ll1mx)/3d0 - (ll2**2*ll1mx**2)/4d0 
     &+ (pi**2*ll2*ll1px)/12d0 -(ll2**3*ll1px)/6d0 
     &+ (pi**2*ll1mx*ll1px)/6d0 - (ll2**2*ll1mx
     &*ll1px)/2d0 + (ll2*ll1mx**2*ll1px)/2d0 
     &- (pi**2*ll1px**2)/12d0 + (ll2**2*ll1px**2)/4d0 
     &- (ll1mx**2*ll1px**2)/4d0 + (ll1mx
     &*ll1px**3)/6d0 - ll1px**4/24d0 + ll1mx*zeta3 
     &- (ll1px*zeta3)/8d0

      case(10)                  !-10-1-1

         ris = -pi**4/30d0 + 3*basis4(x) + 2*basis3_4(x)
     &*ll1px + (pi**2*ll1px**2)/12d0 + (basis2_2(x)
     &*ll1px**2)/2d0 + (llx*ll1px**3)/2d0 
     &- (5*ll1px**4)/24d0 + ll1px*zeta3

      case(11)                  !-10-10

         ris = -pi**4/45d0 + basis2_2(x)**2/2d0 + 2*basis2(x) 
     &+ 2*basis4(x) + 2*basis6(x) + 2*basis3_4(x)
     &*llx + (pi**2*llx*ll1px)/3d0 + basis2_2(x)*llx
     &*ll1px - (pi**2*ll1px**2)/6d0 + llx**2
     &*ll1px**2 - (2*llx*ll1px**3)/3d0 + ll1px**4
     &/6d0 - 2*llx*zeta3 + 2*ll1px*zeta3

      case(12)                  !-10-11

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag = bcflag_save

         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = (-11*pi**4)/720d0 + (pi**2*basis2_2(x))/12d0 
     &- basis2_3(x)*basis2_2(x) + basis2_2(x)**2/2d0 
     &+ b18 
     &- 2*basis3(x) - (3*basis2(x))/2d0 
     &- basis1(x)/2d0 - basis15(x)/4d0 + basis4(x) 
     &+ basis9(x) - basis10(x) + basis3(x) 
     &+ basis11(x) + basis13(x)/2d0 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &- 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 
     &+ (basis2_2(x)*ll2**2)/2d0
     &+ basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0
     &- (basis2_2(x)*ll2**2)/2d0 - 2*basis3_4(x)*ll1mx 
     &+ basis2_2(x)*ll2*ll1mx + 2*basis3_1(x)*ll1px 
     &+ 2*basis3_8(x)*ll1px + 2*basis3_7(x)
     &*ll1px - (pi**2*ll1mx*ll1px)/2d0 - basis2_2(x)
     &*ll1mx*ll1px + ll2*ll1mx**2*ll1px 
     &- (ll1mx**3*ll1px)/3d0 + ll1mx**2*llx
     &*ll1px + (pi**2*ll1px**2)/8d0 - (basis2_3(x)
     &*ll1px**2)/2d0 + (basis2_2(x)*ll1px**2)/2d0
     &-(basis2_1(x)
     &*ll1px**2)/2d0 - (ll2**2*ll1px**2)/4d0 
     &- (3*ll2*ll1mx*ll1px**2)/2d0 - 2*ll1mx
     &*llx*ll1px**2 + ll2*ll1px**3 + (5*ll1mx
     &*ll1px**3)/6d0 + (5*llx*ll1px**3)/6d0 
     &- (11*ll1px**4)/24d0 + (7*ll1mx*zeta3)/4d0 
     &+ (ll1px*zeta3)/4d0

      case(13)                  !-100-1

         ris = -basis2_2(x)**2/2d0 - basis3_2(x)*ll1px

      case(14)                  !-1000

         ris = basis2(x) - basis3_2(x)*llx + (basis2_2(x)
     &*llx**2)/2d0 + (llx**3*ll1px)/6d0

      case(15)                  !-1001

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag = bcflag_save

         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/360d0 + basis2_2(x)*basis2_1(x) + b16 
     &- basis3(x) - basis2(x) - basis1(x) + basis5(x) 
     &+ basis4(x) + basis6(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_2(x)*ll1mx 
     &+ (pi**2*ll1mx**2)/16d0 + ll1mx**4/32d0 
     &- (ll1mx**3*llx)/12d0 + 2*basis3_1(x)*ll1px 
     &- (pi**2*ll1mx*ll1px)/24d0 - (ll1mx**3
     &*ll1px)/24d0 + (ll1mx**2*llx*ll1px)/4d0 
     &- (5*pi**2*ll1px**2)/48d0 - (ll1mx**2*ll1px**2)
     &/16d0 + (ll1mx*llx*ll1px**2)/4d0 - (ll1mx
     &*ll1px**3)/24d0 - (llx*ll1px**3)/12d0 
     &+ (7*ll1px**4)/96d0 + (3*ll1mx*zeta3)/4d0 
     &+ (3*ll1px*zeta3)/4d0

      case(16)                  !-101-1

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/90d0 - (pi**2*basis2_2(x))/12d0 + basis2_3(x)
     &*basis2_2(x) - basis2_2(x)**2/2d0 
     &- b18 
     &- (pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0)
     &+ (3*basis2(x))/2d0 + basis1(x)/2d0 
     &+ basis15(x)/4d0 - basis4(x) 
     &+ 2*basis6(x) - basis11(x) + (basis2_2(x)
     &*ll2**2)/2d0 - basis2_2(x)*ll2*ll1mx 
     &+ basis3_6(x)*ll1px - basis3_3(x)*ll1px 
     &- 2*basis3_1(x)*ll1px - 2*basis3_8(x)*ll1px 
     &- basis3_4(x)*ll1px - basis3_7(x)
     &*ll1px + basis3_5(x)*ll1px + (pi**2*ll2
     &*ll1px)/6d0 - (ll2**3*ll1px)/3d0 
     &+ (pi**2*ll1mx*ll1px)/4d0 + (ll2**2*ll1mx
     &*ll1px)/2d0 - ll2*ll1mx**2*ll1px 
     &+ (ll1mx**3*ll1px)/3d0 - ll1mx**2*llx
     &*ll1px - (3*pi**2*ll1px**2)/8d0 + (basis2_3(x)
     &*ll1px**2)/2d0 - (basis2_2(x)*ll1px**2)/2d0
     &+(basis2_1(x)
     &*ll1px**2)/2d0 + (3*ll2**2*ll1px**2)/4d0 
     &+ (ll2*ll1mx*ll1px**2)/2d0 + ll1mx*llx
     &*ll1px**2 - ll2*ll1px**3 - (5*llx
     &*ll1px**3)/6d0 + (11*ll1px**4)/24d0 
     &+ (ll1px*zeta3)/4d0

      case(17)                  !-1010

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -(basis2_2(x)*basis2_1(x)) - b16 
     &+ basis3_6(x)
     &*llx - basis3_3(x)*llx - basis3_4(x)*llx 
     &+ basis3_7(x)*llx + basis3_5(x)*llx 
     &+ (pi**2*ll2*llx)/6d0 - (ll2**3*llx)/3d0 
     &- (pi**2*ll1mx*llx)/12d0 - basis2_2(x)*ll1mx
     &*llx + (ll2**2*ll1mx*llx)/2d0 - 2*basis3_1(x)
     &*ll1px - (pi**2*llx*ll1px)/12d0 
     &+ (ll2**2*llx*ll1px)/2d0 - ll2*ll1mx
     &*llx*ll1px - ll1mx*llx**2*ll1px 
     &+ (ll1mx*llx*ll1px**2)/2d0 - (3*llx*zeta3)/4d0

      case(18)                  !-1011

         ris = pi**4/288d0 - basis4(x) + basis9(x)
     &/2d0 - basis10(x)/2d0 - basis13(x)/4d0 
     &- basis3_6(x)*ll1mx + basis3_3(x)*ll1mx 
     &+ basis3_4(x)*ll1mx - basis3_7(x)
     &*ll1mx - basis3_5(x)*ll1mx - (pi**2*ll2
     &*ll1mx)/6d0 + (ll2**3*ll1mx)/3d0 + (pi**2
     &*ll1mx**2)/24d0 + (basis2_2(x)*ll1mx**2)/2d0 
     &- (ll2**2*ll1mx**2)/2d0 + (pi**2*ll1mx
     &*ll1px)/4d0 - (ll2**2*ll1mx*ll1px)/2d0 
     &+ ll2*ll1mx**2*ll1px + (ll1mx**2*llx
     &*ll1px)/2d0 + (pi**2*ll1px**2)/24d0 
     &- (ll1mx**2*ll1px**2)/2d0 - ll1px**4/24d0 
     &+ (ll1mx*zeta3)/8d0 + (ll1px*zeta3)/8d0

      case(19)                  !-11-1-1

         ris = -3*cli4pt5 + 3*basis7(x) 
     &- 2*basis3_5(x)*ll1px + (pi**2*ll2*ll1px)
     &/12d0 - (ll2**3*ll1px)/6d0 + (pi**2*ll1px**2)/12d0 
     &- (basis2_3(x)*ll1px**2)/2d0 - (ll2**2
     &*ll1px**2)/2d0 + (ll2*ll1mx*ll1px**2)/2d0 
     &+ (ll2*ll1px**3)/2d0 - (ll1mx*ll1px**3)/2d0 
     &- (7*ll1px*zeta3)/8d0

      case(20)                  !-11-10

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/90d0 - basis2_2(x)**2/2d0 
     &- b18 
     &- (pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0)
     &- 6*cli4pt5 + basis2(x)/2d0 
     &+ (3*basis1(x))/2d0 + (3*basis15(x))/4d0 
     &+ basis4(x) + 2*basis6(x) - 3*basis11(x) 
     &+ 6*basis7(x) - 2*basis3_5(x)*llx 
     &- (pi**2*ll2*llx)/6d0 + (ll2**3*llx)/3d0 
     &- 2*basis3_1(x)*ll1px - 2*basis3_8(x)*ll1px 
     &- 2*basis3_7(x)*ll1px + (pi**2*ll2
     &*ll1px)/2d0 - ll2**3*ll1px 
     &+ (pi**2*ll1mx*ll1px)/3d0 - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 + (pi**2*llx
     &*ll1px)/4d0 - basis2_3(x)*llx*ll1px 
     &- (3*ll2**2*llx*ll1px)/2d0 + ll2*ll1mx
     &*llx*ll1px - ll1mx**2*llx*ll1px 
     &- (17*pi**2*ll1px**2)/24d0 + (basis2_3(x)
     &*ll1px**2)/2d0 - (basis2_2(x)*ll1px**2)/2d0 
     &+ (basis2_1(x)*ll1px**2)/2d0 
     &+ (7*ll2**2*ll1px**2)
     &/4d0 + (3*ll2*ll1mx*ll1px**2)/2d0 + ll2
     &*llx*ll1px**2 + ll1mx*llx*ll1px**2 
     &- 2*ll2*ll1px**3 - (ll1mx*ll1px**3)/2d0 
     &- (7*llx*ll1px**3)/6d0 + (19*ll1px**4)/24d0 
     &+ (7*llx*zeta3)/4d0 - (9*ll1px*zeta3)/4d0

      case(21)                  !-11-11

         ris = (11*pi**4)/480d0 - (pi**2*basis2_3(x))/12d0 
     &+ basis2_3(x)**2/2d0 + 2*basis8(x) + 2
     &*basis10(x) - 2*basis7(x) 
     &- (pi**2*ll2**2)/24d0 + (basis2_3(x)*ll2**2)/2d0 
     &+ ll2**4/8d0 + 2*basis3_5(x)*ll1mx + (5*pi**2
     &*ll2*ll1mx)/12d0 - basis2_3(x)*ll2
     &*ll1mx - (7*ll2**3*ll1mx)/6d0 + ll2**2
     &*ll1mx**2 - (pi**2*ll2*ll1px)/6d0 + (ll2**3
     &*ll1px)/3d0 - (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_3(x)*ll1mx*ll1px + (3*ll2**2
     &*ll1mx*ll1px)/2d0 - 2*ll2*ll1mx**2
     &*ll1px + (pi**2*ll1px**2)/6d0 - (ll2**2
     &*ll1px**2)/2d0 + ll1mx**2*ll1px**2 
     &- (ll1mx*ll1px**3)/3d0 + ll1px**4/12d0 
     &- 2*ll1mx*zeta3 + (ll1px*zeta3)/4d0

      case(22)                  !-110-1

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save         

         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/90d0 + basis2_2(x)**2/2d0 
     &+ b18 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0
     &- (3*basis2(x))/2d0 - basis1(x)/2d0 
     &- basis15(x)/4d0 + basis4(x) 
     &- 2*basis6(x) + basis11(x) + basis3_2(x)
     &*ll1px + basis3_1(x)*ll1px + basis3_8(x)
     &*ll1px + basis3_4(x)*ll1px 
     &+ basis3_7(x)*ll1px - basis3_5(x)
     &*ll1px - (pi**2*ll2*ll1px)/12d0 + (ll2**3
     &*ll1px)/6d0 - (pi**2*ll1mx*ll1px)/6d0 
     &+ (ll2*ll1mx**2*ll1px)/2d0 - (ll1mx**3
     &*ll1px)/6d0 + (ll1mx**2*llx*ll1px)/2d0 
     &+ (3*pi**2*ll1px**2)/8d0 - (basis2_3(x)
     &*ll1px**2)/2d0 + (basis2_2(x)*ll1px**2)/2d0 
     &- (basis2_1(x)*ll1px**2)/2d0 
     &- (3*ll2**2*ll1px**2)
     &/4d0 - (ll2*ll1mx*ll1px**2)/2d0 - ll1mx
     &*llx*ll1px**2 + ll2*ll1px**3 + (5*llx
     &*ll1px**3)/6d0 - (11*ll1px**4)/24d0 
     &- (ll1px*zeta3)/8d0

      case(23)                  !-1100

         ris = -pi**4/90d0 - 2*cli4pt5 - basis2(x)/2d0 
     &+ (3*basis1(x))/2d0 + basis15(x)/4d0 
     &+ basis4(x) + basis6(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis3_2(x)*llx - basis3_1(x)*llx 
     &- basis3_8(x)*llx + basis3_4(x)*llx 
     &- basis3_7(x)*llx - basis3_5(x)*llx 
     &- (pi**2*ll2*llx)/12d0 + (ll2**3*llx)/6d0 
     &+ (pi**2*ll1mx*llx)/6d0 - (ll2*ll1mx**2
     &*llx)/2d0 + (ll1mx**3*llx)/6d0 + (pi**2*llx**2)
     &/24d0 - (basis2_3(x)*llx**2)/2d0 - (ll2**2
     &*llx**2)/4d0 + (ll2*ll1mx*llx**2)/2d0 
     &- (ll1mx**2*llx**2)/2d0 + (pi**2*ll2*ll1px)
     &/6d0 - (ll2**3*ll1px)/3d0 + (pi**2*llx*ll1px)
     &/12d0 - (ll2**2*llx*ll1px)/2d0 +ll2*ll1mx
     &*llx*ll1px + (ll1mx*llx**2*ll1px)/2d0 
     &- (pi**2*ll1px**2)/6d0 + (ll2**2*ll1px**2)/2d0 
     &- (ll1mx*llx*ll1px**2)/2d0 - (ll2
     &*ll1px**3)/3d0 - (llx*ll1px**3)/6d0 
     &+ ll1px**4/6d0 + (7*llx*zeta3)/8d0 
     &        - (3*ll1px*zeta3)/4d0

      case(24)                  !-1101

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/144d0 - basis2_2(x)*basis2_1(x)
     &+ basis2_1(x)**2/2d0 
     &+ b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) 
     &+ basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0)
     &- 4*cli4pt5 + 2*basis8(x) - basis3(x) 
     &+ (3*basis2(x))/2d0 + basis1(x)/2d0 - 2*basis5(x) 
     &- basis12(x) + basis15(x)/4d0 
     &+ 2*basis4(x) - basis9(x) 
     &+ basis10(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis13(x)/4d0 
     &+ basis14(x)/4d0 - basis3_2(x)
     &*ll1mx + basis3_1(x)*ll1mx - basis3_8(x)
     &*ll1mx - basis3_4(x)*ll1mx 
     &+ basis3_7(x)*ll1mx + basis3_5(x)
     &*ll1mx + (pi**2*ll2*ll1mx)/4d0 - (ll2**3
     &*ll1mx)/2d0 - (5*pi**2*ll1mx**2)/48d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (3*ll2**2*ll1mx**2)/4d0-(ll2*ll1mx**3)/3d0 
     &+ (3*ll1mx**4)/32d0 + (ll1mx**3*llx)/4d0 
     &- 2*basis3_1(x)*ll1px + (pi**2*ll2*ll1px)/6d0 
     &- (ll2**3*ll1px)/3d0 - (3*pi**2*ll1mx*ll1px)
     &/8d0 + (ll2**2*ll1mx*ll1px)/2d0 - ll2
     &*ll1mx**2*ll1px + (ll1mx**3*ll1px)/24d0 
     &- (ll1mx**2*llx*ll1px)/4d0 
     &- (7*pi**2*ll1px**2)/48d0 + (ll2**2*ll1px**2)/2d0 
     &+ (9*ll1mx**2*ll1px**2)/16d0 - (ll1mx*llx
     &*ll1px**2)/4d0 - (ll2*ll1px**3)/3d0 
     &+ (ll1mx*ll1px**3)/24d0 - (llx*ll1px**3)
     &/12d0 + (17*ll1px**4)/96d0 - (ll1mx*zeta3)/8d0 
     &- (7*ll1px*zeta3)/4d0

      case(25)                  !-111-1

         ris = -pi**4/288d0 + (pi**2*basis2_3(x))/12d0 
     &- basis2_3(x)**2/2d0 + (pi**2*ll2**2)/24d0 
     &- (basis2_3(x)*ll2**2)/2d0 - ll2**4/8d0 
     &- (pi**2*ll2*ll1mx)/12d0 + basis2_3(x)*ll2
     &*ll1mx + (ll2**3*ll1mx)/2d0 - (ll2**2
     &*ll1mx**2)/2d0 - basis3_6(x)*ll1px - (pi**2
     &*ll2*ll1px)/12d0 + (ll2**3*ll1px)/6d0 
     &+ (pi**2*ll1mx*ll1px)/12d0 - (ll2**2*ll1mx
     &*ll1px)/2d0 + (ll2*ll1mx**2*ll1px)/2d0 
     &+ (7*ll1px*zeta3)/8d0

      case(26)                  !-1110

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/288d0 + basis2_2(x)*basis2_1(x)
     &- basis2_1(x)**2/2d0 
     &- (b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) 
     &+ basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0))
     &+ 5*cli4pt5 
     &- basis8(x) + 2*basis3(x) - basis2(x) - basis1(x) 
     &+ 2*basis12(x) - basis15(x)/2d0 
     &- 3*basis4(x) + basis9(x)/2d0 
     &- basis10(x)/2d0 - 2*basis6(x) 
     &+ 4*basis11(x) - 4*basis7(x) 
     &- basis13(x)/4d0 + 2*basis3_8(x)*ll1mx 
     &- (pi**2*ll2*ll1mx)/12d0 +(ll2**3*ll1mx)/6d0 
     &- (pi**2*ll1mx**2)/8d0 - (basis2_3(x)
     &*ll1mx**2)/2d0 + (basis2_2(x)*ll1mx**2)/2d0 
     &- (basis2_1(x)*ll1mx**2)/2d0
     &-(ll2**2*ll1mx**2)/2d0 
     &+ (2*ll2*ll1mx**3)/3d0 - (7*ll1mx**4)/24d0 
     &- basis3_6(x)*llx - (pi**2*ll2*llx)/12d0 
     &+ (ll2**3*llx)/6d0 + basis2_3(x)*ll1mx
     &*llx - (ll2*ll1mx**2*llx)/2d0 + (ll1mx**3
     &*llx)/3d0 + 2*basis3_1(x)*ll1px - (pi**2*ll2
     &*ll1px)/3d0 + (2*ll2**3*ll1px)/3d0 + (pi**2
     &*ll1mx*ll1px)/6d0 + (3*pi**2*ll1px**2)/8d0 
     &- ll2**2*ll1px**2 + (2*ll2*ll1px**3)/3d0 
     &+ (llx*ll1px**3)/3d0 - (3*ll1px**4)/8d0 
     &- (7*ll1mx*zeta3)/4d0 + (7*llx*zeta3)/8d0 
     &+ (13*ll1px*zeta3)/8d0

      case(27)                  !-1111

         ris = cli4pt5 - basis8(x) 
     &+ basis3_6(x)
     &*ll1mx - (basis2_3(x)*ll1mx**2)/2d0 + (ll2
     &*ll1mx**3)/6d0 - (ll1mx**3*ll1px)/6d0

      case(28)                  !0-1-1-1

         ris = pi**4/90d0 - basis4(x) - basis3_4(x)
     &*ll1px - (pi**2*ll1px**2)/12d0 - (basis2_2(x)
     &*ll1px**2)/2d0 - (llx*ll1px**3)/3d0 
     &+ ll1px**4/8d0

      case(29)                  !0-1-10

         ris = -basis2_2(x)**2/2d0 - basis3_4(x)*llx 
     &- (pi**2*llx*ll1px)/6d0 - basis2_2(x)*llx*ll1px 
     &- (llx**2*ll1px**2)/2d0 + (llx*ll1px**3)/6d0 
     &+ llx*zeta3

      case(30)                  !0-1-11

         bcflag_save=bcflag
         ris = -basis18(x)
         bcflag=bcflag_save

      case(31)                  !0-10-1

         ris = pi**4/45d0 + basis2_2(x)**2/2d0 - 2*basis2(x) 
     &- 2*basis4(x) - 2*basis6(x) + 2*basis3_2(x)
     &*ll1px + (pi**2*ll1px**2)/6d0 + (llx
     &*ll1px**3)/3d0 - ll1px**4/6d0 - 2*ll1px*zeta3

      case(32)                  !0-100

         ris = -3*basis2(x) + 2*basis3_2(x)*llx 
     &- (basis2_2(x)*llx**2)/2d0

      case(33)                  !0-101

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/180d0 - basis2_2(x)*basis2_1(x) - b16 
     &+ 2*basis3(x) + 2*basis2(x) + 2*basis1(x) - 2*basis5(x) 
     &- 2*basis4(x) - 2*basis6(x) - basis13(x)/2d0 
     &+ basis14(x)/2d0 - 2*basis3_2(x)*ll1mx 
     &- (pi**2*ll1mx**2)/8d0 - ll1mx**4/16d0 
     &+ (ll1mx**3*llx)/6d0 - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll1mx*ll1px)/12d0 + (ll1mx**3
     &*ll1px)/12d0 - (ll1mx**2*llx*ll1px)/2d0 
     &+ (5*pi**2*ll1px**2)/24d0 + (ll1mx**2*ll1px**2)
     &/8d0 - (ll1mx*llx*ll1px**2)/2d0 + (ll1mx
     &*ll1px**3)/12d0 + (llx*ll1px**3)/6d0 
     &- (7*ll1px**4)/48d0 - (3*ll1mx*zeta3)/2d0 
     &- (3*ll1px*zeta3)/2d0

      case(34)                  !0-11-1

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/90d0 + (pi**2*basis2_2(x))/12d0 - basis2_3(x)
     &*basis2_2(x) + basis2_2(x)**2/2d0 
     &+ b18 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0
     &+ 6*cli4pt5 - basis2(x)/2d0 
     &- (3*basis1(x))/2d0 - (3*basis15(x))/4d0 
     &- basis4(x) - 2*basis6(x) + 3*basis11(x) 
     &- 6*basis7(x) - (basis2_2(x)*ll2**2)/2d0 + basis2_2(x)
     &*ll2*ll1mx - basis3_6(x)*ll1px 
     &+ basis3_3(x)*ll1px - basis3_2(x)*ll1px 
     &+ 3*basis3_1(x)
     &*ll1px + 3*basis3_8(x)*ll1px 
     &+ 2*basis3_7(x)*ll1px - (7*pi**2*ll2
     &*ll1px)/12d0 + (7*ll2**3*ll1px)/6d0 
     &- (5*pi**2*ll1mx*ll1px)/12d0 - (ll2**2*ll1mx
     &*ll1px)/2d0 + (3*ll2*ll1mx**2*ll1px)/2d0 
     &- (ll1mx**3*ll1px)/2d0 + (3*ll1mx**2*llx
     &*ll1px)/2d0 + (17*pi**2*ll1px**2)/24d0 
     &- (basis2_3(x)*ll1px**2)/2d0 + (basis2_2(x)
     &*ll1px**2)/2d0 - (basis2_1(x)*ll1px**2)/2d0 
     &- (7*ll2**2*ll1px**2)/4d0 - (3*ll2*ll1mx
     &*ll1px**2)/2d0 - 2*ll1mx*llx*ll1px**2 
     &+ 2*ll2*ll1px**3 + (ll1mx*ll1px**3)/2d0 
     &+ (7*llx*ll1px**3)/6d0 - (19*ll1px**4)/24d0 
     &+ (17*ll1px*zeta3)/8d0

      case(35)                  !0-110

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/45d0 + basis2_2(x)*basis2_1(x) + b16 
     &+ 4*cli4pt5 + basis2(x) - 3*basis1(x) 
     &- basis15(x)/2d0 - 2*basis4(x) 
     &- 2*basis6(x) + 4*basis11(x) 
     &- 4*basis7(x) - basis3_6(x)*llx 
     &+ basis3_3(x)*llx - basis3_2(x)*llx 
     &+ basis3_1(x)*llx 
     &+ basis3_8(x)*llx - (pi**2*ll2*llx)/12d0 
     &+ (ll2**3*llx)/6d0 - (pi**2*ll1mx*llx)/12d0 
     &+ basis2_2(x)*ll1mx*llx 
     &- (ll2**2*ll1mx*llx)
     &/2d0 + (ll2*ll1mx**2*llx)/2d0 - (ll1mx**3
     &*llx)/6d0 + (ll1mx**2*llx**2)/2d0 + 2*basis3_1(x)
     &*ll1px - (pi**2*ll2*ll1px)/3d0 + (2*ll2**3
     &*ll1px)/3d0 + (pi**2*ll1px**2)/3d0 - ll2**2
     &*ll1px**2 + (2*ll2*ll1px**3)/3d0 + (llx
     &*ll1px**3)/3d0 - ll1px**4/3d0 - (llx*zeta3)/8d0 
     &+ (3*ll1px*zeta3)/2d0

      case(36)                  !0-111

         ris = -pi**4/72d0 - cli4pt5 - basis8(x) 
     &- basis3(x) - basis2(x)/2d0 + basis1(x)/2d0 
     &+ 2*basis5(x) - basis12(x) 
     &+ basis15(x)/4d0 + 2*basis4(x) 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_6(x)*ll1mx 
     &- basis3_3(x)*ll1mx + basis3_2(x)*ll1mx 
     &- basis3_1(x)*ll1mx - basis3_8(x)*ll1mx 
     &+ (3*pi**2*ll1mx**2)/16d0 - (basis2_2(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/4d0 - (ll2*ll1mx**3)/3d0 
     &+ (19*ll1mx**4)/96d0 - (7*ll1mx**3*llx)/12d0 
     &+ (pi**2*ll2*ll1px)/6d0 - (ll2**3*ll1px)/3d0 
     &- (pi**2*ll1mx*ll1px)/24d0 - (ll1mx**3
     &*ll1px)/24d0 + (ll1mx**2*llx*ll1px)/4d0 
     &- (13*pi**2*ll1px**2)/48d0 + (ll2**2*ll1px**2)/2d0 
     &- (ll1mx**2*ll1px**2)/16d0 + (ll1mx*llx
     &*ll1px**2)/4d0 - (ll2*ll1px**3)/3d0 - (ll1mx
     &*ll1px**3)/24d0 - (llx*ll1px**3)/4d0 
     &+ (23*ll1px**4)/96d0 + (7*ll1mx*zeta3)/4d0

      case(37)                  !00-1-1

         ris = -pi**4/90d0 + basis2(x) + basis4(x) 
     &+ basis6(x) - basis3_2(x)*ll1px 
     &- (pi**2*ll1px**2)/12d0 - (llx*ll1px**3)/6d0 
     &+ ll1px**4/12d0 + ll1px*zeta3

      case(38)                  !00-10

         ris = 3*basis2(x) - basis3_2(x)*llx

      case(39)                  !00-11

         ris = -pi**4/72d0 - 2*cli4pt5 - basis3(x) 
     &- (3*basis2(x))/2d0 + basis1(x)/2d0 + basis5(x) 
     &+ basis15(x)/4d0 + 2*basis4(x) 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 2*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_2(x)*ll1mx 
     &+ (pi**2*ll1mx**2)/16d0 + ll1mx**4/32d0 
     &- (ll1mx**3*llx)/12d0 + (pi**2*ll2*ll1px)/6d0 
     &- (ll2**3*ll1px)/3d0 - (pi**2*ll1mx
     &*ll1px)/24d0 - (ll1mx**3*ll1px)/24d0 
     &+ (ll1mx**2*llx*ll1px)/4d0 - (13*pi**2
     &*ll1px**2)/48d0 + (ll2**2*ll1px**2)/2d0 
     &- (ll1mx**2*ll1px**2)/16d0 + (ll1mx*llx
     &*ll1px**2)/4d0 - (ll2*ll1px**3)/3d0 
     &- (ll1mx*ll1px**3)/24d0 - (llx*ll1px**3)/4d0 
     &+ (23*ll1px**4)/96d0 + (3*ll1mx*zeta3)/4d0

      case(40)                  !000-1

         ris = -basis2(x)

      case(41)                  !0000

         ris = llx**4/24d0

      case(42)                  !0001

         ris = basis1(x)

      case(43)                  !001-1

         ris = pi**4/90d0 + 2*cli4pt5 + basis2(x)/2d0 
     &- (3*basis1(x))/2d0 - basis15(x)/4d0 
     &- basis4(x) - basis6(x) + 2*basis11(x) 
     &- 2*basis7(x) + basis3_1(x)*ll1px - (pi**2*ll2
     &*ll1px)/6d0 + (ll2**3*ll1px)/3d0 
     &+ (pi**2*ll1px**2)/6d0 - (ll2**2*ll1px**2)/2d0 
     &+ (ll2*ll1px**3)/3d0 + (llx*ll1px**3)/6d0 
     &- ll1px**4/6d0 + (3*ll1px*zeta3)/4d0

      case(44)                  !0010

         ris = -3*basis1(x) + basis3_1(x)*llx

      case(45)                  !0011
         
         ris = pi**4/90d0 - basis3(x) + basis1(x) + basis5(x) 
     &- basis3_1(x)*ll1mx + (pi**2*ll1mx**2)/12d0 
     &+ ll1mx**4/24d0 - (ll1mx**3*llx)/6d0 
     &+ ll1mx*zeta3

      case(46)                  !01-1-1
         ris = -pi**4/90d0 - 3*cli4pt5 - basis2(x)/2d0 
     &+ basis1(x)/2d0 + basis15(x)/4d0 + basis4(x) 
     &- basis11(x) + 3*basis7(x) + basis3_2(x)
     &*ll1px - basis3_1(x)*ll1px - basis3_8(x)
     &*ll1px + basis3_4(x)*ll1px 
     &- basis3_7(x)*ll1px - basis3_5(x)
     &*ll1px + (pi**2*ll2*ll1px)/6d0 - (ll2**3
     &*ll1px)/3d0 + (pi**2*ll1mx*ll1px)/6d0 
     &- (ll2*ll1mx**2*ll1px)/2d0 + (ll1mx**3
     &*ll1px)/6d0 - (ll1mx**2*llx*ll1px)/2d0 
     &- (pi**2*ll1px**2)/8d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ (ll2**2*ll1px**2)/4d0 + ll2*ll1mx
     &*ll1px**2 + ll1mx*llx*ll1px**2 
     &- (ll2*ll1px**3)/2d0 - (ll1mx*ll1px**3)/2d0 
     &- (llx*ll1px**3)/6d0 + ll1px**4/6d0 
     &- (3*ll1px*zeta3)/4d0

      case(47)                  !01-10

         bcflag_save=bcflag         
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/45d0 - b16 - 4*cli4pt5 
     &- basis2(x) 
     &+ 3*basis1(x) + basis15(x)/2d0 + 2*basis4(x) 
     &+ 2*basis6(x) - 4*basis11(x) 
     &+ 4*basis7(x) + basis3_2(x)*llx - basis3_1(x)*llx 
     &- basis3_8(x)*llx + basis3_4(x)*llx 
     &- basis3_7(x)*llx - basis3_5(x)*llx 
     &- (pi**2*ll2*llx)/12d0 + (ll2**3*llx)/6d0 
     &+ (pi**2*ll1mx*llx)/6d0 - (ll2*ll1mx**2
     &*llx)/2d0 + (ll1mx**3*llx)/6d0 - (ll1mx**2
     &*llx**2)/2d0 - 2*basis3_1(x)*ll1px + (pi**2*ll2
     &*ll1px)/3d0 - (2*ll2**3*ll1px)/3d0 + (pi**2
     &*llx*ll1px)/12d0 + basis2_1(x)*llx*ll1px 
     &- (ll2**2*llx*ll1px)/2d0 + ll2*ll1mx
     &*llx*ll1px + ll1mx*llx**2*ll1px 
     &- (pi**2*ll1px**2)/3d0 + ll2**2*ll1px**2 
     &- (ll1mx*llx*ll1px**2)/2d0 - (2*ll2
     &*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ ll1px**4/3d0 + (7*llx*zeta3)/8d0 
     &- (3*ll1px*zeta3)/2d0

      case(48)                  !01-11

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)


         ris = pi**4/72d0 + (pi**2*basis2_1(x))/12d0 - basis2_3(x)
     &*basis2_1(x) + basis2_2(x)*basis2_1(x) -basis2_1(x)**2/2d0 
     &- (b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0))
     &+ 6*cli4pt5 + 3*basis3(x) 
     &- basis2(x)/2d0 - (3*basis1(x))/2d0 - 2*basis5(x) 
     &+ 3*basis12(x) - (3*basis15(x))/4d0 
     &- 4*basis4(x) - 4*basis6(x) + 6*basis11(x) 
     &- 6*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 - (basis2_1(x)*ll2**2)/2d0 
     &- basis3_2(x)*ll1mx + basis3_1(x)*ll1mx 
     &+ 3*basis3_8(x)*ll1mx - basis3_4(x)
     &*ll1mx + basis3_7(x)*ll1mx 
     &+ basis3_5(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 + basis2_1(x)*ll2*ll1mx
     &-(ll2**3
     &*ll1mx)/6d0 - (17*pi**2*ll1mx**2)/48d0 
     &- (basis2_3(x)*ll1mx**2)/2d0 
     &+ (basis2_2(x)*ll1mx**2)/2d0 
     &-(basis2_1(x)*ll1mx**2)/2d0 
     &- (ll2**2*ll1mx**2)/4d0 + ll2*ll1mx**3 
     &- (47*ll1mx**4)/96d0 + (11*ll1mx**3*llx)/12d0 
     &+ 2*basis3_1(x)*ll1px - (pi**2*ll2*ll1px)/2d0 
     &+ ll2**3*ll1px - (pi**2*ll1mx*ll1px)/24d0 
     &- basis2_1(x)*ll1mx*ll1px + (ll2**2*ll1mx
     &*ll1px)/2d0 - ll2*ll1mx**2*ll1px 
     &+ (ll1mx**3*ll1px)/24d0 - (5*ll1mx**2*llx
     &*ll1px)/4d0 + (29*pi**2*ll1px**2)/48d0 
     &- (3*ll2**2*ll1px**2)/2d0 + (9*ll1mx**2
     &*ll1px**2)/16d0 - (ll1mx*llx*ll1px**2)/4d0 
     &+ ll2*ll1px**3 + (ll1mx*ll1px**3)/24d0 
     &+ (7*llx*ll1px**3)/12d0 - (55*ll1px**4)/96d0 
     &- (29*ll1mx*zeta3)/8d0 + (3*ll1px*zeta3)/2d0

      case(49)                  !010-1

         bcflag_save=bcflag
         ris = basis16(x)
         bcflag=bcflag_save

      case(50)                  !0100

         ris = 3*basis1(x) - 2*basis3_1(x)*llx 
     &+ (basis2_1(x)*llx**2)/2d0

      case(51)                  !0101

         ris = -pi**4/45d0 + basis2_1(x)**2/2d0 + 2*basis3(x) 
     &- 2*basis1(x) - 2*basis5(x) + 2*basis3_1(x)*ll1mx 
     &- (pi**2*ll1mx**2)/6d0 - ll1mx**4/12d0 
     &+ (ll1mx**3*llx)/3d0 - 2*ll1mx*zeta3

      case(52)                  !011-1

         bcflag_save=bcflag
         ris = basis17(x)
         bcflag=bcflag_save

      case(53)                  !0110

         ris = -basis2_1(x)**2/2d0 - basis3_3(x)*llx 
     &+ (pi**2*ll1mx*llx)/6d0 - basis2_1(x)*ll1mx*llx 
     &- (ll1mx**2*llx**2)/2d0 + llx*zeta3

      case(54)                  !0111

         ris = pi**4/90d0 - basis3(x) + basis3_3(x)*ll1mx 
     &- (pi**2*ll1mx**2)/12d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll1mx**3*llx)/3d0

      case(55)                  !1-1-1-1

         ris = cli4pt5 - basis7(x) 
     &+ basis3_5(x)
     &*ll1px - (pi**2*ll1px**2)/12d0 + (basis2_3(x)
     &*ll1px**2)/2d0 +(ll2**2*ll1px**2)/2d0 - (ll2
     &*ll1mx*ll1px**2)/2d0 - (ll2*ll1px**3)/3d0 
     &+ (ll1mx*ll1px**3)/3d0

      case(56)                  !1-1-10

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/480d0 + basis2_2(x)**2/2d0 
     &+ b18 
     &+ pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0
     &+ 3*cli4pt5 - basis3(x) 
     &- basis2(x) - basis1(x) - basis15(x)/2d0 
     &+ basis9(x)/2d0 - basis10(x)/2d0 
     &- 2*basis6(x) + 2*basis11(x) 
     &- 3*basis7(x) + basis13(x)/4d0 
     &+ basis3_5(x)*llx + (pi**2*ll2*llx)/12d0 
     &- (ll2**3*llx)/6d0 + 2*basis3_1(x)*ll1px 
     &+ 2*basis3_8(x)*ll1px + 2*basis3_7(x)
     &*ll1px - (pi**2*ll2*ll1px)/4d0 + (ll2**3
     &*ll1px)/2d0 - (pi**2*ll1mx*ll1px)/4d0 
     &+ ll2*ll1mx**2*ll1px - (ll1mx**3
     &*ll1px)/3d0 - (pi**2*llx*ll1px)/6d0 
     &+ basis2_3(x)*llx*ll1px + ll2**2*llx
     &*ll1px - ll2*ll1mx*llx*ll1px 
     &+ ll1mx**2*llx*ll1px + (5*pi**2*ll1px**2)
     &/12d0 - (basis2_3(x)*ll1px**2)/2d0 + (basis2_2(x)
     &*ll1px**2)/2d0 - (basis2_1(x)*ll1px**2)/2d0
     &-ll2**2
     &*ll1px**2 - (3*ll2*ll1mx*ll1px**2)/2d0 
     &- (ll2*llx*ll1px**2)/2d0 - ll1mx*llx
     &*ll1px**2 + (3*ll2*ll1px**3)/2d0 + (ll1mx
     &*ll1px**3)/2d0 + llx*ll1px**3 - (5*ll1px**4)
     &/8d0 - (ll1mx*zeta3)/8d0 - (7*llx*zeta3)/8d0 
     &+ (5*ll1px*zeta3)/4d0

      case(57)                  !1-1-11

         ris = -pi**4/288d0 + (pi**2*basis2_3(x))/12d0 
     &- basis2_3(x)**2/2d0 + (pi**2*ll2**2)/24d0 
     &- (basis2_3(x)*ll2**2)/2d0 - ll2**4/8d0 
     &- basis3_5(x)*ll1mx - (pi**2*ll2*ll1mx)/6d0 
     &+ basis2_3(x)*ll2*ll1mx + (2*ll2**3
     &*ll1mx)/3d0 - (ll2**2*ll1mx**2)/2d0 
     &+ (pi**2*ll1mx*ll1px)/6d0 - basis2_3(x)
     &*ll1mx*ll1px - ll2**2*ll1mx*ll1px 
     &+ ll2*ll1mx**2*ll1px + (ll2*ll1mx
     &*ll1px**2)/2d0 - (ll1mx**2*ll1px**2)/2d0 
     &+ (7*ll1mx*zeta3)/8d0

      case(58)                  !1-10-1

         bcflag_save=bcflag
         b18 = basis18(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = (11*pi**4)/720d0 - basis2_2(x)**2/2d0 
     &- b18 
     &- (pi**4/480d0 - (pi**2*basis2_2(x))/12d0 
     &+ basis2_3(x)*basis2_2(x) - basis2_2(x)**2/2d0 
     &- 3*cli4pt5 + basis3(x) 
     &+ basis2(x) + basis1(x) + basis15(x)/2d0 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ 2*basis6(x) - 2*basis11(x) 
     &+ 3*basis7(x) - basis13(x)/4d0 + (basis2_2(x)
     &*ll2**2)/2d0 + basis3_4(x)*ll1mx - basis2_2(x)
     &*ll2*ll1mx - 2*basis3_1(x)*ll1px 
     &- 2*basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px + (pi**2*ll2*ll1px)/4d0 - (ll2**3
     &*ll1px)/2d0 + (5*pi**2*ll1mx*ll1px)/12d0 
     &+ basis2_2(x)*ll1mx*ll1px - ll2*ll1mx**2
     &*ll1px + (ll1mx**3*ll1px)/3d0 - ll1mx**2
     &*llx*ll1px - (5*pi**2*ll1px**2)/12d0 
     &+ (basis2_3(x)*ll1px**2)/2d0 - (basis2_2(x)
     &*ll1px**2)/2d0 + (basis2_1(x)*ll1px**2)/2d0 
     &+ ll2**2*ll1px**2 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + 2*ll1mx*llx*ll1px**2 
     &- (3*ll2*ll1px**3)/2d0 - (2*ll1mx
     &*ll1px**3)/3d0 - llx*ll1px**3 
     &+ (5*ll1px**4)/8d0 - (7*ll1mx*zeta3)/8d0 
     &- (5*ll1px*zeta3)/4d0)
     &+ 2*basis3(x) + (3*basis2(x))/2d0 
     &+ basis1(x)/2d0 + basis15(x)/4d0 - basis4(x) 
     &- basis9(x) + basis10(x) 
     &+ 2*basis6(x) - basis11(x) - basis13(x)/2d0 
     &- basis3_6(x)*ll1px + basis3_3(x)*ll1px 
     &- basis3_2(x)*ll1px - basis3_1(x)*ll1px 
     &- basis3_8(x)*ll1px - 2*basis3_7(x)
     &*ll1px - (pi**2*ll2*ll1px)/12d0 + (ll2**3
     &*ll1px)/6d0 + (pi**2*ll1mx*ll1px)/12d0 
     &- (ll2**2*ll1mx*ll1px)/2d0 - (ll2
     &*ll1mx**2*ll1px)/2d0 + (ll1mx**3*ll1px)
     &/6d0 - (ll1mx**2*llx*ll1px)/2d0 - (pi**2
     &*ll1px**2)/8d0 + (basis2_3(x)*ll1px**2)/2d0 
     &- (basis2_2(x)*ll1px**2)/2d0 
     &+(basis2_1(x)*ll1px**2)/2d0 
     &+ (ll2**2*ll1px**2)/4d0 + (3*ll2*ll1mx
     &*ll1px**2)/2d0 + ll1mx*llx*ll1px**2 
     &- ll2*ll1px**3 - (ll1mx*ll1px**3)/2d0 
     &- (5*llx*ll1px**3)/6d0 + (11*ll1px**4)/24d0 
     &+ (ll1mx*zeta3)/4d0 - (3*ll1px*zeta3)/8d0

      case(59)                  !1-100

         ris = pi**4/72d0 + 2*cli4pt5 + basis3(x) 
     &+ (3*basis2(x))/2d0 - basis1(x)/2d0 - basis5(x) 
     &- basis15(x)/4d0 - 2*basis4(x) 
     &- 2*basis6(x) + 2*basis11(x) 
     &- 2*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 - (pi**2*ll1mx**2)/16d0 
     &- ll1mx**4/32d0 - basis3_6(x)*llx + basis3_3(x)
     &*llx - basis3_2(x)*llx + basis3_1(x)*llx 
     &+ basis3_8(x)*llx - (pi**2*ll2*llx)/12d0 
     &+ (ll2**3*llx)/6d0 - (pi**2*ll1mx*llx)/12d0 
     &- (ll2**2*ll1mx*llx)/2d0 + (ll2*ll1mx**2
     &*llx)/2d0 - (ll1mx**3*llx)/12d0 - (pi**2*llx**2)
     &/24d0 + (basis2_3(x)*llx**2)/2d0 + (ll2**2
     &*llx**2)/4d0 - (ll2*ll1mx*llx**2)/2d0 
     &+ (ll1mx**2*llx**2)/2d0 - (pi**2*ll2*ll1px)
     &/6d0 + (ll2**3*ll1px)/3d0 + (pi**2*ll1mx
     &*ll1px)/24d0 + (ll1mx**3*ll1px)/24d0 
     &- (ll1mx**2*llx*ll1px)/4d0 + (13*pi**2
     &*ll1px**2)/48d0 - (ll2**2*ll1px**2)/2d0 
     &+ (ll1mx**2*ll1px**2)/16d0 - (ll1mx*llx
     &*ll1px**2)/4d0 + (ll2*ll1px**3)/3d0 + (ll1mx
     &*ll1px**3)/24d0 + (llx*ll1px**3)/4d0 
     &- (23*ll1px**4)/96d0 - (3*ll1mx*zeta3)/4d0 
     &- (llx*zeta3)/8d0

      case(60)                  !1-101

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)


         ris = -pi**4/72d0 + basis2_2(x)*basis2_1(x)
     &- basis2_1(x)**2/2d0 
     &-(b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) 
     &+ basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0))
     &+ 4*cli4pt5 
     &- 2*basis8(x) + basis3(x) - (3*basis2(x))/2d0 
     &- basis1(x)/2d0 + 2*basis5(x) + basis12(x) 
     &- basis15(x)/4d0 + 2*basis11(x) 
     &- 2*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 + basis3_6(x)*ll1mx 
     &- basis3_3(x)*ll1mx + basis3_2(x)*ll1mx 
     &- basis3_1(x)
     &*ll1mx + basis3_8(x)*ll1mx - (pi**2*ll2
     &*ll1mx)/12d0 + (ll2**3*ll1mx)/6d0 + (5*pi**2
     &*ll1mx**2)/48d0 - (basis2_3(x)*ll1mx**2)/2d0 
     &+ (basis2_2(x)*ll1mx**2)/2d0 
     &-(basis2_1(x)*ll1mx**2)/2d0 
     &- (ll2**2*ll1mx**2)/4d0 + (ll2*ll1mx**3)/3d0 
     &- (3*ll1mx**4)/32d0 - (ll1mx**3*llx)/4d0 
     &+ 2*basis3_1(x)*ll1px - (pi**2*ll2*ll1px)/6d0 
     &+ (ll2**3*ll1px)/3d0 - (pi**2*ll1mx*ll1px)
     &/24d0 - (ll1mx**3*ll1px)/24d0 + (ll1mx**2
     &*llx*ll1px)/4d0 + (pi**2*ll1px**2)/16d0 
     &- (ll2**2*ll1px**2)/2d0 - (ll1mx**2
     &*ll1px**2)/16d0 + (ll1mx*llx*ll1px**2)/4d0 
     &+ (ll2*ll1px**3)/3d0 - (ll1mx*ll1px**3)/24d0 
     &+ (llx*ll1px**3)/12d0 - (3*ll1px**4)/32d0 
     &+ (5*ll1mx*zeta3)/8d0 + (3*ll1px*zeta3)/2d0

      case(61)                  !1-11-1

         ris = (-23*pi**4)/1440d0 - (pi**2*basis2_3(x))/12d0 
     &+ basis2_3(x)**2/2d0 - 2*basis8(x) 
     &- 2*basis10(x) + 2*basis7(x) 
     &- (pi**2*ll2**2)/24d0 + (basis2_3(x)*ll2**2)/2d0 
     &+ ll2**4/8d0 - (pi**2*ll2*ll1mx)/12d0 
     &- basis2_3(x)*ll2*ll1mx - (ll2**3
     &*ll1mx)/6d0 + 2*basis3_6(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &+ ll2**2*ll1mx*ll1px - (pi**2*ll1px**2)/6d0 
     &+ (ll2**2*ll1px**2)/2d0 - ll2*ll1mx
     &*ll1px**2 + (ll1mx*ll1px**3)/3d0 
     &- ll1px**4/12d0 + (ll1mx*zeta3)/4d0 
     &- 2*ll1px*zeta3

      case(62)                  !1-110

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/72d0 - basis2_2(x)*basis2_1(x)
     &+ basis2_1(x)**2/2d0 
     &+ b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0) - 6*cli4pt5 
     &- 3*basis3(x) + basis2(x)/2d0 + (3*basis1(x))/2d0 
     &+ 2*basis5(x) - 3*basis12(x) 
     &+ (3*basis15(x))/4d0 + 4*basis4(x) 
     &+ 4*basis6(x) - 6*basis11(x) 
     &+ 6*basis7(x) + basis13(x)/4d0 
     &- basis14(x)/4d0 - 2*basis3_8(x)
     &*ll1mx + (3*pi**2*ll1mx**2)/16d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/4d0 - (ll2*ll1mx**3)/2d0 
     &+ (31*ll1mx**4)/96d0 + 2*basis3_6(x)*llx 
     &+ (pi**2*ll2*llx)/6d0 - (ll2**3*llx)/3d0 
     &- (pi**2*ll1mx*llx)/12d0 - basis2_3(x)
     &*ll1mx*llx + (ll2**2*ll1mx*llx)/2d0 
     &- (5*ll1mx**3*llx)/12d0 - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/2d0 - ll2**3*ll1px 
     &- (pi**2*ll1mx*ll1px)/24d0 - (ll1mx**3
     &*ll1px)/24d0 + (ll1mx**2*llx*ll1px)/4d0 
     &- (29*pi**2*ll1px**2)/48d0 + (3*ll2**2
     &*ll1px**2)/2d0 - (ll1mx**2*ll1px**2)/16d0 
     &+ (ll1mx*llx*ll1px**2)/4d0 - ll2
     &*ll1px**3 - (ll1mx*ll1px**3)/24d0 
     &- (7*llx*ll1px**3)/12d0 + (55*ll1px**4)/96d0 
     &+ (11*ll1mx*zeta3)/4d0 - (7*llx*zeta3)/4d0 
     &- (3*ll1px*zeta3)/2d0

      case(63)                  !1-111

         ris = -3*cli4pt5 + 3*basis8(x) 
     &- 2*basis3_6(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - (ll2**3*ll1mx)/6d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 
     &- (7*ll1mx*zeta3)/8d0

      case(64)                  !10-1-1

         ris = -pi**4/480d0 - basis3(x) + basis9(x)/2d0 
     &- basis10(x)/2d0 + basis13(x)/4d0 
     &+ basis3_6(x)*ll1px - basis3_3(x)*ll1px 
     &- basis3_4(x)*ll1px + basis3_7(x)
     &*ll1px + basis3_5(x)*ll1px + (pi**2*ll2
     &*ll1px)/6d0 - (ll2**3*ll1px)/3d0 + (ll2**2
     &*ll1mx*ll1px)/2d0 - (pi**2*ll1px**2)/6d0 
     &- (basis2_1(x)*ll1px**2)/2d0
     &+(ll2**2*ll1px**2)/2d0 
     &- ll2*ll1mx*ll1px**2 - (ll1mx*llx
     &*ll1px**2)/2d0 + (ll1mx*ll1px**3)/2d0 
     &- (ll1mx*zeta3)/8d0 - (ll1px*zeta3)/8d0

      case(65)                  !10-10

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/180d0 + b16 - 2*basis3(x) 
     &- 2*basis2(x) - 2*basis1(x) + 2*basis5(x) 
     &+ 2*basis4(x) + 2*basis6(x) + basis13(x)/2d0 
     &- basis14(x)/2d0 + (pi**2*ll1mx**2)/8d0 
     &+ ll1mx**4/16d0 + basis3_6(x)*llx - basis3_3(x)
     &*llx - basis3_4(x)*llx + basis3_7(x)
     &*llx + basis3_5(x)*llx + (pi**2*ll2*llx)
     &/6d0 - (ll2**3*llx)/3d0 - (pi**2*ll1mx*llx)/12d0 
     &+ (ll2**2*ll1mx*llx)/2d0 - (ll1mx**3
     &*llx)/6d0 + 2*basis3_1(x)*ll1px - (pi**2*ll1mx
     &*ll1px)/12d0 - (ll1mx**3*ll1px)/12d0 
     &- (pi**2*llx*ll1px)/12d0 - basis2_1(x)*llx*ll1px 
     &+ (ll2**2*llx*ll1px)/2d0 - ll2*ll1mx
     &*llx*ll1px + (ll1mx**2*llx*ll1px)/2d0 
     &- ll1mx*llx**2*ll1px - (5*pi**2*ll1px**2)
     &/24d0 - (ll1mx**2*ll1px**2)/8d0 + ll1mx*llx
     &*ll1px**2 - (ll1mx*ll1px**3)/12d0 - (llx
     &*ll1px**3)/6d0 + (7*ll1px**4)/48d0 + (3*ll1mx
     &*zeta3)/2d0 - (3*llx*zeta3)/4d0 + (3*ll1px*zeta3)/2d0

      case(66)                  !10-11

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/72d0 - (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &+ b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0)
     &- 4*cli4pt5 + 2*basis8(x)
     & - basis3(x) + (3*basis2(x))/2d0 + basis1(x)/2d0 
     &- 2*basis5(x) - basis12(x) 
     &+ basis15(x)/4d0 - 2*basis11(x) 
     &+ 2*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- basis3_6(x)*ll1mx + basis3_3(x)*ll1mx 
     &- 2*basis3_8(x)*ll1mx + basis3_4(x)
     &*ll1mx - basis3_7(x)*ll1mx 
     &- basis3_5(x)*ll1mx - basis2_1(x)*ll2*ll1mx 
     &+ (pi**2*ll1mx**2)/16d0 + (basis2_3(x)
     &*ll1mx**2)/2d0 - (basis2_2(x)*ll1mx**2)/2d0 
     &+ (basis2_1(x)*ll1mx**2)/2d0
     &+(ll2**2*ll1mx**2)/4d0 
     &- (5*ll2*ll1mx**3)/6d0 + (25*ll1mx**4)/96d0 
     &- (ll1mx**3*llx)/4d0 - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/6d0 - (ll2**3*ll1px)/3d0 
     &+ (pi**2*ll1mx*ll1px)/8d0 + basis2_1(x)*ll1mx
     &*ll1px - (ll2**2*ll1mx*ll1px)/2d0 
     &+ ll2*ll1mx**2*ll1px + (ll1mx**3
     &*ll1px)/24d0 + (3*ll1mx**2*llx*ll1px)/4d0 
     &- (pi**2*ll1px**2)/16d0 + (ll2**2*ll1px**2)/2d0 
     &- (7*ll1mx**2*ll1px**2)/16d0 - (ll1mx*llx
     &*ll1px**2)/4d0 - (ll2*ll1px**3)/3d0 
     &+ (ll1mx*ll1px**3)/24d0 - (llx*ll1px**3)
     &/12d0 + (3*ll1px**4)/32d0 + (ll1mx*zeta3)/4d0 
     &- (3*ll1px*zeta3)/2d0

      case(67)                  !100-1

         bcflag_save=bcflag
         b16 = basis16(x)
         bcflag=bcflag_save         

         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = pi**4/360d0 - b16 + basis3(x) + basis2(x) 
     &+ basis1(x) - basis5(x) - basis4(x) 
     &- basis6(x) - basis13(x)/4d0 + basis14(x)
     &/4d0 - (pi**2*ll1mx**2)/16d0 - ll1mx**4/32d0 
     &+ (ll1mx**3*llx)/12d0 - basis3_1(x)*ll1px 
     &+ (pi**2*ll1mx*ll1px)/24d0 + (ll1mx**3
     &*ll1px)/24d0 - (ll1mx**2*llx*ll1px)/4d0 
     &+ (5*pi**2*ll1px**2)/48d0 + (ll1mx**2*ll1px**2)
     &/16d0 - (ll1mx*llx*ll1px**2)/4d0 + (ll1mx
     &*ll1px**3)/24d0 + (llx*ll1px**3)/12d0 
     &- (7*ll1px**4)/96d0 - (3*ll1mx*zeta3)/4d0 
     &- (3*ll1px*zeta3)/4d0

      case(68)                  !1000

         ris = -basis1(x)+basis3_1(x)*llx
     &-(basis2_1(x)*llx**2)/2d0 
     &- (ll1mx*llx**3)/6d0

      case(69)                  !1001

         ris = -basis2_1(x)**2/2d0 - basis3_1(x)*ll1mx

      case(70)                  !101-1

         bcflag_save=bcflag
         b17= basis17(x)
         bcflag=bcflag_save
         
         llx = log(x)
         ll1mx = log(1d0-x)
         ll1px = log(1d0+x)

         ris = -pi**4/144d0 + (pi**2*basis2_1(x))/12d0 - basis2_3(x)
     &*basis2_1(x) + basis2_2(x)*basis2_1(x) -basis2_1(x)**2/2d0 
     &+ 4*cli4pt5 
     &- (b17-( -pi**4/288d0 
     &- (pi**2*basis2_1(x))/12d0 + basis2_3(x)
     &*basis2_1(x) - basis2_2(x)*basis2_1(x) +basis2_1(x)**2/2d0 
     &- 5*cli4pt5 
     &+ basis8(x) - 2*basis3(x) + basis2(x) + basis1(x) 
     &- 2*basis12(x) + basis15(x)/2d0 
     &+ 3*basis4(x) - basis9(x)/2d0 
     &+ basis10(x)/2d0 + 2*basis6(x) 
     &- 4*basis11(x) + 4*basis7(x) 
     &+ basis13(x)/4d0 + (basis2_1(x)*ll2**2)/2d0 
     &- 2*basis3_8(x)*ll1mx + (pi**2*ll2
     &*ll1mx)/12d0 - basis2_1(x)*ll2*ll1mx 
     &- (ll2**3*ll1mx)/6d0 + (pi**2*ll1mx**2)/8d0 
     &+ (basis2_3(x)*ll1mx**2)/2d0 - (basis2_2(x)
     &*ll1mx**2)/2d0 + (basis2_1(x)*ll1mx**2)/2d0 
     &+ (ll2**2*ll1mx**2)/2d0-(2*ll2*ll1mx**3)/3d0 
     &+ (7*ll1mx**4)/24d0 - (ll1mx**3*llx)/3d0 
     &- basis3_3(x)*ll1px - 2*basis3_1(x)*ll1px 
     &+ (pi**2*ll2*ll1px)/3d0-(2*ll2**3*ll1px)/3d0 
     &- (3*pi**2*ll1px**2)/8d0 + ll2**2*ll1px**2 
     &- (2*ll2*ll1px**3)/3d0 - (llx*ll1px**3)/3d0 
     &+ (3*ll1px**4)/8d0 + (7*ll1mx*zeta3)/4d0 
     &- (5*ll1px*zeta3)/8d0))
     &- 2*basis8(x) + basis3(x) - (3*basis2(x))/2d0 
     &- basis1(x)/2d0 + 2*basis5(x) + basis12(x) 
     &- basis15(x)/4d0 - 2*basis4(x) 
     &+ basis9(x) - basis10(x) 
     &+ 2*basis11(x) - 2*basis7(x) 
     &- basis13(x)/4d0 - basis14(x)/4d0 
     &- (basis2_1(x)*ll2**2)/2d0 + 2*basis3_8(x)*ll1mx 
     &-(pi**2*ll2*ll1mx)/6d0
     &+basis2_1(x)*ll2*ll1mx 
     &+ (ll2**3*ll1mx)/3d0 - (pi**2*ll1mx**2)/16d0 
     &- (basis2_3(x)*ll1mx**2)/2d0 + (basis2_2(x)
     &*ll1mx**2)/2d0 - (basis2_1(x)*ll1mx**2)/2d0 
     &- (3*ll2**2*ll1mx**2)/4d0 + (5*ll2*ll1mx**3)
     &/6d0 - (25*ll1mx**4)/96d0 + (ll1mx**3*llx)/4d0 
     &+ 2*basis3_3(x)*ll1px + 2*basis3_1(x)*ll1px - (pi**2
     &*ll2*ll1px)/6d0 + (ll2**3*ll1px)/3d0 
     &- (pi**2*ll1mx*ll1px)/24d0 - (ll1mx**3
     &*ll1px)/24d0 + (ll1mx**2*llx*ll1px)/4d0 
     &+ (7*pi**2*ll1px**2)/48d0 - (ll2**2*ll1px**2)/2d0 
     &- (ll1mx**2*ll1px**2)/16d0 + (ll1mx*llx
     &*ll1px**2)/4d0 + (ll2*ll1px**3)/3d0 - (ll1mx
     &*ll1px**3)/24d0 + (llx*ll1px**3)/12d0 
     &- (17*ll1px**4)/96d0 - (3*ll1mx*zeta3)/4d0 
     &- (ll1px*zeta3)/4d0

      case(71)                  !1010

         ris = pi**4/45d0 + basis2_1(x)**2/2d0 - 2*basis3(x) 
     &+ 2*basis1(x) + 2*basis5(x) + (pi**2*ll1mx**2)/6d0 
     &+ ll1mx**4/12d0 + 2*basis3_3(x)*llx - (pi**2
     &*ll1mx*llx)/3d0 + basis2_1(x)*ll1mx*llx 
     &- (ll1mx**3*llx)/3d0 + ll1mx**2*llx**2 
     &+ 2*ll1mx*zeta3 - 2*llx*zeta3

      case(72)                  !1011

         ris = -pi**4/30d0 + 3*basis3(x) - 2*basis3_3(x)*ll1mx 
     &+ (pi**2*ll1mx**2)/12d0 - (basis2_1(x)*ll1mx**2)/2d0 
     &- (ll1mx**3*llx)/2d0 - ll1mx*zeta3

      case(73)                  !11-1-1

         ris = (7*pi**4)/720d0 + basis8(x) 
     &+ basis10(x) - basis7(x) + (pi**2*ll2
     &*ll1mx)/12d0 - (ll2**3*ll1mx)/6d0 + (ll2**2
     &*ll1mx**2)/4d0 - basis3_6(x)*ll1px 
     &- (pi**2*ll2*ll1px)/6d0 + (ll2**3*ll1px)/3d0 
     &- (ll2**2*ll1mx*ll1px)/2d0 + (pi**2
     &*ll1px**2)/12d0 - (ll2**2*ll1px**2)/4d0 
     &+ (ll2*ll1mx*ll1px**2)/2d0 - (ll1mx
     &*ll1px**3)/6d0 + ll1px**4/24d0 - (ll1mx*zeta3)
     &/8d0 + ll1px*zeta3

      case(74)                  !11-10

         ris = pi**4/72d0 + cli4pt5 + basis8(x) 
     &+ basis3(x) + basis2(x)/2d0 - basis1(x)/2d0 
     &- 2*basis5(x) + basis12(x) 
     &- basis15(x)/4d0 - 2*basis4(x) 
     &- 2*basis6(x) + 2*basis11(x) 
     &- 2*basis7(x) - basis13(x)/4d0 
     &+ basis14(x)/4d0 + (pi**2*ll2*ll1mx)/12d0 
     &- (ll2**3*ll1mx)/6d0 - (5*pi**2*ll1mx**2)/48d0 
     &+ (ll2**2*ll1mx**2)/4d0 - (ll2*ll1mx**3)/6d0 
     &- ll1mx**4/32d0 - basis3_6(x)*llx - (pi**2
     &*ll2*llx)/12d0 + (ll2**3*llx)/6d0 + (pi**2
     &*ll1mx*llx)/12d0 - (ll2**2*ll1mx*llx)/2d0 
     &+ (ll2*ll1mx**2*llx)/2d0 + (ll1mx**3*llx)
     &/12d0 - (pi**2*ll2*ll1px)/6d0+(ll2**3*ll1px)
     &/3d0 + (pi**2*ll1mx*ll1px)/24d0 + (ll1mx**3
     &*ll1px)/24d0 - (ll1mx**2*llx*ll1px)/4d0 
     &+ (13*pi**2*ll1px**2)/48d0 - (ll2**2*ll1px**2)/2d0 
     &+ (ll1mx**2*ll1px**2)/16d0 - (ll1mx*llx
     &*ll1px**2)/4d0 + (ll2*ll1px**3)/3d0 + (ll1mx
     &*ll1px**3)/24d0 + (llx*ll1px**3)/4d0 
     &- (23*ll1px**4)/96d0 - (13*ll1mx*zeta3)/8d0 
     &+ (7*llx*zeta3)/8d0

      case(75)                  !11-11
         
         ris = 3*cli4pt5 - 3*basis8(x) 
     &+ basis3_6(x)*ll1mx - (pi**2*ll2*ll1mx)/6d0 
     &+ (ll2**3*ll1mx)/3d0 + (pi**2*ll1mx**2)/24d0 
     &- (ll2**2*ll1mx**2)/4d0 + (7*ll1mx*zeta3)/4d0

      case(76)                  !110-1

         ris = -pi**4/288d0 + basis4(x) 
     &- basis9(x)/2d0 + basis10(x)/2d0 
     &+ basis13(x)/4d0 + (pi**2*ll1mx**2)/24d0 
     &- basis3_3(x)*ll1px - (pi**2*ll1px**2)/24d0 
     &+ ll1px**4/24d0 + (5*ll1mx*zeta3)/8d0 
     &+ (7*ll1px*zeta3)/8d0

      case(77)                  !1100

         ris = -pi**4/90d0 + basis3(x) - basis1(x) - basis5(x) 
     &- (pi**2*ll1mx**2)/12d0 - ll1mx**4/24d0 - basis3_3(x)
     &*llx + (pi**2*ll1mx*llx)/6d0 + (ll1mx**3
     &*llx)/6d0 - (ll1mx**2*llx**2)/4d0 
     &- ll1mx*zeta3 + llx*zeta3

      case(78)                  !1101

         ris = pi**4/30d0 - 3*basis3(x) + basis3_3(x)*ll1mx 
     &+ (pi**2*ll1mx**2)/12d0 + 2*ll1mx*zeta3

      case(79)                  !111-1

         ris = -cli4pt5 + basis8(x) 
     &+ (pi**2*ll2*ll1mx)/12d0 -(ll2**3*ll1mx)/6d0 
     &- (pi**2*ll1mx**2)/24d0 + (ll2**2*ll1mx**2)/4d0 
     &- (ll2*ll1mx**3)/6d0 - (7*ll1mx*zeta3)/8d0

      case(80)                  !1110

         ris = -pi**4/90d0 + basis3(x) - (pi**2*ll1mx**2)/12d0 
     &- ll1mx*zeta3

      case(81)                  !1111

         ris = ll1mx**4/24d0
         
      end select
      
c --- set the imaginary part back to zero if it has been modified to
c --- get the branch cuts right (and should be zero). Also, set imaginary
c --- part of result to zero if x is between 0 and 1.

      if (bcflag.eq.1) then
         x = x - dcmplx(0d0,1d-60)
         xre = dreal(x)
         if (xre.ge.0d0.and.xre.le.1d0) then
            ris = dcmplx(dreal(ris),0d0)
         endif
      endif

      HPL4else=ris
      return
      end function
