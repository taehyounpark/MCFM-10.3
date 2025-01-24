      double complex function HPL3else(n1, n2, n3, x)
      implicit none
      double precision pi, zeta2, zeta3,ll2,xre
      double complex x, ris
      double complex basis3_1,basis3_2,basis3_3,basis3_4
      double complex basis3_5,basis3_6,basis3_7,basis3_8
      double complex basis2_1,basis2_2,basis2_3
      double complex ccli2,cli3
      double complex ll1px,ll1mx,llx
      integer n1,n2,n3,j,bcflag

      pi=3.1415926535897932385D0
      zeta3=1.20205690315959428539973816151d0
      zeta2=pi**2/6d0
      ll2 = dlog(2d0)

      j=1+(n3+1)*1+(n2+1)*3+(n1+1)*9

      ris=dcmplx(0d0,0d0)
      bcflag = 0

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
c     #####################################################
      ll1px = log(1d0+x)
      ll1mx = log(1d0-x)
      llx = log(x)
      basis3_1 = cli3(x) 
      basis3_2 = cli3(-x)
      basis3_3 = cli3(1d0-x)
      basis3_4 = cli3(1d0/(1d0+x)) 
      basis3_5 = cli3((1d0+x)/2d0) 
      basis3_6 = cli3((1d0-x)/2d0) 
      basis3_7 = cli3((1d0-x)/(1d0+x)) 
      basis3_8 = cli3(2d0*x/(x-1d0))
      basis2_1 = ccli2(x)
      basis2_2 = ccli2(-x) 
      basis2_3 = ccli2((1d0-x)/2d0)
c     #####################################################

      select case(j)
      case(1)
         ris = ll1px**3/6d0
      case(2)
         ris = -(pi**2*ll1px)/6d0 + ll1px**3/6d0 
     &- basis3_4 
     &        + zeta3
      case(3)
         ris = (pi**2*ll2)/12d0 - ll2**3/6d0 
     &- (pi**2*ll1px)/12d0 + (ll2**2*ll1px)/2d0 
     &- (ll2*ll1px**2)/2d0 + basis3_5 - (7*zeta3)/8d0
      case(4)
         ris = (pi**2*ll1px)/3d0 + llx*ll1px**2 
     &- ll1px**3/3d0 + ll1px*basis2_2+2*basis3_4
     &-2*zeta3
      case(5)
         ris = (llx**2*ll1px)/2d0+llx*basis2_2
     &-basis3_2
      case(6)
         ris = (pi**2*ll2)/6d0 - ll2**3/3d0 
     &- (pi**2*ll1mx)/12d0 + (ll2**2*ll1mx)/2d0 
     &- (pi**2*ll1px)/12d0 + (ll2**2*ll1px)/2d0 
     &- ll2*ll1mx*ll1px - ll1mx*llx*ll1px 
     &+(ll1mx*ll1px**2)/2d0-ll1mx*basis2_2
     &+basis3_6 
     &- basis3_3-basis3_4+ basis3_7 + basis3_5
     &-(3*zeta3)/4d0
      case(7)
         ris = -(pi**2*ll2)/6d0 + ll2**3/3d0 
     &+ (pi**2*ll1px)/4d0 - (3*ll2**2*ll1px)/2d0 
     &+ ll2*ll1mx*ll1px + ll2*ll1px**2 
     &- ll1mx*ll1px**2 - ll1px*basis2_3 
     &- 2*basis3_5 + (7*zeta3)/4d0
      case(8)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &+ (pi**2*ll1mx)/6d0 - (ll2*ll1mx**2)/2d0 
     &+ ll1mx**3/6d0+(pi**2*llx)/12d0- (ll2**2*llx)/2d0 
     &+ ll2*ll1mx*llx - (ll1mx**2*llx)/2d0 
     &+ (pi**2*ll1px)/12d0 - (ll2**2*ll1px)/2d0 
     &+ ll2*ll1mx*ll1px - (ll1mx*ll1px**2)/2d0 
     &- llx*basis2_3 + basis3_2-basis3_1-basis3_8 
     &+ basis3_4 - basis3_7 - basis3_5 + (7*zeta3)/8d0
      case(9)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &- (ll2*ll1mx**2)/2d0 + (ll1mx**2*ll1px)/2d0 
     &+ ll1mx*basis2_3 - basis3_6 + (7*zeta3)/8d0
      case(10)
         ris = -(pi**2*ll1px)/6d0 - (llx*ll1px**2)/2d0 
     &+ ll1px**3/6d0 - ll1px*basis2_2 - basis3_4 
     &+ zeta3
      case(11)
         ris = -(llx*basis2_2) + 2*basis3_2
      case(12)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &- (pi**2*ll1mx)/12d0 - (ll2**2*ll1mx)/2d0 
     &+ (ll2*ll1mx**2)/2d0 - ll1mx**3/6d0 
     &+ (ll1mx**2*llx)/2d0 + ll1mx*basis2_2 
     &- basis3_6 + basis3_3 - basis3_2 +basis3_1
     &+basis3_8 
     &- zeta3/8d0
      case(13)
         ris = -basis3_2
      case(14)
         ris = llx**3/6d0
      case(15)
         ris = basis3_1
      case(16)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &+ (pi**2*ll1mx)/6d0 - (ll2*ll1mx**2)/2d0 
     &+ ll1mx**3/6d0 - (ll1mx**2*llx)/2d0 
     &+ (pi**2*ll1px)/12d0 - (ll2**2*ll1px)/2d0 
     &+ ll2*ll1mx*ll1px + ll1mx*llx*ll1px 
     &-(ll1mx*ll1px**2)/2d0+ll1px*basis2_1
     &+basis3_2 
     &- basis3_1 - basis3_8 + basis3_4 -basis3_7
     &-basis3_5 
     &+ (7*zeta3)/8d0
      case(17)
         ris = llx*basis2_1 - 2*basis3_1
      case(18)
         ris = (pi**2*ll1mx)/6d0 - (ll1mx**2*llx)/2d0 
     &- ll1mx*basis2_1 - basis3_3 + zeta3
      case(19)
         ris = (pi**2*ll2)/12d0 - ll2**3/6d0 
     &- (pi**2*ll1px)/6d0 + ll2**2*ll1px 
     &- ll2*ll1mx*ll1px - (ll2*ll1px**2)/2d0 
     &+ (ll1mx*ll1px**2)/2d0 + ll1px*basis2_3 
     &+ basis3_5 - (7*zeta3)/8d0
      case(20)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &- (pi**2*ll1mx)/12d0 - (ll2**2*ll1mx)/2d0 
     &+ (ll2*ll1mx**2)/2d0 - ll1mx**3/6d0 
     &- (pi**2*llx)/12d0 + (ll2**2*llx)/2d0 
     &- ll2*ll1mx*llx + (ll1mx**2*llx)/2d0 
     &+ llx*basis2_3 - basis3_6 + basis3_3 
     &- basis3_2 + basis3_1 + basis3_8 - zeta3/8d0
      case(21)
         ris = (pi**2*ll2)/6d0 - ll2**3/3d0 
     &- (pi**2*ll1mx)/12d0 + (ll2**2*ll1mx)/2d0 
     &- ll1mx*basis2_3 + 2*basis3_6 - (7*zeta3)/4d0
      case(22)
         ris = (pi**2*ll2)/6d0 - ll2**3/3d0 
     &- (pi**2*ll1mx)/12d0 + (ll2**2*ll1mx)/2d0 
     &- (pi**2*ll1px)/12d0 + (ll2**2*ll1px)/2d0 
     &- ll2*ll1mx*ll1px - ll1mx*llx*ll1px 
     &+(ll1mx*ll1px**2)/2d0-ll1px*basis2_1
     &+basis3_6 
     &- basis3_3 - basis3_4+basis3_7+basis3_5
     &- (3*zeta3)/4d0
      case(23)
         ris =-(ll1mx*llx**2)/2d0-llx*basis2_1
     &+basis3_1
      case(24)
         ris = -(pi**2*ll1mx)/3d0 + ll1mx**2*llx 
     &+ ll1mx*basis2_1 + 2*basis3_3 - 2*zeta3
      case(25)
         ris = -(pi**2*ll2)/12d0 + ll2**3/6d0 
     &+ (pi**2*ll1mx)/12d0 - (ll2**2*ll1mx)/2d0 
     &+ (ll2*ll1mx**2)/2d0 - basis3_6 + (7*zeta3)/8d0
      case(26)
         ris = (pi**2*ll1mx)/6d0 - basis3_3 + zeta3
      case(27)
         ris = -ll1mx**3/6d0
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
      
      HPL3else=ris
      return
      end
