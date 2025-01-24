C=============================================================================
C---  basis functions
C=============================================================================
c---  Li2

      double complex  function ccli2(z)
      implicit none
      double complex ris, z, bsli2_inside,bsli2_outside, wcli2
      double complex zlocal
      double precision zabs, pi, zeta2, border, tiny, arg

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0

      border = 0.3d0 
      tiny = 1d-14
      zabs = abs(z)
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=-wcli2(1d0/z)-zeta2-0.5d0*log(-z)**2
      elseif (zabs.le.border) then 
         ris=bsli2_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli2_outside(zlocal)
      endif

      ccli2=ris
      return
      end
      
c---  recursion
      
      double complex  function wcli2(z)
      implicit none
      double complex z, ccli2
      wcli2 =  ccli2(z)
      return
      end

c--- Li3

      double complex  function cli3(z)
      implicit none
      double complex ris, z, bsli3_inside,bsli3_outside, wcli3
      double complex zlocal
      double precision zabs,border, pi, zeta2, zeta3,tiny,arg
      
      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
     
      border = 0.3d0
      zabs = abs(z)
      tiny = 1d-14
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=wcli3(1d0/z)-log(-z)**3/6d0-zeta2*log(-z)
      elseif (zabs.le.border) then 
         ris=bsli3_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli3_outside(zlocal)
      endif

      cli3=ris
      return
      end
      
c---  recursion

      double complex  function wcli3(z)
      implicit none
      double complex z, cli3
      wcli3 =  cli3(z)
      return
      end

c--- Li4

      double complex  function cli4(z)
      implicit none
      double complex ris, z, bsli4_outside, bsli4_inside, wcli4
      double complex zlocal
      double precision zabs, pi, zeta2, zeta3, zeta4, border,tiny,arg
     
      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      
      border = 0.3d0
      zabs = abs(z)
      tiny = 1d-14
      zlocal=z

      if (zabs.gt.1d0+tiny) then
         ris=-wcli4(1d0/z) -log(-z)**4/24d0 - 7d0*zeta4/4d0 
     &           - zeta2*log(-z)**2/2d0
      elseif (zabs.le.border) then 
         ris=bsli4_inside(z)
      else
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         ris=bsli4_outside(zlocal)
      endif

      
      cli4=ris
      return
      end
      
c     --- recursion for li4
      
      double complex  function wcli4(z)
      implicit none
      double complex z, cli4
      wcli4 = cli4(z)
      return
      end
c --- the case Li4(1-z^2) needs some special treatment because of its branch cut structure
c --- (that's what 'sbc' stands for: special branch cut)

      double complex function cli4_sbc(z)
      implicit none
      double complex ris, z, cli4, myi,basis14
      double complex ll1,ll2,ll3
      double precision pi,zabs,zreal
      integer s
      
      pi=3.1415926535897932385D0
      zabs = abs(z)
      zreal = dreal(z)
      myi = dcmplx(0d0,1d0)
           
      if (zabs.le.1d0) then !normal li4
         if (zreal.gt.0d0) then
            ris = cli4(1d0 - z**2)
         else if (zreal.eq.0d0 .and. s(z).eq.1) then !also normal li4
            ris = cli4(1d0 - z**2 - dcmplx(0d0,1d-60))
         else                   ! special branch cut configuration
            ris = cli4(1d0 - z**2- dcmplx(0d0,1d-60))
     &           - myi*pi*s(z)/3d0*(log(1d0 - z)+log(1d0+z))**3 
         endif
      else 
         ll1=log(1d0/z)
         ll2=log(1d0 - 1d0/z)
         ll3=log(1d0 + 1d0/z)
         ris = -2d0/3d0*ll1**4 + 4d0/3d0*ll2
     &        *ll1**3 + 4d0/3d0*ll3*ll1**3 
     &        - ll2**2*ll1**2 - ll3**2
     &        *ll1**2 - 2d0*ll2*ll3
     &        *ll1**2 - pi**2/3d0*ll1**2 + 1d0/3d0
     &        *ll2**3*ll1 + 1d0/3d0
     &        *ll3**3*ll1 + ll2
     &        *ll3**2*ll1 + pi**2/3d0*ll1
     &        *ll2 + ll3*ll2**2
     &        *ll1 + pi**2/3d0*ll1*ll3 
     &        - 1d0/24d0*ll2**4 - 1d0/24d0
     &        *ll3**4 - 1d0/6d0*ll2
     &        *ll3**3 - pi**2/12d0*ll2**2 
     &        - pi**2/12d0*ll3**2 
     &        - 1d0/4d0*ll2**2*ll3**2 
     &        - 1d0/6d0*ll2**3*ll3 
     &        - pi**2/6d0*ll2*ll3 
     &        - 7*pi**4/360d0 
     &        -basis14(1d0/z)
      endif
      
      cli4_sbc = ris
      return
      end

c --- the case Li4(4z/(1+z)^2) also needs some special treatment because of its branch cut structure

      double complex function cli4_sbc_2(z)
      implicit none
      double complex ris, z, cli4, myi, wcli4_sbc_2
      double complex arg,llx,ll1px,wcli4sbc2,zlocal
      double precision pi,zabs,zreal,ll2,tiny,arg2
      integer s
      
      pi=3.1415926535897932385D0
      ll2 = dlog(2d0)
      zabs = abs(z)
      zreal = dreal(z)
      myi = dcmplx(0d0,1d0)
      llx = 1
      ll1px = 1
      wcli4sbc2 = 1
      tiny = 1d-14
      zlocal=z
     
      ris = dcmplx(0d0,0d0)

      if (zabs.lt.1d0) then
         ris = cli4(4d0*z/(1d0+z)**2)
      elseif (zabs.lt.1d0+tiny) then 
         arg2=atan2(dimag(zlocal),dreal(zlocal))
         zlocal=dcmplx(cos(arg2),sin(arg2))
         arg = dcmplx(dreal(4d0*zlocal/(1d0+zlocal)**2),s(zlocal)*1d-60)
         ris = cli4(arg)
      else    
         wcli4sbc2 =  wcli4_sbc_2(1d0/z)
         llx = log(1d0/z)    
         ll1px = log(1d0+1d0/z)          
         ris = wcli4sbc2 + myi*pi*s(z)*
     &(4d0*ll2**2*llx - 8d0*ll2**2* ll1px 
     &+ 2d0*ll2*llx**2 - 8d0*ll2*llx
     &* ll1px + 8d0*ll2*ll1px**2 
     &+ 1d0/3d0*llx**3 - 2d0* ll1px*llx**2 
     &+ 4d0* ll1px**2*llx - 8d0/3d0* ll1px**3 
     &+ 8d0/3d0*ll2**3)


c             ris = wcli4_sbc_2(1d0/z) + myi*pi*s(z)*
c     &(4d0*ll2**2*log(1d0/z) - 8d0*ll2**2*log(1d0+1d0/z) 
c     &+ 2d0*ll2*log(1d0/z)**2 - 8d0*ll2*log(1d0/z)
c     &*log(1d0+1d0/z) + 8d0*ll2*log(1d0+1d0/z)**2 
c     &+ 1d0/3d0*log(1d0/z)**3 - 2d0*log(1d0+1d0/z)*log(1d0/z)**2 
c     &+ 4d0*log(1d0+1d0/z)**2*log(1d0/z) - 8d0/3d0*log(1d0+1d0/z)**3 
c     &+ 8d0/3d0*ll2**3)
      endif
      
      cli4_sbc_2 = ris
      return
      end

c     --- recursion for cli4_sbc_2
      
      double complex  function wcli4_sbc_2(z)
      implicit none
      double complex z, cli4_sbc_2
      wcli4_sbc_2 =  cli4_sbc_2(z)
      return
      end

C-----------------------------------------------------------------------
C     mapping of H_2-2(z) into convergent region
      
      double complex  function ch2m2(z)
      implicit none
      double complex ris,z,bsh2m2_inside,bsh2m2_outside,cli4,ccli2
      double complex HPL4,wch2m2,myi,zlocal !,cli4_sbc,cli3
      double precision pi,zeta2,zeta3,zeta4,zabs,zreal,border,tiny,arg
      integer s

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      myi = dcmplx(0d0,1d0)
      tiny=1d-14
      
      border = 0.3d0
      zabs = abs(z)
      zreal = dreal(z)
      zlocal=z


      if (zabs.lt.border) then ! inside circle of |z| = 0.3, we employ the log(1+z) expansion
         ris = bsh2m2_inside(z)            
      elseif (zabs.lt.1d0+tiny) then
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         if (zreal.ge.0d0) then ! on the half annulus 0.3 < |z| < 1 ; Re(z) >= 0, we have the log(x) exp.
            ris = bsh2m2_outside(zlocal)
         else                ! for Re(z) < 0, we map back to Re(z) > 0 by using the fact that HPL4(n1,n2,n3,n4,z) = (+-) HPL4(-n1,-n2,-n3,-n4,-z) (if n4 =/= 0):
            ris = HPL4(0,-1,0,1,-zlocal) 
         endif
      else                      ! For |z| > 1, we use the inversion formula to map into the unit circle. 
         ris = dcmplx(0d0,0d0) +
     &        wch2m2(1d0/z) 
     &        + 37d0*pi**4/720d0 
     &        - HPL4(0,1,0,0,1d0/z) 
     &        - log(1d0/z)**4/24d0 
     &        - pi**2/12d0*log(1d0/z)**2 
     &        - pi**2/6d0*ccli2(1d0/z) 
     &        - cli4(-1d0/z) 
     &        + 3d0*zeta3*log(1d0/z)/2d0 
     &        - pi**3*myi*s(z)*log(1d0/z)/12d0
      endif
      
      ch2m2=ris
      return
      end
      
      
c     --- recursion for H_2-2(z)
      
      double complex  function wch2m2(z)
      implicit none
      double complex z, ch2m2
           
      wch2m2 = ch2m2(z)
      return
      end

C------------------------------------------------------------------------------
C     mapping of H21-1(z) into convergent region 
      
      double complex  function ch21m1(z)
      implicit none
      double complex ris,z,bsh21m1_inside,bsh21m1_outside_1,zlocal
      double complex bsh21m1_outside_2,cli4,ccli2,HPL4,wch21m1,myi,ch2m2
      double precision pi,zeta2,zeta3,zeta4,border,zreal,zabs,ll2,tiny
      double precision arg
      integer s

      pi=3.1415926535897932385D0
      zeta2=pi**2/6d0
      zeta3=1.20205690315959428539973816151d0
      zeta4=pi**4/90d0
      ll2 = dlog(2d0)
      border = 0.3d0
      myi = dcmplx(0d0,1d0)
      tiny=1d-14

      zabs = abs(z)
      zreal = dreal(z)
      zlocal=z


      if (zabs.lt.border) then ! inside circle of |z| = 0.3, we employ the log(1+z) expansion
         ris = bsh21m1_inside(z)           
      elseif (zabs.lt.1d0+tiny) then
         if (zabs.gt.1d0) then
            arg=atan2(dimag(zlocal),dreal(zlocal))
            zlocal=dcmplx(cos(arg),sin(arg))
         endif
         if (zreal.ge.0d0) then ! on the half annulus 0.3 < |z| < 1 ; Re(z) >= 0, we have the log(x) exp.
            ris = bsh21m1_outside_1(zlocal)
         else                ! for Re(z) < 0, we map back to Re(z) > 0 by using the fact that HPL4(n1,n2,n3,n4,z) = (+-) HPL4(-n1,-n2,-n3,-n4,-z) (if n4 =/= 0):
            ris = bsh21m1_outside_2(zlocal)
         endif
      else                      ! For |z| > 1, we use the inversion formula to map into the unit circle. 
         ris = -wch21m1(1d0/z)-pi**4/144d0 -ch2m2(1d0/z) 
     &        + log(1d0/z)**4/24d0 
     &        + pi**2*ll2**2/3d0 - ll2**4/12d0 
     &        + 3d0*pi**2*ll2*log(1d0/z)/4d0 
     &        + pi**2*log(1d0/z)**2/8d0 
     &        + pi**2*ccli2(1d0/z)/4d0 
     &        - 2*cli4(dcmplx(0.5d0,0d0)) 
     &        + cli4(-1d0/z) 
     &        - 7d0*zeta3*log(1d0/z)/8d0 
     &        + myi*s(z)*(pi**3*ll2/6d0 
     &        + pi**3*log(1d0/z)/12d0 
     &        - 0.5d0*pi*ll2**2*log(1d0/z) 
     &        - 0.5d0*pi*ll2*log(1d0/z)**2 
     &        - pi*ll2*ccli2(1d0/z))
     &        - HPL4(0,0,1,-1,1d0/z) 
     &        + HPL4(0,0,1,0,1d0/z) 
     &        + HPL4(0,1,0,0,1d0/z) 
     &        + HPL4(0,1,1,0,1d0/z)
      endif

      ch21m1=ris
      return
      end
      

c     --- recursion for H21-1(z)
      
      double complex  function wch21m1(z)
      implicit none
      double complex z, ch21m1
            
      wch21m1 =  ch21m1(z)
      return
      end
