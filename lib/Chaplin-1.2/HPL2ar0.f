c ------------------------------------------------
      double complex function HPL2ar0(n1,n2,x)
      implicit none
      integer n1,n2,j,bcflag
      double complex x,ris,myi,llx
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


         ris = (x**2)/2d0 - (x**3)/2d0 + (11d0*x**4)
     &        /24d0 - (5d0*x**5)/12d0 + (137d0*x**6)/360d0 - (7d0*x**
     &        7)/20d0 + (363d0*x**8)/1120d0 - (761d0*x**9)/2520d0 + (
     &        7129d0*x**10)/25200d0
         
      case(2)                   !-10
         
         llx = log(x)
         
         ris = -x + (x**2)/4d0 - (x**3)/9d0 + (x**4)
     &        /16d0 - (x**5)/25d0 + (x**6)/36d0 - (x**7)/49d0 + (x**8
     &        )/64d0 - (x**9)/81d0 + (x**10)/100d0 + x*llx - (x**2*ll
     &        x)/2d0 + (x**3*llx)/3d0 - (x**4*llx)/4d0 + (x**5*llx)/5
     &        d0 - (x**6*llx)/6d0 + (x**7*llx)/7d0 - (x**8*llx)/8d0 +
     &        (x**9*llx)/9d0 - (x**10*llx)/10d0
         
      case(3)                   !-11
         
         
         ris = (x**2)/2d0 - (x**3)/6d0 + (5d0*x**4)/
     &        24d0 - (7d0*x**5)/60d0 + (47d0*x**6)/360d0 - (37d0*x**7
     &        )/420d0 + (319d0*x**8)/3360d0 - (533d0*x**9)/7560d0 + (
     &        1879d0*x**10)/25200d0
         
      case(4)                   !0-1
         
         
         ris = x - (x**2)/4d0 + (x**3)/9d0 - (x**4)/
     &        16d0 + (x**5)/25d0 - (x**6)/36d0 + (x**7)/49d0 - (x**8)
     &        /64d0 + (x**9)/81d0 - (x**10)/100d0
         
      case(5)                   !00
         
         llx = log(x)
         
         ris = (llx**2)/2d0
         
      case(6)                   !01
         
         
         ris = x + (x**2)/4d0 + (x**3)/9d0 + (x**4)/
     &        16d0 + (x**5)/25d0 + (x**6)/36d0 + (x**7)/49d0 + (x**8)
     &        /64d0 + (x**9)/81d0 + (x**10)/100d0
         
      case(7)                   !1-1
         
         
         ris = (x**2)/2d0 + (x**3)/6d0 + (5d0*x**4)/
     &        24d0 + (7d0*x**5)/60d0 + (47d0*x**6)/360d0 + (37d0*x**7
     &        )/420d0 + (319d0*x**8)/3360d0 + (533d0*x**9)/7560d0 + (
     &        1879d0*x**10)/25200d0
         
      case(8)                   !10
         
         llx = log(x)
         
         ris = -x - (x**2)/4d0 - (x**3)/9d0 - (x**4)
     &        /16d0 - (x**5)/25d0 - (x**6)/36d0 - (x**7)/49d0 - (x**8
     &        )/64d0 - (x**9)/81d0 - (x**10)/100d0 + x*llx + (x**2*ll
     &        x)/2d0 + (x**3*llx)/3d0 + (x**4*llx)/4d0 + (x**5*llx)/5
     &        d0 + (x**6*llx)/6d0 + (x**7*llx)/7d0 + (x**8*llx)/8d0 +
     &        (x**9*llx)/9d0 + (x**10*llx)/10d0
         
      case(9)                   !11
         
         
         ris = (x**2)/2d0 + (x**3)/2d0 + (11d0*x**4)
     &        /24d0 + (5d0*x**5)/12d0 + (137d0*x**6)/360d0 + (7d0*x**
     &        7)/20d0 + (363d0*x**8)/1120d0 + (761d0*x**9)/2520d0 + (
     &        7129d0*x**10)/25200d0
c     End of expansions around x = 0
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
      
      HPL2ar0=ris
      return
      end function
