C=============================================================================
C---  HPLs  of Rank 2  
C=============================================================================

      double  complex function HPL2else(n1,n2,x)
      implicit none 
      double precision pi,ll2,xre
      double complex x, ris,myi,ll1x,ll1mx,llx
      double complex ccli2
      integer n1,n2,j,bcflag

      pi=3.1415926535897932385D0
      myi=dcmplx(0d0,1d0)

      ll2 = dlog(2d0)
      j = 3*(n1+1) + (n2+1) +1

      ris=dcmplx(0d0,0d0)
      bcflag = 0

c---  +i*epsilon to get branch cuts right ---
      if (dimag(x).eq.0d0) then
         x = x + dcmplx(0d0,1d-60)
         bcflag = 1
      endif
c---  
      ll1x = log(1d0+x)
      ll1mx = log(1d0-x)
      llx = log(x)
      
      select case(j)
      case(1)
         ris=ll1x**2/2d0
      case(2)
         ris=ccli2(-x) + llx*ll1x
      case(3)
         ris= pi**2/12d0 - ll2**2/2d0 + ll2*ll1mx 
     &        - ll1mx*ll1x - ccli2((1d0-x)/2d0)
      case(4)
         ris=-ccli2(-x)
      case(5)
         ris=llx**2/2d0
      case(6)
         ris=ccli2(x)
      case(7)
         ris=-pi**2/12d0 + ccli2((1d0-x)/2d0) + ll2**2/2d0 
     &        -ll2*ll1mx
      case(8)
         ris=-ccli2(x)-ll1mx*llx
      case(9)
         ris=ll1mx**2/2d0
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

      HPL2else=ris 
      return
      end
