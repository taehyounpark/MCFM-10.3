C=============================================================================
C---  Auxiliary functions
C=============================================================================

c ------------------------------------------------
c --- some auxiliary functions:
      integer function s(x)
      implicit none
      double complex x
      integer ris
      ris = 1
      if (dimag(x).lt.0d0) then
         ris = -1
      endif
      s=ris
      return
      end function
c ------------------------------------------------
c$$$c s1(x) = s(x^2)...if the imaginary part of x^2 is zero, it depends
c$$$c on the value of the real part of x if we return -1 or 1 because:
c$$$c (x+I*eps) = x^2 + 2*x*I*eps - eps^2
c$$$      integer function s1(x)
c$$$      implicit none
c$$$      double complex x
c$$$      integer ris
c$$$      ris = 1
c$$$      if (dimag(x**2).lt.0d0) then
c$$$         ris = -1
c$$$      else if (dimag(x**2).eq.0d0.and.dreal(x).lt.0d0) then
c$$$         ris = -1
c$$$      endif
c$$$      s1=ris
c$$$      return
c$$$      end function
c$$$c ------------------------------------------------
c$$$c s2(x) = s(4x/(1+x)^2)...if the imaginary part of the arg is zero, it depends
c$$$c on the value of the real part of x if we return -1 or 1 because:
c$$$c 4*(x+I*eps)/(1+x+I*eps)^2 = 4x/(1+x)^2 + 4*(1-x)*I*eps/(1+x)^3
c$$$      integer function s2(x)
c$$$      implicit none
c$$$      double complex x,arg
c$$$      integer ris
c$$$      arg = 4d0*x/(1d0+x)**2
c$$$      ris = 1
c$$$      if (dimag(arg).lt.0d0) then
c$$$         ris = -1
c$$$      else if (dimag(arg).eq.0d0.and.dreal(1d0-x).lt.0d0) then
c$$$         ris = -1
c$$$      endif
c$$$      s2=ris
c$$$      return
c$$$      end function
c$$$c ------------------------------------------------
