C=============================================================================
C---  basis mappings
C=============================================================================

c ---------------------------------------------------------
      double complex function basis2_1(x)
      implicit none
      double complex x,ccli2            
      basis2_1=ccli2(x)
      return
      end
c ---------------------------------------------------------
      double complex function basis2_2(x)
      implicit none
      double complex x,ccli2
      basis2_2=ccli2(-x)
      return
      end
c ---------------------------------------------------------
      double complex function basis2_3(x)
      implicit none
      double complex x,ccli2
      basis2_3=ccli2((1d0-x)/2d0)
      return
      end
c ---------------------------------------------------------
c ---------------------------------------------------------
c     basis3_1(z) = cli3(z) 
      double complex function basis3_1(x)
      implicit none
      double complex x,cli3           
      basis3_1=cli3(x)
      return
      end
c ---------------------------------------------------------
c     basis3_2(z) = cli3(-z)
      double complex function basis3_2(x)
      implicit none
      double complex x,cli3
      basis3_2=cli3(-x)
      return
      end
c ---------------------------------------------------------
c     basis3_3(z) = cli3(1-z)
      double complex function basis3_3(x)
      implicit none
      double complex x,cli3
      basis3_3 = cli3(1d0-x)
      return
      end
c ---------------------------------------------------------
c     basis3_4(z) = cli3(1/(1+z)) 
      double complex function basis3_4(x)
      implicit none
      double complex x,cli3
      basis3_4 = cli3(1d0/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis3_5(z) = cli3((1+z)/2) 
      double complex function basis3_5(x)
      implicit none
      double complex x,cli3
      basis3_5 = cli3((1d0+x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis3_6(z) = cli3((1-z)/2) 
      double complex function basis3_6(x)
      implicit none
      double complex x,cli3
      basis3_6 = cli3((1d0-x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis3_7(z) = cli3((1-z)/(1+z)) 
      double complex function basis3_7(x)
      implicit none
      double complex x,cli3
      basis3_7 = cli3((1d0-x)/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis3_8(z) = cli3(2z/(z-1))
      double complex function basis3_8(x)
      implicit none
      double complex x,cli3
      basis3_8 = cli3(2d0*x/(x-1d0))
      return
      end
c ---------------------------------------------------------
c ---------------------------------------------------------
c     basis1(x) = cli4(x) 
      double complex function basis1(x)
      implicit none
      double complex x,cli4
      basis1=cli4(x)
      return
      end
c ---------------------------------------------------------
c     basis2(x) = cli4(-x)
      double complex function basis2(x)
      implicit none
      double complex x,cli4
      basis2=cli4(-x)
      return
      end
c ---------------------------------------------------------
c     basis3(x) = cli4(1-x)
      double complex function basis3(x)
      implicit none
      double complex x,cli4
      basis3 = cli4(1d0-x)
      return
      end
c ---------------------------------------------------------
c     basis4(x) = cli4(1/(1+x)) 
      double complex function basis4(x)
      implicit none
      double complex x,cli4
      basis4 = cli4(1d0/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis5(x) = cli4(x/(x-1))
      double complex function basis5(x)
      implicit none
      double complex x,cli4
      basis5 = cli4(x/(x-1d0))
      return
      end
c ---------------------------------------------------------
c     basis6(x) = cli4(x/(x+1))
      double complex function basis6(x)
      implicit none
      double complex x,cli4
      basis6 = cli4(x/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis7(x) = cli4((1+x)/2) 
      double complex function basis7(x)
      implicit none
      double complex x,cli4
      basis7 = cli4((1d0+x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis8(x) = cli4((1-x)/2)
      double complex function basis8(x)
      implicit none
      double complex x,cli4
      basis8 = cli4((1d0-x)/2d0)
      return
      end
c ---------------------------------------------------------
c     basis9(x) = cli4((1-x)/(1+x))
      double complex function basis9(x)
      implicit none
      double complex x,cli4
      basis9 = cli4((1d0-x)/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis10(x) = cli4((x-1)/(x+1))
      double complex function basis10(x)
      implicit none
      double complex x,cli4
      basis10 = cli4((x-1d0)/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis11(x) = cli4(2x/(1+x))
      double complex function basis11(x)
      implicit none
      double complex x,cli4
      basis11 = cli4(2d0*x/(1d0+x))
      return
      end
c ---------------------------------------------------------
c     basis12(x) = cli4(2x/(x-1)) 
      double complex function basis12(x)
      implicit none
      double complex x,cli4
      basis12 = cli4(2d0*x/(x-1d0))
      return
      end
c ---------------------------------------------------------
c     basis13(x) = cli4(1-x^2) = cli4_sbc 
      double complex function basis13(x)
      implicit none
      double complex x,cli4_sbc !,cli4
      basis13=cli4_sbc(x)
      return
      end
c ---------------------------------------------------------
c     basis14(x) = cli4(x^2/(x^2-1)) 
      double complex function basis14(x)
      implicit none
      double complex x,cli4
      basis14 = cli4(x**2/(x**2-1d0))
      return
      end
c ---------------------------------------------------------
c     basis15(x) = cli4(4x/(1+x)^2) = cli4_sbc_2  
      double complex function basis15(x)
      implicit none      
      double complex x,cli4_sbc_2
      basis15=cli4_sbc_2(x)
      return
      end
c ---------------------------------------------------------
c     basis16(x) = ch2m2(x) 
      double complex function basis16(x)
      implicit none
      double complex x,ch2m2
      basis16=ch2m2(x)
      return
      end
c ---------------------------------------------------------
c     basis17(x) = ch21m1(x) 
      double complex function basis17(x)
      implicit none
      double complex x,ch21m1
      basis17=ch21m1(x)
      return
      end
c ---------------------------------------------------------
c     basis18(x) = ch21m1(-x)
      double complex function basis18(x)
      implicit none
      double complex x,ch21m1
      basis18=ch21m1(-x)
      return
      end
c ---------------------------------------------------------
