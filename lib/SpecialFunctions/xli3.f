 
c================================================================ 
      FUNCTION xli3(z)
c
c T. RIEMANN, 06/2008, see also Mathematica notebook Li3.nb
c TRANSFORMS THE ARGUMENT OF trilogarithm WITH ARBITRARY COMPLEX ARGUMENT 
c TO THE REGION |Z|<1, Real(Z)<0.5 AND CALLS THEN A SERIES
c BASED ON SOME FORMULAE OF LAST CENTURIES TAKEN FROM MATHWORLD:
c http://mathworld.wolfram.com/Trilogarithm.html
c  
c THE Real ARGUMENT Z WITH Z>1 IS FORBIDDEN - THERE IS THE CUT
c (AND WE WOULD HAVE LOG OF NEGATIVE REAL ARGUMENT HERE)
c
c CALLS THE SERIES xcli3 
c 

      implicit none
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      complex(dp)::xli3
      integer::Num
      real(dp)::REZ,AIZ,AAZ,AAS,AIS,RE1,RES
      complex(dp)::Z,xcli3,cli3
      EXTERNAL xcli3
c
c ... Num = number of terms in sum xcli3(Num,z) for Li_3(z)
c     Num=12 gives about 12 safe digits
c     Num=20 gives about 18 safe digits
c
      Num=20      
c
c      F1=zeta2
c
      REZ=REAL(Z)
      AIZ=AIMAG(Z)
      AAZ=ABS(Z)
c
      IF(AIZ) 70,71,70       
 71   IF(REZ-1._dp) 70,70,72
 72   WRITE(*,*) 'FROM xli3(z): z is real with z>1, forbidden'
      stop 
 70   CONTINUE       
      IF(AAZ) 11,9,11
 9    cli3=0
      GOTO 99
 11   IF(AAZ-1._dp) 6,4,1
 1    RE1=REAL(1._dp/Z)
      IF(RE1-0.5_dp) 3,3,2
 2    CONTINUE
c real(1/z)>1/2, |z|>1, this is my case 22
      cli3=-xcli3(Num,1._dp-1._dp/z)-xcli3(Num,1._dp-z)
     & +LOG(1._dp/z)*zeta2 
     &   - 0.5_dp*(LOG(1._dp/z))**2 * LOG(1._dp-1._dp/z)
     &   + (LOG(1._dp/z))**3/6._dp - LOG(-z)*zeta2 - (LOG(-z))**3/6._dp 
     &      + zeta3
      GOTO 99
 3    CONTINUE
      cli3=xcli3(Num,1._dp/z) - zeta2 *LOG(-z) -(LOG(-z))**3/6._dp
      GOTO 99
 4    IF(REZ-1._dp) 6,5,1
 5    cli3=CMPLX(zeta3,0._dp,kind=dp)
      GOTO 99
 6    IF(REZ-0.5_dp) 7,7,8
 7    CONTINUE
      cli3=xcli3(Num,z)
      GOTO 99
 8    CONTINUE
c |z|<1, Real(z)>0.5, this is my case 12 
      cli3=-xcli3(Num,1._dp-1._dp/z) - xcli3(Num,1._dp-z) + LOG(z)*zeta2 
     & - 0.5_dp*(LOG(z))**2 * LOG(1._dp-z) + (LOG(z))**3/6._dp 
     &      + zeta3

 99   CONTINUE
      AAS= ABS(cli3)
      RES=REAL(cli3)
      AIS=AIMAG(cli3)
c      IF(JPRINT) 97,97,98
c 98   CONTINUE
c 97   CONTINUE
      xli3=cli3
      END

c================================================================ 

