      subroutine runF_0(m0sq,Gr,Czero0,N0)
C---  Expression for rearrangement of Eq. 5.65
C---  Calculates D0, requires D00
C---  Small terms of order Gr(i,j)*Dij
C---  Denominator m0sq
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
      integer ep,N0,m,n,np
      parameter(np=3)
      real(dp):: m0sq,Gr(np,np)
      complex(dp):: Czero0(-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(dii(z2(n,m))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dd00+N0,ep-1)
      
      Dv(dd0+N0,ep)=
     . (8d0*Dv(dd00+N0,ep)
     . -pole
     . -2d0*Czero0(ep)
     . +bit)/(2d0*m0sq)

      enddo
      
      return
      end
