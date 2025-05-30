      subroutine runF_iiii(i1,i2,i3,i4,m0sq,Gr,Czero4,N0)
C---  Expression for rearrangement of extension of Eq. 5.69
C---  Calculates Diiii, requires D00iiii
C---  Small terms of order Gr(i,j)*Dijklmn
C---  Denominator m0sq
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
      integer ep,N0,i1,i2,i3,i4,m,n,np
      parameter(np=3)
      real(dp):: m0sq,Gr(np,np)
      complex(dp):: Czero4(z4max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(diiiiii(z6(n,m,i1,i2,i3,i4))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dzziiii(z4(i1,i2,i3,i4))+N0,ep-1)
      
      Dv(diiii(z4(i1,i2,i3,i4))+N0,ep)=
     . (24d0*Dv(dzziiii(z4(i1,i2,i3,i4))+N0,ep)
     . -pole
     . -2d0*Czero4(z4(i1,i2,i3,i4),ep)
     . +bit)/(2d0*m0sq)

      enddo
      
      return
      end
