      subroutine runF_iii(i1,i2,i3,m0sq,Gr,Czero3,N0)
C---  Expression for rearrangement of extension of Eq. 5.69
C---  Calculates Diii, requires D00iii
C---  Small terms of order Gr(i,j)*Dijklm
C---  Denominator m0sq
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
      integer ep,N0,i1,i2,i3,m,n,np
      parameter(np=3)
      real(dp):: m0sq,Gr(np,np)
      complex(dp):: Czero3(z3max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(diiiii(z5(n,m,i1,i2,i3))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dzziii(z3(i1,i2,i3))+N0,ep-1)
      
      Dv(diii(z3(i1,i2,i3))+N0,ep)=
     . (20d0*Dv(dzziii(z3(i1,i2,i3))+N0,ep)
     . -pole
     . -2d0*Czero3(z3(i1,i2,i3),ep)
     . +bit)/(2d0*m0sq)

      enddo
      
      return
      end
