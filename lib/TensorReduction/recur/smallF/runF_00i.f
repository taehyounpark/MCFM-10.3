      subroutine runF_00i(i1,f,Gr,Shat3,N0)
C---  Expression for rearrangement of Eq. 5.68
C---  Calculates D00i
C---  Small terms of order f(i)*Dij,Gr(i,j)*Dijk
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
      integer ep,N0,i1,np
      parameter(np=3)
      real(dp):: f(np),Gr(np,np)
      complex(dp):: Shat3(np,z2max,-2:0)
       
      do ep=-2,0
      Dv(dzzi(i1)+N0,ep)=
     . (Shat3(i1,z2(i1,i1),ep)
     . -f(i1)*Dv(dii(z2(i1,i1))+N0,ep)
     . -Gr(i1,1)*Dv(diii(z3(1,i1,i1))+N0,ep) 
     . -Gr(i1,2)*Dv(diii(z3(2,i1,i1))+N0,ep)
     . -Gr(i1,3)*Dv(diii(z3(3,i1,i1))+N0,ep))/4d0

      enddo
      
      return
      end
