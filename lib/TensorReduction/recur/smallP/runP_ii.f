      subroutine runP_ii(k,i1,i2,f,Gr,Shat3,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,k,i1,i2,np
      parameter(np=3)
      real(dp):: f(np),Gr(np,np)
      complex(dp):: Shat3(np,z2max,-2:0)
       
      do ep=-2,0
      Dv(dii(z2(i1,i2))+N0,ep)=
     . (Shat3(k,z2(i1,i2),ep)
     . -2d0*delta(k,i1)*Dv(dzzi(i2)+N0,ep)
     . -2d0*delta(k,i2)*Dv(dzzi(i1)+N0,ep)
     . -Gr(k,1)*Dv(diii(z3(1,i1,i2))+N0,ep) 
     . -Gr(k,2)*Dv(diii(z3(2,i1,i2))+N0,ep) 
     . -Gr(k,3)*Dv(diii(z3(3,i1,i2))+N0,ep))/f(k) 
      enddo
      
      return
      end
