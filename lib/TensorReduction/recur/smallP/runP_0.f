      subroutine runP_0(k,f,Gr,Shat1,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,k,np
      parameter(np=3)
      real(dp):: f(np),Gr(np,np)
      complex(dp):: Shat1(np,-2:0)
       
      do ep=-2,0
      Dv(dd0+N0,ep)=
     . (Shat1(k,ep)
     . -Gr(k,1)*Dv(di(1)+N0,ep) 
     . -Gr(k,2)*Dv(di(2)+N0,ep) 
     . -Gr(k,3)*Dv(di(3)+N0,ep))/f(k) 
      enddo
      
      return
      end
