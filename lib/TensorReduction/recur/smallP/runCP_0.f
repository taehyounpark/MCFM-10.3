      subroutine runCP_0(k,f,Gr,Shat1,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      integer ep,N0,k,np
      parameter(np=2)
      real(dp):: f(np),Gr(np,np)
      complex(dp):: Shat1(np,-2:0)
       
      do ep=-2,0
      Cv(cc0+N0,ep)=
     . (Shat1(k,ep)
     . -Gr(k,1)*Cv(ci(1)+N0,ep) 
     . -Gr(k,2)*Cv(ci(2)+N0,ep))/f(k) 
      enddo
      
      return
      end
