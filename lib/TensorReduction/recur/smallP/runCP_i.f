      subroutine runCP_i(k,i1,f,Gr,Shat2,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      integer ep,N0,k,i1,np
      parameter(np=2)
      real(dp):: f(np),Gr(np,np)
      complex(dp):: Shat2(np,np,-2:0)
       
      do ep=-2,0
      Cv(ci(i1)+N0,ep)=
     . (Shat2(k,i1,ep)
     . -2d0*delta(k,i1)*Cv(cc00+N0,ep)
     . -Gr(k,1)*Cv(cii(z2(1,i1))+N0,ep) 
     . -Gr(k,2)*Cv(cii(z2(2,i1))+N0,ep))/f(k) 
      enddo
      
      return
      end
