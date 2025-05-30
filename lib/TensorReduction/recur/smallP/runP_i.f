      subroutine runP_i(k,i1,f,Gr,Shat2,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,k,i1,np
      parameter(np=3)
      real(dp):: f(np),Gr(np,np)
      complex(dp):: Shat2(np,np,-2:0)
       
      do ep=-2,0
      Dv(di(i1)+N0,ep)=
     . (Shat2(k,i1,ep)
     . -2d0*delta(k,i1)*Dv(dd00+N0,ep)
     . -Gr(k,1)*Dv(dii(z2(1,i1))+N0,ep) 
     . -Gr(k,2)*Dv(dii(z2(2,i1))+N0,ep)
     . -Gr(k,3)*Dv(dii(z2(3,i1))+N0,ep))/f(k) 
      enddo
      
      return
      end
