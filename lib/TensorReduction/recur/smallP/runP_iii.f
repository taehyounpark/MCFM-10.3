      subroutine runP_iii(k,i1,i2,i3,f,Gr,Shat4,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,k,i1,i2,i3,np
      parameter(np=3)
      real(dp):: f(np),Gr(np,np)
      complex(dp):: Shat4(np,z3max,-2:0)
       
      do ep=-2,0
      Dv(diii(z3(i1,i2,i3))+N0,ep)=
     . (Shat4(k,z3(i1,i2,i3),ep)
     . -2d0*delta(k,i1)*Dv(dzzii(z2(i2,i3))+N0,ep)
     . -2d0*delta(k,i2)*Dv(dzzii(z2(i1,i3))+N0,ep)
     . -2d0*delta(k,i3)*Dv(dzzii(z2(i1,i2))+N0,ep)
     . -Gr(k,1)*Dv(diiii(z4(1,i1,i2,i3))+N0,ep) 
     . -Gr(k,2)*Dv(diiii(z4(2,i1,i2,i3))+N0,ep)
     . -Gr(k,3)*Dv(diiii(z4(3,i1,i2,i3))+N0,ep))/f(k) 
      enddo
      
      return
      end
