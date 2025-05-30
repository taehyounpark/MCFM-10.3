      subroutine runCP_iiii(k,i1,i2,i3,i4,f,Gr,Shat5,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      integer ep,N0,k,i1,i2,i3,i4,np
      parameter(np=2)
      real(dp):: f(np),Gr(np,np)
      complex(dp):: Shat5(np,z4max,-2:0)
       
      do ep=-2,0
      Cv(ciiii(z4(i1,i2,i3,i4))+N0,ep)=
     . (Shat5(k,z4(i1,i2,i3,i4),ep)
     . -2d0*delta(k,i1)*Cv(czziii(z3(i2,i3,i4))+N0,ep)
     . -2d0*delta(k,i2)*Cv(czziii(z3(i1,i3,i4))+N0,ep)
     . -2d0*delta(k,i3)*Cv(czziii(z3(i1,i2,i4))+N0,ep)
     . -2d0*delta(k,i4)*Cv(czziii(z3(i1,i2,i3))+N0,ep)
     . -Gr(k,1)*Cv(ciiiii(z5(1,i1,i2,i3,i4))+N0,ep) 
     . -Gr(k,2)*Cv(ciiiii(z5(2,i1,i2,i3,i4))+N0,ep))/f(k) 
      enddo
      
      return
      end
