      subroutine runCF_00(i1,f,Gr,Shat2,N0)
C---  Expression for rearrangement of Eq. 5.66
C---  Calculates C00
C---  Small terms of order f(i)*Ci,Gr(i,j)*Cij
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
      integer ep,N0,i1,np
      parameter(np=2)
      real(dp):: f(np),Gr(np,np)
      complex(dp):: Shat2(np,np,-2:0)
       
      do ep=-2,0
      Cv(cc00+N0,ep)=
     . (Shat2(i1,i1,ep)
     . -f(i1)*Cv(ci(i1)+N0,ep)
     . -Gr(i1,1)*Cv(cii(z2(1,i1))+N0,ep) 
     . -Gr(i1,2)*Cv(cii(z2(2,i1))+N0,ep))/2d0

      enddo
      
      return
      end
