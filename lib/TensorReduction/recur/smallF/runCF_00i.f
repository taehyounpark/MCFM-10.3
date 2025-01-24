      subroutine runCF_00i(i1,f,Gr,Shat3,N0)
C---  Expression for rearrangement of Eq. 5.68
C---  Calculates C00i
C---  Small terms of order f(i)*Cij,Gr(i,j)*Cijk
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
      complex(dp):: Shat3(np,z2max,-2:0)
       
      do ep=-2,0
      Cv(czzi(i1)+N0,ep)=
     . (Shat3(i1,z2(i1,i1),ep)
     . -f(i1)*Cv(cii(z2(i1,i1))+N0,ep)
     . -Gr(i1,1)*Cv(ciii(z3(1,i1,i1))+N0,ep) 
     . -Gr(i1,2)*Cv(ciii(z3(2,i1,i1))+N0,ep))/4d0

      enddo
      
      return
      end
