      subroutine runCY_0000l(k,l,Xtwiddle,Gtwiddle,Shat5zz,N0)
      implicit none
C---  Expression for C0000 obtained from 5.50, following the comment after 
C     5.50 on how to add adcing adcitional "00" pairs
C---  (similar to Eq. 5.56a but with "00" added) 
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f' 
      include 'lib/TensorReduction/Include/pvCv.f' 
      include 'lib/TensorReduction/recur/Include/Carraydef.f' 
      include 'lib/TensorReduction/recur/Include/Carrays.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      real(dp):: Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      complex(dp):: Shat5zz(np,z2max,-2:0)

      do ep=-2,0

      Cv(czzzzi(l)+N0,ep) = -(
     . -(Gtwiddle(k,1)*Shat5zz(1,z2(l,l),ep)
     .  +Gtwiddle(k,2)*Shat5zz(2,z2(l,l),ep)
     .  +Xtwiddle(0,k)*Cv(czzii(z2(l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Cv(czziii(z3(k,l,l))+N0,ep)))/(4d0*Gtwiddle(k,l))

      enddo

      return
      end
  

