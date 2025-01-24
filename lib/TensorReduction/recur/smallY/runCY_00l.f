      subroutine runCY_00l(k,l,Xtwiddle,Gtwiddle,Shat3,N0)
      implicit none
C---  Expression for Eq. 5.56a
C---  Calculates C00l
C---  Small terms of order Xtwiddle(0,k)*Cii,Xtwiddle(0,0)*Ciii
C---  Denominator Gtwiddle(k,l)
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f' 
      include 'lib/TensorReduction/Include/pvCv.f' 
      include 'lib/TensorReduction/recur/Include/Carraydef.f' 
      include 'lib/TensorReduction/recur/Include/Carrays.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      real(dp):: Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      complex(dp):: Shat3(np,z2max,-2:0)

      do ep=-2,0
      Cv(czzi(l)+N0,ep)=
     . (Gtwiddle(k,1)*Shat3(1,z2(l,l),ep)
     .  +Gtwiddle(k,2)*Shat3(2,z2(l,l),ep)
     .  +Xtwiddle(0,k)*Cv(cii(z2(l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Cv(ciii(z3(k,l,l))+N0,ep))/(4*Gtwiddle(k,l))
      enddo


      return
      end
  



