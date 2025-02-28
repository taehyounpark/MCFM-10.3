      subroutine runCY_00ll(k,l,Xtwiddle,Gtwiddle,Shat4,N0)
      implicit none
C---  Expression for Eq. 5.58a
C---  Calculates C00ll
C---  Small terms of order Xtwiddle(0,k)*Ciii,Xtwiddle(0,0)*Ciiii
C---  Denominator Gtwiddle(k,l)
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f' 
      include 'lib/TensorReduction/Include/pvCv.f' 
      include 'lib/TensorReduction/recur/Include/Carraydef.f' 
      include 'lib/TensorReduction/recur/Include/Carrays.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      real(dp):: Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      complex(dp):: Shat4(np,z3max,-2:0)

      do ep=-2,0
      Cv(czzii(z2(l,l))+N0,ep)=
     .  (Gtwiddle(k,1)*Shat4(1,z3(l,l,l),ep)
     .  +Gtwiddle(k,2)*Shat4(2,z3(l,l,l),ep)
     .  +Xtwiddle(0,k)*Cv(ciii(z3(l,l,l))+N0,ep)
     .  -Xtwiddle(0,0)*Cv(ciiii(z4(k,l,l,l))+N0,ep))/(6*Gtwiddle(k,l))
      enddo


      return
      end
  



