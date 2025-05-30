      subroutine runCY_00lll(k,l,Xtwiddle,Gtwiddle,Shat5,N0)
      implicit none
C---  Expression for Eq. 5.60a
C---  Calculates C00lll
C---  Small terms of order Xtwiddle(0,k)*Ciiii,Xtwiddle(0,0)*Ciiiii
C---  Denominator Gtwiddle(k,l)
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f' 
      include 'lib/TensorReduction/Include/pvCv.f' 
      include 'lib/TensorReduction/recur/Include/Carraydef.f' 
      include 'lib/TensorReduction/recur/Include/Carrays.f' 
      integer ep,N0,k,l,np
      parameter(np=2)
      real(dp):: Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      complex(dp):: Shat5(np,z4max,-2:0)

      do ep=-2,0

      Cv(czziii(z3(l,l,l))+N0,ep)=
     . (Gtwiddle(k,1)*Shat5(1,z4(l,l,l,l),ep)
     . +Gtwiddle(k,2)*Shat5(2,z4(l,l,l,l),ep)
     . +Xtwiddle(0,k)*Cv(ciiii(z4(l,l,l,l))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciiiii(z5(k,l,l,l,l))+N0,ep))/(8*Gtwiddle(k,l))
      enddo


      return
      end
  



