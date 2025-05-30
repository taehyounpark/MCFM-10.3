      subroutine runCY_00i(k,l,i1,Xtwiddle,Gtwiddle,Shat3,N0)
      implicit none
C---  Expression for Eq. 5.56b
C---  Calculates C00i1
C---  Requires C00l
C---  Small terms of order Xtwiddle(0,k)*Dli1,Xtwiddle(0,0)*Dkli1
C---  Denominator Gtwiddle(k,l)
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f' 
      include 'lib/TensorReduction/Include/pvCv.f' 
      include 'lib/TensorReduction/recur/Include/Carraydef.f' 
      include 'lib/TensorReduction/recur/Include/Carrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=2)
      real(dp):: Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      complex(dp):: Shat3(np,z2max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Cv(czzi(i1)+N0,ep)=
     . (-2*Gtwiddle(k,i1)*Cv(czzi(l)+N0,ep)
     . +Gtwiddle(k,1)*Shat3(1,z2(l,i1),ep)
     . +Gtwiddle(k,2)*Shat3(2,z2(l,i1),ep)
     . +Xtwiddle(0,k)*Cv(cii(z2(l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Cv(ciii(z3(k,l,i1))+N0,ep))/(2*Gtwiddle(k,l))

      enddo

      return
      end
  



