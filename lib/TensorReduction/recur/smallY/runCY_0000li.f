      subroutine runCY_0000li(k,l,i1,Xtwiddle,Gtwiddle,Shat6zz,N0)
      implicit none
C---  Expression for Eq. 5.58b
C---  Calculates C0000li, requires C0000ll
C---  Small terms of order Xtwiddle(0,k)*C00iii,Xtwiddle(0,0)*C00iiii
C---  Denominator Gtwiddle(k,l)
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f' 
      include 'lib/TensorReduction/Include/pvCv.f' 
      include 'lib/TensorReduction/recur/Include/Carraydef.f' 
      include 'lib/TensorReduction/recur/Include/Carrays.f' 
      integer ep,N0,k,l,i1,np
      parameter(np=2)
      real(dp):: Xtwiddle(0:np,0:np),Gtwiddle(np,np)
      complex(dp):: Shat6zz(np,z3max,-2:0)

      if ((i1 .eq. l) .or. (i1 .eq. 0)) then
      return
      endif

      do ep=-2,0

      Cv(czzzzii(z2(l,i1))+N0,ep)=
     . (-2d0*Gtwiddle(k,i1)*Cv(czzzzii(z2(l,l))+N0,ep)
     . +Gtwiddle(k,1)*Shat6zz(1,z3(l,l,i1),ep)
     . +Gtwiddle(k,2)*Shat6zz(2,z3(l,l,i1),ep)
     . +Xtwiddle(0,k)*Cv(czziii(z3(l,l,i1))+N0,ep)
     . -Xtwiddle(0,0)*Cv(czziiii(z4(k,l,l,i1))+N0,ep)
     . )/(4d0*Gtwiddle(k,l))
 
      enddo


      return
      end
  



