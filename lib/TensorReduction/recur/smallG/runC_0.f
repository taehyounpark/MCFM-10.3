      subroutine runC_0(j,DetGr,Xtwiddle0,Gtwiddle,Shat1,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      integer ep,N0,j,np
      parameter(np=2)
      real(dp):: DetGr,Xtwiddle0(np),Gtwiddle(np,np)
      complex(dp):: Shat1(np,-2:0)
       
      do ep=-2,0
      Cv(cc0+N0,ep)=
     . -(Gtwiddle(j,1)*Shat1(1,ep)
     .  +Gtwiddle(j,2)*Shat1(2,ep)
     . -DetGr*Cv(ci(j)+N0,ep))/Xtwiddle0(j)
      enddo
      return
      end
