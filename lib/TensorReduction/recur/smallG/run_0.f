      subroutine run_0(j,DetGr,Xtwiddle0,Gtwiddle,Shat1,N0)
      implicit none
C---Fixes D0 according to 5.41 Denner-Dittmaier
C---with corrections of order Delta Di
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,j
      real(dp):: DetGr,Xtwiddle0(3),Gtwiddle(3,3)
      complex(dp):: Shat1(3,-2:0)
       
      do ep=-2,0
      Dv(dd0+N0,ep)=-(
     . +Gtwiddle(j,1)*Shat1(1,ep)
     . +Gtwiddle(j,2)*Shat1(2,ep)
     . +Gtwiddle(j,3)*Shat1(3,ep)
     . -DetGr*Dv(di(j)+N0,ep))/Xtwiddle0(j)
      enddo

      return
      end
