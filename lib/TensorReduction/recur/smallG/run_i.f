      subroutine run_i(j,i1,DetGr,Xtwiddle0,Gtwiddle,Shat2,N0)
      implicit none
C---Fixes Di according to 5.43 Denner-Dittmaier
C---knowing D00 with corrections of order Delta_3 Dii
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,j,i1,n
      real(dp):: DetGr,Xtwiddle0(3),Gtwiddle(3,3)
      complex(dp):: Shat2(3,3,-2:0),Shat2s(3,3,-2:0)
       

      do ep=-2,0
      do n=1,3
         Shat2s(n,i1,ep)=Shat2(n,i1,ep)-2d0*delta(n,i1)*Dv(N0+dd00,ep)
      enddo 
      Dv(di(i1)+N0,ep)=-(
     . +Gtwiddle(j,1)*Shat2s(1,i1,ep)
     . +Gtwiddle(j,2)*Shat2s(2,i1,ep)
     . +Gtwiddle(j,3)*Shat2s(3,i1,ep)
     . -DetGr*Dv(dii(z2(j,i1))+N0,ep))/Xtwiddle0(j)

      enddo
      return
      end
