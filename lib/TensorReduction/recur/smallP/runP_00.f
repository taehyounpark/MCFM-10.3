      subroutine runP_00(m0sq,Gr,Czero0,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,m,n,np
      parameter(np=3)
      real(dp):: m0sq,Gr(np,np)
      complex(dp):: Czero0(-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(dii(z2(n,m))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dd00+N0,ep-1)
      
      Dv(dd00+N0,ep)=
     . (pole
     . +2d0*Czero0(ep)
     . +2d0*m0sq*Dv(dd0+N0,ep)
     . -bit)/8d0
      enddo
      
      return
      end
