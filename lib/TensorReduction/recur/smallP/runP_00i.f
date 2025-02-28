      subroutine runP_00i(i1,m0sq,Gr,Czero1,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,i1,m,n,np
      parameter(np=3)
      real(dp):: m0sq,Gr(np,np)
      complex(dp):: Czero1(z1max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(diii(z3(n,m,i1))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dzzi(i1)+N0,ep-1)
      
      Dv(dzzi(i1)+N0,ep)=
     . (pole
     . +2d0*Czero1(i1,ep)
     . +2d0*m0sq*Dv(di(i1)+N0,ep)
     . -bit)/12d0
      enddo
      
      return
      end
