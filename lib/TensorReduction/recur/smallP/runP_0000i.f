      subroutine runP_0000i(i1,Gr,S0000i,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,i1,m,n,np
      parameter(np=3)
      real(dp):: Gr(np,np)
      complex(dp):: S0000i(np,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Dv(dzziii(z3(n,m,i1))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Dv(dzzzzi(i1)+N0,ep-1)

c--- note: we have simplified the recursion relation by using DD (5.9)
c--- so that S0000i = 2C00i(0)+2*m0sq*D00i      
      Dv(dzzzzi(i1)+N0,ep)=
     . (pole
     . +S0000i(i1,ep)
     . -bit)/16d0
      enddo
      
      return
      end
