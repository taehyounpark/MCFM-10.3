      subroutine runF_i(i1,m0sq,Gr,Czero1,N0)
C---  Expression for rearrangement of Eq. 5.67
C---  Calculates Di, requires D00i
C---  Small terms of order Gr(i,j)*Dijk
C---  Denominator m0sq
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
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
      
      Dv(di(i1)+N0,ep)=    
     . (12d0*Dv(dzzi(i1)+N0,ep)
     . -pole
     . -2d0*Czero1(i1,ep)
     . +bit)/(2d0*m0sq)

      enddo
      
      return
      end
