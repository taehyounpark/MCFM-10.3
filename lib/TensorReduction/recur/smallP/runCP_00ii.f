      subroutine runCP_00ii(i1,i2,m0sq,Gr,Bzero2,N0)
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      integer ep,N0,i1,i2,m,n,np
      parameter(np=2)
      real(dp):: m0sq,Gr(np,np)
      complex(dp):: Bzero2(z2max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciiii(z4(n,m,i1,i2))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czzii(z2(i1,i2))+N0,ep-1)
      
      Cv(czzii(z2(i1,i2))+N0,ep)=
     . (pole
     . +2d0*Bzero2(z2(i1,i2),ep)
     . +2d0*m0sq*Cv(cii(z2(i1,i2))+N0,ep)
     . -bit)/16d0
      enddo
      
      return
      end
