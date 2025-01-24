      subroutine runCF_ii(i1,i2,m0sq,Gr,Bzero2,N0)
C---  Expression for rearrangement of Eq. 5.69
C---  Calculates Cii, requires C00ii
C---  Small terms of order Gr(i,j)*Cijkl
C---  Denominator m0sq
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
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
      
      Cv(cii(z2(i1,i2))+N0,ep)=
     . (16d0*Cv(czzii(z2(i1,i2))+N0,ep)
     . -pole
     . -2d0*Bzero2(z2(i1,i2),ep)
     . +bit)/(2d0*m0sq)

      enddo
      
      return
      end
