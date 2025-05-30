      subroutine runCF_iiiii(i1,i2,i3,i4,i5,m0sq,Gr,Bzero5,N0)
C---  Expression for rearrangement of extension of Eq. 5.69
C---  Calculates Ciiii, requires C00iiiii
C---  Small terms of order Gr(i,j)*Cijklmno
C---  Denominator m0sq
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvCnames.f'
      include 'lib/TensorReduction/Include/pvCv.f'
      include 'lib/TensorReduction/recur/Include/Carraydef.f'
      include 'lib/TensorReduction/recur/Include/Carrays.f'
      include 'lib/TensorReduction/Include/pvweenumber.f' 
      integer ep,N0,i1,i2,i3,i4,i5,m,n,np
      parameter(np=2)
      real(dp):: m0sq,Gr(np,np)
      complex(dp):: Bzero5(z5max,-2:0),bit,pole
       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gr(n,m)*Cv(ciiiiiii(z7(n,m,i1,i2,i3,i4,i5))+N0,ep)
      enddo
      enddo
      
      pole=czip
      if (ep .gt. -2) pole=4d0*Cv(czziiiii(z5(i1,i2,i3,i4,i5))+N0,ep-1)
      
      Cv(ciiiii(z5(i1,i2,i3,i4,i5))+N0,ep)=
     . (28d0*Cv(czziiiii(z5(i1,i2,i3,i4,i5))+N0,ep)
     . -pole
     . -2d0*Bzero5(z5(i1,i2,i3,i4,i5),ep)
     . +bit)/(2d0*m0sq)

      enddo
      
      return
      end
