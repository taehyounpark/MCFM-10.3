      subroutine runCY_i1i2(i,j,i1,i2,f,Xtwiddle,Gtt,Gtwiddle,Shat3,
     . Bzero2,N0)
      implicit none
C---  Expression for Eq. 5.59
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f' 
      include 'lib/TensorReduction/Include/pvCnames.f' 
      include 'lib/TensorReduction/Include/pvCv.f' 
      include 'lib/TensorReduction/recur/Include/Carraydef.f' 
      include 'lib/TensorReduction/recur/Include/Carrays.f' 
      integer ep,N0,i,j,i1,i2,n,m,np
      parameter(np=2)
      real(dp):: Xtwiddle(0:np,0:np),Gtwiddle(np,np),f(np),
     . Gtt(np,np,np,np)
      complex(dp):: Shat3(np,z2max,-2:0),Bzero2(z2max,-2:0),
     . bit,pole


      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit+Gtt(i,n,j,m)*f(n)
     . *(Shat3(m,z2(i1,i2),ep)
     . -2d0*delta(m,i1)*Cv(czzi(i2)+N0,ep)
     . -2d0*delta(m,i2)*Cv(czzi(i1)+N0,ep))
      enddo
      enddo

      pole=czip
      if (ep .gt. -2) pole=-4*Cv(czzii(z2(i1,i2))+N0,ep-1)

      Cv(cii(z2(i1,i2))+N0,ep)=
     .  (Gtwiddle(i,j)
     . *(8d0*Cv(czzii(z2(i1,i2))+N0,ep)+pole-Bzero2(z2(i1,i2),ep))
     . +bit+Xtwiddle(0,j)*Cv(ciii(z3(i,i1,i2))+N0,ep))/Xtwiddle(i,j)

      enddo
      return
      end
  



