      subroutine run_000000(k,l,DetGr,f,Gtwiddle,Gtt,
     . Shat5zzzz,Shat6zzzz,S000000,N0)
      implicit none
C------Fixes D000000 
C-----knowing D0000 with correction of order Delta D0000ii
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'  
      include 'lib/TensorReduction/Include/pvDnames.f'  
      include 'lib/TensorReduction/Include/pvDv.f'  
      include 'lib/TensorReduction/recur/Include/Darraydef.f'  
      include 'lib/TensorReduction/recur/Include/Darrays.f'  
      integer ep,N0,k,l,n,m,np
      parameter(np=3)
      real(dp):: DetGr,Gtwiddle(np,np),Gtt(np,np,np,np),f(np)
      complex(dp):: Shat5zzzz(np,-2:0),S000000(-2:0),
     . Shat6zzzz(np,z1max,-2:0),bit,pole

       
      do ep=-2,0
      bit=czip
      do n=1,np
      do m=1,np
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shat5zzzz(m,ep)-f(n)*f(m)*Dv(dd0000+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2) pole=-4d0*Gtwiddle(k,l)*Dv(dd000000+N0,ep-1)

      Dv(dd000000+N0,ep)=
     . -(pole
     . +DetGr*Dv(dzzzzii(z2(k,l))+N0,ep)
     . -Gtwiddle(k,l)*S000000(ep)
     . -Gtwiddle(1,l)*Shat6zzzz(1,k,ep)
     . -Gtwiddle(2,l)*Shat6zzzz(2,k,ep)
     . -Gtwiddle(3,l)*Shat6zzzz(3,k,ep)
     . +Gtwiddle(k,l)
     . *(Shat6zzzz(1,1,ep)+Shat6zzzz(2,2,ep)+Shat6zzzz(3,3,ep))
     . +bit)/(12d0*Gtwiddle(k,l))
      enddo


      return
      end
  



