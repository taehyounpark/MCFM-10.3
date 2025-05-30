      subroutine run_0000(k,l,DetGr,f,Gtwiddle,Gtt,
     . Shati00,Si00m,S0000,N0)
      implicit none
C------Fixes D0000 using 5.46
C-----knowing D00 with correction of order Delta D00ii
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRconstants.f'
      include 'lib/TensorReduction/Include/pvDnames.f'
      include 'lib/TensorReduction/Include/pvDv.f'
      include 'lib/TensorReduction/recur/Include/Darraydef.f'
      include 'lib/TensorReduction/recur/Include/Darrays.f'
      integer ep,N0,k,l,n,m
      real(dp):: DetGr,Gtwiddle(3,3),Gtt(3,3,3,3),f(3)
      complex(dp):: Shati00(3,-2:0),S0000(-2:0),Si00m(3,3,-2:0),bit,pole

       
      do ep=-2,0
      bit=czip
      do n=1,3
      do m=1,3
      bit=bit
     . +Gtt(k,n,l,m)*(f(n)*Shati00(m,ep)-f(n)*f(m)*Dv(dd00+N0,ep))
      enddo
      enddo
      pole=czip
      if (ep .gt. -2) pole=-4d0*Gtwiddle(k,l)*Dv(dd0000+N0,ep-1)

      Dv(dd0000+N0,ep)=-(
     . +pole
     . +DetGr*Dv(dzzii(z2(k,l))+N0,ep)
     . -Gtwiddle(k,l)*S0000(ep)
     . -Gtwiddle(1,l)*Si00m(1,k,ep)
     . -Gtwiddle(2,l)*Si00m(2,k,ep)
     . -Gtwiddle(3,l)*Si00m(3,k,ep)
     . +Gtwiddle(k,l)*(Si00m(1,1,ep)+Si00m(2,2,ep)+Si00m(3,3,ep))
     . +bit)/(8d0*Gtwiddle(k,l))
      enddo


      return
      end
  



