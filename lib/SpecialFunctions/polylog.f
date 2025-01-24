      function polylog(n,x)
      implicit none
      include 'types.f'
      real(dp):: x
      real(dp):: polylog
      complex(dp):: Hc1(-1:1),Hc2(-1:1,-1:1),Hc3(-1:1,-1:1,-1:1),Hc4(-1:1,-1:1,-1:1,-1:1)
      real(dp):: Hr1(-1:1),Hr2(-1:1,-1:1),Hr3(-1:1,-1:1,-1:1),Hr4(-1:1,-1:1,-1:1,-1:1)
      real(dp):: Hi1(-1:1),Hi2(-1:1,-1:1),Hi3(-1:1,-1:1,-1:1),Hi4(-1:1,-1:1,-1:1,-1:1)
      integer:: n
      if (n==2) then
         call hplog(x,2,Hc1,Hc2,Hc3,Hc4,
     &                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1)
         polylog=real(Hc2(0,1),kind=dp)
      elseif (n==3) then
         call hplog(x,3,Hc1,Hc2,Hc3,Hc4,
     &                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1)
         polylog=real(Hc3(0,0,1),kind=dp)
      elseif (n==4) then
         call hplog(x,4,Hc1,Hc2,Hc3,Hc4,
     &                  Hr1,Hr2,Hr3,Hr4,Hi1,Hi2,Hi3,Hi4,-1,1)
         polylog=real(Hc4(0,0,0,1),kind=dp)
      else
         write(6,*) 'polylog(n,x) is only defined for n=2,3,4'
         stop
      endif
      return
      end

