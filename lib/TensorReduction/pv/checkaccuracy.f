      subroutine checkaccuracy(trhs,tq,prec,failed) 
      implicit none
      include 'lib/TensorReduction/Include/types.f'
c      include 'lib/TensorReduction/Include/pvverbose.f'
      real(dp)::prec
      complex(dp):: trhs,tq 
      logical:: failed

c--- equation we are testing is tq+trhs=0

      if (abs(trhs-tq) < 1d-15) return

c--- if sum is small relative to difference, equation is satisfied
      if (abs((trhs+tq)/(trhs-tq)) .lt. prec) then
        return
      endif

c--- if both are small, no problem
      if ((abs(trhs) .lt. 1d-6) .and. (abs(tq) .lt. 1d-6)) then
        return
      endif

c--- otherwise, equation is violated
      failed=.true.

      return

c          write(6,*) 'checkacc: ',ep,trhs/tq,tq,prec
c      if ( (abs(tq) .gt. prec .and. abs(trhs/tq) .gt. prec) .or. 
c     .     (abs(tq) .lt. prec .and. abs(trhs) .gt. prec)) then 
c          if (pvverbose) write(6,*) 'checkacc: ',ep,trhs/tq,tq,prec
c          failed=.true.
c          pause 
c      endif

      end

      
