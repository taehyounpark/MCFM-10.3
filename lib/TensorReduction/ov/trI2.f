      function trI2(p1,m1,m2,mu2,ep) 
c--- this is a switchyard routine: depending on value of
c--- integer variable, this routine either:
c---  TRscalarselect = 1    calls QCDLoop for scalar integral
c---  TRscalarselect = 2    calls OneLOop for scalar integral
c---  TRscalarselect = 3    calls both routines and compares results
c      use avh_olo
      use mod_qcdloop_c
      implicit none
      include 'lib/TensorReduction/Include/types.f'
      include 'lib/TensorReduction/Include/TRscalarselect.f'
      real(dp):: p1,m1,m2,mu2
      complex(dp):: trI2,resQCDLoop,resOneLOop,result(0:2)
      integer ep

      TRscalarselect=1

c--- call to QCDLoop if necessary      
      if ((TRscalarselect .eq. 1) .or.(TRscalarselect .eq. 3)) then
        resQCDLoop=qlI2(p1,m1,m2,mu2,ep)
      endif
      
      if (TRscalarselect .eq. 1) then
        trI2=resQCDLoop
        return
      endif
      
      call olo_b0(result,p1,m1,m2,dsqrt(mu2))
      resOneLOop=result(abs(ep))
      trI2=resOneLOop
      
      if (TRscalarselect .eq. 3) then
        if ((cdabs(resOneLOop) .gt. 1d-12) .and.
     &      (abs(resQCDLoop/resOneLOop-1d0) .gt. 1d-12)) then
          write(6,*) 'trI2: ',p1,m1,m2,mu2,ep
          write(6,*) 'QCDLoop:',resQCDLoop
          write(6,*) 'OneLOop:',resOneLOop
          write(6,*) '->ratio:',resQCDLoop/resOneLOop
        endif
      endif
      
      return
      end
      
