      function trI1(m1,mu2,ep) 
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
      real(dp):: m1,mu2
      complex(dp)::trI1,resQCDLoop,resOneLOop,result(0:2)
      integer ep

      TRscalarselect=1

c--- call to QCDLoop if necessary      
      if ((TRscalarselect .eq. 1) .or.(TRscalarselect .eq. 3)) then
        resQCDLoop=qlI1(m1,mu2,ep)
      endif
      
      if (TRscalarselect .eq. 1) then
        trI1=resQCDLoop
        return
      endif

c--- initializations should be moved elsewhere and called only once
      call olo_onshell(1d-8)
      call olo_unit(-1,'error')

      call olo_a0(result,m1,dsqrt(mu2))
      resOneLOop=result(abs(ep))
      trI1=resOneLOop
      
      if (TRscalarselect .eq. 3) then
        if ((cdabs(resOneLOop) .gt. 1d-12) .and.
     &      (abs(resQCDLoop/resOneLOop-1d0) .gt. 1d-12)) then
          write(6,*) 'trI1: ',m1,ep
          write(6,*) 'QCDLoop:',resQCDLoop
          write(6,*) 'OneLOop:',resOneLOop
          write(6,*) '->ratio:',resQCDLoop/resOneLOop
        endif
      endif
      
      return
      end
      
