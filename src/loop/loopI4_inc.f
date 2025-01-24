!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      real(dp):: mu2
      complex(dp):: loopI4_res,resQCDLoop,resOneLOop,result(0:2)
      integer:: ep
c--- this is a switchyard routine: depending on value of
c--- integer variable, this routine either:
c---  scalarselect = 1    calls QCDLoop for scalar integral
c---  scalarselect = 2    calls OneLOop for scalar integral
c---  scalarselect = 3    calls both routines and compares results

c--- call to QCDLoop if necessary      
      if ((scalarselect .eq. 1) .or.(scalarselect .eq. 3)) then
        resQCDLoop=qlI4(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
      endif
      
      if (scalarselect .eq. 1) then
        loopI4_res=resQCDLoop
        return
      endif

      call olo_d0(result,p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,sqrt(mu2))
      resOneLOop=result(abs(ep))
      loopI4_res=resOneLOop
      
      if (scalarselect .eq. 3) then
        if ((cdabs(resOneLOop) .gt. 1d-9) .and.
     &      (abs(resQCDLoop/resOneLOop-1d0) .gt. 1d-9)) then
          write(6,*) 'loopI4: ',p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep
          write(6,*) 'QCDLoop:',resQCDLoop
          write(6,*) 'OneLOop:',resOneLOop
          write(6,*) '->ratio:',resQCDLoop/resOneLOop
        endif
      endif
      return
      
