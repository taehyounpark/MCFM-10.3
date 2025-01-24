!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine hard_DY(p,order,msq0,msq1,msq2)
      implicit none
      include 'types.f'
c---- LO * (hard function) for Drell-Yan processes, in units of as/4/pi
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'kprocess.f'
      real(dp),intent(out):: msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
      integer:: order
      real(dp):: p(mxpart,4),s12,msq0(-nf:nf,-nf:nf),hard(2),msqextra(-nf:nf,-nf:nf)

      s12=two*(p(1,4)*p(2,4)-p(1,1)*p(2,1)
     &        -p(1,2)*p(2,2)-p(1,3)*p(2,3))

      call hardqq(s12,musq,hard)
c get back to coefficients of (as/4/pi) from (as/2/pi) returned by hardqq
      hard(1)=hard(1)*2._dp
      hard(2)=hard(2)*4._dp

      if     (kcase == kW_only) then
        call qqb_w(p,msq0)
      elseif (kcase == kZ_only) then
        call qqb_z(p,msq0)
      elseif ((kcase==kWHbbar) .or. (kcase==kWHgaga) .or. (kcase==kWH__WW)) then
        call qqb_wh(p,msq0)
      elseif ((kcase==kZHbbar) .or. (kcase==kZHgaga) .or. (kcase==kZH__WW)) then
        call qqb_zh(p,msq0)
      else
        write(6,*) 'Unrecognized case in hard_DY: kcase = ',kcase
        stop

      endif

      msq1(:,:)=msq0(:,:)*hard(1)
      msq2(:,:)=msq0(:,:)*hard(2)

c--- add contributions that do not factorize over LO
      if     ((kcase==kWHbbar) .or. (kcase==kWHgaga) .or. (kcase==kWH__WW)) then
        call qqb_whas2(p,msqextra)
        msq2(:,:)=msq2(:,:)+msqextra(:,:)/ason4pi**2
      elseif ((kcase==kZHbbar) .or. (kcase==kZHgaga) .or. (kcase==kZH__WW)) then
        call qqb_zhas2(p,msqextra)
        msq2(:,:)=msq2(:,:)+msqextra(:,:)/ason4pi**2
      endif

      return
      end


