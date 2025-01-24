!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine hard_H(p,order,msq0,msq1,msq2)
      implicit none
      include 'types.f'
c---- LO * (hard function) for Drell-Yan processes, in units of as/4/pi
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'scale.f'
      include 'kprocess.f'
      real(dp),intent(out):: msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
      integer order
      real(dp):: p(mxpart,4),s12,msq0(-nf:nf,-nf:nf),hard(2)

      s12=two*(p(1,4)*p(2,4)-p(1,1)*p(2,1)
     &        -p(1,2)*p(2,2)-p(1,3)*p(2,3))

      call hardgg(s12,musq,hard)
c get back to coefficients of (as/4/pi) from (as/2/pi) returned by hardqq
      hard(1)=hard(1)*2._dp
      hard(2)=hard(2)*4._dp

      if     (kcase == kggfus0) then
        call gg_h(p,msq0)
      else
        write(6,*) 'Unrecognized case in hard_H: kcase = ',kcase
        stop

      endif

      msq1(:,:)=msq0(:,:)*hard(1)
      msq2(:,:)=msq0(:,:)*hard(2)

      return
      end


