!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine hard_VV(p,order,msq0,msq1,msq2)
      implicit none
      include 'types.f'
c---- LO * (hard function) for diboson processes, in units of as/4/pi
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'kprocess.f'
      include 'qcdcouple.f'
      include 'ipsgen.f'
      real(dp),intent(out):: msq0(-nf:nf,-nf:nf),msq1(-nf:nf,-nf:nf),msq2(-nf:nf,-nf:nf)
      integer order
      real(dp):: p(mxpart,4),msqfull(-nf:nf,-nf:nf,0:2),msqgg(-nf:nf,-nf:nf),pttwo

      msqfull(:,:,:)=zip

! back-to-back check
      if (pttwo(3,4,p) < 1.e-3_dp) then
        msq0=zip
        msq1=zip
        msq2=zip
        return
      endif

c At NNLO this routine is also called with ipsgen = 2 to fill the gg contribution
      if ((order == 2) .and. (ipsgen == 1)) then
c Numerical safety cut for 2-loop amplitudes
        if (pttwo(3,4,p) > 2.5_dp) then
          call vv_amps(p,msqfull,order)
        endif

      endif

c Use existing BDK-style MCFM results up to NLO
      epinv=0._dp
      epinv2=0._dp

      msq0=zip
      msq1=zip
      if (kcase == kWWqqbr) then
        if ((order == 2) .and. (ipsgen == 2)) then
          call gg_WW_int(p,msqgg) ! additional gluon-gluon contribution
          msqfull(0,0,2)=msqgg(0,0)
        else
          call qqb_ww(p,msq0)
          if (order > 0) call qqb_ww_v(p,msq1)
        endif
      elseif (kcase == kWZbbar) then
        call qqb_wz(p,msq0)
        if (order > 0) call qqb_wz_v(p,msq1)
      elseif (kcase == kZZlept) then
        if ((order == 2) .and. (ipsgen == 2)) then
          call gg_ZZ_all(p,msqgg) ! additional gluon-gluon contribution
          msqfull(0,0,2)=msqgg(0,0)
        else
          call qqb_zz(p,msq0)
          if (order > 0) call qqb_zz_v(p,msq1)
        endif
      else
        write(6,*) 'Unexpected case in hard_VV: kcase = ',kcase
        stop
      endif
      msqfull(:,:,0)=msq0(:,:)
      msqfull(:,:,1)=msq1(:,:)-ason2pi*CF*(1d0-pisq/6._dp)*msq0(:,:)

c change to contain just coefficients of (as/4/pi) and (as/4/pi^2)
      msq1(:,:)=msqfull(:,:,1)/ason4pi
      msq2(:,:)=msqfull(:,:,2)/ason4pi**2

      return
      end

c---          printed=.false.
c---          do j=-nf,nf
c---          do k=-nf,nf
c---!          if (abs(msq(j,k,1)/(msqv(j,k)-ason2pi*CF*(1d0-pisq/6._dp)*msq0(j,k))-1._dp) > 1.e-4_dp) then
c---          if (abs(msq(j,k,2)/msq(j,k,0)) > 1.e3_dp) then
c---!          if ((msq(j,k,0) /= 0._dp) .or. (msq(j,k,1) /= 0._dp) .or. (msq(j,k,2) /= 0._dp)) then
c---            write(6,*) j,k,msq0(j,k),msq(j,k,:)
c---            printed=.true.
c---          endif
c---          enddo
c---          enddo
c---!          pause

c          do j=1,100
c!          call qqb_ww_v(p,msqv)
c          call vv_amps(p,msq,1)
c          enddo
c          stop

c---! Check agreement
c---          call vv_amps(p,msq,order)
c---          write(6,*) 'j,k,rat,GvMT,MCFM'
c---          do j=-nf,nf
c---          do k=-nf,nf
c---          if (msq(j,k,0) /= 0._dp) then
c---!            write(6,*) 'MCFM',j,k,msq0(j,k)
c---!            write(6,*) 'GvMT',j,k,msq(j,k,0)
c---            write(6,*) '0L',j,k,msq(j,k,0)/msq0(j,k),msq(j,k,0),msq0(j,k)
c---          endif
c---          enddo
c---          enddo
c---!! translation to hard function, c.f. src/Diboson2/source/test1.f
c---          msqv(:,:)=msqv(:,:)-ason2pi*CF*(1d0-pisq/6._dp)*msq0(:,:)
c---          do j=-nf,nf
c---          do k=-nf,nf
c---          if ((msq(j,k,1) /= 0._dp) .or. (msqv(j,k) /= 0._dp)) then
c---            write(6,*) '1L',j,k,msq(j,k,1)/msqv(j,k),msq(j,k,1),msqv(j,k)
c---          endif
c---          enddo
c---          enddo
c---          pause

