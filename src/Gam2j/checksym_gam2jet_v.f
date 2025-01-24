!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine checksym_gam2jet_v(p)
c Routine for checking symmetry between beams in the photon+2jet virtual processes
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'flags.f'
      include 'nf.f'
      include 'ptilde.f'
      real(dp):: p(mxpart,4),p12(mxpart,4),p1245(mxpart,4)
      real(dp):: msq(-nf:nf,-nf:nf),msq12(-nf:nf,-nf:nf),msq1245(-nf:nf,-nf:nf)
      real(dp):: msqv(-nf:nf,-nf:nf),msqv12(-nf:nf,-nf:nf),msqv1245(-nf:nf,-nf:nf)
      integer j,k

c Momentum array with p1 and p2 interchanged
      p12(:,:)=p(:,:)
      p12(1,:)=p(2,:)
      p12(2,:)=p(1,:)

c Momentum array with p1 and p2, p4 and p5 interchanged
      p1245(:,:)=p12(:,:)
      p1245(4,:)=p12(5,:)
      p1245(5,:)=p12(4,:)

c First do Gflag checks, which require only p12 throughout
      Gflag=.true.
      Qflag=.false.

      call qqb_gam2j(p,msq)
      call qqb_gam2j(p12,msq12)

      call qqb_gam2j_v(p,msqv)
      call qqb_gam2j_v(p12,msqv12)

      write(6,*) '******************** Gflag ********************'
      do j=-nf,nf
      do k=-nf,nf
        if (abs(msq(j,k)) > 1.e-25_dp) then
          write(6,*) 'check LO  ',j,k,msq(j,k),msq12(k,j),msq(j,k)/msq12(k,j)
        endif
      enddo
      enddo
      write(6,*)

      do j=-nf,nf
      do k=-nf,nf
        if (abs(msq(j,k)) > 1.e-25_dp) then
          write(6,*) 'check virt',j,k,msqv(j,k),msqv12(k,j),msqv(j,k)/msqv12(k,j)
        endif
      enddo
      enddo
      write(6,*)

c Now do Qflag checks, which require both p12 and p1245 depending on parton flavors
      Gflag=.false.
      Qflag=.true.

      call qqb_gam2j(p,msq)
      call qqb_gam2j(p12,msq12)
      call qqb_gam2j(p1245,msq1245)

      call qqb_gam2j_v(p,msqv)
      call qqb_gam2j_v(p12,msqv12)
      call qqb_gam2j_v(p1245,msqv1245)

      write(6,*) '******************** Qflag ********************'
      do j=-nf,nf
      do k=-nf,nf
        if (abs(msq(j,k)) > 1.e-25_dp) then
          if (j*k <= 0) then
            write(6,*) 'check LO  ',j,k,msq(j,k),msq12(k,j),msq(j,k)/msq12(k,j)
          else
            write(6,*) 'check LO  ',j,k,msq(j,k),msq1245(k,j),msq(j,k)/msq1245(k,j)
          endif
        endif
      enddo
      enddo
      write(6,*)

      do j=-nf,nf
      do k=-nf,nf
        if (abs(msq(j,k)) > 1.e-25_dp) then
          if (j*k <= 0) then
            write(6,*) 'check virt',j,k,msqv(j,k),msqv12(k,j),msqv(j,k)/msqv12(k,j)
          else
            write(6,*) 'check virt',j,k,msqv(j,k),msqv1245(k,j),msqv(j,k)/msqv1245(k,j)
          endif
        endif
      enddo
      enddo
      write(6,*)

      stop

      return
      end

