!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine softcheck(p)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      integer j,k
      real(dp):: p(mxpart,4),msqg(-nf:nf,-nf:nf),
     &     msqs(-nf:nf,-nf:nf)

      call writeout(p)

      write(6,*) 'j, k, msqg, msqs, g/s: '

      call qqb_QQb_mix_g(p,msqg)
      call qqb_QQb_mix_sft(p,msqs)


      do j=-nf,nf
         k = -j
         if (j  /=  0) then
            write(6,*) j, k,
     &           msqg(j,k),msqs(j,k),msqg(j,k)/msqs(j,k)
         end if
      end do
      write(6,*) 'pause in softcheck.f :push enter to continue'
      read(5,*)

      end subroutine softcheck
