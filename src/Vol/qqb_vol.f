!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_vol(P,msq)
      implicit none
      include 'types.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'sprods_com.f'

      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      call dotem(mxpart,p,s)
      msq(:,:)=0._dp

      msq(2,-1)=1._dp
      return
      end

c      subroutine qqb_vol_g(P,msq)
c      implicit none
c      include 'types.f'

c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'masses.f'
c      include 'sprods_com.f'
c      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
c      integer:: j,k,N
c      N=7
c      call dotem(N,p,s)
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(j,k)=0._dp
c      enddo
c      enddo

c      msq(2,-1)=1._dp
c      return
c      end

c      subroutine qqb_vol_gs(P,msq)
c      implicit none
c      include 'types.f'

c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'masses.f'
c      include 'ptilde.f'
c      include 'sprods_com.f'
c      real(dp):: P(mxpart,4),msq(maxd,-nf:nf,-nf:nf)
c      integer:: j,k,nd,N
c      N=7
c      call dotem(N,p,s)
c      do nd=1,maxd
c      do j=-nf,nf
c      do k=-nf,nf
c      msq(nd,j,k)=0._dp
c      enddo
c      enddo
c      enddo
c      ndmax=0

c      return
c      end
