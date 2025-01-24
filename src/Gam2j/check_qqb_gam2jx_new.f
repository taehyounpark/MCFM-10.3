!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine check_qqb_gam2jx_new(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'pp.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: mqq(0:2,fn:nf,fn:nf),ppmsqx(0:2,ppmax),
     & msqx_cs(0:2,-nf:nf,-nf:nf)
      logical:: swap56
      common/swap56/swap56
!$omp threadprivate(/swap56/)

      call qqb_gam2jx_new(p,msq,mqq,ppmsqx,msqx_cs)

c swap56 is a variable which, if true, returns the matrix elements
c with 5 and 6 in the 'unnatural' ordering not usually provided by msqx_cs

c--- Only pick out Gflag pieces
      do j=-4,4
      do k=-4,4

      if     ((j==0) .and. (k==0)) then
        if (swap56) then
          msqx_cs(:,j,k)=3._dp*ppmsqx(:,pp(j,k,1,-1))
     &                  +2._dp*ppmsqx(:,pp(j,k,2,-2))
        else
          ! msqx_cs already filled
        endif
      elseif ((j>0) .and. (k<0)) then
        ! msqx_cs already filled and only one ordering is necessary
      elseif ((j<0) .and. (k>0)) then
        ! msqx_cs already filled and only one ordering is necessary
      elseif ((j>0) .and. (k==0)) then
        if (swap56) then
          msqx_cs(:,j,k)=ppmsqx(:,pp(j,k,k,j))
        else
          ! msqx_cs already filled
        endif
      elseif ((j<0) .and. (k==0)) then
        if (swap56) then
          msqx_cs(:,j,k)=ppmsqx(:,pp(j,k,k,j))
        else
          ! msqx_cs already filled
        endif
      elseif ((j==0) .and. (k>0)) then
        if (swap56) then
          msqx_cs(:,j,k)=ppmsqx(:,pp(j,k,j,k))
        else
          ! msqx_cs already filled
        endif
      elseif ((j==0) .and. (k<0)) then
        if (swap56) then
          msqx_cs(:,j,k)=ppmsqx(:,pp(j,k,j,k))
        else
          ! msqx_cs already filled
        endif
      endif

      enddo
      enddo

      if (swap56) then
c Only the sum of color structures is returned in gvecx routine
      else
c DEBUG: only keep one color ordering
        msqx_cs(0,:,:)=zip
        msqx_cs(2,:,:)=zip
      endif
c return msq array in usual format by performing sum
      msq(:,:)=msqx_cs(0,:,:)+msqx_cs(1,:,:)+msqx_cs(2,:,:)

      return
      end

