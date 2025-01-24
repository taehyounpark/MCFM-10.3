!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine check_qqb_gam2j_gvecx_new(p,n,in,msqv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'pp.f'
      integer:: in,j,k
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),n(4)
      real(dp):: msqv_cs(0:2,-nf:nf,-nf:nf),ppmsqvx(ppmax)
      logical:: swap56
      common/swap56/swap56
!$omp threadprivate(/swap56/)

      call qqb_gam2j_gvecx_new(p,n,in,msqv,msqv_cs,ppmsqvx)

c swap56 is a variable which, if true, returns the matrix elements
c with 5 and 6 in the 'unnatural' ordering not usually provided by msqv_cs

c The detailed flavor matrix ppmsqvx is only returned as the sum of all
c color structures in gvecx routine, so if we're doing swap56 then simply
c put result in color structure 1, with 0 and 2 being zeroed later anyway

c--- Only pick out Gflag pieces
      do j=-4,4
      do k=-4,4

      if     ((j==0) .and. (k==0)) then
        if (swap56) then
          msqv_cs(1,j,k)=3._dp*ppmsqvx(pp(j,k,1,-1))
     &                  +2._dp*ppmsqvx(pp(j,k,2,-2))
        else
          ! msqv_cs already filled
        endif
      elseif ((j>0) .and. (k<0)) then
        if (swap56) then
          msqv_cs(1,j,k)=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
        else
          ! msqv_cs already filled and only one ordering is necessary
        endif
      elseif ((j<0) .and. (k>0)) then
        if (swap56) then
          msqv_cs(1,j,k)=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
        else
          ! msqv_cs already filled and only one ordering is necessary
        endif
      elseif ((j>0) .and. (k==0)) then
        if (swap56) then
          msqv_cs(1,j,k)=ppmsqvx(pp(j,k,k,j))
        else
          ! msqv_cs already filled
        endif
      elseif ((j<0) .and. (k==0)) then
        if (swap56) then
          msqv_cs(1,j,k)=ppmsqvx(pp(j,k,k,j))
        else
          ! msqv_cs already filled
        endif
      elseif ((j==0) .and. (k>0)) then
        if (swap56) then
          msqv_cs(1,j,k)=ppmsqvx(pp(j,k,j,k))
        else
          ! msqv_cs already filled
        endif
      elseif ((j==0) .and. (k<0)) then
        if (swap56) then
          msqv_cs(1,j,k)=ppmsqvx(pp(j,k,j,k))
        else
          ! msqv_cs already filled
        endif
      endif

      enddo
      enddo

c DEBUG: only keep one color ordering
      msqv_cs(0,:,:)=zip
      msqv_cs(2,:,:)=zip

c return msqv array in usual format by performing sum
      msqv(:,:)=msqv_cs(0,:,:)+msqv_cs(1,:,:)+msqv_cs(2,:,:)

      return
      end

