!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine dkqqb_hww_gs(p,msqc)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     May, 2012.                                                       *
c     calculate the subtraction term for radiation in the              *
c     W^- decay for the process                                        *
c                                                                      *
c     q(-p1)+qbar(-p2) = H-> W^+(nu(p3)+e+(p4))                        *
c                       +W^-(e-(p5)+nubar(p6)+g(p7))                   *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ptilde.f'
      include 'qqgg.f'
      include 'incldip.f'
      real(dp):: msqc(maxd,-nf:nf,-nf:nf),p(mxpart,4)
      real(dp)::
     & msq57_6(-nf:nf,-nf:nf),msq67_5(-nf:nf,-nf:nf),
     & dummyv(-nf:nf,-nf:nf),sub57_6(4),sub67_5(4),dsubv
      integer:: j,k,nd
      external qqb_hww,donothing_gvec

      ndmax=2

      do j=-nf,nf
      do k=-nf,nf
      do nd=1,ndmax
        msqc(nd,j,k)=0._dp
        incldip(nd)=.true.
      enddo
      enddo
      enddo

      call dips(1,p,5,7,6,sub57_6,dsubv,msq57_6,dummyv,
     & qqb_hww,donothing_gvec)
      call dips(2,p,6,7,5,sub67_5,dsubv,msq67_5,dummyv,
     & qqb_hww,donothing_gvec)

      msqc(1,0,0)=sub57_6(qq)*msq57_6(0,0)*2._dp*cf
      msqc(2,0,0)=sub67_5(qq)*msq67_5(0,0)*2._dp*cf


      return
      end

