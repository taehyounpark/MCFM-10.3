!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_gamgam_g(p,msq)
      implicit none
      include 'types.f'

c***********************************************************************
c     Authors: R.K. Ellis and John M. Campbell                         *
c     December, 2010.                                                  *
c***********************************************************************
c                                                                      *
c     Matrix element for gamma + gamma + 1 parton production           *
c     in order alpha_s, averaged over initial colours and spins        *
c                                                                      *
c     q(-p1)+qbar(-p2) --> gamma(p3)+gamma(p4)+gluon(p5)               *
c                                                                      *
c     (and all crossings)                                              *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     & Cgamgam,qa_gagag(2),
     & ag_gagaa(2),qg_gagaq(2),ga_gagaa(2),gq_gagaq(2),gg_gaga
      integer:: j,k
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      call spinoru(5,p,za,zb)

      do j=1,2
      ag_gagaa(j)=-aveqg*Cgamgam(5,1,3,4,2,j)
      qg_gagaq(j)=-aveqg*Cgamgam(1,5,3,4,2,j)
      ga_gagaa(j)=-aveqg*Cgamgam(5,2,3,4,1,j)
      gq_gagaq(j)=-aveqg*Cgamgam(2,5,3,4,1,j)

      qa_gagag(j)=+aveqq*Cgamgam(1,2,3,4,5,j)
      enddo

c--- averaging factor already included in ggtogagag

      gg_gaga=zip
c JC removing gg contribution
c      if (omitgg) then
c        gg_gaga=0._dp
c      else
c        gg_gaga=ggtogagag()
c      endif

      do j=-nf,nf
      do k=-nf,nf

c--set msq=0 to initalize
      msq(j,k)=0._dp
c--qa
      if ((j > 0) .and. (k < 0)) then
          if (j == -k) then
          msq(j,k)=qa_gagag(jj(j))
          endif
c--aq
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) then
          msq(j,k)=qa_gagag(kk(k))
          endif
c--qg
      elseif ((j > 0) .and. (k == 0)) then
            msq(j,k)=qg_gagaq(jj(j))
c--ag
      elseif ((j < 0) .and. (k == 0)) then
            msq(j,k)=ag_gagaa(-jj(j))
c--gq
      elseif ((j == 0) .and. (k > 0)) then
            msq(j,k)=gq_gagaq(kk(k))
c--ga
      elseif ((j == 0) .and. (k < 0)) then
            msq(j,k)=ga_gagaa(-kk(k))

c--gg
      elseif ((j == 0) .and. (k == 0)) then
            msq(j,k)=gg_gaga
      endif

      enddo
      enddo

      return
      end
