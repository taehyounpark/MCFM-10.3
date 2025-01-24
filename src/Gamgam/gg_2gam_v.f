!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_2gam_v(p,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Authors: R.K. Ellis and John M. Campbell                         *
c     December, 2010.                                                  *
c***********************************************************************
c                                                                      *
c     Matrix element for gamma + gamma production,                     *
c     averaged over initial colours and spins                          *
c                                                                      *
c     g(-p1)+g(-p2) --> gamma(p3)+gamma(p4)                            *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),facgg,Qsum,
     & virtgamgam
      real(dp),parameter::statfac=0.5_dp

c--set msq=0 to initalize
      msq(:,:)=0._dp

c      scheme='dred'
c This bug pointed out by Xiaran Zhao, Aug 21 2018
      scheme='tH-V'

      call dotem(3,p,s)

c--- initialize gg 2-loop matrix elements
      Qsum=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
      facgg=4._dp*esq*gsq/(16._dp*pisq)*Qsum

      msq(0,0)=avegg*V*facgg**2*statfac*virtgamgam(s(1,2),s(1,3),s(2,3))

      return
      end


      function gg_2gam_msbar(p)
      implicit none
      include 'types.f'
c***********************************************************************
c     Authors: R.K. Ellis and John M. Campbell                         *
c     December, 2010.                                                  *
c***********************************************************************
c                                                                      *
c     Matrix element for gamma + gamma production,                     *
c     averaged over initial colours and spins                          *
c                                                                      *
c     g(-p1)+g(-p2) --> gamma(p3)+gamma(p4)                            *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      real(dp) :: gg_2gam_msbar
      real(dp), intent(in) :: p(mxpart,4)
      real(dp) :: facgg,Qsum, virtgamgam_msbar
      real(dp),parameter::statfac=0.5_dp

      call dotem(3,p,s)

      Qsum=Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
      facgg=4._dp*esq*gsq/(16._dp*pisq)*Qsum

      gg_2gam_msbar = avegg*V*facgg**2*statfac
     & *virtgamgam_msbar(s(1,2),s(1,3),s(2,3))

      end


