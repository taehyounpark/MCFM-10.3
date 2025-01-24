!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_gamgam_v(p,msq)
      implicit none
      include 'types.f'

c***********************************************************************
c     Authors: R.K. Ellis and John M. Campbell                         *
c     December, 2010.                                                  *
c***********************************************************************
c                   &                                                  *
c     Matrix element for gamma + gamma production,                     *
c     averaged over initial colours and spins                          *
c                   &                                                  *
c     q(-p1)+qbar(-p2) --> gamma(p3)+gamma(p4)                         *
c                   &                                                  *
c***********************************************************************
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zcouple_cms.f'
      include 'scheme.f'
      integer:: j,k
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: qa,aq,gg,fac,qagamgam
      real(dp), parameter:: statfac=0.5_dp
c      real(dp):: facgg,Qsum

c     scheme='tH-V'
      scheme='dred'

      fac=8._dp*xn*abs(zesq)**2*ason2pi*statfac

      call dotem(4,p,s)
      qa=+qagamgam(1,2,3)*fac*aveqq
      aq=qa


c JC removing gg contribution
c--- initialize gg 2-loop matrix elements
c      Qsum=+Q(1)**2+Q(2)**2+Q(3)**2+Q(4)**2+Q(5)**2
c      facgg=4._dp*abs(zesq)*gsq/(16._dp*pisq)*Qsum
c      gg=avegg*V*facgg**2*statfac*virtgamgam(s(1,2),s(1,3),s(2,3))
      gg=zip

      do j=-nf,nf
      do k=-nf,nf
c--set msq=0 to initalize
      msq(j,k)=0._dp
c--qa
      if ((j > 0) .and. (k < 0)) then
          if (j == -k) msq(j,k)=Q(j)**4*qa
c--aq
      elseif ((j < 0) .and. (k > 0)) then
          if (j == -k) msq(j,k)=Q(k)**4*aq
c--gg
      elseif ((j == 0) .and. (k == 0)) then
          msq(j,k)=gg
      endif

      enddo
      enddo

      return
      end


      function qagamgam_cw(i1,i2,i3)
c===== corrected version of q(i1)+qb(i2)=>gamma(i3)+gamma(i4)
c====== C. Williams July 2015
      implicit none
      include 'types.f'
      real(dp):: qagamgam_cw

      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      integer:: i1,i2,i3
      complex(dp):: lnrat,l12,l13,l23
      real(dp):: s12,s13,s23,T0
      s12=s(i1,i2)
      s13=s(i1,i3)
      s23=s(i2,i3)
      T0=s13/s23+s23/s13
      l12=lnrat(musq,-s12)
      l13=lnrat(musq,-s13)
      l23=lnrat(musq,-s23)

c==== divergent parts
      qagamgam_cw=(-2*epinv*epinv2-(3._dp+2*real(l12))*epinv)*T0

c======= finite part
      qagamgam_cw=qagamgam_cw+
     & real((2*l12**2 + l13**2 + l23**2 - 2*l12*(l13 + l23) +
     &       2*Pi**2)*s12**2 +
     &    (-6 + l13**2 - 2*l12*(1 + l13) - l23 + Pi**2)*s13**2 +
     &    (-6 - l13 + l23**2 - 2*l12*(1 + l23) + Pi**2)*s23**2 +
     &    2*s12*((-l12 + l23)*s13 + (-l12 + l13)*s23))/(s13*s23)

c===== normalization
      qagamgam_cw=cf*qagamgam_cw
      return
      end

      function qagamgam(i1,i2,i3)
      implicit none
      include 'types.f'
      real(dp):: qagamgam

c----Matrix element for gamma + gamma production
c----in order alpha_s
c---
c---  0 -> gamma(p1) + gamma(p2) + q(p3) + qb(p4)
c---
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scale.f'
      include 'scheme.f'
      integer:: i1,i2,i3
      complex(dp):: lnrat,l12,l13,l23
      real(dp):: s12,s13,s23,T0,deltar
      s12=s(i1,i2)
      s13=s(i1,i3)
      s23=s(i2,i3)
      T0=s13/s23+s23/s13
      l12=lnrat(-s12,musq)
      l13=lnrat(-s13,musq)
      l23=lnrat(-s23,musq)

      if     (scheme == 'tH-V') then
        deltar=0._dp
      elseif (scheme == 'dred') then
        deltar=1._dp
      else
        write(6,*) 'Invalid scheme in qqb_gamgam_v.f'
        stop
      endif

cId,anscdr=cf*((-2*epinv*epinv2-epinv*(3-2*[ln(-s12)])-7-[ln(-s12)]**2)*T0
c    +[ln(-s13)]*(s23-2*s12)/s13
c    +[ln(-s23)]*(s13-2*s12)/s23
c    +(s12**2+s23**2)/s13/s23*(([ln(-s12)]-[ln(-s23)])**2+pisq)
c    +(s12**2+s13**2)/s13/s23*(([ln(-s12)]-[ln(-s13)])**2+pisq)
c    -4*[ln(-s12)]);


      qagamgam=cf*((-2._dp*epinv*epinv2-epinv*(3._dp-2._dp*real(l12))
     & -7._dp+deltar-real(l12**2))*T0
     & +real(l13)*(s23-2._dp*s12)/s13
     & +real(l23)*(s13-2._dp*s12)/s23
     & +(s12**2+s23**2)/s13/s23*(real((l12-l23)**2)+pisq)
     & +(s12**2+s13**2)/s13/s23*(real((l12-l13)**2)+pisq)
     & -4._dp*real(l12))

      return
      end


