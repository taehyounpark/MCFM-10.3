!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_hgagag_v(p,msq)
      implicit none
      include 'types.f'
c----Author: R.K. Ellis August 2011
c----Virtual corrections matrix element for H production
c----in the heavy quark (mt=Infinity) limit.
c----averaged over initial colours and spins
c----f(p1)+f(p2) --> gamma(p3)+gamma(p4)+g(p5)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
c     (Taken from Ravindran, Smith, van Neerven hep-ph/0201114)
c     Modified by overall factors
      integer:: iglue,j,k
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf)
      real(dp):: ss,tt,uu,msqgamgam,
     & virtgg,virtqa,virtaq,virtqg,virtgq,hdecay,Asq,fac
      parameter(iglue=5)

      scheme='tH-V'

      call dotem(iglue,p,s)
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)

      Asq=(as/(3._dp*pi))**2/vevsq

c   Deal with Higgs decay
      hdecay=msqgamgam(hmass)/((s(3,4)-hmass**2)**2+(hmass*hwidth)**2)

      fac=ason2pi*Asq*gsq*hdecay
      call hjetfill(ss,tt,uu,virtgg,virtqa,virtaq,virtqg,virtgq)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0._dp
      if ((j==0).and.(k==0)) msq(j,k)=avegg*fac*virtgg
      if ((j>0).and.(k==-j)) msq(j,k)=aveqq*fac*virtqa
      if ((j<0).and.(k==-j)) msq(j,k)=aveqq*fac*virtaq
      if ((j==0).and.(k /= 0)) msq(j,k)=aveqg*fac*virtgq
      if ((j /= 0).and.(k==0)) msq(j,k)=aveqg*fac*virtqg
      enddo
      enddo

      return
      end
