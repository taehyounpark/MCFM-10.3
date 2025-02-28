!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_dm_monojet_v(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'scheme.f'
      include 'dm_params.f'
c     (Taken from Ravindran, Smith, van Neerven hep-ph/0201114)
c     Modified by overall factors
      integer:: iglue,j,k
      real(dp):: p(mxpart,4),msq(fn:nf,fn:nf)
      real(dp):: ss,tt,uu,
     & virtgg,virtqa,virtaq,virtqg,virtgq,hdecay,Asq,fac
      parameter(iglue=5)

      scheme='tH-V'

      call dotem(iglue,p,s)
      ss=s(1,2)
      tt=s(1,iglue)
      uu=s(2,iglue)

c      Asq=(as/(3d0*pi))**2/vevsq
      Asq=one/dm_lam**6*as**2*16d0

      call dmsdecay(p,3,4,hdecay)
c      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      fac=ason2pi*Asq*gsq*hdecay
      call hjetfill(ss,tt,uu,virtgg,virtqa,virtaq,virtqg,virtgq)

      do j=-nf,nf
      do k=-nf,nf
      msq(j,k)=0d0
      if ((j==0).and.(k==0)) msq(j,k)=avegg*fac*virtgg
      if ((j>0).and.(k==-j)) msq(j,k)=aveqq*fac*virtqa
      if ((j<0).and.(k==-j)) msq(j,k)=aveqq*fac*virtaq
      if ((j==0).and.(k /= 0)) msq(j,k)=aveqg*fac*virtgq
      if ((j /= 0).and.(k==0)) msq(j,k)=aveqg*fac*virtqg
      enddo
      enddo

      return
      end
