!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine dkqqb_wh_v_massless(p,msq)
 !==== C.Williams Sept 2015
c===== routine which calculates qqb => W (H=bb) with Hbb at NLO (virtual)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'epinv.f'
      include 'epinv2.f'
      include 'scheme.f'
      include 'scale.f'
      include 'alfacut.f'
      include 'nf.f'
      include 'qcdcouple.f'
      real(kind=dp) p(mxpart,4),msq(-nf:nf,-nf:nf),msq0(-nf:nf,-nf:nf)
      real(kind=dp) ct,ls,L,ff
      real(kind=dp) s12,subuv
      integer j,k

      msq(:,:)=zip
      s12=two*(p(1,4)*p(2,4)
     &        -p(1,3)*p(2,3)-p(1,2)*p(2,2)-p(1,1)*p(2,1))

      scheme='dred'
c==== LO
      call qqb_wh(p,msq0)
      ls=log(musq/s12)
      L=-ls
c==== virtual form factor  (in units of ason4pi)
      ff=(-2._dp - 2._dp*epinv*epinv2 - (2._dp*ls)*epinv
     &     - ls**2 + 7.d0*Pi**2/6._dp)
      ff=ff*cf

c==== final final integrated counter term  (in units of ason4pi)
      ct=epinv*(epinv2-L)+half*L**2+1.5_dp*(epinv-L)+5._dp-half*pisq
      ct=ct+1.5_dp*(aff-1._dp-log(aff))-log(aff)**2

c==== factor of 2 for each dipole and scheme def
      if(scheme=='dred') ct=ct-half
      ct=2._dp*ct*cf

c===== UV renormalization of vertex (in units of ason4pi)
      subuv=-3._dp*cf*epinv

c==== total factor (ason2pi since 2*Re)
      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=ason2pi*(ff+ct+subuv)*msq0(j,k)
         enddo
      enddo


      return
      end


