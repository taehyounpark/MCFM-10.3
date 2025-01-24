!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_gamgamas2(p,tlrem,msq)
      implicit none
      include 'types.f'
c---
c--- O(as^2) contributions that do not factorize over LO matrix element
c---
c---Matrix element squared averaged over initial colors and spins
c---for diphoton production
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'facscale.f'
      include 'ckm.f'
      include 'scet_const.f'
      include 'hbbparams.f'
      include 'noglue.f'
      include 'ewcharge.f'
      include 'zcouple_cms.f'
      include 'sprods_com.f'
      integer:: j,k
      real(dp):: p(mxpart,4),fac,msq(-nf:nf,-nf:nf),
     & gggaga,msqgggaga,facgg,qsum,tlrem
      real(dp),parameter:: statfac=0.5_dp

      fac=xn*abs(zesq)**2*statfac

      qsum=zip
      do j=1,nf
         qsum=qsum+Q(j)**2
      enddo
      facgg=4._dp*abs(zesq)*gsq/(16._dp*pisq)*Qsum
      gggaga=avegg*V*facgg**2*msqgggaga(s(1,2),s(1,3),s(2,3))*statfac

      msq(:,:)=zip
      do j=-nf,nf
         k=-j

         if (j*k > 0) cycle    ! skip, qq, aa

         if (j > 0) then
c======add on two-loop A functions which go like sum over quark charges
           msq(j,k)=fac*aveqq*tlrem*ason2pi**2*qsum*Q(j)**2
         elseif (j < 0) then
c======add on two-loop A functions which go like sum over quark charges
           msq(j,k)=fac*aveqq*tlrem*ason2pi**2*qsum*Q(k)**2
         elseif (j==0) then
           msq(j,k)=gggaga
         endif


      enddo

      return
      end


