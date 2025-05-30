!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c---------------------------------------------------------------
c   This subroutine checks the number of external dipoles----
c---absorbing the correct number into the fragmenation functions
c   it then returns the finite (msq_qcd*dip) ---------------------
c---------------------------------------------------------------
c--- Author C. Williams Feb 2011
c-----------------------------------------------------------------

c--- Passed in
c---   p:      array of momenta to evaluate matrix elements
c---   p_phys: array of momenta to evaluate integrated dipoles

      subroutine qqb_gamgam_fragdips(p,p_phys,qcd_tree,msq_out)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcharge.f'
      include 'frag.f'
      include 'lastphot.f'
      include 'zcouple_cms.f'
      real(dp):: p(mxpart,4),p_phys(mxpart,4)
      real(dp):: msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      integer:: j,k
      real(dp):: virt_dips,xl,dot,fsq
      real(dp):: aewo2pi,fi_gaq
      external qcd_tree


      aewo2pi=abs(zesq)/(fourpi*twopi)
      fsq=frag_scale**2

      xl=log(-two*dot(p_phys,2,lastphot)/fsq)
      virt_dips=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,lastphot,2,2))

      msq_out(:,:)=zip
c--- fill underlying QCD matrix elements
      call qcd_tree(p,msq_qcd)

c--- fill output array
      do j=-nf,nf
         do k=-nf,nf

c   factor of two cancelled by statistical factor because two photons
            if((j==0).and.(k /= 0)) then
                  msq_out(j,k)=msq_qcd(j,k)*Q(k)**2*virt_dips
            elseif((j /= 0).and.(k==0)) then
                 msq_out(j,k)=msq_qcd(j,k)*Q(j)**2*virt_dips
            endif

         enddo
      enddo


      return
      end subroutine



