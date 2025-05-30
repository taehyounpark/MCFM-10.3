!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c---------------------------------------------------------------
c   This subroutine checks the number of external dipoles----
c---absorbing the correct number into the fragmenation functions
c   it then returns the finite (msq_qcd*dip) ---------------------
c---------------------------------------------------------------
c--- Author C. Williams Nov 2012
c-----------------------------------------------------------------
c modified under frag re-write March 2013

      subroutine qqb_dm_monophot_fragdips(p,p_phys,qcd_tree,msq_out)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'frag.f'
      real(dp):: p(mxpart,4),p_phys(mxpart,4)
      real(dp):: msq_qcd(-nf:nf,-nf:nf),msq_out(-nf:nf,-nf:nf)
      integer:: j,k
      real(dp):: virt_dips,xl,dot,fsq
      real(dp):: aewo2pi,fi_gaq
      external qcd_tree

      aewo2pi=esq/(fourpi*twopi)

      fsq=frag_scale**2

      xl=log(-two*dot(p_phys,2,5)/fsq)

      virt_dips=+aewo2pi*(fi_gaq(z_frag,p_phys,xl,5,2,2))

      do j=-nf,nf
         do k=-nf,nf
            msq_qcd(j,k)=0d0
            msq_out(j,k)=0d0
         enddo
      enddo

      call qcd_tree(p,msq_qcd)

      do j=-nf,nf
         do k=-nf,nf

           if((j==0).and.(k /= 0)) then
              msq_out(j,k)=Q(k)**2*msq_qcd(j,k)*virt_dips
           elseif((j /= 0).and.(k==0)) then
              msq_out(j,k)=Q(j)**2*msq_qcd(j,k)*virt_dips
           endif

         enddo
      enddo


      return
      end subroutine



