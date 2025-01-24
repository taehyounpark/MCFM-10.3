!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine gg_2gam_z(p,z)
      implicit none
      include 'types.f'
c***********************************************************************
c     John M. Campbell                                                 *
c     February, 2016.                                                  *
c***********************************************************************
      include 'constants.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'PR_new.f'
      include 'agq.f'
      integer:: is
      real(dp):: z,xl12,p(mxpart,4),dot,ii_gg,ii_gq,tempgg,tempgq

      xl12=log(two*dot(p,1,2)/musq)

      do is=1,3
      tempgg=ason2pi*xn*ii_gg(z,xl12,is)
      Q1(g,g,g,is)=tempgg
      Q2(g,g,g,is)=tempgg

      tempgq=ason4pi*two*cf*ii_gq(z,xl12,is)
      Q1(g,q,g,is)=tempgq
      Q1(g,a,g,is)=tempgq
      Q2(g,q,g,is)=tempgq
      Q2(g,a,g,is)=tempgq
      enddo

      return
      end
