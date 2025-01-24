!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine hbbdecay_g(p,ib,ibb,ig,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: J.M. Campbell, June 2012                                 *
c                                                                      *
c     matrix element squared for the process of                        *
c     Higgs decay  H --> b(ib)+b~(ibb)+g(ig)                           *
c     with bottom mass included                                        *
c***********************************************************************
      include 'mxpart.f'
      include 'masses.f'
      integer:: ib,ibb,ig,j,k
      real(dp):: p(mxpart,4),s,s56,s57,s67,msq,msqhbbg

      s(j,k)=2*(p(j,4)*p(k,4)-p(j,1)*p(k,1)-p(j,2)*p(k,2)-p(j,3)*p(k,3))
      s56=s(ib,ibb)+2._dp*mb**2
      s57=s(ib,ig)
      s67=s(ibb,ig)

      msq=msqhbbg(s56,s57,s67)
      return
      end

      function msqhbbg(s56,s57,s67)
      implicit none
      include 'types.f'
      real(dp):: msqhbbg
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'hbbparams.f'
      real(dp):: s56,s57,s67

      msqhbbg=xn*gwsq*mb_eff**2/(4._dp*wmass**2)*gsq*Cf
     &*(2._dp*(s56-4._dp*mb**2)
     &*(-4._dp*mb**2/s57**2-4._dp*mb**2/s67**2
     &  +4._dp*(s56-2._dp*mb**2)/s57/s67)
     & +4._dp*(2._dp*s56+s57+s67-6._dp*mb**2)/s57
     & +4._dp*(2._dp*s56+s57+s67-6._dp*mb**2)/s67
     & -8._dp*mb**2*(s67/s57**2+s57/s67**2))

      return
      end
