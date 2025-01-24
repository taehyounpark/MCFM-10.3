!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine htautaudecay(p,jm,jp,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: J.M. Campbell, August 2012                               *
c                                                                      *
c     matrix element squared for the process of                        *
c     Higgs decay  H --> tau^-(jm)+tau^+(jp)                           *
c     with tau mass included                                           *
c***********************************************************************

      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      integer:: jm,jp
      real(dp):: p(mxpart,4),s56,msq,msqhtautau

      s56=two*(p(jm,4)*p(jp,4)-p(jm,1)*p(jp,1)
     &        -p(jm,2)*p(jp,2)-p(jm,3)*p(jp,3))+two*mtau**2

      msq=msqhtautau(s56)

      return
      end



      function msqhtautau(s)
      implicit none
      include 'types.f'
      real(dp):: msqhtautau
      real(dp):: s
      include 'constants.f'
      include 'masses.f'
      include 'ewcouple.f'

      msqhtautau=gwsq*mtausq/(four*wmass**2)*two*(s-four*mtau**2)

      return
      end


