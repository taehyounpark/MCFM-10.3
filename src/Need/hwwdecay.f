!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine hwwdecay(p,j1,j2,j3,j4,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: R.K. Ellis, September 2012                               *
c                                                                      *
c     matrix element squared for the process of                        *
c     Higgs decay  H --> WW --> ne(j1)+e+(j2)+e-(j3)+nu~(j4)           *
c***********************************************************************
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      integer:: j1,j2,j3,j4
      real(dp):: p(mxpart,4),msq,dot,hdecay,s12,s13,s24,s34

      s12=2._dp*dot(p,j1,j2)
      s13=2._dp*dot(p,j1,j3)
      s24=2._dp*dot(p,j2,j4)
      s34=2._dp*dot(p,j3,j4)

      hdecay=gwsq**3*wmass**2*s13*s24
      hdecay=hdecay
     &  /(((s12-wmass**2)**2+(wmass*wwidth)**2)
     &   *((s34-wmass**2)**2+(wmass*wwidth)**2))
      msq=hdecay

      return
      end


