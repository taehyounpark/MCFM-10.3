!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine setupqt()
      implicit none
      include 'types.f'
      include 'taucut.f'
c      include 'nproc.f'

      qtcut = 0.0015_dp
c      qtcut = 0.006_dp

      return

c      select case (nproc)
c          case (1,6)
c              qtcut = 6.e-3_dp
c          case (31,32)
c              qtcut = 6.e-3_dp
c          case (61,81,82)
c              qtcut = 3.e-4_dp
c          case (71,76)
c              qtcut = 1.5e-4_dp
c          case (91:110)
c              qtcut = 3.e-3_dp
c          case (111,112,119)
c              qtcut = 4.e-3_dp
c          case (203:204)
c              qtcut = 1.e-4_dp
c          case (285)
c              qtcut = 1.e-4_dp
c          case (290,295,300,305)
c              qtcut = 3.e-4_dp
c          case default
c              write(6,*) 'Pre-determined value of qtcut not available'
c              write(6,*) 'for nproc = ',nproc
c              stop
c      end select

      end

