!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine setuptau()
      implicit none
      include 'types.f'
      include 'taucut.f'
      include 'nproc.f'

      if (nproc == 1610) then
         taucut = 0.001_dp
      else
         taucut = 0.0015_dp**sqrt(2._dp)
c         taucut = 0.006_dp**sqrt(2._dp)
      endif

      return

      select case (nproc)
          case (1,6)
              taucut = 6.e-3_dp
          case (31,32)
              taucut = 6.e-3_dp
          case (61,81,82)
              taucut = 3.e-4_dp
          case (71,76)
              taucut = 1.5e-4_dp
          case (91:110)
              taucut = 3.e-3_dp
          case (111,112,119)
              taucut = 4.e-3_dp
          case (1610)
              taucut = 0.001_dp
          case (203:204)
              taucut = 1.e-4_dp
          case (285)
              taucut = 1.e-4_dp
          case (290,295,300,305)
              taucut = 3.e-4_dp
          case default
             taucut = 0.0015_dp**sqrt(2._dp)
             write(6,*) 'Pre-determined value of taucut not available'
             write(6,*) 'for nproc = ',nproc
             stop
      end select

      end
