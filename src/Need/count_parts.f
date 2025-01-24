!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      ! functions for counting particles
      function count_photo()
      implicit none
      integer:: count_photo
      include 'mxpart.f'
      integer:: j
      logical:: is_photon
      external is_photon

      count_photo=0
      do j=1,mxpart
         if (is_photon(j)) then
            count_photo=count_photo+1
         endif
      enddo

      return
      end


      function count_jets()
      implicit none
      integer:: count_jets
      include 'npart.f'
      integer:: j
      logical:: is_hadronic
      external is_hadronic

c---- Count final state jets
      count_jets=0

      do j=3,npart+2
         if (is_hadronic(j)) then
            count_jets=count_jets+1
         endif
      enddo

      return
      end
