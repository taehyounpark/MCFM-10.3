!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function ubdggmsq(p1,p2,p3,p4,p5,p6,p7,za,zb)
c     Matrix element for
c     ub(-p1)+d(-p2)=e-(p3)+nu~(p4)+gamma(p5)+g(p6)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer:: p1,p2,p3,p4,p5,p6,p7
      real(dp):: ubdggmsq
      complex(dp):: AB(2,2,2),BA(2,2,2)

      call a7Wgamg_g(p1,p2,p3,p4,p5,p6,p7,za,zb,AB,BA)

      ubdggmsq=
     & +cdabs(AB(1,1,1))**2+cdabs(AB(1,1,2))**2
     & +cdabs(AB(1,2,1))**2+cdabs(AB(1,2,2))**2
     & +cdabs(AB(2,2,2))**2+cdabs(AB(2,2,1))**2
     & +cdabs(AB(2,1,2))**2+cdabs(AB(2,1,1))**2
     & +cdabs(BA(1,1,1))**2+cdabs(BA(1,1,2))**2
     & +cdabs(BA(1,2,1))**2+cdabs(BA(1,2,2))**2
     & +cdabs(BA(2,2,2))**2+cdabs(BA(2,2,1))**2
     & +cdabs(BA(2,1,2))**2+cdabs(BA(2,1,1))**2
     & -(
     & +cdabs(AB(1,1,1)+BA(1,1,1))**2
     & +cdabs(AB(1,1,2)+BA(1,1,2))**2
     & +cdabs(AB(1,2,1)+BA(1,2,1))**2
     & +cdabs(AB(1,2,2)+BA(1,2,2))**2
     & +cdabs(AB(2,2,2)+BA(2,2,2))**2
     & +cdabs(AB(2,2,1)+BA(2,2,1))**2
     & +cdabs(AB(2,1,2)+BA(2,1,2))**2
     & +cdabs(AB(2,1,1)+BA(2,1,1))**2)/xn**2


      return
      end

