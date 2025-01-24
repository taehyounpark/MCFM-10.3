!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine A7Wpgamgamg(p1,p2,p3,p4,p5,p6,p7,za,zb,b)
      implicit none
c---- Matrix element for Wgamgamma radiation from line 12
c  u(-p1) +dbar(-p2) --> ve(p3)+e^+(p4)+gamgam(p5)+g(p6)
c  including radiation from decay electron

c Overall factor of  i_*e^2*gw^2*gs
c     including QED type propagators
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5,p6,p7,b5,b6,b7
      complex(dp)::b(2,2,2)
c      amplitude b(h5,h6,h7)

      b5=p1
      b6=p1
      b7=p1
c      do b5=1,1!2
c      do b6=1,1!2
c      do b7=1,1!2
      call A7Wpgamgamg_mmm(p1,p2,p3,p4,p5,p6,p7,b5,b6,b7,za,zb,b(1,1,1))
      call A7Wpgamgamg_mmp(p1,p2,p3,p4,p5,p6,p7,b5,b6,b7,za,zb,b(1,1,2))
      call A7Wpgamgamg_mpm(p1,p2,p3,p4,p5,p6,p7,b5,b6,b7,za,zb,b(1,2,1))
      call A7Wpgamgamg_mpp(p1,p2,p3,p4,p5,p6,p7,b5,b6,b7,za,zb,b(1,2,2))
      call A7Wpgamgamg_pmm(p1,p2,p3,p4,p5,p6,p7,b5,b6,b7,za,zb,b(2,1,1))
      call A7Wpgamgamg_pmp(p1,p2,p3,p4,p5,p6,p7,b5,b6,b7,za,zb,b(2,1,2))
      call A7Wpgamgamg_ppm(p1,p2,p3,p4,p5,p6,p7,b5,b6,b7,za,zb,b(2,2,1))
      call A7Wpgamgamg_ppp(p1,p2,p3,p4,p5,p6,p7,b5,b6,b7,za,zb,b(2,2,2))
c      write(6,*) b5,b6,b7
c      write(6,*) 'b(1,1,1)',b(1,1,1)
c      write(6,*) 'b(1,1,2)',b(1,1,2)
c      write(6,*) 'b(1,2,1)',b(1,2,1)
c      write(6,*) 'b(2,1,1)',b(2,1,1)
c      write(6,*) 'b(2,2,2)',b(2,2,2)
c      write(6,*) 'b(2,2,1)',b(2,2,1)
c      write(6,*) 'b(2,1,2)',b(2,1,2)
c      write(6,*) 'b(1,2,2)',b(1,2,2)
c      enddo
c      enddo
c      enddo
c      pause
      return
      end
