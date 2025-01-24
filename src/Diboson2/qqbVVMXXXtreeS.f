!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqbVVMXXXtreeS(switch,p1,p2,p3,p4,p5,p6,qqb1)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      integer p1,p2,p3,p4,p5,p6,n1,n2,n3
      complex(dp)::qqb1(2,2,2),tmp(2,2,2),A6trees
      logical:: switch
c     A6 spinor content A6trees(1,2,3,4,5,6) --> <1><3><6>[2],[4],[5]
      tmp(2,1,1)=+A6trees(p3,p4,p1,p2,p6,p5,za,zb) !<3,1,5>[4,2,6]
      tmp(2,2,1)=-A6trees(p4,p3,p1,p2,p6,p5,za,zb) !<4,1,5>[3,2,6]
      tmp(2,1,2)=+A6trees(p3,p4,p1,p2,p5,p6,za,zb) !<3,1,6>[4,2,5]
      tmp(2,2,2)=-A6trees(p4,p3,p1,p2,p5,p6,za,zb) !<4,1,6>[3,2,5]

      tmp(1,1,1)=+A6trees(p3,p4,p2,p1,p6,p5,za,zb) ! <3><2><5> [4],[1],[6]
      tmp(1,2,1)=-A6trees(p4,p3,p2,p1,p6,p5,za,zb) ! <4><2><5> [3],[1],[6]
      tmp(1,1,2)=+A6trees(p3,p4,p2,p1,p5,p6,za,zb) ! <3><2><6> [4],[1],[5]
      tmp(1,2,2)=-A6trees(p4,p3,p2,p1,p5,p6,za,zb) ! <4><2><6> [3],[1],[5]

      if (switch) then ! if we are switching (56)<-->(78)
         do n3=1,2; do n2=1,2; do n1=1,2
         qqb1(n1,n2,n3)=tmp(n1,n3,n2)
         enddo; enddo; enddo
      else
         qqb1(:,:,:)=tmp(:,:,:)
      endif

      return
      end
