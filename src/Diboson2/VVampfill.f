!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine VVampfill(order,ABCF,qtype,p1,p2,p5,p6,p7,p8,qqbAj,amp)
c     applies coupling to amplitudes and removes factor of four
      implicit none
      include 'types.f'
      include 'coupfac_vv.f'
      include 'ABCF.f'
      complex(dp):: qqbVVMLLL,qqbVVMRLL,amp(2,2,2),qqbAj(0:2,4,10)
      integer p1,p2,p5,p6,p7,p8
      integer:: order,qtype
      amp(1,1,1)=coupfac(ABCF,qtype,1,1,1)
     & *qqbVVMLLL(order,ABCF,qtype,p1,p2,p5,p6,p7,p8,qqbAj) ! LLL
      amp(1,1,2)=coupfac(ABCF,qtype,1,1,2)
     & *qqbVVMLLL(order,ABCF,qtype,p1,p2,p5,p6,p8,p7,qqbAj) ! LLR
      amp(1,2,1)=coupfac(ABCF,qtype,1,2,1)
     & *qqbVVMLLL(order,ABCF,qtype,p1,p2,p6,p5,p7,p8,qqbAj) ! LRL
      amp(1,2,2)=coupfac(ABCF,qtype,1,2,2)
     & *qqbVVMLLL(order,ABCF,qtype,p1,p2,p6,p5,p8,p7,qqbAj) ! LRR

      amp(2,1,1)=coupfac(ABCF,qtype,2,1,1)
     & *qqbVVMRLL(order,ABCF,qtype,p1,p2,p5,p6,p7,p8,qqbAj) ! RLL
      amp(2,2,1)=coupfac(ABCF,qtype,2,2,1)
     & *qqbVVMRLL(order,ABCF,qtype,p1,p2,p6,p5,p7,p8,qqbAj) ! RRL
      amp(2,2,2)=coupfac(ABCF,qtype,2,2,2)
     & *qqbVVMRLL(order,ABCF,qtype,p1,p2,p6,p5,p8,p7,qqbAj) ! RRR
      amp(2,1,2)=coupfac(ABCF,qtype,2,1,2)
     & *qqbVVMRLL(order,ABCF,qtype,p1,p2,p5,p6,p8,p7,qqbAj) ! RLR
c     Divide by a factor of 4 in order to agree with normalization
c     used e.g. in qq_zz.f
      amp(:,:,:)=0.25_dp*amp(:,:,:)
      return
      end
