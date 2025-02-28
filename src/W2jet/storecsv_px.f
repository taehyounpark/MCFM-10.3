!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine storecsv_px(i,j)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the W2jet_gvec matrix elements into elements (..,i,j) of p1p2
      include 'mmsqv_cs.f'
      integer:: i,j,k
      real(dp):: p1p2(0:2,-1:1,-1:1)
      common/p1p2/p1p2
!$omp threadprivate(/p1p2/)

      do k=0,2
        p1p2(k,i,j)=mmsqv_cs(k,+1,+1)
      enddo

      return
      end


      subroutine storecsv_px_swap(i,j)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the W2jet_gvec matrix elements into elements (..,i,j) of p1p2
      include 'mmsqv_cs.f'
      integer:: i,j,k
      integer, save:: swap(0:2) = (/ 0, 2, 1 /)
      real(dp):: p1p2(0:2,-1:1,-1:1)
      common/p1p2/p1p2
!$omp threadprivate(/p1p2/)

      do k=0,2
        p1p2(k,i,j)=mmsqv_cs(swap(k),+1,+1)
      enddo

      return
      end

