!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine swapjet(pjet,jetindex,i,j)
        use types
        use jettagging
        use SCET, only: currentNd
      implicit none
c--- swaps jets i..j in pjet

      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'cplx.h'
      include 'jetlabel.f'
      integer:: i,j,k,jetindex(mxpart)
      real(dp):: pjet(mxpart,4),tmp
      integer :: itmp
      character(len=2):: chartmp

c--- escape if we're trying to swap the same jets
      if (i == j) return

      do k=1,4
        tmp=pjet(i,k)
        pjet(i,k)=pjet(j,k)
        pjet(j,k)=tmp
      enddo

      chartmp=jetlabel(i)
      jetlabel(i)=jetlabel(j)
      jetlabel(j)=chartmp

      itmp = jetcontent(i,currentNd)
      jetcontent(i,currentNd) = jetcontent(j,currentNd)
      jetcontent(j,currentNd) = itmp

      return
      end
