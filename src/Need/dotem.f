!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine dotem(N,p,s)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      real(dp):: p(mxpart,4),s(mxpart,mxpart),dot
      integer:: j,k,N
c---returns 2*piDpj for massless particles
      do j=1,N
      s(j,j)=0._dp
      do k=j+1,N
      s(j,k)=2*dot(p,j,k)
      s(k,j)=s(j,k)
      enddo
      enddo
      return
      end
