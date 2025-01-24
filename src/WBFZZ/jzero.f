!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine jzero(n2,n1,zab,zba,j0)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      complex(dp):: zab(mxpart,4,mxpart),zba(mxpart,4,mxpart),
     & j0(4,2)
c---The one Z-current multiplied by i
c---order of indices Lorentz,quark-line helicity

      integer:: n1,n2,nu
      do nu=1,4
      j0(nu,1)=zab(n2,nu,n1)
      j0(nu,2)=zba(n2,nu,n1)
      enddo
      return
      end
