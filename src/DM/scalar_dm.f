!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine scalar_dm(i1,i2,za,zb,amp)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'dm_params.f'
      include 'zprods_decl.f'
      integer:: i1,i2
      complex(dp):: amp(2,2)
      real(dp):: Bp,Beta,sab

      amp(1,2)=czip
      amp(2,1)=czip

      sab=Dble(za(i1,i2)*zb(i2,i1))
c      write(6,*) 'in scalar dm ',sab
c      pause
      Beta=one-4d0*xmass**2/sab
      Beta=sqrt(beta)
      Bp=one+Beta

      amp(1,1)=-za(i1,i2)*Bp+xmass**2/zb(i2,i1)/Bp
      amp(2,2)=zb(i2,i1)*Bp-xmass**2/za(i1,i2)/Bp

      return
      end
