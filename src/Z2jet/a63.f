!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function a63(st,j1,j2,j3,j4,j5,j6,za,zb)
      implicit none
c     implementation of last parts of
c     Eqs. (2.15) and (2.16 ) of Nucl Phys. 513, 3 (1998)
      include 'types.f'
      complex(dp):: a63
      include 'mxpart.f'
      include 'zprods_decl.f'
      complex(dp):: a6ax
      integer:: j1,j2,j3,j4,j5,j6
      character(len=2):: st

      if (st == 'pp') then
      a63=-a6ax(j1,j4,j2,j3,j5,j6,za,zb)
      elseif (st == 'pm') then
      a63=+a6ax(j1,j4,j3,j2,j5,j6,za,zb)
      else
      write(6,*) 'Unimplemented st in a63',st
      stop
      endif
      return
      end

