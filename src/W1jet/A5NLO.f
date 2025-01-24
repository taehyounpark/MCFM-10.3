!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine A5NLO(j1,j2,j3,j4,j5,za,zb,includeaxial,A5LOm,A5NLOm,A5ax)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      integer:: j1,j2,j3,j4,j5
      logical includeaxial
      complex(dp):: A51,A52,A53,A5NLOm,A5LOm,A5ax,zab2
      real(dp):: s125,mt2

c As originally written, the functions A51, A52 correspond to
c 0 --> q_R(1)+qb_L(3)+g_R(2)+ebar_L(4)+e_R(5)
c with all RH couplings
c However we want it in our
c standard form
c       0--> qb_R(1)+q_L(2)++e_L(3)+ebar_R(4)+g_L(5)
c with all LH couplings

c so we have made the changes

c                    'q+g+qb-'   (A51)
c                   (1 ---> 2)
c                   (2 ---> 5)
c                   (3 ---> 1)
c                   (4 ---> 4)
c                   (5 ---> 3)

c                    'q+qb-g+'   (A52)
c                   (1 ---> 2)
c                   (2 ---> 1)
c                   (3 ---> 5)
c                   (4 ---> 4)
c                   (5 ---> 3)

c  and also exchanged za and zb.


c--- corresponds to (1V.1) times minus i, with the (A51) change
      zab2(j1,j2,j3,j4)=za(j1,j2)*zb(j2,j4)+za(j1,j3)*zb(j3,j4)

      A5ax=czip

c      A5LOm=-zb(j1,j4)**2/(zb(j2,j5)*zb(j5,j1)*zb(j4,j3))
      s125=s(j1,j2)+s(j1,j5)+s(j2,j5)
      A5LOm=-zb(j1,j4)*zab2(j3,j2,j5,j1)/(zb(j2,j5)*zb(j5,j1)*s125)
      A5NLOm=A51(j2,j5,j1,j4,j3,zb,za)+A52(j2,j1,j5,j4,j3,zb,za)/xnsq

      if (includeaxial) then
        mt2=mt**2
        A5ax=A53(j2,j1,j5,j4,j3,zb,za,mt2)/xn
      endif

      return
      end
