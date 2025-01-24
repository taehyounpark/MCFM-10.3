!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      use scpmpmC12x34m0_unsym_generic
      implicit none
      include 'Inc/zprods_decl.f'
      integer p1,p2,p3,p4
      complex(dp)::scpmpmC12x34m0_res
      scpmpmC12x34m0_res=
     &  scpmpmC12x34m0_unsym(p1,p2,p3,p4,za,zb)
     & +scpmpmC12x34m0_unsym(p3,p4,p1,p2,za,zb)
     & +scpmpmC12x34m0_unsym(p2,p1,p4,p3,zb,za)
     & +scpmpmC12x34m0_unsym(p4,p3,p2,p1,zb,za)
      return


