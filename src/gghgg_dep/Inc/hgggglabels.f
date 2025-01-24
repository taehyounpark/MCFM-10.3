!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      integer, parameter :: dmax=16, cmax=18, bmax=9
c Boxes
      integer, parameter ::
     & d1x2x34=1, d1x4x32=2, d2x1x43=3, d2x3x41=4,
     & d3x2x14=5, d3x4x12=6, d4x1x23=7, d4x3x21=8,
     & d2x34x1=9, d1x23x4=10, d3x41x2=11, d4x12x3=12,
     & d4x1x2=13, d1x2x3=14, d3x4x1=15, d2x3x4=16

c Triangles
      integer, parameter ::
     & c1x234=1, c2x341=2, c3x412=3, c4x123=4,
     & c12x34=5, c23x41=6,
     & c1x23=7, c1x43=8, c2x34=9, c2x14=10,
     & c3x21=11, c3x41=12, c4x12=13, c4x32=14,
     & c1x2=15, c2x3=16, c3x4=17, c4x1=18

c Bubbles
      integer, parameter ::
     & b123=1, b234=2, b341=3, b412=4,
     & b12=5, b23=6, b34=7, b41=8,
     & b1234=9

