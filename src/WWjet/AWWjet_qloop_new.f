!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWWjet_qloop_new(helname,swapVV,j1,j2,j3,j4,j5,j6,j7,za,zb,scints,A7q)
c--- Note: array of scalar integrals (scints) is passed in
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'WWjetlabels.f'
      include 'zprods_decl.f'
c      include 'KCdef.f'
c      include 'verbose.f'
      integer j1,j2,j3,j4,j5,j6,j7
      complex(dp)::
     & dd12x34x56m0,
     & dd12x34x56m4,
     & dd12x56x34m0,
c     & dd12x56x34m4,
     & dd56x12x34m0,
     & dd56x12x34m4,
     & cc7x12x3456m2,cc7x34x1256m2,cc7x56x1234m2,
     & cc34x56m0,cc12x56m0,cc12x34m0,
     & cc34x56m2,cc12x56m2,cc12x34m2,
     & bb12,bb34,bb56,bb347,bb567!,bb127
      complex(dp):: A7q,rat,
     & qloop_b12,qloop_b127,
     & qloop_c34x56m2,
     & qloop_c12x56m2,
     & qloop_c7x12m2,
     & qloop_c7x34m2,
     & qloop_c34x56m0,
     & qloop_d12x34x56m0,
     & qloop_d12x34x56m4
      logical:: swapVV
      character(len=2):: helname

c      call parser('runWW','QL',helname,KCD,KCC,KCB,rattot,tot)
c      include 'runQL.f'

c New implementation of boxes
      dd12x34x56m0=qloop_d12x34x56m0(j1,j2,j3,j4,j5,j6,j7,za,zb)
      dd56x12x34m0=qloop_d12x34x56m0(j5,j6,j1,j2,j3,j4,j7,za,zb)
      dd12x56x34m0=qloop_d12x34x56m0(j1,j2,j5,j6,j3,j4,j7,za,zb)

      dd56x12x34m4=qloop_d12x34x56m4(j1,j2,j3,j4,j5,j6,j7,za,zb)
c      dd12x56x34m4=dd56x12x34m4
      dd12x34x56m4=dd56x12x34m4

c 3-mass triangles
       cc34x56m0=qloop_c34x56m0(j1,j2,j3,j4,j5,j6,j7,za,zb)
       cc12x56m0=qloop_c34x56m0(j3,j4,j1,j2,j5,j6,j7,za,zb)
       cc12x34m0=qloop_c34x56m0(j5,j6,j1,j2,j3,j4,j7,za,zb)

c Massive parts of triangles
       cc7x12x3456m2=qloop_c7x12m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
       cc7x34x1256m2=+qloop_c7x34m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
       cc7x56x1234m2=+qloop_c7x34m2(j1,j2,j5,j6,j3,j4,j7,za,zb)
       cc34x56m2=+qloop_c34x56m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
       cc12x56m2=+qloop_c12x56m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
       cc12x34m2=+qloop_c12x56m2(j1,j2,j5,j6,j3,j4,j7,za,zb)

c Bubbles (independent of mass)
c      bb127=qloop_b127(j1,j2,j3,j4,j5,j6,j7,za,zb)
      bb567=qloop_b127(j5,j6,j3,j4,j1,j2,j7,za,zb)
      bb347=qloop_b127(j3,j4,j1,j2,j5,j6,j7,za,zb)

      bb12=qloop_b12(j1,j2,j3,j4,j5,j6,j7,za,zb)
      bb56=qloop_b12(j5,j6,j3,j4,j1,j2,j7,za,zb)
      bb34=qloop_b12(j3,j4,j1,j2,j5,j6,j7,za,zb)

c Rational term (independent of mass)
c--- each triangle contributes (+1/2), each box (-1/6) and there
c--- are three boxes, each with the same rational contribution
      rat=0.5d0*(
     &  cc7x12x3456m2+cc7x34x1256m2+cc7x56x1234m2
     & +cc34x56m2+cc12x56m2+cc12x34m2
     & -dd12x34x56m4)

c--- translation from notation of WWjetcomputescalars.f to here
c--- (for positive helicity gluon)
c---
c---  box1 -> d12x56x34sl
c---  box2 -> d34x12x56slb
c---  box3 -> d34x12x56sla
c---
c---  tri56x34 -> c56x34
c---  tri12x56 -> c12x56sl
c---  tri12x34 -> c12x34sl
c---
c---  bub12 -> b12sl
c---  bub34 -> b34
c---  bub56 -> b56
c---  bub127 -> b127
c---  bub347 -> b347sl
c---  bub567 -> b567sl

      if (swapVV) then
        A7q=
     &  -dd12x56x34m0*scints(4,d7x12x34sl,0)
     &  -dd12x34x56m0*scints(4,d7x12x56sl,0)
     &  -dd56x12x34m0*scints(4,d7x34x12sl,0)
     &  -cc34x56m0*scints(3,c56x34,0)
     &  -cc12x56m0*scints(3,c12x34sl,0)
     &  -cc12x34m0*scints(3,c12x56sl,0)
     &  -bb12*(scints(2,b12sl,0)-scints(2,b127,0))
     &  -bb34*(scints(2,b56,0)-scints(2,b127,0))
     &  -bb56*(scints(2,b34,0)-scints(2,b127,0))
     &  -bb347*(scints(2,b567sl,0)-scints(2,b127,0))
     &  -bb567*(scints(2,b347sl,0)-scints(2,b127,0))
c     &  -bb127*(scints(2,b127,0)-scints(2,b127,0))
     &  -rat
c     &  -bb12*scints(2,b12sl,0)
c     &  -bb34*scints(2,b56,0)
c     &  -bb56*scints(2,b34,0)
c     &  -bb127*scints(2,b127,0)
c     &  -bb347*scints(2,b567sl,0)
c     &  -bb567*scints(2,b347sl,0)
c Use the fact that the sum of bubbles is zero
      else
        A7q=
     &  +dd12x56x34m0*scints(4,d7x12x56sl,0)
     &  +dd12x34x56m0*scints(4,d7x12x34sl,0)
     &  +dd56x12x34m0*scints(4,d7x34x12sl,0)
     &  +cc34x56m0*scints(3,c56x34,0)
     &  +cc12x56m0*scints(3,c12x56sl,0)
     &  +cc12x34m0*scints(3,c12x34sl,0)
     &  +bb12*(scints(2,b12sl,0)-scints(2,b127,0))
     &  +bb34*(scints(2,b34,0)-scints(2,b127,0))
     &  +bb56*(scints(2,b56,0)-scints(2,b127,0))
     &  +bb347*(scints(2,b347sl,0)-scints(2,b127,0))
     &  +bb567*(scints(2,b567sl,0)-scints(2,b127,0))
c     &  +bb127*(scints(2,b127,0)-scints(2,b127,0))
     &  +rat
c     &  +bb12*scints(2,b12sl,0)
c     &  +bb34*scints(2,b34,0)
c     &  +bb56*scints(2,b56,0)
c     &  +bb127*scints(2,b127,0)
c     &  +bb347*scints(2,b347sl,0)
c     &  +bb567*scints(2,b567sl,0)
c Use the fact that the sum of bubbles is zero
      endif

c      if (verbose) then
c      write(6,*)
c      write(6,*) 'QL:Comparison to KCheck'
c      if     (swapVV .eqv. .false.) then
c        write(6,*) 'KCD(d23x1x45q)',dd12x56x34/KCD(d23x1x45q)
c        write(6,*) 'KCD(d23x1x67q)',dd12x34x56/KCD(d23x1x67q)
c        write(6,*) 'KCD(d45x1x67q)',dd56x12x34/KCD(d45x1x67q)
c        write(6,*) 'KCC(c1x23q)   ',cc7x12x3456/KCC(c1x23q)
c        write(6,*) 'KCC(c1x45q)   ',cc7x34x1256/KCC(c1x45q)
c        write(6,*) 'KCC(c1x67q)   ',cc7x56x1234/KCC(c1x67q)
c        write(6,*) 'KCC(c45x67q)',cc34x56/KCC(c45x67q)
c        write(6,*) 'KCC(c23x67q)',cc12x56/KCC(c23x67q)
c        write(6,*) 'KCC(c23x45q)',cc12x34/KCC(c23x45q)
c        write(6,*) 'KCB(b23q)     ',bb12/KCB(b23q)
c        write(6,*) 'KCB(b123q)    ',bb127/KCB(b123q)
c        write(6,*) 'KCB(b45q)     ',bb34/KCB(b45q)
c        write(6,*) 'KCB(b145q)    ',bb347/KCB(b145q)
c        write(6,*) 'KCB(b67q)     ',bb56/KCB(b67q)
c        write(6,*) 'KCB(b167q)    ',bb567/KCB(b167q)
c        write(6,*) 'rattot        ',rat/rattot
c        write(6,*) 'TOTAL         ',A7q/tot
c      elseif (swapVV .eqv. .true.) then
c        write(6,*) 'KCD(d23x1x45q)',dd12x56x34/KCD(d23x1x67q)
c        write(6,*) 'KCD(d23x1x67q)',dd12x34x56/KCD(d23x1x45q)
c        write(6,*) 'KCD(d45x1x67q)',dd56x12x34/KCD(d45x1x67q)
c        write(6,*) 'KCC(c1x23q)   ',cc7x12x3456/KCC(c1x23q)
c        write(6,*) 'KCC(c1x45q)   ',cc7x34x1256/KCC(c1x67q)
c        write(6,*) 'KCC(c1x67q)   ',cc7x56x1234/KCC(c1x45q)
c        write(6,*) 'KCC(c45x67q)',cc34x56/KCC(c45x67q)
c        write(6,*) 'KCC(c23x67q)',cc12x56/KCC(c23x45q)
c        write(6,*) 'KCC(c23x45q)',cc12x34/KCC(c23x67q)
c        write(6,*) 'KCB(b23q)     ',bb12/KCB(b23q)
c        write(6,*) 'KCB(b123q)    ',bb127/KCB(b123q)
c        write(6,*) 'KCB(b45q)     ',bb34/KCB(b67q)
c        write(6,*) 'KCB(b145q)    ',bb347/KCB(b167q)
c        write(6,*) 'KCB(b67q)     ',bb56/KCB(b45q)
c        write(6,*) 'KCB(b167q)    ',bb567/KCB(b145q)
c        write(6,*) 'rattot        ',rat/rattot
c        write(6,*) 'TOTAL         ',A7q/tot
c      endif
c      endif

      return
      end

