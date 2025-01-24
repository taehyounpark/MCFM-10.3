!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AWWjet_qloop(helname,swapVV,p,j1,j2,j3,j4,j5,j6,j7,za,zb,scints,A7q)
c--- Note: array of scalar integrals (scints) is passed in
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'WWjetlabels.f'
      include 'zprods_decl.f'
c      include 'KCdef.f'
c      include 'verbose.f'
      integer:: j1,j2,j3,j4,j5,j6,j7
      real(dp):: p(mxpart,4)
      complex(dp):: A7q,
     & qloop_box,qloop_b12,qloop_b127,
     & box1,box2,box3,tri56x34,tri12x56,tri12x34,
     & bub12,bub34,bub56,bub127,bub347,bub567,rat,qloop_c34x56m0
      logical:: swapVV
      character(len=2):: helname

c      call parser('runWW','QL',helname,KCD,KCC,KCB,rattot,tot)
c      include 'runQL.f'

c--- boxes
      box1=qloop_box(j5,j6,j1,j2,j3,j4,j7,za,zb)
      box2=qloop_box(j4,j3,j1,j2,j6,j5,j7,za,zb)
      box3=qloop_box(j2,j1,j6,j5,j3,j4,j7,za,zb)

c      if (verbose) then
c      write(6,*)
c      write(6,*) 'QL:Comparison to KCheck'
c      if     (swapVV .eqv. .false.) then
c        write(6,*) 'KCD(d23x1x45q)',box1/KCD(d23x1x45q)
c        write(6,*) 'KCD(d23x1x67q)',box2/KCD(d23x1x67q)
c        write(6,*) 'KCD(d45x1x67q)',box3/KCD(d45x1x67q)
c      elseif (swapVV .eqv. .true.) then
c        write(6,*) 'KCD(d23x1x45q)',box1/KCD(d23x1x67q)
c        write(6,*) 'KCD(d23x1x67q)',box2/KCD(d23x1x45q)
c        write(6,*) 'KCD(d45x1x67q)',box3/KCD(d45x1x67q)
c      endif
c      endif

c--- triangles
      tri56x34=qloop_c34x56m0(j1,j2,j6,j5,j3,j4,j7,za,zb)
      tri12x34=qloop_c34x56m0(j6,j5,j1,j2,j3,j4,j7,za,zb)
      tri12x56=qloop_c34x56m0(j3,j4,j6,j5,j1,j2,j7,za,zb)
c      write(6,*) 'tri56x34',tri56x34/
c     &      (qloop_tri(j1,j2,j6,j5,j3,j4,j7,za,zb)
c     &      +qloop_tri(j1,j2,j3,j4,j6,j5,j7,za,zb))
c      write(6,*) 'tri12x34',tri12x34/
c     &      (qloop_tri(j6,j5,j1,j2,j3,j4,j7,za,zb)
c     &      +qloop_tri(j6,j5,j3,j4,j1,j2,j7,za,zb))
c      write(6,*) 'tri12x56',tri12x56/
c     &      (qloop_tri(j3,j4,j6,j5,j1,j2,j7,za,zb)
c     &      +qloop_tri(j3,j4,j1,j2,j6,j5,j7,za,zb))
c      pause

c      if (verbose) then
c      write(6,*)
c      if     (swapVV .eqv. .false.) then
c        write(6,*) 'KCC(d45x1x67q)',tri56x34/KCC(c45x67q)
c        write(6,*) 'KCC(d23x1x67q)',tri12x56/KCC(c23x67q)
c        write(6,*) 'KCC(d45x1x67q)',tri12x34/KCC(c23x45q)
c      elseif (swapVV .eqv. .true.) then
c        write(6,*) 'KCC(d45x1x67q)',tri56x34/KCC(c45x67q)
c        write(6,*) 'KCC(d23x1x67q)',tri12x56/KCC(c23x45q)
c        write(6,*) 'KCC(d45x1x67q)',tri12x34/KCC(c23x67q)
c      endif
c      endif

c--- bubbles

c      bub12=qloop_bub12(j1,j2,j6,j5,j3,j4,j7,za,zb)
c     &     +qloop_bub12(j1,j2,j3,j4,j6,j5,j7,za,zb)
c      bub56=qloop_bub12(j6,j5,j1,j2,j3,j4,j7,za,zb)
c     &     +qloop_bub12(j6,j5,j3,j4,j1,j2,j7,za,zb)
c      bub34=qloop_bub12(j3,j4,j6,j5,j1,j2,j7,za,zb)
c     &     +qloop_bub12(j3,j4,j1,j2,j6,j5,j7,za,zb)

      bub12=qloop_b12(j1,j2,j3,j4,j6,j5,j7,za,zb)
      bub56=qloop_b12(j6,j5,j3,j4,j1,j2,j7,za,zb)
      bub34=qloop_b12(j3,j4,j1,j2,j6,j5,j7,za,zb)

      bub127=qloop_b127(j1,j2,j3,j4,j6,j5,j7,za,zb)
      bub567=qloop_b127(j6,j5,j3,j4,j1,j2,j7,za,zb)
      bub347=qloop_b127(j3,j4,j1,j2,j6,j5,j7,za,zb)

c      if (verbose) then
c      write(6,*)
c      if     (swapVV .eqv. .false.) then
c        write(6,*) 'KCB(b23q)',bub12/KCB(b23q)
c        write(6,*) 'KCB(b123q)',bub127/KCB(b123q)
c        write(6,*) 'KCB(b45q)',bub34/KCB(b45q)
c        write(6,*) 'KCB(b145q)',bub347/KCB(b145q)
c        write(6,*) 'KCB(b67q)',bub56/KCB(b67q)
c        write(6,*) 'KCB(b167q)',bub567/KCB(b167q)
c      elseif (swapVV .eqv. .true.) then
c        write(6,*) 'KCB(b23q)',bub12/KCB(b23q)
c        write(6,*) 'KCB(b123q)',bub127/KCB(b123q)
c        write(6,*) 'KCB(b45q)',bub34/KCB(b67q)
c        write(6,*) 'KCB(b145q)',bub347/KCB(b167q)
c        write(6,*) 'KCB(b67q)',bub56/KCB(b45q)
c        write(6,*) 'KCB(b167q)',bub567/KCB(b145q)
c      endif
c      endif

c--- rational piece
      call qloop_rat(j7,j1,j2,j6,j5,j3,j4,za,zb,rat)

c      if (verbose) then
c      write(6,*)
c      write(6,*) 'rational',rat/rattot
c      endif

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
     &  -box1*scints(4,d34x12x56slb,0)
     &  -box2*scints(4,d12x56x34sl,0)
     &  -box3*scints(4,d34x12x56sla,0)
     &  -tri56x34*scints(3,c56x34,0)
     &  -tri12x56*scints(3,c12x34sl,0)
     &  -tri12x34*scints(3,c12x56sl,0)
     &  -bub12*scints(2,b12sl,0)
     &  -bub34*scints(2,b56,0)
     &  -bub56*scints(2,b34,0)
     &  -bub127*scints(2,b127,0)
     &  -bub347*scints(2,b567sl,0)
     &  -bub567*scints(2,b347sl,0)
     &  -rat
      else
        A7q=
     &  +box1*scints(4,d7x12x56sl,0)
     &  +box2*scints(4,d7x12x34sl,0)
     &  +box3*scints(4,d7x34x12sl,0)
     &  +tri56x34*scints(3,c56x34,0)
     &  +tri12x56*scints(3,c12x56sl,0)
     &  +tri12x34*scints(3,c12x34sl,0)
     &  +bub12*scints(2,b12sl,0)
     &  +bub34*scints(2,b34,0)
     &  +bub56*scints(2,b56,0)
     &  +bub127*scints(2,b127,0)
     &  +bub347*scints(2,b347sl,0)
     &  +bub567*scints(2,b567sl,0)
     &  +rat
      endif

c      if (verbose) then
c      write(6,*)
c      write(6,*) 'total',A7q/tot
c      endif


      return
      end

