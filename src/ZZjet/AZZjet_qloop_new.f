!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine AZZjet_qloop_new(helname,swapVV,j1,j2,j3,j4,j5,j6,j7,za,zb,mq,scintsm0,scintsmq,A7q)
c---  Note: array of scalar integrals (scints) is passed in
c---  computed in new notation where basic amplitude is q1- qb2+ l3- lb4+ l5- lb6+ g7+
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'WWjetlabels.f'
      include 'KCdef.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'verbose.f'
      integer j1,j2,j3,j4,j5,j6,j7,icontrib
      real(dp):: msq,msq2,s12,s34,s56,s127,s347,s567,mq
      complex(dp):: A7q(0:2),rat,scintsm0(4,60,-2:0),scintsmq(4,60,-2:0)
      complex(dp)::
     & dd12x34x56,dd12x34x56m0L2,dd12x34x56m0LR,
     & dd12x34x56m2L2,dd12x34x56m2LR,
     & dd12x34x56m4L2,dd12x34x56m4LR,
     & dd12x56x34,dd12x56x34m0L2,dd12x56x34m0LR,
     & dd12x56x34m2L2,dd12x56x34m2LR,
     & dd12x56x34m4L2,dd12x56x34m4LR,
     & dd56x12x34,dd56x12x34m0L2,dd56x12x34m0LR,
     & dd56x12x34m2L2,dd56x12x34m4L2,
     & dd56x12x34m2LR,dd56x12x34m4LR,
     & cc7x12x3456,cc7x34x1256,cc7x56x1234,
     & cc7x12x3456m0,cc7x34x1256m0,cc7x56x1234m0,
     & cc7x12x3456m2,cc7x34x1256m2,cc7x56x1234m2,
     & cc34x56,cc12x56,cc12x34,
     & cc34x56m0,cc12x56m0,cc12x34m0,
     & cc34x56m2,cc12x56m2,cc12x34m2,
     & cc7x12x3456LRm2,cc7x34x1256LRm2,cc7x56x1234LRm2,
     & cc34x56LRm2,cc12x56LRm2,cc12x34LRm2,
     & bb12,bb34,bb56,bb127,bb347,bb567
      complex(dp)::
     & qloop_b12,qloop_b127,
     & qloop_c34x56m2,qloop_ct34x56m2,
     & qloop_c12x56m2,qloop_c7x12m2,
     & qloop_c7x34m2,qloop_ct7x34m2,
     & qloop_ct7x12m2,qloop_ct12x56m2,
     & qloop_c34x56m0,qloop_dt56x12x34m2,
     & qloop_d12x34x56m0,
     & qloop_dt12x34x56m2,qloop_d56x12x34m2,qloop_d12x34x56m2,
     & qloop_dt12x34x56m4,qloop_d12x34x56m4
      character(4):: helname
      logical:: swapVV


c New implementation of boxes
      dd12x34x56m0L2=qloop_d12x34x56m0(j1,j2,j3,j4,j5,j6,j7,za,zb)
      dd56x12x34m0L2=qloop_d12x34x56m0(j5,j6,j1,j2,j3,j4,j7,za,zb)
      dd12x56x34m0L2=qloop_d12x34x56m0(j1,j2,j5,j6,j3,j4,j7,za,zb)

      dd12x34x56m2L2=qloop_d12x34x56m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
      dd12x56x34m2L2=qloop_d12x34x56m2(j1,j2,j5,j6,j3,j4,j7,za,zb)

      dd12x56x34m2LR=qloop_dt12x34x56m2(j1,j2,j5,j6,j3,j4,j7,za,zb)
      dd56x12x34m4LR=qloop_dt12x34x56m4(j1,j2,j3,j4,j5,j6,j7,za,zb)
      dd12x56x34m4LR=dd56x12x34m4LR

      dd56x12x34m4L2=qloop_d12x34x56m4(j1,j2,j3,j4,j5,j6,j7,za,zb)
      dd12x56x34m4L2=dd56x12x34m4L2
      dd12x34x56m4L2=dd56x12x34m4L2

      dd56x12x34m2L2=qloop_d56x12x34m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
      dd12x56x34m0LR=czip
      dd12x34x56m0LR=czip
      dd56x12x34m0LR=czip

      dd56x12x34m2LR=qloop_dt56x12x34m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
      dd12x34x56m2LR=qloop_dt12x34x56m2(j1,j2,j3,j4,j5,j6,j7,za,zb)

      dd12x34x56m4LR=dd56x12x34m4LR

      s12=s(j1,j2)
      s34=s(j3,j4)
      s56=s(j5,j6)
      s127=s12+s(j1,j7)+s(j2,j7)
      s347=s34+s(j3,j7)+s(j4,j7)
      s567=s56+s(j5,j7)+s(j6,j7)
c Infrared relations for m^0 parts of triangles
       cc7x12x3456m0=(dd12x56x34m0L2/(s127*s347-s12*s34)
     &              +dd12x34x56m0L2/(s127*s567-s12*s56))*(s12-s127)
       cc7x34x1256m0=(dd12x56x34m0L2/(s127*s347-s12*s34)
     &              +dd56x12x34m0L2 /(s347*s567-s34*s56))*(s34-s347)
       cc7x56x1234m0=(dd12x34x56m0L2/(s127*s567-s12*s56)
     &              +dd56x12x34m0L2 /(s347*s567-s34*s56))*(s56-s567)

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

      bb127=qloop_b127(j1,j2,j3,j4,j5,j6,j7,za,zb)
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
     & -dd12x34x56m4L2)

      cc7x12x3456LRm2=qloop_ct7x12m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
      cc7x34x1256LRm2=qloop_ct7x34m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
      cc7x56x1234LRm2=qloop_ct7x34m2(j1,j2,j5,j6,j3,j4,j7,za,zb)
      cc34x56LRm2=qloop_ct34x56m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
      cc12x56LRm2=qloop_ct12x56m2(j1,j2,j3,j4,j5,j6,j7,za,zb)
      cc12x34LRm2=qloop_ct12x56m2(j1,j2,j5,j6,j3,j4,j7,za,zb)

c loop over masses
      do icontrib=0,2
      if (icontrib == 0) then
        msq = 0d0
        msq2 = 0d0
        scints(:,:,:)=scintsm0(:,:,:)
        if (verbose) call parser('runZZ','QL_mb0',helname,KCD,KCC,KCB,rattot,tot)
      else
        msq = mq**2
        msq2 = mq**4
        scints(:,:,:)=scintsmq(:,:,:)
        if (icontrib == 1) then
          if (verbose) call parser('runZZ','QL_mbvLsq',helname,KCD,KCC,KCB,rattot,tot)
        else
          if (verbose) call parser('runZZ','QL_mbvLvR',helname,KCD,KCC,KCB,rattot,tot)
        endif
      endif

      if (icontrib < 2) then
        dd12x56x34=dd12x56x34m0L2+msq*dd12x56x34m2L2+msq2*dd12x56x34m4L2
        dd12x34x56=dd12x34x56m0L2+msq*dd12x34x56m2L2+msq2*dd12x34x56m4L2
        dd56x12x34=dd56x12x34m0L2+msq*dd56x12x34m2L2+msq2*dd56x12x34m4L2

        cc34x56=cc34x56m0+msq*cc34x56m2
        cc12x56=cc12x56m0+msq*cc12x56m2
        cc12x34=cc12x34m0+msq*cc12x34m2

c In massless case, do not include triangles with poles (net contribution zero)
        if (icontrib == 0) then
          cc7x12x3456=czip
          cc7x34x1256=czip
          cc7x56x1234=czip
        else
          cc7x12x3456=cc7x12x3456m0+msq*cc7x12x3456m2
          cc7x34x1256=cc7x34x1256m0+msq*cc7x34x1256m2
          cc7x56x1234=cc7x56x1234m0+msq*cc7x56x1234m2
        endif
      else
        dd12x56x34=dd12x56x34m0LR+msq*dd12x56x34m2LR+msq2*dd12x56x34m4LR
        dd12x34x56=dd12x34x56m0LR+msq*dd12x34x56m2LR+msq2*dd12x34x56m4LR
        dd56x12x34=dd56x12x34m0LR+msq*dd56x12x34m2LR+msq2*dd56x12x34m4LR
        cc7x12x3456=msq*cc7x12x3456LRm2
        cc7x34x1256=msq*cc7x34x1256LRm2
        cc7x56x1234=msq*cc7x56x1234LRm2
        cc34x56=msq*cc34x56LRm2
        cc12x56=msq*cc12x56LRm2
        cc12x34=msq*cc12x34LRm2
        bb12=czip
        bb34=czip
        bb56=czip
        bb127=czip
        bb347=czip
        bb567=czip
        rat=czip
      endif

c--- translation from notation of WWjetcomputescalars.f to here
c--- (for positive helicity gluon)
c---
c---  box1 -> d12x56x34sl
c---  box2 -> d34x12x56slb
c---  box3 -> d34x12x56sla
c---
c---  tri123 -> c56x34
c---  tri145 -> c12x56sl
c---  tri167 -> c12x34sl
c---
c---  bub12 -> b12sl
c---  bub34 -> b34
c---  bub56 -> b56
c---  bub127 -> b127
c---  bub347 -> b347sl
c---  bub567 -> b567sl

      if (swapVV) then
        A7q(icontrib)=
     &  +dd12x56x34*scints(4,d34x12x56slb,0)
     &  +dd12x34x56*scints(4,d12x56x34sl,0)
     &  +dd56x12x34*scints(4,d34x12x56sla,0)
     &  +cc34x56*scints(3,c56x34,0)
     &  +cc12x56*scints(3,c12x34sl,0)
     &  +cc12x34*scints(3,c12x56sl,0)
     &  +cc7x12x3456*scints(3,c7x12sl,0)
     &  +cc7x34x1256*scints(3,c7x56sl,0)
     &  +cc7x56x1234*scints(3,c7x34sl,0)
     &  +bb12*scints(2,b12sl,0)
     &  +bb34*scints(2,b56,0)
     &  +bb56*scints(2,b34,0)
     &  +bb127*scints(2,b127,0)
     &  +bb347*scints(2,b567sl,0)
     &  +bb567*scints(2,b347sl,0)
     &  +rat
      else
        A7q(icontrib)=
     &  +dd12x56x34*scints(4,d7x12x56sl,0)
     &  +dd12x34x56*scints(4,d7x12x34sl,0)
     &  +dd56x12x34*scints(4,d7x34x12sl,0)
     &  +cc34x56*scints(3,c56x34,0)
     &  +cc12x56*scints(3,c12x56sl,0)
     &  +cc12x34*scints(3,c12x34sl,0)
     &  +cc7x12x3456*scints(3,c7x12sl,0)
     &  +cc7x34x1256*scints(3,c7x34sl,0)
     &  +cc7x56x1234*scints(3,c7x56sl,0)
     &  +bb12*(scints(2,b12sl,0)-scints(2,b127,0))
     &  +bb34*(scints(2,b34,0)-scints(2,b127,0))
     &  +bb56*(scints(2,b56,0)-scints(2,b127,0))
c     &  +bb127*(scints(2,b127,0)-scints(2,b127,0))
     &  +bb347*(scints(2,b347sl,0)-scints(2,b127,0))
     &  +bb567*(scints(2,b567sl,0)-scints(2,b127,0))
     &  +rat
c Use the fact that the sum of bubbles is zero
c     &  +bb12*scints(2,b12sl,0)
c     &  +bb34*scints(2,b34,0)
c     &  +bb56*scints(2,b56,0)
c     &  +bb127*scints(2,b127,0)
c     &  +bb347*scints(2,b347sl,0)
c     &  +bb567*scints(2,b567sl,0)
      endif

      if (verbose) then
      write(6,*)
      write(6,*) 'QL:Comparison to KCheck, msq =',msq
      if     (swapVV .eqv. .false.) then
        write(6,*) 'KCD(d23x1x45q)',dd12x56x34/KCD(d23x1x45q)
        write(6,*) 'KCD(d23x1x67q)',dd12x34x56/KCD(d23x1x67q)
        write(6,*) 'KCD(d45x1x67q)',dd56x12x34/KCD(d45x1x67q)
        write(6,*) 'KCC(c1x23q)   ',cc7x12x3456/KCC(c1x23q)
        write(6,*) 'KCC(c1x45q)   ',cc7x34x1256/KCC(c1x45q)
        write(6,*) 'KCC(c1x67q)   ',cc7x56x1234/KCC(c1x67q)
        write(6,*) 'KCC(c45x67q)',cc34x56/KCC(c45x67q)
        write(6,*) 'KCC(c23x67q)',cc12x56/KCC(c23x67q)
        write(6,*) 'KCC(c23x45q)',cc12x34/KCC(c23x45q)
        write(6,*) 'KCB(b23q)     ',bb12/KCB(b23q)
        write(6,*) 'KCB(b123q)    ',bb127/KCB(b123q)
        write(6,*) 'KCB(b45q)     ',bb34/KCB(b45q)
        write(6,*) 'KCB(b145q)    ',bb347/KCB(b145q)
        write(6,*) 'KCB(b67q)     ',bb56/KCB(b67q)
        write(6,*) 'KCB(b167q)    ',bb567/KCB(b167q)
        write(6,*) 'rattot        ',rat/rattot
        write(6,*) 'TOTAL         ',A7q(icontrib)/tot
      elseif (swapVV .eqv. .true.) then
        write(6,*) 'KCD(d23x1x45q)',dd12x56x34/KCD(d23x1x67q)
        write(6,*) 'KCD(d23x1x67q)',dd12x34x56/KCD(d23x1x45q)
        write(6,*) 'KCD(d45x1x67q)',dd56x12x34/KCD(d45x1x67q)
        write(6,*) 'KCC(c1x23q)   ',cc7x12x3456/KCC(c1x23q)
        write(6,*) 'KCC(c1x45q)   ',cc7x34x1256/KCC(c1x67q)
        write(6,*) 'KCC(c1x67q)   ',cc7x56x1234/KCC(c1x45q)
        write(6,*) 'KCC(c45x67q)',cc34x56/KCC(c45x67q)
        write(6,*) 'KCC(c23x67q)',cc12x56/KCC(c23x45q)
        write(6,*) 'KCC(c23x45q)',cc12x34/KCC(c23x67q)
        write(6,*) 'KCB(b23q)     ',bb12/KCB(b23q)
        write(6,*) 'KCB(b123q)    ',bb127/KCB(b123q)
        write(6,*) 'KCB(b45q)     ',bb34/KCB(b67q)
        write(6,*) 'KCB(b145q)    ',bb347/KCB(b167q)
        write(6,*) 'KCB(b67q)     ',bb56/KCB(b45q)
        write(6,*) 'KCB(b167q)    ',bb567/KCB(b145q)
        write(6,*) 'rattot        ',rat/rattot
        write(6,*) 'TOTAL         ',A7q(icontrib)/tot
      endif
      endif

      enddo

      return
      end

