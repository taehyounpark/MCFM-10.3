!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine ZZjetqloopcomputescalars(mq,p1,p2,p3,p4,p5,p6,p7,scints)
c--- routine to compute all scalar integrals used in the ZZ+jet calculation
c--- of the qloop contribution
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'WWjetlabels.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'scalarselect.f'
      include 'first.f'
      logical:: writescalars
      integer p1,p2,p3,p4,p5,p6,p7,iep
      real(dp):: t,s12,s56,s34,s127,s567,s347,mq,msq

c--- statement function
      t(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)

c--- initialize QCDLoop, if necessary
      if (first) then
c        call qlinit
        first=.false.
        writescalars=.false.
        if (writescalars) then
          open(unit=67,file='scalars.out',status='unknown')
        endif
      endif

      s12=s(p1,p2)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s127=t(p1,p2,p7)
      s347=t(p3,p4,p7)
      s567=t(p5,p6,p7)

      msq=mq**2

c---  NB: only need finite pieces, poles will be handled separately
      iep=0

      scints(4,d7x34x12sl,iep)
     & =loopI4(s56,s12,s34,zip,s347,s567,msq,msq,msq,msq,musq,iep)

      scints(4,d7x12x56sl,iep)
     & =loopI4(s56,s12,zip,s34,s347,s127,msq,msq,msq,msq,musq,iep)

      scints(4,d7x12x34sl,iep)
     & =loopI4(zip,s12,s34,s56,s127,s567,msq,msq,msq,msq,musq,iep)

      scints(3,c56x34,iep)=loopI3(s34,s56,s127,msq,msq,msq,musq,iep)
      scints(3,c12x56sl,iep)=loopI3(s56,s347,s12,msq,msq,msq,musq,iep)
      scints(3,c12x34sl,iep)=loopI3(s34,s567,s12,msq,msq,msq,musq,iep)

      scints(3,c7x12sl,iep)=loopI3(zip,s12,s127,msq,msq,msq,musq,iep)
      scints(3,c7x34sl,iep)=loopI3(zip,s34,s347,msq,msq,msq,musq,iep)
      scints(3,c7x56sl,iep)=loopI3(zip,s56,s567,msq,msq,msq,musq,iep)

      scints(2,b34,iep)=loopI2(s34,msq,msq,musq,iep)
      scints(2,b56,iep)=loopI2(s56,msq,msq,musq,iep)
      scints(2,b127,iep)=loopI2(s127,msq,msq,musq,iep)

      scints(2,b347sl,iep)=loopI2(s347,msq,msq,musq,iep)
      scints(2,b567sl,iep)=loopI2(s567,msq,msq,musq,iep)
      scints(2,b12sl,iep)=loopI2(s12,msq,msq,musq,iep)

      return
      end

