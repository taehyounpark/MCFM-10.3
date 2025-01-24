!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine WWjetcomputescalars(p1,p2,p3,p4,p5,p6,p7,scints)
c--- routine to compute all scalar integrals used in the WW+jet calculation
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
      logical:: writescalars
      integer:: p1,p2,p3,p4,p5,p6,p7,iep
      real(dp):: t,s12,s17,s27,s56,s34,s127,s134,s256,s567,s347,s156,s234
c      complex(dp):: loopI4,loopI3,loopI2
      complex(dp):: onemassbox,easy,hard,threemassbox

c--- statement function
      t(p1,p2,p3)=s(p1,p2)+s(p1,p3)+s(p2,p3)


      s12=s(p1,p2)
      s17=s(p1,p7)
      s27=s(p2,p7)
      s34=s(p3,p4)
      s56=s(p5,p6)
      s127=t(p1,p2,p7)
      s134=t(p1,p3,p4)
      s256=t(p2,p5,p6)
      s347=t(p3,p4,p7)
      s567=t(p5,p6,p7)
      s234=t(p2,p3,p4)
      s156=t(p1,p5,p6)

c---  NB: only need finite pieces, poles will be handled separately
c      do iep=-2,0
      iep=0
c      do iep=0,0

c--- integrals appearing at leading color
      scints(4,d7x1x34,iep)
     & =hard(zip,zip,s34,s256,s17,s134,zip,zip,zip,zip)
c     & =loopI4(zip,zip,s34,s256,s17,s134,zip,zip,zip,zip,musq,iep)

      scints(4,d27x1x34,iep)
     & =threemassbox(zip,s27,s56,s34,s127,s134,zip,zip,zip,zip)
c     & =loopI4(s56,s34,zip,s27,s127,s134,zip,zip,zip,zip,musq,iep)

      scints(4,d7x2x56,iep)
     & =hard(zip,zip,s56,s134,s27,s256,zip,zip,zip,zip)
c     & =loopI4(zip,zip,s56,s134,s27,s256,zip,zip,zip,zip,musq,iep)

      scints(4,d17x2x56,iep)
     & =threemassbox(zip,s17,s34,s56,s127,s256,zip,zip,zip,zip)
c     & =loopI4(s34,s56,zip,s17,s127,s256,zip,zip,zip,zip,musq,iep)

      scints(4,d1x7x2,iep)
     & =onemassbox(zip,zip,zip,s127,s17,s27,zip,zip,zip,zip)
c     & =loopI4(zip,zip,zip,s127,s17,s27,zip,zip,zip,zip,musq,iep)

      scints(3,c17x34,iep)=loopI3(s34,s17,s256,zip,zip,zip,musq,iep)
      scints(3,c27x56,iep)=loopI3(s56,s27,s134,zip,zip,zip,musq,iep)
      scints(3,c56x34,iep)=loopI3(s34,s56,s127,zip,zip,zip,musq,iep)

      scints(2,b134,iep)=loopI2(s134,zip,zip,musq,iep)
      scints(2,b256,iep)=loopI2(s256,zip,zip,musq,iep)
      scints(2,b34,iep)=loopI2(s34,zip,zip,musq,iep)
      scints(2,b17,iep)=loopI2(s17,zip,zip,musq,iep)
      scints(2,b56,iep)=loopI2(s56,zip,zip,musq,iep)
      scints(2,b127,iep)=loopI2(s127,zip,zip,musq,iep)
      scints(2,b27,iep)=loopI2(s27,zip,zip,musq,iep)

c      return

c---integrals only appearing at subleading color
      scints(4,d1x34x7sl,iep)
c     & =loopI4(zip,s34,zip,s256,s347,s134,zip,zip,zip,zip,musq,iep)
     & =easy(zip,s34,zip,s256,s347,s134,zip,zip,zip,zip)

      scints(4,d1x7x34sl,iep)
c     & =loopI4(s34,zip,zip,s256,s347,s17,zip,zip,zip,zip,musq,iep)
     & =hard(s34,zip,zip,s256,s347,s17,zip,zip,zip,zip)

      scints(4,d2x56x7sl,iep)
c     & =loopI4(zip,s56,zip,s134,s567,s256,zip,zip,zip,zip,musq,iep)
     & =easy(zip,s56,zip,s134,s567,s256,zip,zip,zip,zip)

      scints(4,d2x7x56sl,iep)
c     & =loopI4(s56,zip,zip,s134,s567,s27,zip,zip,zip,zip,musq,iep)
     & =hard(s56,zip,zip,s134,s567,s27,zip,zip,zip,zip)

      scints(4,d17x2x56sl,iep)
c     & =loopI4(s34,s56,zip,s17,s127,s256,zip,zip,zip,zip,musq,iep)
     & =threemassbox(zip,s17,s34,s56,s127,s256,zip,zip,zip,zip)

      scints(4,d1x2x7sl,iep)
c     & =loopI4(zip,zip,zip,s127,s27,s12,zip,zip,zip,zip,musq,iep)
     & =onemassbox(zip,zip,zip,s127,s27,s12,zip,zip,zip,zip)

      scints(4,d2x1x7sl,iep)
c     & =loopI4(zip,zip,zip,s127,s17,s12,zip,zip,zip,zip,musq,iep)
     & =onemassbox(zip,zip,zip,s127,s17,s12,zip,zip,zip,zip)

      scints(4,d1x2x56sl,iep)
c     & =loopI4(s56,zip,zip,s347,s256,s12,zip,zip,zip,zip,musq,iep)
     & =hard(s56,zip,zip,s347,s256,s12,zip,zip,zip,zip)

      scints(4,d7x34x12sl,iep)
c     & =loopI4(s56,s12,s34,zip,s347,s567,zip,zip,zip,zip,musq,iep)
     & =threemassbox(zip,s56,s12,s34,s567,s347,zip,zip,zip,zip)

      scints(4,d2x1x34sl,iep)
c     & =loopI4(zip,zip,s34,s567,s12,s134,zip,zip,zip,zip,musq,iep)
     & =hard(zip,zip,s34,s567,s12,s134,zip,zip,zip,zip)

      scints(4,d7x12x56sl,iep)
c     & =loopI4(s56,s12,zip,s34,s347,s127,zip,zip,zip,zip,musq,iep)
     & =threemassbox(zip,s34,s56,s12,s347,s127,zip,zip,zip,zip)

      scints(4,d7x12x34sl,iep)
c     & =loopI4(zip,s12,s34,s56,s127,s567,zip,zip,zip,zip,musq,iep)
     & =threemassbox(zip,s12,s34,s56,s127,s567,zip,zip,zip,zip)

      scints(4,d27x1x34sl,iep)
     & =threemassbox(zip,s27,s56,s34,s127,s134,zip,zip,zip,zip)
c     & =loopI4(s56,s34,zip,s27,s127,s134,zip,zip,zip,zip,musq,iep)

      scints(3,c12x56sl,iep)=loopI3(s56,s347,s12,zip,zip,zip,musq,iep)
      scints(3,c12x34sl,iep)=loopI3(s34,s567,s12,zip,zip,zip,musq,iep)

      scints(2,b347sl,iep)=loopI2(s347,zip,zip,musq,iep)
      scints(2,b567sl,iep)=loopI2(s567,zip,zip,musq,iep)
      scints(2,b12sl,iep)=loopI2(s12,zip,zip,musq,iep)

c--- identify subleading integrals that have already been computed
      scints(3,c17x34sl,iep)=scints(3,c17x34,iep)
      scints(3,c27x56sl,iep)=scints(3,c27x56,iep)
      scints(3,c56x34sl,iep)=scints(3,c56x34,iep)
      scints(3,c34x56sl,iep)=scints(3,c56x34sl,iep)

      scints(2,b134sl,iep)=scints(2,b134,iep)
      scints(2,b256sl,iep)=scints(2,b256,iep)
      scints(2,b34sl,iep) =scints(2,b34,iep)
      scints(2,b17sl,iep) =scints(2,b17,iep)
      scints(2,b56sl,iep) =scints(2,b56,iep)
      scints(2,b27sl,iep)=scints(2,b27,iep)
      scints(2,b127sl,iep)=scints(2,b127,iep)
c      enddo

c      if (writescalars) then
c      write(67,*) 'scints(4,d2x3x4,0) ',scints(4,d2x3x4,0)
c      write(67,*) 'scints(4,d1x2x3,0) ',scints(4,d1x2x3,0)
c      write(67,*) 'scints(4,d1x23x4,0)',scints(4,d1x23x4,0)
c      write(67,*) 'scints(4,d1x2x56,0)',scints(4,d1x2x56,0)
c      write(67,*) 'scints(4,d12x3x4,0)',scints(4,d12x3x4,0)
c      write(67,*) 'done boxes'
c      write(67,*) 'scints(3,c23x4,0)  ',scints(3,c23x4,0)
c      write(67,*) 'scints(3,c1x23,0)  ',scints(3,c1x23,0)
c      write(67,*) 'scints(3,c2x3,0)   ',scints(3,c2x3,0)
c      write(67,*) 'scints(3,c2x56,0)  ',scints(3,c2x56,0)
c      write(67,*) 'scints(3,c3x4,0)   ',scints(3,c3x4,0)
c      write(67,*) 'scints(3,c12x3,0)  ',scints(3,c12x3,0)
c      write(67,*) 'scints(3,c1x2,0)   ',scints(3,c1x2,0)
c      write(67,*) 'scints(3,c1x256,0) ',scints(3,c1x256,0)
c      write(67,*) 'scints(3,c123x4,0) ',scints(3,c123x4,0)
c      write(67,*) 'scints(3,c12x56,0) ',scints(3,c12x56,0)
c      write(67,*) 'done triangles'
c      write(67,*) 'scints(2,b23,0)    ',scints(2,b23,0)
c      write(67,*) 'scints(2,b256,0)   ',scints(2,b256,0)
c      write(67,*) 'scints(2,b123,0)   ',scints(2,b123,0)
c      write(67,*) 'scints(2,b2x1m,0)  ',scints(2,b2x1m,0)
c      write(67,*) 'scints(2,b56,0)    ',scints(2,b56,0)
c      write(67,*) 'scints(2,b12,0)    ',scints(2,b12,0)
c      write(67,*) 'scints(2,b1256,0)  ',scints(2,b1256,0)
c      write(67,*) 'done bubbles'
c      write(67,*) 'scints(1,a0m,0)    ',scints(1,a0m,0)
c      write(67,*) 'done tadpole'
c      write(67,*)
c      endif

c--- triangles appearing in qloop
c      scints(3,c34x56,iep)=loopI3(s34,s56,s127,zip,zip,zip,musq,iep)
c      scints(3,c12x56sl,iep)=loopI3(s56,s347,s12,zip,zip,zip,musq,iep)
c      scints(3,c12x34sl,iep)=loopI3(s34,s567,s12,zip,zip,zip,musq,iep)
c      tmp(1)=s34**2
c      tmp(2)=s56**2
c      tmp(3)=s127**2
c      tmp(4)=abs(2d0*s56*s34)
c      tmp(5)=abs(2d0*s56*s127)
c      tmp(6)=abs(2d0*s34*s127)
c      write(6,*) 'tri1',((s34-s56-s127)**2-4d0*s56*s127)
c     &/max(tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6))
c      tmp(1)=s56**2
c      tmp(2)=s347**2
c      tmp(3)=s12**2
c      tmp(4)=abs(2d0*s56*s347)
c      tmp(5)=abs(2d0*s56*s12)
c      tmp(6)=abs(2d0*s12*s347)
c      write(6,*) 'tri2',((s56-s347-s12)**2-4d0*s347*s12)
c     &/max(tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6))
c      tmp(1)=s34**2
c      tmp(2)=s567**2
c      tmp(3)=s12**2
c      tmp(4)=abs(2d0*s34*s567)
c      tmp(5)=abs(2d0*s34*s12)
c      tmp(6)=abs(2d0*s12*s567)
c      write(6,*) 'tri3',((s34-s567-s12)**2-4d0*s567*s12)
c     &/max(tmp(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6))

c Extra integrals needed for up-type W-Z (34<->56)
c--- integrals appearing at leading color
      scints(4,d7x1x56,iep)
     & =hard(zip,zip,s56,s234,s17,s156,zip,zip,zip,zip)

      scints(4,d1x27x34,iep)
     & =threemassbox(zip,s27,s34,s56,s127,s156,zip,zip,zip,zip)

      scints(4,d7x2x34,iep)
     & =hard(zip,zip,s34,s156,s27,s234,zip,zip,zip,zip)

      scints(4,d2x17x56,iep)
     & =threemassbox(zip,s17,s56,s34,s127,s234,zip,zip,zip,zip)

      scints(3,c56x17,iep)=loopI3(s56,s17,s234,zip,zip,zip,musq,iep)
      scints(3,c34x27,iep)=loopI3(s34,s27,s156,zip,zip,zip,musq,iep)

      scints(2,b156,iep)=loopI2(s156,zip,zip,musq,iep)
      scints(2,b234,iep)=loopI2(s234,zip,zip,musq,iep)


      scints(4,d2x34x7sl,iep)
     & =easy(zip,s34,zip,s156,s347,s234,zip,zip,zip,zip)

      scints(4,d2x7x34sl,iep)
     & =hard(s34,zip,zip,s156,s347,s27,zip,zip,zip,zip)

      scints(4,d1x56x7sl,iep)
     & =easy(zip,s56,zip,s234,s567,s156,zip,zip,zip,zip)

      scints(4,d1x7x56sl,iep)
     & =hard(s56,zip,zip,s234,s567,s17,zip,zip,zip,zip)

      scints(4,d1x27x34sl,iep)
     & =threemassbox(zip,s27,s34,s56,s127,s156,zip,zip,zip,zip)

      scints(4,d2x1x7sl,iep)
     & =onemassbox(zip,zip,zip,s127,s17,s12,zip,zip,zip,zip)

      scints(4,d1x2x7sl,iep)
     & =onemassbox(zip,zip,zip,s127,s27,s12,zip,zip,zip,zip)

      scints(4,d2x1x56sl,iep)
     & =hard(s56,zip,zip,s347,s156,s12,zip,zip,zip,zip)

      scints(4,d7x34x12sl,iep)
     & =threemassbox(zip,s56,s12,s34,s567,s347,zip,zip,zip,zip)

      scints(4,d1x2x34sl,iep)
     & =hard(zip,zip,s34,s567,s12,s234,zip,zip,zip,zip)

      scints(4,d7x12x56sl,iep)
     & =threemassbox(zip,s34,s56,s12,s347,s127,zip,zip,zip,zip)
      scints(4,d12x7x34sl,iep)=scints(4,d7x12x56sl,iep)

      scints(4,d7x12x34sl,iep)
     & =threemassbox(zip,s12,s34,s56,s127,s567,zip,zip,zip,zip)

      scints(4,d17x2x34sl,iep)
     & =threemassbox(zip,s17,s56,s34,s127,s234,zip,zip,zip,zip)

      scints(3,c56x17sl,iep)=scints(3,c56x17,iep)
      scints(3,c17x56sl,iep)=scints(3,c56x17sl,iep)
      scints(3,c34x27sl,iep)=scints(3,c34x27,iep)

      scints(2,b156sl,iep)=scints(2,b156,iep)
      scints(2,b234sl,iep)=scints(2,b234,iep)

      return
      end

