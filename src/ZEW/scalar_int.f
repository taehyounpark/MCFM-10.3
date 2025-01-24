!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function xI1(msq,mu2,ep)
      use loopI1_generic
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: msq,mu2,xI1
      integer ep

      xI1=real(loopI1(msq,mu2,ep),dp)

      end function xI1


      function xI2(psq,m1sq,m2sq,mu2,ep)
      use loopI1_generic
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: psq,m1sq,m2sq,mu2,xI2
      integer ep

      xI2=real(loopI2(psq,m1sq,m2sq,mu2,ep),dp)

      end function xI2


      function xI3(p1sq,p2sq,p3sq,m1sq,m2sq,m3sq,mu2,ep)
      use loopI1_generic
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: p1sq,p2sq,p3sq,m1sq,m2sq,m3sq,mu2,xI3
      integer ep

      xI3=real(loopI3(p1sq,p2sq,p3sq,m1sq,m2sq,m3sq,mu2,ep),dp)

      end function xI3


      function xI4(p1sq,p2sq,p3sq,p4sq,s12,s23,
     &     m1sq,m2sq,m3sq,m4sq,mu2,ep)
      use loopI1_generic
      use loopI2_generic
      use loopI3_generic
      use loopI4_generic
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: p1sq,p2sq,p3sq,p4sq,s12,s23,
     &     m1sq,m2sq,m3sq,m4sq,mu2,xI4
      integer ep

      xI4=real(loopI4(p1sq,p2sq,p3sq,p4sq,s12,s23,m1sq,m2sq,m3sq,m4sq,
     &     mu2,ep),dp)

      end function xI4


      function db0(psq,m1sq,m2sq)
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: psq,m1sq,m2sq,db0
      complex(dp):: mcfmccdb0

c FF call
c      call ffxdb0(cdb0,cdb0p,psq,m1sq,m2sq,ier)
c      db0=real(cdb0,dp)
c MCFM implementation
      db0=real(mcfmccdb0(psq,m1sq,m2sq),dp)

c      write(6,*) 'db0 comparison',psq,m1sq,m2sq,db0/real(cdb0,dp)

      end function db0


      function ccdb0(psq,m1sq,m2sq)
      implicit none
      include 'types.f'
      include 'scalarselect.f'
      real(dp):: psq,m1sq,m2sq
      complex(dp):: ccdb0,mcfmccdb0

c FF call
c      call ffxdb0(cdb0,cdb0p,psq,m1sq,m2sq,ier)
c      ccdb0 = cdb0
c MCFM implementation
      ccdb0 = mcfmccdb0(psq,m1sq,m2sq)

c      write(6,*) 'ccdb0 comparison',psq,m1sq,m2sq,ccdb0/cdb0

      end function ccdb0
