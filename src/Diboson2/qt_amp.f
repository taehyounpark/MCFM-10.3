!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      module qT_amp
      implicit none
      private
      public:: qTsum    !qT00,qT10,qT20,qT11
      public:: DelI1,DelI2
      include 'types.f'
      include 'constants.f'
      include 'zeta.f'
      include 'nf.f'
      complex(dp),parameter::
     & DelI1=CF*(-half*pisq+impi*three/two),
     & DelI2=CF*CA*(-607./162-1181./432*pisq+187./72*zeta3+7./96*pi**4
     & +impi*(961./216+11./72*pisq-13./2*zeta3))
     & +CF**2*(-9./8*pisq+1./8*pi**4+impi*(3./8-5./4*pisq+6d0*zeta3))
     & +NF*CF*(41./81+97./216*pisq-17./36*zeta3
     & +impi*(-65./108-1./36*pisq))

      contains

      function qTsum(amp0,amp1)
c     returns helicity-summed interference of amp0 and amp1, defined in MXXX.f
c      include 'types.f'
      integer:: n1,n2,n3
      complex(dp):: qTsum,amp0(2,2,2),amp1(2,2,2)
      qTsum=czip
      do n3=1,2
      do n2=1,2
      do n1=1,2
      qTsum=qTsum+amp0(n1,n2,n3)*conjg(amp1(n1,n2,n3))
      enddo
      enddo
      enddo
      end function qTsum

c---      function qT00(qqb,qType,p1,p2,p5,p6,p7,p8)
c---!     returns tree level amplitude squared, using amplitudes defined in MXXX.f
c---      include 'types.f'
c---      integer:: p1,p2,p5,p6,p7,p8,n1,n2,n3
c---      integer:: qType
c---      complex(dp):: qT00,amp0(2,2,2)
c---      logical::qqb
c---      call MXXX(0,qqb,qType,p1,p2,p5,p6,p7,p8,amp0)
c---      qT00=czip
c---      do n3=1,2
c---      do n2=1,2
c---      do n1=1,2
c---      qT00=qT00+abs(amp0(n1,n2,n3))**2
c---      enddo
c---      enddo
c---      enddo
c---      end function qT00
c---
c---      function qT10(qqb,qType,p1,p2,p5,p6,p7,p8)
c---      include 'types.f'
c---      integer:: p1,p2,p5,p6,p7,p8,n1,n2,n3
c---      integer:: qType
c---      complex(dp):: qT10,amp0(2,2,2),amp1(2,2,2)
c---      logical::qqb
c---
c---      call MXXX(0,qqb,qType,p1,p2,p5,p6,p7,p8,amp0)
c---      call MXXX(1,qqb,qType,p1,p2,p5,p6,p7,p8,amp1)
c---      qT10=czip
c---      do n3=1,2
c---      do n2=1,2
c---      do n1=1,2
c---      qT10=qT10+two*amp1(n1,n2,n3)*conjg(amp0(n1,n2,n3))
c---      enddo
c---      enddo
c---      enddo
c---      end function qT10
c---
c---      function qT20(qqb,qType,p1,p2,p5,p6,p7,p8)
c---      include 'types.f'
c---      integer:: p1,p2,p5,p6,p7,p8,n1,n2,n3
c---      integer:: qType
c---      complex(dp):: qT20,amp0(2,2,2),amp2(2,2,2)
c---      logical::qqb
c---      call MXXX(0,qqb,qType,p1,p2,p5,p6,p7,p8,amp0)
c---      call MXXX(2,qqb,qType,p1,p2,p5,p6,p7,p8,amp2)
c---      qT20=czip
c---      do n3=1,2
c---      do n2=1,2
c---      do n1=1,2
c---      qT20=qT20+two*amp2(n1,n2,n3)*conjg(amp0(n1,n2,n3))
c---      enddo
c---      enddo
c---      enddo
c---      end function qT20
c---
c---
c---      function qT11(qqb,qType,p1,p2,p5,p6,p7,p8)
c---      include 'types.f'
c---      integer:: p1,p2,p5,p6,p7,p8,n1,n2,n3
c---      integer:: qType
c---      complex(dp):: qT11,amp1(2,2,2)
c---      logical::qqb
c---      call MXXX(1,qqb,qType,p1,p2,p5,p6,p7,p8,amp1)
c---      qT11=czip
c---      do n3=1,2
c---      do n2=1,2
c---      do n1=1,2
c---      qT11=qT11+abs(amp1(n1,n2,n3))**2
c---      enddo
c---      enddo
c---      enddo
c---      end function qT11

      end module qT_amp

