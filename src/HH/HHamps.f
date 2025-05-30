!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

       subroutine HHamps(p1,p2,p3,p4,gauge)
      use loopI3_generic
      use loopI4_generic
      implicit none
      include 'scalarselect.f'
      include 'types.f'
c--    Formula taken from Glover and van der Bij
c--    Nucl. Phys. B309 (1988) 202

      include 'constants.f'
      include 'cplx.h'
      include 'masses.f'
      include 'scale.f'
      include 'yukawas.f'
      integer n,nmin
      real(dp):: p1(4),p2(4),p3(4),p4(4),ss,tt,uu,
     & p3sq,p4sq,mQsq,mhsq,mQsq_yuk
      complex(dp):: triangle(1,2),box(1,2),gauge(1,2),
     & D123,D213,D132,C12,C23,C13,C34

c----Signs for momenta chosen such that p1+p2+p3+p4=0
      ss=(p1(4)+p2(4))**2
     &  -(p1(1)+p2(1))**2-(p1(2)+p2(2))**2-(p1(3)+p2(3))**2
      tt=(p1(4)+p3(4))**2
     &  -(p1(1)+p3(1))**2-(p1(2)+p3(2))**2-(p1(3)+p3(3))**2
      uu=(p2(4)+p3(4))**2
     &  -(p2(1)+p3(1))**2-(p2(2)+p3(2))**2-(p2(3)+p3(3))**2
      p3sq=p3(4)**2-p3(1)**2-p3(2)**2-p3(3)**2
      p4sq=p4(4)**2-p4(1)**2-p4(2)**2-p4(3)**2
      mhsq=hmass**2

      gauge(:,:)=czip
      nmin=4
      if (mc_yuk == 0) nmin=5
      if (mb_yuk == 0) nmin=6
      do n=nmin,6
      if (n == 4) then
        mqsq=mc**2
        mqsq_yuk=mc_yuk**2
      endif
      if (n == 5) then
        mqsq=mb**2
        mqsq_yuk=mb_yuk**2
      endif
      if (n == 6) then
        mqsq=mt**2
        mqsq_yuk=mt_yuk**2
      endif

      C12=loopI3(0._dp,0._dp ,ss,mQsq,mQsq,mQsq,musq,0)
      C23=loopI3(0._dp,p3sq,uu,mQsq,mQsq,mQsq,musq,0)
      C13=loopI3(0._dp,p3sq,tt,mQsq,mQsq,mQsq,musq,0)
      C34=loopI3(p3sq,p4sq,ss,mQsq,mQsq,mQsq,musq,0)
      D123=loopI4(0._dp,0._dp,p3sq,p4sq,ss,uu,mQsq,mQsq,mQsq,mQsq,musq,0)
      D213=loopI4(0._dp,0._dp,p3sq,p4sq,ss,tt,mQsq,mQsq,mQsq,mQsq,musq,0)
      D132=loopI4(0._dp,p3sq,0._dp,p4sq,tt,uu,mQsq,mQsq,mQsq,mQsq,musq,0)

      triangle(1,1)=-12._dp*mhsq*sqrt(mQsq*mQsq_yuk)/cplx2(ss-mhsq,hmass*hwidth)
     & *(2._dp+(4._dp*mQsq-ss)*C12)
      triangle(1,2)=czip

      box(1,1)=-4._dp*mQsq_yuk*(
     & mQsq*(8._dp*mQsq-ss-2._dp*mhsq)*(D123+D213+D132)
     & +(uu*tt-mhsq**2)/ss*(4._dp*mQsq-mhsq)*D132+2._dp+4._dp*mQsq*C12
     & +2._dp/ss*(mhsq-4._dp*mQsq)*((tt-mhsq)*C13+(uu-mhsq)*C23))

      box(1,2)=-2._dp*mQsq_yuk*(
     & 2._dp*(8._dp*mQsq+ss-2._dp*mhsq)*(mQsq*(D123+D213+D132)-C34)
     & -2._dp*(ss*C12+(tt-mhsq)*C13+(uu-mhsq)*C23)
     & +1._dp/(uu*tt-mhsq**2)*(ss*uu*(8._dp*uu*mQsq-uu**2-mhsq**2)*D123
     & +ss*tt*(8._dp*tt*mQsq-tt**2-mhsq**2)*D213
     & +(8._dp*mQsq+ss-2._dp*mhsq)*(ss*(ss-2._dp*mhsq)*C12
     & +ss*(ss-4._dp*mhsq)*C34
     & +2._dp*tt*(mhsq-tt)*C13+2._dp*uu*(mhsq-uu)*C23)))

      gauge(1,1)=gauge(1,1)+triangle(1,1)+box(1,1)
      gauge(1,2)=gauge(1,2)+triangle(1,2)+box(1,2)

      enddo

      return
      end
