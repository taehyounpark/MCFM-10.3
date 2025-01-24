!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qqbVVMRLL(order,ABCF,qtype,p1,p2,p5,p6,p7,p8,qqbAj)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'kprocess.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'ABCF.f'
      include 'first.f'
      integer p1,p2,p5,p6,p7,p8,order,qtype,ii
      integer,parameter:: dntype=1,uptype=2
      real(dp):: s12,t,u,p3sq,p4sq
      complex(dp):: qqbVVMRLL,qqbAj(0:2,4,10),
     & zab2,zba2_3_1,E(9),stringRLL(9),AjF(1:10)

      zab2(p1,p5,p6,p2)=za(p1,p5)*zb(p5,p2)+za(p1,p6)*zb(p6,p2)

      if (first) then
          call Qformfill
          first=.false.
      endif

      qqbVVMRLL=czip   !  initialize to zero
      if (ABCF == A) then
          if ((kcase == kWWqqbr).and.(qtype == dntype)) return
      elseif (ABCF == B) then
          if ((kcase == kWWqqbr).and.(qtype == uptype)) return
      elseif (ABCF==C) then
          ! no C-type diagrams for WZ or in 0 or 1 loops
          if ((kcase==kWZbbar) .or. (order < 2)) return
      elseif (ABCF==F) then ! no F-type diagrams for ZZ
          if (kcase==kZZlept) return
      endif

      s12=s(p1,p2)
      t=s(p1,p5)+s(p1,p6)+s(p5,p6)
      u=s(p2,p5)+s(p2,p6)+s(p5,p6)
      p3sq=s(p5,p6)
      p4sq=s(p7,p8)

      IF (ABCF==F) then
      call qqbAjF(order,AjF)
      qqbAj(order,F,:)=AjF(:)
      endif

c     1503.04812v2, Eq.3.21
      E(1)=qqbAj(order,ABCF,1)
      E(2)=qqbAj(order,ABCF,2)
     & +2._dp/s12*(qqbAj(order,ABCF,9)-qqbAj(order,ABCF,10))
      E(3)=qqbAj(order,ABCF,3)
     & -2._dp/s12*(qqbAj(order,ABCF,9)-qqbAj(order,ABCF,10))
      E(4)=qqbAj(order,ABCF,4)
      E(5)=2._dp*(qqbAj(order,ABCF,9)+qqbAj(order,ABCF,10))
      E(6)=2._dp*qqbAj(order,ABCF,7)
     & +2._dp*(u-p3sq)/s12*(qqbAj(order,ABCF,9)-qqbAj(order,ABCF,10))
      E(7)=2._dp*qqbAj(order,ABCF,8)
     &-2._dp*(t-p3sq)/s12*(qqbAj(order,ABCF,9)-qqbAj(order,ABCF,10))
      E(8)=2._dp*qqbAj(order,ABCF,5)-2._dp/s12
     & *((u-s12-p3sq)*qqbAj(order,ABCF,9)+(t-p4sq)*qqbAj(order,ABCF,10))
      E(9)=2._dp*qqbAj(order,ABCF,6)-2._dp/s12
     & *((t-s12-p3sq)*qqbAj(order,ABCF,10)+(u-p4sq)*qqbAj(order,ABCF,9))


c     1503.04812v2, Eq.3.20, signs adjusted to account for -p1 and -p2
      zba2_3_1=zab2(p1,p5,p6,p2)
      stringRLL(1)= -za(p1,p5)*za(p1,p7)*zb(p1,p6)*zb(p1,p8)*zba2_3_1
      stringRLL(2)= -za(p1,p5)*za(p2,p7)*zb(p1,p6)*zb(p2,p8)*zba2_3_1
      stringRLL(3)= -za(p2,p5)*za(p1,p7)*zb(p2,p6)*zb(p1,p8)*zba2_3_1
      stringRLL(4)= -za(p2,p5)*za(p2,p7)*zb(p2,p6)*zb(p2,p8)*zba2_3_1
      stringRLL(5)= -za(p5,p7)*zb(p6,p8)*zba2_3_1
      stringRLL(6)= za(p1,p5)*za(p1,p7)*zb(p1,p6)*zb(p2,p8)
      stringRLL(7)= za(p2,p5)*za(p1,p7)*zb(p2,p6)*zb(p2,p8)
      stringRLL(8)= za(p1,p5)*za(p1,p7)*zb(p2,p6)*zb(p1,p8)
      stringRLL(9)= za(p1,p5)*za(p2,p7)*zb(p2,p6)*zb(p2,p8)

      do ii=1,9
         qqbVVMRLL=qqbVVMRLL+stringRLL(ii)*E(ii)
      enddo

      return
      end
