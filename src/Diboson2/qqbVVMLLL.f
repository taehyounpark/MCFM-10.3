!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qqbVVMLLL(order,ABCF,qtype,p1,p2,p5,p6,p7,p8,qqbAj)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'kprocess.f'
      include 'mxpart.f'
      include 'zprods_com.f'
      include 'sprods_com.f'
      include 'ABCF.f'
      include 'first.f'
      complex(dp):: qqbVVMLLL,qqbAj(0:2,4,10),
     & zab2,zba1_3_2,E(9),stringLLL(9),AjF(1:10)
      integer p1,p2,p5,p6,p7,p8,ii,order,qtype
      integer,parameter:: dntype=1,uptype=2
      real(dp):: s12,t,u,p3sq,p4sq

      zab2(p1,p5,p6,p2)=za(p1,p5)*zb(p5,p2)+za(p1,p6)*zb(p6,p2)

      if (first) then
          call Qformfill
          first=.false.
      endif

      qqbVVMLLL=czip   !  initialize to zero
      if (ABCF == A) then
          if ((kcase == kWWqqbr).and.(qtype == dntype)) return
      elseif (ABCF == B) then
          if ((kcase == kWWqqbr).and.(qtype == uptype)) return
      elseif (ABCF==C) then
c no C-type diagrams for WZ or in 0 or 1 loops
          if ((kcase==kWZbbar) .or. (order < 2)) return
      elseif (ABCF==F) then
c no F-type diagrams for ZZ
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


c     1503.04812v2, Eq.3.19, signs adjusted to account for -p1 and -p2

c     spinor content <5><7><2> [6],[8],[1] -->
c     spinor content <3><5><2> [4],[6],[1] in 123456 notation

      zba1_3_2=zab2(p2,p5,p6,p1)
      stringLLL(1)= -za(p1,p5)*za(p1,p7)*zb(p1,p6)*zb(p1,p8)*zba1_3_2
      stringLLL(2)= -za(p1,p5)*za(p2,p7)*zb(p1,p6)*zb(p2,p8)*zba1_3_2
      stringLLL(3)= -za(p2,p5)*za(p1,p7)*zb(p2,p6)*zb(p1,p8)*zba1_3_2
      stringLLL(4)= -za(p2,p5)*za(p2,p7)*zb(p2,p6)*zb(p2,p8)*zba1_3_2
      stringLLL(5)= -za(p5,p7)*zb(p6,p8)*zba1_3_2
      stringLLL(6)= za(p1,p5)*za(p2,p7)*zb(p1,p6)*zb(p1,p8)
      stringLLL(7)= za(p2,p5)*za(p2,p7)*zb(p2,p6)*zb(p1,p8)
      stringLLL(8)= za(p2,p5)*za(p1,p7)*zb(p1,p6)*zb(p1,p8)
      stringLLL(9)= za(p2,p5)*za(p2,p7)*zb(p1,p6)*zb(p2,p8)


      do ii=1,9
         qqbVVMLLL=qqbVVMLLL+stringLLL(ii)*E(ii)
      enddo

      end
