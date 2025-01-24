!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qtassemble(order,qtcut,tBa0,tBb0,tBa1,tBb1,
     & tBa2,tBb2,tSb1,tSb2,H)
      implicit none
      include 'types.f'
c---- Given beam, soft and hard functions, calculates
c---- O(alphas^order) correction after all convolutions
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'zeta.f'
      include 'kpart.f'
      integer:: order
      real(dp)::qtassemble,qtcut,
     & tBa0,tBb0,tBa1(0:1),tBb1(0:1),tBa2(0:2),tBb2(0:2),
     & tSb1(0:2),tSb2(0:4),H(2),
     & Lqt2cut,full1(0:2),full2(0:4)

      Lqt2cut = 2*log(qtcut/scale)

      if (coeffonly) then
        qtassemble=zip
      else
        qtassemble=tBa0*tBb0
      endif
c     integral of -L0(qt2/mu2) from mu2 to qtcut**2 =   -Lqt2cut
c     integral of +2*L1(qt2/mu2) from mu2 to qtcut**2 = +Lqt2cut^2
c     integral of -3*L2(qt2/mu2) from mu2 to qtcut**2 = -Lqt2cut^3
c     integral of +4*L3(qt2/mu2) from mu2 to qtcut**2 = +Lqt2cut^4
      if ( (order ==1) .or.
     &    ((order ==2) .and. (coeffonly .eqv. .false.)) ) then

      full1(0)=H(1)*tBa0*tBb0
     & +tBa0*tBb0*tSb1(0)+tBa0*tBb1(0)+tBa1(0)*tBb0
      full1(1)=tBa0*tBb0*tSb1(1)+tBa0*tBb1(1)+tBa1(1)*tBb0
      full1(2)=tBa0*tBb0*tSb1(2)
      qtassemble=qtassemble+ason4pi
     & *(full1(0)+Lqt2cut*(-full1(1)+Lqt2cut*full1(2)))
c      write(6,*) 'ason4pi*full1(2)',ason4pi*full1(2)
c      write(6,*) 'ason4pi*full1(1)',ason4pi*full1(1)
c      write(6,*) 'ason4pi*full1(0)',ason4pi*full1(0)
c      pause

c      write(6,*) 'QT:  full1(2)',full1(2)
c      write(6,*) 'QT: -full1(1)',-full1(1)
c      write(6,*) 'QT:  full1(0)',full1(0)
c      write(6,*)

      endif

      if (order > 1) then
      full2(0)=H(1)*tBa0*tBb0*tSb1(0)+H(1)*tBa0*tBb1(0)
     & +H(1)*tBa1(0)*tBb0+H(2)*tBa0*tBb0+tBa0*tBb0*tSb2(0)
     & -4*tBa0*tBb0*tSb2(3)*zeta3
     & +tBa0*tBb1(0)*tSb1(0)-4*tBa0*tBb1(1)*tSb1(2)*zeta3+tBa0*tBb2(0)
     & +tBa1(0)*tBb0*tSb1(0)+tBa1(0)*tBb1(0)
     & -4*tBa1(1)*tBb0*tSb1(2)*zeta3+tBa2(0)*tBb0

      full2(1)=H(1)*tBa0*tBb0*tSb1(1)+H(1)*tBa0*tBb1(1)
     & +H(1)*tBa1(1)*tBb0+tBa0*tBb0*tSb2(1)+tBa0*tBb1(0)*tSb1(1)
     & +tBa0*tBb1(1)*tSb1(0)+tBa0*tBb2(1)+tBa1(0)*tBb0*tSb1(1)
     & +tBa1(0)*tBb1(1)+tBa1(1)*tBb0*tSb1(0)+tBa1(1)*tBb1(0)
     & -16*tBa0*tBb0*tSb2(4)*zeta3+tBa2(1)*tBb0

      full2(2)=H(1)*tBa0*tBb0*tSb1(2)+tBa0*tBb0*tSb2(2)
     & +tBa0*tBb1(0)*tSb1(2)+tBa0*tBb1(1)*tSb1(1)+tBa0*tBb2(2)
     & +tBa1(0)*tBb0*tSb1(2)+tBa1(1)*tBb0*tSb1(1)
     & +tBa1(1)*tBb1(1)+tBa2(2)*tBb0

      full2(3)=tBa0*tBb0*tSb2(3)+tBa0*tBb1(1)*tSb1(2)
     & +tBa1(1)*tBb0*tSb1(2)

      full2(4)=tBa0*tBb0*tSb2(4)

c      write(6,*) 'QT:  full2(4)',full2(4)
c      write(6,*) 'QT: -full2(3)',-full2(3)
c      write(6,*) 'QT:  full2(2)',full2(2)
c      write(6,*) 'QT: -full2(1)',-full2(1)
c      write(6,*) 'QT:  full2(0)',full2(0)
c      write(6,*)

c      write(6,*) '      ason4pi*H(1)',ason4pi*H(1)
c      write(6,*) '   ason4pi**2*H(2)',ason4pi**2*H(2)
c      write(6,*) '   ason4pi*tSb1(0)',ason4pi*tSb1(0)
c      write(6,*) '   ason4pi*tSb1(1)',ason4pi*tSb1(1)
c      write(6,*) '   ason4pi*tSb1(2)',ason4pi*tSb1(2)

c      write(6,*) 'ason4pi**2*tSb2(0)',ason4pi**2*tSb2(0)
c      write(6,*) 'ason4pi**2*tSb2(1)',ason4pi**2*tSb2(1)
c      write(6,*) 'ason4pi**2*tSb2(2)',ason4pi**2*tSb2(2)
c      write(6,*) 'ason4pi**2*tSb2(3)',ason4pi**2*tSb2(3)
c      write(6,*) 'ason4pi**2*tSb2(4)',ason4pi**2*tSb2(4)

c      write(6,*) '   ason4pi*tBa1(0)',ason4pi*tBa1(0)
c      write(6,*) '   ason4pi*tBa1(1)',ason4pi*tBa1(1)


c      write(6,*) 'ason4pi**2*tBa2(0)',ason4pi**2*tBa2(0)
c      write(6,*) 'ason4pi**2*tBa2(1)',ason4pi**2*tBa2(1)
c      write(6,*) 'ason4pi**2*tBa2(2)',ason4pi**2*tBa2(2)

c      write(6,*) '   ason4pi*tBb1(0)',ason4pi*tBb1(0)
c      write(6,*) '   ason4pi*tBb1(1)',ason4pi*tBb1(1)

c      write(6,*) 'ason4pi**2*tBb2(0)',ason4pi**2*tBb2(0)
c      write(6,*) 'ason4pi**2*tBb2(1)',ason4pi**2*tBb2(1)
c      write(6,*) 'ason4pi**2*tBb2(2)',ason4pi**2*tBb2(2)

c      write(6,*) 'ason4pi**2*full2(4)',ason4pi**2*full2(4)
c      write(6,*) 'ason4pi**2*full2(3)',ason4pi**2*full2(3)
c      write(6,*) 'ason4pi**2*full2(2)',ason4pi**2*full2(2)
c      write(6,*) 'ason4pi**2*full2(1)',ason4pi**2*full2(1)
c      write(6,*) 'ason4pi**2*full2(0)',ason4pi**2*full2(0)
c      pause

      qtassemble=qtassemble+ason4pi**2
     & *(full2(0)+Lqt2cut*(-full2(1)+Lqt2cut*(full2(2)
     & +Lqt2cut*(-full2(3)+Lqt2cut*full2(4)))))

c      qtassemble=qtassemble+ason4pi**2
c     & *(full2(0)+Lqt2cut*(-0*full2(1)+Lqt2cut*(0*full2(2)
c     & +Lqt2cut*(-0*full2(3)+Lqt2cut*0*full2(4)))))
      endif

      return
      end
