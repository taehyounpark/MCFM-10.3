!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function GLYassemble(order,qtcut,pp,Ba0,Bb0,Ba1,Bb1,
     & Ba2,Bb2,F1,F2,H)
      implicit none
      include 'types.f'
c---- Given beam, soft and hard functions, calculates
c---- O(alphas^order) correction after all convolutions
c----
c---- pp is an integer label (= qq or gg) that specifies
c---- the identities of the partons in the LO process
      include 'constants.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'zeta.f'
      include 'kpart.f'
      integer:: order,pp
      real(dp)::GLYassemble,qtcut,
     & Ba0,Bb0,Ba1(0:2),Bb1(0:2),Ba2(0:4),Bb2(0:4),
     & F1(0:1,0:2),F2(0:1,0:4),H(2),
     & Lqt2cut,full1(0:2),full2(0:4)

      Lqt2cut = 2*log(qtcut/scale)

      if (coeffonly) then
        GLYassemble=zip
      else
        GLYassemble=Ba0*Bb0
      endif
c     integral of -L0(qt2/mu2) from mu2 to qtcut**2 =   -Lqt2cut
c     integral of +2*L1(qt2/mu2) from mu2 to qtcut**2 = +Lqt2cut^2
c     integral of -3*L2(qt2/mu2) from mu2 to qtcut**2 = -Lqt2cut^3
c     integral of +4*L3(qt2/mu2) from mu2 to qtcut**2 = +Lqt2cut^4
      if ( (order ==1) .or.
     &    ((order ==2) .and. (coeffonly .eqv. .false.)) ) then
      full1(2)=+Ba0*Bb0*F1(pp,2)+Ba0*Bb1(2)+Ba1(2)*Bb0;;
      full1(1)=+Ba0*Bb0*F1(pp,1)+Ba0*Bb1(1)+Ba1(1)*Bb0;
      full1(0)=H(1)*Ba0*Bb0+Ba0*Bb0*F1(pp,0)+Ba0*Bb1(0)+Ba1(0)*Bb0;

      GLYassemble=GLYassemble+ason4pi
     & *(full1(0)+Lqt2cut*(-full1(1)+Lqt2cut*full1(2)))
c      write(6,*) 'Ba0',Ba0
c      write(6,*) 'Bb0',Bb0
c      write(6,*) F1(pp,2),Ba1(2)/Ba0,Bb1(2)/Bb0
c      write(6,*) 'ason4pi',ason4pi
c      write(6,*) 'Lqt2cut',Lqt2cut
c      write(6,*) 'scale',scale
c      write(6,*) 'full1(2)',full1(2)
c      write(6,*) 'full1(1)',full1(1)
c      write(6,*) 'full1(0)',full1(0)
c      pause

c      write(6,*) 'GLY: full1(2)',full1(2)
c      write(6,*) 'GLY:-full1(1)',-full1(1)
c      write(6,*) 'GLY: full1(0)',full1(0)

      endif

      if (order > 1) then
      full2(4)=+Ba0*Bb0*F2(pp,4)+Ba0*Bb1(2)*F1(pp,2)+Ba0*Bb2(4)
     & +Ba1(2)*Bb0*F1(pp,2)+Ba1(2)*Bb1(2)+Ba2(4)*Bb0;

      full2(3)=+Ba0*Bb0*F2(pp,3)+Ba0*Bb1(1)*F1(pp,2)+Ba0*Bb1(2)*F1(pp,1)
     & +Ba0*Bb2(3)+Ba1(1)*Bb0*F1(pp,2)+Ba1(1)*Bb1(2)+Ba1(2)*Bb0*F1(pp,1)
     & +Ba1(2)*Bb1(1)+Ba2(3)*Bb0

      full2(2)=+H(1)*Ba0*Bb0*F1(pp,2)+H(1)*Ba0*Bb1(2)+H(1)*Ba1(2)*Bb0
     & +Ba0*Bb0*F2(pp,2)+Ba0*Bb1(0)*F1(pp,2)+Ba0*Bb1(1)*F1(pp,1)
     & +Ba0*Bb1(2)*F1(pp,0)+Ba0*Bb2(2)+Ba1(0)*Bb0*F1(pp,2)+Ba1(0)*Bb1(2)
     & +Ba1(1)*Bb0*F1(pp,1)+Ba1(1)*Bb1(1)+Ba1(2)*Bb0*F1(pp,0)
     & +Ba1(2)*Bb1(0)+Ba2(2)*Bb0

      full2(1)=-(-H(1)*Ba0*Bb0*F1(pp,1)-H(1)*Ba0*Bb1(1)-H(1)*Ba1(1)*Bb0
     & -Ba0*Bb0*F2(pp,1)+16*Ba0*Bb0*F2(pp,4)*zeta3-Ba0*Bb1(0)*F1(pp,1)
     & -Ba0*Bb1(1)*F1(pp,0)+16*Ba0*Bb1(2)*F1(pp,2)*zeta3
     & -Ba0*Bb2(1)+16*Ba0*Bb2(4)*zeta3-Ba1(0)*Bb0*F1(pp,1)
     & -Ba1(0)*Bb1(1)-Ba1(1)*Bb0*F1(pp,0)-Ba1(1)*Bb1(0)
     & +16*Ba1(2)*Bb0*F1(pp,2)*zeta3
     & +16*Ba1(2)*Bb1(2)*zeta3-Ba2(1)*Bb0+16*Ba2(4)*Bb0*zeta3)

      full2(0)=H(1)*Ba0*Bb0*F1(pp,0)+H(1)*Ba0*Bb1(0)+H(1)*Ba1(0)*Bb0
     & +H(2)*Ba0*Bb0+Ba0*Bb0*F2(pp,0)-4*Ba0*Bb0*F2(pp,3)*zeta3
     & +Ba0*Bb1(0)*F1(pp,0)-4*Ba0*Bb1(1)*F1(pp,2)*zeta3
     & -4*Ba0*Bb1(2)*F1(pp,1)*zeta3+Ba0*Bb2(0)-4*Ba0*Bb2(3)*zeta3
     & +Ba1(0)*Bb0*F1(pp,0)+Ba1(0)*Bb1(0)
     & -4*Ba1(1)*Bb0*F1(pp,2)*zeta3-4*Ba1(1)*Bb1(2)*zeta3
     & -4*Ba1(2)*Bb0*F1(pp,1)*zeta3
     & -4*Ba1(2)*Bb1(1)*zeta3+Ba2(0)*Bb0-4*Ba2(3)*Bb0*zeta3

c      write(6,*) '      ason4pi*H(1)',ason4pi*H(1)
c      write(6,*) '   ason4pi**2*H(2)',ason4pi**2*H(2)
c      write(6,*) '   ason4pi*F1(pp,0)',ason4pi*F1(pp,0)
c      write(6,*) '   ason4pi*F1(pp,1)',ason4pi*F1(pp,1)
c      write(6,*) '   ason4pi*F1(pp,2)',ason4pi*F1(pp,2)

c      write(6,*) 'ason4pi**2*F2(pp,0)',ason4pi**2*F2(pp,0)
c      write(6,*) 'ason4pi**2*F2(pp,1)',ason4pi**2*F2(pp,1)
c      write(6,*) 'ason4pi**2*F2(pp,2)',ason4pi**2*F2(pp,2)
c      write(6,*) 'ason4pi**2*F2(pp,3)',ason4pi**2*F2(pp,3)
c      write(6,*) 'ason4pi**2*F2(pp,4)',ason4pi**2*F2(pp,4)

c      write(6,*) '   ason4pi*Ba1(0)',ason4pi*Ba1(0)
c      write(6,*) '   ason4pi*Ba1(1)',ason4pi*Ba1(1)
c      write(6,*) '   ason4pi*Ba1(2)',ason4pi*Ba1(2)

c      write(6,*) 'ason4pi**2*Ba2(0)',ason4pi**2*Ba2(0)
c      write(6,*) 'ason4pi**2*Ba2(1)',ason4pi**2*Ba2(1)
c      write(6,*) 'ason4pi**2*Ba2(2)',ason4pi**2*Ba2(2)
c      write(6,*) 'ason4pi**2*Ba2(3)',ason4pi**2*Ba2(3)
c      write(6,*) 'ason4pi**2*Ba2(4)',ason4pi**2*Ba2(4)

c      write(6,*) '   ason4pi*Bb1(0)',ason4pi*Bb1(0)
c      write(6,*) '   ason4pi*Bb1(1)',ason4pi*Bb1(1)

c      write(6,*) 'ason4pi**2*Bb2(0)',ason4pi**2*Bb2(0)
c      write(6,*) 'ason4pi**2*Bb2(1)',ason4pi**2*Bb2(1)
c      write(6,*) 'ason4pi**2*Bb2(2)',ason4pi**2*Bb2(2)
c      write(6,*) 'ason4pi**2*Bb2(3)',ason4pi**2*Bb2(3)
c      write(6,*) 'ason4pi**2*Bb2(4)',ason4pi**2*Bb2(4)

c      write(6,*) 'ason4pi**2*full2(4)',ason4pi**2*full2(4)
c      write(6,*) 'ason4pi**2*full2(3)',ason4pi**2*full2(3)
c      write(6,*) 'ason4pi**2*full2(2)',ason4pi**2*full2(2)
c      write(6,*) 'ason4pi**2*full2(1)',ason4pi**2*full2(1)
c      write(6,*) 'ason4pi**2*full2(0)',ason4pi**2*full2(0)
c      stop

c      write(6,*) 'GLY: full2(4)',full2(4)
c      write(6,*) 'GLY:-full2(3)',-full2(3)
c      write(6,*) 'GLY: full2(2)',full2(2)
c      write(6,*) 'GLY:-full2(1)',-full2(1)
c      write(6,*) 'GLY: full2(0)',full2(0)

      GLYassemble=GLYassemble+ason4pi**2
     & *(full2(0)+Lqt2cut*(-full2(1)+Lqt2cut*(full2(2)
     & +Lqt2cut*(-full2(3)+Lqt2cut*full2(4)))))
      endif

      return
      end
