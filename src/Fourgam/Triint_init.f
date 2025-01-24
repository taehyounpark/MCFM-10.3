!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c===== T. Dennen, May 2014
c===== Initialise array of triangle integral values.
c===== For use with qqb->4gamma and qqb->2j2gamma
c===== coefficients by same author.
      subroutine Triint_init(i1,i2,i3,i4,i5,i6,Triint,ord)
      use loopI3_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'scalarselect.f'
      integer:: i1,i2,i3,i4,i5,i6
      integer:: ord
      complex(dp):: Triint(90)
      real(dp):: t

      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)

      Triint(:) = czip

      Triint(1) = loopI3(0d0,0d0,s(i1,i2),0d0,0d0,0d0,musq,ord)
      Triint(2) = loopI3(0d0,s(i1,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(3) = loopI3(s(i5,i6),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(4) = loopI3(0d0,s(i1,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(5) = loopI3(s(i4,i6),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(6) = loopI3(0d0,s(i1,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(7) = loopI3(s(i4,i5),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(8) = loopI3(s(i3,i6),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(9) = loopI3(0d0,s(i1,i3),0d0,0d0,0d0,0d0,musq,ord)
      Triint(10) = loopI3(s(i3,i5),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(11) = loopI3(s(i3,i4),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(12) = loopI3(s(i2,i6),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(13) = loopI3(s(i2,i5),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(14) = loopI3(s(i2,i4),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(15) = loopI3(s(i2,i3),0d0,0d0,0d0,0d0,0d0,musq,ord)
      Triint(16) = loopI3(0d0,s(i2,i3),t(i4,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(17) = loopI3(0d0,t(i2,i3,i4),s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(18) = loopI3(s(i1,i2),0d0,t(i4,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(19) = loopI3(s(i1,i2),t(i3,i4,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(20) = loopI3(t(i1,i2,i3),0d0,s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(21) = loopI3(t(i1,i2,i3),s(i4,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(22) = loopI3(s(i1,i2),t(i3,i4,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(23) = loopI3(t(i1,i2,i3),s(i4,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(24) = loopI3(0d0,t(i2,i3,i5),s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(25) = loopI3(s(i1,i2),t(i3,i5,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(26) = loopI3(0d0,t(i2,i3,i6),s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(27) = loopI3(0d0,s(i2,i4),t(i3,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(28) = loopI3(t(i1,i2,i4),0d0,s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(29) = loopI3(t(i1,i2,i4),s(i3,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(30) = loopI3(t(i1,i2,i4),s(i3,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(31) = loopI3(0d0,t(i2,i4,i5),s(i3,i6),0d0,0d0,0d0,musq,ord)
      Triint(32) = loopI3(0d0,t(i2,i4,i6),s(i3,i5),0d0,0d0,0d0,musq,ord)
      Triint(33) = loopI3(0d0,s(i2,i5),t(i3,i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(34) = loopI3(t(i1,i2,i5),0d0,s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(35) = loopI3(t(i1,i2,i5),s(i3,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(36) = loopI3(t(i1,i2,i5),s(i3,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(37) = loopI3(0d0,t(i2,i5,i6),s(i3,i4),0d0,0d0,0d0,musq,ord)
      Triint(38) = loopI3(0d0,s(i2,i6),t(i3,i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(39) = loopI3(t(i1,i2,i6),0d0,s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(40) = loopI3(t(i1,i2,i6),s(i3,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(41) = loopI3(t(i1,i2,i6),s(i3,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(42) = loopI3(s(i1,i3),0d0,t(i4,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(43) = loopI3(s(i1,i3),t(i2,i4,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(44) = loopI3(s(i1,i3),t(i2,i4,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(45) = loopI3(s(i1,i3),t(i2,i5,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(46) = loopI3(t(i1,i3,i4),0d0,s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(47) = loopI3(t(i1,i3,i4),s(i2,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(48) = loopI3(t(i1,i3,i4),s(i2,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(49) = loopI3(t(i1,i3,i5),0d0,s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(50) = loopI3(t(i1,i3,i5),s(i2,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(51) = loopI3(t(i1,i3,i5),s(i2,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(52) = loopI3(t(i1,i3,i6),0d0,s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(53) = loopI3(t(i1,i3,i6),s(i2,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(54) = loopI3(t(i1,i3,i6),s(i2,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(55) = loopI3(s(i1,i4),0d0,t(i3,i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(56) = loopI3(s(i1,i4),t(i2,i3,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(57) = loopI3(s(i1,i4),t(i2,i3,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(58) = loopI3(s(i1,i4),t(i2,i5,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(59) = loopI3(t(i1,i4,i5),0d0,s(i3,i6),0d0,0d0,0d0,musq,ord)
      Triint(60) = loopI3(t(i1,i4,i5),s(i2,i3),0d0,0d0,0d0,0d0,musq,ord)
      Triint(61) = loopI3(t(i1,i4,i5),s(i2,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(62) = loopI3(t(i1,i4,i6),0d0,s(i3,i5),0d0,0d0,0d0,musq,ord)
      Triint(63) = loopI3(t(i1,i4,i6),s(i2,i3),0d0,0d0,0d0,0d0,musq,ord)
      Triint(64) = loopI3(t(i1,i4,i6),s(i2,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(65) = loopI3(s(i1,i5),0d0,t(i3,i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(66) = loopI3(s(i1,i5),t(i2,i3,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(67) = loopI3(s(i1,i5),t(i2,i3,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(68) = loopI3(s(i1,i5),t(i2,i4,i6),0d0,0d0,0d0,0d0,musq,ord)
      Triint(69) = loopI3(t(i1,i5,i6),0d0,s(i3,i4),0d0,0d0,0d0,musq,ord)
      Triint(70) = loopI3(t(i1,i5,i6),s(i2,i3),0d0,0d0,0d0,0d0,musq,ord)
      Triint(71) = loopI3(t(i1,i5,i6),s(i2,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(72) = loopI3(s(i1,i6),0d0,t(i3,i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(73) = loopI3(s(i1,i6),t(i2,i3,i4),0d0,0d0,0d0,0d0,musq,ord)
      Triint(74) = loopI3(s(i1,i6),t(i2,i3,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(75) = loopI3(s(i1,i6),t(i2,i4,i5),0d0,0d0,0d0,0d0,musq,ord)
      Triint(76) = loopI3(s(i1,i2),s(i3,i4),s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(77) = loopI3(s(i1,i2),s(i3,i5),s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(78) = loopI3(s(i1,i2),s(i3,i6),s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(79) = loopI3(s(i1,i3),s(i2,i4),s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(80) = loopI3(s(i1,i3),s(i2,i5),s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(81) = loopI3(s(i1,i3),s(i2,i6),s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(82) = loopI3(s(i1,i4),s(i2,i3),s(i5,i6),0d0,0d0,0d0,musq,ord)
      Triint(83) = loopI3(s(i1,i4),s(i2,i5),s(i3,i6),0d0,0d0,0d0,musq,ord)
      Triint(84) = loopI3(s(i1,i4),s(i2,i6),s(i3,i5),0d0,0d0,0d0,musq,ord)
      Triint(85) = loopI3(s(i1,i5),s(i2,i3),s(i4,i6),0d0,0d0,0d0,musq,ord)
      Triint(86) = loopI3(s(i1,i5),s(i2,i4),s(i3,i6),0d0,0d0,0d0,musq,ord)
      Triint(87) = loopI3(s(i1,i5),s(i2,i6),s(i3,i4),0d0,0d0,0d0,musq,ord)
      Triint(88) = loopI3(s(i1,i6),s(i2,i3),s(i4,i5),0d0,0d0,0d0,musq,ord)
      Triint(89) = loopI3(s(i1,i6),s(i2,i4),s(i3,i5),0d0,0d0,0d0,musq,ord)
      Triint(90) = loopI3(s(i1,i6),s(i2,i5),s(i3,i4),0d0,0d0,0d0,musq,ord)

      return
      end
