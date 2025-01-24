!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c===== T. Dennen, May 2014
c===== Initialise array of box integral values.
c===== For use with qqb->4gamma and qqb->2j2gamma
c===== coefficients by same author.
      subroutine Boxint_init(i1,i2,i3,i4,i5,i6,Boxint,ord)
      use loopI4_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'scalarselect.f'
      integer:: i1,i2,i3,i4,i5,i6
      integer:: ord
      complex(dp):: Boxint(195)
      real(dp):: t

      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)


      Boxint(:) = czip

      Boxint(1) =
     &  loopI4(0d0,0d0,0d0,t(i4,i5,i6),s(i1,i2),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

c      write(6,*) Boxint(1),musq,t(i4,i5,i6)
 !     pause
      Boxint(2) =
     &  loopI4(0d0,0d0,t(i3,i4,i5),0d0,s(i1,i2),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(3) =
     &  loopI4(0d0,t(i2,i3,i4),0d0,0d0,s(i5,i6),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(4) =
     &  loopI4(t(i1,i2,i3),0d0,0d0,0d0,s(i5,i6),s(i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(5) =
     &  loopI4(0d0,0d0,t(i3,i4,i6),0d0,s(i1,i2),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(6) =
     &  loopI4(0d0,t(i2,i3,i4),0d0,0d0,s(i5,i6),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(7) =
     &  loopI4(t(i1,i2,i3),0d0,0d0,0d0,s(i5,i6),s(i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(8) =
     &  loopI4(0d0,t(i2,i3,i5),0d0,0d0,s(i4,i6),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(9) =
     &  loopI4(t(i1,i2,i3),0d0,0d0,0d0,s(i4,i6),s(i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(10) =
     &  loopI4(0d0,0d0,t(i3,i5,i6),0d0,s(i1,i2),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(11) =
     &  loopI4(0d0,t(i2,i3,i5),0d0,0d0,s(i4,i6),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(12) =
     &  loopI4(0d0,t(i2,i3,i6),0d0,0d0,s(i4,i5),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(13) =
     &  loopI4(0d0,t(i2,i3,i6),0d0,0d0,s(i4,i5),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(14) =
     &  loopI4(0d0,0d0,0d0,t(i3,i5,i6),s(i1,i2),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(15) =
     &  loopI4(t(i1,i2,i4),0d0,0d0,0d0,s(i5,i6),s(i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(16) =
     &  loopI4(t(i1,i2,i4),0d0,0d0,0d0,s(i5,i6),s(i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(17) =
     &  loopI4(0d0,t(i2,i4,i5),0d0,0d0,s(i3,i6),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(18) =
     &  loopI4(t(i1,i2,i4),0d0,0d0,0d0,s(i3,i6),s(i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(19) =
     &  loopI4(0d0,0d0,t(i4,i5,i6),0d0,s(i1,i2),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(20) =
     &  loopI4(0d0,t(i2,i4,i5),0d0,0d0,s(i3,i6),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(21) =
     &  loopI4(0d0,t(i2,i4,i6),0d0,0d0,s(i3,i5),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(22) =
     &  loopI4(0d0,t(i2,i4,i6),0d0,0d0,s(i3,i5),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(23) =
     &  loopI4(0d0,0d0,0d0,t(i3,i4,i6),s(i1,i2),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(24) =
     &  loopI4(t(i1,i2,i5),0d0,0d0,0d0,s(i4,i6),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(25) =
     &  loopI4(t(i1,i2,i5),0d0,0d0,0d0,s(i4,i6),s(i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(26) =
     &  loopI4(t(i1,i2,i5),0d0,0d0,0d0,s(i3,i6),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(27) =
     &  loopI4(0d0,t(i2,i5,i6),0d0,0d0,s(i3,i4),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(28) =
     &  loopI4(0d0,t(i2,i5,i6),0d0,0d0,s(i3,i4),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(29) =
     &  loopI4(0d0,0d0,0d0,t(i3,i4,i5),s(i1,i2),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(30) =
     &  loopI4(t(i1,i2,i6),0d0,0d0,0d0,s(i4,i5),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(31) =
     &  loopI4(t(i1,i2,i6),0d0,0d0,0d0,s(i4,i5),s(i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(32) =
     &  loopI4(t(i1,i2,i6),0d0,0d0,0d0,s(i3,i5),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(33) =
     &  loopI4(0d0,0d0,0d0,t(i4,i5,i6),s(i1,i3),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(34) =
     &  loopI4(0d0,0d0,t(i2,i4,i5),0d0,s(i1,i3),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(35) =
     &  loopI4(0d0,0d0,t(i2,i4,i6),0d0,s(i1,i3),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(36) =
     &  loopI4(0d0,0d0,t(i2,i5,i6),0d0,s(i1,i3),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(37) =
     &  loopI4(t(i1,i3,i4),0d0,0d0,0d0,s(i5,i6),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(38) =
     &  loopI4(t(i1,i3,i4),0d0,0d0,0d0,s(i5,i6),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(39) =
     &  loopI4(0d0,t(i3,i4,i5),0d0,0d0,s(i2,i6),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(40) =
     &  loopI4(t(i1,i3,i4),0d0,0d0,0d0,s(i2,i6),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(41) =
     &  loopI4(0d0,t(i3,i4,i6),0d0,0d0,s(i2,i5),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(42) =
     &  loopI4(t(i1,i3,i5),0d0,0d0,0d0,s(i4,i6),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(43) =
     &  loopI4(t(i1,i3,i5),0d0,0d0,0d0,s(i4,i6),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(44) =
     &  loopI4(t(i1,i3,i5),0d0,0d0,0d0,s(i2,i6),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(45) =
     &  loopI4(0d0,t(i3,i5,i6),0d0,0d0,s(i2,i4),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(46) =
     &  loopI4(t(i1,i3,i6),0d0,0d0,0d0,s(i4,i5),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(47) =
     &  loopI4(t(i1,i3,i6),0d0,0d0,0d0,s(i4,i5),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(48) =
     &  loopI4(t(i1,i3,i6),0d0,0d0,0d0,s(i2,i5),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(49) =
     &  loopI4(0d0,0d0,t(i2,i3,i5),0d0,s(i1,i4),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(50) =
     &  loopI4(0d0,0d0,t(i2,i3,i6),0d0,s(i1,i4),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(51) =
     &  loopI4(t(i1,i4,i5),0d0,0d0,0d0,s(i3,i6),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(52) =
     &  loopI4(t(i1,i4,i5),0d0,0d0,0d0,s(i3,i6),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(53) =
     &  loopI4(t(i1,i4,i5),0d0,0d0,0d0,s(i2,i6),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(54) =
     &  loopI4(t(i1,i4,i6),0d0,0d0,0d0,s(i3,i5),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(55) =
     &  loopI4(t(i1,i4,i6),0d0,0d0,0d0,s(i3,i5),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(56) =
     &  loopI4(t(i1,i4,i6),0d0,0d0,0d0,s(i2,i5),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(57) =
     &  loopI4(0d0,0d0,t(i2,i3,i4),0d0,s(i1,i5),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(58) =
     &  loopI4(t(i1,i5,i6),0d0,0d0,0d0,s(i3,i4),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(59) =
     &  loopI4(t(i1,i5,i6),0d0,0d0,0d0,s(i3,i4),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(60) =
     &  loopI4(t(i1,i5,i6),0d0,0d0,0d0,s(i2,i4),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(61) =
     &  loopI4(0d0,s(i2,i3),0d0,s(i5,i6),t(i1,i2,i3),t(i2,i3,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(62) =
     &  loopI4(s(i1,i2),0d0,s(i4,i5),0d0,t(i1,i2,i3),t(i3,i4,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(63) =
     &  loopI4(s(i1,i2),0d0,s(i4,i6),0d0,t(i1,i2,i3),t(i3,i4,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(64) =
     &  loopI4(0d0,s(i2,i3),0d0,s(i4,i6),t(i1,i2,i3),t(i2,i3,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(65) =
     &  loopI4(s(i1,i2),0d0,s(i5,i6),0d0,t(i1,i2,i3),t(i3,i5,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(66) =
     &  loopI4(0d0,s(i2,i3),0d0,s(i4,i5),t(i1,i2,i3),t(i2,i3,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(67) =
     &  loopI4(0d0,s(i2,i4),0d0,s(i5,i6),t(i1,i2,i4),t(i2,i4,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(68) =
     &  loopI4(s(i1,i2),0d0,s(i3,i5),0d0,t(i1,i2,i4),t(i4,i3,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(69) =
     &  loopI4(s(i1,i2),0d0,s(i3,i6),0d0,t(i1,i2,i4),t(i4,i3,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(70) =
     &  loopI4(0d0,s(i2,i4),0d0,s(i3,i6),t(i1,i2,i4),t(i2,i4,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(71) =
     &  loopI4(0d0,s(i2,i4),0d0,s(i3,i5),t(i1,i2,i4),t(i2,i4,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(72) =
     &  loopI4(0d0,s(i2,i5),0d0,s(i4,i6),t(i1,i2,i5),t(i2,i5,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(73) =
     &  loopI4(s(i1,i2),0d0,s(i3,i4),0d0,t(i1,i2,i5),t(i5,i3,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(74) =
     &  loopI4(0d0,s(i2,i5),0d0,s(i3,i6),t(i1,i2,i5),t(i2,i5,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(75) =
     &  loopI4(0d0,s(i2,i5),0d0,s(i3,i4),t(i1,i2,i5),t(i2,i5,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(76) =
     &  loopI4(0d0,s(i2,i6),0d0,s(i4,i5),t(i1,i2,i6),t(i2,i6,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(77) =
     &  loopI4(0d0,s(i2,i6),0d0,s(i3,i5),t(i1,i2,i6),t(i2,i6,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(78) =
     &  loopI4(0d0,s(i2,i6),0d0,s(i3,i4),t(i1,i2,i6),t(i2,i6,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(79) =
     &  loopI4(s(i1,i3),0d0,s(i4,i5),0d0,t(i1,i3,i2),t(i2,i4,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(80) =
     &  loopI4(s(i1,i3),0d0,s(i4,i6),0d0,t(i1,i3,i2),t(i2,i4,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(81) =
     &  loopI4(s(i1,i3),0d0,s(i5,i6),0d0,t(i1,i3,i2),t(i2,i5,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(82) =
     &  loopI4(0d0,s(i3,i4),0d0,s(i5,i6),t(i1,i3,i4),t(i3,i4,i2),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(83) =
     &  loopI4(s(i1,i3),0d0,s(i2,i5),0d0,t(i1,i3,i4),t(i4,i2,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(84) =
     &  loopI4(s(i1,i3),0d0,s(i2,i6),0d0,t(i1,i3,i4),t(i4,i2,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(85) =
     &  loopI4(0d0,s(i3,i5),0d0,s(i4,i6),t(i1,i3,i5),t(i3,i5,i2),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(86) =
     &  loopI4(s(i1,i3),0d0,s(i2,i4),0d0,t(i1,i3,i5),t(i5,i2,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(87) =
     &  loopI4(0d0,s(i3,i6),0d0,s(i4,i5),t(i1,i3,i6),t(i3,i6,i2),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(88) =
     &  loopI4(s(i1,i4),0d0,s(i3,i5),0d0,t(i1,i4,i2),t(i2,i3,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(89) =
     &  loopI4(s(i1,i4),0d0,s(i3,i6),0d0,t(i1,i4,i2),t(i2,i3,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(90) =
     &  loopI4(s(i1,i4),0d0,s(i5,i6),0d0,t(i1,i4,i2),t(i2,i5,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(91) =
     &  loopI4(s(i1,i4),0d0,s(i2,i5),0d0,t(i1,i4,i3),t(i3,i2,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(92) =
     &  loopI4(s(i1,i4),0d0,s(i2,i6),0d0,t(i1,i4,i3),t(i3,i2,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(93) =
     &  loopI4(s(i1,i4),0d0,s(i2,i3),0d0,t(i1,i4,i5),t(i5,i2,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(94) =
     &  loopI4(s(i1,i5),0d0,s(i3,i4),0d0,t(i1,i5,i2),t(i2,i3,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(95) =
     &  loopI4(s(i1,i5),0d0,s(i3,i6),0d0,t(i1,i5,i2),t(i2,i3,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(96) =
     &  loopI4(s(i1,i5),0d0,s(i4,i6),0d0,t(i1,i5,i2),t(i2,i4,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(97) =
     &  loopI4(s(i1,i5),0d0,s(i2,i4),0d0,t(i1,i5,i3),t(i3,i2,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(98) =
     &  loopI4(s(i1,i5),0d0,s(i2,i6),0d0,t(i1,i5,i3),t(i3,i2,i6),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(99) =
     &  loopI4(s(i1,i5),0d0,s(i2,i3),0d0,t(i1,i5,i4),t(i4,i2,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(100) =
     &  loopI4(s(i1,i6),0d0,s(i3,i4),0d0,t(i1,i6,i2),t(i2,i3,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(101) =
     &  loopI4(s(i1,i6),0d0,s(i3,i5),0d0,t(i1,i6,i2),t(i2,i3,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(102) =
     &  loopI4(s(i1,i6),0d0,s(i4,i5),0d0,t(i1,i6,i2),t(i2,i4,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(103) =
     &  loopI4(s(i1,i6),0d0,s(i2,i4),0d0,t(i1,i6,i3),t(i3,i2,i4),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(104) =
     &  loopI4(s(i1,i6),0d0,s(i2,i5),0d0,t(i1,i6,i3),t(i3,i2,i5),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(105) =
     &  loopI4(s(i1,i6),0d0,s(i2,i3),0d0,t(i1,i6,i4),t(i4,i2,i3),
     & 0d0,0d0,0d0,0d0,musq,ord)

      Boxint(106) =
     &  loopI4(0d0,0d0,s(i3,i4),s(i5,i6),s(i1,i2),t(i2,i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(107) =
     &  loopI4(0d0,s(i2,i3),s(i4,i5),0d0,t(i1,i2,i3),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(108) =
     &  loopI4(s(i1,i2),0d0,0d0,s(i5,i6),t(i1,i2,i3),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(109) =
     &  loopI4(s(i1,i2),s(i3,i4),0d0,0d0,s(i5,i6),t(i3,i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(110) =
     &  loopI4(0d0,s(i2,i3),s(i4,i6),0d0,t(i1,i2,i3),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(111) =
     &  loopI4(s(i1,i2),s(i3,i4),0d0,0d0,s(i5,i6),t(i3,i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(112) =
     &  loopI4(0d0,0d0,s(i3,i5),s(i4,i6),s(i1,i2),t(i2,i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(113) =
     &  loopI4(s(i1,i2),0d0,0d0,s(i4,i6),t(i1,i2,i3),s(i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(114) =
     &  loopI4(s(i1,i2),s(i3,i5),0d0,0d0,s(i4,i6),t(i3,i5,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(115) =
     &  loopI4(0d0,s(i2,i3),s(i5,i6),0d0,t(i1,i2,i3),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(116) =
     &  loopI4(s(i1,i2),s(i3,i5),0d0,0d0,s(i4,i6),t(i3,i5,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(117) =
     &  loopI4(0d0,0d0,s(i3,i6),s(i4,i5),s(i1,i2),t(i2,i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(118) =
     &  loopI4(s(i1,i2),0d0,0d0,s(i4,i5),t(i1,i2,i3),s(i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(119) =
     &  loopI4(s(i1,i2),s(i3,i6),0d0,0d0,s(i4,i5),t(i3,i6,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(120) =
     &  loopI4(s(i1,i2),s(i3,i6),0d0,0d0,s(i4,i5),t(i3,i6,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(121) =
     &  loopI4(0d0,s(i2,i4),s(i3,i5),0d0,t(i1,i2,i4),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(122) =
     &  loopI4(s(i1,i2),0d0,0d0,s(i5,i6),t(i1,i2,i4),s(i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(123) =
     &  loopI4(0d0,s(i2,i4),s(i3,i6),0d0,t(i1,i2,i4),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(124) =
     &  loopI4(0d0,0d0,s(i4,i5),s(i3,i6),s(i1,i2),t(i2,i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(125) =
     &  loopI4(s(i1,i2),s(i4,i5),0d0,0d0,s(i3,i6),t(i4,i5,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(126) =
     &  loopI4(0d0,s(i2,i4),s(i5,i6),0d0,t(i1,i2,i4),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(127) =
     &  loopI4(0d0,0d0,s(i4,i6),s(i3,i5),s(i1,i2),t(i2,i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(128) =
     &  loopI4(s(i1,i2),s(i4,i6),0d0,0d0,s(i3,i5),t(i4,i6,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(129) =
     &  loopI4(0d0,s(i2,i5),s(i3,i4),0d0,t(i1,i2,i5),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(130) =
     &  loopI4(0d0,s(i2,i5),s(i3,i6),0d0,t(i1,i2,i5),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(131) =
     &  loopI4(0d0,s(i2,i5),s(i4,i6),0d0,t(i1,i2,i5),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(132) =
     &  loopI4(0d0,0d0,s(i5,i6),s(i3,i4),s(i1,i2),t(i2,i5,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(133) =
     &  loopI4(0d0,s(i2,i6),s(i3,i4),0d0,t(i1,i2,i6),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(134) =
     &  loopI4(0d0,s(i2,i6),s(i3,i5),0d0,t(i1,i2,i6),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(135) =
     &  loopI4(0d0,s(i2,i6),s(i4,i5),0d0,t(i1,i2,i6),s(i1,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(136) =
     &  loopI4(0d0,0d0,s(i2,i4),s(i5,i6),s(i1,i3),t(i3,i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(137) =
     &  loopI4(s(i1,i3),0d0,0d0,s(i5,i6),t(i1,i3,i2),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(138) =
     &  loopI4(s(i1,i3),s(i2,i4),0d0,0d0,s(i5,i6),t(i2,i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(139) =
     &  loopI4(s(i1,i3),s(i2,i4),0d0,0d0,s(i5,i6),t(i2,i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(140) =
     &  loopI4(0d0,0d0,s(i2,i5),s(i4,i6),s(i1,i3),t(i3,i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(141) =
     &  loopI4(s(i1,i3),0d0,0d0,s(i4,i6),t(i1,i3,i2),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(142) =
     &  loopI4(s(i1,i3),s(i2,i5),0d0,0d0,s(i4,i6),t(i2,i5,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(143) =
     &  loopI4(s(i1,i3),s(i2,i5),0d0,0d0,s(i4,i6),t(i2,i5,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(144) =
     &  loopI4(0d0,0d0,s(i2,i6),s(i4,i5),s(i1,i3),t(i3,i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(145) =
     &  loopI4(s(i1,i3),0d0,0d0,s(i4,i5),t(i1,i3,i2),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(146) =
     &  loopI4(s(i1,i3),s(i2,i6),0d0,0d0,s(i4,i5),t(i2,i6,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(147) =
     &  loopI4(s(i1,i3),s(i2,i6),0d0,0d0,s(i4,i5),t(i2,i6,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(148) =
     &  loopI4(0d0,s(i3,i4),s(i2,i5),0d0,t(i1,i3,i4),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(149) =
     &  loopI4(s(i1,i3),0d0,0d0,s(i5,i6),t(i1,i3,i4),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(150) =
     &  loopI4(0d0,s(i3,i4),s(i2,i6),0d0,t(i1,i3,i4),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(151) =
     &  loopI4(s(i1,i3),s(i4,i5),0d0,0d0,s(i2,i6),t(i4,i5,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(152) =
     &  loopI4(s(i1,i3),s(i4,i6),0d0,0d0,s(i2,i5),t(i4,i6,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(153) =
     &  loopI4(0d0,s(i3,i5),s(i2,i4),0d0,t(i1,i3,i5),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(154) =
     &  loopI4(0d0,s(i3,i5),s(i2,i6),0d0,t(i1,i3,i5),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(155) =
     &  loopI4(0d0,s(i3,i6),s(i2,i4),0d0,t(i1,i3,i6),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(156) =
     &  loopI4(0d0,s(i3,i6),s(i2,i5),0d0,t(i1,i3,i6),s(i1,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(157) =
     &  loopI4(0d0,0d0,s(i2,i3),s(i5,i6),s(i1,i4),t(i4,i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(158) =
     &  loopI4(s(i1,i4),0d0,0d0,s(i5,i6),t(i1,i4,i2),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(159) =
     &  loopI4(s(i1,i4),s(i2,i3),0d0,0d0,s(i5,i6),t(i2,i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(160) =
     &  loopI4(s(i1,i4),s(i2,i3),0d0,0d0,s(i5,i6),t(i2,i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(161) =
     &  loopI4(s(i1,i4),0d0,0d0,s(i3,i6),t(i1,i4,i2),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(162) =
     &  loopI4(s(i1,i4),s(i2,i5),0d0,0d0,s(i3,i6),t(i2,i5,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(163) =
     &  loopI4(s(i1,i4),s(i2,i5),0d0,0d0,s(i3,i6),t(i2,i5,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(164) =
     &  loopI4(s(i1,i4),0d0,0d0,s(i3,i5),t(i1,i4,i2),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(165) =
     &  loopI4(s(i1,i4),s(i2,i6),0d0,0d0,s(i3,i5),t(i2,i6,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(166) =
     &  loopI4(s(i1,i4),s(i2,i6),0d0,0d0,s(i3,i5),t(i2,i6,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(167) =
     &  loopI4(s(i1,i4),0d0,0d0,s(i5,i6),t(i1,i4,i3),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(168) =
     &  loopI4(s(i1,i4),s(i3,i5),0d0,0d0,s(i2,i6),t(i3,i5,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(169) =
     &  loopI4(s(i1,i4),s(i3,i6),0d0,0d0,s(i2,i5),t(i3,i6,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(170) =
     &  loopI4(0d0,s(i4,i5),s(i2,i3),0d0,t(i1,i4,i5),s(i1,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(171) =
     &  loopI4(0d0,s(i4,i6),s(i2,i3),0d0,t(i1,i4,i6),s(i1,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(172) =
     &  loopI4(s(i1,i5),0d0,0d0,s(i4,i6),t(i1,i5,i2),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(173) =
     &  loopI4(s(i1,i5),s(i2,i3),0d0,0d0,s(i4,i6),t(i2,i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(174) =
     &  loopI4(s(i1,i5),s(i2,i3),0d0,0d0,s(i4,i6),t(i2,i3,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(175) =
     &  loopI4(s(i1,i5),0d0,0d0,s(i3,i6),t(i1,i5,i2),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(176) =
     &  loopI4(s(i1,i5),s(i2,i4),0d0,0d0,s(i3,i6),t(i2,i4,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(177) =
     &  loopI4(s(i1,i5),s(i2,i4),0d0,0d0,s(i3,i6),t(i2,i4,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(178) =
     &  loopI4(s(i1,i5),0d0,0d0,s(i3,i4),t(i1,i5,i2),s(i2,i6),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(179) =
     &  loopI4(s(i1,i5),s(i2,i6),0d0,0d0,s(i3,i4),t(i2,i6,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(180) =
     &  loopI4(s(i1,i5),s(i2,i6),0d0,0d0,s(i3,i4),t(i2,i6,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(181) =
     &  loopI4(s(i1,i5),0d0,0d0,s(i4,i6),t(i1,i5,i3),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(182) =
     &  loopI4(s(i1,i5),s(i3,i4),0d0,0d0,s(i2,i6),t(i3,i4,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(183) =
     &  loopI4(s(i1,i5),s(i3,i6),0d0,0d0,s(i2,i4),t(i3,i6,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(184) =
     &  loopI4(s(i1,i6),0d0,0d0,s(i4,i5),t(i1,i6,i2),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(185) =
     &  loopI4(s(i1,i6),s(i2,i3),0d0,0d0,s(i4,i5),t(i2,i3,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(186) =
     &  loopI4(s(i1,i6),s(i2,i3),0d0,0d0,s(i4,i5),t(i2,i3,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(187) =
     &  loopI4(s(i1,i6),0d0,0d0,s(i3,i5),t(i1,i6,i2),s(i2,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(188) =
     &  loopI4(s(i1,i6),s(i2,i4),0d0,0d0,s(i3,i5),t(i2,i4,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(189) =
     &  loopI4(s(i1,i6),s(i2,i4),0d0,0d0,s(i3,i5),t(i2,i4,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(190) =
     &  loopI4(s(i1,i6),0d0,0d0,s(i3,i4),t(i1,i6,i2),s(i2,i5),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(191) =
     &  loopI4(s(i1,i6),s(i2,i5),0d0,0d0,s(i3,i4),t(i2,i5,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(192) =
     &  loopI4(s(i1,i6),s(i2,i5),0d0,0d0,s(i3,i4),t(i2,i5,i4),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(193) =
     &  loopI4(s(i1,i6),0d0,0d0,s(i4,i5),t(i1,i6,i3),s(i2,i3),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(194) =
     &  loopI4(s(i1,i6),s(i3,i4),0d0,0d0,s(i2,i5),t(i3,i4,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)

      Boxint(195) =
     &  loopI4(s(i1,i6),s(i3,i5),0d0,0d0,s(i2,i4),t(i3,i5,i2),0d0,
     & 0d0,0d0,0d0,musq,ord)


      return
      end
