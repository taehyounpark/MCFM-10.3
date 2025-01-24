!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c===== T. Dennen, May 2014
c===== Initialise array of bubble integral values.
c===== For use with qqb->4gamma and qqb->2j2gamma
c===== coefficients by same author.
      subroutine Bubint_init(i1,i2,i3,i4,i5,i6,Bubint,ord)
      use loopI2_generic
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'scale.f'
      include 'scalarselect.f'
      integer:: i1,i2,i3,i4,i5,i6
      integer:: ord
      complex(dp):: Bubint(25)
      real(dp):: t

      t(i1,i2,i3) = s(i1,i2) + s(i2,i3) + s(i3,i1)


      Bubint(1) = loopI2(t(i1,i2,i3),0d0,0d0,musq,ord)
      Bubint(2) = loopI2(t(i1,i2,i4),0d0,0d0,musq,ord)
      Bubint(3) = loopI2(t(i1,i2,i5),0d0,0d0,musq,ord)
      Bubint(4) = loopI2(t(i1,i2,i6),0d0,0d0,musq,ord)
      Bubint(5) = loopI2(t(i1,i3,i4),0d0,0d0,musq,ord)
      Bubint(6) = loopI2(t(i1,i3,i5),0d0,0d0,musq,ord)
      Bubint(7) = loopI2(t(i1,i3,i6),0d0,0d0,musq,ord)
      Bubint(8) = loopI2(t(i1,i4,i5),0d0,0d0,musq,ord)
      Bubint(9) = loopI2(t(i1,i4,i6),0d0,0d0,musq,ord)
      Bubint(10) = loopI2(t(i1,i5,i6),0d0,0d0,musq,ord)
      Bubint(11) = loopI2(s(i1,i2),0d0,0d0,musq,ord)
      Bubint(12) = loopI2(s(i5,i6),0d0,0d0,musq,ord)
      Bubint(13) = loopI2(s(i4,i6),0d0,0d0,musq,ord)
      Bubint(14) = loopI2(s(i4,i5),0d0,0d0,musq,ord)
      Bubint(15) = loopI2(s(i3,i6),0d0,0d0,musq,ord)
      Bubint(16) = loopI2(s(i3,i5),0d0,0d0,musq,ord)
      Bubint(17) = loopI2(s(i3,i4),0d0,0d0,musq,ord)
      Bubint(18) = loopI2(s(i1,i3),0d0,0d0,musq,ord)
      Bubint(19) = loopI2(s(i2,i6),0d0,0d0,musq,ord)
      Bubint(20) = loopI2(s(i2,i5),0d0,0d0,musq,ord)
      Bubint(21) = loopI2(s(i2,i4),0d0,0d0,musq,ord)
      Bubint(22) = loopI2(s(i1,i4),0d0,0d0,musq,ord)
      Bubint(23) = loopI2(s(i2,i3),0d0,0d0,musq,ord)
      Bubint(24) = loopI2(s(i1,i5),0d0,0d0,musq,ord)
      Bubint(25) = loopI2(s(i1,i6),0d0,0d0,musq,ord)

      return
      end
