!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function qqbgWH_topamp(i1,i2,i3,i4,i5,za,zb,mt2)
c===== C.Williams July 2015
c====== amplitude for q(i1)^-+qb(i2)^+ + ell^-(i3)+ell^+(i4)+g(i5) + Higgs
c===== where the Higgs is radiated from the top loop
      use loopI2_generic
      use loopI3_generic
      implicit none
      include 'types.f'
      complex(dp):: qqbgWH_topamp
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'mxpart.f'
      include 'scale.f'
      include 'masses.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'scalarselect.f'
      integer:: i1,i2,i3,i4,i5
c====== basis integrals B0mt(1) = Bub(s12345,mt) B0mt(2)=Bub(s1234,mt)
      real(dp)::s12345,s1234,s134,s234,s34,xd
      complex(dp):: zab2,zab3
      complex(dp):: helpart,intfunc,fac,prop34
      complex(dp):: tri,bmH,bmgH
      real(dp):: mt2
      logical:: useeft_wh
      common/useeft_wh/useeft_wh
!$omp threadprivate(/useeft_wh/)

 !====== statement functions
      include 'cplx.h'

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
      zab3(i1,i2,i3,i5,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)
     &     +za(i1,i5)*zb(i5,i4)

c====== useful invariants
      s34=s(i3,i4)
      s134=s(i1,i3)+s(i3,i4)+s(i1,i4)
      s234=s(i2,i3)+s(i2,i4)+s(i3,i4)
      s1234=s134+s234+s(i1,i2)-s(i3,i4)
c====== s12345 is the Higgs mass
      s12345=s(i1,i5)+s(i2,i5)+s(i3,i5)+s(i4,i5)+s1234

      prop34=cplx1(s34)/cplx2(s34-wmass**2,wmass*wwidth)

      if(useeft_wh) then
      fac=-im*sqrt(gsq)*gwsq*prop34*as/(3._dp*pi*sqrt(vevsq))/rt2
      else
      fac=im*sqrt(gsq)*ason4pi*gwsq*prop34*mt2/sqrt(wmass**2-im*wmass*wwidth)*gw/rt2/2._dp
      endif

c=== Spin independent function this is just the form factor
c==== and can be written in terms of existing MCFM functions

 !     intfunc=-ason4pi*sqrt(esq/xw)/wmass
 !    &     *ffDDHK_ql(s12345,s1234,mt2)*s1234*(four*xn*mt2/s1234)

c===== helicity dependent part (includes 1/s1234 from gluon propagator)

      helpart=((s134*zab2(i3,i2,i4,i5)*zab3(i1,i2,i3,i4,i5)*zb(i4,i2)
     &     +s234*za(i1,i3)*zb(i5,i2)
     & *(za(i1,i2)*zb(i4,i1)*zb(i5,i2)-za(i2,i3)*zb(i4,i3)*zb(i5,i2)
     &  +s134*zb(i5,i4))))/(s1234*s134*s234*s34)

      xd=one/(s12345-s1234)

      if(useeft_wh) then
      intfunc=1._dp
      else
         tri=loopI3(s12345,s1234,0._dp,mt2,mt2,mt2,musq,0)
         bmgH=loopI2(s1234,mt2,mt2,musq,0)
         bmH=loopI2(s12345,mt2,mt2,musq,0)
         intfunc=8._dp*((0.5_dp*(one-4._dp*mt2*xd)*tri
     &        + s1234*xd**2*(bmgH-bmH)-xd))
      endif
c      write(6,*) 'tri coefficient ',0.5_dp*(one-4_dp*mt2*xd)*helpart*8
c      write(6,*) 'bub(mH) coefficient ',-s1234*xd**2*helpart*8_dp
c      write(6,*) 'bub(mHg) coefficient ',s1234*xd**2*helpart*8_dp
c      write(6,*) 'rational piece ',-xd*helpart*8_dp
c      write(6,*) 'total ',intfunc*helpart
c      stop

c===== total
      qqbgWH_topamp=helpart*intfunc*fac

      return
      end

