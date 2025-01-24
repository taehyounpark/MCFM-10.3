!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c===== C.Williams July 15, routine to construct interference between
c===== top-loop mediated diagrams and tree level in VH +jet production
      subroutine A5NLO_VHtop(i1,i2,i3,i4,i5,za,zb,A5LO,A5NLOt)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'zprods_decl.f'
c-==== this routine is really just a wrapper
      integer:: i1,i2,i3,i4,i5
      complex(dp):: A5LO,A5NLOt
      complex(dp)::qqbgWH_topamp,Alo5_VHt
      external qqbgWH_topamp,Alo5_VHt
      real(dp):: mt2
      logical:: useeft_wh
      common/useeft_wh/useeft_wh
!$omp threadprivate(/useeft_wh/)

c      useeft_wh=.true. ! for comparison with vh@nnlo
      useeft_wh=.false. ! for best prediction

      mt2=mt**2
c      if(useeft_wh) then
c      mrun2=mt**2
c      else
c      mrun2=massfrun(mt,scale,amz,nloop)**2
c      endif

      A5LO=Alo5_VHt(i1,i2,i3,i4,i5,za,zb)
      A5NLOt=qqbgWH_topamp(i1,i2,i3,i4,i5,za,zb,mt2)
c===== add in bottom loops too
      if(useeft_wh) return
      mt2=mb**2
c      mrun2=massfrun(mb,scale,amz,nloop)**2
      A5NLOt=A5NLOt+qqbgWH_topamp(i1,i2,i3,i4,i5,za,zb,mt2)

      return
      end

      function Alo5_VHt(i1,i2,i3,i4,i5,za,zb)
      implicit none
      include 'types.f'
      complex(dp):: Alo5_VHt
c======= this is the tree-level amplitude again, defined for a positive
c======= helicity gluon (i5), only reason this is used is to make factors
c+====== extracted for interference more convienent compared to those
c======= defined in the _v.f routine

c===== i1 is the negative helicity quark, i3 is the negative helicity lepton
      include 'constants.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'masses.f'
      include 'qcdcouple.f'
      integer:: i1,i2,i3,i4,i5
      complex(dp):: prop125,prop34
      complex(dp):: fac,amp
      complex(dp)::zab2
      real(dp):: s34,s125
c=====begin statement function
      include 'cplx.h'
c=====end statement function

      zab2(i1,i2,i3,i4)=za(i1,i2)*zb(i2,i4)+za(i1,i3)*zb(i3,i4)

      s125=s(i1,i2)+s(i2,i5)+s(i1,i5)
      s34=s(i3,i4)

      prop125=cplx1(s125)/cplx2(s125-wmass**2,wmass*wwidth)
      prop34=cplx1(s34)/cplx2(s34-wmass**2,wmass*wwidth)


c====== overall pre-factor
      fac=prop34*prop125*sqrt(gsq*gwsq)*gwsq*sqrt(wmass**2-im*wmass*wwidth)*rt2*im

c====== amplitude
      amp=-((za(i1,i3)*zab2(i1,i2,i5,i4))
     &     /(s125*s34*za(i1,i5)*za(i2,i5)))

      Alo5_VHt=amp*fac
      return
      end



