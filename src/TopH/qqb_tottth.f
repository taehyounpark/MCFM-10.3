!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_tottth(p,msq)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     December, 1999.                                                  *
c     calculate the element squared
c     for the process                                                  *
c----My notation                                                       *
c     q(-p1) +qbar(-p2)=t(p3)+t(p4)+h(p5)

c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'kpart.f'
      include 'masses.f'
      include 'msbarmasses.f'
      include 'msq_cs.f'
      include 'first.f'

      integer:: j
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),
     & wtqqb,wtgg1,wtgg2,wtgg0,massfrun,mt_eff,facqq,facgg,ytsq
      save mt_eff
!$omp threadprivate(mt_eff)
      if (first) then
c--- run mt to appropriate scale
        if (kpart==klord) then
          mt_eff=massfrun(mt_msbar,hmass,amz,1)
        else
          mt_eff=massfrun(mt_msbar,hmass,amz,2)
        endif
        first=.false.
      endif

      call ttqqHampsq(p,3,4,1,2,wtqqb)
      call ttggHdriver(p,wtgg1,wtgg2,wtgg0)

      ytsq=0.5d0*gwsq*(mt_eff/wmass)**2
      facqq=V/8d0*gsq**2*ytsq
      facgg=facqq*xn
c----set all elements to zero
      msq(:,:)=0._dp

c---fill qb-q, gg and q-qb elements
      do j=-nf,nf
      if (j < 0) then
          msq(j,-j)=aveqq*facqq*wtqqb
      elseif (j == 0) then
          msq(j,j)=avegg*facgg*(wtgg1+wtgg2+wtgg0)
          msq_cs(1,j,j)=avegg*facgg*wtgg1
          msq_cs(2,j,j)=avegg*facgg*wtgg2
          msq_cs(0,j,j)=avegg*facgg*wtgg0
      elseif (j > 0) then
          msq(j,-j)=aveqq*facqq*wtqqb
      endif
      enddo
      return
      end

