!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_tbbdk(p,msq)
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: R.K. Ellis                                               *
c     January, 2012.                                                   *
c     calculate the LO element squared and subtraction terms           *
c     for the process                                                  *
c                                                                      *
c     u(-p1) +dbar(-p2)=t(nu(p3)+e+(p4)+b(p5))+bbar(p6)                *
c     or                                                               *
c     d(-p1) +ubar(-p2)=t~(e-(p3)+nu~(p4)+bb(p5))+b(p6)                *
c                                                                      *
c     Top (antitop) is kept strictly on-shell                          *
c     although all spin correlations are retained.                     *
c                                                                      *
c     NOTE: this routine is a replacement for qqb_tbb.f that includes  *
c           the effect of the b-quark mass. In the massless case it is *
c           approximately 10% slower than that routine                 *
c                                                                      *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'nwz.f'
      integer:: j,k,hb,hc,jmax,jmin
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,qqb,qbq
      complex(dp)::  prop
      complex(dp)::  mdecay(2,2),mprodqa(2,2),mprodaq(2,2),
     & mtotqa(2,2),mtotaq(2,2)
      parameter(jmax=1,jmin=2)

      prop=cplx2(zip,mt*twidth)
      fac=aveqq*xn**2*gwsq**4/abs(prop)**2

c----set all elements to zero
      msq(:,:)=0._dp

      if (nwz == +1) then
        call schantoponshell(1,2,p,0,mprodqa)
        call schantoponshell(2,1,p,0,mprodaq)
        call tdecay(p,3,4,5,mdecay)
      elseif (nwz == -1) then
        call schanatoponshell(1,2,p,0,mprodqa)
        call schanatoponshell(2,1,p,0,mprodaq)
        call adecay(p,3,4,5,mdecay)
      endif


c--- Calculate complete amplitudes and square
      qqb=0._dp
      qbq=0._dp
      do hb=1,2
      do hc=1,2
      mtotqa(hb,hc)=czip
      mtotaq(hb,hc)=czip
      if (nwz == +1) then
          do j=1,jmax
          mtotqa(hb,hc)=mtotqa(hb,hc)+mdecay(hb,j)*mprodqa(j,hc)
          mtotaq(hb,hc)=mtotaq(hb,hc)+mdecay(hb,j)*mprodaq(j,hc)
          enddo
      elseif (nwz == -1) then
          do j=jmin,2
          mtotqa(hb,hc)=mtotqa(hb,hc)+mprodqa(hb,j)*mdecay(j,hc)
          mtotaq(hb,hc)=mtotaq(hb,hc)+mprodaq(hb,j)*mdecay(j,hc)
          enddo
      endif
      qqb=qqb+abs(mtotqa(hb,hc))**2
      qbq=qbq+abs(mtotaq(hb,hc))**2
      enddo
      enddo

c---fill qb-q and q-qb elements
      do j=-nf,nf
      do k=-nf,nf
      if     ((j > 0) .and. (k < 0)) then
      msq(j,k)=Vsq(j,k)*fac*qqb
      elseif ((j < 0) .and. (k > 0)) then
      msq(j,k)=Vsq(j,k)*fac*qbq
      endif
      enddo
      enddo

      return
      end
