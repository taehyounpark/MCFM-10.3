!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine powcorr_gg(order,tauc,x1,x2,Q,beama,beamb,powc)
c--- Implementation of power corrections for gg processes according to:
c---
c---   M.~A.~Ebert, I.~Moult, I.~W.~Stewart, F.~J.~Tackmann, G.~Vita and H.~X.~Zhu,
c---   ``Power Corrections for N-Jettiness Subtractions at ${\cal O}(\alpha_s)$,''
c---   JHEP {\bf 1812}, 084 (2018),  doi:10.1007/JHEP12(2018)084
c---   [arXiv:1807.10764 [hep-ph]].
c---
c---   NLO results, with both values of keepSubL, for gg channel
c---    checked against Fig. 1 of this reference (hadronic definition, tauboost = F)
c---
c---    also checked gg and qg+gq against Fig. 10 for (tauboost = F) and Fig. 8 (tauboost = T)
c---
c--- Extension to NNLO (for keepSubL = .false. only) provided by:
c---
c---   I.~Moult, L.~Rothen, I.~W.~Stewart, F.~J.~Tackmann and H.~X.~Zhu,
c---   ``N -jettiness subtractions for $gg\to H$ at subleading power,''
c---   Phys.\ Rev.\ D {\bf 97}, no. 1, 014013 (2018),  doi:10.1103/PhysRevD.97.014013
c---   [arXiv:1710.03227 [hep-ph]].
c---
c---   Checked NLO and NNLO results for gg and qg+gq, with tauboost = T,
c---    against Figs. 10 and 11 of this reference
c---
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'kpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'taucut.f'
      include 'beamtype.f'
      integer order,j
      real(dp):: x1,x2,Q,beama(-5:5),beamb(-5:5),powc(-nf:nf,-nf:nf),tauc
      real(dp):: dxpdfa(-5:5),dxpdfb(-5:5),err,dx,loga,logb,log3a,log3b,eYonrho
      real(dp):: Cgg1,Cdgg1,Cgdg1,Cqg1,Cgq1,Cgg2,Cdgg2,Cgdg2,Cqg2,Cgq2
      real(dp), parameter:: xfrac=0.05_dp
      logical, parameter:: keepSubL=.false.

      powc(:,:)=zip

c--- compute derivatives of pdfs
      dx=min(xfrac*x1,one-x1)      ! keep variation of x1 inside [0,1]
      call dxpdf_dfridr(ih1,1,x1,dx,err,dxpdfa)
      dx=min(xfrac*x2,one-x2)      ! keep variation of x2 inside [0,1]
      call dxpdf_dfridr(ih2,2,x2,dx,err,dxpdfb)

      if (tauboost) then
        eYonrho=one
      else
c note that e^Y is given by Y = log(x1/x2)/two = yraptwo(3,4,p)
        eYonrho=sqrt(x1/x2)
      endif

      loga=log(tauc/Q)
      logb=loga
      if ((keepSubL) .and. (order < 2)) then
        loga=loga+log(eYonrho)
        logb=logb-log(eYonrho)
      endif

      log3a=loga**3/eYonrho/Q
      log3b=logb**3*eYonrho/Q

      loga=loga/eYonrho/Q
      logb=logb*eYonrho/Q

c NLO coefficients
      if ( (order == 1) .or.
     &    ((order == 2) .and. (coeffonly .eqv. .false.)) ) then

        Cgg1=four*CA*(loga+logb)
        Cdgg1=-four*CA*loga
        Cgdg1=-four*CA*logb

        Cqg1=-two*CF*loga
        Cgq1=-two*CF*logb

      else

        Cgg1=zip; Cdgg1=zip; Cgdg1=zip; Cqg1=zip; Cgq1=zip

      endif

c Note that the NNLO coefficients are computed only for keepSubL = F
      if (order == 2) then

        Cgg2=-16._dp*CA**2*(log3a+log3b)
        Cdgg2=16._dp*CA**2*log3a
        Cgdg2=16._dp*CA**2*log3b

        Cqg2=four*CF*(CF+CA)*log3a
        Cgq2=four*CF*(CF+CA)*log3b

      else

        Cgg2=zip; Cdgg2=zip; Cgdg2=zip; Cqg2=zip; Cgq2=zip

      endif

c Fill array of power corrections
      powc(0,0)=ason4pi*tauc*(
     &            (Cgg1+ason4pi*Cgg2)*beama(0)*beamb(0)
     &           +(Cdgg1+ason4pi*Cdgg2)*x1*dxpdfa(0)*beamb(0)
     &           +(Cgdg1+ason4pi*Cgdg2)*x2*beama(0)*dxpdfb(0))
      do j=-nf,nf
        if (j == 0) cycle
        powc(0,j)=ason4pi*(Cgq1+ason4pi*Cgq2)*tauc*beama(0)*beamb(j)
        powc(j,0)=ason4pi*(Cqg1+ason4pi*Cqg2)*tauc*beama(j)*beamb(0)
      enddo

      return
      end
