!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine powcorr_qa(order,tauc,x1,x2,Q,beama,beamb,powc)
c--- Implementation of power corrections for DY processes according to:
c---
c---   M.~A.~Ebert, I.~Moult, I.~W.~Stewart, F.~J.~Tackmann, G.~Vita and H.~X.~Zhu,
c---   ``Power Corrections for N-Jettiness Subtractions at ${\cal O}(\alpha_s)$,''
c---   JHEP {\bf 1812}, 084 (2018),  doi:10.1007/JHEP12(2018)084
c---   [arXiv:1807.10764 [hep-ph]].
c---
c---   Checked NLO results for qa and qg+gq against
c---    Fig. 5 (tauboost = F) and Fig. 3 (tauboost = T)
c---
c--- Extension to NNLO (for keepSubL = .false. only) provided by:
c---
c---   I.~Moult, L.~Rothen, I.~W.~Stewart, F.~J.~Tackmann and H.~X.~Zhu,
c---   ``Subleading Power Corrections for N-Jettiness Subtractions,''
c---   Phys.\ Rev.\ D {\bf 95}, no. 7, 074023 (2017), doi:10.1103/PhysRevD.95.074023
c---   [arXiv:1612.00450 [hep-ph]].
c---
c---   Checked NLO and NNLO results for qa and qg+gq, with tauboost = T,
c---    against Figs. 8 and 9 of this reference
c---
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'kpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'taucut.f'
      include 'beamtype.f'
      integer order,j,k
      real(dp):: x1,x2,Q,beama(-5:5),beamb(-5:5),powc(-nf:nf,-nf:nf),tauc
      real(dp):: dxpdfa(-5:5),dxpdfb(-5:5),err,dx,loga,logb,log3a,log3b,eYonrho
      real(dp):: Cqa1,Cdqa1,Cqda1,Cqg1,Cgq1,Cqa2,Cdqa2,Cqda2,Cqg2,Cgq2
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

        Cqa1=four*CF*(loga+logb)
        Cdqa1=-four*CF*loga
        Cqda1=-four*CF*logb

        Cqg1=-two*TR*logb
        Cgq1=-two*TR*loga

      else

        Cqa1=zip; Cdqa1=zip; Cqda1=zip; Cqg1=zip; Cgq1=zip

      endif

c Note that the NNLO coefficients are computed only for keepSubL = F
      if (order == 2) then

        Cqa2=-16._dp*CF**2*(log3a+log3b)
        Cdqa2=16._dp*CF**2*log3a
        Cqda2=16._dp*CF**2*log3b

        Cqg2=four*TR*(CF+CA)*log3b
        Cgq2=four*TR*(CF+CA)*log3a

      else

        Cqa2=zip; Cdqa2=zip; Cqda2=zip; Cqg2=zip; Cgq2=zip

      endif

c Fill array of power corrections
      do j=-nf,nf
        if (j == 0) cycle
        powc(0,j)=ason4pi*(Cgq1+ason4pi*Cgq2)*tauc*beama(0)*beamb(j)
        powc(j,0)=ason4pi*(Cqg1+ason4pi*Cqg2)*tauc*beama(j)*beamb(0)
        do k=-nf,nf
          if (j*k >= 0) cycle
          powc(j,k)=ason4pi*tauc*(
     &               (Cqa1+ason4pi*Cqa2)*beama(j)*beamb(k)
     &              +(Cdqa1+ason4pi*Cdqa2)*x1*dxpdfa(j)*beamb(k)
     &              +(Cqda1+ason4pi*Cqda2)*x2*beama(j)*dxpdfb(k))
        enddo
      enddo

      return
      end
