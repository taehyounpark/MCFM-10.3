!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module topwidth
          use types
          implicit none

      private

      public :: lotopdecaywidth
      public :: nloratiotopdecay
      public :: nnlotopdecay

      real(dp), private, save ::  mt1, besq, xi, ga
!$omp threadprivate(mt1,besq,xi,ga)

      real(dp), private, save :: ason2pi
!$omp threadprivate(ason2pi)

      contains

      function lotopdecaywidth(mt,mb,mw,gamw)
      implicit none
      include 'types.f'
      real(dp):: lotopdecaywidth
c***********************************************************************
c     Authors: R.K. Ellis and J. Campbell, February 2012               *
c                                                                      *
c     LO width of the top quark, including the effect of the           *
c     bottom quark mass and off shellness                              *
c                                                                      *
c     Formula is taken from Eq. (27) of                                *
c                                                                      *
c     \bibitem{Czarnecki:1990kv}                                       *
c     A.~Czarnecki,                                                    *
c     ``QCD corrections to the decay t ---> W b                        *
c       in dimensional regularization,''                               *
c     Phys.\ Lett.\ B {\bf 252}, 467 (1990).                           *
c     %%CITATION = PHLTA,B252,467;%%                                   *
c                                                                      *
c                                                                      *
c***********************************************************************

      include 'zerowidth.f'
      real(dp):: mb,mt,mw,om,omsq,be,dgauss,xlo,xhi,gamw

      real(dp), parameter :: tiny = 1d-8

      mt1=mt
      om=mw/mt
      be=mb/mt
      besq=be**2

      if (zerowidth) then
          omsq=om**2
          lotopdecaywidth=Gamma0(mt,besq,omsq)
      else
          xlo=0._dp
          xhi=(1._dp-be)**2
          ga=gamw/mw
          xi=(mt/mw)**2
          lotopdecaywidth=dgauss(Gamma0int,xlo,xhi,tiny)
      endif


      end function lotopdecaywidth

      function nloratiotopdecay(mt,mb,mw,gamw,scale)
          use constants
      implicit none
      real(dp):: nloratiotopdecay
c***********************************************************************
c     Authors: R.K. Ellis and J. Campbell, February 2012               *
c                                                                      *
c     ratio NLO/LO for the width of the top quark, including the       *
c     effect of the bottom quark mass.                                 *
c                                                                      *
c     Formula is taken from Eq. (27) of                                *
c                                                                      *
c     \bibitem{Czarnecki:1990kv}                                       *
c     A.~Czarnecki,                                                    *
c     ``QCD corrections to the decay t ---> W b                        *
c       in dimensional regularization,''                               *
c     Phys.\ Lett.\ B {\bf 252}, 467 (1990).                           *
c     %%CITATION = PHLTA,B252,467;%%                                   *
c                                                                      *
c     Formula has been improved in order to have a                     *
c     smooth mb -> 0 limit (in agreement with topwidth.f)              *
c                                                                      *
c***********************************************************************

c formula taken from Eq.(27) of
c---  %\cite{Czarnecki:1990kv}
c---  \bibitem{Czarnecki:1990kv}
c---  A.~Czarnecki,
c---  %``QCD corrections to the decay t ---> W b in dimensional regularization,''
c---  Phys.\ Lett.\ B {\bf 252}, 467 (1990).
c---- %%CITATION = PHLTA,B252,467;%%

      include 'zerowidth.f'
      include 'nlooprun.f'! nlooprun
      include 'couple.f'! amz
      real(dp), intent(in) :: mt,mb,mw,gamw,scale

      real(dp):: om,omsq,lo,ho,be,xlo,xhi

      real(dp), parameter :: tiny = 1d-8

      real(dp) :: alphas, dgauss

      ason2pi = alphas(scale,amz,nlooprun)/2._dp/pi

      mt1=mt
      om=mw/mt
      be=mb/mt
      ga=gamw/mw
      xi=(mt/mw)**2
      besq=be**2
      omsq=om**2

      if (zerowidth) then
        lo=Gamma0(mt,besq,omsq)
        ho=asGamma1(mt,besq,omsq)
      else
        xlo=0._dp
        xhi=(1._dp-be)**2
        lo=dgauss(Gamma0int,xlo,xhi,tiny)
        ho=dgauss(asGamma1int,xlo,xhi,tiny)
      endif

      nloratiotopdecay=ho/lo

      end function nloratiotopdecay


      function nnlotopdecay(mt,mw,scale)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'scet_const.f'
      include 'nlooprun.f'! nlooprun
      include 'couple.f'! amz
      real(dp):: nnlotopdecay
      real(dp):: w,logw,log2,XA,XH,XL,XNA
      real(dp) :: alphas
      real(dp), intent(in) :: mt,mw,scale
      include 'scet_beta.f'

c Formulae taken from section V of
c---  @article{Blokland:2005vq,
c---      author = "Blokland, Ian Richard and Czarnecki, Andrzej and Slusarczyk, Maciej and Tkachov, Fyodor",
c---      title = "{Next-to-next-to-leading order calculations for heavy-to-light decays}",
c---      eprint = "hep-ph/0503039",
c---      archivePrefix = "arXiv",
c---      reportNumber = "ALBERTA-THY-14-04, UVIC-TH-04-08, UVIC-TH-04-07",
c---      doi = "10.1103/PhysRevD.71.054004",
c---      journal = "Phys. Rev. D",
c---      volume = "71",
c---      pages = "054004",
c---      year = "2005",
c---      note = "[Erratum: Phys.Rev.D 79, 019901 (2009)]" }

      w=(mw/mt)**2
      logw=log(w)
      log2=log(2d0)

      XL = (-4d0/9d0 + 23d0*pisq/108d0 + zeta3 )
     & + w * ( -19d0/6d0 + 2d0*pisq/9d0 )
     & + w**2 * ( 745d0/72d0 - 31d0*pisq/36d0 - 3d0*zeta3 - 7d0/4d0*logw )
     & + w**3 * ( -5839d0/648d0 + 7d0*pisq/27d0 + 2d0*zeta3 + 5d0/3d0*logw )
     & + w**4 * ( 4253d0/8640d0 + pisq/4d0 - 17d0/144d0*logw )
     & + w**5 * ( -689d0/27000d0 + pisq/15d0 - 7d0/900d0*logw )
     & + w**6 * ( -13187d0/181440d0 + pisq/36d0 + 1d0/48d0*logw )
     & + w**7 * ( -2282381d0/37044000d0 + pisq/70d0 + 2263d0/88200d0*logw )

      XH = ( 12991d0/1296d0 - 53d0*pisq/54d0 - zeta3/3d0 )
     & + w * ( -35d0/108d0 - 4d0*pisq/9d0 + 4d0*zeta3 )
     & + w**2 * ( - 6377d0/432d0 + 25d0*pisq/18d0 + zeta3 )
     & + w**3 * ( 319d0/27d0 - 31d0*pisq/27d0 - 2d0*zeta3/3d0 )
     & + w**4 * ( 76873d0/8640d0 - 8d0*pisq/9d0 )
     & + w**5 * ( 237107d0/27000d0 - 8d0*pisq/9d0 )

      XA = ( 5d0 - 119d0*pisq/48d0 - 53d0*zeta3/8d0 - 11d0*pisq**2/720d0 + 19d0/4d0*pisq*log2 )
     & + w * ( -73d0/8d0 + 41d0*pisq/8d0 - 41d0*pisq**2/90d0 )
     & + w**2 * ( -7537d0/288d0 + 523d0*pisq/96d0 + 295d0*zeta3/32d0 - 191d0*pisq**2/720d0
     &            - 27d0/16d0*pisq*log2 + ( 115d0/48d0 - 5d0*pisq/16d0)*logw )
     & + w**3 * ( 16499d0/864d0 - 407d0*pisq/216d0 - 7d0*zeta3/2d0 + 7d0*pisq**2/120d0
     &            - pisq*log2 + ( -367d0/144d0 + 5d0*pisq/9d0 )*logw )
     & + w**4 * ( -1586479d0/259200d0 + 2951d0*pisq/6912d0 + 9d0*zeta3/2d0
     &            + ( 31979d0/17280d0 - pisq/16d0 )*logw )
     & + w**5 * ( -11808733d0/6480000d0 + 37d0*pisq/2400d0 + 6d0*zeta3/5d0
     &            + ( 13589d0/27000d0 - pisq/60d0 )*logw )

      XNA = ( 521d0/576d0 + 505d0*pisq/864d0 + 9d0*zeta3/16d0 + 11d0*pisq**2/1440d0
     &        - 19d0/8d0*pisq*log2 )
     & + w * ( 91d0/48d0 + 329d0*pisq/144d0 - 13d0*pisq**2/60d0 )
     & + w**2 * ( -12169d0/576d0 + 2171d0*pisq/576d0 + 377d0*zeta3/64d0 - 77d0*pisq**2/288d0
     &            + 27d0/32d0*pisq*log2 + ( 73d0/16d0 - 3d0*pisq/32d0 )*logw )
     & + w**3 * ( 13685d0/864d0 - 47d0*pisq/72d0 - 19d0*zeta3/4d0 + 43d0*pisq**2/720d0
     &            + 1d0/2d0*pisq*log2 + ( -1121d0/432d0 - pisq/6d0 )*logw )
     & + w**4 * ( -420749d0/103680d0 - 3263d0*pisq/13824d0 - 9d0*zeta3/8d0
     &            + ( 11941d0/6912d0 - 3d0*pisq/32d0 )*logw )
     & + w**5 * ( -4868261d0/12960000d0 - 557d0*pisq/4800d0 - 3d0*zeta3/10d0
     &            + ( 153397d0/216000d0 - pisq/40d0 )*logw )

c Change normalization to be coefficients of (as/2/pi)^2*Gamma0(mw -> 0)
      XL  = 4d0*CF * XL * TR
      XH  = 4d0*CF * XH * TR
      XA  = 4d0*CF * XA * CF
      XNA = 4d0*CF * XNA * CA

c Factor in front is no longer Gamma0(mw -> 0) but Gamma0(exact)
      XL  = XL / ((1d0-w)**2*(1d0+2d0*w))
      XH  = XH / ((1d0-w)**2*(1d0+2d0*w))
      XA  = XA / ((1d0-w)**2*(1d0+2d0*w))
      XNA = XNA / ((1d0-w)**2*(1d0+2d0*w))

      ason2pi = alphas(scale,amz,nlooprun)/2._dp/pi

c Remaining prefactor: Gamma0(exact)
      nnlotopdecay = ason2pi**2 * (XL*nf + XH + XA + XNA)

c Extra term to reinstate correct scale-dependence
      nnlotopdecay = nnlotopdecay
     & + be0 * ason2pi * nloratiotopdecay(mt, 0d0, mw, 0d0, scale) * log(scale/mt)

c      write(6,*) 'XL',XL
c      write(6,*) 'XH',XH
c      write(6,*) 'XA',XA
c      write(6,*) 'XNA',XNA

c Targets taken from 1301.713 equation (6.2) for mt = 172.85 and mw = 80.419
c      write(6,*) 'XL/target',XL/7.978d0
c      write(6,*) 'XH/target',XH/(-0.1166d0)
c      write(6,*) 'XA/target',XA/23.939d0
c      write(6,*) 'XNA/target',XNA/(-134.24d0)

      return
      end function nnlotopdecay

      function Gamma0(mt,besq,omsq)
      implicit none
      include 'types.f'
      real(dp):: Gamma0
c--   Author: John M. Campbell and R.K. Ellis, January 2012
c--   Taken from formula (2) of
c--   Fermilab-PUB-12-078-T

      include 'constants.f'
      include 'ewcouple.f'
      real(dp):: mt,omsq,besq,Gammainfty,f,P3b
      Gammainfty=GF*mt**3/(8._dp*rt2*pi)
      P3b=0.5_dp*sqrt(1._dp+omsq**2+besq**2-2._dp*(omsq+besq+omsq*besq))
      f=(1._dp-besq)**2+omsq*(1._dp+besq)-2._dp*omsq**2
      Gamma0=Gammainfty*2._dp*P3b*f
      return
      end function Gamma0

      function asGamma1(mt,besq,omsq)
      implicit none
      include 'types.f'
      real(dp):: asGamma1
c--   Author: John M. Campbell and R.K. Ellis, January 2012
c-----Taken from formula (5) of
c-----Fermilab-PUB-12-078-T

      include 'constants.f'
      include 'ewcouple.f'
      real(dp):: mt,P0,P3,PP,PM,W0,WP,wm,YW,z,omsq,om,be,
     & P0b,P3b,Pmb,PPb,Ywb,Wmb,f,besq,ddilog,GammaInfty,term4,
     & term7,term9

c     Statement functions.
      P0(z)=0.5_dp*(1._dp-omsq+z)
      P3(z)=0.5_dp*sqrt(1._dp+omsq**2+z**2-2._dp*(omsq+z+omsq*z))
      PP(z)=P0(z)+P3(z)
      PM(z)=P0(z)-P3(z)
      W0(z)=0.5_dp*(1._dp+omsq-z)
      WP(z)=W0(z)+P3(z)
      WM(z)=W0(z)-P3(z)
      YW(z)=0.5_dp*log(WP(z)/WM(z))
c      YP(z)=0.5_dp*log(PP(z)/PM(z))
c     End statement functions.

      f=(1._dp-besq)**2+omsq*(1._dp+besq)-2._dp*omsq**2
      om=sqrt(omsq)
      be=sqrt(besq)
      P0b=P0(besq)
      P3b=P3(besq)
      Pmb=PM(besq)
      Wmb=WM(besq)
c      WPb=WP(besq)
      PPb=PP(besq)
c      Ypb=YP(besq)
      Ywb=YW(besq)

      GammaInfty=Gf*mt**3/8._dp/pi/rt2
      if (besq > 0._dp) then
c      term4=Ypb*log(4._dp*P3b**2/Ppb**2/Wpb)
      term4=(log(PPb)-log(be))*log(4._dp*P3b**2*Wmb/(omsq*PPb**2))
      term7=
     & +(3._dp-besq+11._dp*besq**2-besq**3+omsq*(6._dp-12._dp*besq+2._dp*besq**2)
     & -omsq**2*(21._dp+5._dp*besq)+12._dp*omsq**3)*log(Ppb)
     & -(-besq+11._dp*besq**2-besq**3+omsq*(-12._dp*besq+2._dp*besq**2)
     & -omsq**2*(5._dp*besq))*log(be)
c     & -(3._dp+6._dp*omsq-omsq**2*21._dp+12._dp*omsq**3)*log(be)
      term9=
     & +6._dp*(1._dp-4._dp*besq+3._dp*besq**2+omsq*(3._dp+besq)-4._dp*omsq**2)
     & *(P3b-0.5_dp*(1._dp-omsq))*log(be)
     & +3._dp*(1._dp-omsq)*(-4._dp*besq+3._dp*besq**2+omsq*(besq))*log(be)
c     & +3._dp*(1._dp-omsq)*(1._dp+3._dp*omsq-4._dp*omsq**2)*log(be)
      else
      term4=log(PPb)*log(4._dp*P3b**2*Wmb/(omsq*PPb**2))
      term7=
     & +(3._dp-besq+11._dp*besq**2-besq**3+omsq*(6._dp-12._dp*besq+2._dp*besq**2)
     & -omsq**2*(21._dp+5._dp*besq)+12._dp*omsq**3)*log(Ppb)
      term9=0._dp
      endif


c--- equation for alphas*Gamma1
      asGamma1=GammaInfty*ason2pi*Cf*(
     & 8._dp*f*P0b*(ddilog(1._dp-Pmb)-ddilog(1._dp-Ppb)
     &  -2._dp*ddilog(1._dp-Pmb/Ppb)+term4
     &  +Ywb*log(Ppb))
     & +4._dp*(1._dp-besq)*((1._dp-besq)**2+omsq*(1._dp+besq)-4._dp*omsq**2)*Ywb
     & +term7
     & +8._dp*f*P3b*log(om/4._dp/P3b**2)+term9
     & +(5._dp-22._dp*besq+5._dp*besq**2+9._dp*omsq*(1._dp+besq)-6._dp*omsq**2)*P3b)

      return
      end function asGamma1

      function asGamma1int(omsq)
      implicit none
      include 'types.f'
      real(dp):: asGamma1int

c--   Author R.K. Ellis April 2012
c--   Integrand for NLO width with W-offshell
      include 'constants.f'
      real(dp):: omsq,asGamma1

      asGamma1int=ga*xi/pi
     & /((1._dp-xi*omsq)**2+ga**2)*asGamma1(mt1,besq,omsq)

      end function asGamma1int

      function Gamma0int(omsq)
      implicit none
      include 'types.f'
      real(dp):: Gamma0int

c--   Author R.K. Ellis April 2012
c--   Integrand for width with W-offshell
      include 'constants.f'
      real(dp):: omsq,Gamma0

      Gamma0int=(ga*xi/pi)/((1._dp-xi*omsq)**2+ga**2)*Gamma0(mt1,besq,omsq)

      end function Gamma0int

      end module

