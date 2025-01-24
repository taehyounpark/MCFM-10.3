!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function nloratiotopdecay(mt,mb,mw,gamw)
      implicit none
      include 'types.f'
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
      real(dp):: mb,mt,mt1,mw,om,omsq,lo,ho,be,besq,
     & Gamma0,Gamma0int,asGamma1,asGamma1int,dgauss,
     & xlo,xhi,xi,gamw,ga
      real(dp):: cachemass,cacheratio,tiny
      data cachemass,cacheratio,tiny/0._dp,0._dp,1.e-8_dp/
      common/transfer/mt1,besq,xi,ga
      save cachemass,cacheratio,tiny
!$omp threadprivate(cachemass,cacheratio,tiny,/transfer/)
      external Gamma0int,asGamma1int

c--- check to see if result has already been computed
c     if (abs(mt*mw-mb-cachemass) < tiny) then
c       nloratiotopdecay=cacheratio
c       return
c     endif

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


c--- set-up caching variables
      cachemass=mt*mw-mb
      cacheratio=nloratiotopdecay

      return
      end

