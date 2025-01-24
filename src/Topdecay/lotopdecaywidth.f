!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

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
      real(dp):: mb,mt,mt1,mw,om,omsq,be,besq,dgauss,Gamma0,
     & xlo,xhi,xi,gamw,ga,Gamma0int
      real(dp):: cachemass,cachewidth,tiny
      data cachemass,cachewidth,tiny/0._dp,0._dp,1.e-8_dp/
      common/transfer/mt1,besq,xi,ga
      save cachemass,cachewidth,tiny
!$omp threadprivate(cachemass,cachewidth,tiny,/transfer/)
      external Gamma0int


c--- check to see if result has already been computed
c     if (abs(mt*mw-mb-cachemass) < tiny) then
c       lotopdecaywidth=cachewidth
c       return
c     endif

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

c--- set-up caching variables
      cachemass=mt*mw-mb
      cachewidth=lotopdecaywidth

      return
      end



