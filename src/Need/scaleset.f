!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine scaleset(rscalestart,fscalestart,p)
          use Scalevar
          use ptveto, only: jetptveto
      implicit none
      include 'types.f'
c--- wrapper routine to set a dynamic scale; please refer to individual
c--- routines for exact definitions of the scales;
c--- upgraded 3/2016 to convert string to integer on first call, for speed
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'scale.f'
      include 'facscale.f'
      include 'dynamicscale.f'
      include 'frag.f'
      include 'nlooprun.f'
      include 'qcdcouple.f'
      include 'couple.f'
      include 'stopscales.f'
      include 'masses.f'
      include 'first.f'
      include 'mpicommon.f'
      include 'taucut.f'
      real(dp):: rscalestart,fscalestart,p(mxpart,4),mu0
      integer, parameter ::
     & khmass=1,kwmass=2,kzmass=3,kmt=4,km34=5,km345=6,km3456=7,
     & kMsqpt34sq=8,kMsqpt345sq=9,kMsqpt5sq=10,kMsqptj1sq=11,kMsqsumptjsq=12,
     & km34sqsumptjsq=13,kptphoton=14,kHT=15,kddis=16,kmVpmH=17,kshat=18,kptj1=19,
     & kfixed=20,kHTprime=21, km34sqpt34sq=22, kpt3pt4d2=23, kHTos=24, kMsqptgamm1sq=25,
     & kMsqpt3456sq=26,kptveto=27
      integer, save :: scaleindex
!$omp threadprivate(scaleindex)

      ! functions
      real(dp) :: massvec, pttwo

c Initialization and write-out
      if (first) then
!$omp master
        if (rank == 0) then
        write(6,*)
        write(6,*)'************** Dynamic scale choice ****************'
        write(6,*)'*                                                  *'
        write(6,*)'*                 RENORMALIZATION                  *'
        write(6,45) ' mu_ren  =',rscalestart,dynstring
        write(6,*)'*                                                  *'
        write(6,*)'*                  FACTORIZATION                   *'
        write(6,45) ' mu_fac  =',fscalestart,dynstring
        if (frag) then
        write(6,*)'*                                                  *'
        write(6,*)'*                  FRAGMENTATION                   *'
        write(6,45) ' mu_frag =',frag_scalestart,dynstring
        endif
        write(6,*)'*                                                  *'
        write(6,*)'****************************************************'
        endif
!$omp end master
        first=.false.
        if     ((dynstring == 'mh') .or. (dynstring == 'mH')
     &     .or. (dynstring == 'Mh') .or. (dynstring == 'MH')) then
          scaleindex=khmass
        elseif ((dynstring == 'mw') .or. (dynstring == 'mW')
     &     .or. (dynstring == 'Mw') .or. (dynstring == 'MW')) then
          scaleindex=kwmass
        elseif ((dynstring == 'mz') .or. (dynstring == 'mZ')
     &     .or. (dynstring == 'Mz') .or. (dynstring == 'MZ')) then
          scaleindex=kzmass
        elseif ((dynstring == 'mt') .or. (dynstring == 'mT')
     &     .or. (dynstring == 'Mt') .or. (dynstring == 'MT')) then
          scaleindex=kmt
        elseif (dynstring == 'm(34)') then
          scaleindex=km34
        elseif (dynstring == 'm(345)') then
          scaleindex=km345
        elseif (dynstring == 'm(3456)') then
          scaleindex=km3456
        elseif (dynstring == 'sqrt(M^2+pt34^2)') then
          scaleindex=kMsqpt34sq
        elseif (dynstring == 'sqrt(M^2+pt345^2)') then
          scaleindex=kMsqpt345sq
        elseif (dynstring == 'sqrt(M^2+pt3456^2)') then
          scaleindex=kMsqpt3456sq
        elseif (dynstring == 'sqrt(M^2+pt5^2)') then
          scaleindex=kMsqpt5sq
        elseif (dynstring == 'sqrt(M^2+ptj1^2)') then
          scaleindex=kMsqptj1sq
        elseif (dynstring == 'sqrt(M^2+sumptj^2)') then
          scaleindex=kMsqsumptjsq
        elseif (dynstring == 'sqrt(m(34)^2+pt34^2)') then
          scaleindex=km34sqpt34sq
        elseif (dynstring == 'sqrt(m(34)^2+sumptj^2)') then
          scaleindex=km34sqsumptjsq
        elseif (dynstring == 'pt(j1)') then
          scaleindex=kptj1
        elseif (dynstring == 'pt(photon)') then
          scaleindex=kptphoton
        elseif (dynstring == 'HT') then
          scaleindex=kHT
        elseif (dynstring == 'DDIS') then
          scaleindex=kddis
        elseif (dynstring == 'mV+mH') then
          scaleindex=kmVpmH
        elseif (dynstring == 's-hat') then
          scaleindex=kshat
        elseif (dynstring == 'fixed') then
          scaleindex=kfixed
        elseif (dynstring == '(pt3+pt4)/2') then
          scaleindex=kpt3pt4d2
         elseif (dynstring == 'HTprime') then
           scaleindex=kHTprime
         elseif (dynstring == 'HTos') then
           scaleindex=kHTos
         elseif (dynstring == 'sqrt(M^2+ptgamm1^2)') then
           scaleindex=kMsqptgamm1sq
         elseif (dynstring == 'pt(veto)') then
           scaleindex=kptveto
        else
          write(6,*) 'Dynamic scale choice not recognized'
          write(6,*) '   dynamicscale = ',dynstring
          stop
        endif
      endif

c Set scale
      if     (scaleindex == khmass) then
        mu0=hmass
      elseif (scaleindex == kwmass) then
        mu0=wmass
      elseif (scaleindex == kzmass) then
        mu0=zmass
      elseif (scaleindex == kmt) then
        mu0=mt
      elseif (scaleindex == km34) then
        call scaleset_m34(p,mu0)
      elseif (scaleindex == km345) then
        call scaleset_m345(p,mu0)
      elseif (scaleindex == km3456) then
        call scaleset_m3456(p,mu0)
      elseif (scaleindex == kMsqpt34sq) then
        call scaleset_Msqpt34sq(p,mu0)
      elseif (scaleindex == kMsqpt345sq) then
        call scaleset_Msqpt345sq(p,mu0)
      elseif (scaleindex == kMsqpt3456sq) then
        call scaleset_Msqpt3456sq(p,mu0)
      elseif (scaleindex == kMsqpt5sq) then
        call scaleset_Msqpt5sq(p,mu0)
      elseif (scaleindex == kMsqptj1sq) then
        call scaleset_Msqptj1sq(p,mu0)
      elseif (scaleindex == kMsqsumptjsq) then
        call scaleset_Msqsumptjsq(p,mu0)
      elseif (scaleindex == km34sqsumptjsq) then
        call scaleset_m34sqsumptjsq(p,mu0)
      elseif (scaleindex == km34sqpt34sq) then
        mu0 = max(2d0, sqrt(massvec(p(3,:)+p(4,:)) + pttwo(3,4,p)**2))
      elseif (scaleindex == kptj1) then
        call scaleset_ptj1(p,mu0)
      elseif (scaleindex == kptphoton) then
        call scaleset_ptphoton(p,mu0)
      elseif (scaleindex == kHT) then
        call scaleset_HT(p,mu0)
      elseif (scaleindex == kDDIS) then
        !call scaleset_ddis(p,mu0)
        ! all scale dependence is handled elsewhere
        ! initialize with something meaningful
        mu0 = mt
      elseif (scaleindex == kmVpmH) then
        call scaleset_mVpmH(p,mu0)
      elseif (scaleindex == kshat) then
        call scaleset_shat(p,mu0)
      elseif (scaleindex == kpt3pt4d2) then
        mu0 = scaleset_pt3pt4d2(p)
       elseif (scaleindex == kHTprime) then
         call scaleset_HTprime(p,mu0)
       elseif (scaleindex == kHTos) then
         call scaleset_HTos(p,mu0)
       elseif (scaleindex == kMsqptgamm1sq) then
         call scaleset_Msqptgamm1sq(p,mu0)
       elseif (scaleindex == kptveto) then
         mu0 = jetptveto
      elseif (scaleindex == kfixed) then
        ! just keep fixed scales, for debugging purposes
        mu0 = 1._dp
      else
        write(6,*) 'Dynamic scale choice not recognized'
        write(6,*) '   scaleindex = ',scaleindex
        stop
      endif

      frag_scale=frag_scalestart*mu0
      if  (frag_scale > 900._dp) frag_scale=900._dp
      if  (frag_scale < 1._dp) frag_scale=1._dp


      ! renormalization and factorization scales now set in Mod/mod_Scalevar.f90
      call usescales(rscalestart*mu0, fscalestart*mu0)

      return

 45   format(1x,'* ',a15,f6.2,' x ',a24,' *')

      end

