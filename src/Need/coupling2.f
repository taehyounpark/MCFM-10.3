!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine coupling2
          use types
          use Scalevar, only: doScalevar, maxscalevar
      implicit none
c--- this routine calculates alpha-s using the now-determined
c--- value of nflav and makes CKM matrix diagonal if necessary

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
c      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'scale.f'
      include 'verbose.f'
      include 'nlooprun.f'
c      include 'kprocess.f'
c      include 'ewinput.f'
      include 'nflav.f'
      include 'b0.f'
      include 'dynamicscale.f'
      include 'stopscales.f'
      include 'fourthgen.f'
      include 'couple.f'
      include 'kpart.f'
      include 'nproc.f'
      include 'taucut.f'
      include 'scalevar.f'
      include 'mpicommon.f'
      include 'blha.f'
      real(dp):: alphas,cmass,bmass
      real(dp):: Vud,Vus,Vub,Vcd,Vcs,Vcb
      common/cabib/Vud,Vus,Vub,
     &             Vcd,Vcs,Vcb
      common/qmass/cmass,bmass
c      common/em/aemmz
c      common/mykpart/mykpart
      logical, save :: first = .true.

c--- set up the beta-function
      b0=(xn*11._dp-2._dp*nflav)/6._dp

c--- initialize the pdf set
      nlooprun=0
      if (useblha==0) then
         call pdfwrap
      endif

      musq=scale**2
c--- set up masses used in running of alphas (alfamz.f)
      cmass=sqrt(mcsq)
      if (fourthgen) then
        continue ! normal b mass already set in chooser.f
      else
        bmass=sqrt(mbsq) ! use the mass specified in the input file
      endif

      if (nflav < 5) then
        bmass=1001._dp
      endif
      if (nflav < 4) then
        cmass=1000._dp
      endif

c--- set the number of loops to use in the running of alpha_s
c--- if it hasn't been set by pdfwrap already
      if (nlooprun == 0) then
        if (kpart==klord) then
          nlooprun=1
        else
          nlooprun=2
        endif
      endif

c--- initialize alpha_s
      as=alphas(abs(scale),amz,nlooprun)

      ason2pi=as/twopi
      ason4pi=as/fourpi
      gsq=fourpi*as

      if (useblha == 0) then
c--- if we're doing W + jets, automatically make the CKM matrix
c--- diagonal since we're not interested in these small effects
      if ((nproc == 22) .or. (nproc == 27)) then
c      if ((nproc == 11) .or. (nproc == 16) .or.
c     &    (nproc == 22) .or. (nproc == 27)) then
c        if ( ((nproc == 11) .or. (nproc == 16))
c     &        .and. usescet .and. abovecut ) then
c! for W+jet process in above-cut SCET calculation, leave as-is
c          continue
c        else
          Vud=1._dp
          Vus=0._dp
          Vub=0._dp
          Vcd=0._dp
          Vcs=1._dp
          Vcb=0._dp
c        endif
      endif

      if (first .and. rank == 0) then
      write(6,*)
      write(6,*) '***************** CKM mixing matrix ****************'
      write(6,*) '*                                                  *'
      write(6,47) Vud,Vus,Vub
      write(6,48) Vcd,Vcs,Vcb
      write(6,*) '****************************************************'
      if ((nproc == 22) .or. (nproc == 27)) then
c      if ((nproc == 11) .or. (nproc == 16) .or.
c     &    (nproc == 22) .or. (nproc == 27)) then
c        if ( ((nproc == 11) .or. (nproc == 16))
c     &        .and. usescet .and. abovecut ) then
c for W+jet process in above-cut SCET calculation, leave as-is
c          continue
c        else
          write(6,*) '* Forced to be diagonal for simplicity in W + 2j   *'
          write(6,*) '****************************************************'
c        endif
      endif
 47   format(' *      Vud=',g10.5,'Vus=',g10.5,'Vub=',g10.5,'  *')
 48   format(' *      Vcd=',g10.5,'Vcs=',g10.5,'Vcb=',g10.5,'  *')
      endif
      endif

c--- special write-out for stop+b case
      if ((verbose .and. first .and. rank == 0) .and. (initrenscale_L > 0._dp)) then
      write(6,*)
      write(6,*) '************* Strong coupling, alpha_s  ************'
      write(6,*) '*                                                  *'
      write(6,49) 'alpha_s (zmass)    ',amz
      write(6,49) 'alpha_s (hvy scale)',as_H
      write(6,49) 'alpha_s (lgt scale)',as_L
      write(6,50) ' (using ',nlooprun,'-loop running of alpha_s)'
      write(6,*) '****************************************************'
      else if (verbose .and. first .and. rank ==0 .and. scale > 0._dp) then
      write(6,*)
      write(6,*) '************* Strong coupling, alpha_s  ************'
      write(6,*) '*                                                  *'
      if (dynamicscale .eqv. .false.) then
      write(6,49) 'alpha_s (scale)',gsq/fourpi
      write(6,49) 'alpha_s (zmass)',amz
      else
      write(6,*) '*  Dynamic scale - alpha_s changed event-by-event  *'
      write(6,49) 'alpha_s (zmass)',amz
      endif
      write(6,50) ' (using ',nlooprun,'-loop running of alpha_s)'
      write(6,*) '****************************************************'
 49   format(' *  ',a20,f12.8,16x,'*')
 50   format(' *  ',6x,a8,i1,a25,8x,'*')
      endif

      if (verbose .and. first .and. rank == 0 .and. doscalevar) then
      write(6,*)
      write(6,*) '****************************************************'
      write(6,*) '*                                                  *'
      write(6,*) '*     Computing scale variation in histograms      *'
      write(6,51) maxscalevar
      write(6,*) '*                                                  *'
      write(6,*) '****************************************************'
 51   format(' *    using ',i1,'-point variation by a factor of two    *')
      endif

      first = .false.

      return
      end
