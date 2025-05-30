!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_4ftwdk(p,wt,wt2,switch)
      implicit none
      include 'types.f'
c--- Variable passed in to this routine:

c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)

c---     wt:  weight of this event

c---    wt2:  weight^2 of this event

c--- switch:  an integer:: equal to 0 or 1, depending on the type of event
c---                0  --> lowest order, virtual or real radiation
c---                1  --> counterterm for real radiation

      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'nqcdjets.f'
      real(dp):: p(mxpart,4),wt,wt2,pt
      real(dp):: tiny,swap,costheta,ylight
      integer:: switch,n,nplotmax,j
      integer tag
      logical:: failed
      parameter(tiny=1.e-8_dp)
      include 'first.f'
      common/nplotmax/nplotmax

c***********************************************************************
c                                                                      *
c     INITIAL BOOKKEEPING                                              *
c                                                                      *
c***********************************************************************

      if (first) then
c--- Initialize histograms only
        tag=tagbook
        goto 99
      else
c--- Add event in histograms
        tag=tagplot
      endif

c***********************************************************************
c                                                                      *
c     DEFINITIONS OF QUANTITIES TO PLOT                                *
c                                                                      *
c***********************************************************************

c--- Jets have already been reordered so that bottom, anti-bottom
c--- quarks are in the correct positions; just have to order jets in
c--- positions 7 and 8 according to pt
      if (p(8,4) > tiny) then
        if (pt(8,p) > pt(7,p)) then
          do j=1,4
          swap=p(7,j)
          p(7,j)=p(8,j)
          p(8,j)=swap
          enddo
        endif
      endif

c***********************************************************************
c                                                                      *
c     FILL HISTOGRAMS                                                  *
c                                                                      *
c***********************************************************************

c--- Call histogram routines
   99 continue

c--- "n" will count the number of histograms
      n=nextnplot

c--- Syntax of "bookplot" routine is:

c---   call bookplot(n,tag,titlex,var,wt,wt2,xmin,xmax,dx,llplot)

c---        n:  internal number of histogram
c---      tag:  "book" to initialize histogram, "plot" to fill
c---   titlex:  title of histogram
c---      var:  value of quantity being plotted
c---       wt:  weight of this event (passed in)
c---      wt2:  weight of this event (passed in)
c---     xmin:  lowest value to bin
c---     xmax:  highest value to bin
c---       dx:  bin width
c---   llplot:  equal to "lin"/"log" for linear/log scale

c--- fill variables that require top reconstruction
      call singletopreconstruct(p,failed,costheta,ylight)

      if ((failed) .and. (first .eqv. .false.)) then
c--- make sure to increment by the number of plots in the 'else' section
      n=n+2

      else

c--- cos(theta*)
      call bookplot(n,tag,'cos(theta*)',costheta,wt,wt2,
     & -1._dp,1._dp,0.1_dp,'lin')
      n=n+1

      call bookplot(n,tag,'|y(light)|',abs(ylight),wt,wt2,
     & 0._dp,5._dp,0.2_dp,'lin')
      n=n+1

      endif

c--- single-particle plots
      do j=3,7
        if ((abs(p(j,4)) > tiny) .or. (first)) then
          call genplot1(p,j,tag,wt,wt2,n)
        else
          n=n+2
        endif
      enddo
c--- two-particle plots
      call genplot2(p,3,4,tag,wt,wt2,n)
      if (((abs(p(5,4)) > tiny) .and. (abs(p(6,4)) > tiny))
     &    .or. (first)) then
        call genplot2(p,5,6,tag,wt,wt2,n)
      else
        n=n+3
      endif
c--- three-particle plots
      if ((abs(p(5,4)) > tiny) .or. (first)) then
        call genplot3(p,3,4,5,tag,wt,wt2,n)
      else
        n=n+3
      endif

c--- additional plots that may be present at NLO
      if ((abs(p(8,4)) > tiny) .or. (first)) then
        call genplot1(p,8,tag,wt,wt2,n)
      else
        n=n+2
      endif

c***********************************************************************
c                                                                      *
c     FINAL BOOKKEEPING                                                *
c                                                                      *
c***********************************************************************

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1


c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif

      return
      end

      subroutine singletopreconstruct(p,failed,costheta,ylight)
      implicit none
      include 'types.f'
c--- Given the usual momentum array, try to reconstruct which
c--- set of momenta reconstruct the top antitop quark;
c--- given those assignments, the routine returns the following quantities:
c---
c---    costheta      the angle between the lepton and the light jet,
c---                    computed in the top rest frame
c---    ylight            the rapidity of the light jet
c---
c--- For real radiation events, the jet algorithm may result in
c--- events for which the invariant mass is not exactly equal to mt,
c--- despite the phase space producing tops exactly on-shell. In that
c--- case, allow |s-mt| <= 'toler' [GeV]
c---
c--- If the event does not contain exactly one b-jet and one light jet
c--- then failed is equal to .true.

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'plabel.f'
      logical:: failed
      integer:: i
      real(dp):: tiny,p(mxpart,4),
     & ptop(4),ptoprest(4),plep(4),plight(4),pleprest(4),plightrest(4),
     & mtcand,costheta,ylight,yrap,ptopalt(4),mtcandalt
      parameter (tiny=1.e-8_dp)
c      real(dp):: small,toler
c      parameter (small=1.e-4_dp,toler=5d9)

c--- default: reconstruction is okay
      failed=.false.

c--- note: this routine was written for the case 'notag=1' in chooser.f
c--- (the default is notag=0) and inclusive=F, with the result that
c--- only 2 jets are ever accepted in the event
      if (p(5,4) < tiny) then
c        write(6,*) 'No b-quark present'
      failed=.true.
      return
      endif

      if (p(7,4) < tiny) then
c        write(6,*) 'No light jet present'
      failed=.true.
      return
      endif

c--- this requires that two jets are b and light-jet
      do i=1,4
        ptop(i)=p(3,i)+p(4,i)+p(5,i)
        ptopalt(i)=p(3,i)+p(4,i)+p(5,i)+p(7,i)
        plight(i)=p(7,i)
      if (plabel(3) == 'el') then
        plep(i)=p(3,i)
      else
        plep(i)=p(4,i)
      endif
      enddo
      mtcand=sqrt(ptop(4)**2-ptop(1)**2-ptop(2)**2-ptop(3)**2)
      mtcandalt=sqrt(ptopalt(4)**2-ptopalt(1)**2
     &               -ptopalt(2)**2-ptopalt(3)**2)
c      write(6,*) 'Event accepted, top candidate mass=',mtcand,mtcandalt

c--- check to see if 345+7 is a better top candidate than 345:
c--- if so, assume light get radiated in decay, so no angle to construct
      if (abs(mtcandalt-mt) < abs(mtcand-mt)) then
c        write(6,*) 'Radiation in decay: ',mtcandalt,' vs ',mtcand
      failed=.true.
      return
      endif

      ptoprest(4)=mtcand
      do i=1,3
      ptoprest(i)=0._dp
      enddo

      call boostx(plep,ptop,ptoprest,pleprest)
      call boostx(plight,ptop,ptoprest,plightrest)

c--- cos(theta*)
      costheta=(pleprest(1)*plightrest(1)
     &         +pleprest(2)*plightrest(2)
     &         +pleprest(3)*plightrest(3))
     &        /sqrt(pleprest(1)**2+pleprest(2)**2+pleprest(3)**2)
     &        /sqrt(plightrest(1)**2+plightrest(2)**2+plightrest(3)**2)
c      write(6,*) 'cos(theta*)=',costheta

c--- rapidity of light jet
      ylight=yrap(7,p)

      return
      end

