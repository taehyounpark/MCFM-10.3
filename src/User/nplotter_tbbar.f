!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_tbbar(p,wt,wt2,switch)
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
      real(dp):: p(mxpart,4),wt,wt2,tiny
      integer:: switch,n,nplotmax,j
      integer tag
      include 'first.f'
      common/nplotmax/nplotmax
      parameter(tiny=1.e-8_dp)

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

c--- single-particle plots
      do j=3,6
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
      if ((abs(p(7,4)) > tiny) .or. (first)) then
        call genplot1(p,7,tag,wt,wt2,n)
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





