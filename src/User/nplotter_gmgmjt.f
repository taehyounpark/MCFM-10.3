!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_gmgmjt(p,wt,wt2,switch)
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
      real(dp):: p(mxpart,4),wt,wt2,pt,pord(mxpart,4),pt3,pt4,
     & pt5,pt6
      integer:: switch,n,nplotmax
      integer tag
      common/nplotmax/nplotmax
      include 'first.f'


c***********************************************************************
c                                                                      *
c     INITIAL BOOKKEEPING                                              *
c                                                                      *
c***********************************************************************

      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
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

c--- order photons by pt
      pt3=pt(3,p)
      pt4=pt(4,p)

      pord(:,:)=p(:,:)

      if (pt4 > pt3) then
        pord(3,:)=p(4,:)
        pord(4,:)=p(3,:)
      endif

c--- order jets by pt if there are two present
      if (abs(p(6,4)) > 1.e-8_dp) then
        pt5=pt(5,p)
        pt6=pt(6,p)
        if (pt6 > pt5) then
          pord(5,:)=p(6,:)
          pord(6,:)=p(5,:)
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

c--- usual plots for highest pt photon
      call autoplot1(pord,3,tag,wt,wt2,n)

c--- usual plots for 2n.e-_dphighest pt photon
      call autoplot1(pord,4,tag,wt,wt2,n)

c--- usual plots for highest pt jet
      call autoplot1(pord,5,tag,wt,wt2,n)

c--- additional plots that may be present at NLO
      if (abs(p(6,4)) > 1.e-8_dp) then
        call autoplot1(pord,6,tag,wt,wt2,n)
        call autoplot2(pord,56,5,6,tag,wt,wt2,n)
      else
        n=n+6
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
