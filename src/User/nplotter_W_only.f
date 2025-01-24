!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_W_only(p,wt,wt2,switch)
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
      include 'jetlabel.f'
      include 'nwz.f'
      real(dp):: p(mxpart,4),wt,wt2,yrap,pt,r,yraptwo,etaraptwo,
     & y5,pt5,Re5,y34,eta34
      integer:: switch,n,nplotmax
      integer tag
      include 'first.f'
      common/nplotmax/nplotmax
      real(dp) :: yel, ptel

c***********************************************************************
c                                                                      *
c     INITIAL BOOKKEEPING                                              *
c                                                                      *
c***********************************************************************


      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        yel = 1d3
c--- If there is no NLO jet, these initial y5, pt5 will not pass the cut
        y5=1d3
        pt5=1d3
c--- If Re5 is not changed by the NLO value, it will be out of
c--- the plotting range
        Re5=1d3
        jets=1
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

c--- W rapidity and pseudorapidity
      y34=yraptwo(3,4,p)
      eta34=etaraptwo(3,4,p)

c--- For W+ processes plot e^+(4).  Otherwise (W-) plot e^-(3).
      if(nwz == +1) then
         yel = yrap(4,p)
         ptel = pt(4,p)
      else
         yel = yrap(3,p)
         ptel = pt(3,p)
      endif
c---      eventpart=4+jets
      if(jets > 0) then
         pt5=pt(5,p)
         y5=yrap(5,p)
         if(nwz == +1) then
            Re5=R(p,4,5)
         else
            Re5=R(p,3,5)
         endif
      else
        pt5=-1._dp
        y5=1d3
        Re5=1d3
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

       call bookplot(n,tag,'W rapidity',y34,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
       n=n+1
       call bookplot(n,tag,'W ps-rap',eta34,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
       n=n+1

       call bookplot(n,tag,'y_lep',yel,wt,wt2,-6._dp,6._dp,0.25_dp,'lin')
       n=n+1
       call bookplot(n,tag,'pt_lep',ptel,wt,wt2,0._dp,100._dp,2._dp,'lin')
       n=n+1

      call bookplot(n,tag,'DeltaRe5',Re5,wt,wt2,0._dp,5._dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y5',y5,wt,wt2,-3._dp,3._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,wt2,0._dp,80._dp,2._dp,'lin')
      n=n+1


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

