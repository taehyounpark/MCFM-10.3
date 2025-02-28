!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_Z_only(p,wt,wt2,switch)
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
      include 'ewcorr.f'
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: yrap,pt,yraptwo,pttwo,r
c---  Z->e+e-(31) or b bbar(33): both measured, rapidities and momenta of 3 and 4 can
c---  be calculated, also the invariant mass m34
      real(dp):: y3,y4,y5,y34,pt3,pt4,pt5,pt34,m34,r35
      real(dp):: costheta,p3(4),p4(4),p34(4),wt_ew
      integer:: switch,n,nplotmax
      integer tag
      include 'first.f'
      common/nplotmax/nplotmax

c***********************************************************************
c                                                                      *
c     INITIAL BOOKKEEPING                                              *
c                                                                      *
c***********************************************************************
      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        y3=1d3
        y4=1d3
        y5=1d3
        y34=1d3
        pt3=0._dp
        pt4=0._dp
        pt5=1d3
        pt34=0._dp
        m34=0._dp
        r35=1d3
        jets=1
        costheta=0._dp
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

      y3=yrap(3,p)
      y4=yrap(4,p)
      y34=yraptwo(3,4,p)
      pt3=pt(3,p)
      pt4=pt(4,p)
      pt34=pttwo(3,4,p)
      m34=sqrt((p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &         -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2)
      p3(:)=p(3,:)
      p4(:)=p(4,:)
      p34(:)=p(3,:)+p(4,:)
!    compute lepton asymmetry as a function of m34
!    (see for example Eq.(3) of PLB718 (2013) 752)
      costheta=p34(3)/abs(p34(3))
     & *((p3(4)+p3(3))*(p4(4)-p4(3))-(p3(4)-p3(3))*(p4(4)+p4(3)))
     & /m34/sqrt(m34**2+p34(1)**2+p34(2)**2)

      if(jets > 0) then
         pt5=pt(5,p)
         y5=yrap(5,p)
         r35=R(p,3,5)
      else
         pt5=-1._dp
         y5=1d3
         r35=1d3
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

      call bookplot(n,tag,'m34',m34,wt,wt2,0._dp,8000._dp,100._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt3',pt3,wt,wt2,0._dp,4000._dp,50._dp,'lin')
      n=n+1
c-- for EW corrections, this is the histogram without them
      if (kewcorr /= knone) then
        call bookplot(n,tag,'m34 - no EW',m34,wt_noew,wt_noew**2,
     &   0._dp,8000._dp,100._dp,'lin')
        n=n+1
        wt_ew=wt_noew+wt
        call bookplot(n,tag,'m34 - with EW',m34,wt_ew,wt_ew**2,
     &   0._dp,8000._dp,100._dp,'lin')
        n=n+1
        call bookplot(n,tag,'m34 - +RELEW+',m34,wt_ew,wt_ew**2,
     &   0._dp,8000._dp,100._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt3 - no EW',pt3,wt_noew,wt_noew**2,
     &   0._dp,4000._dp,50._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt3 - with EW',pt3,wt_ew,wt_ew**2,
     &   0._dp,4000._dp,50._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt3 - +RELEW+',pt3,wt_ew,wt_ew**2,
     &   0._dp,4000._dp,50._dp,'lin')
        n=n+1
      else
        wt_ew=wt
      endif

      call bookplot(n,tag,'y3',y3,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y4',y4,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y34',y34,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt3',pt3,wt,wt2,0._dp,80._dp,2._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt4',pt4,wt,wt2,0._dp,80._dp,2._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt34',pt34,wt,wt2,0._dp,50._dp,2._dp,'lin')
      n=n+1
      call bookplot(n, tag,'m34',m34,wt,wt2,70._dp,110._dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DeltaR35',r35,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y5',y5,wt,wt2,-3.2_dp,3.2_dp,0.5_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,wt2,0._dp,100._dp,2._dp,'lin')
      n=n+1

!--- these histograms must be kept
      if ((costheta > 0._dp) .or. (tag == tagbook)) then
        call bookplot(n, tag,'m34 forward lepton',
     &   m34,wt_ew,wt_ew**2,40._dp,200._dp,5._dp,'lin')
      endif
      n=n+1
      if ((costheta <= 0._dp) .or. (tag == tagbook)) then
        call bookplot(n, tag,'m34 backward lepton',
     &   m34,wt_ew,wt_ew**2,40._dp,200._dp,5._dp,'lin')
      endif
      n=n+1
c--- placeholder for lepton FB asymmetry (will be filled properly in histofin)
      call bookplot(n, tag,'lepton +FB+ asymmetry',
     & m34,wt_ew,wt_ew**2,40._dp,200._dp,5._dp,'lin')
      n=n+1

c--- these histograms must be kept
      if ((costheta > 0._dp) .or. (tag == tagbook)) then
        call bookplot(n, tag,'m34 forward lepton',
     &   m34,wt_ew,wt_ew**2,200._dp,8000._dp,200._dp,'lin')
      endif
      n=n+1
      if ((costheta <= 0._dp) .or. (tag == tagbook)) then
        call bookplot(n, tag,'m34 backward lepton',
     &   m34,wt_ew,wt_ew**2,200._dp,8000._dp,200._dp,'lin')
      endif
      n=n+1
c--- placeholder for lepton FB asymmetry (will be filled properly in histofin)
      call bookplot(n, tag,'lepton +FB+ asymmetry',
     & m34,wt_ew,wt_ew**2,200._dp,8000._dp,200._dp,'lin')
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

