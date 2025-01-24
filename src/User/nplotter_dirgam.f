!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_dirgam(p,wt,wt2,switch,nd)
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
c-- nd determines whether we are analysing a photon dipole and hence have
c---to rescale accordingly

      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'jetlabel.f'
      include 'frag.f'
      include 'phot_dip.f'
      real(dp):: p(mxpart,4),wt,wt2,wtj
      real(dp):: yrap,pt,r,mjj
      real(dp):: y3,y4,y5,pt3,pt4,pt5,r45,r35,r34,s34,m34
      integer:: switch,n,nplotmax
      integer tag
      integer:: nd
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
        pt3=1d3
        pt4=1d3
        m34=-1d3
c---Initialise jet
        y5=1d3
        pt5=1d3
c--- If re5 is not changed by the NLO value, it will be out of
c--- the plotting range
        r34=1d3
        r35=1d3
        r45=1d3
        mjj=1d7
        jets=2
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
c--- Photons order based on pt

      if(jets<=1) then
         pt3=pt(3,p)
         pt4=pt(4,p)
         y3=yrap(3,p)
         y4=yrap(4,p)
         r34=R(p,3,4)
         s34=2._dp*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &           -p(4,3)*p(3,3))
         pt5=1d7
         y5=1d7
         r45=1d7
         r35=1d7
      else
         pt3=pt(3,p)
         y3=yrap(3,p)
         if(pt(4,p)>pt(5,p)) then
            pt4=pt(4,p)
            y4=yrap(4,p)
            pt5=pt(5,p)
            y5=yrap(5,p)
            r34=R(p,3,4)
            r35=R(p,3,5)
            s34=2._dp*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &              -p(4,3)*p(3,3))
         else
            pt4=pt(5,p)
            y4=yrap(5,p)
            pt5=pt(4,p)
            y5=yrap(4,p)
            r34=R(p,3,5)
            r35=R(p,3,4)
            s34=2._dp*(p(5,4)*p(3,4)-p(5,1)*p(3,1)-p(5,2)*p(3,2)
     &              -p(5,3)*p(3,3))
         endif
         r45=R(p,4,5)
         mjj=2._dp*(p(4,4)*p(5,4)-p(4,1)*p(5,1)-p(4,2)*p(5,2)
     &           -p(4,3)*p(5,3))

         mjj=sqrt(mjj)
      endif

      m34=sqrt(s34)

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

      call bookplot(n,tag,'pt_phot',pt3,wt,wt2,0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'y_phot',y3,wt,wt2,-3._dp,3._dp,0.1_dp,'lin')
      n=n+1

c quantities only defined for 1 or more jets
      if(jets > 0) then
        wtj=wt
      else
        wtj=zip
      endif
      call bookplot(n,tag,'pt_j1',pt4,wtj,wtj**2,0._dp,200._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'y_j1',y4,wtj,wtj**2,-5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DeltaR34',r34,wtj,wtj**2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'m_gamj',m34,wtj,wtj**2,0._dp,300._dp,20._dp,'log')
      n=n+1

c quantities only defined for 2 or more jets
      if(jets > 1) then
        wtj=wt
      else
        wtj=zip
      endif
      call bookplot(n,tag,'y_j2',y5,wtj,wtj**2,-3._dp,3._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt_j2',pt5,wtj,wtj**2,0._dp,200._dp,2._dp,'log')
      n=n+1
      call bookplot(n,tag,'mjj',mjj,wtj,wtj**2,0._dp,150._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'DeltaR35',r35,wtj,wtj**2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DeltaR45',r45,wtj,wtj**2,0._dp,5._dp,0.1_dp,'lin')
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



