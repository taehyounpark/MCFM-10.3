!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_Hjets(p,wt,wt2,switch)
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
      include 'histo.f'
      include 'jetlabel.f'
      include 'nproc.f'
      include 'first.f'
      integer:: i5,i6,i7,nu
      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: yrap,pt,
     & yj1,ptj1,yj2,ptj2,yj3,ptj3,getet,y3,y4,pt3,pt4,
     & pt5,pt6,pt7,tmp5(4),tmp6(4),tmp7(4),oldpt(5:7),
     & sumy,y34,eta34,pt34,yraptwo,etaraptwo,pttwo,ay345,yrapthree
      integer:: switch,n,nplotmax
      integer tag
      common/nplotmax/nplotmax

c***********************************************************************
c                                                                      *
c     INITIAL BOOKKEEPING                                              *
c                                                                      *
c***********************************************************************

c--- Initialize dummy values for all quantities that could be plotted
      ptj1=-1._dp
      ptj2=-1._dp
      ptj3=-1._dp
      yj1=99._dp
      yj2=99._dp
      yj3=99._dp
      y34=99._dp
      eta34=99._dp
      ay345=-1._dp

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

c--- H rapidity and pseudorapidity
      y34=yraptwo(3,4,p)
      eta34=etaraptwo(3,4,p)
      pt34=pttwo(3,4,p)

c--- decay products pt and rapidity
      y3=yrap(3,p)
      pt3=pt(3,p)
      y4=yrap(4,p)
      pt4=pt(4,p)

c--- BEGIN: order jets according to pt
      if (jets > 0) then
      pt5=getet(p(5,4),p(5,1),p(5,2),p(5,3))
      if (jets > 1) pt6=getet(p(6,4),p(6,1),p(6,2),p(6,3))
      if (jets > 2) pt7=getet(p(7,4),p(7,1),p(7,2),p(7,3))
      i5=5
      i6=6
      i7=7
      oldpt(5)=pt5
      if (jets > 1) oldpt(6)=pt6
      if (jets > 2) oldpt(7)=pt7
c--- sort for 2 jets
      if (jets == 2) then
        if (pt6 > pt5) then
          i5=6
          i6=5
        endif
      endif
c--- sort for 3 jets
      if (jets == 3) then
        if ((pt5 > pt6) .and. (pt5 > pt7)) then
           i5=5
          if (pt6 > pt7) then
            i6=6
            i7=7
          else
            i6=7
            i7=6
          endif
        endif
        if ((pt6 > pt5) .and. (pt6 > pt7)) then
           i5=6
          if (pt5 > pt7) then
            i6=5
            i7=7
          else
            i6=7
            i7=5
          endif
        endif
        if ((pt7 > pt5) .and. (pt7 > pt6)) then
           i5=7
          if (pt5 > pt6) then
            i6=5
            i7=6
          else
            i6=6
            i7=5
          endif
        endif
      endif
c--- perform exchange
      do nu=1,4
           tmp5(nu)=p(i5,nu)
           tmp6(nu)=p(i6,nu)
           tmp7(nu)=p(i7,nu)
      enddo
      do nu=1,4
           p(5,nu)=tmp5(nu)
           p(6,nu)=tmp6(nu)
           p(7,nu)=tmp7(nu)
      enddo
      pt5=oldpt(i5)
      if (jets > 1) pt6=oldpt(i6)
      if (jets > 2) pt7=oldpt(i7)
c--- END: ordering
      endif

c--- Calculate quantities to plot
      if (jets > 0) then
        yj1=yrap(5,p)
        ptj1=pt(5,p)
      else
        yj1=99d0
        ptj1=-1d0
      endif

      if (jets >= 2) then
        yj2=yrap(6,p)
        ptj2=pt(6,p)
c        deleta=abs(yj1-yj2)
      endif

      if (jets >= 3) then
        yj3=yrap(7,p)
        ptj3=pt(7,p)
      if     (abs(yj1-yj3) >
     &      max(abs(yj1-yj2),abs(yj2-yj3))) then
c          deleta=abs(yj1-yj3)
        sumy=yj1+yj3
        yj3=yj2
      elseif (abs(yj2-yj3) >
     &      max(abs(yj1-yj2),abs(yj1-yj3))) then
c          deleta=abs(yj2-yj3)
        sumy=yj2+yj3
        yj3=yj1
      else
c        deleta=abs(yj1-yj2)
        sumy=yj1+yj2
      endif
      endif

      if (jets > 0) then
        ay345=abs(yrapthree(3,4,5,p))
      else
        ay345=99d0
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

      call bookplot(n,tag,'H rapidity',y34,wt,wt2,-5._dp,5._dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'H pt',pt34,wt,wt2,0._dp,1200._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'H abs y',abs(y34),wt,wt2,0._dp,4.8_dp,0.3_dp,'lin')
      n=n+1
      call bookplot(n,tag,'H abs eta',abs(eta34),wt,wt2,0._dp,4.8_dp,0.3_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y3',y3,wt,wt2,-5._dp,5._dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt3',pt3,wt,wt2,0._dp,200._dp,5._dp,'lin')
      n=n+1
      call bookplot(n,tag,'y4',y4,wt,wt2,-5._dp,5._dp,0.4_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt4',pt4,wt,wt2,0._dp,200._dp,5._dp,'lin')
      n=n+1

      call bookplot(n,tag,'Jet 1 pt',ptj1,wt,wt2,
     &              0._dp,1200._dp,50._dp,'log')
      n=n+1
      call bookplot(n,tag,'Jet 1 abs eta',abs(yj1),wt,wt2,
     &              0._dp,10._dp,0.25_dp,'lin')
      n=n+1

      call bookplot(n,tag,'abs y(345)',ay345,wt,wt2,
     &              0._dp,10._dp,0.5_dp,'lin')
      n=n+1

c      call bookplot(n,tag,'Jet 1 pt low',ptj1,wt,wt2,
c     &              0._dp,2000._dp,50._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'Jet 1 pt high',ptj1,wt,wt2,
c     &              0._dp,40000._dp,1000._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'Jet 1 eta',yj1,wt,wt2,
c     &              -5._dp,5._dp,0.2_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'Jet 2 pt log',ptj2,wt,wt2,
c     &              20._dp,260._dp,5._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'Jet 2 pt lin',ptj2,wt,wt2,
c     &              20._dp,260._dp,5._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'Jet 2 eta',yj2,wt,wt2,
c     &              -5._dp,5._dp,0.2_dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'Jet 3 pt log',ptj3,wt,wt2,
c     &              20._dp,170._dp,5._dp,'log')
c      n=n+1
c      call bookplot(n,tag,'Jet 3 pt lin',ptj3,wt,wt2,
c     &              20._dp,170._dp,5._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'Jet 3 eta',yj3,wt,wt2,
c     &              -5._dp,5._dp,0.2_dp,'lin')
c      n=n+1

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

