!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_auto(p,wt,wt2)
          use types
          use Scalevar
      implicit none
c--- Variable passed in to this routine:

c---      p:  4-momenta of particles in the format p(i,4)
c---          with the particles numbered according to the input file
c---          and components labelled by (px,py,pz,E)

c---     wt:  weight of this event

c---    wt2:  weight^2 of this event


      include 'vegas_common.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'histo.f'
      include 'mcfmplotinfo.f'
      include 'scalevar.f'
      include 'kpart.f'
      real(dp):: p(mxpart,4),wt,wt2,tiny
      integer:: n,nplotmax,j,m,i1,i2,i3,i4,ilomomenta
      integer tag
      common/ilomomenta/ilomomenta
      common/nplotmax/nplotmax
      parameter(tiny=1.e-8_dp)
      integer, save::imaxmom
      include 'first.f'

c***********************************************************************
c                                                                      *
c     INITIAL BOOKKEEPING                                              *
c                                                                      *
c***********************************************************************

      if (first) then
c---- Initialize histograms and fix the maximum number of
c---- momenta for making single-particle plots
        tag=tagbook
        imaxmom=ilomomenta
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
      do j=3,imaxmom
        call autoplot1(p,j,tag,wt,wt2,n)
      enddo
c--- two-particle plots

        if (origkpart == kreal .or. origkpart == ktota .or. origkpart == ktodk
     &          .or. origkpart == ksnlo .or. first) then
            if ((p(imaxmom+1,4) /= 0._dp)) then
                call autoplot1(p,imaxmom+1,tag,wt,wt2,n)
            else
                n=n+2
            endif
        elseif (origkpart == knnlo .or. first) then
            if ((p(imaxmom+1,4) /= 0._dp)) then
                call autoplot1(p,imaxmom+1,tag,wt,wt2,n)
            else
                n=n+2
            endif
            if ((p(imaxmom+2,4) /= 0._dp)) then
                call autoplot1(p,imaxmom+2,tag,wt,wt2,n)
            else
                n=n+2
            endif
        endif

      j=1
      do while (mcfmplotinfo(j) > 0)
        m=mcfmplotinfo(j)
        if     (m < 10) then
c--- one-particle plots
          i1=m
          if (i1 == 0) i1=10   ! special code: 0 -> 10
          call autoplot1(p,i1,tag,wt,wt2,n)
        elseif (m < 100) then
c--- two-particle plots
          i1=m/10
          i2=mod(m,10)
          if (i2 == 0) i2=10   ! special code: 0 -> 10
          call autoplot2(p,m,i1,i2,tag,wt,wt2,n)
        elseif (m < 1000) then
c--- three-particle plots
          i1=m/100
          i2=(m-i1*100)/10
          i3=mod(m,10)
          if (i3 == 0) i3=10   ! special code: 0 -> 10
          call autoplot3(p,m,i1,i2,i3,tag,wt,wt2,n)
        elseif (m < 10000) then
c--- four-particle plots
          i1=m/1000
          i2=(m-i1*1000)/100
          i3=(m-i1*1000-i2*100)/10
          i4=mod(m,10)
          if (i4 == 0) i4=10   ! special code: 0 -> 10
          call autoplot4(p,m,i1,i2,i3,i4,tag,wt,wt2,n)
        else
          write(6,*) 'Unforeseen plot in nplotter_auto.f'
          stop
        endif
        j=j+1
      enddo

c--- We have over-counted the number of histograms by 1 at this point
      n=n-1

c--- Set the maximum number of plots, on the first call
      if (first) then
        first=.false.
        nplotmax=n
      endif

      return
      end
