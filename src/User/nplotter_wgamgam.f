!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_wgamgam(p,wt,wt2,switch,nd)
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
      real(dp):: p(mxpart,4),wt,wt2,etmiss
      real(dp):: s34,m34
      real(dp):: ptgam5,ptgam6,ptgam
      real(dp):: m3456,s56,m56
      real(dp):: ptmiss,etvec(4)
      integer:: switch,n,nplotmax
      integer tag
      integer:: nd
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
        ptgam=0._dp
        ptmiss=0._dp
        m34=0._dp
        m3456=0._dp
        m56=0._dp
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
c-----m(l,l)
      s34=2._dp*(p(4,4)*p(3,4)-p(4,1)*p(3,1)-p(4,2)*p(3,2)
     &        -p(4,3)*p(3,3))
      m34=sqrt(s34)
c-----m3456
      m3456=sqrt((p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     &           -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &           -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     &           -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2)
c-----m(gam,gam)
      s56=2._dp*(p(5,4)*p(6,4)-p(5,1)*p(6,1)-p(5,2)*p(6,2)
     &        -p(5,3)*p(6,3))
      m56=sqrt(s56)
c-----pT(photon)-hardest
      ptgam5 = 0._dp
      ptgam5 = 0._dp
      ptgam5 = sqrt(p(5,1)**2+p(5,2)**2)
      ptgam6 = sqrt(p(6,1)**2+p(6,2)**2)
      ptgam  = max(ptgam5,ptgam6)
c-----missing ET
      ptmiss=etmiss(p,etvec)

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

c-----m34
      call bookplot(n,tag,'m(l,l)',m34,wt,wt2,0._dp,200._dp,2._dp,'lin')
      n=n+1
c-----m34
      call bookplot(n,tag,'m(l,l)',m34,wt,wt2,0._dp,500._dp,5._dp,'lin')
      n=n+1
c-----m3456
      call bookplot(n,tag,'m(l,l,gam,gam)',
     & m3456,wt,wt2,0._dp,500._dp,5._dp,'lin')
      n=n+1
c-----m56
      call bookplot(n,tag,'m(gam,gam)',m56,wt,wt2,0._dp,500._dp,5._dp,'lin')
      n=n+1
c-----hardest photon pT
      call bookplot(n,tag,'pT(gam)',ptgam,wt,wt2,0._dp,500._dp,10._dp,'lin')
      n=n+1
c-----missing transverse momentum
      call bookplot(n,tag,'pT(miss)',ptmiss,wt,wt2,0._dp,500._dp,5._dp,'lin')
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





c-----m(l,l,gamma)
c-----m345
c     m345=sqrt((p(3,4)+p(4,4)+p(5,4))**2-(p(3,1)+p(4,1)+p(5,1))**2
c    &          -(p(3,2)+p(4,2)+p(5,2))**2-(p(3,3)+p(4,3)+p(5,3))**2)
c-----m346
c     m346=sqrt((p(3,4)+p(4,4)+p(6,4))**2-(p(3,1)+p(4,1)+p(6,1))**2
c    &          -(p(3,2)+p(4,2)+p(6,2))**2-(p(3,3)+p(4,3)+p(6,3))**2)
c------photon-photon separation
c       Rgamgam=R(p,5,6)
c-----m345
c      call bookplot(n,tag,'m(l,l,gam1)',m345,wt,wt2,0._dp,200._dp,5._dp,'lin')
c      n=n+1
c-----m346
c      call bookplot(n,tag,'m(l,l,gam2)',m346,wt,wt2,0._dp,200._dp,5._dp,'lin')
c      n=n+1
c-----R(gamma,gamma)
c      call bookplot(n,tag,'R(gam,gam)',Rgamgam,wt,wt2,
c     &0_dp,5._dp,0.1_dp,'lin')
c      n=n+1

