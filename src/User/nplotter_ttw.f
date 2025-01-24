!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_ttw(p,wt,wt2,switch)
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
      include 'plabel.f'
      include 'npart.f'
      real(dp):: p(mxpart,4),wt,wt2,tiny,ptthree,yrapthree,pttwo,yraptwo
      integer:: switch,n,nplotmax
      integer tag
      parameter(tiny=1.e-8_dp)
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

cc--- missing ET
c      MET=etmiss(p,vecmet)
cc--- HT(jet) - scalar sum of jet pt
c      HTjet=0._dp
c      do j=3,npart+2-switch
c        if ((plabel(j) == 'pp') .or. (plabel(j) == 'bq')
c     & .or. (plabel(j) == 'ba')) then
c          HTjet=HTjet+pt(j,p)
c        endif
c      enddo

c***********************************************************************
c                                                                      *
c     FILL HISTOGRAMS                                                  *
c                                                                      *
c***********************************************************************

c--- Call histogram routines
   99 continue

c--- "n" will count the number of histograms
      n=nextnplot

      call bookplot(n,tag,'pt top',ptthree(3,4,5,p),wt,wt2,0._dp,500._dp,25._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt antitop',ptthree(6,7,8,p),wt,wt2,0._dp,500._dp,25._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt W',pttwo(9,10,p),wt,wt2,0._dp,500._dp,25._dp,'lin')
      n=n+1

      call bookplot(n,tag,'y top',yrapthree(3,4,5,p),wt,wt2,-5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y antitop',yrapthree(6,7,8,p),wt,wt2,-5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y W',yraptwo(9,10,p),wt,wt2,-5._dp,5._dp,0.2_dp,'lin')
      n=n+1

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

c      call bookplot(n,tag,'HTjet',HTjet,wt,wt2,0._dp,800._dp,20._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'MET',MET,wt,wt2,0._dp,200._dp,10._dp,'lin')
c      n=n+1

cc--- This section tailored for comparison with CMS PAS SUS-11-020

cc--- only fill histograms if both b-quark and anti-b-quark are present
c      if ((first) .or. ((p(5,4) > tiny).and.(p(6,4) > tiny))) then

c      call bookplot(n,tag,'HTjet with 2b',HTjet,wt,wt2,
c     & 0._dp,800._dp,20._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'MET with 2b',MET,wt,wt2,0._dp,200._dp,10._dp,'lin')
c      n=n+1

c      if ((first) .or. ((MET > 30._dp) .and. (HTjet > 80._dp))) then
c      call bookplot(n,tag,'SR1',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
c      endif
c      n=n+1

c      if ((first) .or. ((MET > 120._dp) .and. (HTjet > 200._dp))) then
c      call bookplot(n,tag,'SR3',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
c      endif
c      n=n+1

c      if ((first) .or. ((MET > 50._dp) .and. (HTjet > 200._dp))) then
c      call bookplot(n,tag,'SR4',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
c      endif
c      n=n+1

c      if ((first) .or. ((MET > 50._dp) .and. (HTjet > 320._dp))) then
c      call bookplot(n,tag,'SR5',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
c      endif
c      n=n+1

c      if ((first) .or. ((MET > 120._dp) .and. (HTjet > 320._dp))) then
c      call bookplot(n,tag,'SR6',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
c      endif
c      n=n+1

cc--- double-binned histogram in MET and HT
c!      binindex=METHTdoublebin(MET,HTjet)
c!      call bookplot(n,tag,'MET-HT double bin',binindex,wt,wt2,
c!     & -0.5_dp,99.5_dp,1._dp,'lin')
c!      n=n+1

c      else

cc--- increment n by the # of histograms skipped
c!      n=n+8
c      n=n+7

c      endif

cc--- This section tailored for comparison with CMS PAS SUS-11-020

cc--- This section tailored for comparison with CMS PAS SUS-11-010

c      if ((first) .or. ((MET > 120._dp) .and. (HTjet > 400._dp))) then
c      call bookplot(n,tag,'Region 1',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
c      endif
c      n=n+1

c      if ((first) .or. ((MET > 50._dp) .and. (HTjet > 400._dp))) then
c      call bookplot(n,tag,'Region 2',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
c      endif
c      n=n+1

c      if ((first) .or. ((MET > 120._dp) .and. (HTjet > 200._dp))) then
c      call bookplot(n,tag,'Region 3',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
c      endif
c      n=n+1

c      if ((first) .or. ((MET > 100._dp) .and. (HTjet > 80._dp))) then
c      call bookplot(n,tag,'Region 4',0.5_dp,wt,wt2,0._dp,1._dp,1._dp,'lin')
c      endif
c      n=n+1


c--- This section tailored for comparison with CMS PAS SUS-11-010


cc--- single-particle plots
c      do j=3,10
c        call genplot1(p,j,tag,wt,wt2,n)
c      enddo
cc--- two-particle plots
c      call genplot2(p,3,4,tag,wt,wt2,n)
c      call genplot2(p,7,8,tag,wt,wt2,n)
c      call genplot2(p,9,10,tag,wt,wt2,n)
c      call genplot2(p,5,6,tag,wt,wt2,n)
cc--- three-particle plots
c      call genplot3(p,3,4,5,tag,wt,wt2,n)
c      call genplot3(p,6,7,8,tag,wt,wt2,n)

cc--- additional plots that may be present at NLO
c      if (abs(p(11,4)) > 1.e-8_dp) then
c        call genplot1(p,11,tag,wt,wt2,n)
c      else
c        n=n+1
c      endif

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

      function METHTdoublebin(MET,HT)
      implicit none
      include 'types.f'
      real(dp):: METHTdoublebin

c--- returns a real(dp):: number that indicates which bin
c--- the double-binned MET and HT fall into
      real(dp):: MET,HT,boundMET(11),boundHT(9)
      integer:: j,bin(2),maxbinMET,maxbinHT
      data boundMET/30._dp,40._dp,50._dp,60._dp,70._dp,80._dp,90._dp,100._dp,
     &  110._dp,120._dp,130._dp/
      data maxbinMET/11/
      data boundHT/80._dp,140._dp,200._dp,260._dp,320._dp,380._dp,440._dp,500._dp,560._dp/
      data maxbinHT/9/
      save maxbinMET,maxbinHT,boundMET,boundHT

c--- see which bin MET falls into
      bin(1)=0
      do j=1,maxbinMET-1
        if ((MET >= boundMET(j)) .and. (MET < boundMET(j+1)))
     &    bin(1)=j
      enddo
      if (MET >= boundMET(maxbinMET)) bin(1)=maxbinMET

      bin(2)=0
      do j=1,maxbinHT-1
        if ((HT >= boundHT(j)) .and. (HT < boundHT(j+1)))
     &    bin(2)=j
      enddo
      if (HT >= boundHT(maxbinHT)) bin(2)=maxbinHT

      METHTdoublebin=real((bin(1)-1)*maxbinHT+bin(2),dp)

c--- underflow
      if ((bin(1) == 0) .or. (bin(2) == 0)) then
        METHTdoublebin=0._dp
      endif

      return
      end
