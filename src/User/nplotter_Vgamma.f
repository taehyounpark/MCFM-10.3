!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_Vgamma(p,wt,wt2,switch,nd)
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
      include 'mxpart.f'
      include 'histo.f'
      include 'jetlabel.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      include 'kprocess.f'

      real(dp):: p(mxpart,4),wt,wt2
      real(dp):: etarap,pt,r
      real(dp):: y3,y4,y5,y6,pt3,pt4,pt5,pt6,re5,rea5,re6,rea6
      real(dp):: yWgam,pt34,pttwo
      real(dp):: yraptwo,m34,m345,yellgam,dot
      real(dp) :: yrapthree, y345,wttag,mTlnu
      integer:: switch,n,nplotmax,nd
      integer tag
      integer:: j,ilep,igam,inu,isojets
      include 'first.f'
      common/nplotmax/nplotmax
      logical:: dummy,passedcuts_wgamma_ew,passedcuts_wgammajet_ew
      real(dp):: ptgamm1,ptlep,mtranselgam,ygamm1,ylepgamm1,Rlepgamm1
      common/ew_observables/ptgamm1,ptlep,mtranselgam
      common/ew_observables_extra/ygamm1,ylepgamm1,Rlepgamm1
!$omp threadprivate(/ew_observables/)
!$omp threadprivate(/ew_observables_extra/)
ccccc!$omp threadprivate(first,/nplotmax/)
      logical, parameter :: CMS_SMP_20_005 = .false.


c***********************************************************************
c                                                                      *
c     INITIAL BOOKKEEPING                                              *
c                                                                      *
c***********************************************************************


      if (first) then
c--- Initialize histograms, without computing any quantities; instead
c--- set them to dummy values
        tag=tagbook
        y3=1.e3_dp
        y4=1.e3_dp
        pt3=1.e7_dp
        pt4=1.e7_dp
        pt34=1.e7_dp
c---Initialise photon
        y5=1.e3_dp
        pt5=1.e7_dp
        yWgam=1.e3_dp
        yellgam=1.e8_dp
        m34=-1._dp
        m345=-1._dp
c----Initialise jet values will not pass cuts in there is an NLO jet
        y6=1.e3_dp
        pt6=1.e7_dp
c--- If re5 is not changed by the NLO value, it will be out of
c--- the plotting range
        re5=1.e3_dp
        rea5=1.e3_dp
        re6=1.e3_dp
        rea6=1.e3_dp
        y345=1.e8_dp
        mTlnu=-1._dp
        jets=nqcdjets
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
c--- If nproc=290, plot e^+(4). If nproc=295, plot e^-(3).
c--- For nproc=300 plot e^+(4) and e^-(3)
      if((nprocbelow == 290) .or. (nprocbelow == 291) .or. (nprocbelow == 294)) then
         y4=etarap(4,p)
         pt4=pt(4,p)
      elseif((nprocbelow == 295) .or. (nprocbelow == 296) .or. (nprocbelow == 299)) then
         y3=etarap(3,p)
         pt3=pt(3,p)
      elseif(nprocbelow >= 300) then
         if(pt(3,p)>pt(4,p)) then
            y3=etarap(3,p)
            pt3=pt(3,p)
            y4=etarap(4,p)
            pt4=pt(4,p)
         else
            y3=etarap(4,p)
            pt3=pt(4,p)
            y4=etarap(3,p)
            pt4=pt(3,p)
         endif
      endif

      pt34=pttwo(3,4,p)
      m34=sqrt(max(two*dot(p,3,4),zip))
      pt5=pt(5,p)
      y5=etarap(5,p)
      re5=-1._dp
      rea5=-1._dp
      if((nprocbelow == 290) .or. (nprocbelow == 291) .or. (nprocbelow == 294)) then
         re5=R(p,4,5)
      elseif((nprocbelow == 295) .or. (nprocbelow == 296) .or. (nprocbelow == 299)) then
         re5=R(p,3,5)
      elseif((nprocbelow == 300) .or. (nprocbelow == 302)) then
         re5=R(p,3,5)
         rea5=R(p,4,5)
      endif

c---- Radiation Zero plot for all processes
      yWgam=y5-yraptwo(3,4,p)
      yellgam=1.e8_dp
      if((nprocbelow==290) .or. (nprocbelow == 291) .or. (nprocbelow == 294)) then
         yellgam=y5-y4
      elseif((nprocbelow==295) .or. (nprocbelow == 296) .or. (nprocbelow == 299)) then
         yellgam=-(y5-y3)
      endif

      mTlnu=(p(3,1)*p(4,1)+p(3,2)*p(4,2))/sqrt((p(3,1)**2+p(3,2)**2)
     &      *(p(4,1)**2+p(4,2)**2))
      mTlnu=2._dp*sqrt(p(3,1)**2+p(3,2)**2)
     &      *sqrt(p(4,1)**2+p(4,2)**2)*(1._dp-mTlnu)
      mTlnu=sqrt(max(mTlnu,zip))

      if(jets > 0) then
         pt6=pt(6,p)
         y6=etarap(6,p)
         re6=-1._dp
         rea6=-1._dp
         if((nprocbelow == 290) .or. (nprocbelow == 291) .or. (nprocbelow == 294)) then
            re6=R(p,4,6)
         elseif((nprocbelow == 295) .or. (nprocbelow == 296) .or. (nprocbelow == 299)) then
            re6=R(p,3,6)
         elseif((nprocbelow == 300) .or. (nprocbelow == 302)) then
            re6=R(p,3,6)
            rea6=R(p,4,6)
         endif
      else ! put out of range of plotting
         pt6=1.e7_dp
         y6=1.e7_dp
         re6=1.e7_dp
         rea6=1.e7_dp
      endif


      if(nprocbelow<300) then
c---  transverse mass of (e-gam,nu) system for Wgamma
        if ((nprocbelow == 290) .or. (nprocbelow == 291)
     & .or. (nprocbelow == 292) .or. (nprocbelow == 294)) then
          inu=3
          ilep=4
        else
          inu=4
          ilep=3
        endif
        igam=5
        if ((nprocbelow == 291) .or. (nprocbelow == 296)) then
          if (p(6,4) > 1.e-8_dp) then
            if (pt(6,p) > pt(5,p)) then
              pt5=pt(6,p)
              igam=6
            endif
          endif
        endif
        m345=(p(ilep,4)+p(igam,4))**2-(p(ilep,1)+p(igam,1))**2
     &      -(p(ilep,2)+p(igam,2))**2-(p(ilep,3)+p(igam,3))**2
        m345=m345+(p(ilep,1)+p(igam,1))**2
     &           +(p(ilep,2)+p(igam,2))**2
        m345=sqrt(max(m345,zip))+sqrt(p(inu,1)**2+p(inu,2)**2)
        m345=m345**2
        do j=1,2
           m345=m345-(p(3,j)+p(4,j)+p(5,j))**2
        enddo
        m345=sqrt(max(m345,zip))

c---  invariant mass of (Z,gam) system for Zgamma
      else
         m345=(p(3,4)+p(4,4)+p(5,4))**2
         do j=1,3
           m345=m345-(p(3,j)+p(4,j)+p(5,j))**2
         enddo
         m345=sqrt(max(m345,zip))
      endif

      y345 = yrapthree(3,4,5,p)

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

      if (CMS_SMP_20_005) then

        if (kcase == kWga_ew) then
          pt5=ptgamm1
          y5=ygamm1
          Re5=Rlepgamm1
          yellgam=ylepgamm1
          m345=mtranselgam
        endif

c Table 2
        call bookplot(n,tag,'pt5a',pt5,wt,wt2,30._dp,50._dp,20._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt5b',pt5,wt,wt2,50._dp,70._dp,20._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt5c',pt5,wt,wt2,70._dp,100._dp,30._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt5d',pt5,wt,wt2,100._dp,150._dp,50._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt5e',pt5,wt,wt2,150._dp,200._dp,50._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt5f',pt5,wt,wt2,200._dp,300._dp,100._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt5g',pt5,wt,wt2,300._dp,500._dp,200._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt5h',pt5,wt,wt2,500._dp,800._dp,300._dp,'lin')
        n=n+1
        call bookplot(n,tag,'pt5i',pt5,wt,wt2,800._dp,1500._dp,700._dp,'lin')
        n=n+1
c Table 3
        call bookplot(n,tag,'y5a',y5,wt,wt2,-2.5_dp,-1.9_dp,0.6_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5b',y5,wt,wt2,-1.9_dp,-1.5_dp,0.4_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5c',y5,wt,wt2,-1.5_dp,-1.1_dp,0.4_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5d',y5,wt,wt2,-1.1_dp,-0.75_dp,0.35_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5e',y5,wt,wt2,-0.75_dp,-0.45_dp,0.3_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5f',y5,wt,wt2,-0.45_dp,-0.15_dp,0.3_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5g',y5,wt,wt2,-0.15_dp,0.15_dp,0.3_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5h',y5,wt,wt2,0.15_dp,0.45_dp,0.3_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5i',y5,wt,wt2,0.45_dp,0.75_dp,0.3_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5j',y5,wt,wt2,0.75_dp,1.1_dp,0.35_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5k',y5,wt,wt2,1.1_dp,1.5_dp,0.4_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5l',y5,wt,wt2,1.5_dp,1.9_dp,0.4_dp,'lin')
        n=n+1
        call bookplot(n,tag,'y5m',y5,wt,wt2,1.9_dp,2.5_dp,0.6_dp,'lin')
        n=n+1
c Table 4
        call bookplot(n,tag,'Re5a',Re5,wt,wt2,0.7_dp,4._dp,0.3_dp,'lin')
        n=n+1
        call bookplot(n,tag,'Re5b',Re5,wt,wt2,4._dp,5._dp,0.5_dp,'lin')
        n=n+1
c Table 5
        call bookplot(n,tag,'yellgama',yellgam,wt,wt2,-5._dp,-3.4_dp,1.6_dp,'lin')
        n=n+1
        call bookplot(n,tag,'yellgamb',yellgam,wt,wt2,-3.4_dp,-1.8_dp,0.8_dp,'lin')
        n=n+1
        call bookplot(n,tag,'yellgamc',yellgam,wt,wt2,-1.8_dp,1.8_dp,0.4_dp,'lin')
        n=n+1
        call bookplot(n,tag,'yellgamd',yellgam,wt,wt2,1.8_dp,3.4_dp,0.8_dp,'lin')
        n=n+1
        call bookplot(n,tag,'yellgame',yellgam,wt,wt2,3.4_dp,5._dp,1.6_dp,'lin')
        n=n+1
c Table 6
        call bookplot(n,tag,'mTclustera',m345,wt,wt2,0._dp,100._dp,100._dp,'lin')
        n=n+1
        call bookplot(n,tag,'mTclusterb',m345,wt,wt2,100._dp,300._dp,50._dp,'lin')
        n=n+1
        call bookplot(n,tag,'mTclusterc',m345,wt,wt2,300._dp,500._dp,100._dp,'lin')
        n=n+1
        call bookplot(n,tag,'mTclusterd',m345,wt,wt2,500._dp,800._dp,150._dp,'lin')
        n=n+1
        call bookplot(n,tag,'mTclustere',m345,wt,wt2,800._dp,1200._dp,200._dp,'lin')
        n=n+1
        call bookplot(n,tag,'mTclusterf',m345,wt,wt2,1200._dp,1700._dp,500._dp,'lin')
        n=n+1
        call bookplot(n,tag,'mTclusterg',m345,wt,wt2,1700._dp,2500._dp,800._dp,'lin')
        n=n+1
c Table 7
        call bookplot(n,tag,'jts',real(jets,dp),wt,wt2,-0.5_dp,2.5_dp,1._dp,'lin')
        n=n+1
        isojets=jets
        if (jets == 2) then
          if ((R(p,ilep,7) < 0.4_dp) .or. (R(p,5,7) < 0.4_dp)) isojets=isojets-1
        endif
        if (jets >= 1) then
          if ((R(p,ilep,6) < 0.4_dp) .or. (R(p,5,6) < 0.4_dp)) isojets=isojets-1
        endif
        call bookplot(n,tag,'isojets',real(isojets,dp),wt,wt2,-0.5_dp,2.5_dp,1._dp,'lin')
        n=n+1
c Table 8
        if ((isojets == 0) .and. (m345 > 150._dp)) then
          wttag=1._dp
        else
          wttag=0._dp
        endif
        call bookplot(n,tag,'yellgamvetoa',yellgam,wttag*wt,wttag*wt2,-5._dp,-3.4_dp,1.6_dp,'lin')
        n=n+1
        call bookplot(n,tag,'yellgamvetob',yellgam,wttag*wt,wttag*wt2,-3.4_dp,-1.8_dp,0.8_dp,'lin')
        n=n+1
        call bookplot(n,tag,'yellgamvetoc',yellgam,wttag*wt,wttag*wt2,-1.8_dp,1.8_dp,0.4_dp,'lin')
        n=n+1
        call bookplot(n,tag,'yellgamvetod',yellgam,wttag*wt,wttag*wt2,1.8_dp,3.4_dp,0.8_dp,'lin')
        n=n+1
        call bookplot(n,tag,'yellgamvetoe',yellgam,wttag*wt,wttag*wt2,3.4_dp,5._dp,1.6_dp,'lin')
        n=n+1

      else

c--- cross-section with pt(photon)>60
      if (pt5 > 60._dp) then
        wttag=wt
      else
        wttag=0._dp
      endif
      call bookplot(n,tag,'xsecpt60',0.5_dp,wttag,wttag**2,zip,1._dp,1._dp,'lin')
      n=n+1

c--- cross-section with pt(photon)>90
      if (pt5 > 90._dp) then
        wttag=wt
      else
        wttag=0._dp
      endif
      call bookplot(n,tag,'xsecpt90',0.5_dp,wttag,wttag**2,zip,1._dp,1._dp,'lin')
      n=n+1

      if((nprocbelow == 290).or.(nprocbelow == 291).or.(nprocbelow==294).or.(nprocbelow>=300)) then
         call bookplot(n,tag,'y4',y4,wt,wt2,-5._dp,5._dp,0.2_dp,'lin')
         n=n+1
         call bookplot(n,tag,'pt4',pt4,wt,wt2,zip,100._dp,two,'lin')
         n=n+1
      endif
      if((nprocbelow==295).or.(nprocbelow == 296).or.(nprocbelow==299).or.(nprocbelow>=300)) then
         call bookplot(n,tag,'y3',y3,wt,wt2,-5._dp,5._dp,0.2_dp,'lin')
         n=n+1
         call bookplot(n,tag,'pt3',pt3,wt,wt2,zip,100._dp,two,'lin')
         n=n+1
      endif
      call bookplot(n,tag,'boson invariant mass',
     & m34,wt,wt2,zip,200._dp,5._dp,'lin')
      n=n+1

      call bookplot(n,tag,'mT(l,nu)',mTlnu,wt,wt2,zip,150._dp,5._dp,'lin')
      n=n+1

      call bookplot(n,tag,'DeltaRe5',re5,wt,wt2,zip,5._dp,0.1_dp,'lin')
      n=n+1
      if(nprocbelow==300) then
      call bookplot(n,tag,'DeltaRea5',rea5,wt,wt2,
     & zip,5._dp,0.1_dp,'lin')
      n=n+1
      endif
      call bookplot(n,tag,'pt34',pt34,wt,wt2,zip,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'y5',y5,wt,wt2,-4.5_dp,4.5_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt5',pt5,wt,wt2,zip,1500._dp,30._dp,'log')
      n=n+1
      if (jets == 0) then
        wttag=wt
      else
        wttag=0._dp
      endif
      call bookplot(n,tag,'vetopt5',pt5,wttag,wttag**2,zip,1500._dp,30._dp,'log')
      n=n+1
      if(nprocbelow<300) then
      call bookplot(n,tag,'transverse cluster mass, m(e-gam,nu)',
     & m345,wt,wt2,zip,2500._dp,100._dp,'lin')
      n=n+1
      if (jets == 0) then
        wttag=wt
      else
        wttag=0._dp
      endif
      call bookplot(n,tag,'vetotransverse cluster mass, m(e-gam,nu)',
     & m345,wttag,wttag**2,zip,2500._dp,100._dp,'lin')
      n=n+1
      call bookplot(n,tag,'ydiff(ellgam)',yellgam,wt,wt2,-5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      if ((jets == 0) .and. (m345 > 150._dp)) then
        wttag=wt
      else
        wttag=0._dp
      endif
      call bookplot(n,tag,'vetoydiff(ellgam)',yellgam,wttag,wttag**2,-5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      if ((jets == 0) .and. (m345 > 110._dp)) then
        wttag=wt
      else
        wttag=0._dp
      endif
      call bookplot(n,tag,'veto7ydiff(ellgam)',yellgam,wttag,wttag**2,-5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      else
      call bookplot(n,tag,'(Z,gam) invariant mass',
     & m345,wt,wt2,zip,200._dp,5._dp,'lin')
      n=n+1
      call bookplot(n,tag,'(Z,gam) invariant mass',
     & m345,wt,wt2,zip,4000._dp,50._dp,'log')
      n=n+1
      endif
      call bookplot(n,tag,'ydiff(Vgam)',yWgam,wt,wt2,-5._dp,5._dp,
     & 0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt6',pt6,wt,wt2,zip,200._dp,two,'lin')
      n=n+1
      call bookplot(n,tag,'y6',y6,wt,wt2,-5._dp,5._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DeltaRe6',re6,wt,wt2,zip,5._dp,0.1_dp,'lin')
      n=n+1
c      if(nprocbelow < 305) then
      call bookplot(n,tag,'DeltaRea6',rea6,wt,wt2,
     & zip,5._dp,0.1_dp,'lin')
      n=n+1

      call bookplot(n,tag,'y345',y345,wt,wt2,-5._dp,5._dp,0.2_dp,'lin')
      n=n+1

c get photon pt that passes more generic EW  cuts (fills common block)
      if (kcase == kWgajew) then
        dummy=passedcuts_wgammajet_ew(switch,p)
      else
        dummy=passedcuts_wgamma_ew(switch,p)
      endif
      call bookplot(n,tag,'ptgamm1',ptgamm1,wt,wt2,zip,1500._dp,30._dp,'log')
      n=n+1
      call bookplot(n,tag,'ptlep',ptlep,wt,wt2,zip,1500._dp,30._dp,'log')
      n=n+1
      call bookplot(n,tag,'mtranselgam',mtranselgam,wt,wt2,zip,2500._dp,100._dp,'lin')
      n=n+1

c--- cross-section with pt(photon)>60
      if (ptgamm1> 60._dp) then
        wttag=wt
      else
        wttag=0._dp
      endif
      call bookplot(n,tag,'xsecptgamm160',0.5_dp,wttag,wttag**2,zip,1._dp,1._dp,'lin')
      n=n+1

c--- cross-section with pt(photon)>90
      if (ptgamm1 > 90._dp) then
        wttag=wt
      else
        wttag=0._dp
      endif
      call bookplot(n,tag,'xsecptgamm190',0.5_dp,wttag,wttag**2,zip,1._dp,1._dp,'lin')
      n=n+1

      endif

      call bookplot(n,tag,'total cross', 0.5d0,wt,wt2,0d0,100d0,1d0,'lin')
      n=n+1
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
