!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine nplotter_Ztjdk(p,wt,wt2,switch)
      implicit none
      include 'types.f'
c -- Taylored for process pp -> Z(mu-(p3)+mu+(p4))+t(-->nu(p5)+e^+(p6)+b(p7))+q(p8) + g(p9)
c -- Because the momentum assignments for charged leptons and neutrinos are different, this routine
c -- should be modified when considering tbar-Z production.
c -- R. Rontsch 2013-03-01

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
      include 'jetlabel.f'
      real(dp):: p(mxpart,4),wt,wt2,pt,
     & pttwo,yraptwo,etaraptwo,etarapthree,yrapthree,
     & tiny,etarap,r
      real(dp):: wtjet(6),wtjet2(6)
      real(dp):: ptl,pta,ptep,ptmiss,ptj1,ptj2,ptb,ptt,
     & yl,ya,yla,ylaj,yjep,ybep,etab,etala,
     & etalaj,etalab,etaep,etaW,etaWb,etaWj,
     & mla,mlaj,mlab,mWb,mWj,mZbep,mZjep,etaj1,etaj2,
     & Dylaj,Detalaj,ptla,Detajb,DetaZt,DetaZj,DetaZb,Detatj,DetaZep,
     & DetaWb,DetaWj,DRbl,DRba,DRbep,DRjl,DRja,DRjep,HTTOT,dphibep,
     & dphijep,dphila,dphimissep,mtransW,fphi
      real(dp):: twomass,threemass,fourmass
      integer:: switch,n,nplotmax,j,iorderl(mxpart),
     & jetindex(mxpart),ijet
      integer:: countljet,countbjet
      real(dp):: pljet(mxpart,4),pbjet(mxpart,4),pbl(mxpart,4),
     &     ptljet(mxpart)
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
c---  Add event in histograms
        tag=tagplot
      endif

c***********************************************************************
c                                                                      *
c     DEFINITIONS OF QUANTITIES TO PLOT                                *
c                                                                      *
c***********************************************************************

c --  jet-binned weights
       wtjet(:)=0._dp
       wtjet2(:)=0._dp
      if (jets == 0) then
         wtjet(1)=wt
         wtjet2(1)=wt2
      elseif (jets == 1) then
         wtjet(2)=wt
         wtjet2(2)=wt2
      elseif (jets == 2) then
         wtjet(3)=wt
         wtjet2(3)=wt2
      elseif (jets == 3) then
         wtjet(4)=wt
         wtjet2(4)=wt2
      elseif (jets == 4) then
         wtjet(5)=wt
         wtjet2(5)=wt2
      elseif (jets == 5) then
         wtjet(6)=wt
         wtjet2(6)=wt2
      endif

c -- start with jet definitions
      ijet=0
       do j=3,9
          if ((plabel(j) == 'pp') .or. (plabel(j) == 'bq')
     &         .or.(plabel(j) == 'ba')) then
             if (p(j,4) > tiny) then
                ijet=ijet+1
                jetindex(ijet)=j
             endif
          endif
       enddo

c -- identify jets as either light or b-tagged. Their are countljet light jets and countbjet b-tagged jets.
       call idjet(p,jetindex,countljet,countbjet,
     &     pljet,pbjet)

c -- pt-order the light jets
       ptljet=0._dp
       do j=1,countljet
          ptljet(j)=pt(j,pljet)
       enddo
       call arraysort(countljet,ptljet,iorderl)

c -- now I construct a set of momentum where 7 is defined as the bottom jet,
c     8,9,..., are defined as the light jets, and 3,4 are the leptons which reproduce the Z-mass
c -- this will be used to define observables

c -- first, all the leptonic momenta
       pbl =0._dp
       pbl(1,:)=p(1,:)
       pbl(2,:)=p(2,:)
       pbl(3,:)=p(3,:)
       pbl(4,:)=p(4,:)
       pbl(5,:)=p(5,:)
       pbl(6,:)=p(6,:)

c -- now jet momenta
       pbl(7,:)=pbjet(1,:)
       do j=1,countljet
          pbl(7+j,:)=pljet(iorderl(j),:)
       enddo

c -- define observables
       ptl=pt(3,pbl)
       pta=pt(4,pbl)
       ptla=pttwo(3,4,pbl)
       ptep=pt(6,pbl)
       ptmiss=pt(5,pbl)
       if (countljet > 0) ptj1=pt(8,pbl)
       if (countbjet > 0) ptb=pt(7,pbl)

       yl=etarap(3,pbl)
       ya=etarap(4,pbl)
       yla=yraptwo(3,4,pbl)
       if (countljet > 0) ylaj=yrapthree(3,4,8,pbl)
       if (countljet > 0) etalaj=etarapthree(3,4,8,pbl)
       if (countbjet > 0) etalab=etarapthree(3,4,7,pbl)
       etaep=etarap(6,pbl)
       if (countljet > 0) yjep=yraptwo(6,8,pbl)
       if (countbjet > 0) ybep=yraptwo(6,7,pbl)
       if (countljet > 0) etaj1=etarap(8,pbl)
       if (countbjet > 0) etab=etarap(7,pbl)
       etala=etaraptwo(3,4,pbl)
       etaW=etaraptwo(5,6,pbl)
       if (countbjet > 0) etaWb=etarapthree(5,6,7,pbl)
       if (countljet > 0) etaWj=etarapthree(5,6,8,pbl)

c --  for NLO, second light jet
       if (countljet >= 2) then
         ptj2=pt(9,pbl)
         etaj2=etarap(9,pbl)
       else
         ptj2=-1._dp
         etaj2=1d5
       endif

       mla=twomass(3,4,pbl)
       if (countljet > 0) mlaj=threemass(3,4,8,pbl)
       if (countbjet > 0) mlab=threemass(3,4,7,pbl)
       if (countbjet > 0) mWb=threemass(5,6,7,pbl)
       if (countljet > 0) mWj=threemass(5,6,8,pbl)
       if (countbjet > 0) mZbep=fourmass(3,4,6,7,pbl)
       if (countljet > 0) mZjep=fourmass(3,4,6,8,pbl)
       if ((countbjet > 0) .and. (countljet > 0)) Detajb=etaj1-etab
       if (countbjet > 0) DetaZt=etala-etaWb
       if (countljet > 0) DetaZj=etala-etaj1
       if (countbjet > 0) DetaZb=etala-etab
       DetaZep=etala-etaep
       if ((countbjet > 0) .and. (countljet > 0)) Detatj=etaj1-etaWb
       if (countbjet > 0) DetaWb=etaW-etab
       if (countljet > 0) DetaWj=etaW-etaj1
       ptt=ptmiss+ptb+ptep

       if (countljet > 0) Dylaj=yla-etaj1
       if (countljet > 0) Detalaj=etala-etaj1

       if (countbjet > 0) DRbl=R(pbl,3,7)
       if (countbjet > 0) DRba=R(pbl,4,7)
       if (countbjet > 0) DRbep=R(pbl,6,7)
       if (countljet > 0) DRjl=R(pbl,3,8)
       if (countljet > 0) DRja=R(pbl,4,8)
       if (countljet > 0) DRjep=R(pbl,6,8)

       HTTOT=ptl+pta+ptep+ptb+ptj1+ptj2+ptmiss

       if (countbjet > 0) dphibep=fphi(6,7,pbl)
       if (countljet > 0) dphijep=fphi(6,8,pbl)
       dphila=fphi(3,4,pbl)
       dphimissep=fphi(5,6,pbl)
       mtransW=sqrt(2._dp*ptmiss*ptep*(1._dp-cos(dphimissep)))
       
       if (countbjet == 0) then
         ptb=0._dp
         etalab=99._dp
         ybep=99._dp
         etab=99._dp
         etaWb=99._dp
         mlab=0._dp
         mWb=0._dp
         mZbep=0._dp
         DetaZt=99._dp
         DetaZb=99._dp
         DetaWb=99._dp
         DRbl=99._dp
         DRba=99._dp
         DRbep=99._dp
         dphibep=99._dp
       endif
       if (countljet == 0) then
         ptj1=0._dp
         ylaj=99._dp
         etalaj=99._dp
         yjep=99._dp
         etaj1=99._dp
         etaWj=99._dp
         mlaj=0._dp
         mWj=0._dp
         mZjep=0._dp
         DetaZj=99._dp
         DetaWj=99._dp
         Dylaj=99._dp
         Detalaj=99._dp
         DRjl=99._dp
         DRja=99._dp
         DRjep=99._dp
         dphijep=99._dp
       endif
       if ((countbjet == 0) .and. (countljet == 0)) then
         Detajb=99._dp
         Detatj=99._dp
       endif

c**************** END DEFINITIONS **************




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

c -- xsecs with different numbers of jets
c      call bookplot(n,tag,'0j',1._dp,nojetwt,nojetwt2,0._dp,2._dp,1._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'1j',1._dp,onejetwt,onejetwt2,0._dp,2._dp,1._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'2j',1._dp,twojetwt,twojetwt2,0._dp,2._dp,1._dp,'lin')
c      n=n+1
c      call bookplot(n,tag,'3j',1._dp,thrjetwt,thrjetwt2,0._dp,2._dp,1._dp,'lin')
c      n=n+1

c -- One-particle plots
c -- lepton plots
c -- plot jet-binned cross-sections
      call bookplot(n,tag,'0jet',0.5_dp,wtjet(1),wtjet2(1),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'1jet',0.5_dp,wtjet(2),wtjet2(2),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'2jet',0.5_dp,wtjet(3),wtjet2(3),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'3jet',0.5_dp,wtjet(4),wtjet2(4),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'4jet',0.5_dp,wtjet(5),wtjet2(5),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'5jet',0.5_dp,wtjet(6),wtjet2(6),
     &     0._dp,1._dp,1._dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt(mu-)',ptl,wt,wt2,
     & 0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(mu-)',yl,wt,wt2,
     & -2.2_dp,2.2_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt(e+)',pta,wt,wt2,
     & 0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(mu+)',ya,wt,wt2,
     & -2.2_dp,2.2_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt(e+)',ptep,wt,wt2,
     & 0._dp,200._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(e+)',etaep,wt,wt2,
     & -2.2_dp,2.2_dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt,miss',ptmiss,wt,wt2,
     & 0._dp,200._dp,10._dp,'log')
      n=n+1
c -- jet plots
      call bookplot(n,tag,'pt(j1)',ptj1,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(j1)',etaj1,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt(j2)',ptj2,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(j2)',etaj2,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1
      call bookplot(n,tag,'pt(b)',ptb,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'eta(b)',etab,wt,wt2,-6._dp,6._dp,0.2_dp,'lin')
      n=n+1

c -- Two particle plots
      call bookplot(n,tag,'pt(mu-mu+)',ptla,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'m(mu-mu+)',mla,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'m(W,trans)',mtransW,wt,wt2,0._dp,500._dp,
     &     10._dp,'log')
      n=n+1
      call bookplot(n,tag,'y(mu-mu+)',yla,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y(b,e+)',ybep,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'y(j,e+)',yjep,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(mu-mu+)',etala,wt,wt2,-5._dp,5._dp,
     &     0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(W)',etaW,wt,wt2,-5._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Deta(j,b)',Detajb,wt,wt2,
     &     -5._dp,5._dp,0.1_dp,'lin')
      n=n+1

      call bookplot(n,tag,'DR(b,l)',DRbl,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DR(b,a)',DRba,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DR(b,e+)',DRbep,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DR(j,l)',DRjl,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DR(j,mu+)',DRja,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'DR(j,e+)',DRjep,wt,wt2,0._dp,5._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dphi(mu-mu+)',dphila,wt,wt2,0._dp,3.14_dp,
     &     0.157_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dphi(b,e+)',dphibep,wt,wt2,0._dp,3.14_dp,
     &     0.157_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dphi(j,e+)',dphijep,wt,wt2,0._dp,3.14_dp,
     &     0.157_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Dphi(mu-mu+)',dphila,wt,wt2,0._dp,3.14_dp,
     &     0.157_dp,'lin')
      n=n+1

c -- Three particle plots
      call bookplot(n,tag,'m(mu-mu+b)',mlab,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'m(Wb)',mWb,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'m(Wj)',mWj,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'m(mu-mu+j)',mlaj,wt,wt2,
     &     100._dp,300._dp,5._dp,'log')
      n=n+1
      call bookplot(n,tag,'y(mu-mu+j)',ylaj,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(mu-mu+j)',etalaj,
     &     wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(mu-mu+b)',etalab,
     &     wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(Wb)',etaWb,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'eta(Wj)',etaWj,wt,wt2,-4._dp,4._dp,0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Delta y(mu-mu+,j)',Dylaj,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Delta eta(mu-mu+,j)',Detalaj,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Delta eta(Z,j)',DetaZj,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Delta eta(Z,b)',DetaZb,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Delta eta(W,j)',DetaWj,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Delta eta(W,b)',DetaWb,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1
      call bookplot(n,tag,'Delta eta(Z,e+)',DetaZep,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
c -- four particle plots
      n=n+1
      call bookplot(n,tag,'m(Zbe+)',mZbep,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'m(Zje+)',mZjep,wt,wt2,0._dp,500._dp,10._dp,'log')
c -- involving t
       n=n+1
      call bookplot(n,tag,'Delta eta(t,j)',Detatj,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      call bookplot(n,tag,'Delta eta(t,Z)',DetaZt,wt,wt2,-5._dp,5._dp,
     & 0.1_dp,'lin')
      n=n+1
c--total momentum plots
      call bookplot(n,tag,'pt(t)',ptt,wt,wt2,0._dp,500._dp,10._dp,'log')
      n=n+1
      call bookplot(n,tag,'HTTOT',HTTOT,wt,wt2,0._dp,500._dp,10._dp,'log')
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
