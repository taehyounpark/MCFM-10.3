!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine gen_stop(r,njets,p,wt,*)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'mxdim.f'
      include 'limits.f'
      include 'zerowidth.f'
      include 'kprocess.f'
      include 'reset.f'
      include 'kpart.f'
      include 'x1x2.f'
      include 'energy.f'
      include 'notag.f'
      include 'first.f'
      include 'ipsgen.f'
      include 'jetcuts.f'
c---- Generate phase space for 2-->2+n process
c---- with (345) being a top and 6,..,5+n the jets
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign
c---- from physical values
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(3n-4), where n is the number
c----  of final state particles)
c---- This routine has a minimum of 4 final state particles, hence
c---- the twopi**2 correction factor is given by the ratio of
c---- (1/twopi)**(3n-4) present in the phase space and the factor
c---- of [(1/twopi)**2]**(n-1) from the number of branchings
c---- For the specific case 'ttdkay' where one of the jets is
c---- associated with the top quark decay, we must add an extra
c---- factor of (1/twopi) since the number of jets generated is
c---- larger than the value of 'njets' passed
      real(dp):: r(mxdim)
      real(dp):: p(mxpart,4),psumjet(4),pcm(4),Q(4)
      real(dp):: wt,wt0,wtbg
      real(dp):: pt,etamax,etamin
      real(dp):: y,sinhy,coshy,phi,mv2,wtbw,mjets
      real(dp):: ybar,ptsumjet2,ycm,sumpst,q0st,rshat,dely,xjac
      real(dp), save :: ptjetmin_use,etajetmax_use,ptbreak
!$omp threadprivate(ptjetmin_use,etajetmax_use,ptbreak)
      real(dp):: plstar,estar,plstarsq,y5starmax,y5starmin,mtrans
      real(dp):: bm(4),wp(4),nn(4),ep(4),pbg(4),g(4),wtwp,wtepnn
      integer:: j,nu,njets,ijet,in
      logical:: oldzerowidth
      parameter(wt0=1._dp/twopi**2)

      if (first .or. reset) then
        first=.false.
        reset=.false.
        ptjetmin_use = ptjetmin
        etajetmax_use = etajetmax
        if ((kpart==kreal) .or. (notag > 0)) then
          etajetmax_use=100._dp
        endif
      endif

      ptbreak=ptjetmin_use

      do nu=1,4
        do j=1,5+njets
          p(j,nu)=0._dp
        enddo
        psumjet(nu)=0._dp
        pcm(nu)=0._dp
      enddo

      wt=2._dp*pi

      do ijet=1,njets
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
        wt=wt/16._dp/pi**3

c       ! some other ways to generate in pT
c       ! directly in pT
c       pt = sqrts*r(ijet)/2._dp
c       xjac = r(ijet)*sqrts**2/4._dp
c       ! 1/pt
c       pt = sqrts*r(ijet)/(sqrts+2*r(ijet)-sqrts*r(ijet))
c       xjac = r(ijet)*sqrts**3 / (2*r(ijet) - sqrts*(r(ijet)-1))**3
c       ! -log(pT)
c       pt = exp(log(sqrts/2._dp) - (1-r(ijet))/r(ijet))
c       xjac = exp(2._dp*log(sqrts/2._dp) - 2._dp*(1-r(ijet))/r(ijet))/r(ijet)**2

        ! for nnlo double real ipsgen=1
        !if ((kcase == kbq_tpq_jet) .and. (kpart == kreal) .and. usescet .and. (ipsgen == 1)) then
c           if (ijet == 1) then
                ! in 1/sqrt(pT)
                pt = r(ijet)**2*sqrts/2._dp
                xjac = pt*sqrts/2._dp*2._dp*r(ijet)
c          else
c              ! 1/pt^2
c              pt = 1._dp / sqrt(4._dp/sqrts**2 + (1-r(ijet))/r(ijet))
c              xjac = sqrts**4 / 2._dp / (sqrts**2*(r(ijet)-1._dp) - 4*r(ijet))**2
c           endif
        ! for nlo real, ipsgen == 2
        !elseif ((kcase == kbq_tpq) .and. (kpart == kreal) .and. (.not. usescet) .and. (ipsgen == 2)) then
            !error stop "todo"
        !else
            !error stop "please explicitly setup gen_stop for this case"
            !call genpt(r(ijet),ptbreak,.true.,pt,xjac)
        !endif
        wt=wt*xjac

        etamax=sqrts/2._dp/pt
        if (etamax**2 <= 1._dp) then
c            write(6,*) 'etamax**2 <= 1._dp in gen_stop.f',etamax**2
            wt=0._dp
            return 1
        endif

        etamax=log(etamax+sqrt(etamax**2-1._dp))

        etamax=min(etamax,etajetmax_use)
        y=etamax*(2._dp*r(njets+ijet)-1._dp)
        wt=wt*2._dp*etamax

        sinhy=sinh(y)
        coshy=sqrt(1._dp+sinhy**2)

        p(5+ijet,4)=pt*coshy

        phi=2._dp*pi*r(2*njets+ijet)
        wt=wt*2._dp*pi

        p(5+ijet,1)=pt*cos(phi)
        p(5+ijet,2)=pt*sin(phi)
        p(5+ijet,3)=pt*sinhy

c--- generate a b-quark as particle 6 for s-channel processes
        if (((kcase==kt_bbar) .or. (kcase==ktdecay))
     &   .and. (ijet == 1)) then
          mtrans=sqrt(pt**2+mb**2)
          p(5+ijet,4)=mtrans*coshy
          p(5+ijet,3)=mtrans*sinhy
      endif

        do nu=1,4
          psumjet(nu)=psumjet(nu)+p(5+ijet,nu)
        enddo
      enddo

      if (kcase == ktopanom) then
          call breitw(r(3*njets+7),wsqmin,wsqmax,mt,twidth,mv2,wtbw)
      else
c--- now generate Breit-Wigner, but always with zero width
          oldzerowidth=zerowidth
          zerowidth=.true.
          call breitw(one,wsqmin,wsqmax,mt,twidth,mv2,wtbw)
          zerowidth=oldzerowidth
      endif
      wt=wt*wtbw
c--- for one jet, mjets must be exactly zero
      if (njets == 1) then
        mjets=0._dp
      else
c--- invariant mass of jets
        mjets=psumjet(4)**2-psumjet(1)**2-psumjet(2)**2-psumjet(3)**2
c--- check that argument of upcoming sqrt is not negative
        if (mjets < 0._dp) then
          wt=0._dp
          return 1
        endif
        mjets=sqrt(mjets)
      endif

      if (psumjet(4)-psumjet(3) == 0._dp) then
        wt=0._dp
        return 1
      endif
      ybar=(psumjet(4)+psumjet(3))/(psumjet(4)-psumjet(3))
c--- check that argument of upcoming log is not negative or infinite
      if (ybar <= 0._dp) then
        wt=0._dp
        return 1
      endif
      ybar=0.5_dp*log(ybar)

      ptsumjet2=psumjet(1)**2+psumjet(2)**2
      plstarsq=((sqrts**2-mv2-mjets**2)**2
     & -4._dp*(mjets**2*mv2+ptsumjet2*sqrts**2))/(4._dp*sqrts**2)
c--- check that argument of upcoming sqrt is not negative
      if (plstarsq < 0._dp) then
        wt=0._dp
        return 1
      endif
      plstar=sqrt(plstarsq)
      Estar=plstarsq+ptsumjet2+mjets**2
c--- check that argument of upcoming sqrt is not negative
      if (Estar < 0._dp) then
        wt=0._dp
        return 1
      endif
      Estar=sqrt(Estar)
      if (Estar-plstar == 0._dp) then
        wt=0._dp
        return 1
      endif
      y5starmax=(Estar+plstar)/(Estar-plstar)
c--- check that argument of upcoming log is not negative or infinite
      if (y5starmax <= 0._dp) then
        wt=0._dp
        return 1
      endif
      y5starmax=0.5_dp*log(y5starmax)
      y5starmin=-y5starmax

      etamax=ybar-y5starmin
      etamin=ybar-y5starmax
      dely=etamax-etamin
      ycm=etamin+r(3*njets+1)*dely
      sinhy=sinh(ycm)
      coshy=sqrt(1._dp+sinhy**2)

c--- now make the initial state momenta
      sumpst=ptsumjet2+(psumjet(3)*coshy-psumjet(4)*sinhy)**2
      q0st=mv2+sumpst
c--- check that argument of upcoming sqrt is not negative
      if (q0st < 0._dp) then
        wt=0._dp
        return 1
      endif
      q0st=sqrt(q0st)
      rshat=mjets**2+sumpst
c--- check that argument of upcoming sqrt is not negative
      if (rshat < 0._dp) then
        wt=0._dp
        return 1
      endif
      rshat=q0st+sqrt(rshat)
      pcm(4)=rshat*coshy
      pcm(3)=rshat*sinhy

      xx(1)=(pcm(4)+pcm(3))/sqrts
      xx(2)=(pcm(4)-pcm(3))/sqrts

      if   ((xx(1)*xx(2) > 1._dp)) then
c        write(6,*) 'gen_stop: xx(1)*xx(2),xx(1),xx(2)',
c     &   xx(1)*xx(2),xx(1),xx(2)
      endif

      if   ((xx(1) > 1._dp) .or. (xx(2) > 1._dp)) then
         wt=0._dp
         return 1
      endif

      wt=wt*dely
      do j=1,4
        Q(j)=pcm(j)-psumjet(j)
      enddo

      p(1,4)=-xx(1)*sqrts/2._dp
      p(1,3)=p(1,4)
      p(2,4)=-xx(2)*sqrts/2._dp
      p(2,3)=-p(2,4)

      wt=wt*rshat/(sqrts**2*q0st)

c--- If we're calculating top decay then generate the additional jet
c--- for the real contribution here, after the decay
      if ( ((kcase==kttdkay) .or. (kcase==ktdecay))
     &     .and. (kpart==kreal) ) then
        in=3*njets+2
        call phi1_2(r(in),r(in+1),r(in+2),r(in+3),Q,pbg,wp,wtwp,*999)
        in=in+4
          call phi3m(r(in),r(in+1),pbg,bm,g,mb,zip,wtbg,*999)
          call phi3m0(r(in+2),r(in+3),wp,nn,ep,wtepnn,*999)
            wt=wt0*wt*wtwp*wtbg*wtepnn/twopi
        do nu=1,4
          p(7,nu)=g(nu)
        enddo
      elseif (kcase == ktopanom .and. ipsgen == 1) then
        ! this is for radiation in production, when this routine is called with njets=2
        ! or for born topology pieces with njets=1
        call phi1_2m(0._dp,r(3*njets+2),r(3*njets+3),r(3*njets+4),zip,
     &  Q,bm,wp,wtwp,*999)
        call phi3m0(r(3*njets+5),r(3*njets+6),wp,nn,ep,wtepnn,*999)
        wt=wt0*wt*wtwp*wtepnn
      elseif (kcase == ktopanom .and. ipsgen == 2) then
        ! this is for radiation in decay, when this routine is called with njets=1
        in=3*njets+2
        call phi1_2(r(in),r(in+1),r(in+2),r(in+3),Q,pbg,wp,wtwp,*999)
        in=in+4
          call phi3m0(r(in),r(in+1),pbg,bm,g,wtbg,*999)
          call phi3m0(r(in+2),r(in+3),wp,nn,ep,wtepnn,*999)
            wt=wt0*wt*wtwp*wtbg*wtepnn/twopi
        do nu=1,4
          p(7,nu)=g(nu)
        enddo
      else
        call phi1_2m(mb,r(3*njets+2),r(3*njets+3),r(3*njets+4),zip,
     &  Q,bm,wp,wtwp,*999)
        call phi3m0(r(3*njets+5),r(3*njets+6),wp,nn,ep,wtepnn,*999)
        wt=wt0*wt*wtwp*wtepnn
      endif

      do nu=1,4
        p(3,nu)=nn(nu)
        p(4,nu)=ep(nu)
        p(5,nu)=bm(nu)
      enddo


      return

  999 wt=0._dp
      return 1

      end

