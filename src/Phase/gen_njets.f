!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gen_njets(r,njets,p,wt,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'mxdim.f'
      include 'limits.f'
      include 'nodecay.f'
      include 'reset.f'
      include 'hdecaymode.f'
      include 'masses.f'
      include 'breit.f'
      include 'kpart.f'
      include 'x1x2.f'
      include 'notag.f'
      include 'energy.f'
      include 'first.f'
      include 'taucut.f'
      include 'nproc.f'
      include 'jetcuts.f'
c---- generate phase space for 2-->2+n process
c---- with (34) being a vector boson and 5,..,4+n the jets
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign
c---- from physical values
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
c----
c---- if 'nodecay' is true, then the vector boson decay into massless
c---- particles is not included and 2 less integration variables
c---- are required
      real(dp):: r(mxdim),rdk1,rdk2
      real(dp):: p(mxpart,4),p3(4),p34(4),psumjet(4),pcm(4),Q(4)
      real(dp):: wt
      real(dp):: pt,etamax,etamin
      real(dp):: y,sinhy,coshy,phi,mv2,wtbw,mjets
      real(dp):: ybar,ptsumjet2,ycm,sumpst,q0st,rshat
      real(dp):: costh,sinth,dely,xjac
      real(dp):: ptjetmin_save
      real(dp):: plstar,estar,plstarsq,y5starmax,y5starmin,mf,beta
      integer:: j,nu,njets,ijet
      real(dp) :: etajetmax_use

      etajetmax_use = 1.e10_dp

      do nu=1,4
        do j=1,4+njets
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

c--- relax cuts for real and calculations with SCET above taucut
        if (origkpart == kresummed) then
            pt = r(ijet)**2*sqrts/2._dp
            xjac = pt*sqrts/2._dp*2._dp*r(ijet)
        elseif (((usescet .eqv. .false.) .and. (kpart == kvirt)) .or.
     &          ((usescet) .and. (abovecut .eqv. .false.) .and. (origkpart /= kn3lo))) then
            call genpt(r(ijet),ptjetmin,.true.,pt,xjac)
        else
           ! linear
            !pt = sqrts*r(ijet)/2._dp
            !xjac = r(ijet)*sqrts**2/4._dp
           ! sqrt
            !pt = r(ijet)**2*sqrts/2._dp
            !xjac = pt*sqrts/2._dp*2._dp*r(ijet)
           ! 1/pT
            pt = sqrts*r(ijet)/(sqrts+2*r(ijet)-sqrts*r(ijet))
            xjac = r(ijet)*sqrts**3 / (2*r(ijet) - sqrts*(r(ijet)-1))**3
           ! -log(pT)
            !pt = exp(log(sqrts/2._dp) - (1-r(ijet))/r(ijet))
            !xjac = exp(2._dp*log(sqrts/2._dp) - 2._dp*(1-r(ijet))/r(ijet))/r(ijet)**2
        !elseif ((kpart .eq. kreal) .or. (usescet .and. abovecut)) then
          !call genpt(r(ijet),ptjetmin_use,.false.,pt,xjac)
        !else
          !call genpt(r(ijet),ptjetmin_use,.true.,pt,xjac)
        endif
        wt=wt*xjac
        etamax=sqrts/2._dp/pt
        if (etamax**2 <= 1._dp) then
!          write(6,*) 'etamax**2 <= 1._dp in gen_njets.f',etamax**2
          wt=0._dp
          return 1
        endif
        etamax=log(etamax+sqrt(etamax**2-1._dp))

c--- only generate up to maximum jet rapidity, if not using SCET
c        if ((usescet .eqv. .false.) .and. (new_pspace .eqv. .false.)) then
           etamax=min(etamax,etajetmax_use)
c        endif
        y=etamax*(2._dp*r(njets+ijet)-1._dp)
        wt=wt*2._dp*etamax

        sinhy=sinh(y)
        coshy=sqrt(1._dp+sinhy**2)

        p(4+ijet,4)=pt*coshy

        phi=2._dp*pi*r(2*njets+ijet)
        wt=wt*2._dp*pi

        p(4+ijet,1)=pt*cos(phi)
        p(4+ijet,2)=pt*sin(phi)
        p(4+ijet,3)=pt*sinhy

        do nu=1,4
          psumjet(nu)=psumjet(nu)+p(4+ijet,nu)
        enddo
      enddo

c--- now generate Breit-Wigner
      call breitw(r(3*njets+1),wsqmin,wsqmax,mass3,width3,mv2,wtbw)
      wt=wt*wtbw/2._dp/pi
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

c--- check that upcoming denominator does not vanish
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
c--- check that upcoming denominator does not vanish
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
      ycm=etamin+r(3*njets+2)*dely
      sinhy=sinh(ycm)
      coshy=sqrt(1._dp+sinhy**2)

c--- now make the initial state momenta
c      write(6,*) ptsumjet2,psumjet(3),coshy,psumjet(4),sinhy
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

c     if   ((xx(1)*xx(2) > 1._dp) .and. (xxerror .eqv. .false.)) then
c       xxerror=.true.
c        write(6,*) 'gen_njets: xx(1)*xx(2),xx(1),xx(2)',
c     &   xx(1)*xx(2),xx(1),xx(2)
c     endif

      if   ((xx(1) > 1._dp) .or. (xx(2) > 1._dp)) then
c     & .or. (xx(1) < xmin ) .or. (xx(2) < xmin )) then
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

c      if (p(1,4) == p(1,4)) then
c        continue
c      else
c        write(6,*) 'q0st,ptsumjet2,psumjet',q0st,ptsumjet2,psumjet
c        write(6,*) 'ycm,ybar,coshy,sinhy',ycm,ybar,coshy,sinhy
c        stop
c      endif

      wt=wt*rshat/(sqrts**2*q0st)

c--- dummy values if there's no decay
      if (nodecay) then
        rdk1=0.5_dp
        rdk2=0.5_dp
      else
        rdk1=r(3*njets+3)
        rdk2=r(3*njets+4)
      endif
      if (hdecaymode == 'tlta') then
        mf=mtau
      elseif (hdecaymode == 'bqba') then
        mf=mb
      else
        mf=0._dp
      endif

c--- decay boson into leptons, in boson rest frame
      if (mv2 > 4._dp*mf**2) then
        beta=sqrt(1._dp-4._dp*mf**2/mv2)
      else
        beta=0._dp
      endif
      costh=2._dp*rdk1-1._dp
      sinth=sqrt(1._dp-costh**2)
      phi=2._dp*pi*rdk2
      p34(4)=sqrt(mv2)/2._dp
      p34(1)=beta*p34(4)*sinth*cos(phi)
      p34(2)=beta*p34(4)*sinth*sin(phi)
      p34(3)=beta*p34(4)*costh

c--- boost into lab frame
      call boost(sqrt(mv2),Q,p34,p3)
      do j=1,4
      p(3,j)=p3(j)
      p(4,j)=Q(j)-p(3,j)
      enddo

      wt=wt*beta/8._dp/pi

      return
      end











