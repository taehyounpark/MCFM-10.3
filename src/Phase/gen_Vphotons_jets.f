!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gen_Vphotons_jets(r,nphots,njets,p,wt,*)
      implicit none
      include 'types.f'
c---- generate phase space for 2-->2+nphots+njets process
c----   with (34) being a vector boson
c----   and 5,..,4+nphots the photons
c----   and 5+nphots,..,4+nphots+njets the jets
c----
c---- This routine is adapted from gen_njets.f
c----
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign
c---- from physical values
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
c----
c---- if 'nodecay' is true, then the vector boson decay into massless
c---- particles is not included and 2 less integration variables
c---- are required

      include 'constants.f'
      include 'mxpart.f'
      include 'mxdim.f'
      include 'limits.f'
      include 'nodecay.f'
      include 'leptcuts.f'
      include 'reset.f'
      include 'breit.f'
      include 'kpart.f'
      include 'x1x2.f'
      include 'energy.f'
      include 'first.f'
      include 'jetcuts.f'
      real(dp):: r(mxdim),rdk1,rdk2
      real(dp):: p(mxpart,4),p3(4),p34(4),psumjet(4),pcm(4),Q(4),wt
      real(dp):: pt,etamax,etamin
      real(dp):: y,sinhy,coshy,phi,mv2,wtbw,mjets
      real(dp):: ybar,ptsumjet2,ycm,sumpst,q0st,rshat
      real(dp):: costh,sinth,dely,xjac
      real(dp), save :: ptjetmin_use,etajetmax_use,pbreak
!$omp threadprivate(ptjetmin_use,etajetmax_use,pbreak)
      real(dp):: ptmin_part,etamax_part
      real(dp):: plstar,estar,plstarsq,y5starmax,y5starmin
      integer:: j,nu,nphots,njets,nphotsjets,ijet

      if (first .or. reset) then
        first=.false.
        reset=.false.

        ptjetmin_use = ptjetmin
        etajetmax_use = etajetmax

        if (kpart==kreal) then
c--- if we're generating phase space for real emissions, then we need
c--- to produce partons spanning the whole phase space pt>0,eta<10;
c--- in this case, pbreak=ptjetmin simply means that we
c--- generate pt approx. 1/x for pt > pbreak and
c--- pt approx. uniformly for pt < pbreak
c          pbreak=ptjetmin
c          ptjetmin=0._dp
          etajetmax_use = 20._dp
        else
c--- for lord and virt, the partons produced here can be generated
c--- right up to the jet cut boundaries and there is no need for pbreak
          pbreak=0._dp
        endif
c--- in case this routine is used for very small values of ptjetmin
        if ((ptjetmin_use < 5._dp) .and. (kpart /= kreal)) pbreak=5._dp
c--- for processes in which it is safe to jet ptmin to zero at NLO
      if ((kpart==kreal) .and. (pbreak < 1.e-8_dp)) pbreak=5._dp
      endif

c--- total number of photons and jets
      nphotsjets=nphots+njets

      p(:,:)=0._dp
      Q(:)=0._dp
      psumjet(:)=0._dp
      pcm(:)=0._dp

      wt=2._dp*pi

      do ijet=1,nphotsjets
c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
        wt=wt/16._dp/pi**3

        if (ijet <= nphots) then
          if     (nphots == 3) then
            ptmin_part=min(gammptmin,gammpt2,gammpt3)
          elseif (nphots == 2) then
            ptmin_part=min(gammptmin,gammpt2)
          elseif (nphots == 1) then
            ptmin_part=gammptmin
          else
            write(6,*) 'Unexpected # photons in gen_Vphotons_jets: ',
     &                  nphots
            stop
          endif
          etamax_part=gammrapmax
c          pbreak_part=0._dp
          if (kpart==kreal) then
c--- cannot generate exactly to match, since dipoles transform photon
          etamax_part=20._dp
c          ptmin_part=0._dp
c          pbreak_part=gammpt
          endif
        else
          ptmin_part=ptjetmin_use
          etamax_part=etajetmax_use
c          pbreak_part=pbreak
          if (kpart==kreal) then
c--- cannot generate exactly to match, since dipoles transform jet
          etamax_part=20._dp
c          ptmin_part=0._dp
c          pbreak_part=gammpt
          endif
        endif

c        if ((flatreal) .and. (kpart==kreal)) then
cc--- generate flat pt for the real contribution
c          pt=r(ijet)*((sqrts/2._dp)-ptmin_part)+ptmin_part
c          wt=wt*((sqrts/2._dp)-ptmin_part)*pt
c        else
cc--- favour small pt region
c          hmin=1._dp/sqrt((sqrts/2._dp)**2+pbreak_part**2)
c          hmax=1._dp/sqrt(ptmin_part**2+pbreak_part**2)
c          delh=hmax-hmin
c          h=hmin+r(ijet)*delh
c          pt=sqrt(1._dp/h**2-pbreak_part**2)
c          wt=wt*delh/h**3
c        endif

        if (kpart==kreal) then
          call genpt(r(ijet),ptmin_part,.false.,pt,xjac)
        else
          call genpt(r(ijet),ptmin_part,.true.,pt,xjac)
        endif
        wt=wt*xjac

        etamax=sqrts/2._dp/pt
        if (etamax**2 <= 1._dp) then
c          write(6,*) 'etamax**2 <= 1._dp in gen_Vphotons_jets.f',etamax**2
          wt=0._dp
          return 1
        endif
        etamax=log(etamax+sqrt(etamax**2-1._dp))

        etamax=min(etamax,etamax_part)
        y=etamax*(2._dp*r(nphotsjets+ijet)-1._dp)
        wt=wt*2._dp*etamax

        sinhy=sinh(y)
        coshy=sqrt(1._dp+sinhy**2)

        p(4+ijet,4)=pt*coshy

        phi=2._dp*pi*r(2*nphotsjets+ijet)
        wt=wt*2._dp*pi

        p(4+ijet,1)=pt*cos(phi)
        p(4+ijet,2)=pt*sin(phi)
        p(4+ijet,3)=pt*sinhy

        do nu=1,4
          psumjet(nu)=psumjet(nu)+p(4+ijet,nu)
        enddo
      enddo


c      if (kcase==kW_2gam) then
cc--- this modified version of the Breit-Wigner form allows for some fraction of
cc--- events to be generated more efficiently below mflatsq
cc--- [currently used only for W+2 photons]
c        if (gammpt > mass3) then
c           mflatsq=wsqmin
c        else
c          mflatsq=(mass3-gammpt/2._dp)**2
c          if (mflatsq < wsqmin) mflatsq=wsqmin
c        endif
c        call breitw_mod(r(3*nphotsjets+1),wsqmin,wsqmax,
c     &                  mflatsq,mass3,width3,mv2,wtbw)
c      else
c--- [currently used only for everything else]
c--- now generate Breit-Wigner
        call breitw(r(3*nphotsjets+1),wsqmin,wsqmax,
     &              mass3,width3,mv2,wtbw)
c        write(6,*) sqrt(mv2),wtbw
c      endif
      wt=wt*wtbw/2._dp/pi
c--- for one jet, mjets must be exactly zero
      if (nphotsjets == 1) then
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
      if (plstarsq <= 0._dp) then
        wt=0._dp
        return 1
      endif
      plstar=sqrt(plstarsq)
      Estar=sqrt(plstarsq+ptsumjet2+mjets**2)
c--- check that upcoming denominator does not vanish
      if (Estar-plstar == 0._dp) then
        wt=0._dp
        return 1
      endif
      y5starmax=0.5_dp*log((Estar+plstar)/(Estar-plstar))
      y5starmin=-y5starmax

      etamax=ybar-y5starmin
      etamin=ybar-y5starmax
      dely=etamax-etamin
      ycm=etamin+r(3*nphotsjets+2)*dely
      sinhy=sinh(ycm)
      coshy=sqrt(1._dp+sinhy**2)

c--- now make the initial state momenta
      sumpst=ptsumjet2+(psumjet(3)*coshy-psumjet(4)*sinhy)**2
      q0st=sqrt(mv2+sumpst)
      rshat=q0st+sqrt(mjets**2+sumpst)
      pcm(4)=rshat*coshy
      pcm(3)=rshat*sinhy

      xx(1)=(pcm(4)+pcm(3))/sqrts
      xx(2)=(pcm(4)-pcm(3))/sqrts

!     if   ((xx(1)*xx(2) > 1._dp) .and. (xxerror .eqv. .false.)) then
!       xxerror=.true.
!        write(6,*) 'gen_Vphotons_jets: xx(1)*xx(2),xx(1),xx(2)',
!     &   xx(1)*xx(2),xx(1),xx(2)
!     endif

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

c--- dummy values if there's no decay
      if (nodecay) then
        rdk1=0.5_dp
        rdk2=0.5_dp
      else
        rdk1=r(3*nphotsjets+3)
        rdk2=r(3*nphotsjets+4)
      endif

c--- decay boson into leptons, in boson rest frame
      costh=2._dp*rdk1-1._dp
      sinth=sqrt(1._dp-costh**2)
      phi=2._dp*pi*rdk2
      p34(4)=sqrt(mv2)/2._dp
      p34(1)=p34(4)*sinth*cos(phi)
      p34(2)=p34(4)*sinth*sin(phi)
      p34(3)=p34(4)*costh

c--- boost into lab frame
      call boost(sqrt(mv2),Q,p34,p3)
      do j=1,4
      p(3,j)=p3(j)
      p(4,j)=Q(j)-p(3,j)
      enddo

c--- veto PS regions for W+2 photons
c      if (kcase==kW_2gam) then
c        if (ipsgen == 1) then
c          s345=2._dp*(dot(p,3,4)+dot(p,3,5)+dot(p,4,5))
c          if (abs(sqrt(s345)-mass3) < 5._dp*width3) return 1
c          s346=2._dp*(dot(p,3,4)+dot(p,3,6)+dot(p,4,6))
c          if (abs(sqrt(s346)-mass3) < 5._dp*width3) return 1
c        endif
c      endif

      wt=wt/8._dp/pi

      return
      end











