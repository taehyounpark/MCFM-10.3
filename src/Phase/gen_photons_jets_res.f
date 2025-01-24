!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gen_photons_jets_res(r,nphots,njets,p,wt,*)
      implicit none
      include 'types.f'
      ! basically an unconstrained version of gen_photons_jets
c----  WARNING: although mostly written in a generic manner, this
c----           routine has been tailored for direct photon (dirgam)
c----           and photon+two jets (gamjet) production

c---- generate phase space for 2-->nphots+njets process
c----   and 3,..,2+nphots the photons
c----   and 2+nphots,..,2+nphots+njets the jets
c----
c---- This routine is adapted from gen_Vphotons_jets.f
c----
c---- r(mxdim),p1(4),p2(4) are inputs reversed in sign
c---- from physical values
c---- phase space for -p1-p2 --> p3+p4+p5+p6
c---- with all 2 pi's (ie 1/(2*pi)^(4+2n))
c----
      include 'constants.f'
      include 'mxpart.f'
      include 'mxdim.f'
      include 'leptcuts.f'
      include 'x1x2.f'
      include 'energy.f'
      real(dp):: r(mxdim),p(mxpart,4),psumjet(4),wt
      real(dp):: xpt,etamax
      real(dp):: y,sinhy,coshy,phi,xjac,tmp(mxpart,4)
c      real(dp):: ptjetmin,etajetmin,etajetmax,pbreak
      real(dp):: ptmin_part,etamax_part
      integer:: j,nu,nphots,njets,nphotsjets,ijet,icount
      integer,save:: perm=0
      integer,parameter:: perm2_3(0:1)=(/3,4/)
      integer,parameter:: perm2_4(0:1)=(/4,3/)
      integer,parameter:: perm3_3(0:5)=(/3,3,4,4,5,5/)
      integer,parameter:: perm3_4(0:5)=(/4,5,3,5,3,4/)
      integer,parameter:: perm3_5(0:5)=(/5,4,5,3,4,3/)
      integer,parameter:: perm4_3(0:23)=(/3,3,3,3,3,3,4,4,4,4,4,4,
     & 5,5,5,5,5,5,6,6,6,6,6,6/)
      integer,parameter:: perm4_4(0:23)=(/4,4,5,5,6,6,3,3,5,5,6,6,
     & 3,3,4,4,6,6,3,3,4,4,5,5/)
      integer,parameter:: perm4_5(0:23)=(/5,6,4,6,4,5,5,6,3,6,3,5,
     & 4,6,3,6,3,4,4,5,3,5,3,4/)
      integer,parameter:: perm4_6(0:23)=(/6,5,6,4,5,4,6,5,6,3,5,3,
     & 6,4,6,3,4,3,5,4,5,3,4,3/)
c      save ptjetmin,etajetmin,etajetmax,pbreak
comp threadprivate(ptjetmin,etajetmin,etajetmax,pbreak,perm)

c---  total number of photons and jets
      nphotsjets=nphots+njets

      do nu=1,4
        do j=1,4+nphotsjets
          p(j,nu)=0._dp
        enddo
        psumjet(nu)=0._dp
      enddo

      wt=2._dp*pi


      icount=0
c--- loop over all the jets and photons except the last one
      do ijet=1,nphotsjets-1

c--- generate the pt of jet number ijet
c--- rapidity limited by E=pT*coshy
      wt=wt/16._dp/pi**3

      ptmin_part = 0d0
      etamax_part = 20d0
      if (ijet <= nphots) then
          if (nphots == 2) then
            ptmin_part=min(gammptmin,gammpt2)
          else
            write(6,*) 'Unexpected # photons in gen_photons_jets: ',nphots
            error stop
          endif
      endif

      if (ptmin_part > 0._dp) then
          call genpt(r(icount+1),ptmin_part,.false.,xpt,xjac)
      else
          xpt = r(icount+1)**2*sqrts/2._dp
          xjac = xpt*sqrts/2._dp*2._dp*r(icount+1)
      endif

      wt=wt*xjac

      etamax=sqrts/2._dp/xpt
      if (etamax**2 <= 1._dp) then
        !write(6,*) 'etamax**2 <= 1._dp in gen_photons_jets.f',etamax**2
        wt=0._dp
        return 1
      endif
      etamax=log(etamax+sqrt(etamax**2-1._dp))

      etamax=min(etamax,etamax_part)
      y=etamax*(2._dp*r(icount+2)-1._dp)
      wt=wt*2._dp*etamax

      sinhy=sinh(y)
      coshy=sqrt(1._dp+sinhy**2)

      p(2+ijet,4)=xpt*coshy

      phi=2._dp*pi*r(icount+3)
      wt=wt*2._dp*pi

      p(2+ijet,1)=xpt*cos(phi)
      p(2+ijet,2)=xpt*sin(phi)
      p(2+ijet,3)=xpt*sinhy

      do nu=1,4
        psumjet(nu)=psumjet(nu)+p(2+ijet,nu)
      enddo

      icount=icount+3

c--- end of loop over photons and jets
      enddo

c---- have generated nphotsjets-1 massless momenta, last pt is fixed but rapidity is unconstrained
      xpt=sqrt(psumjet(1)**2+psumjet(2)**2)
      etamax=sqrts/2._dp/xpt
      if (etamax**2 <= 1._dp) then
c         write(6,*) 'etamax**2 <= 1._dp in gen_photons_jets.f',etamax**2
         wt=0._dp
         return 1
      endif
      etamax=log(etamax+sqrt(etamax**2-1._dp))

      etamax = min(etamax,15d0)

      y=etamax*(2._dp*r(3*nphotsjets-2)-1._dp)
      wt=wt*2._dp*etamax

      sinhy=sinh(y)
      coshy=sqrt(1._dp+sinhy**2)

      p(2+nphotsjets,4)=xpt*coshy

      p(2+nphotsjets,1)=-psumjet(1)
      p(2+nphotsjets,2)=-psumjet(2)
      p(2+nphotsjets,3)=xpt*sinhy

      do nu=1,4
      psumjet(nu)=psumjet(nu)+p(2+nphotsjets,nu)
      enddo

c--- now make the initial state momenta
      xx(1)=(psumjet(4)+psumjet(3))/sqrts
      xx(2)=(psumjet(4)-psumjet(3))/sqrts
      if (psumjet(3) == 0._dp) then
        wt=0._dp
        return 1
      endif
      if (abs(psumjet(4)/psumjet(3)-1._dp) < 1.e-12_dp) then
        wt=0._dp
        return 1
      endif

      if   ((xx(1) > 1._dp) .or. (xx(2) > 1._dp)) then
         wt=0._dp
         return 1
      endif

      p(1,4)=-0.5_dp*xx(1)*sqrts
      p(1,1)=0._dp
      p(1,2)=0._dp
      p(1,3)=-0.5_dp*xx(1)*sqrts

      p(2,4)=-0.5_dp*xx(2)*sqrts
      p(2,1)=0._dp
      p(2,2)=0._dp
      p(2,3)=+0.5_dp*xx(2)*sqrts

      wt=wt/sqrts**2

      do j=3+nphotsjets,mxpart
      do nu=1,4
      p(j,nu)=0._dp
      enddo
      enddo

c--- symmetrize phase-space by permuting photon momenta appropriately
      if     (nphots == 2) then
        perm=int(r(3*nphotsjets-1)*2._dp)
        if (perm < 0) perm=0
        if (perm > 1) perm=1
        tmp(3,:)=p(3,:)
        tmp(4,:)=p(4,:)
        p(3,:)=tmp(perm2_3(perm),:)
        p(4,:)=tmp(perm2_4(perm),:)
c        perm=mod(perm+1,2)
      elseif (nphots == 3) then
        perm=int(r(3*nphotsjets-1)*6._dp)
        if (perm < 0) perm=0
        if (perm > 5) perm=5
        tmp(3,:)=p(3,:)
        tmp(4,:)=p(4,:)
        tmp(5,:)=p(5,:)
        p(3,:)=tmp(perm3_3(perm),:)
        p(4,:)=tmp(perm3_4(perm),:)
        p(5,:)=tmp(perm3_5(perm),:)
c        perm=mod(perm+1,6)
      elseif (nphots==4) then
        perm=int(r(3*nphotsjets-1)*24._dp)
        if (perm < 0) perm=0
        if (perm > 23) perm=23
        tmp(3,:)=p(3,:)
        tmp(4,:)=p(4,:)
        tmp(5,:)=p(5,:)
        tmp(6,:)=p(6,:)
        p(3,:)=tmp(perm4_3(perm),:)
        p(4,:)=tmp(perm4_4(perm),:)
        p(5,:)=tmp(perm4_5(perm),:)
        p(6,:)=tmp(perm4_6(perm),:)
c        perm=mod(perm+1,24)
      elseif (nphots >= 5) then
        write(6,*) 'Symmetrize gen_photons_jets.f before proceeding'
        write(6,*) 'nphots = ',nphots
        stop
      endif

      return
      end











