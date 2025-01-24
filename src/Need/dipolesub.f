!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c***********************************************************************
c     Author: J. M. Campbell                                           *
c     August, 1999                                                     *
c     Calculates the nj-jet subtraction term corresponding to dipole   *
c     nd with momentum p and dipole kinematics (ip,jp) wrt kp          *
c     Automatically chooses dipole kind                                *
c     Returns the dipoles in sub,subv and matrix elements in msq,msqv  *
c     nd labels the dipole configurations                              *
c     ip labels the emitter parton                                     *
c     jp labels the emitted parton                                     *
c     kp labels the spectator parton                                   *
c     subr_born is the subroutine which call the born process          *
c     subr_corr is the subroutine which call the born process dotted   *
c     with vec for an emitted gluon only                               *
c***********************************************************************

      subroutine dips(nd,p,ip,jp,kp,sub,subv,msq,msqv,
     & subr_born,subr_corr)
        use singletop2_scale_m
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'kprocess.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      include 'incldip.f'
      include 'msq_cs.f'
      include 'msqv_cs.f'
      include 'msq_struc.f'
      real(dp):: p(mxpart,4),ptrans(mxpart,4),sub(4),subv,vecsq
      real(dp):: x,omx,z,omz,y,omy,u,omu,sij,sik,sjk,dot,vec(4)
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),vtilde
      integer:: nd,ip,jp,kp,nu,j,k,ipt
c--      logical:: includedipole
      external subr_born,subr_corr

c---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0._dp
      enddo
      subv=0._dp
      msq = 0._dp
      msqv = 0._dp
      msq_cs = 0._dp
      msqv_cs = 0._dp
      msq_struc = 0._dp
      msq_strucv = 0._dp
      incldip(nd)=.true.

      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp)
      sjk=two*dot(p,jp,kp)

c**********************************************************************
c************************** INITIAL-INITIAL ***************************
c**********************************************************************
      if ((ip <= 2) .and. (kp <= 2)) then
        omx=-(sij+sjk)/sik
        x=one-omx
        vtilde=sij/sik

c---Modification so that only close to singular subtracted
        if (-vtilde > aii) then
           incldip(nd)=.false.
           return
        endif

        call transform(p,ptrans,x,ip,jp,kp)

        if ((kcase==ktwo_ew) .and.
     &      (ip == 1) .and. (jp < 5)) then ! (13) or (14) singularity
c------ swap particles 3 and 4
          vec(:)=ptrans(4,:)
          ptrans(4,:)=ptrans(3,:)
          ptrans(3,:)=vec(:)
        endif

        if ((kcase==kWgajew) .and. (jp == 5)) then
c------ swap particles 5 and 6 to get photon and jet in correct positions
          vec(:)=ptrans(6,:)
          ptrans(6,:)=ptrans(5,:)
          ptrans(5,:)=vec(:)
        endif

        call storeptilde(nd,ptrans)

c-- Check to see if this dipole will be included
c        incldip(nd)=includedipole(nd,ptrans)
c--if not return
c        if (incldip(nd) .eqv. .false.) return

        vecsq=-sij*sjk/sik
        do nu=1,4
          vec(nu)=(p(jp,nu)-vtilde*p(kp,nu))/sqrt(-vecsq)
        enddo

c     write (*,*) "initial initial"
c--- if using a dynamic scale, set that scale with dipole kinematics
      if (dynamicscale) then
        if (kcase == ktopanom .and. use_DDIS) then
          call singletop2_set_dipscale(nd, ptrans, gsq)
        else
          call scaleset(initscale,initfacscale,ptrans)
          dipscale(nd)=facscale
        endif
      endif

c--- for "gamgam" process, gvec contribution is not written in
c--- a canonical way; instead, pass emitted vector directly as "vec"
        if ((kcase==kgamgam) .or. (kcase==kgg2gam)) then
          do nu=1,4
            vec(nu)=p(jp,nu)
          enddo
      endif

        call subr_born(ptrans,msq)
        call subr_corr(ptrans,vec,ip,msqv)

        sub(qq)=-gsq/x/sij*(two/omx-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/sij*(x/omx+x*omx)
        subv   =-4._dp*gsq/x/sij*omx/x

c**********************************************************************
c************************** INITIAL-FINAL *****************************
c**********************************************************************
      elseif ((ip <= 2) .and. (kp > 2)) then
        u=sij/(sij+sik)

        omx=-sjk/(sij+sik)
        x=one-omx
        omu=sik/(sij+sik)
c---npart is the number of particles in the final state
c---transform the momenta so that only the first npart+1 are filled
        call transform(p,ptrans,x,ip,jp,kp)

        if ((kcase==kWgajew) .and. (jp == 5)) then
c------ swap particles 5 and 6 to get photon and jet in correct positions
          vec(:)=ptrans(6,:)
          ptrans(6,:)=ptrans(5,:)
          ptrans(5,:)=vec(:)
        endif

        call storeptilde(nd,ptrans)

c-- Check to see if this dipole will be included
c        incldip(nd)=includedipole(nd,ptrans)
c-- if not return
c        if (incldip(nd) .eqv. .false.) return

c     write (*,*) "initial final"
c--- if using a dynamic scale, set that scale with dipole kinematics
      if (dynamicscale) then
        if (kcase == ktopanom .and. use_DDIS) then
          call singletop2_set_dipscale(nd, ptrans, gsq)
        else
          call scaleset(initscale,initfacscale,ptrans)
          dipscale(nd)=facscale
        endif
      endif

c--- Calculate the matrix element now because it might be needed
c--- in the final-initial segment, regardless of whether or not the
c--- alfa cut fails here
        call subr_born(ptrans,msq)
c---Modification so that only close to singular subtracted
c---Do not set incldip because initial-final can fail
c---but final initial needs still to be tested

        if ((kcase==kH_1jet) .and. (ip==1) .and. (jp==6)) then
c--- do nothing
        else
          if (u > aif) return
        endif

        do nu=1,4
           vec(nu)=(p(jp,nu)/u-p(kp,nu)/omu)/sqrt(sjk)
        enddo

        call subr_corr(ptrans,vec,ip,msqv)
        sub(qq)=-gsq/x/sij*(two/(omx+u)-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/sij*(one/(omx+u)-one+x*omx)
        subv   =-4._dp*gsq/x/sij*(omx/x*u*(one-u))

c**********************************************************************
c************************** FINAL-INITIAL *****************************
c**********************************************************************
      elseif ((ip > 2) .and. (kp <= 2)) then
c-- Check to see if this dipole will be included - should have been
c-- already determined at this point in the initial-final phase
c        if (incldip(nd) .eqv. .false.) return

c--- note, here we assume that msq kinematics are already taken care of
c--- for msq, although msqv must be recalculated each time
        omx=-sij/(sjk+sik)
c---Modification so that only close to singular subtracted
        if (omx > afi) return

        x=one-omx
        z=sik/(sik+sjk)
        omz=sjk/(sik+sjk)

        call transform(p,ptrans,x,ip,jp,kp)

        if ((kcase==kWgajew) .and. (jp == 5)) then
c------ swap particles 5 and 6 to get photon and jet in correct positions
          vec(:)=ptrans(6,:)
          ptrans(6,:)=ptrans(5,:)
          ptrans(5,:)=vec(:)
        endif

        call storeptilde(nd,ptrans)

        do nu=1,4
          vec(nu)=(z*p(ip,nu)-omz*p(jp,nu))/sqrt(sij)
        enddo
c---call msqv again because vec has changed
        do j=1,mxpart
        do k=1,4
          ptrans(j,k)=ptilde(nd,j,k)
        enddo
        enddo

c       write (*,*) "final initial"
c--- if using a dynamic scale, set that scale with dipole kinematics
      if (dynamicscale) then
        if (kcase == ktopanom .and. use_DDIS) then
          call singletop2_set_dipscale(nd, ptrans, gsq)
        else
          call scaleset(initscale,initfacscale,ptrans)
          dipscale(nd)=facscale
        endif
      endif

c--- do something special if we're doing W+2,Z+2jet (jp  /=  7)
        if ((jp  /=  7) .and. (kcase /= kHWWjet)
     &      .and. (kcase /= kHZZjet)
     &      .and. (kcase /= kWH1jet)
     &      .and. (kcase /= kZH1jet)
     &      .and. (kcase /= kqq_HWW)
     &      .and. (kcase /= kqq_HZZ)) then
          if (ip < 7) then
c ie for cases 56_i,65_i
            ipt=5
          else
c ie for cases 76_i,75_i
            ipt=6
          endif
        else
c ie for cases 57_i,67_i
          ipt=ip
        endif

c--- do something special for H(->Za)+jet
        if (kcase == kHi_Zaj) ipt=6

c--- do something special for HWW/HZZ+2 jets
        if ((kcase==kHWW2jt) .or. (kcase==kHZZ2jt)) then
          if (jp  /=  9) then
            if (ip < 9) then
c ie for cases 78_i,87_i
              ipt=7
            else
c ie for cases 98_i,97_i
              ipt=8
            endif
          else
c ie for cases 79_i,89_i
            ipt=ip
          endif
        endif

c--- special case for VV+jet
        if ((kcase==kWW_jet) .or. (kcase==kWZ_jet) .or. (kcase==kZZ_jet)) ipt=7
c--- special case for Wgam+jet
        if ((kcase==kWgajet) .or. (kcase==kWgajew)) ipt=6

c--- do something special for direct photon production
        if (kcase==kdirgam) ipt=4

        call subr_born(ptrans,msq)
        call subr_corr(ptrans,vec,ipt,msqv)

        sub(qq)=+gsq/x/sij*(two/(omz+omx)-one-z)
        sub(gq)=+gsq/x/sij
        sub(gg)=+2._dp*gsq/x/sij*(one/(omz+omx)+one/(z+omx)-two)
        subv   =+4._dp*gsq/x/sij

c**********************************************************************
c*************************** FINAL-FINAL ******************************
c**********************************************************************
      elseif ((ip > 2) .and. (kp > 2)) then
c------Eq-(5.2)
        y=sij/(sij+sjk+sik)

c---Modification so that only close to singular subtracted
        if (y > aff) then
          incldip(nd)=.false.
          return
        endif

        z=sik/(sjk+sik)
        omz=one-z
        omy=one-y
c---calculate the ptrans-momenta

        call transform(p,ptrans,y,ip,jp,kp)

        if ((kcase==kWgajew) .and. (jp == 5)) then
c------ swap particles 5 and 6 to get photon and jet in correct positions
          vec(:)=ptrans(6,:)
          ptrans(6,:)=ptrans(5,:)
          ptrans(5,:)=vec(:)
        endif

        call storeptilde(nd,ptrans)

c-- Check to see if this dipole will be included
c        incldip(nd)=includedipole(nd,ptrans)
c        if (incldip(nd) .eqv. .false.) return

        do nu=1,4
          vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
        enddo

c        write (*,*) "final final"
c--- if using a dynamic scale, set that scale with dipole kinematics
        if (dynamicscale) then
         if (kcase == ktopanom .and. use_DDIS) then
           call singletop2_set_dipscale(nd, ptrans, gsq)
         else
           call scaleset(initscale,initfacscale,ptrans)
           dipscale(nd)=facscale
         endif
        endif

        call subr_born(ptrans,msq)
        if (kcase==kepem3j) then
          ipt=5
        else
          if (ip < kp) then
            ipt=5
          else
            ipt=6
          endif
        endif

c--- do something special for HWW/HZZ+2 jets
        if ((kcase==kHWW2jt) .or. (kcase==kHZZ2jt)) then
          if (ip < kp) then
            ipt=7
          else
            ipt=8
          endif
      endif

c--- we only need this for the (473) (476) config
c--- for the identical qa qa final state
        if (kcase == kZgajet .and. ip == 4) then
            ipt = 4
        endif

        call subr_corr(ptrans,vec,ipt,msqv)

        sub(qq)=gsq/sij*(two/(one-z*omy)-one-z)
        sub(gq)=gsq/sij
        sub(gg)=gsq/sij*(two/(one-z*omy)+two/(one-omz*omy)-four)
        subv   =+4._dp*gsq/sij/sij

      endif

      return
      end

