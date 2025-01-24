!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine dips_mass(nd,p,ip,jp,kp,sub,subv,msq,msqv,
     & subr_born,subr_corr)
          use singletop2_scale_m, only: singletop2_set_dipscale
          use singletop2_nnlo_vars, only: currentContrib
          use singletop_interf_hxd, only: extend_trans_st => extend_trans
      implicit none
      include 'types.f'

c***********************************************************************
c     Author: Keith Ellis                                              *
c     June 2002                                                        *
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
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'kprocess.f'
      include 'alfacut.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      include 'breit.f'
      include 'incldip.f'
      include 'debug.f'
      include 'ipsgen.f'
      include 'nproc.f'
      real(dp):: p(mxpart,4),ptrans(mxpart,4),sub(4),subv,vecsq,
     & x,omx,z,omz,y,omy,u,omu,pij,pik,pjk,dot,q(4),qsq,qij(4),qijsq,
     & vec(4),root,vtilde,pold(mxpart,4),pext(mxpart,4)
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),zp,zm
      real(dp):: mksq,misq,mjsq,mijsq,muisq,mujsq,muksq,
     & muijsq,kappa,vijk,vtijk,viji,ztmi,ztmj,muk,mqsq,subv_gg,subv_gq
      real(dp):: yp
      integer:: nd,ip,jp,kp,nu,j,jproc,ipt
c--- common block to handle the possibility of multiple definitions
c--- of subv in the final-final section
      common/subv_ff/subv_gg,subv_gq
      external subr_born,subr_corr
      parameter(kappa=zip)
      logical, save :: first = .TRUE.
!$omp threadprivate(first,/subv_ff/)

      if (first) then
        first=.false.
c        if (debug) write(6,*) 'dips_mass:mass2',mass2
      endif
      mqsq=mass2**2

      do nu=1,4
      do j=1,mxpart
      ptrans(j,nu)=zip
      pold(j,nu)=p(j,nu)
      enddo
      enddo

      if (nproc == 1610 .and. currentContrib == 6 .and. ipsgen == 8) then
          p(3,:)=pold(3,:)+pold(4,:)+pold(5,:)+pold(7,:)
          p(4,:)=pold(6,:)
          p(5,:)=pold(8,:)
          p(7,:)=zip
      elseif ((kcase==kt_bbar) .or. (kcase==kbq_tpq)) then
c--- if we're doing single-top, reduce # of momenta from 7 to 5
        do nu=1,4
          p(3,nu)=pold(3,nu)+pold(4,nu)+pold(5,nu)
          p(4,nu)=pold(6,nu)
          p(5,nu)=pold(7,nu)
          p(6,nu)=zip
          p(7,nu)=zip
        enddo
      endif

      if ((kcase==ktt_bbl) .or.
     &    (kcase==ktt_bbh) .or.
     &    (kcase==ktt_bbu)) then
c--- if we're doing ttb case, reduce # of momenta from 9 to 5
        do nu=1,4
          p(3,nu)=pold(3,nu)+pold(4,nu)+pold(5,nu)
          p(4,nu)=pold(6,nu)+pold(7,nu)+pold(8,nu)
          p(5,nu)=pold(9,nu)
          p(6,nu)=zip
          p(7,nu)=zip
        enddo
      endif

      if ((kcase==kW_twdk) .or. (kcase==kW_cwdk)) then
c--- if we're doing W+t, reduce # of momenta from 8 to 6
        do nu=1,4
          p(3,nu)=pold(3,nu)
          p(4,nu)=pold(4,nu)
          p(5,nu)=pold(5,nu)+pold(6,nu)+pold(7,nu)
          p(6,nu)=pold(8,nu)
          p(7,nu)=zip
          p(8,nu)=zip
        enddo
      endif

      if ((kcase==kZ_tdkj) .or. (kcase==kH_tdkj)) then
c--- if we're doing Z+t or H+t, reduce # of momenta from 9 to 7
        do nu=1,4
          p(3,nu)=pold(3,nu)
          p(4,nu)=pold(4,nu)
          p(5,nu)=pold(5,nu)+pold(6,nu)+pold(7,nu)
          p(6,nu)=pold(8,nu)
          p(7,nu)=pold(9,nu)
          p(8,nu)=zip
        enddo
      endif

      if ((kcase==k4ftwdk) .or. (kcase==kbq_tpq_jet)) then
c--- if we're doing (4F) t-channel single top with decay,
c--- reduce # of momenta from 8 to 6
        do nu=1,4
          p(3,nu)=pold(3,nu)+pold(4,nu)+pold(5,nu)
          p(4,nu)=pold(6,nu)
          p(5,nu)=pold(7,nu)
          p(6,nu)=pold(8,nu)
          p(7,nu)=zip
          p(8,nu)=zip
        enddo
      endif

      if (kcase==kqq_ttw) then
c--- reduce # of momenta from 11 to 7
        do nu=1,4
          p(3,nu)=pold(9,nu)
          p(4,nu)=pold(10,nu)
          p(5,nu)=pold(3,nu)+pold(4,nu)+pold(5,nu)
          p(6,nu)=pold(6,nu)+pold(7,nu)+pold(8,nu)
          p(7,nu)=pold(11,nu)
          p(8,nu)=zip
        enddo
      endif

c---Initialize the dipoles to zero
      sub=zip
      subv=zip
      msq = 0._dp
      msqv = 0._dp
      incldip(nd)=.true.

c--- default is all particles massless
      misq=zip
      mjsq=zip
      mksq=zip
      mijsq=zip

      pij=two*dot(p,ip,jp)
      pik=two*dot(p,ip,kp)
      pjk=two*dot(p,jp,kp)

      if ((ip <= 2) .and. (kp <= 2)) then
c**********************************************************************
c************************** INITIAL-INITIAL ***************************
c**********************************************************************
        omx=-(pij+pjk)/pik
        x=one-omx
        vtilde=pij/pik

c---Modification so that only close to singular subtracted
        if (-vtilde > aii) then
           incldip(nd)=.false.
           goto 99
        endif

        call transform_mass(p,ptrans,x,ip,jp,kp,misq,mjsq,mksq,mijsq)

        if ((kcase==kt_bbar) .or. (kcase==kbq_tpq)
     &  .or.(kcase==kW_twdk) .or. (kcase==kW_cwdk)
     &  .or.(kcase==kZ_tdkj) .or. (kcase==kH_tdkj)
     &  .or.(kcase==ktt_bbl) .or. (kcase==ktt_bbh)
     &  .or.(kcase==ktt_bbu) .or. (kcase==k4ftwdk)
     &  .or.(kcase==kqq_ttw) .or. (kcase==kbq_tpq_jet)) then
          if     ((kcase==kW_twdk) .or. (kcase==kW_cwdk)) then
            call extend_trans_wt(pold,p,ptrans,pext)
          elseif ((kcase==kZ_tdkj) .or. (kcase==kH_tdkj)) then
            call extend_trans_ztj(pold,p,ptrans,pext)
          elseif ((kcase==ktt_bbl)
     &       .or. (kcase==ktt_bbh)
     &       .or. (kcase==ktt_bbu)) then
            call extend_trans_ttb(pold,p,ptrans,pext)
          elseif ((kcase==k4ftwdk) .or. (kcase==kbq_tpq_jet)) then
            call extend_trans_stopb(pold,p,ptrans,pext)
          elseif ((kcase==kqq_ttw)) then
            call extend_trans_ttw(pold,p,ptrans,pext)
          elseif (kcase == kbq_tpq .and. currentContrib == 6) then
            call extend_trans_st(pold,p,ptrans,pext)
          else
            call extend_trans(pold,p,ptrans,pext)
          endif
          do j=1,mxpart
          do nu=1,4
            ptrans(j,nu)=pext(j,nu)
          enddo
          enddo
        endif

        call storeptilde(nd,ptrans)

        do nu=1,4
          vec(nu)=p(jp,nu)-pij/pik*p(kp,nu)
        enddo
        vecsq=-pij*pjk/pik

c--- if using a dynamic scale, set that scale with dipole kinematics
      if (nproc == 1610 .or. nproc == 1650 ) then
        call singletop2_set_dipscale(nd, ptrans, gsq)
      elseif (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif

        call subr_born(ptrans,msq)
        call subr_corr(ptrans,vec,ip,msqv)

        sub(qq)=-gsq/x/pij*(two/omx-one-x)
        sub(gq)=-gsq/pij
        sub(qg)=-gsq/x/pij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/pij*(x/omx+x*omx)
        subv   =+4._dp*gsq/x/pij*omx/x/vecsq

c**********************************************************************
c************************** INITIAL-FINAL *****************************
c**********************************************************************
      elseif ((ip <= 2) .and. (kp > 2)) then
        omx=-pjk/(pij+pik)
        x=one-omx
        u=pij/(pij+pik)
        omu=pik/(pij+pik)

c---determine mass of spectator
        mksq=max(p(kp,4)**2-p(kp,1)**2-p(kp,2)**2-p(kp,3)**2,zip)
        if (mksq>zip) then
          muksq=mksq/2._dp/(-dot(p,ip,jp)-dot(p,ip,kp))
          zp=omx/(omx+muksq)
        else
          zp=one
        endif

c---npart is the number of particles in the final state
c---transform the momenta so that only the first npart+1 are filled
        call transform_mass(p,ptrans,x,ip,jp,kp,misq,mjsq,mksq,mijsq)

        if ((kcase==kt_bbar) .or. (kcase==kbq_tpq)
     &  .or.(kcase==kW_twdk) .or. (kcase==kW_cwdk)
     &  .or.(kcase==kZ_tdkj) .or. (kcase==kH_tdkj)
     &  .or.(kcase==ktt_bbl) .or. (kcase==ktt_bbh)
     &  .or.(kcase==ktt_bbu) .or. (kcase==k4ftwdk)
     &  .or.(kcase==kqq_ttw) .or. (kcase==kbq_tpq_jet)) then
          if     ((kcase==kW_twdk) .or. (kcase==kW_cwdk)) then
            call extend_trans_wt(pold,p,ptrans,pext)
          elseif ((kcase==kZ_tdkj) .or. (kcase==kH_tdkj)) then
            call extend_trans_ztj(pold,p,ptrans,pext)
          elseif ((kcase==ktt_bbl)
     &       .or. (kcase==ktt_bbh)
     &       .or. (kcase==ktt_bbu)) then
            call extend_trans_ttb(pold,p,ptrans,pext)
          elseif ((kcase==k4ftwdk) .or. (kcase==kbq_tpq_jet)) then
            call extend_trans_stopb(pold,p,ptrans,pext)
          elseif ((kcase==kqq_ttw)) then
            call extend_trans_ttw(pold,p,ptrans,pext)
          elseif (kcase == kbq_tpq .and. currentContrib == 6) then
            call extend_trans_st(pold,p,ptrans,pext)
          else
            call extend_trans(pold,p,ptrans,pext)
          endif
          do j=1,mxpart
          do nu=1,4
            ptrans(j,nu)=pext(j,nu)
          enddo
          enddo
        endif

        call storeptilde(nd,ptrans)

c--- if using a dynamic scale, set that scale with dipole kinematics
      if (nproc == 1610 .or. nproc == 1650) then
        call singletop2_set_dipscale(nd, ptrans, gsq)
      elseif (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif

c--- Calculate the matrix element now because it might be needed
c--- in the final-initial segment, regardless of whether or not the
c--- alfa cut fails here
        call subr_born(ptrans,msq)
c---Modification so that only close to singular subtracted
c---Do not set incldip because initial-final can fail
c---but final initial still needs to be tested
c--- [note: for massive dipoles initial-final =/= final-initial,
c---        but for the case 4ftwdk we call this routine with all
c---        masses = 0, so the special case is handled gracefully]

c---Modification so that only close to singular subtracted
        if (u > aif) goto 99

        do nu=1,4
           vec(nu)=(p(jp,nu)/u-p(kp,nu)/omu)/sqrt(pjk)
        enddo

        call subr_corr(ptrans,vec,ip,msqv)
        sub(qq)=-gsq/x/pij*(two/(omx+u)-one-x)
        sub(gq)=-gsq/pij
        sub(qg)=-gsq/x/pij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/pij*(one/(omx+u)-one+x*omx)
        subv   =-4._dp*gsq/x/pij*(omx/x*u*(one-u))

c**********************************************************************
c************************** FINAL-INITIAL *****************************
c**********************************************************************
      elseif ((ip > 2) .and. (kp <= 2)) then
        do jproc=1,4
        if ((jproc==qq) .and. (qqproc .eqv. .false.)) goto 79
        if ((jproc==gq) .and. (gqproc .eqv. .false.)) goto 79
        if ((jproc==qg) .and. (qgproc .eqv. .false.)) goto 79
        if ((jproc==gg) .and. (ggproc .eqv. .false.)) goto 79

        if (jproc==qq) then
        mijsq=mqsq
c--- the masses of i and j have been switched
        misq=mqsq
        mjsq=zip
        elseif (jproc==qg) then
        go to 79
        elseif (jproc==gq) then
        mijsq=zip
        misq=mqsq
        mjsq=mqsq
        elseif (jproc==gg) then
        goto 79
        endif
        omx=(mijsq-misq-mjsq-pij)/(pjk+pik)
        x=one-omx

        do nu=1,4
          qij(nu)=p(ip,nu)+p(jp,nu)
          q(nu)=qij(nu)+p(kp,nu)
        enddo
        qsq=q(4)**2-q(1)**2-q(2)**2-q(3)**2
        qijsq=qij(4)**2-qij(1)**2-qij(2)**2-qij(3)**2
        call transform_mass(p,ptrans,x,ip,jp,kp,misq,mjsq,mksq,mijsq)

        if ((kcase==kt_bbar) .or. (kcase==kbq_tpq)
     &  .or.(kcase==kW_twdk) .or. (kcase==kW_cwdk)
     &  .or.(kcase==kZ_tdkj) .or. (kcase==kH_tdkj)
     &  .or.(kcase==ktt_bbl) .or. (kcase==ktt_bbh)
     &  .or.(kcase==ktt_bbu) .or. (kcase==k4ftwdk)
     &  .or.(kcase==kqq_ttw) .or. (kcase==kbq_tpq_jet)) then
          if     ((kcase==kW_twdk) .or. (kcase==kW_cwdk)) then
            call extend_trans_wt(pold,p,ptrans,pext)
          elseif ((kcase==kZ_tdkj) .or. (kcase==kH_tdkj)) then
            call extend_trans_ztj(pold,p,ptrans,pext)
          elseif ((kcase==ktt_bbl)
     &       .or. (kcase==ktt_bbh)
     &       .or. (kcase==ktt_bbu)) then
            call extend_trans_ttb(pold,p,ptrans,pext)
          elseif ((kcase==k4ftwdk) .or. (kcase==kbq_tpq_jet)) then
            call extend_trans_stopb(pold,p,ptrans,pext)
          elseif ((kcase==kqq_ttw)) then
            call extend_trans_ttw(pold,p,ptrans,pext)
          elseif (kcase == kbq_tpq .and. currentContrib == 6) then
            call extend_trans_st(pold,p,ptrans,pext)
          else
            call extend_trans(pold,p,ptrans,pext)
          endif
          do j=1,mxpart
          do nu=1,4
            ptrans(j,nu)=pext(j,nu)
          enddo
          enddo
        endif

        call storeptilde(nd,ptrans)
        z=pik/(pik+pjk)
        omz=pjk/(pik+pjk)
c--- note that musq is related to msq by musq = msq/(2pij_tilde.pa)
c--- and 2pij_tilde.pa = (Qsq-mijsq)/x
        zp=omx*(Qsq-mijsq)/x+mijsq+misq-mjsq
c---    quadratic
        root=omx*(Qsq-mijsq)/x+mijsq-misq-mjsq
        root=sqrt(root**2-4._dp*misq*mjsq)
        zm=(zp-root)/(2._dp*(omx*(Qsq-mijsq)/x+mijsq))
        zp=(zp+root)/(2._dp*(omx*(Qsq-mijsq)/x+mijsq))

c--- if using a dynamic scale, set that scale with dipole kinematics
      if (nproc == 1610 .or. nproc == 1650) then
        call singletop2_set_dipscale(nd, ptrans, gsq)
      elseif (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif

        call subr_born(ptrans,msq)

c---Modification so that only close to singular subtracted
        if (omx > afi) goto 99

        do nu=1,4
          vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
        enddo

        call subr_corr(ptrans,vec,ip,msqv)

        if (jproc == qq) then
        sub(qq)=+gsq/x/(qijsq-mijsq)*(two/(omz+omx)-one-z-2._dp*mqsq/pij)
        elseif (jproc == gq) then
        sub(gq)=+gsq/x/(qijsq-mijsq)
        subv   =+4._dp*gsq/x/qijsq/(qijsq-mijsq)
        endif
 79     continue
        enddo

c**********************************************************************
c*************************** FINAL-FINAL ******************************
c**********************************************************************
      elseif ((ip > 2) .and. (kp > 2)) then
c------Eq-(5.2)

c----Form momentum vectors
      y=pij/(pij+pjk+pik)
      z=pik/(pjk+pik)
      omz=one-z
      omy=one-y
       do nu=1,4
c  create 4-momentum of sum emitter+emittee
          qij(nu)=p(ip,nu)+p(jp,nu)
c  create 4-momentum of sum emitter+emittee+spectator
          q(nu)=qij(nu)+p(kp,nu)
       enddo
c  and square them
      Qsq=q(4)**2-q(1)**2-q(2)**2-q(3)**2
      qijsq=qij(4)**2-qij(1)**2-qij(2)**2-qij(3)**2


c---loop over the different possibilities
      do jproc=1,4
      if ((jproc==qq) .and. (qqproc .eqv. .false.)) goto 80
      if ((jproc==gq) .and. (gqproc .eqv. .false.)) goto 80
      if ((jproc==qg) .and. (qgproc .eqv. .false.)) goto 80
      if ((jproc==gg) .and. (ggproc .eqv. .false.)) goto 80

      if (jproc==qq) then
c q->qg
      mijsq=mqsq
c--- the masses of i and j have been switched
      misq=mqsq
      mjsq=zip
      elseif (jproc==qg) then
c q->gq
      go to 80
      elseif (jproc==gq) then
c g->qqbar
      mijsq=zip
      misq=mqsq
      mjsq=mqsq
      elseif (jproc==gg) then
c g->gg
      mijsq=zip
      misq=zip
      mjsq=zip
      endif

      muisq=misq/Qsq
      mujsq=mjsq/Qsq
      muijsq=mijsq/Qsq

c---determine mass of spectator
      mksq=max(p(kp,4)**2-p(kp,1)**2-p(kp,2)**2-p(kp,3)**2,zip)
      if (mksq>zip) then
        muksq=mksq/Qsq
      muk=sqrt(muksq)
        yp=one-2._dp*muk*(one-muk)/(one-muisq-mujsq-muksq)
      else
        muksq=zip
        muk=zip
        yp=one
      endif

      if (y > yp) then
         write(6,*) 'Problems with phase space in dips_mass.f'
         write(6,*) "y, yp = ", y, yp
         incldip(nd)=.false.
         go to 99
         !stop
      endif

c---Modification so that only close to singular subtracted
       if (y > aff*yp) then
         incldip(nd)=.false.
         go to 99
       endif


c      viji=sqrt((one-muijsq-muisq)**2-4._dp*mijsq*muisq)
c     & /(one-muijsq-muisq)
c      write(6,*) 'viji',viji
c      vijk=sqrt((one-qijsq/Qsq-muksq)**2-4._dp*qijsq/Qsq*muksq)
c     & /(one-qijsq/Qsq-muksq)
c      write(6,*) vijk

c--- Note that identities of i and j have been exchanged
      viji=sqrt(((one-mujsq-muisq-muksq)*y)**2-4._dp*muisq*mujsq)
     & /((one-mujsq-muisq-muksq)*y+2._dp*mujsq)

      vijk=sqrt((2._dp*muksq+(one-mujsq-muisq-muksq)*omy)**2-4._dp*muksq)
     & /((one-mujsq-muisq-muksq)*omy)

c      ym=2._dp*mui*muj/(one-muisq-mujsq-muksq)

      zp=(2._dp*mujsq+(one-muisq-mujsq-muksq)*y)
     & /(2._dp*(muisq+mujsq+(one-muisq-mujsq-muksq)*y))
      zm=zp*(one-viji*vijk)
      zp=zp*(one+viji*vijk)
c---calculate the ptrans-momenta
       call transform_mass(p,ptrans,y,ip,jp,kp,misq,mjsq,mksq,mijsq)

        if ((kcase==kt_bbar) .or. (kcase==kbq_tpq)
     &  .or.(kcase==kW_twdk) .or. (kcase==kW_cwdk)
     &  .or.(kcase==kZ_tdkj) .or. (kcase==kH_tdkj)
     &  .or.(kcase==ktt_bbl) .or. (kcase==ktt_bbh)
     &  .or.(kcase==ktt_bbu) .or. (kcase==k4ftwdk)
     &  .or.(kcase==kqq_ttw) .or. (kcase==kbq_tpq_jet)) then
          if     ((kcase==kW_twdk) .or. (kcase==kW_cwdk)) then
            call extend_trans_wt(pold,p,ptrans,pext)
          elseif ((kcase==kZ_tdkj) .or. (kcase==kH_tdkj)) then
            call extend_trans_ztj(pold,p,ptrans,pext)
          elseif ((kcase==ktt_bbl)
     &       .or. (kcase==ktt_bbh)
     &       .or. (kcase==ktt_bbu)) then
            call extend_trans_ttb(pold,p,ptrans,pext)
          elseif ((kcase==k4ftwdk) .or. (kcase==kbq_tpq_jet)) then
            call extend_trans_stopb(pold,p,ptrans,pext)
          elseif ((kcase==kqq_ttw)) then
            call extend_trans_ttw(pold,p,ptrans,pext)
          elseif (kcase == kbq_tpq .and. currentContrib == 6) then
            call extend_trans_st(pold,p,ptrans,pext)
          else
            call extend_trans(pold,p,ptrans,pext)
          endif
          do j=1,mxpart
          do nu=1,4
            ptrans(j,nu)=pext(j,nu)
          enddo
          enddo
        endif

c       write(6,*) 'Dipole ',nd, 'ptrans'
c       call writeout(ptrans)

       call storeptilde(nd,ptrans)

       ztmi=z-0.5_dp+0.5_dp*vijk
       ztmj=omz-0.5_dp+0.5_dp*vijk

c--- if using a dynamic scale, set that scale with dipole kinematics
        if (nproc == 1610 .or. nproc == 1650) then
          call singletop2_set_dipscale(nd, ptrans, gsq)
        elseif (dynamicscale) then
          call scaleset(initscale,initfacscale,ptrans)
          dipscale(nd)=facscale
        endif

       call subr_born(ptrans,msq)

       do nu=1,4
         vec(nu)=ztmi*p(ip,nu)-ztmj*p(jp,nu)
       enddo

       if (kcase==kqq_tbg) then
         ipt=5
       elseif (kcase==kbq_tpq_jet) then
         ipt=7
       else
         if (ip < kp) then
           ipt=5
         else
           ipt=6
         endif
       endif
       call subr_corr(ptrans,vec,ipt,msqv)

      if     (jproc == qq) then
        vtijk=sqrt((one-muijsq-muksq)**2-4._dp*muijsq*muksq)
     &  /(one-muijsq-muksq)

        sub(qq)=gsq/(qijsq-mijsq)*(two/(one-z*omy)
     &  -vtijk/vijk*(one+z+2._dp*mqsq/pij))

      elseif (jproc == gq) then
        sub(gq)=gsq/(qijsq-mijsq)/vijk*(
     &          one-two*kappa*(zp*zm-mqsq/qijsq))
        subv   =+4._dp*gsq/(qijsq-mijsq)/qijsq/vijk
        subv_gq=subv ! put in common block

      elseif (jproc == gg) then
        sub(gg)=two*gsq/(qijsq-mijsq)*(one/(one-z*omy)+one/(one-omz*omy)
     &    -(two-kappa*zp*zm)/vijk)
        subv   =+4._dp*gsq/(qijsq-mijsq)/pij/vijk
        subv_gg=subv ! put in common block
      endif

 80   continue
      enddo
      endif

c--- fall through to here, so that p retains the value it entered with
   99 continue

      do j=1,mxpart
      do nu=1,4
        p(j,nu)=pold(j,nu)
      enddo
      enddo

      return
      end

