!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c***********************************************************************
c     Author: J. M. Campbell                                           *
c     August, 2001                                                     *
c                                                                      *
c     Replica of dipolesub.f, except for the fact that extra matrix    *
c     element arrays are called in the Born term                       *
c                                                                      *
c     Calculates the nj-jet subtraction term corresponding to dipole   *
c     nd with momentum p and dipole kinematics (ip,jp) wrt kp          *
c     Automatically chooses dipole kind                                *
c     Returns the dipoles in sub,subv and matrix elements in msq,msqv  *
c     nd labels the dipole configurations                              *
c     ip labels the emitter parton                                     *
c     jp labels the emitted parton                                     *
c     kp labels the spectator parton                                   *
c     msq - lowest order matrix elements at rescaled momentum, msq(j,k)*
c     msqv -  lowest order matrix elements at rescaled momentum        *
c      with emitter contracted with appropriate vector, msqv(j,k)      *
c     subr_born is the subroutine which call the born process          *
c     subr_corr is the subroutine which call the born process dotted   *
c      with vec for an emitted gluon only                              *
c     mqq - 4-quark contribution to lowest order matrix elements sqd.  *
c     msqx - lowest order matrix elements with 4 indices, msqx(j,k,l,m)*
c            Sum_{l,m} msqx(j,k,l,m) = msq(j,k)                        *
c     mg - 2-quark contribution to lowest order matrix elements sqd,   *
c           separated by colours                                       *
c     mvg - 2-quark contribution to lowest order matrix elements sqd,  *
c           separated by colours, contracted with appropriate vector   *
c     mvxg - lowest order matrix elements with 4 indices and           *
c        contracted with appropriate vector, msqvx(j,k,l,m)            *
c        Sum_{l,m} msqvx(j,k,l,m) = msqv(j,k)                          *
c***********************************************************************

      subroutine dipsx(nd,p,ip,jp,kp,sub,subv,msq,msqv,
     & subr_born,subr_corr,mqq,msqx,mg,mvg,mvxg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      include 'alfacut.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      include 'incldip.f'
      real(dp):: p(mxpart,4),ptrans(mxpart,4),sub(4),subv,vecsq
      real(dp):: x,omx,z,omz,y,omy,u,omu,sij,sik,sjk,dot,vec(4)
      real(dp):: msq(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf),vtilde
      real(dp):: mqq(0:2,-nf:nf,-nf:nf)
      real(dp):: msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: mg(0:2,-nf:nf,-nf:nf)
      real(dp):: mvg(0:2,-nf:nf,-nf:nf)
      real(dp):: mvxg(-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      integer:: nd,ip,jp,kp,nu,j,k
c--      logical:: includedipole
      external subr_born,subr_corr

c---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0._dp
      enddo
      subv=0._dp
      msq = 0._dp
      msqv = 0._dp
      mqq=0._dp
      msqx=0._dp
      mg=0._dp
      mvg=0._dp
      mvxg=0._dp
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
        call storeptilde(nd,ptrans)

        vecsq=-sij*sjk/sik
        do nu=1,4
          vec(nu)=p(jp,nu)-vtilde*p(kp,nu)
        enddo

c--- if using a dynamic scale, set that scale with dipole kinematics
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif

        call subr_born(ptrans,msq,mqq,msqx,mg)
        call subr_corr(ptrans,vec,ip,msqv,mvg,mvxg)

        sub(qq)=-gsq/x/sij*(two/omx-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/sij*(x/omx+x*omx)
        subv   =+4._dp*gsq/x/sij*omx/x/vecsq

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
        call storeptilde(nd,ptrans)

c--- if using a dynamic scale, set that scale with dipole kinematics
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif

c--- Calculate the matrix element now because it might be needed
c--- in the final-initial segment, regardless of whether or not the
c--- alfa cut fails here
        call subr_born(ptrans,msq,mqq,msqx,mg)
c---Modification so that only close to singular subtracted
c---Do not set incldip because initial-final can fail
c---but final initial needs still to be tested
        if (u > aif) return

        do nu=1,4
           vec(nu)=p(jp,nu)/u-p(kp,nu)/omu
        enddo

        call subr_corr(ptrans,vec,ip,msqv,mvg,mvxg)
        sub(qq)=-gsq/x/sij*(two/(omx+u)-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/sij*(one/(omx+u)-one+x*omx)
        subv   =-4._dp*gsq/x/sij*(omx/x*u*(one-u)/sjk)

c**********************************************************************
c************************** FINAL-INITIAL *****************************
c**********************************************************************
      elseif ((ip > 2) .and. (kp <= 2)) then
c--- note, here we assume that msq kinematics are already taken care of
c--- for msq, although msqv must be recalculated each time
        omx=-sij/(sjk+sik)
c---Modification so that only close to singular subtracted
        if (omx > afi) return

        x=one-omx
        z=sik/(sik+sjk)
        omz=sjk/(sik+sjk)
        do nu=1,4
          vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
        enddo
c---call msqv again because vec has changed
        do j=1,mxpart
        do k=1,4
          ptrans(j,k)=ptilde(nd,j,k)
        enddo
        enddo

c--- if using a dynamic scale, set that scale with dipole kinematics
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif

c--- do something special if we're doing W+2,Z+2jet (jp  /=  7)
        if (jp  /= 7) then
          if (ip < 7) then
c ie for cases 56_i,65_i
          call subr_corr(ptrans,vec,5,msqv,mvg,mvxg)
          else
c ie for cases 76_i,75_i
          call subr_corr(ptrans,vec,6,msqv,mvg,mvxg)
          endif
        else
c ie for cases 57_i,67_i
          call subr_corr(ptrans,vec,ip,msqv,mvg,mvxg)
        endif

        sub(qq)=+gsq/x/sij*(two/(omz+omx)-one-z)
        sub(gq)=+gsq/x/sij
        sub(gg)=+2._dp*gsq/x/sij*(one/(omz+omx)+one/(z+omx)-two)
        subv   =+4._dp*gsq/x/sij/sij


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
       call storeptilde(nd,ptrans)

       do nu=1,4
         vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
       enddo

c--- if using a dynamic scale, set that scale with dipole kinematics
      if (dynamicscale) then
        call scaleset(initscale,initfacscale,ptrans)
        dipscale(nd)=facscale
      endif

       call subr_born(ptrans,msq,mqq,msqx,mg)
       if (ip < kp) then
         call subr_corr(ptrans,vec,5,msqv,mvg,mvxg)
       else
         call subr_corr(ptrans,vec,6,msqv,mvg,mvxg)
       endif

       sub(qq)=gsq/sij*(two/(one-z*omy)-one-z)
       sub(gq)=gsq/sij
       sub(gg)=gsq/sij*(two/(one-z*omy)+two/(one-omz*omy)-four)
       subv   =+4._dp*gsq/sij/sij

      endif

      return
      end

