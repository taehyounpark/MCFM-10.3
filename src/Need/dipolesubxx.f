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
c     subr_born is the subroutine which call the born process          *
c     subr_corr is the subroutine which call the born process dotted   *
c      with vec for an emitted gluon only                              *
c     msqx - lowest order matrix elements with 4 indices, msqx(j,k,l,m)*
c            Sum_{l,m} msqx(j,k,l,m) = msq(j,k)                        *
c     mvxg - lowest order matrix elements with 4 indices and           *
c        contracted with appropriate vector, msqvx(j,k,l,m)            *
c        Sum_{l,m} msqvx(j,k,l,m) = msqv(j,k)                          *
c***********************************************************************

      subroutine dipsxx(nd,p,ip,jp,kp,sub,subv,
     & subr_born,subr_corr,msqx,msqvx)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'qqgg.f'
      include 'ptilde.f'
      real(dp):: p(mxpart,4),ptrans(mxpart,4),sub(4),subv,vecsq
      real(dp):: x,omx,z,omz,y,omy,u,omu,sij,sik,sjk,dot,vec(4)
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: msqx(0:2,-nf:nf,-nf:nf,-nf:nf,-nf:nf)
      real(dp):: msqvx(0:2,-1:1,-1:1,-1:1,-1:1)
      integer:: nd,ip,jp,kp,nu,j,k
      external subr_born,subr_corr

c---Initialize the dipoles to zero
      do j=1,4
      sub(j)=0._dp
      enddo

      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp)
      sjk=two*dot(p,jp,kp)

      if ((ip <= 2) .and. (kp <= 2)) then
c**********************************************************************
c************************** INITIAL-INITIAL ***************************
c**********************************************************************
        omx=-(sij+sjk)/sik
        x=one-omx

        call transform(p,ptrans,x,ip,jp,kp)
        call storeptilde(nd,ptrans)
        do nu=1,4
          vec(nu)=p(jp,nu)-sij/sik*p(kp,nu)
        enddo
        vecsq=-sij*sjk/sik
        call subr_born(ptrans,msq,msqx)
        call subr_corr(ptrans,vec,ip,msqvx)

        sub(qq)=-gsq/x/sij*(two/omx-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/sij*(x/omx+x*omx)
        subv   =+4._dp*gsq/x/sij*omx/x/vecsq

c**********************************************************************
c************************** INITIAL-FINAL *****************************
c**********************************************************************
      elseif ((ip <= 2) .and. (kp > 2)) then

        omx=-sjk/(sij+sik)
        x=one-omx
        u=sij/(sij+sik)
        omu=sik/(sij+sik)
c---npart is the number of particles in the final state
c---transform the momenta so that only the first npart+1 are filled
        call transform(p,ptrans,x,ip,jp,kp)
        call storeptilde(nd,ptrans)
        do nu=1,4
           vec(nu)=p(jp,nu)/u-p(kp,nu)/omu
        enddo
        call subr_born(ptrans,msq,msqx)
        call subr_corr(ptrans,vec,ip,msqvx)
        sub(qq)=-gsq/x/sij*(two/(omx+u)-one-x)
        sub(gq)=-gsq/sij
        sub(qg)=-gsq/x/sij*(one-two*x*omx)
        sub(gg)=-2._dp*gsq/x/sij*(one/(omx+u)-one+x*omx)
        subv   =-4._dp*gsq/x/sij*(omx/x*u*(one-u)/sjk)

      elseif ((ip > 2) .and. (kp <= 2)) then
c**********************************************************************
c************************** FINAL-INITIAL *****************************
c**********************************************************************
c--- note, here we assume that msq kinematics are already taken care of
c--- for msq, although msqv must be recalculated each time
        omx=-sij/(sjk+sik)
        x=one-omx
        z=sik/(sik+sjk)
        omz=sjk/(sik+sjk)
        do nu=1,4
          vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
        enddo
c---call again because vec has changed
        do j=1,mxpart
        do k=1,4
          ptrans(j,k)=ptilde(nd,j,k)
        enddo
        enddo
c--- do something special if (jp  /=  5)
        if (jp  /=  5) then
          if (ip < 5) then
c ie for cases 34_i,43_i
          call subr_corr(ptrans,vec,3,msqvx)
          else
c ie for cases 54_i,53_i
          call subr_corr(ptrans,vec,4,msqvx)
          endif
        else
c ie for cases 35_i,45_i
          call subr_corr(ptrans,vec,ip,msqvx)
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
       z=sik/(sjk+sik)
       omz=one-z
       omy=one-y
c---calculate the ptrans-momenta

       call transform(p,ptrans,y,ip,jp,kp)
       call storeptilde(nd,ptrans)
       do nu=1,4
         vec(nu)=z*p(ip,nu)-omz*p(jp,nu)
       enddo
       call subr_born(ptrans,msq,msqx)
       if (ip < kp) then
         call subr_corr(ptrans,vec,3,msqvx)
       else
         call subr_corr(ptrans,vec,4,msqvx)
       endif

       sub(qq)=gsq/sij*(two/(one-z*omy)-one-z)
       sub(gq)=gsq/sij
       sub(gg)=gsq/sij*(two/(one-z*omy)+two/(one-omz*omy)-four)
       subv   =+4._dp*gsq/sij/sij

      endif

      return
      end

