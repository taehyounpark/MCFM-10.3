!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c***********************************************************************
c     This subroutine calculates dipoles with an identified            *
c     Final state photon                                               *
c     C. Williams Dec 2010                                             *
c     Returns the dipoles in sub,subv and matrix elements in msq,msqv  *
c     nd labels the dipole configurations                              *
c     ip labels the emitter PHOTON                                     *
c     jp labels the emitted parton                                     *
c     kp labels the spectator parton                                   *
c     subr_born is the subroutine which call the born process          *
c     subr_corr is the subroutine which call the born process dotted   *
c     with vec for an emitted gluon only                               *
c***********************************************************************

      subroutine dipsfrag(nd,p,ip,jp,kp,sub,msq,subr_born)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'ewcouple.f'
      include 'ptilde.f'
      include 'dynamicscale.f'
      include 'initialscales.f'
      include 'dipolescale.f'
      include 'facscale.f'
      include 'betacut.f'
      include 'lastphot.f'
      include 'incldip.f'
      real(dp):: p(mxpart,4),ptrans(mxpart,4),sub
      real(dp):: z,omz,sij,sik,sjk,dot,u,p_phys(mxpart,4)
      real(dp):: msq(-nf:nf,-nf:nf),tmp
      integer:: nd,ip,jp,kp,j,k,ipt
      logical:: check_nv
      external subr_born
      z=0._dp
      omz=1._dp
      u=0._dp
      sub=0._dp
      do j=-nf,nf
         do k=-nf,nf
            msq(j,k)=0._dp
         enddo
      enddo

      incldip(nd)=.true.

      sij=two*dot(p,ip,jp)
      sik=two*dot(p,ip,kp)
      sjk=two*dot(p,jp,kp)


c******************************************************************************
c*********************** INITIAL - INITIAL ************************************
c******************************************************************************

c*** I === I not implemented (need photon PDFS rather than frags) need rapidity cuts to remove collinear sing
      if ((ip <= 2) .and. (kp <= 2)) return

c******************************************************************************
c*********************** INITIAL - FINAL **************************************
c******************************************************************************

c*** I === F not implemented (need photon PDFS rather than frags) need rapidity cuts to remove collinear sing
      if ((ip <= 2) .and. (kp > 2)) return

c******************************************************************************
c*********************** FINAL - INITIAL  *************************************
c******************************************************************************


      if ((ip > 2) .and. (kp <= 2)) then


         if(check_nv(p,ip,jp,kp).eqv..false.) then
            incldip(nd)=.false.
            return
         endif
         call transformfrag(p,ptrans,z,ip,jp,kp)
         ipt=ip
         if (ip < lastphot) then
             do j=1,4
             tmp=ptrans(ip,j)
             ptrans(ip,j)=ptrans(lastphot,j)
             ptrans(lastphot,j)=tmp
             enddo
             ipt=lastphot
         endif

        if (dynamicscale) then
c--- rescale momentum (c.f. code below) in order to obtain physical
c--- momentum for setting dynamic scale
           p_phys(:,:)=ptrans(:,:)
           p_phys(lastphot,:)=z*ptrans(lastphot,:)
           call scaleset(initscale,initfacscale,p_phys)
           dipscale(nd)=facscale
        endif

c-- here, ptrans contains the right momentum for evaluating
c-- the LO matrix element with the last photon replaced by a jet
         call subr_born(ptrans,msq)
         omz=one-z
         sub=two*(esq/sij)*((one+omz**2)/z)

c--- now rescale momentum of last photon entry to represent
c--- the observed photon momentum
c--- NB: momentum will no longer be conserved in ptrans
         ptrans(lastphot,:)=z*ptrans(lastphot,:)

         call storeptilde(nd,ptrans)
c--- store z for use in isolation routine
         call store_zdip(nd,z)

c----  Check that photon will still pass cuts
c         if (phot_pass(ptrans,ipt,z) .eqv. .false.) return
c--- if using a dynamic scale, set that scale with dipole kinematics
c        if (dynamicscale) then
c           call rescale_z_dip(ptrans,nd,ipt)
c           call scaleset(initscale,initfacscale,ptrans)
c           call return_z_dip(ptrans,nd,ipt)
c           dipscale(nd)=facscale
c        endif



c******************************************************************************
c*********************** FINAL - FINAL ****************************************
c******************************************************************************

      elseif ((ip > 2) .and. (kp > 2)) then

         write(6,*) 'Final-final fragmentation dipole not implemented.'
         stop

         z=(sij+sik)/(sik+sjk+sij)
         omz=one-z

         u=sij/(sij+sik)

         if(u>bff) then
            incldip(nd)=.false.
            return
         endif

c--- Calculate the ptrans-momenta
         call transformfrag(p,ptrans,z,ip,jp,kp)
         ipt=ip
         if (ip < lastphot) then
             do j=1,4
             tmp=ptrans(ip,j)
             ptrans(ip,j)=ptrans(lastphot,j)
             ptrans(lastphot,j)=tmp
             enddo
             ipt=lastphot
         endif

        if (dynamicscale) then
c--- rescale momentum (c.f. code below) in order to obtain physical
c--- momentum for setting dynamic scale
           p_phys(:,:)=ptrans(:,:)
           p_phys(lastphot,:)=z*ptrans(lastphot,:)
           call scaleset(initscale,initfacscale,p_phys)
           dipscale(nd)=facscale
        endif


         call subr_born(ptrans,msq)
         omz=one-z
         sub=two*(esq/sij)*((one+omz**2)/z)

c--- now rescale momentum of last photon entry to represent
c--- the observed photon momentum
c--- NB: momentum will no longer be conserved in ptrans
         ptrans(lastphot,:)=z*ptrans(lastphot,:)

         call storeptilde(nd,ptrans)
c--- store z for use in isolation routine
         call store_zdip(nd,z)

      endif

      return
      end

