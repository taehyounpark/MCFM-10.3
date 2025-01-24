!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wgam_v(p,msqv)
      implicit none
      include 'types.f'

c----Author R.K.Ellis August 2002
c====Virtual corrections to
c     q(-p1)+qbar(-p2)-->e^-(p3)+nu(p4)+gamma(p5)
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'scheme.f'
      include 'qcdcouple.f'
      include 'ckm.f'
      include 'zprods_decl.f'
      include 'nwz.f'
      include 'blha.f'
      include 'zcouple_cms.f'
      integer:: j,k
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),qbq,qqb
      real(dp):: fac
      complex(dp):: agamtree,agamvirt

      call spinoru(5,p,za,zb)

      fac=ason2pi*cf*aveqq*2._dp*xn*abs((zesq/zxw)**2*zesq)

      scheme='dred'
      if (nwz == -1) then
c ie ub-d

      if (useblha == 0) then
      qbq=+fac*2._dp*real(
     &conjg(agamtree(1,2,3,4,5,za,zb,+1))*agamvirt(1,2,3,4,5,za,zb,+1))
     &    +fac*2._dp*real(
     &conjg(agamtree(1,2,3,4,5,za,zb,-1))*agamvirt(1,2,3,4,5,za,zb,-1))
      endif
c ie d-ub
      qqb=+fac*2._dp*real(
     &conjg(agamtree(2,1,3,4,5,za,zb,+1))*agamvirt(2,1,3,4,5,za,zb,+1))
     &    +fac*2._dp*real(
     &conjg(agamtree(2,1,3,4,5,za,zb,-1))*agamvirt(2,1,3,4,5,za,zb,-1))

      elseif (nwz == +1) then
c ie db-u

      if (useblha == 0) then
      qbq=+fac*2._dp*real(
     &conjg(agamtree(2,1,4,3,5,zb,za,+1))*agamvirt(2,1,4,3,5,zb,za,+1))
     &    +fac*2._dp*real(
     &conjg(agamtree(2,1,4,3,5,zb,za,-1))*agamvirt(2,1,4,3,5,zb,za,-1))
      endif
c ie u-db

      qqb=+fac*2._dp*real(
     &conjg(agamtree(1,2,4,3,5,zb,za,+1))*agamvirt(1,2,4,3,5,zb,za,+1))
     &    +fac*2._dp*real(
     &conjg(agamtree(1,2,4,3,5,zb,za,-1))*agamvirt(1,2,4,3,5,zb,za,-1))

      endif

      do j=-nf,nf
      do k=-nf,nf
c--- set msqv=0 to initalize
      msqv(j,k)=0._dp
          if ((j > 0) .and. (k < 0)) then
            msqv(j,k)=Vsq(j,k)*qqb
          elseif ((j < 0) .and. (k > 0)) then
            msqv(j,k)=Vsq(j,k)*qbq
          endif
      enddo
      enddo

      return
      end




      function agamvirt(p1,p2,p3,p4,p5,za,zb,hgamma)
      implicit none
      include 'types.f'
      complex(dp):: agamvirt

c----Author R.K.Ellis August 2002
c     q(-p1)+qbar(-p2)-->e^-(p3)+nu(p4)+gamma(p5)
      include 'constants.f'
      include 'mxpart.f'
      include 'zprods_decl.f'
      include 'masses.f'
      complex(dp):: agamtree,vpole,fagamma,fbgamma,prop,vpl
      real(dp):: s34,s12
      integer:: p1,p2,p3,p4,p5,hgamma


      s34=real(za(p3,p4)*zb(p4,p3))
      s12=real(za(p1,p2)*zb(p2,p1))
      prop=s34/cmplx(s34-wmass**2,wmass*wwidth,kind=dp)
      vpl=vpole(s12)
      if    (hgamma == +1) then
        agamvirt=prop*(
     &  +Qd*fagamma(p1,p2,p3,p4,p5,za,zb)
     &  +Qu*fbgamma(p1,p2,p3,p4,p5,za,zb))
     &          +vpl*agamtree(p1,p2,p3,p4,p5,za,zb,+1)
      elseif (hgamma == -1) then
        agamvirt=prop*(
     &  +Qu*fagamma(p2,p1,p4,p3,p5,zb,za)
     &  +Qd*fbgamma(p2,p1,p4,p3,p5,zb,za))
     &          +vpl*agamtree(p1,p2,p3,p4,p5,za,zb,-1)
      endif

      return
      end
