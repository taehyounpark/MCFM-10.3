!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_dm_monojet_v(p,msq)
      implicit none
      include 'types.f'
c----- routine for calcuating the virtual corrections to monojet production

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'dm_params.f'
      include 'scheme.f'
      include 'qcdcouple.f'
      include 'epinv.f'
      include 'nflav.f'
      real(dp)::  qqbg(2),qbqg(2)
      real(dp):: qbgq(2),qgqb(2)
      real(dp):: gqqb(2),gqbq(2)
      real(dp)::  qqbg_nf(2),qbqg_nf(2)
      real(dp):: qbgq_nf(2),qgqb_nf(2)
      real(dp):: gqqb_nf(2),gqbq_nf(2)
      real(dp):: qqbg_sum(nf),qbqg_sum(nf)
      real(dp):: qgqb_sum(nf),qbgq_sum(nf)
      real(dp):: gqqb_sum(nf),gqbq_sum(nf)
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      integer:: j,k
      real(dp):: msq0(-nf:nf,-nf:nf),fac,subuv
c      logical:: check_QED
c      common/check_QED/check_QED
      real(dp):: fac_dm,s34
      complex(dp):: cprop
      real(dp):: propsq

c      check_QED=.false.

c      if(check_QED) then
c------ make factor into photon propogator for testing
c         s34=two*(p(3,4)*p(4,4)-p(3,3)*p(4,3)
c     &        -p(3,2)*p(4,2)-p(3,1)*p(4,1))
c         fac_dm=esq**2/s34**2
c         do j=1,nf
c            dmL(j)=Q(j)
c            dmR(j)=Q(j)
c         enddo
c      else
c         fac_dm=one/dm_lam**4
c      endif



c------ Fix Med-width
c      if((medwidth==1d0).and.(first)) then
c         medwidth=medmass/8d0/pi
c      elseif((medwidth==0d0).and.(first)) then
c         medwidth=medmass/3d0
c      endif
c      if(first.and.(effective_th.eqv..false.)) then
c         write(6,*) 'automatically set medwidth ',medwidth
c      endif

      if(effective_th) then
c--------- effective theory pre-factor
         fac_dm=one/dm_lam**4
      else
c-------- full theory => for V,A and PS is simply g_dmx**2*g_dmq**2/prop(34)**2
c-------- for S and G need to be more careful, but G happens elsewhere and S can be done
c-------- in special routine
         s34=(p(3,4)+p(4,4))**2
     &        -(p(3,1)+p(4,1))**2
     &        -(p(3,2)+p(4,2))**2
     &        -(p(3,3)+p(4,3))**2
         cprop=cone/cplx2((s34-medmass**2),medmass*medwidth)
         propsq=abs(cprop)**2
         fac_dm=propsq*g_dmq**2*g_dmx**2
      endif



      scheme='dred'
      msq(:,:)=0d0
      qqbg_sum(:)=0d0
      qbqg_sum(:)=0d0
      qgqb_sum(:)=0d0
      qbgq_sum(:)=0d0
      gqqb_sum(:)=0d0
      gqbq_sum(:)=0d0
      qgqb_nf(:)=0d0
      qbgq_nf(:)=0d0
      gqqb_nf(:)=0d0
      gqbq_nf(:)=0d0
      qqbg_nf(:)=0d0
      qbqg_nf(:)=0d0
c------ lowest order

      call qqb_dm_monojet(p,msq0)
c------ UV counterterm contains the finite renormalization to arrive at the MS bar scheme (V,A check for S)

      subuv=ason2pi*xn
     & *(epinv*(11._dp-two*real(nflav,dp)/xn)-1._dp)/6._dp

      fac=8d0*cf*xnsq*gsq*fac_dm


      if(dm_mediator=='vector') then
         call qqb_dm_monojet_v_Vamps(p,1,2,5,3,4,qgqb)
         call qqb_dm_monojet_v_Vamps(p,5,2,1,3,4,qbgq)
         call qqb_dm_monojet_v_Vamps(p,1,5,2,3,4,qqbg)
         call qqb_dm_monojet_v_Vamps(p,2,5,1,3,4,qbqg)
         call qqb_dm_monojet_v_Vamps(p,2,1,5,3,4,gqqb)
         call qqb_dm_monojet_v_Vamps(p,5,1,2,3,4,gqbq)
      elseif(dm_mediator=='axvect') then
         call qqb_dm_monojet_v_Axamps(p,1,2,5,3,4,qgqb)
         call qqb_dm_monojet_v_Axamps(p,5,2,1,3,4,qbgq)
         call qqb_dm_monojet_v_Axamps(p,1,5,2,3,4,qqbg)
         call qqb_dm_monojet_v_Axamps(p,2,5,1,3,4,qbqg)
         call qqb_dm_monojet_v_Axamps(p,2,1,5,3,4,gqqb)
         call qqb_dm_monojet_v_Axamps(p,5,1,2,3,4,gqbq)
c------ nf axial pieces
         call qqb_dm_monojet_nf_ax(p,1,2,5,3,4,qgqb_nf)
         call qqb_dm_monojet_nf_ax(p,5,2,1,3,4,qbgq_nf)
         call qqb_dm_monojet_nf_ax(p,1,5,2,3,4,qqbg_nf)
         call qqb_dm_monojet_nf_ax(p,2,5,1,3,4,qbqg_nf)
         call qqb_dm_monojet_nf_ax(p,2,1,5,3,4,gqqb_nf)
         call qqb_dm_monojet_nf_ax(p,5,1,2,3,4,gqbq_nf)
      elseif(dm_mediator=='scalar') then
         call qqb_dm_monojet_v_Samps(p,1,2,5,3,4,qgqb)
         call qqb_dm_monojet_v_Samps(p,5,2,1,3,4,qbgq)
         call qqb_dm_monojet_v_Samps(p,1,5,2,3,4,qqbg)
         call qqb_dm_monojet_v_Samps(p,2,5,1,3,4,qbqg)
         call qqb_dm_monojet_v_Samps(p,2,1,5,3,4,gqqb)
         call qqb_dm_monojet_v_Samps(p,5,1,2,3,4,gqbq)
         fac=fac/4d0
      elseif(dm_mediator=='pseudo') then
         call qqb_dm_monojet_v_PSamps(p,1,2,5,3,4,qgqb)
         call qqb_dm_monojet_v_PSamps(p,5,2,1,3,4,qbgq)
         call qqb_dm_monojet_v_PSamps(p,1,5,2,3,4,qqbg)
         call qqb_dm_monojet_v_PSamps(p,2,5,1,3,4,qbqg)
         call qqb_dm_monojet_v_PSamps(p,2,1,5,3,4,gqqb)
         call qqb_dm_monojet_v_PSamps(p,5,1,2,3,4,gqbq)
         fac=fac/4d0
      elseif(dm_mediator=='gluonO') then
         call gg_dm_monojet_v(p,msq)
         return
      endif




c----- note nf pieces are linear in couplings since
c----- only come from tree

      do j=1,nf
         qqbg_sum(j)=qqbg_sum(j)
     &        +dmL(j)**2*qqbg(1)
     &        +dmR(j)**2*qqbg(2)
     &        +dmL(j)*qqbg_nf(1)+dmR(j)*qqbg_nf(2)
         qbqg_sum(j)=qbqg_sum(j)
     &        +dmL(j)**2*qbqg(1)
     &        +dmR(j)**2*qbqg(2)
     &        +dmL(j)*qbqg_nf(1)+dmR(j)*qbqg_nf(2)
         qgqb_sum(j)=qgqb_sum(j)
     &        +dmL(j)**2*qgqb(1)
     &        +dmR(j)**2*qgqb(2)
     &        +dmL(j)*qgqb_nf(1)+dmR(j)*qgqb_nf(2)
         qbgq_sum(j)=qbgq_sum(j)
     &        +dmL(j)**2*qbgq(1)
     &        +dmR(j)**2*qbgq(2)
     &        +dmL(j)*qbgq_nf(1)+dmR(j)*qbgq_nf(2)
         gqqb_sum(j)=gqqb_sum(j)
     &        +dmL(j)**2*gqqb(1)
     &        +dmR(j)**2*gqqb(2)
     &        +dmL(j)*gqqb_nf(1)+dmR(j)*gqqb_nf(2)
         gqbq_sum(j)=gqbq_sum(j)
     &        +dmL(j)**2*gqbq(1)
     &        +dmR(j)**2*gqbq(2)
     &        +dmL(j)*gqqb_nf(1)+dmR(j)*gqqb_nf(2)
      enddo



      do j=-nf,nf
         do k=-nf,nf
           if( (j  /=  0) .and. (k  /=  0) .and. (j  /=  -k)) goto 20

           if((j==0).and.(k==0)) then
              msq(j,k)=0d0
           elseif((j>0).and.(k<0)) then
              msq(j,k)=aveqq*qqbg_sum(j)*fac
     &            -subuv*msq0(j,k)
           elseif((j>0).and.(k==0)) then
              msq(j,k)=aveqg*qgqb_sum(j)*fac
     &            -subuv*msq0(j,k)
           elseif((j==0).and.(k>0)) then
              msq(j,k)=aveqg*gqqb_sum(k)*fac
     &            -subuv*msq0(j,k)
           elseif((j==0).and.(k<0)) then
              msq(j,k)=aveqg*gqbq_sum(abs(k))*fac
     &            -subuv*msq0(j,k)
           elseif((j<0).and.(k>0)) then
              msq(j,k)=aveqq*qbqg_sum(k)*fac
     &            -subuv*msq0(j,k)
           elseif((j<0).and.(k==0)) then
              msq(j,k)=aveqg*qbgq_sum(abs(j))*fac
     &            -subuv*msq0(j,k)
           endif


 20        continue
        enddo
      enddo


      return
      end


