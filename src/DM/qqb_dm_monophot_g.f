!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_dm_monophot_g(p,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'dm_params.f'
      include 'qcdcouple.f'
      include 'ewcharge.f'
      include 'ewcouple.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      complex(dp):: qqbg(2,2,2,2,2),qbqg(2,2,2,2,2)
      complex(dp):: qgqb(2,2,2,2,2),qbgq(2,2,2,2,2)
      complex(dp):: gqqb(2,2,2,2,2),gqbq(2,2,2,2,2)
      real(dp):: qqbg_sum(nf),qbqg_sum(nf)
      real(dp):: qgqb_sum(nf),qbgq_sum(nf)
      real(dp):: gqbq_sum(nf),gqqb_sum(nf)
      real(dp):: fac
      integer:: h2,h3,h4,h5
      integer:: j,k

      real(dp):: s34
      complex(dp):: cprop
      real(dp):: propsq

      logical:: check_QED
      common/check_QED/check_QED
!$omp threadprivate(/check_QED/)



      check_QED=.false.

      if(check_QED) then
c------ make factor into photon propogator for testing
         s34=2d0*(p(3,4)*p(4,4)-p(3,3)*p(4,3)
     &        -p(3,2)*p(4,2)-p(3,1)*p(4,1))
         dm_lam=sqrt(s34/esq)
         do j=1,nf
            dmL(j)=Q(j)
            dmR(j)=Q(j)
         enddo
      endif

      if(effective_th) then
         fac=16d0*esq*cf*xn*gsq/dm_lam**4
      else
c--------full theory => for V,A and PS is simply g_dmx**2*g_dmq**2/prop(34)**2
c--------for S and G need to be more careful, but G happens elsewhere and S can be done
c--------in special routine
         s34=(p(3,4)+p(4,4))**2
     &        -(p(3,1)+p(4,1))**2
     &        -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
         cprop=cone/cplx2((s34-medmass**2),medmass*medwidth)
         propsq=abs(cprop)**2
         fac=16d0*esq*cf*xn*gsq*propsq*g_dmq**2*g_dmx**2
      endif


      qqbg_sum(:)=0d0
      qbqg_sum(:)=0d0
      qgqb_sum(:)=0d0
      qbgq_sum(:)=0d0
      gqqb_sum(:)=0d0
      gqbq_sum(:)=0d0

      msq(:,:)=0
c------- wrapper for separte operators
      if(dm_mediator=='vector') then
         call qqb_dm_monophot_g_Vamps(p,1,2,6,5,3,4,qgqb)
         call qqb_dm_monophot_g_Vamps(p,6,2,1,5,3,4,qbgq)
         call qqb_dm_monophot_g_Vamps(p,1,6,2,5,3,4,qqbg)
         call qqb_dm_monophot_g_Vamps(p,2,6,1,5,3,4,qbqg)
         call qqb_dm_monophot_g_Vamps(p,2,1,6,5,3,4,gqqb)
         call qqb_dm_monophot_g_Vamps(p,6,1,2,5,3,4,gqbq)
      elseif(dm_mediator=='axvect') then
         call qqb_dm_monophot_g_Axamps(p,1,2,6,5,3,4,qgqb)
         call qqb_dm_monophot_g_Axamps(p,6,2,1,5,3,4,qbgq)
         call qqb_dm_monophot_g_Axamps(p,1,6,2,5,3,4,qqbg)
         call qqb_dm_monophot_g_Axamps(p,2,6,1,5,3,4,qbqg)
         call qqb_dm_monophot_g_Axamps(p,2,1,6,5,3,4,gqqb)
         call qqb_dm_monophot_g_Axamps(p,6,1,2,5,3,4,gqbq)
      elseif(dm_mediator=='scalar') then
         call qqb_dm_monophot_g_Samps(p,1,2,6,5,3,4,qgqb)
         call qqb_dm_monophot_g_Samps(p,6,2,1,5,3,4,qbgq)
         call qqb_dm_monophot_g_Samps(p,1,6,2,5,3,4,qqbg)
         call qqb_dm_monophot_g_Samps(p,2,6,1,5,3,4,qbqg)
         call qqb_dm_monophot_g_Samps(p,2,1,6,5,3,4,gqqb)
         call qqb_dm_monophot_g_Samps(p,6,1,2,5,3,4,gqbq)
         fac=fac/4d0
      elseif(dm_mediator=='pseudo') then
         call qqb_dm_monophot_g_PSamps(p,1,2,6,5,3,4,qgqb)
         call qqb_dm_monophot_g_PSamps(p,6,2,1,5,3,4,qbgq)
         call qqb_dm_monophot_g_PSamps(p,1,6,2,5,3,4,qqbg)
         call qqb_dm_monophot_g_PSamps(p,2,6,1,5,3,4,qbqg)
         call qqb_dm_monophot_g_PSamps(p,2,1,6,5,3,4,gqqb)
         call qqb_dm_monophot_g_PSamps(p,6,1,2,5,3,4,gqbq)
         fac=fac/4d0
      endif



      do j=1,nf
      do h2=1,2
         do h3=1,2
            do h4=1,2
               do h5=1,2
               qqbg_sum(j)=qqbg_sum(j)
     &              +abs(dmL(j)*qqbg(1,h2,h3,h4,h5))**2
     &              +abs(dmR(j)*qqbg(2,h2,h3,h4,h5))**2
               qbqg_sum(j)=qbqg_sum(j)
     &              +abs(dmL(j)*qbqg(1,h2,h3,h4,h5))**2
     &              +abs(dmR(j)*qbqg(2,h2,h3,h4,h5))**2
               qgqb_sum(j)=qgqb_sum(j)
     &              +abs(dmL(j)*qgqb(1,h2,h3,h4,h5))**2
     &              +abs(dmR(j)*qgqb(2,h2,h3,h4,h5))**2
               qbgq_sum(j)=qbgq_sum(j)
     &              +abs(dmL(j)*qbgq(1,h2,h3,h4,h5))**2
     &              +abs(dmR(j)*qbgq(2,h2,h3,h4,h5))**2
               gqqb_sum(j)=gqqb_sum(j)
     &              +abs(dmL(j)*gqqb(1,h2,h3,h4,h5))**2
     &              +abs(dmR(j)*gqqb(2,h2,h3,h4,h5))**2
               gqbq_sum(j)=gqbq_sum(j)
     &              +abs(dmL(j)*gqbq(1,h2,h3,h4,h5))**2
     &              +abs(dmR(j)*gqbq(2,h2,h3,h4,h5))**2
            enddo
         enddo
         enddo
      enddo
      enddo


      do j=-nf,nf
         do k=-nf,nf
           if( (j  /=  0) .and. (k  /=  0) .and. (j  /=  -k)) goto 20

           if((j==0).and.(k==0)) then
              msq(j,k)=0d0
           elseif((j>0).and.(k<0)) then
              msq(j,k)=aveqq*qqbg_sum(j)*fac*Q(j)**2
           elseif((j>0).and.(k==0)) then
              msq(j,k)=aveqg*qgqb_sum(j)*fac*Q(j)**2
           elseif((j==0).and.(k>0)) then
              msq(j,k)=aveqg*gqqb_sum(k)*fac*Q(k)**2
           elseif((j==0).and.(k<0)) then
              msq(j,k)=aveqg*gqbq_sum(abs(k))*fac*Q(-k)**2
           elseif((j<0).and.(k>0)) then
              msq(j,k)=aveqq*qbqg_sum(k)*fac*Q(k)**2
           elseif((j<0).and.(k==0)) then
              msq(j,k)=aveqg*qbgq_sum(abs(j))*fac*Q(-j)**2
           endif


 20        continue
        enddo
      enddo


      return
      end
