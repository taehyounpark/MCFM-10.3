!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_dm_monojet(p,msq)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'dm_params.f'
      include 'qcdcouple.f'
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf)
      complex(dp):: qqbg(2,2,2,2),qbqg(2,2,2,2)
      complex(dp):: qgqb(2,2,2,2),qbgq(2,2,2,2)
      complex(dp):: gqqb(2,2,2,2),gqbq(2,2,2,2)
      real(dp):: qqbg_sum(nf),qbqg_sum(nf)
      real(dp):: qgqb_sum(nf),qbgq_sum(nf)
      real(dp):: gqbq_sum(nf),gqqb_sum(nf)
      real(dp):: fac
      integer:: j,k
      integer:: h2,h3,h4
      complex(dp):: cprop
      real(dp):: propsq
      real(dp):: s34
c      logical:: check_QED
c      common/check_QED/check_QED

c      check_QED=.false.

c      if(check_QED) then
c------ make factor into photon propogator for testing
c         s34=2d0*(p(3,4)*p(4,4)-p(3,3)*p(4,3)
c     &        -p(3,2)*p(4,2)-p(3,1)*p(4,1))
c         fac=gsq*V*esq**2/s34**2*4d0
c         do j=1,nf
c            dmL(j)=Q(j)
c            dmR(j)=Q(j)
c         enddo
c      else
c         fac=4d0/dm_lam**4*gsq*V
c      endif

c------ Fix Med-width
c      if((medwidth==1d0).and.(first)) then
c         medwidth=medmass/8d0/pi
c        elseif((medwidth==0d0).and.(first)) then
c       medwidth=medmass/3d0
c      endif
c      if(first.and.(effective_th.eqv..false.)) then
c     write(6,*) 'automatically set medwidth ',medwidth
c      endif

      if(effective_th) then
c--------- effective theory pre-factor
         fac=4d0/dm_lam**4*gsq*V
      else
c-------- full theory => for V,A and PS is simply g_dmx**2*g_dmq**2/prop(34)**2
c-------- for S and G need to be more careful, but G happens elsewhere and S can be done
c-------- in special routine
         s34=(p(3,4)+p(4,4))**2
     &        -(p(3,1)+p(4,1))**2
     &        -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
         cprop=cone/cplx2((s34-medmass**2),medmass*medwidth)
         propsq=abs(cprop)**2
         fac=4d0*gsq*V*propsq*g_dmq**2*g_dmx**2
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
         call qqb_dm_monojet_Vamps(p,1,2,5,3,4,qgqb)
         call qqb_dm_monojet_Vamps(p,5,2,1,3,4,qbgq)
         call qqb_dm_monojet_Vamps(p,1,5,2,3,4,qqbg)
         call qqb_dm_monojet_Vamps(p,2,5,1,3,4,qbqg)
         call qqb_dm_monojet_Vamps(p,2,1,5,3,4,gqqb)
         call qqb_dm_monojet_Vamps(p,5,1,2,3,4,gqbq)
      elseif(dm_mediator=='axvect') then
         call qqb_dm_monojet_Axamps(p,1,2,5,3,4,qgqb)
         call qqb_dm_monojet_Axamps(p,5,2,1,3,4,qbgq)
         call qqb_dm_monojet_Axamps(p,1,5,2,3,4,qqbg)
         call qqb_dm_monojet_Axamps(p,2,5,1,3,4,qbqg)
         call qqb_dm_monojet_Axamps(p,2,1,5,3,4,gqqb)
         call qqb_dm_monojet_Axamps(p,5,1,2,3,4,gqbq)
      elseif(dm_mediator=='scalar') then
         call qqb_dm_monojet_Samps(p,1,2,5,3,4,qgqb)
         call qqb_dm_monojet_Samps(p,5,2,1,3,4,qbgq)
         call qqb_dm_monojet_Samps(p,1,5,2,3,4,qqbg)
         call qqb_dm_monojet_Samps(p,2,5,1,3,4,qbqg)
         call qqb_dm_monojet_Samps(p,2,1,5,3,4,gqqb)
         call qqb_dm_monojet_Samps(p,5,1,2,3,4,gqbq)
         fac=fac/4d0
      elseif(dm_mediator=='pseudo') then
         call qqb_dm_monojet_PSamps(p,1,2,5,3,4,qgqb)
         call qqb_dm_monojet_PSamps(p,5,2,1,3,4,qbgq)
         call qqb_dm_monojet_PSamps(p,1,5,2,3,4,qqbg)
         call qqb_dm_monojet_PSamps(p,2,5,1,3,4,qbqg)
         call qqb_dm_monojet_PSamps(p,2,1,5,3,4,gqqb)
         call qqb_dm_monojet_PSamps(p,5,1,2,3,4,gqbq)
         fac=fac/4d0
      elseif(dm_mediator=='gluonO') then
c---------- Fill msq elsewhere
         call gg_dm_monojet(p,msq)
         return
      elseif(dm_mediator=='scalmt') then
c---------- Fill msq elsewhere
         call gg_dm_top(p,msq)
         return
      endif



      do j=1,nf
      do h2=1,2
         do h3=1,2
            do h4=1,2
               qqbg_sum(j)=qqbg_sum(j)
     &              +abs(dmL(j)*qqbg(1,h2,h3,h4))**2
     &              +abs(dmR(j)*qqbg(2,h2,h3,h4))**2
               qbqg_sum(j)=qbqg_sum(j)
     &              +abs(dmL(j)*qbqg(1,h2,h3,h4))**2
     &              +abs(dmR(j)*qbqg(2,h2,h3,h4))**2
               qgqb_sum(j)=qgqb_sum(j)
     &              +abs(dmL(j)*qgqb(1,h2,h3,h4))**2
     &              +abs(dmR(j)*qgqb(2,h2,h3,h4))**2
               qbgq_sum(j)=qbgq_sum(j)
     &              +abs(dmL(j)*qbgq(1,h2,h3,h4))**2
     &              +abs(dmR(j)*qbgq(2,h2,h3,h4))**2
               gqqb_sum(j)=gqqb_sum(j)
     &              +abs(dmL(j)*gqqb(1,h2,h3,h4))**2
     &              +abs(dmR(j)*gqqb(2,h2,h3,h4))**2
               gqbq_sum(j)=gqbq_sum(j)
     &              +abs(dmL(j)*gqbq(1,h2,h3,h4))**2
     &              +abs(dmR(j)*gqbq(2,h2,h3,h4))**2
            enddo
         enddo
      enddo
      enddo


c---- compare to z
c      call spinoru(5,p,za,zb)
c      call zgamps2(1,2,4,3,5,za,zb,ztes)
c      write(6,*) ztes
c      zt2=0d0
c      do j=1,2
c         do k=1,2
c            zt2=zt2+ztes(j,k)
c         enddo
c      enddo
c      write(6,*) zt2,qqbg
c      pause

c     call writeout(p)
c--- compare to paddy
c      call pjfox_dm(2,p,pjf_msq)
c      write(6,*) 'my qbqg ',qbqg
c      write(6,*) 'my qqbg_sum', qqbg
c      write(6,*) 'paddys ',pjf_msq
c      write(6,*) 'ratio ',qqbg/pjf_msq
c      pause

      do j=-nf,nf
         do k=-nf,nf
           if( (j  /=  0) .and. (k  /=  0) .and. (j  /=  -k)) goto 20

           if((j==0).and.(k==0)) then
              msq(j,k)=0d0
           elseif((j>0).and.(k<0)) then
              msq(j,k)=aveqq*qqbg_sum(j)*fac
           elseif((j>0).and.(k==0)) then
              msq(j,k)=aveqg*qgqb_sum(j)*fac
           elseif((j==0).and.(k>0)) then
              msq(j,k)=aveqg*gqqb_sum(k)*fac
           elseif((j==0).and.(k<0)) then
              msq(j,k)=aveqg*gqbq_sum(abs(k))*fac
           elseif((j<0).and.(k>0)) then
              msq(j,k)=aveqq*qbqg_sum(k)*fac
           elseif((j<0).and.(k==0)) then
              msq(j,k)=aveqg*qbgq_sum(abs(j))*fac
           endif


 20        continue
        enddo
      enddo


      return
      end


