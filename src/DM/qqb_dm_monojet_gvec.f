!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_dm_monojet_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: R.K. Ellis                                               *
c     September, 1999.                                                 *
c     Matrix element for Z production                                  *
c     averaged over initial colours and spins                          *
c     contracted with the vector n(mu) (orthogonal to p5)              *
c     u(-p1)+dbar(-p2)--> g(p5)+ Z^+(l(p3)+a(p4))
c----- CW MODIFIED FOR DM                    *
c***********************************************************************
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'nflav.f'
      include 'dm_params.f'
      integer:: j,k,in
c--in is the label of the parton dotted with n
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: fac,n(4)
c      complex(dp):: prop
c      real(dp):: ans1,ans2
      real(dp):: qqbg(2),qbqg(2)
      real(dp):: gqqb(2),gqbq(2)
      real(dp):: qgqb(2),qbgq(2)
      complex(dp):: cprop
      real(dp):: propsq,s34,fac_dm

      msq(:,:)=0._dp

      if(dm_mediator=='gluonO') then
         call gg_dm_monojet_gvec(p,n,in,msq)
         return
      endif

      call dotem(5,p,s)


c-----Protect from photon pole by cutting off at some value about 10 GeV
c      if (s(3,4) < 4d0*mbsq) return

c      if(check_QED) then
c         dm_lam=sqrt(s(3,4)/esq)
c         do j=1,nf
c            dmL(j)=Q(j)
c            dmR(j)=Q(j)
c         enddo
c      endif

      if(effective_th) then
c--------- effective theory pre-factor
         fac_dm=one/dm_lam**4
      else
c-------- full theory => for V,A and PS is simply g_dmx**2*g_dmq**2/prop(34)**2
c-------- for S and G need to be more careful, but G happens elsewhere and S can be done
c-------- in special routine
         s34=(p(3,4)+p(4,4))**2
     &        -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
         cprop=cone/cmplx((s34-medmass**2),medmass*medwidth,kind=dp)
         propsq=abs(cprop)**2
         fac_dm=propsq*g_dmq**2*g_dmx**2
      endif



      fac=16d0*cf*xn*gsq*fac_dm

c      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)

c      p1p2(:,:)=0d0

      gqqb(:)=0d0
      gqbq(:)=0d0
      qgqb(:)=0d0
      qbgq(:)=0d0
      qqbg(:)=0d0
      qbqg(:)=0d0


      if(dm_mediator=='vector') then
         if (in == 1) then
            call dmmonojn(5,2,1,p,n,gqqb)
            call dmmonojn(2,5,1,p,n,gqbq)
c     p1p2(0,-1)=-aveqg*fac*dmmonojn(5,2,1,p,n)
c     p1p2(0,+1)=-aveqg*fac*dmmonojn(2,5,1,p,n)
         elseif (in == 2) then
            call dmmonojn(1,5,2,p,n,qgqb)
            call dmmonojn(5,1,2,p,n,qbgq)
c     p1p2(+1,0)=-aveqg*fac*dmmonojn(1,5,2,p,n)
c     p1p2(-1,0)=-aveqg*fac*dmmonojn(5,1,2,p,n)
         elseif (in == 5) then
            call dmmonojn(2,1,5,p,n,qbqg)
            call dmmonojn(1,2,5,p,n,qqbg)
c     p1p2(-1,1)=+aveqq*fac*dmmonojn(2,1,5,p,n)
c     p1p2(1,-1)=+aveqq*fac*dmmonojn(1,2,5,p,n)
         endif
      elseif(dm_mediator=='axvect') then
         if (in == 1) then
            call dmmonojn_ax(5,2,1,p,n,gqqb)
            call dmmonojn_ax(2,5,1,p,n,gqbq)
c     p1p2(0,-1)=-aveqg*fac*dmmonojn(5,2,1,p,n)
c     p1p2(0,+1)=-aveqg*fac*dmmonojn(2,5,1,p,n)
         elseif (in == 2) then
            call dmmonojn_ax(1,5,2,p,n,qgqb)
            call dmmonojn_ax(5,1,2,p,n,qbgq)
c     p1p2(+1,0)=-aveqg*fac*dmmonojn(1,5,2,p,n)
c     p1p2(-1,0)=-aveqg*fac*dmmonojn(5,1,2,p,n)
         elseif (in == 5) then
            call dmmonojn_ax(2,1,5,p,n,qbqg)
            call dmmonojn_ax(1,2,5,p,n,qqbg)
c     p1p2(-1,1)=+aveqq*fac*dmmonojn(2,1,5,p,n)
c      p1p2(1,-1)=+aveqq*fac*dmmonojn(1,2,5,p,n)
         endif
      elseif(dm_mediator=='scalar') then
         if (in == 1) then
            call dmmonojn_scal(5,2,1,p,n,gqqb)
            call dmmonojn_scal(2,5,1,p,n,gqbq)
c     p1p2(0,-1)=-aveqg*fac*dmmonojn(5,2,1,p,n)
c     p1p2(0,+1)=-aveqg*fac*dmmonojn(2,5,1,p,n)
         elseif (in == 2) then
            call dmmonojn_scal(1,5,2,p,n,qgqb)
            call dmmonojn_scal(5,1,2,p,n,qbgq)
c     p1p2(+1,0)=-aveqg*fac*dmmonojn(1,5,2,p,n)
c     p1p2(-1,0)=-aveqg*fac*dmmonojn(5,1,2,p,n)
         elseif (in == 5) then
            call dmmonojn_scal(2,1,5,p,n,qbqg)
            call dmmonojn_scal(1,2,5,p,n,qqbg)
c     p1p2(-1,1)=+aveqq*fac*dmmonojn(2,1,5,p,n)
c     p1p2(1,-1)=+aveqq*fac*dmmonojn(1,2,5,p,n)
         endif
         fac=fac/4d0
      elseif(dm_mediator=='pseudo') then
         if (in == 1) then
            call dmmonojn_Pscal(5,2,1,p,n,gqqb)
            call dmmonojn_Pscal(2,5,1,p,n,gqbq)
c     p1p2(0,-1)=-aveqg*fac*dmmonojn(5,2,1,p,n)
c     p1p2(0,+1)=-aveqg*fac*dmmonojn(2,5,1,p,n)
         elseif (in == 2) then
            call dmmonojn_Pscal(1,5,2,p,n,qgqb)
            call dmmonojn_Pscal(5,1,2,p,n,qbgq)
c     p1p2(+1,0)=-aveqg*fac*dmmonojn(1,5,2,p,n)
c     p1p2(-1,0)=-aveqg*fac*dmmonojn(5,1,2,p,n)
         elseif (in == 5) then
            call dmmonojn_Pscal(2,1,5,p,n,qbqg)
            call dmmonojn_Pscal(1,2,5,p,n,qqbg)
c     p1p2(-1,1)=+aveqq*fac*dmmonojn(2,1,5,p,n)
c     p1p2(1,-1)=+aveqq*fac*dmmonojn(1,2,5,p,n)
         endif
         fac=fac/4d0
      endif



      do j=-nflav,nflav
      do k=-nflav,nflav
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) goto 19

      if     ((j == 0) .and. (k == 0)) then
          msq(j,k)=0d0
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=+aveqq*fac*(abs(dmL(j))**2*qqbg(1)
     &              +abs(dmR(j))**2*qqbg(2))
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=+aveqq*fac*(abs(dmL(-j))**2*qbqg(1)
     &              +abs(dmR(-j))**2*qbqg(2))
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=+aveqg*fac*(abs(dmL(j))**2*qgqb(1)
     &              +abs(dmR(j))**2*qgqb(2))
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=+aveqg*fac*(abs(dmL(-j))**2*qbgq(1)
     &              +abs(dmR(-j))**2*qbgq(2))
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=+aveqg*fac*(abs(dmL(k))**2*gqqb(1)
     &              +abs(dmR(k))**2*gqqb(2))
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=+aveqg*fac*(abs(dmL(-k))**2*gqbq(1)
     &              +abs(dmR(-k))**2*gqbq(2))
      endif


   19 continue
      enddo
      enddo

      return
      end


      subroutine dmmonojn(j1,j2,j5,p,n1,qqbg)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'dm_params.f'
      include 'zprods_decl.f'
c----- ME squared for q(j1)+qb(j2)+x(j3)+x~(j4)+g(j5) X n_nu
      integer:: j1,j2,j5,j3,j4
      real(dp):: n1(4),p(mxpart,4),q(mxpart,4)
      complex(dp)::amp(2,2,2),zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: qqbg(2)
c      real(dp):: z1jetn
c----- return as funciton of quark line helicity
      integer:: h1,h2,h3

      j3=3
      j4=4

      call checkndotp(p,n1,j5)

      if(xmass>1d-8) then
c---------generate massless phase space
         call gen_masslessvecs(p,q,j3,j4)
c---------generate spinors
         call spinoru(5,q,za,zb)
         call spinork(5,q,zab,zba,n1)
      else
c--------massless dm can use usual spinoru
         call spinoru(5,p,za,zb)
         call spinork(5,p,zab,zba,n1)

      endif


      do h1=1,5
         do h2=1,5
            s(h1,h2)=Dble(za(h1,h2)*zb(h2,h1))
         enddo
      enddo


c------ bulid amplitudes

c------- helicity conserving

      amp(1,1,2)= (-2*za(j2,j3)*zab(j1,j1)*zb(j4,j1))/s(j1,j5) +
     &  (2*za(j2,j3)*zab(j2,j2)*zb(j4,j1))/s(j2,j5) -
     &  (2*za(j3,j5)*zab(j2,j5)*zb(j4,j1))/s(j2,j5) +
     &  (2*za(j2,j3)*zab(j5,j1)*zb(j5,j4))/s(j1,j5)
      amp(2,1,2)=(2*za(j1,j3)*zab(j1,j1)*zb(j4,j2))/s(j1,j5) -
     &  (2*za(j3,j5)*zab(j1,j5)*zb(j4,j2))/s(j1,j5) -
     &  (2*za(j1,j3)*zab(j2,j2)*zb(j4,j2))/s(j2,j5) +
     &  (2*za(j1,j3)*zab(j5,j2)*zb(j5,j4))/s(j2,j5)
      amp(1,2,1)=(-2*za(j2,j4)*zab(j1,j1)*zb(j3,j1))/s(j1,j5) +
     &  (2*za(j2,j4)*zab(j2,j2)*zb(j3,j1))/s(j2,j5) -
     &  (2*za(j4,j5)*zab(j2,j5)*zb(j3,j1))/s(j2,j5) +
     &  (2*za(j2,j4)*zab(j5,j1)*zb(j5,j3))/s(j1,j5)
      amp(2,2,1)=(2*za(j1,j4)*zab(j1,j1)*zb(j3,j2))/s(j1,j5) -
     &  (2*za(j4,j5)*zab(j1,j5)*zb(j3,j2))/s(j1,j5) -
     &  (2*za(j1,j4)*zab(j2,j2)*zb(j3,j2))/s(j2,j5) +
     &  (2*za(j1,j4)*zab(j5,j2)*zb(j5,j3))/s(j2,j5)

c------ helicity violating
      amp(1,1,1)=(2*xmass*za(j2,j3)*zab(j1,j1)*zb(j3,j1))/
     &   (s(j1,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j2,j3)*zab(j2,j2)*zb(j3,j1))/
     &   (s(j2,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j3,j5)*zab(j2,j5)*zb(j3,j1))/
     &   (s(j2,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j2,j4)*zab(j1,j1)*zb(j4,j1))/
     &   (s(j1,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j2,j4)*zab(j2,j2)*zb(j4,j1))/
     &   (s(j2,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j4,j5)*zab(j2,j5)*zb(j4,j1))/
     &   (s(j2,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j2,j3)*zab(j5,j1)*zb(j5,j3))/
     &   (s(j1,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j2,j4)*zab(j5,j1)*zb(j5,j4))/
     &   (s(j1,j5)*zb(j4,j3))
      amp(2,1,1)= (-2*xmass*za(j1,j3)*zab(j1,j1)*zb(j3,j2))/
     &   (s(j1,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j3,j5)*zab(j1,j5)*zb(j3,j2))/
     &   (s(j1,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j1,j3)*zab(j2,j2)*zb(j3,j2))/
     &   (s(j2,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j1,j4)*zab(j1,j1)*zb(j4,j2))/
     &   (s(j1,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j4,j5)*zab(j1,j5)*zb(j4,j2))/
     &   (s(j1,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j1,j4)*zab(j2,j2)*zb(j4,j2))/
     &   (s(j2,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j1,j3)*zab(j5,j2)*zb(j5,j3))/
     &   (s(j2,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j1,j4)*zab(j5,j2)*zb(j5,j4))/
     &   (s(j2,j5)*zb(j4,j3))
      amp(1,2,2)=(-2*xmass*za(j2,j3)*zab(j1,j1)*zb(j3,j1))/
     &   (s(j1,j5)*za(j3,j4)) +
     &  (2*xmass*za(j2,j3)*zab(j2,j2)*zb(j3,j1))/
     &   (s(j2,j5)*za(j3,j4)) -
     &  (2*xmass*za(j3,j5)*zab(j2,j5)*zb(j3,j1))/
     &   (s(j2,j5)*za(j3,j4)) +
     &  (2*xmass*za(j2,j4)*zab(j1,j1)*zb(j4,j1))/
     &   (s(j1,j5)*za(j3,j4)) -
     &  (2*xmass*za(j2,j4)*zab(j2,j2)*zb(j4,j1))/
     &   (s(j2,j5)*za(j3,j4)) +
     &  (2*xmass*za(j4,j5)*zab(j2,j5)*zb(j4,j1))/
     &   (s(j2,j5)*za(j3,j4)) +
     &  (2*xmass*za(j2,j3)*zab(j5,j1)*zb(j5,j3))/
     &   (s(j1,j5)*za(j3,j4)) -
     &  (2*xmass*za(j2,j4)*zab(j5,j1)*zb(j5,j4))/
     &   (s(j1,j5)*za(j3,j4))
      amp(2,2,2)=(2*xmass*za(j1,j3)*zab(j1,j1)*zb(j3,j2))/
     &   (s(j1,j5)*za(j3,j4)) -
     &  (2*xmass*za(j3,j5)*zab(j1,j5)*zb(j3,j2))/
     &   (s(j1,j5)*za(j3,j4)) -
     &  (2*xmass*za(j1,j3)*zab(j2,j2)*zb(j3,j2))/
     &   (s(j2,j5)*za(j3,j4)) -
     &  (2*xmass*za(j1,j4)*zab(j1,j1)*zb(j4,j2))/
     &   (s(j1,j5)*za(j3,j4)) +
     &  (2*xmass*za(j4,j5)*zab(j1,j5)*zb(j4,j2))/
     &   (s(j1,j5)*za(j3,j4)) +
     &  (2*xmass*za(j1,j4)*zab(j2,j2)*zb(j4,j2))/
     &   (s(j2,j5)*za(j3,j4)) +
     &  (2*xmass*za(j1,j3)*zab(j5,j2)*zb(j5,j3))/
     &   (s(j2,j5)*za(j3,j4)) -
     &  (2*xmass*za(j1,j4)*zab(j5,j2)*zb(j5,j4))/
     &   (s(j2,j5)*za(j3,j4))


       qqbg(:)=0d0
       do h1=1,2
          do h2=1,2
             do h3=1,2
                 qqbg(h1)=qqbg(h1)+abs(0.25d0*amp(h1,h2,h3))**2
              enddo
           enddo
        enddo

c        write(6,*) j1,j2,j5
c        write(6,*) one/16d0*abs(amp(2,2,1))**2/s(j3,j4)**2
c     &/z1jetn(j1,j2,j5,p,n1)
c        write(6,*) one/16d0*abs(amp(2,1,2))**2/s(j3,j4)**2
c     &/z1jetn(j2,j1,j5,p,n1)
c        pause

      return
      end


c      function z1jetn(j1,j2,j5,p,n)
c      implicit none
c      include 'types.f'
c     real(dp):: z1jetn

c---calculates the amplitude squared for the process
c   q(p1)+qbar(p2) --> Z(l(p3)+a(p4))+g(p5)
c   contracted with the vector n(mu)
c   before spin/color average
c---overall factor of 16 gs**2*gw**4*xw**2*CF*xn removed
c--note QED propagator included.
c      include 'constants.f'
c      include 'nf.f'
c      include 'mxpart.f'
c      include 'cplx.h'
c      include 'sprods_com.f'

c      integer:: j1,j2,j3,j4,j5
c      real(dp):: n(4),p(mxpart,4),nDn,nDp1,nDp2,nDp3,test
c      j3=3
c      j4=4

c      nDp1=n(4)*p(j1,4)-n(3)*p(j1,3)-n(2)*p(j1,2)-n(1)*p(j1,1)
c      nDp2=n(4)*p(j2,4)-n(3)*p(j2,3)-n(2)*p(j2,2)-n(1)*p(j2,1)
c      nDp3=n(4)*p(j3,4)-n(3)*p(j3,3)-n(2)*p(j3,2)-n(1)*p(j3,1)
c      nDn=n(4)**2-n(3)**2-n(2)**2-n(1)**2

c      call checkndotp(p,n,j5)

c      z1jetn=((nDp1*s(j2,j3)/s(j1,j5)-nDp2*s(j1,j4)/s(j2,j5))**2
c     & +two*(s(j2,j3)*nDp1/s(j1,j5)-s(j1,j4)*nDp2/s(j2,j5))*(nDp2+nDp3)
c     & -(s(j1,4)-s(j2,3))**2*s(j3,j4)*nDn/(four*s(j1,j5)*s(j2,j5))
c     & +(nDp2+nDp3)**2)/s(j3,j4)**2

c      return
c      end




      subroutine dmmonojn_ax(j1,j2,j5,p,n1,qqbg)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'dm_params.f'
      include 'zprods_decl.f'
c----- ME squared for q(j1)+qb(j2)+x(j3)+x~(j4)+g(j5) X n_nu
      integer:: j1,j2,j5,j3,j4
      real(dp):: n1(4),p(mxpart,4),q(mxpart,4)
      complex(dp)::amp(2,2,2),zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: qqbg(2)
c      real(dp):: z1jetn
c----- return as funciton of quark line helicity
      integer:: h1,h2,h3
      real(dp):: bp,beta

      j3=3
      j4=4

      call checkndotp(p,n1,j5)

      if(xmass>1d-8) then
c---------generate massless phase space
         call gen_masslessvecs(p,q,j3,j4)
c---------generate spinors
         call spinoru(5,q,za,zb)
         call spinork(5,q,zab,zba,n1)
      else
c--------massless dm can use usual spinoru
         call spinoru(5,p,za,zb)
         call spinork(5,p,zab,zba,n1)

      endif


      do h1=1,5
         do h2=1,5
            s(h1,h2)=Dble(za(h1,h2)*zb(h2,h1))
         enddo
      enddo

cbp = 1/2(1+beta)
c beta = sqrt(1-4xmass**2/s34)
      beta=sqrt(1d0-4d0*xmass**2/s(j3,j4))
      bp=0.5d0*(one+beta)



c------ bulid amplitudes

c------- helicity conserving

      amp(1,1,2)= (-2*za(j2,j4)*zab(j1,j1)*zb(j3,j1))/s(j1,j5) +
     &  (4*bp*za(j2,j4)*zab(j1,j1)*zb(j3,j1))/s(j1,j5) +
     &  (2*za(j2,j4)*zab(j2,j2)*zb(j3,j1))/s(j2,j5) -
     &  (4*bp*za(j2,j4)*zab(j2,j2)*zb(j3,j1))/s(j2,j5) -
     &  (2*za(j4,j5)*zab(j2,j5)*zb(j3,j1))/s(j2,j5) +
     &  (4*bp*za(j4,j5)*zab(j2,j5)*zb(j3,j1))/s(j2,j5) +
     &  (2*za(j2,j4)*zab(j5,j1)*zb(j5,j3))/s(j1,j5) -
     &  (4*bp*za(j2,j4)*zab(j5,j1)*zb(j5,j3))/s(j1,j5)
      amp(2,1,2)=(2*za(j1,j4)*zab(j1,j1)*zb(j3,j2))/s(j1,j5) -
     &  (4*bp*za(j1,j4)*zab(j1,j1)*zb(j3,j2))/s(j1,j5) -
     &  (2*za(j4,j5)*zab(j1,j5)*zb(j3,j2))/s(j1,j5) +
     &  (4*bp*za(j4,j5)*zab(j1,j5)*zb(j3,j2))/s(j1,j5) -
     &  (2*za(j1,j4)*zab(j2,j2)*zb(j3,j2))/s(j2,j5) +
     &  (4*bp*za(j1,j4)*zab(j2,j2)*zb(j3,j2))/s(j2,j5) +
     &  (2*za(j1,j4)*zab(j5,j2)*zb(j5,j3))/s(j2,j5) -
     &  (4*bp*za(j1,j4)*zab(j5,j2)*zb(j5,j3))/s(j2,j5)
      amp(1,2,1)=(2*za(j2,j3)*zab(j1,j1)*zb(j4,j1))/s(j1,j5) -
     &  (4*bp*za(j2,j3)*zab(j1,j1)*zb(j4,j1))/s(j1,j5) -
     &  (2*za(j2,j3)*zab(j2,j2)*zb(j4,j1))/s(j2,j5) +
     &  (4*bp*za(j2,j3)*zab(j2,j2)*zb(j4,j1))/s(j2,j5) +
     &  (2*za(j3,j5)*zab(j2,j5)*zb(j4,j1))/s(j2,j5) -
     &  (4*bp*za(j3,j5)*zab(j2,j5)*zb(j4,j1))/s(j2,j5) -
     &  (2*za(j2,j3)*zab(j5,j1)*zb(j5,j4))/s(j1,j5) +
     &  (4*bp*za(j2,j3)*zab(j5,j1)*zb(j5,j4))/s(j1,j5)
      amp(2,2,1)=(-2*za(j1,j3)*zab(j1,j1)*zb(j4,j2))/s(j1,j5) +
     &  (4*bp*za(j1,j3)*zab(j1,j1)*zb(j4,j2))/s(j1,j5) +
     &  (2*za(j3,j5)*zab(j1,j5)*zb(j4,j2))/s(j1,j5) -
     &  (4*bp*za(j3,j5)*zab(j1,j5)*zb(j4,j2))/s(j1,j5) +
     &  (2*za(j1,j3)*zab(j2,j2)*zb(j4,j2))/s(j2,j5) -
     &  (4*bp*za(j1,j3)*zab(j2,j2)*zb(j4,j2))/s(j2,j5) -
     &  (2*za(j1,j3)*zab(j5,j2)*zb(j5,j4))/s(j2,j5) +
     &  (4*bp*za(j1,j3)*zab(j5,j2)*zb(j5,j4))/s(j2,j5)

c------ helicity violating
      amp(1,1,1)=(2*xmass*za(j2,j3)*zab(j1,j1)*zb(j3,j1))/(s(j1,j5)*
     &     zb(j4,j3)) -
     &  (2*xmass*za(j2,j3)*zab(j2,j2)*zb(j3,j1))/(s(j2,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j3,j5)*zab(j2,j5)*zb(j3,j1))/(s(j2,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j2,j4)*zab(j1,j1)*zb(j4,j1))/(s(j1,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j2,j4)*zab(j2,j2)*zb(j4,j1))/(s(j2,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j4,j5)*zab(j2,j5)*zb(j4,j1))/(s(j2,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j2,j3)*zab(j5,j1)*zb(j5,j3))/(s(j1,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j2,j4)*zab(j5,j1)*zb(j5,j4))/(s(j1,j5)*zb(j4,j3))
      amp(2,1,1)=(-2*xmass*za(j1,j3)*zab(j1,j1)*zb(j3,j2))/
     &   (s(j1,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j3,j5)*zab(j1,j5)*zb(j3,j2))/
     &   (s(j1,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j1,j3)*zab(j2,j2)*zb(j3,j2))/
     &   (s(j2,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j1,j4)*zab(j1,j1)*zb(j4,j2))/
     &   (s(j1,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j4,j5)*zab(j1,j5)*zb(j4,j2))/
     &   (s(j1,j5)*zb(j4,j3)) +
     &  (2*xmass*za(j1,j4)*zab(j2,j2)*zb(j4,j2))/
     &   (s(j2,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j1,j3)*zab(j5,j2)*zb(j5,j3))/
     &   (s(j2,j5)*zb(j4,j3)) -
     &  (2*xmass*za(j1,j4)*zab(j5,j2)*zb(j5,j4))/
     &   (s(j2,j5)*zb(j4,j3))
      amp(1,2,2)=(2*xmass*za(j2,j3)*zab(j1,j1)*zb(j3,j1))/
     &   (s(j1,j5)*za(j3,j4)) -
     &  (2*xmass*za(j2,j3)*zab(j2,j2)*zb(j3,j1))/
     &   (s(j2,j5)*za(j3,j4)) +
     &  (2*xmass*za(j3,j5)*zab(j2,j5)*zb(j3,j1))/
     &   (s(j2,j5)*za(j3,j4)) +
     &  (2*xmass*za(j2,j4)*zab(j1,j1)*zb(j4,j1))/
     &   (s(j1,j5)*za(j3,j4)) -
     &  (2*xmass*za(j2,j4)*zab(j2,j2)*zb(j4,j1))/
     &   (s(j2,j5)*za(j3,j4)) +
     &  (2*xmass*za(j4,j5)*zab(j2,j5)*zb(j4,j1))/
     &   (s(j2,j5)*za(j3,j4)) -
     &  (2*xmass*za(j2,j3)*zab(j5,j1)*zb(j5,j3))/
     &   (s(j1,j5)*za(j3,j4)) -
     &  (2*xmass*za(j2,j4)*zab(j5,j1)*zb(j5,j4))/
     &   (s(j1,j5)*za(j3,j4))
      amp(2,2,2)=(-2*xmass*za(j1,j3)*zab(j1,j1)*zb(j3,j2))/
     &   (s(j1,j5)*za(j3,j4)) +
     &  (2*xmass*za(j3,j5)*zab(j1,j5)*zb(j3,j2))/
     &   (s(j1,j5)*za(j3,j4)) +
     &  (2*xmass*za(j1,j3)*zab(j2,j2)*zb(j3,j2))/
     &   (s(j2,j5)*za(j3,j4)) -
     &  (2*xmass*za(j1,j4)*zab(j1,j1)*zb(j4,j2))/
     &   (s(j1,j5)*za(j3,j4)) +
     &  (2*xmass*za(j4,j5)*zab(j1,j5)*zb(j4,j2))/
     &   (s(j1,j5)*za(j3,j4)) +
     &  (2*xmass*za(j1,j4)*zab(j2,j2)*zb(j4,j2))/
     &   (s(j2,j5)*za(j3,j4)) -
     &  (2*xmass*za(j1,j3)*zab(j5,j2)*zb(j5,j3))/
     &   (s(j2,j5)*za(j3,j4)) -
     &  (2*xmass*za(j1,j4)*zab(j5,j2)*zb(j5,j4))/
     &   (s(j2,j5)*za(j3,j4))


       qqbg(:)=0d0
       do h1=1,2
          do h2=1,2
             do h3=1,2
                 qqbg(h1)=qqbg(h1)+abs(0.25d0*amp(h1,h2,h3))**2
              enddo
           enddo
        enddo

c        write(6,*) j1,j2,j5
c        write(6,*) one/16d0*abs(amp(2,2,1))**2/s(j3,j4)**2
c     &/z1jetn(j1,j2,j5,p,n1)
c        write(6,*) one/16d0*abs(amp(2,1,2))**2/s(j3,j4)**2
c     &/z1jetn(j2,j1,j5,p,n1)
c        pause

      return
      end
