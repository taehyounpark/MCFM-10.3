!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wwg_gvec(p,n,in,msq)
      implicit none
c----Matrix element for WW production
c----averaged over initial colours and spins
c    contracted with the vector v(mu)
c For nwz=+1
c     u(-p1)+dbar(-p2)--> g(p7)+ W^+(n(p3)+e^+(p4))+ W^-(mu(p5)+nu~(p6))
c---
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      include 'pchoice.f'
      integer in
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4)
      real(dp):: n(4)
      integer:: jk,tjk,polq,minus,mplus,jtype
      real(dp):: qdks(mxpart,4),ave,s34,s56,s127,fac,fac1
      complex(dp):: C(2,2)
      complex(dp):: u_ub(3,2),d_db(3,2),ub_u(3,2),db_d(3,2),
     &               u_g(3,2), d_g(3,2), g_ub(3,2),g_db(3,2),
     &               ub_g(3,2),db_g(3,2),g_u(3,2),g_d(3,2),
     &               amp(5),propwp,propwm,prop12,cprop,A(2)
      parameter(minus=1,mplus=2)
      real(dp), parameter :: mp(nf)=(/-1d0,+1d0,-1d0,+1d0,-1d0/)

      msq(:,:)=zip

c      fac=gw**4
      fac=esq**2/abs(zxw)**2
      fac1=gsq*cf*xn

c----Change the momenta to DKS notation
c   We have --- f(p1) + f'(p2)-->mu^-(p5)+nubar(p6)+e^+(p4)+nu(p3)+g(p7)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)+g(p7)

      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(5,j)
      qdks(4,j)=p(6,j)
      qdks(5,j)=p(4,j)
      qdks(6,j)=p(3,j)
      qdks(7,j)=p(7,j)
      enddo

      call spinoru(7,qdks,za,zb)

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7);s34=s(3,4);s56=s(5,6)
      propwp=cmplx(s(3,4),kind=dp)/cmplx(s34-wmass**2,wmass*wwidth,dp)
      propwm=cmplx(s(5,6),kind=dp)/cmplx(s56-wmass**2,wmass*wwidth,dp)
      prop12=cmplx(s127,kind=dp)/cmplx(s127-zmass**2,zmass*zwidth,dp)
      cprop=propwp*propwm

c-- couplings according to hep-ph/9803250 Eqs. 3.4 and 3.6
c-- First argument is left or right, second is d(1) or u(2)
      do j=1,2
      C(minus,j)=mp(j)*2d0*Q(j)*zxw+prop12*(1d0-2d0*mp(j)*Q(j)*zxw)
      C(mplus,j)=mp(j)*2d0*Q(j)*zxw*(1d0-prop12)
      enddo

      if     (in  ==  1) then
        call a7treen(2,7,3,4,5,6,1,qdks,n,za,zb,2,g_ub)
        call a7treen(7,2,3,4,5,6,1,qdks,n,za,zb,2,g_u)
        call a7treen(2,7,6,5,4,3,1,qdks,n,za,zb,1,g_db)
        call a7treen(7,2,6,5,4,3,1,qdks,n,za,zb,1,g_d)
      elseif (in  ==  2) then
        call a7treen(1,7,3,4,5,6,2,qdks,n,za,zb,2,ub_g)
        call a7treen(7,1,3,4,5,6,2,qdks,n,za,zb,2,u_g)
        call a7treen(1,7,6,5,4,3,2,qdks,n,za,zb,1,db_g)
        call a7treen(7,1,6,5,4,3,2,qdks,n,za,zb,1,d_g)
      elseif (in  ==  7) then
        call a7treen(1,2,3,4,5,6,7,qdks,n,za,zb,2,ub_u)
        call a7treen(2,1,3,4,5,6,7,qdks,n,za,zb,2,u_ub)
        call a7treen(1,2,6,5,4,3,7,qdks,n,za,zb,1,db_d)
        call a7treen(2,1,6,5,4,3,7,qdks,n,za,zb,1,d_db)
      else
        write(6,*) 'Invalid value of in=',in,' in qqb_wwg_gvec.f!'
        stop
      endif

      do j=-4,4
      do k=-4,4

c-- skip gluon-gluon case
      if((j  ==  0).and.(k  ==  0)) goto 19

c-- skip non-diagonal quark flavors except gluon
      if ((j  /=  0 .and. k  /=  0) .and. (j  /=  -k)) goto 19

c--- only do g-q and g-qb for in=1
      if ((in  ==  1) .and. (j  /=  0)) goto 19

c--- only do q-g and qb-g for in=2
      if ((in  ==  2) .and. (k  /=  0)) goto 19

c--- only do q-qb and qb-q for in=7
      if ((in  ==  7) .and. ((j  ==  0) .or. (k  ==  0))) goto 19

      jk=max(j,k)
      ave=aveqq
      if (j  ==  0 .or. k  ==  0) then
          jk=j+k
          ave=aveqg
      endif

      do polq=1,2

c---sum is over diagram type (a) and (b)
      do jtype=1,3
          if    (j  <  0 .and. tau(jk)  ==  -1d0 .and. k  /=  0) then
            amp(jtype)=db_d(jtype,polq)
          elseif(j  <  0 .and. tau(jk)  ==   1d0 .and. k  /=  0) then
            amp(jtype)=ub_u(jtype,polq)
          elseif(j  >  0 .and. tau(jk)  ==  -1d0 .and. k  /=  0) then
            amp(jtype)=d_db(jtype,polq)
          elseif(j  >  0 .and. tau(jk)  ==   1d0 .and. k  /=  0) then
            amp(jtype)=u_ub(jtype,polq)
          elseif(j  ==  0 .and. tau(jk)  ==   1d0 .and. jk  >  0) then
            amp(jtype)=g_u(jtype,polq)
          elseif(j  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  >  0) then
            amp(jtype)=g_d(jtype,polq)
          elseif(j  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  <  0) then
            amp(jtype)=g_ub(jtype,polq)
          elseif(j  ==  0 .and. tau(jk)  ==   1d0 .and. jk  <  0) then
            amp(jtype)=g_db(jtype,polq)
          elseif(k  ==  0 .and. tau(jk)  ==   1d0 .and. jk  >  0) then
            amp(jtype)=u_g(jtype,polq)
          elseif(k  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  >  0) then
            amp(jtype)=d_g(jtype,polq)
          elseif(k  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  <  0) then
            amp(jtype)=ub_g(jtype,polq)
          elseif(k  ==  0 .and. tau(jk)  ==   1d0 .and. jk  <  0) then
            amp(jtype)=db_g(jtype,polq)
          endif
      enddo

c---tjk is equal to 2 (u,c) or 1 (d,s,b)
      tjk=2-mod(abs(jk),2)

c-- include coupling in l1 to account for non-leptonic W decays
      A(polq)=l1*cmplx(fac,kind=dp)*(
     & cprop*(amp(1)+C(polq,tjk)*amp(2))+amp(3)*zxw)
      enddo

      msq(j,k)=fac1*ave*
     & (cdabs(A(minus))**2+cdabs(A(mplus))**2)

   19 continue
      enddo
      enddo

      return
      end

