!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine hjetmass_gvec(p,n,in,msq)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'constants.f'
      include 'masses.f'
      include 'hdecaymode.f'
      include 'first.f'
c  in is the label of the momentum contracted with n
      integer j,k,in
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: n(4),p(mxpart,4),hdecay,s34,fac,
     & p1p2(-1:1,-1:1),msqgamgam
      complex(dp):: hjetmass_gg_gvec_n1
      complex(dp):: hjetmass_gg_gvec_n2
      complex(dp):: hjetmass_gg_gvec_n3
      real(dp):: hjetmass_qqghn
      msq(:,:)=0._dp

      s34=(p(3,4)+p(4,4))**2
     & -(p(3,1)+p(4,1))**2-(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2

c   Deal with Higgs decay
      if (hdecaymode == 'tlta') then
          call htautaudecay(p,3,4,hdecay)
      elseif (hdecaymode == 'bqba') then
          call hbbdecay(p,3,4,hdecay)
      elseif (hdecaymode == 'gaga') then
           hdecay=msqgamgam(sqrt(s34))
      else
      write(6,*) 'Unimplemented process in gg_hgg_gvec'
      stop
      endif
      hdecay=hdecay/((s34-hmass**2)**2+(hmass*hwidth)**2)

      fac=hdecay

      do j=-1,+1
      do k=-1,+1
      p1p2(j,k)=0d0
      enddo
      enddo

      if (first) then
          call sushi_bernini(18)
          first = .false.
      end if

c     write (*,*) "check subtractions HTL"
      if (in  ==  1) then
      p1p2(0,-1)=-aveqg*fac*hjetmass_qqghn(2,5,1,p,n)
c     write (*,*) p1p2(0,-1)/(-aveqg*fac*qqghn(2,5,1,p,n))
      p1p2(0,+1)=-aveqg*fac*hjetmass_qqghn(2,5,1,p,n)
c     write (*,*) p1p2(0,+1)/(-aveqg*fac*qqghn(2,5,1,p,n))
      p1p2(0,0)=dreal(hjetmass_gg_gvec_n1(p,n))*fac/256d0/16d0
c     write (*,*) p1p2(0,0)/(avegg*fac*ggghn(5,2,1,p,n))
      elseif (in  ==  2) then
      p1p2(+1,0)=-aveqg*fac*hjetmass_qqghn(1,5,2,p,n)
c     write (*,*) p1p2(+1,0)/(-aveqg*fac*qqghn(1,5,2,p,n))
      p1p2(-1,0)=-aveqg*fac*hjetmass_qqghn(5,1,2,p,n)
c     write (*,*) p1p2(-1,0)/(-aveqg*fac*qqghn(5,1,2,p,n))
      p1p2(0,0)=dreal(hjetmass_gg_gvec_n2(p,n))*fac/256d0/16d0
c     write (*,*) p1p2(0,0)/(avegg*fac*ggghn(1,5,2,p,n))
      elseif (in  ==  5) then
      p1p2(1,-1)=+aveqq*fac*hjetmass_qqghn(1,2,5,p,n)
c     write (*,*) p1p2(1,-1)/(aveqq*fac*qqghn(1,2,5,p,n))
      p1p2(-1,1)=+aveqq*fac*hjetmass_qqghn(2,1,5,p,n)
c     write (*,*) p1p2(-1,1)/(aveqq*fac*qqghn(2,1,5,p,n))
      p1p2(0,0)=dreal(hjetmass_gg_gvec_n3(p,n))*fac/256d0/16d0
c     write (*,*) p1p2(0,0)/(avegg*fac*ggghn(1,2,5,p,n))
      endif
c     write (*,*) "end check subtractions HTL"

      do j=-nf,nf
      do k=-nf,nf
      if     ((j  >  0) .and. (k  ==  -j)) then
          msq(j,k)=p1p2(1,-1)
      elseif ((j  <  0) .and. (k  ==  -j)) then
          msq(j,k)=p1p2(-1,1)
      elseif ((j  ==  0) .and. (k  ==  0)) then
          msq(j,k)=p1p2(0,0)
      elseif ((j  >  0) .and. (k  ==  0)) then
          msq(j,k)=
     &    p1p2(+1,0)
      elseif ((j  <  0) .and. (k  ==  0)) then
          msq(j,k)=
     &    p1p2(-1,0)
      elseif ((j  ==  0) .and. (k  >  0)) then
          msq(j,k)=
     &    p1p2(0,+1)
      elseif ((j  ==  0) .and. (k  <  0)) then
          msq(j,k)=
     &    p1p2(0,-1)
      endif
      enddo
      enddo

      return
      end

      ! taken from ggHg/gg_hg_gvec.f and augmented by mass factor
      function hjetmass_qqghn(j1,j2,j5,p,n)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      integer j1,j2,j5
      real(dp):: hjetmass_qqghn
      real(dp):: Asq,p(mxpart,4),n(4),nDn,nDp1,nDp2,dot,s,t,u
      real(dp):: mDot, aex
      complex(dp):: asmh

      nDp1 = mDot(p(j1,:), n)
      nDp2 = mDot(p(j2,:), n)
      nDn  = mDot(n, n)

      call checkndotp(p,n,j5)

      Asq=(as/(3d0*pi))**2/vevsq

      s=2d0*Dot(p,j1,j2)
      t=2d0*Dot(p,j1,j5)
      u=2d0*Dot(p,j2,j5)

      ! this becomes 1 in HTL:
      aex = (abs(asmh(s,t+u,mt**2) * mt**2))**2

      hjetmass_qqghn=-Asq*gsq*V/2d0*(2d0*(nDp1*u-nDp2*t)**2/s**2
     &                    +0.5d0*nDn*(u+t)**2/s) * aex

      return
      end

      function mDot(p,q)
      include 'types.f'
      real(dp):: mDot,p(4), q(4)
      mDot = p(4)*q(4) - sum(p(1:3)*q(1:3))
      end function

      function hjetmass_gq_gvec_n3(p,next)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: p(mxpart,4), next(4)
      complex(dp):: hjetmass_gq_gvec_n3,squaredamp, aex
      real(dp):: sman, tman, uman
      real(dp):: aex_prefactor
      real(dp):: p1_next, p2_next, p3_next, next_next
      real(dp):: mDot, dot

      sman = 2*dot(p,1,2)
      tman = 2*dot(p,1,5)
      uman = 2*dot(p,2,5)

      p1_next = mDot(p(1,:),next)
      p2_next = mDot(p(2,:),next)
      p3_next = -mDot(p(5,:),next)
      next_next = mDot(next,next)

      aex_prefactor = sqrt(gsq)**3 / 32d0 / pi**2 / sqrt(vevsq)
      aex = 1/sman*aex_prefactor*8d0/3d0
      ! this (squared) becomes 1 in HTL:
      !aex = aex * asmh(sman,tman+uman,mt**2) * mt**2

      squaredamp =
     &  + dconjg(aex)*aex * ( 4.D0*p1_next*p2_next*tman*uman - 2.D0*
     &    p1_next*p3_next*sman*uman - 2.D0*p1_next*p3_next*sman*tman -
     &    2.D0*p1_next**2*uman**2 - 2.D0*p2_next*p3_next*sman*uman - 2.D
     &    0*p2_next*p3_next*sman*tman - 2.D0*p2_next**2*tman**2 - 1.D0/
     &    2.D0*next_next*sman*uman**2 - next_next*sman*tman*uman - 1.D0/
     &    2.D0*next_next*sman*tman**2 )

      hjetmass_gq_gvec_n3 = squaredamp * 8d0 / 2d0

      end function

      function hjetmass_gg_gvec_n1(p,next)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: p(mxpart,4), next(4)
      complex(dp):: hjetmass_gg_gvec_n1,squaredamp
      real(dp):: cex_prefactor
      complex(dp):: c1ex, c2ex, c3ex, c4ex
      complex(dp):: c1smh, c2smh
      real(dp):: sman, tman, uman
      real(dp):: p1_next, p2_next, p3_next, next_next
      real(dp):: mDot, dot

      sman = 2*dot(p,1,2)
      tman = 2*dot(p,1,5)
      uman = 2*dot(p,2,5)

      p1_next = mDot(p(1,:),next)
      p2_next = mDot(p(2,:),next)
      p3_next = -mDot(p(5,:),next)
      next_next = mDot(next,next)

      cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

      c1ex = c1smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c2ex = c2smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c3ex = c2smh(tman,sman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c4ex = c2smh(uman,sman,tman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      squaredamp =
     &  + dconjg(c1ex)*c4ex * (  - 16.D0*p1_next*p2_next*sman**(-1) +
     &    16.D0*p1_next*p3_next*tman**(-1) + 8.D0*p1_next**2*sman**(-1)
     &    *tman**(-1)*uman + 32.D0*p2_next*p3_next*uman**(-1) + 16.D0*
     &    p2_next**2*sman**(-1)*tman*uman**(-1) + 16.D0*p3_next**2*sman
     &    *tman**(-1)*uman**(-1) + 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c1ex)*c3ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c2ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c1ex * ( 16.D0*p1_next*
     &    p2_next*sman**(-1) - 16.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c4ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c2ex)*c3ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman + 32.D0*p2_next*p3_next
     &    *uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*uman**(-1) +
     &    16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0*next_next
     &     )
      squaredamp = squaredamp + dconjg(c2ex)*c2ex * ( 32.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next**2*sman**(-1)*tman**(-1)*
     &    uman - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c1ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c4ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c3ex * (  - 32.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p1_next**2*sman**(-1)*tman**(-1)*
     &    uman - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c2ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman + 32.D0*p2_next*p3_next
     &    *uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*uman**(-1) +
     &    16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0*next_next
     &     )
      squaredamp = squaredamp + dconjg(c3ex)*c1ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c4ex * ( 16.D0*p1_next*
     &    p2_next*sman**(-1) - 16.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c3ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c2ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p1_next*p3_next*tman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c1ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 8.D0*
     &    p1_next**2*sman**(-1)*tman**(-1)*uman + 32.D0*p2_next*p3_next
     &    *uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*uman**(-1) +
     &    16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0*next_next
     &     )

       hjetmass_gg_gvec_n1 = squaredamp * 24 * 16 / sman / tman / uman

      end function

      function hjetmass_gg_gvec_n2(p,next)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      real(dp):: p(mxpart,4), next(4)
      complex(dp):: hjetmass_gg_gvec_n2,squaredamp
      real(dp):: cex_prefactor
      complex(dp):: c1ex, c2ex, c3ex, c4ex
      complex(dp):: c1smh, c2smh
      real(dp):: sman, tman, uman
      real(dp):: p1_next, p2_next, p3_next, next_next
      real(dp):: mDot, dot

      sman = 2*dot(p,1,2)
      tman = 2*dot(p,1,5)
      uman = 2*dot(p,2,5)

      p1_next = mDot(p(1,:),next)
      p2_next = mDot(p(2,:),next)
      p3_next = -mDot(p(5,:),next)
      next_next = mDot(next,next)

      cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

      c1ex = c1smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c2ex = c2smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c3ex = c2smh(tman,sman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c4ex = c2smh(uman,sman,tman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      squaredamp =
     &  + dconjg(c1ex)*c4ex * (  - 8.D0*p1_next*p2_next*sman**(-1) - 8.D
     &    0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c3ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 32.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D
     &    0*next_next )
      squaredamp = squaredamp + dconjg(c1ex)*c2ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c1ex * ( 16.D0*p1_next*
     &    p2_next*sman**(-1) - 16.D0*p2_next*p3_next*uman**(-1) - 8.D0*
     &    p2_next**2*sman**(-1)*tman*uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c4ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 32.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) - 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D
     &    0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c3ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c2ex)*c2ex * ( 32.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c1ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c4ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c3ex * ( 16.D0*p1_next*
     &    p2_next*sman**(-1) - 16.D0*p2_next*p3_next*uman**(-1) - 8.D0*
     &    p2_next**2*sman**(-1)*tman*uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c2ex * ( 8.D0*p1_next*
     &    p2_next*sman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c1ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 32.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D
     &    0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c4ex * (  - 32.D0*p2_next*
     &    p3_next*uman**(-1) + 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c3ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c2ex * (  - 16.D0*p1_next*
     &    p2_next*sman**(-1) + 32.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) - 8.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 16.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D
     &    0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c1ex * (  - 8.D0*p1_next*
     &    p2_next*sman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )

       hjetmass_gg_gvec_n2 = squaredamp * 24 * 16 / sman / tman / uman

      end function

      function hjetmass_gg_gvec_n3(p,next)
      implicit none
      include 'types.f'
      include 'mxpart.f'
      include 'constants.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      complex(dp):: hjetmass_gg_gvec_n3
      real(dp):: p(mxpart,4), next(4)
      complex(dp):: squaredamp
      real(dp):: cex_prefactor
      complex(dp):: c1ex, c2ex, c3ex, c4ex
      complex(dp):: c1smh, c2smh
      real(dp):: sman, tman, uman
      real(dp):: p1_next, p2_next, p3_next, next_next
      real(dp):: mDot, dot

      sman = 2*dot(p,1,2)
      tman = 2*dot(p,1,5)
      uman = 2*dot(p,2,5)

      p1_next = mDot(p(1,:),next)
      p2_next = mDot(p(2,:),next)
      p3_next = -mDot(p(5,:),next)
      next_next = mDot(next,next)

      cex_prefactor = sqrt(gsq)**3 / 4d0 / pi**2 / sqrt(vevsq)

      c1ex = c1smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c2ex = c2smh(sman,tman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c3ex = c2smh(tman,sman,uman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      c4ex = c2smh(uman,sman,tman,mt**2) *
     & cex_prefactor * mt**2 / 16d0

      squaredamp =
     &  + dconjg(c1ex)*c4ex * ( 8.D0*p1_next*p3_next*tman**(-1) - 8.D0*
     &    p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c3ex * (  - 8.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c1ex)*c2ex * (  - 32.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 8.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0
     &    *next_next )
      squaredamp = squaredamp + dconjg(c1ex)*c1ex * (  - 16.D0*p1_next*
     &    p3_next*tman**(-1) - 16.D0*p2_next*p3_next*uman**(-1) - 8.D0*
     &    p3_next**2*sman*tman**(-1)*uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c4ex * ( 8.D0*p1_next*
     &    p3_next*tman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c2ex)*c3ex * (  - 8.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c2ex)*c2ex * (  - 16.D0*p1_next*
     &    p3_next*tman**(-1) - 16.D0*p2_next*p3_next*uman**(-1) - 8.D0*
     &    p3_next**2*sman*tman**(-1)*uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c2ex)*c1ex * (  - 32.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) + 8.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0
     &    *next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c4ex * (  - 32.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) - 8.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0
     &    *next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c3ex * (  - 32.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p3_next**2*sman*tman**(-1)*
     &    uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c3ex)*c2ex * (  - 8.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c3ex)*c1ex * (  - 8.D0*p1_next*
     &    p3_next*tman**(-1) + 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c4ex * (  - 32.D0*p2_next*
     &    p3_next*uman**(-1) + 8.D0*p3_next**2*sman*tman**(-1)*
     &    uman**(-1) - 8.D0*next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c3ex * (  - 32.D0*p1_next*
     &    p2_next*sman**(-1) + 16.D0*p1_next*p3_next*tman**(-1) + 16.D0
     &    *p1_next**2*sman**(-1)*tman**(-1)*uman + 16.D0*p2_next*
     &    p3_next*uman**(-1) + 16.D0*p2_next**2*sman**(-1)*tman*
     &    uman**(-1) - 8.D0*p3_next**2*sman*tman**(-1)*uman**(-1) + 8.D0
     &    *next_next )
      squaredamp = squaredamp + dconjg(c4ex)*c2ex * ( 8.D0*p1_next*
     &    p3_next*tman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )
      squaredamp = squaredamp + dconjg(c4ex)*c1ex * ( 8.D0*p1_next*
     &    p3_next*tman**(-1) - 8.D0*p2_next*p3_next*uman**(-1) )

       hjetmass_gg_gvec_n3 = squaredamp * 24 * 16 / sman / tman / uman

      end function

