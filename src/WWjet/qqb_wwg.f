!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wwg(p,msq)
c Note: LO routine for WW+jet process
c Differs from real radiation in WW process in the following ways:
c       (1) Does not include the possibility of anomalous couplings
c       (2) Is computed in the complex mass scheme
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  f'(p5)+bar{f'}(p6) + n(p3)+ebar(p4)+ g(p7)
c   for the moment --- radiation only from initial line

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zerowidth.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      include 'plabel.f'
      include 'pchoice.f'
      integer:: jk,tjk,polg,polq,minus,mplus,jtype
      parameter(minus=1,mplus=2)
      real(dp):: P(mxpart,4),qdks(mxpart,4),msq(-nf:nf,-nf:nf),
     & ave,s127,fac,fac1,xfac
      complex(dp):: ct(2,2),cs_z(2,2),cs_g(2,2),
     & cgamz(2,2),cz(2,2)
      complex(dp)::u_ub(5,2,2),d_db(5,2,2),ub_u(5,2,2),db_d(5,2,2),
     &               u_g(5,2,2), d_g(5,2,2), g_ub(5,2,2),g_db(5,2,2),
     &               ub_g(5,2,2),db_g(5,2,2),g_u(5,2,2),g_d(5,2,2),
     &               amp(5),prop12,cprop,A(2,2)
      real(dp), parameter :: mp(nf)=(/-1._dp,+1._dp,-1._dp,+1._dp,-1._dp/)

      msq(:,:)=0._dp

c      fac=gw**4
      fac=abs(zesq/zxw)**2
      fac1=two*gsq*cf
c---multiply by factor for c-sbar+u-dbar hadronic decay
      if (plabel(5) == 'qj') fac1=2._dp*xn*fac1

c----Change the momenta to DKS notation
c   swapped possibility if we want to swap momenta for hadronic case
c   We have --- f(p1) + f'(p2)-->mu^-(p3)+nubar(p4)+e^+(p6)+nu(p5)+g(p7)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)+g(p7)
c----
c   or normal configuration
c   We have --- f(p1) + f'(p2)-->mu^-(p5)+nubar(p6)+e^+(p4)+nu(p3)+g(p7)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)+g(p7)

      if ((plabel(5) == 'qj') .and. (plabel(3) == 'el')) then
c----swapped case
      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(3,j)
      qdks(4,j)=p(4,j)
      qdks(5,j)=p(6,j)
      qdks(6,j)=p(5,j)
      qdks(7,j)=p(7,j)
      enddo
      else
c----normal case
      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(5,j)
      qdks(4,j)=p(6,j)
      qdks(5,j)=p(4,j)
      qdks(6,j)=p(3,j)
      qdks(7,j)=p(7,j)
      enddo
      endif

c--   s returned from sprodx (common block) is 2*dot product
      call spinoru(7,qdks,za,zb)

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7)
      prop12=s127/cplx2(s127-zmass**2,zmass*zwidth)
      cprop=1._dp

c-- couplings according to 3.4 and 3.6
      do j=1,2
      ct(minus,j)=1._dp
      ct(mplus,j)=0._dp
      cs_z(minus,j)=+mp(j)*zL(j)*zsin2w*prop12
      cs_z(mplus,j)=-mp(j)*2._dp*Q(j)*zxw*prop12
      cs_g(minus,j)=+mp(j)*2._dp*Q(j)*zxw
      cs_g(mplus,j)=+mp(j)*2._dp*Q(j)*zxw
      cz(minus,j)=0._dp
      cz(mplus,j)=0._dp
      cgamz(minus,j)=0._dp
      cgamz(mplus,j)=0._dp
c-- couplings with or without photon pole
      if (zerowidth .neqv. .true.) then
      cz(minus,j)=two*zxw*zln*zL(j)*prop12
      cz(mplus,j)=two*zxw*zln*zR(j)*prop12
      cgamz(minus,j)=two*zxw*(-Q(j)+zle*zL(j)*prop12)
      cgamz(mplus,j)=two*zxw*(-Q(j)+zle*zR(j)*prop12)
      endif
      enddo

c--
c      l(j)=(tau(j)-two*Q(j)*xw)/sin2w ; r(j)=(-two*Q(j)*xw)/sin2w
c      le=(-1._dp-two*(-1._dp)*xw)/sin2w ; re=(-two*(-1._dp)*xw)/sin2w
c      ln=(+1._dp-two*(+0._dp)*xw)/sin2w ; rn=0._dp
c---

c--- apply a dipole form factor to anomalous couplings (only if tevscale > 0)
      if (tevscale > 0._dp) then
        xfac=1._dp/(1._dp+s127/(tevscale*1d3)**2)**2
      else
        xfac=1._dp
      endif
      xdelg1_z=xfac*delg1_z
      xdelg1_g=xfac*delg1_g
      xdelk_z=xfac*delk_z
      xdelk_g=xfac*delk_g
      xlambda_z=xfac*lambda_z
      xlambda_g=xfac*lambda_g

c      write(6,*)
c      write(6,*) 'starting now'

c---remember ub-u is the basic process.
c---case ubar-u
      call wwgamps(1,2,3,4,5,6,7,za,zb,ub_u)
c---case ubar-g
      call wwgamps(1,7,3,4,5,6,2,za,zb,ub_g)
c---case g-ubar
      call wwgamps(2,7,3,4,5,6,1,za,zb,g_ub)

c---case u-ubar
      call wwgamps(2,1,3,4,5,6,7,za,zb,u_ub)
c---case u-g
      call wwgamps(7,1,3,4,5,6,2,za,zb,u_g)
c---case g-u
      call wwgamps(7,2,3,4,5,6,1,za,zb,g_u)

c---case dbar-d
      call wwgamps(1,2,6,5,4,3,7,za,zb,db_d)
c---case dbar-g
      call wwgamps(1,7,6,5,4,3,2,za,zb,db_g)
c---case g-dbar
      call wwgamps(2,7,6,5,4,3,1,za,zb,g_db)

c---case .e-_dpdbar
      call wwgamps(2,1,6,5,4,3,7,za,zb,d_db)
c---case .e-_dpg
      call wwgamps(7,1,6,5,4,3,2,za,zb,d_g)
c---case g-d
      call wwgamps(7,2,6,5,4,3,1,za,zb,g_d)

c--- fill matrix elements (no bottom-initiated contribution)
      do j=-4,4
      do k=-4,4

c-- skip gluon-gluon case
      if((j == 0).and.(k == 0)) goto 19
c-- skip non-diagonal quark flavors except gluon
      if((j  /=  0 .and. k  /=  0) .and. (j  /=  -k)) goto 19
      jk=max(j,k)
      ave=xn*aveqq
      if (j == 0 .or. k == 0) then
          jk=j+k
          ave=xn*aveqg
      endif
      do polg=1,2
      do polq=1,2
c---sum is over diagram type t,s(Z),e,n,s(photon)
      do jtype=1,5
          if    (j < 0 .and. tau(jk) == -1._dp .and. k  /=  0) then
            amp(jtype)=db_d(jtype,polg,polq)
          elseif(j < 0 .and. tau(jk) ==  1._dp .and. k  /=  0) then
            amp(jtype)=ub_u(jtype,polg,polq)
          elseif(j > 0 .and. tau(jk) == -1._dp .and. k  /=  0) then
            amp(jtype)=d_db(jtype,polg,polq)
          elseif(j > 0 .and. tau(jk) ==  1._dp .and. k  /=  0) then
            amp(jtype)=u_ub(jtype,polg,polq)
          elseif(j == 0 .and. tau(jk) ==  1._dp .and. jk > 0) then
            amp(jtype)=g_u(jtype,polg,polq)
          elseif(j == 0 .and. tau(jk) == -1._dp .and. jk > 0) then
            amp(jtype)=g_d(jtype,polg,polq)
          elseif(j == 0 .and. tau(jk) == -1._dp .and. jk < 0) then
            amp(jtype)=g_ub(jtype,polg,polq)
          elseif(j == 0 .and. tau(jk) ==  1._dp .and. jk < 0) then
            amp(jtype)=g_db(jtype,polg,polq)
          elseif(k == 0 .and. tau(jk) ==  1._dp .and. jk > 0) then
            amp(jtype)=u_g(jtype,polg,polq)
          elseif(k == 0 .and. tau(jk) == -1._dp .and. jk > 0) then
            amp(jtype)=d_g(jtype,polg,polq)
          elseif(k == 0 .and. tau(jk) == -1._dp .and. jk < 0) then
            amp(jtype)=ub_g(jtype,polg,polq)
          elseif(k == 0 .and. tau(jk) ==  1._dp .and. jk < 0) then
            amp(jtype)=db_g(jtype,polg,polq)
        endif
      enddo

c---tjk is equal to 2 (u,c) or 1 (d,s,b)
        tjk=2-mod(abs(jk),2)

      A(polg,polq)=fac*cprop
     & *(ct(polq,tjk)*amp(1)+cs_z(polq,tjk)*amp(2)+cs_g(polq,tjk)*amp(5)
     &  +cz(polq,tjk)*amp(3)+cgamz(polq,tjk)*amp(4))

      enddo
      enddo

      msq(j,k)=fac1*ave*
     & (abs(A(mplus,minus))**2+abs(A(minus,minus))**2
     & +abs(A(mplus,mplus))**2+abs(A(minus,mplus))**2)


   19 continue

      enddo
      enddo

      return
      end




