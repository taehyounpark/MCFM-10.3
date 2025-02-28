!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine gg_hZZgg_gvec(p,n,in,msq)
      implicit none
      include 'types.f'


c     g(-p1)+g(-p2)-->H -->Z(e^-(p3)+e^+(p4))+Z(mu^-(p5)+mu^+(p6))
c                                     +g(p_iglue1=7)+g(p_iglue2=8)

c  in is the label of the momentum contracted with n
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'msq_struc.f'
      include 'qcdcouple.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      integer:: in,i5,i6
      real(dp):: msq(-nf:nf,-nf:nf)
      real(dp):: n(4),p(mxpart,4),hdecay,fac,
     & qqgghn_ab,qqgghn_ba,qqgghn_sym,
     & c1234,c1243,c1423,Asq,s3456
c    & ,p1p2(-1:1,-1:1)
      parameter(i5=7,i6=8)
      msq(:,:)=zip

c---fill dot products
      call spinoru(i6,p,za,zb)

c   Deal with Higgs decay to ZZ
      s3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      hdecay=gwsq**3*zmass**2*four*xw**2/(one-xw)*
     & ( ((l1*l2)**2+(r1*r2)**2)*s(3,5)*s(4,6)
     &  +((r1*l2)**2+(r2*l1)**2)*s(3,6)*s(4,5))
      hdecay=hdecay/((s(3,4)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s(5,6)-zmass**2)**2+(zmass*zwidth)**2)
      hdecay=hdecay/((s3456-hmass**2)**2+(hmass*hwidth)**2)

      Asq=(as/(three*pi))**2/vevsq
      fac=Asq*gsq**2*hdecay

c      p1p2(:,:)=zip

c      write(6,*) 'qq_HZZgg_gvec.f: in',in
c--- NOTE: there seems to be some redundancy in the function calls, e.g.
c--- the entries for (0,-1) = (0,+1). Need to check if this is true
c--- and eliminate where possible

c     The function qqgghn calculates q(p1)+qbar(p2)-->H+g(p3)+g(p4)
c     with p4 contracted with n
c     The function gggghn calculates g(p1)+g(p2)-->H+g(p3)+g(p4)
c     with p1 contracted with n

c--- Note that I have removed all references to p1p2(j,k) since
c--- the appropriate terms are actually used only in their
c--- colour-separated forms

      if     (in == 1) then
        call qqgghn(2,i5,i6,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,+1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_ba,0,+1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_sym,0,+1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,-1)=msq_strucv(igg_ab,0,-1)+msq_strucv(igg_ba,0,-1)
c     &            +msq_strucv(igg_sym,0,-1)
c        p1p2(0,+1)=-aveqg*fac*qqgghn_old(i5,2,i6,1,p,n)
        call qqgghn(i5,2,i6,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,-1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_ba,0,-1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_sym,0,-1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,+1)=msq_strucv(igg_ab,0,+1)+msq_strucv(igg_ba,0,+1)
c     &            +msq_strucv(igg_sym,0,+1)
        call qqgghn(i5,i6,2,1,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,0)=avegg*fac*real(nf,dp)*qqgghn_ab
        msq_strucv(igg_ba,0,0)=avegg*fac*real(nf,dp)*qqgghn_ba
        msq_strucv(igg_sym,0,0)=avegg*fac*real(nf,dp)*qqgghn_sym
c        p1p2(0,0)=half*(msq_strucv(igggg_a,0,0)+msq_strucv(igggg_b,0,0)
c     &                  +msq_strucv(igggg_c,0,0))/two
c     &      +real(nf,dp)*(msq_strucv(igg_ab,0,0)+msq_strucv(igg_ba,0,0)
c     &                  +msq_strucv(igg_sym,0,0))
        call gggghn_amp(1,2,i5,i6,p,n,c1234,c1243,c1423)
        msq_strucv(igggg_a,0,0)=avegg*fac*half*c1234
        msq_strucv(igggg_b,0,0)=avegg*fac*half*c1423
        msq_strucv(igggg_c,0,0)=avegg*fac*half*c1243
c        p1p2(0,0)=avegg*fac*(half*(c1234+c1243+c1423)
c     &                 +real(nf,dp)*qqgghn_old(i5,i6,2,1,p,n))

c        p1p2(0,-1)=-aveqg*fac*qqgghn_old(2,i5,i6,1,p,n)

      elseif (in == 2) then
c        p1p2(+1,0)=-aveqg*fac*qqgghn_old(1,i5,i6,2,p,n)
        call qqgghn(1,i5,i6,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,+1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_ba,+1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_sym,+1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(+1,0)=msq_strucv(igg_ab,+1,0)+msq_strucv(igg_ba,+1,0)
c     &            +msq_strucv(igg_sym,+1,0)
c        p1p2(-1,0)=-aveqg*fac*qqgghn_old(i5,1,i6,2,p,n)
        call qqgghn(i5,1,i6,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,-1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_ba,-1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_sym,-1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(-1,0)=msq_strucv(igg_ab,-1,0)+msq_strucv(igg_ba,-1,0)
c     &            +msq_strucv(igg_sym,-1,0)
        call qqgghn(i6,i5,1,2,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,0)=avegg*fac*real(nf,dp)*qqgghn_ab
        msq_strucv(igg_ba,0,0)=avegg*fac*real(nf,dp)*qqgghn_ba
        msq_strucv(igg_sym,0,0)=avegg*fac*real(nf,dp)*qqgghn_sym
c        p1p2(0,0)=half*(msq_strucv(igggg_a,0,0)+msq_strucv(igggg_b,0,0)
c     &                  +msq_strucv(igggg_c,0,0))/two
c     &      +real(nf,dp)*(msq_strucv(igg_ab,0,0)+msq_strucv(igg_ba,0,0)
c     &                  +msq_strucv(igg_sym,0,0))
        call gggghn_amp(2,1,i5,i6,p,n,c1234,c1243,c1423)
        msq_strucv(igggg_a,0,0)=avegg*fac*half*c1243
        msq_strucv(igggg_b,0,0)=avegg*fac*half*c1423
        msq_strucv(igggg_c,0,0)=avegg*fac*half*c1234
c        p1p2(0,0)=+avegg*fac*(half*(c1234+c1243+c1324)
c     &                 +real(nf,dp)*qqgghn_old(i5,i6,1,2,p,n))

      elseif (in == i5) then
c        p1p2(1,-1)=+aveqq*fac*qqgghn_old(1,2,i6,i5,p,n)
        call qqgghn(1,2,i6,i5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,+1,-1)=aveqq*fac*half*qqgghn_ba
        msq_strucv(igg_ba,+1,-1)=aveqq*fac*half*qqgghn_ab
        msq_strucv(igg_sym,+1,-1)=aveqq*fac*half*qqgghn_sym
c        p1p2(+1,-1)=msq_strucv(igg_ab,+1,-1)+msq_strucv(igg_ba,+1,-1)
c     &             +msq_strucv(igg_sym,+1,-1)
c        p1p2(-1,1)=+aveqq*fac*qqgghn_old(2,1,i6,i5,p,n)
        call qqgghn(2,1,i6,i5,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,-1,+1)=aveqq*fac*half*qqgghn_ba
        msq_strucv(igg_ba,-1,+1)=aveqq*fac*half*qqgghn_ab
        msq_strucv(igg_sym,-1,+1)=aveqq*fac*half*qqgghn_sym
c        p1p2(-1,+1)=msq_strucv(igg_ab,-1,+1)+msq_strucv(igg_ba,-1,+1)
c     &             +msq_strucv(igg_sym,-1,+1)
       call gggghn_amp(i5,i6,1,2,p,n,c1234,c1243,c1423)
        msq_strucv(igggg_a,0,0)=avegg*fac*half*c1234
        msq_strucv(igggg_b,0,0)=avegg*fac*half*c1423
        msq_strucv(igggg_c,0,0)=avegg*fac*half*c1243
c        p1p2(0,0)=+half*avegg*fac*(c1234+c1243+c1423)

      elseif (in == i6) then
c        p1p2(1,-1)=+aveqq*fac*qqgghn_old(1,2,i5,i6,p,n)
        call qqgghn(1,2,i5,i6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,+1,-1)=aveqq*fac*half*qqgghn_ab
        msq_strucv(igg_ba,+1,-1)=aveqq*fac*half*qqgghn_ba
        msq_strucv(igg_sym,+1,-1)=aveqq*fac*half*qqgghn_sym
c        p1p2(+1,-1)=msq_strucv(igg_ab,+1,-1)+msq_strucv(igg_ba,+1,-1)
c     &             +msq_strucv(igg_sym,+1,-1)
c        p1p2(-1,1)=+aveqq*fac*qqgghn_old(2,1,i5,i6,p,n)
        call qqgghn(2,1,i5,i6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,-1,+1)=aveqq*fac*half*qqgghn_ab
        msq_strucv(igg_ba,-1,+1)=aveqq*fac*half*qqgghn_ba
        msq_strucv(igg_sym,-1,+1)=aveqq*fac*half*qqgghn_sym
c        p1p2(-1,+1)=msq_strucv(igg_ab,-1,+1)+msq_strucv(igg_ba,-1,+1)
c     &             +msq_strucv(igg_sym,-1,+1)
c--- for the qg, gq pieces, note that qbar-g and g-qbar are never used
        call qqgghn(1,i5,2,i6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,+1,0)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,+1,0)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,+1,0)=-aveqg*fac*qqgghn_sym
c        p1p2(+1,0)=msq_strucv(igg_ab,+1,0)+msq_strucv(igg_ba,+1,0)
c     &            +msq_strucv(igg_sym,+1,0)
        call qqgghn(2,i5,1,i6,p,n,qqgghn_ab,qqgghn_ba,qqgghn_sym)
        msq_strucv(igg_ab,0,+1)=-aveqg*fac*qqgghn_ab
        msq_strucv(igg_ba,0,+1)=-aveqg*fac*qqgghn_ba
        msq_strucv(igg_sym,0,+1)=-aveqg*fac*qqgghn_sym
c        p1p2(0,+1)=msq_strucv(igg_ab,0,+1)+msq_strucv(igg_ba,0,+1)
c     &            +msq_strucv(igg_sym,0,+1)
        call gggghn_amp(i6,1,2,i5,p,n,c1234,c1243,c1423)
        msq_strucv(igggg_a,0,0)=avegg*fac*half*c1234
        msq_strucv(igggg_b,0,0)=avegg*fac*half*c1243
        msq_strucv(igggg_c,0,0)=avegg*fac*half*c1423
c        p1p2(0,0)=+half*avegg*fac*(c1234+c1243+c1423)
      endif

      return
      end


