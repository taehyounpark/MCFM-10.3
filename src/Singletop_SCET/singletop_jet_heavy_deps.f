!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      module singletop_jet2_deps
      use types
      use singletop2_nnlo_vars
      use singletop2_scale_m
      implicit none

      public :: virt_pp, virt_mm, virt_pm, virt_mp
      public :: wtgvecn

      private

      contains

      function wtgvecn(mq,qwidth,ig,is,ie,in,jn,je,jb,p,vec)
      implicit none
      include 'types.f'
      real(dp):: wtgvecn

      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'ewcouple.f'
      include 'masses.f'
      include 'zprods_com.f'
      integer:: ig,is,ie,in,je,jn,jb
      real(dp):: p(mxpart,4),vec(4),prop,fac,mq,qwidth
      complex(dp):: amp
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp) :: gsq
      common/zabprods/zab,zba
!$omp threadprivate(/zabprods/)

      call checkndotp(p,vec,ig)

      call spinoru(7,p,za,zb)
      call spinork(7,p,zab,zba,vec)
      call ampsn(mq,p,ig,is,ie,in,jn,je,jb,amp)

      prop=(real(za(ie,in)*zb(in,ie))-wmass**2)**2+(wmass*wwidth)**2
      prop=prop*((real(za(jn,je)*zb(je,jn))-wmass**2)**2
     & +(wmass*wwidth)**2)
      prop=prop*(
     &(real(za(jn,je)*zb(je,jn)+za(jn,jb)*zb(jb,jn)+za(je,jb)*zb(jb,je))
     & -mq**2)**2+(mq*qwidth)**2)

      if (corr_on_beam == 1) then
          gsq = 4*pi*as_heavy_beam1
      else
          gsq = 4*pi*as_heavy_beam2
      endif

      fac=xn*cf*gsq*gw**8

      wtgvecn=fac*abs(amp)**2/prop

      return
      end

      subroutine ampsn(mq,p,ig,is,ie,in,jn,je,jb,amp)
      implicit none
      include 'types.f'

      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'zprods_com.f'
      complex(dp):: amp
      real(dp):: p(mxpart,4),dot,taugt,taugs,mq
      integer:: ig,is,ie,in,je,jn,jb
      complex(dp):: zab(mxpart,mxpart),zba(mxpart,mxpart)
      common/zabprods/zab,zba
!$omp threadprivate(/zabprods/)

      taugs=two*dot(p,ig,is)
      taugt=two*(dot(p,ig,jn)+dot(p,ig,je)+dot(p,ig,jb))
      amp= za(ig,ie)*za(jn,jb)*zb(is,in)*zb(je,jn)*zab(jn,ig)*
     & taugt**(-1) + za(ig,ie)*za(jn,jb)*zb(is,in)*zb(je,jb)*zab(jb,ig)
     & *taugt**(-1) - za(ie,je)*za(jn,jb)*zb(is,in)*zb(je,jn)*zab(jn,je
     & )*taugt**(-1) - za(ie,je)*za(jn,jb)*zb(is,in)*zb(je,jb)*zab(jb,
     & je)*taugt**(-1) + za(ie,jn)*za(jn,jb)*zb(ig,in)*zb(je,jn)*zab(ig
     & ,is)*taugs**(-1) + za(ie,jn)*za(jn,jb)*zb(is,in)*zb(je,jn)*zab(
     & is,is)*taugs**(-1) - za(ie,jn)*za(jn,jb)*zb(is,in)*zb(je,jn)*
     & zab(jn,jn)*taugt**(-1) - za(ie,jn)*za(jn,jb)*zb(is,in)*zb(je,jb)
     & *zab(jb,jn)*taugt**(-1) + za(ie,jb)*za(jn,jb)*zb(ig,in)*zb(je,jb
     & )*zab(ig,is)*taugs**(-1) - za(ie,jb)*za(jn,jb)*zb(is,in)*zb(je,
     & jn)*zab(jn,jb)*taugt**(-1) + za(ie,jb)*za(jn,jb)*zb(is,in)*zb(je
     & ,jb)*zab(is,is)*taugs**(-1) - za(ie,jb)*za(jn,jb)*zb(is,in)*zb(
     & je,jb)*zab(jb,jb)*taugt**(-1) + za(jn,jb)*zb(is,in)*zab(ie,je)*
     & mq**2*taugt**(-1)

      return
      end

      function virt_pm(mQ,ig,is,ie,in,ic,p,musq)
      implicit none
      include 'types.f'
      complex(dp):: virt_pm

c***********************************************************************
c     Author: F. Tramontano                                            *
c     December, 2008                                                   *
c***********************************************************************
c---- One-loop +- helicity amplitude for W+c production
c---- hQ=+1  with Q(mu)=Q1(mu)+mQ^2/2/Q.g*g(mu)
c---- hg=-1  with s(mu) as gauge vector
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4)) + Qbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ Q(p5)
c----
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'b0.f'
      include 'scheme.f'
      include 'zprods_decl.f'
      integer:: is,ig,ic,ie,in,nu,jpart
      real(dp):: p(mxpart,4),q(mxpart,4)
      real(dp):: dot,taucg,taugs,taucs,msq,mQ
      real(dp):: qsq,tcg,tcs,qsqhat,qgs,qcg,qcs,cgs1
      real(dp):: ddilog,epin,epin2,xlog
      complex(dp):: spm,lopm
      complex(dp):: apm1,apm2,apm3,apm4,apm5
      complex(dp):: fq2m2,fc002,fc00,ffc00,ffc002,fI3me
      complex(dp):: ffcg1,ffcg3,fun1,fun2
      complex(dp):: fL6m1,fL6m2,fL6m3,L6m1,L6m2,L6m3
      complex(dp):: lnrat,C0fa2m,C0fb2m,I3me
      real(dp), intent(in) :: musq
c      scheme='dred'
      msq=mQ**2
      xlog=log(musq/msq)
      epin=epinv+xlog
      epin2=epinv2*epinv+epinv*xlog+half*xlog**2
      taugs=+two*dot(p,is,ig)
      taucs=+two*dot(p,ic,is)
      taucg=+two*dot(p,ic,ig)
      qsqhat=taucs+taucg+taugs
      qsq=qsqhat+msq
      tcg=taucg+msq
      tcs=taucs+msq
      qgs=qsqhat-taugs
      qcg=qsqhat-taucg
      qcs=qsqhat-taucs
      cgs1=taucs-msq*taugs/taucg
      fq2m2=-msq/qsq+lnrat(-qsqhat,msq)*msq*qsqhat/qsq**2
      ffcg1=fun1(tcg,qsq,msq)
      ffcg3=fun2(qsq,tcg,msq)
      fL6m1=L6m1(taucg,taucs,taugs,msq,qsq)
      fL6m2=L6m2(taucg,taucs,taugs,msq,qsq)
      fL6m3=L6m3(taucg,taucs,taugs,msq,qsq)
      fI3me=I3me(msq,taugs,qsq)
      ffc00=fc00(msq,taugs,qsq)
      ffc002=fc002(msq,taugs,qsq)
      do nu=1,4
      do jpart=1,5
      q(jpart,nu)=p(jpart,nu)
      if (jpart==ic) q(ic,nu)=p(ic,nu)-msq/taucg*p(ig,nu)
      enddo
      enddo
      call spinoru(5,q,za,zb)
      apm1=za(ig,ie)*zb(ig,in)/za(is,ic)**2/zb(ig,is)/zb(ig,ic)
      apm2=za(ig,ie)*zb(is,in)*zb(is,ic)/zb(ig,is)
      apm3=za(ig,ie)*zb(is,in)/za(is,ic)/zb(ig,is)
      apm4=za(is,ie)*zb(is,in)/za(is,ic)**2/zb(ig,is)/zb(ig,ic)
      apm5=-apm1-apm4
      spm=  - 3.D0*taucs*taucg**(-1)*taugs*mQ*cf*apm5 - 5.D0/2.D0*taucs
     &    *taucg**(-1)*taugs**2*mQ**(-1)*cf*apm5*fq2m2 - taucs*taucg*
     &    qsq*mQ**(-3)*cf*apm1*fq2m2 - taucs*taucg*qsq*mQ**(-3)*cf*apm5
     &    *fq2m2 - 3.D0/2.D0*taucs*taucg*mQ**(-1)*cf*apm5*fq2m2 - taucs
     &    *taucg*mQ**(-1)*cf*apm5 + taucs*taugs**(-1)*mQ*xn**(-1)*qcs*
     &    apm4*fL6m2 + 1.D0/2.D0*taucs*taugs**(-1)*mQ*xn**(-1)*cgs1*
     &    apm1*fL6m2 + 1.D0/2.D0*taucs*taugs**(-1)*mQ*xn**(-1)*apm3*
     &    fL6m2 - 4.D0*taucs*taugs*mQ**(-1)*cf*apm5*fq2m2 - 3.D0/2.D0*
     &    taucs*mQ*cf*apm5 + 1.D0/2.D0*taucs**2*taucg**(-2)*taugs*mQ*
     &    xn**(-1)*qcs*cgs1**(-1)*apm4*fL6m1 + 1.D0/2.D0*taucs**2*
     &    taucg**(-2)*taugs*mQ*xn**(-1)*qcg*cgs1**(-1)*apm1*fL6m1 + 1.D0
     &    /2.D0*taucs**2*taucg**(-2)*taugs*mQ*xn**(-1)*cgs1**(-1)*apm3*
     &    fL6m1 - 2.D0*taucs**2*taucg**(-1)*taugs*mQ**(-1)*cf*apm5*
     &    fq2m2 + 3.D0/2.D0*taucs**2*taucg**(-1)*mQ*cf*apm5 - 1.D0/2.D0
     &    *taucs**2*taucg**(-1)*mQ*xn**(-1)*apm1*fL6m1 - 1.D0/2.D0*
     &    taucs**2*taugs**(-1)*mQ*xn**(-1)*qcs*cgs1**(-1)*apm4*fL6m2
      spm = spm - 1.D0/2.D0*taucs**2*taugs**(-1)*mQ*xn**(-1)*qcg*
     & cgs1**(-1)*apm1*fL6m2 - 1.D0/2.D0*taucs**2*taugs**(-1)*mQ*
     &    xn**(-1)*cgs1**(-1)*apm3*fL6m2 - 3.D0/2.D0*taucs**2*mQ**(-1)*
     &    cf*apm5*fq2m2 + 1.D0/2.D0*taucs**3*taucg**(-1)*taugs**(-1)*
     &    qsq*mQ**(-1)*cf*apm5*fq2m2 - 1.D0/2.D0*taucs**3*taucg**(-1)*
     &    taugs**(-1)*mQ*cf*apm5*fq2m2 + taucg**(-2)*taugs*tcg*mQ*xn*
     &    apm3*fL6m3 + 2.D0*taucg**(-2)*taugs*qsq*mQ*cf*apm3*ffcg3 - 2.D
     &    0*taucg**(-2)*taugs**2*qsq*mQ*cf*apm4*ffcg3 - 1.D0/2.D0*
     &    taucg**(-1)*taugs*tcg*mQ*xn*apm1*fL6m3 + taucg**(-1)*taugs*
     &    qsq*mQ*xn**(-1)*apm4*ffcg3 - 2.D0*taucg**(-1)*taugs*mQ*cf*qcs
     &    *apm1 + 3.D0*taucg**(-1)*taugs*mQ*cf*qcs*apm4*ffcg3 +
     &    taucg**(-1)*taugs*mQ*cf*apm3*ffcg3 + taucg**(-1)*taugs*mQ**3*
     &    cf*apm5 + taucg**(-1)*taugs**2*mQ*cf*apm1 - taucg**(-1)*
     &    taugs**2*mQ*cf*apm4*ffcg1 - 3.D0/2.D0*taucg**(-1)*taugs**2*mQ
     &    *cf*apm5 - 1.D0/2.D0*taucg**(-1)*taugs**3*mQ**(-1)*cf*apm5*
     &    fq2m2
      spm = spm - 1.D0/2.D0*taucg**(-1)*qsq*mQ*xn**(-1)*apm3*fL6m2 -
     &    taucg**(-1)*xlog*mQ*b0*apm2 + 3.D0/2.D0*taucg**(-1)*mQ*epin*
     &    cf*apm2 + taucg**(-1)*mQ*epin*b0*apm2 - 1.D0/2.D0*taucg**(-1)
     &    *mQ*epin*xn**(-1)*apm2 + 1.D0/2.D0*taucg**(-1)*mQ*epin*xn*
     &    apm2 - 1.D0/2.D0*taucg**(-1)*mQ*epin2*xn**(-1)*apm2 + 3.D0/2.D
     &    0*taucg**(-1)*mQ*epin2*xn*apm2 - taucg**(-1)*mQ*cf*cgs1*apm3
     &     + 5.D0/2.D0*taucg**(-1)*mQ*cf*apm2 + taucg**(-1)*mQ*xn**(-1)
     &    *cgs1*apm3*fL6m2 - 1.D0/2.D0*taucg**(-1)*mQ*xn**(-1)*pisqo6*
     &    apm2 - 1.D0/2.D0*taucg**(-1)*mQ*xn*qcg*apm3*fL6m3 + 1.D0/2.D0
     &    *taucg**(-1)*mQ*xn*pisqo6*apm2 - 1.D0/2.D0*taucg*taugs*
     &    mQ**(-1)*cf*apm5*fq2m2 - taucg*mQ**(-1)*cf*cgs1*apm1 - 1.D0/2.
     &    D0*taugs**(-1)*mQ*xn**(-1)*qcs*cgs1*apm4*fL6m2 - 3.D0/2.D0*
     &    taugs*mQ*cf*apm5 + 1.D0/2.D0*taugs*mQ*xn*qcs*cgs1**(-1)*apm4*
     &    fL6m3 + 1.D0/2.D0*taugs*mQ*xn*qcg*cgs1**(-1)*apm1*fL6m3 + 1.D0
     &    /2.D0*taugs*mQ*xn*cgs1**(-1)*apm3*fL6m3 + 1.D0/2.D0*taugs*mQ*
     &    xn*apm4*fL6m3
      spm = spm - taugs**2*mQ**(-1)*cf*apm5*fq2m2 - ddilog(msq**(-1)*
     &    tcs)*taucg**(-1)*mQ*xn**(-1)*apm2 + ddilog(msq**(-1)*tcg)*
     &    taucg**(-1)*mQ*xn*apm2 + lnrat( - taucs,msq)*taucs*
     &    taucg**(-1)*taugs*tcs**(-1)*mQ**3*xn**(-1)*apm4 + lnrat( -
     &    taucs,msq)*taucs*taucg**(-1)*mQ*xn**(-1)*qcs*apm4 + lnrat( -
     &    taucs,msq)*taucs*taucg**(-1)*mQ*xn**(-1)*apm3 + lnrat( -
     &    taucs,msq)*taucs**2*taucg**(-1)*tcs**(-1)*qsq*mQ*xn**(-1)*
     &    apm1 - lnrat( - taucs,msq)*taucs**2*tcs**(-1)*mQ*xn**(-1)*
     &    apm4 + lnrat( - taucs,msq)*taucg**(-1)*mQ*epin*xn**(-1)*apm2
     &     + lnrat( - taucs,msq)*tcs**(-1)*mQ*xn**(-1)*cgs1*apm3 -
     &    lnrat( - taucs,msq)*tcs**(-1)*mQ**3*xn**(-1)*cgs1*apm4 -
     &    lnrat( - taucs,msq)**2*taucg**(-1)*mQ*xn**(-1)*apm2 - 2.D0*
     &    lnrat( - taucg,msq)*taucg**(-1)*taugs*tcg**(-1)*mQ*cf*qcg*
     &    apm3 - 4.D0*lnrat( - taucg,msq)*taucg**(-1)*taugs*mQ*cf*apm3
     &     + 2.D0*lnrat( - taucg,msq)*taucg**(-1)*taugs**2*tcg**(-1)*mQ
     &    *cf*qcg*apm4
      spm = spm + 4.D0*lnrat( - taucg,msq)*taucg**(-1)*taugs**2*
     & tcg**(-1)*mQ**3*cf*apm4 - lnrat( - taucg,msq)*taucg**(-1)*mQ*
     &    epin*xn*apm2 - lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ*cf*qcg*
     &    apm1 - lnrat( - taucg,msq)*taugs*tcg**(-1)*qsq*mQ*xn**(-1)*
     &    apm4 + 2.D0*lnrat( - taucg,msq)*taugs*tcg**(-1)*mQ**3*cf*apm4
     &     + 3.D0*lnrat( - taucg,msq)*taugs*mQ*cf*apm1 + 2.D0*lnrat( -
     &    taucg,msq)*taugs**2*tcg**(-1)*mQ*cf*apm4 - lnrat( - taucg,msq
     &    )*tcg**(-1)*mQ*cf*qcs*apm3 - 2.D0*lnrat( - taucg,msq)*
     &    tcg**(-1)*mQ**3*cf*cgs1*apm1 - 2.D0*lnrat( - taucg,msq)*
     &    tcg**(-1)*mQ**3*cf*apm3 + lnrat( - taucg,msq)*mQ*xn**(-1)*qcg
     &    *apm1 + lnrat( - taucg,msq)*mQ*xn**(-1)*apm3 + lnrat( - taucg
     &    ,msq)**2*taucg**(-1)*mQ*xn*apm2 + lnrat( - taugs,msq)*taucs*
     &    taucg**(-1)*taugs*mQ*cf*apm1 + 1.D0/2.D0*lnrat( - taugs,msq)*
     &    taucs**2*taucg**(-1)*mQ**(-1)*cf*qcg*apm5 + lnrat( - taugs,
     &    msq)*taucs**2*taucg**(-1)*mQ**(-1)*cf*apm3 - lnrat( - taugs,
     &    msq)*taucs**2*taucg**(-1)*mQ*cf*apm1
      spm = spm + lnrat( - taugs,msq)*taucs**3*taucg**(-1)*mQ**(-1)*cf*
     & apm1 - 1.D0/2.D0*lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*cf*qcg
     &    *apm5 - lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*xn**(-1)*qcs
     &    *apm1 - lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*xn**(-1)*qcs
     &    *apm5 + lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*xn**(-1)*qcg
     &    *apm1 + lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ*xn*apm3 -
     &    lnrat( - taugs,msq)*taucg**(-1)*taugs*mQ**3*cf*apm1 - lnrat(
     &     - taugs,msq)*taucg**(-1)*mQ*epin*xn*apm2 - 1.D0/2.D0*lnrat(
     &     - taugs,msq)*taugs*mQ**(-1)*xn**(-1)*cgs1*apm5 - 2.D0*lnrat(
     &     - taugs,msq)*taugs*mQ*cf*apm1 - 2.D0*lnrat( - taugs,msq)*
     &    taugs*mQ*cf*apm5 + 1.D0/2.D0*lnrat( - taugs,msq)**2*
     &    taucg**(-1)*mQ*xn*apm2 + 2.D0*lnrat( - qsqhat,msq)*taucs*
     &    taucg**(-2)*taugs*qsq**(-1)*mQ*qsqhat*cf*apm3 - lnrat( -
     &    qsqhat,msq)*taucs*taucg**(-2)*taugs**2*qsq**(-1)*mQ*qsqhat*cf
     &    *apm4 - 1.D0/2.D0*lnrat( - qsqhat,msq)*taucs*taucg**(-1)*
     &    taugs*qsq**(-1)*mQ*qsqhat*cf*apm1
      spm = spm - 5.D0/2.D0*lnrat( - qsqhat,msq)*taucs*taucg**(-1)*
     & taugs*qsq**(-1)*mQ*qsqhat*cf*apm4 - lnrat( - qsqhat,msq)*taucs*
     &    taucg**(-1)*tcs*qsq**(-1)*mQ*qsqhat*xn**(-1)*apm1 + lnrat( -
     &    qsqhat,msq)*taucs*taucg**(-1)*qsq**(-1)*mQ*qsqhat*cf*apm3 -
     &    lnrat( - qsqhat,msq)*taucs*taucg**(-1)*qsq**(-1)*mQ*qsqhat*xn
     &    *qcs*apm1 + 2.D0*lnrat( - qsqhat,msq)*taucs*taucg**(-1)*
     &    qsq**(-1)*mQ**3*qsqhat*cf*apm1 - lnrat( - qsqhat,msq)*taucs*
     &    taucg*qsq**(-1)*mQ**(-1)*qsqhat*cf*apm4 + 3.D0/2.D0*lnrat( -
     &    qsqhat,msq)*taucs*qsq**(-1)*mQ*qsqhat*cf*apm1 - 3.D0/2.D0*
     &    lnrat( - qsqhat,msq)*taucs*qsq**(-1)*mQ*qsqhat*cf*apm4 + 1.D0/
     &    2.D0*lnrat( - qsqhat,msq)*taucs**2*taucg**(-1)*taugs*
     &    qsq**(-1)*mQ**(-1)*qsqhat*cf*apm1 + 1.D0/2.D0*lnrat( - qsqhat
     &    ,msq)*taucs**2*taucg**(-1)*taugs*qsq**(-1)*mQ**(-1)*qsqhat*cf
     &    *apm4 - lnrat( - qsqhat,msq)*taucs**2*taucg**(-1)*qsq**(-1)*
     &    mQ**(-1)*qsqhat*cf*apm3 + 3.D0/2.D0*lnrat( - qsqhat,msq)*
     &    taucs**2*taucg**(-1)*qsq**(-1)*mQ*qsqhat*cf*apm1
      spm = spm + 5.D0/2.D0*lnrat( - qsqhat,msq)*taucs**2*taucg**(-1)*
     & qsq**(-1)*mQ*qsqhat*cf*apm4 - 1.D0/2.D0*lnrat( - qsqhat,msq)*
     &    taucs**3*taucg**(-1)*qsq**(-1)*mQ**(-1)*qsqhat*cf*apm1 + 1.D0/
     &    2.D0*lnrat( - qsqhat,msq)*taucs**3*taucg**(-1)*qsq**(-1)*
     &    mQ**(-1)*qsqhat*cf*apm4 + lnrat( - qsqhat,msq)*taucg**(-2)*
     &    taugs*qsq**(-1)*mQ**3*qsqhat*cf*apm3 + 2.D0*lnrat( - qsqhat,
     &    msq)*taucg**(-2)*taugs**2*qsq**(-1)*mQ*qsqhat*cf*apm3 - 2.D0*
     &    lnrat( - qsqhat,msq)*taucg**(-2)*taugs**2*qsq**(-1)*mQ**3*
     &    qsqhat*cf*apm4 - lnrat( - qsqhat,msq)*taucg**(-2)*taugs**3*
     &    qsq**(-1)*mQ*qsqhat*cf*apm4 - 1.D0/2.D0*lnrat( - qsqhat,msq)*
     &    taucg**(-1)*taugs*tcs*tcg*qsq**(-1)*mQ**(-1)*qsqhat*xn**(-1)*
     &    apm1 + 3.D0*lnrat( - qsqhat,msq)*taucg**(-1)*taugs*qsq**(-1)*
     &    mQ*qsqhat*cf*apm3 - 2.D0*lnrat( - qsqhat,msq)*taucg**(-1)*
     &    taugs*qsq**(-1)*mQ**3*qsqhat*cf*apm4 - 1.D0/2.D0*lnrat( -
     &    qsqhat,msq)*taucg**(-1)*taugs*mQ*qsqhat*xn*apm1 - 2.D0*lnrat(
     &     - qsqhat,msq)*taucg**(-1)*taugs**2*qsq**(-1)*mQ*qsqhat*cf*
     &    apm4
      spm = spm + 2.D0*lnrat( - qsqhat,msq)*taucg**(-1)*qsq**(-1)*mQ**3
     & *qsqhat*cf*apm3 - lnrat( - qsqhat,msq)*taucg**(-1)*mQ*qsqhat*
     &    xn**(-1)*apm3 - 1.D0/2.D0*lnrat( - qsqhat,msq)*taugs*
     &    qsq**(-1)*mQ**(-1)*qsqhat*xn**(-1)*cgs1*apm4 - 3.D0/2.D0*
     &    lnrat( - qsqhat,msq)*taugs*qsq**(-1)*mQ*qsqhat*cf*apm1 - 5.D0/
     &    2.D0*lnrat( - qsqhat,msq)*taugs*qsq**(-1)*mQ*qsqhat*cf*apm4
     &     + lnrat( - qsqhat,msq)*qsq**(-1)*mQ*qsqhat*cf*apm3 + 2.D0*
     &    C0fb2m(tcg,msq)*taucs**(-1)*taucg**(-2)*taugs**2*mQ**5*
     &    xn**(-1)*apm3 + C0fb2m(tcg,msq)*taucs**(-1)*taucg**(-1)*taugs
     &    *qsq*mQ**3*xn**(-1)*apm3 - C0fb2m(tcg,msq)*taucs**(-1)*
     &    taucg**(-1)*taugs*mQ**5*xn**(-1)*qcs*apm4 - C0fb2m(tcg,msq)*
     &    taucs*tcg*mQ*xn**(-1)*apm1 - 3.D0*C0fb2m(tcg,msq)*taucg**(-1)
     &    *taugs*mQ**3*xn**(-1)*apm3 - C0fb2m(tcg,msq)*mQ*xn**(-1)*qcs*
     &    apm3 - C0fb2m(tcg,msq)*mQ**3*xn**(-1)*cgs1*apm1 - 2.D0*
     &    C0fb2m(tcg,msq)*mQ**3*xn**(-1)*apm3 + 2.D0*C0fa2m(tcs,qsq,msq
     &    )*taucs**(-1)*taucg**(-2)*taugs*mQ**3*xn**(-1)*qcs*cgs1*apm3
      spm = spm + C0fa2m(tcs,qsq,msq)*taucs**(-1)*taucg**(-2)*taugs*
     & mQ**5*xn**(-1)*qcs**2*apm4 + C0fa2m(tcs,qsq,msq)*taucs**(-1)*
     &    taucg**(-1)*mQ*xn**(-1)*qcs**2*cgs1*apm3 + C0fa2m(tcs,qsq,msq
     &    )*taucs**(-1)*taucg**(-1)*mQ**3*xn**(-1)*qcs*cgs1*apm3 +
     &    C0fa2m(tcs,qsq,msq)*taucs*taucg**(-1)*tcg*mQ*xn**(-1)*qcs*
     &    apm1 + C0fa2m(tcs,qsq,msq)*taucg**(-1)*mQ**3*xn**(-1)*qcs*
     &    cgs1*apm1 + C0fa2m(tcs,qsq,msq)*taucg**(-1)*mQ**3*xn**(-1)*
     &    qcs*apm3 - 36.D0*ffc002*taucs*taucg**(-1)*taugs*mQ*cf*apm5 -
     &    18.D0*ffc002*taucs*taucg**(-1)*mQ*cf*qgs*apm5 + 18.D0*ffc002*
     &    taucs**2*taucg**(-1)*mQ**(-1)*cf*qgs*apm5 + 36.D0*ffc002*
     &    taucs**2*taucg**(-1)*mQ*cf*apm5 - 6.D0*ffc002*taucs**3*
     &    taucg**(-1)*taugs**(-1)*mQ**(-1)*cf*qgs*apm5 - 12.D0*ffc002*
     &    taucs**3*taucg**(-1)*mQ**(-1)*cf*apm5 + 6.D0*ffc002*
     &    taucg**(-1)*taugs*mQ*cf*qgs*apm5 + 12.D0*ffc002*taucg**(-1)*
     &    taugs*mQ**3*cf*apm5 - 12.D0*ffc00*taucs*taucg**(-1)*taugs*mQ*
     &    cf*apm5
      spm = spm + 8.D0*ffc00*taucs*taucg**(-1)*mQ*cf*apm3 - 20.D0*ffc00
     &    *taucs*taucg**(-1)*mQ**3*cf*apm1 + 8.D0*ffc00*taucs*
     &    taugs**(-1)*mQ*cf*apm3 - 8.D0*ffc00*taucs*mQ*cf*apm1 - 8.D0*
     &    ffc00*taucs*mQ*cf*apm5 + 4.D0*ffc00*taucs**2*taucg**(-1)*
     &    taugs**(-1)*mQ*cf*qcg*apm1 + 8.D0*ffc00*taucs**2*taucg**(-1)*
     &    mQ*cf*apm1 + 28.D0*ffc00*taucs**2*taucg**(-1)*mQ*cf*apm5 - 4.D
     &    0*ffc00*taucs**2*taugs**(-1)*mQ**(-1)*cf*apm3 + 12.D0*ffc00*
     &    taucs**2*taugs**(-1)*mQ*cf*apm1 + 4.D0*ffc00*taucs**2*
     &    mQ**(-1)*cf*apm5 - 4.D0*ffc00*taucs**3*taucg**(-1)*
     &    taugs**(-1)*mQ**(-1)*cf*qgs*apm1 - 4.D0*ffc00*taucs**3*
     &    taucg**(-1)*taugs**(-1)*mQ**(-1)*cf*qgs*apm5 - 4.D0*ffc00*
     &    taucs**3*taucg**(-1)*taugs**(-1)*mQ**(-1)*cf*apm3 + 4.D0*
     &    ffc00*taucg**(-1)*taugs*mQ**3*cf*apm1 - 8.D0*ffc00*
     &    taucg**(-1)*mQ**3*cf*apm3 - 4.D0*ffc00*taugs*mQ*cf*apm5 + 2.D0
     &    *ffc00*mQ**(-1)*xn**(-1)*qgs*cgs1*apm5 + 4.D0*ffc00*mQ*
     &    xn**(-1)*cgs1*apm5
      spm = spm + 4.D0*ffc00*mQ**3*cf*apm1 - 2.D0*fI3me*taucs*
     &    taucg**(-1)*taugs*mQ*cf*apm3 - fI3me*taucs*taucg**(-1)*taugs*
     &    mQ**3*xn**(-1)*apm1 - 2.D0*fI3me*taucs*taugs*mQ*cf*apm4 + 2.D0
     &    *fI3me*taucs**2*taucg**(-1)*taugs*mQ*cf*apm4 + 2.D0*fI3me*
     &    taucg**(-2)*taugs*mQ**3*xn*cgs1*apm3 + fI3me*taucg**(-2)*
     &    taugs**2*mQ**5*xn*apm1 + fI3me*taucg**(-1)*taugs*mQ*xn*cgs1*
     &    apm3 + fI3me*taugs*mQ*xn*cgs1*apm4


      lopm =-mQ/taucg*za(ig,ie)*zb(is,in)*zb(is,ic)/zb(ig,is)

c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      spm=spm + lopm * xn/six

c--- include finite counterterm due to FDH scheme
c--- gw = gw * ( 1 - 2*cf * alphas(MSbar) / (2*pi))
      spm=spm - lopm * cf * two

      if (scheme == 'tH-V') then
        spm = spm - lopm * (cf/two + xn/six)
      endif

      virt_pm=spm

      return
      end

      function virt_mp(mQ,ig,is,ie,in,ic,p,musq)
      implicit none
      include 'types.f'
      complex(dp)::virt_mp

c***********************************************************************
c     Author: F. Tramontano                                            *
c     December, 2008                                                   *
c***********************************************************************
c---- One-loop -+ helicity amplitude for W+c production
c---- hQ=-1  with Q(mu)=Q1(mu)+mQ^2/2/Q.g*g(mu)
c---- hg=+1  with s(mu) as gauge vector
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4)) + Qbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ Q(p5)
c----
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'b0.f'
      include 'scheme.f'
      include 'zprods_decl.f'
      integer::is,ig,ic,ie,in,nu,jpart
      real(dp)::p(mxpart,4),q(mxpart,4)
      real(dp)::dot,taucg,taugs,taucs,msq,mQ
      real(dp)::qsq,tcg,tcs,qsqhat,qcg,qcs,cgs1
      real(dp)::ddilog,epin,epin2,xlog
      complex(dp)::smp,lomp
      complex(dp)::amp1,amp2,amp3,amp4
      complex(dp)::ffcg1,ffcg2,ffcg3,ffcs1,ffcs2,
     & fun1,fun2,fun3,fun4
      complex(dp)::fL6m1,fL6m2,fL6m3,L6m1,L6m2,L6m3
      complex(dp)::lnrat,C0fa2m,C0fb2m,I3me
      real(dp), intent(in) :: musq
c      scheme='dred'
      msq=mQ**2
      xlog=log(musq/msq)
      epin=epinv+xlog
      epin2=epinv2*epinv+epinv*xlog+half*xlog**2
      taugs=+two*dot(p,is,ig)
      taucs=+two*dot(p,ic,is)
      taucg=+two*dot(p,ic,ig)
      qsqhat=taucs+taucg+taugs
      qsq=qsqhat+msq
      tcg=taucg+msq
      tcs=taucs+msq
      qcg=qsqhat-taucg
      qcs=qsqhat-taucs
      cgs1=taucs-msq*taugs/taucg
      ffcg2=fun3(tcg,qsq,msq)
      ffcs1=fun4(tcs,qsq,msq)
      ffcg1=fun1(tcg,qsq,msq)
      ffcg3=fun2(qsq,tcg,msq)
      ffcs2=fun2(tcs,qsq,msq)
      fL6m1=L6m1(taucg,taucs,taugs,msq,qsq)
      fL6m2=L6m2(taucg,taucs,taugs,msq,qsq)
      fL6m3=L6m3(taucg,taucs,taugs,msq,qsq)
      do nu=1,4
      do jpart=1,5
      q(jpart,nu)=p(jpart,nu)
      if (jpart==ic) q(ic,nu)=p(ic,nu)-msq/taucg*p(ig,nu)
      enddo
      enddo
      call spinoru(5,q,za,zb)
      amp1=za(ig,ie)*zb(ig,in)/za(ig,is)**2/zb(is,ic)
      amp2=za(is,ie)*zb(is,in)/za(ig,is)**2/zb(is,ic)
      amp3=za(ig,ie)*za(is,ic)*zb(is,in)/za(ig,is)**2
     & /za(ig,ic)/zb(is,ic)
      amp4=za(ie,ic)*zb(ig,in)/za(ig,is) + za(is,ic)*
     & za(ie,ic)*zb(is,in)/za(ig,is)/za(ig,ic)
      smp=  + taucs*taucg**(-1)*taugs*xn**(-1)*amp1 - 1.D0/2.D0*taucs*
     &    taucg*taugs**(-1)*xn**(-1)*amp3*fL6m2 + 1.D0/2.D0*taucs*
     &    taugs**(-1)*xn**(-1)*cgs1*amp1*fL6m2 - 1.D0/2.D0*taucs*
     &    xn**(-1)*amp3*fL6m1 + 1.D0/2.D0*taucs*xn*amp3*fL6m3 - 2.D0*
     &    taucg**(-2)*taugs**2*mQ**2*cf*amp2*ffcg3 + 1.D0/2.D0*
     &    taucg**(-1)*taugs*tcg*xn**(-1)*amp1*fL6m1 + taucg**(-1)*taugs
     &    *tcg*xn**(-1)*amp1 - 1.D0/2.D0*taucg**(-1)*taugs*tcg*xn*amp1*
     &    fL6m3 + taucg**(-1)*taugs*mQ**2*cf*amp2 - 2.D0*taucg**(-1)*
     &    taugs*mQ**2*cf*amp3*ffcg3 - 1.D0/2.D0*taucg**(-1)*taugs*mQ**2
     &    *xn**(-1)*amp2*fL6m2 + taucg**(-1)*taugs*mQ**2*xn*amp2*ffcg3
     &     + taucg**(-1)*taugs*mQ**2*xn*amp2*fL6m3 - taucg**(-1)*taugs*
     &    mQ**2*xn*amp3*fL6m3 - taucg**(-1)*taugs**2*tcs*xn**(-1)*amp1*
     &    ffcs1 - taucg**(-1)*taugs**2*qsq**(-1)*mQ**2*cf*amp2*ffcg1 -
     &    taucg**(-1)*taugs**2*cf*amp2*ffcg2 - 2.D0*taucg**(-1)*
     &    taugs**2*xn**(-1)*amp1*ffcs2 + taucg*cf*amp2 + taucg*xn**(-1)
     &    *amp2
      smp = smp - 1.D0/2.D0*taucg*xn**(-1)*amp3*fL6m1 + 1.D0/2.D0*taucg
     &    *xn*amp3*fL6m3 - 3.D0*taugs*tcg*qsq**(-1)*cf*amp3*ffcg3 -
     &    taugs*tcg*qsq**(-1)*xn**(-1)*amp2*ffcg3 + 2.D0*taugs*tcg*
     &    qsq**(-1)*xn*amp2*ffcg3 + 2.D0*taugs*cf*amp2 + taugs*xn**(-1)
     &    *amp2 + taugs*xn**(-1)*amp3*ffcs2 + taugs*xn**(-1)*amp3 -
     &    taugs**2*qsq**(-1)*cf*amp2*ffcg3 + 1.D0/2.D0*tcs*xn**(-1)*
     &    amp3*fL6m2 + xlog*b0*amp4 - 3.D0/2.D0*epin*cf*amp4 - epin*b0*
     &    amp4 + 1.D0/2.D0*epin*xn**(-1)*amp4 - 1.D0/2.D0*epin*xn*amp4
     &     + 1.D0/2.D0*epin2*xn**(-1)*amp4 - 3.D0/2.D0*epin2*xn*amp4 -
     &    qsqhat*cf*amp2 + 1.D0/2.D0*qsqhat*xn**(-1)*amp2*fL6m1 - 1.D0/
     &    2.D0*qsqhat*xn*amp2*fL6m3 + cf*cgs1*amp3 - 5.D0/2.D0*cf*amp4
     &     - xn**(-1)*qcs*amp3 - xn**(-1)*cgs1*amp3*fL6m2 + 1.D0/2.D0*
     &    xn**(-1)*pisqo6*amp4 - 1.D0/2.D0*xn*pisqo6*amp4 + ddilog(
     &    msq**(-1)*tcs)*xn**(-1)*amp4 - ddilog(msq**(-1)*tcg)*xn*amp4
     &     - I3me(msq,taugs,qsq)*taucs**(-1)*taucg**(-1)*taugs**2*tcg*
     &    mQ**2*xn**(-1)*amp1
      smp = smp + I3me(msq,taugs,qsq)*taucs**(-1)*taucg*taugs*mQ**2*
     & xn**(-1)*amp3 - I3me(msq,taugs,qsq)*taucs**(-1)*taugs*mQ**2*
     &    xn**(-1)*qcs*amp2 + I3me(msq,taugs,qsq)*taucg**(-2)*taugs**2*
     &    tcg*mQ**2*xn*amp1 + I3me(msq,taugs,qsq)*taucg**(-1)*taugs*
     &    mQ**2*xn*qcs*amp2 + 2.D0*I3me(msq,taugs,qsq)*taucg**(-1)*
     &    taugs*mQ**2*xn*cgs1*amp2 - 2.D0*I3me(msq,taugs,qsq)*
     &    taucg**(-1)*taugs*mQ**2*xn*cgs1*amp3 - I3me(msq,taugs,qsq)*
     &    taugs*mQ**2*xn*amp3 + lnrat( - taucs,msq)*taucs*taucg**(-1)*
     &    taugs*xn**(-1)*amp1 + lnrat( - taucs,msq)*taucs*taucg**(-1)*
     &    xn**(-1)*qcg*amp1 - lnrat( - taucs,msq)*taucs*xn**(-1)*amp3
     &     + lnrat( - taucs,msq)*taucg**(-1)*taugs*tcs**(-1)*qsq*mQ**2*
     &    xn**(-1)*amp2 + lnrat( - taucs,msq)*taucg**(-1)*taugs**2*
     &    tcs**(-1)*mQ**2*xn**(-1)*amp1 - lnrat( - taucs,msq)*taugs*
     &    tcs**(-1)*mQ**2*xn**(-1)*amp3 - lnrat( - taucs,msq)*epin*
     &    xn**(-1)*amp4 + lnrat( - taucs,msq)**2*xn**(-1)*amp4 + 2.D0*
     &    lnrat(
     &  - taucg,msq)*taucg**(-1)*taugs*mQ**2*cf*amp1 + 2.D0*lnrat( -
     &    taucg,msq)*taucg**(-1)*taugs*cf*qcs*amp2 + lnrat( - taucg,msq
     &    )*taucg*taugs*tcg**(-1)*xn*amp2 - lnrat( - taucg,msq)*taucg*
     &    xn**(-1)*amp3 + 2.D0*lnrat( - taucg,msq)*taugs*tcg**(-1)*
     &    mQ**2*cf*amp3 + 3.D0*lnrat( - taucg,msq)*taugs*cf*amp1 - 2.D0
     &    *lnrat( - taucg,msq)*taugs**2*tcg**(-1)*cf*amp2 + lnrat( -
     &    taucg,msq)*epin*xn*amp4 + lnrat( - taucg,msq)*xn**(-1)*qcg*
     &    amp1 - lnrat( - taucg,msq)**2*xn*amp4 + lnrat( - taugs,msq)*
     &    epin*xn*amp4 - 1.D0/2.D0*lnrat( - taugs,msq)**2*xn*amp4 -
     &    lnrat( - qsqhat,msq)*taucg**(-2)*tcg*qsq**(-1)*qsqhat*cf*
     &    qcs**2*amp2 + lnrat( - qsqhat,msq)*taucg**(-1)*taugs**2*
     &    qsq**(-1)*qsqhat*cf*amp2 + lnrat( - qsqhat,msq)*taucg**(-1)*
     &    taugs**2*qsq**(-1)*qsqhat*xn**(-1)*amp1 - lnrat( - qsqhat,msq
     &    )*taucg**(-1)*qsqhat*xn**(-1)*qcg*amp1 + lnrat( - qsqhat,msq)
     &    *taucg*qsq**(-1)*qsqhat*cf*amp2 - 3.D0*lnrat( - qsqhat,msq)*
     &    taucg*qsq**(-1)*qsqhat*cf*amp3
      smp = smp + lnrat( - qsqhat,msq)*taucg*qsq**(-1)*qsqhat*xn**(-1)*
     & amp3 - lnrat( - qsqhat,msq)*taugs*qsq**(-1)*qsqhat*xn*amp2 - 2.D0
     &    *lnrat( - qsqhat,msq)*tcs*qsq**(-1)*qsqhat*cf*amp3 + lnrat(
     &     - qsqhat,msq)*tcs*qsq**(-1)*qsqhat*xn**(-1)*amp3 - lnrat( -
     &    qsqhat,msq)*qsq**(-1)*qsqhat*cf*cgs1*amp3 + 3.D0*lnrat( -
     &    qsqhat,msq)*qsqhat*cf*amp2 + C0fb2m(tcg,msq)*taucs**(-1)*
     &    taucg*mQ**2*xn**(-1)*cgs1*amp3 - C0fb2m(tcg,msq)*taucs**(-1)*
     &    taugs*mQ**2*xn**(-1)*cgs1*amp2 + 2.D0*C0fb2m(tcg,msq)*
     &    taucs**(-1)*taugs*mQ**2*xn**(-1)*cgs1*amp3 - C0fb2m(tcg,msq)*
     &    mQ**2*xn**(-1)*cgs1*amp1 - C0fa2m(tcs,qsq,msq)*taucs**(-1)*
     &    taucg**(-2)*taugs**2*mQ**4*xn**(-1)*qcs*amp2 - 2.D0*C0fa2m(
     &    tcs,qsq,msq)*taucs**(-1)*taucg**(-1)*taugs*mQ**2*xn**(-1)*qcs
     &    *cgs1*amp3 - C0fa2m(tcs,qsq,msq)*taucs**(-1)*mQ**2*xn**(-1)*
     &    qcs*cgs1*amp3 - C0fa2m(tcs,qsq,msq)*taucg**(-1)*taugs**2*
     &    mQ**2*xn**(-1)*amp1 + C0fa2m(tcs,qsq,msq)*taucg**(-1)*mQ**2*
     &    xn**(-1)*qcs*cgs1*amp1
      smp = smp + C0fa2m(tcs,qsq,msq)*taugs*mQ**2*xn**(-1)*amp3


      lomp=za(ie,ic)
     & /za(ig,is)/za(ig,ic)*(za(is,ic)*zb(is,in)+za(ig,ic)*zb(ig,in))

c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      smp=smp + lomp * xn/six

c--- include finite counterterm due to FDH scheme
c--- gw = gw * ( 1 - 2*cf * alphas(MSbar) / (2*pi))
      smp=smp - lomp * cf * two

      if (scheme == 'tH-V') then
        smp = smp - lomp * (cf/two + xn/six)
      endif

      virt_mp=smp
      return
      end

      function virt_mm(mQ,ig,is,ie,in,ic,p,musq)
      implicit none
      include 'types.f'
      complex(dp):: virt_mm

c***********************************************************************
c     Author: F. Tramontano                                            *
c     December, 2008                                                   *
c***********************************************************************
c---- One-loop -- helicity amplitude for W+c production
c---- hQ=-1  with Q(mu)=Q1(mu)+mQ^2/2/Q.g*g(mu)
c---- hg=-1  with s(mu) as gauge vector
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4)) + Qbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ Q(p5)
c----
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'b0.f'
      include 'scheme.f'
      include 'zprods_decl.f'
      integer:: is,ig,ic,ie,in,nu,jpart
      real(dp):: p(mxpart,4),q(mxpart,4)
      real(dp):: dot,taucg,taugs,taucs,msq,mQ
      real(dp):: qsq,tcg,tcs,qsqhat,qgs,qcg,qcs,cgs1
      real(dp):: ddilog,epin,epin2,xlog
      complex(dp):: smm,lomm
      complex(dp):: amm1,amm2,amm3,amm4,amm5
      complex(dp):: fq2m2,fc002,fc00,ffc00,ffc002,fI3me
      complex(dp):: ffcg1,ffcg3,ffcs1,ffcs2,fun1,fun2,fun4
      complex(dp):: fL6m1,fL6m2,fL6m3,L6m1,L6m2,L6m3
      complex(dp):: lnrat,C0fa2m,C0fb2m,I3me
      real(dp), intent(in) :: musq
c      scheme='dred'
      msq=mQ**2
      xlog=log(musq/msq)
      epin=epinv+xlog
      epin2=epinv2*epinv+epinv*xlog+half*xlog**2
      taugs=+two*dot(p,is,ig)
      taucs=+two*dot(p,ic,is)
      taucg=+two*dot(p,ic,ig)
      qsqhat=taucs+taucg+taugs
      qsq=qsqhat+msq
      tcg=taucg+msq
      tcs=taucs+msq
      qgs=qsqhat-taugs
      qcg=qsqhat-taucg
      qcs=qsqhat-taucs
      cgs1=taucs-msq*taugs/taucg
      fq2m2=-msq/qsq+lnrat(-qsqhat,msq)*msq*qsqhat/qsq**2
      ffcs1=fun4(tcs,qsq,msq)
      ffcg1=fun1(tcg,qsq,msq)
      ffcg3=fun2(qsq,tcg,msq)
      ffcs2=fun2(tcs,qsq,msq)
      fL6m1=L6m1(taucg,taucs,taugs,msq,qsq)
      fL6m2=L6m2(taucg,taucs,taugs,msq,qsq)
      fL6m3=L6m3(taucg,taucs,taugs,msq,qsq)
      fI3me=I3me(msq,taugs,qsq)
      ffc00=fc00(msq,taugs,qsq)
      ffc002=fc002(msq,taugs,qsq)
      do nu=1,4
      do jpart=1,5
      q(jpart,nu)=p(jpart,nu)
      if (jpart==ic) q(ic,nu)=p(ic,nu)-msq/taucg*p(ig,nu)
      enddo
      enddo
      call spinoru(5,q,za,zb)
      amm1=za(ig,ie)*zb(ig,in)/za(is,ic)/zb(ig,ic)**2
      amm2=za(ig,ie)*zb(is,in)/zb(ig,ic)
      amm3=za(is,ie)*zb(is,in)/za(is,ic)/zb(ig,ic)**2
      amm4=za(ie,ic)*zb(is,in)*zb(is,ic)/zb(ig,is)/zb(ig,ic)
      amm5=-amm1-amm3
      smm=  - 1.D0/2.D0*taucs*taucg**(-2)*mQ**2*xn**(-1)*qcs*qcg*
     &    cgs1**(-1)*amm1*fL6m1 - 1.D0/2.D0*taucs*taucg**(-2)*mQ**2*
     &    xn**(-1)*qcs*cgs1**(-1)*amm2*fL6m1 - 1.D0/2.D0*taucs*
     &    taucg**(-2)*mQ**2*xn**(-1)*qcs**2*cgs1**(-1)*amm3*fL6m1 - 1.D0
     &    /2.D0*taucs*taucg**(-1)*taugs**(-1)*xn**(-1)*qcs*amm2*fL6m1
     &     - 1.D0/2.D0*taucs*taucg**(-1)*taugs**(-1)*xn**(-1)*qcs**2*
     &    amm3*fL6m1 - 1.D0/2.D0*taucs*taucg**(-1)*qsqhat*xn**(-1)*amm1
     &    *fL6m1 + taucs*taucg*taugs**(-1)*cf*amm3 + 1.D0/2.D0*taucs*
     &    taucg*taugs**(-1)*xn*amm3*fL6m3 + 1.D0/2.D0*taucs*taugs**(-1)
     &    *mQ**2*xn**(-1)*amm1*fL6m2 + taucs*taugs**(-1)*cf*qcs*amm5*
     &    fq2m2 + 1.D0/2.D0*taucs*taugs**(-1)*cf*qcg*amm5*fq2m2 - 4.D0*
     &    taucs*qsq**(-1)*mQ**2*cf*amm1 - 2.D0*taucs*qsq**(-1)*mQ**2*
     &    xn**(-1)*amm1 + 2.D0*taucs*qsq**(-1)*mQ**2*xn*amm1 + 4.D0*
     &    taucs*qsq**(-1)*cf*cgs1*amm1 + 2.D0*taucs*qsq**(-1)*xn**(-1)*
     &    cgs1*amm1 - 2.D0*taucs*qsq**(-1)*xn*cgs1*amm1 - 4.D0*taucs*cf
     &    *amm1*fq2m2
      smm = smm + 3.D0*taucs*cf*amm1 + taucs*cf*amm3 - 2.D0*taucs*
     &    xn**(-1)*amm1*ffcs2 - 2.D0*taucs*xn**(-1)*amm1*fq2m2 + 2.D0*
     &    taucs*xn*amm1*fq2m2 - taucs*xn*amm1 - 1.D0/2.D0*taucs**2*
     &    taugs**(-1)*cf*amm1 - 1.D0/2.D0*taucs**2*taugs**(-1)*cf*amm3
     &     + 2.D0*taucg**(-1)*taugs*mQ**2*cf*amm3*ffcg3 - 2.D0*
     &    taucg**(-1)*mQ**2*cf*amm2*ffcg3 - taucg**(-1)*mQ**2*xn*amm2*
     &    fL6m3 - 1.D0/2.D0*taucg*taugs**(-1)*tcs*xn**(-1)*amm3*fL6m2
     &     - 1.D0/2.D0*taucg*taugs**(-1)*xn**(-1)*amm2*fL6m2 - taucg*
     &    qsq**(-1)*mQ**2*cf*amm3*ffcg3 + 1.D0/2.D0*taugs**(-1)*qsq*
     &    xn**(-1)*amm2*fL6m2 + 1.D0/2.D0*taugs**(-1)*mQ**2*xn**(-1)*
     &    qcs*amm3*fL6m2 + taugs**(-1)*cf*cgs1*amm2 - taugs**(-1)*
     &    xn**(-1)*cgs1*amm2*fL6m2 + 1.D0/2.D0*taugs**(-1)*xn*qcg*amm2*
     &    fL6m3 + taugs*tcs*xn**(-1)*amm1*ffcs1 + taugs*qsq**(-1)*mQ**2
     &    *cf*amm3*ffcg1 - 2.D0*taugs*qsq**(-1)*mQ**2*cf*amm3*ffcg3 +
     &    taugs*xn**(-1)*amm1*ffcs2 + tcg**(-1)*mQ**4*cf*amm3 -
     &    qsq**(-1)*mQ**2*cf*amm2*ffcg3
      smm = smm - mQ**2*cf*amm1 - 2.D0*mQ**2*cf*amm3*ffcg3 - 3.D0*mQ**2
     &    *cf*amm3 - mQ**2*xn**(-1)*amm1*ffcs2 - mQ**2*xn**(-1)*amm3*
     &    ffcg3 - 1.D0/2.D0*mQ**2*xn*qcs*cgs1**(-1)*amm3*fL6m3 - 1.D0/2.
     &    D0*mQ**2*xn*qcg*cgs1**(-1)*amm1*fL6m3 - 1.D0/2.D0*mQ**2*xn*
     &    cgs1**(-1)*amm2*fL6m3 + 1.D0/2.D0*mQ**2*xn*amm1*fL6m3 - mQ**2
     &    *xn*amm3*fL6m3 - xlog*b0*amm2 - xlog*b0*amm4 + 3.D0/2.D0*epin
     &    *cf*amm2 + 3.D0/2.D0*epin*cf*amm4 + epin*b0*amm2 + epin*b0*
     &    amm4 - 1.D0/2.D0*epin*xn**(-1)*amm2 - 1.D0/2.D0*epin*xn**(-1)
     &    *amm4 + 1.D0/2.D0*epin*xn*amm2 + 1.D0/2.D0*epin*xn*amm4 - 1.D0
     &    /2.D0*epin2*xn**(-1)*amm2 - 1.D0/2.D0*epin2*xn**(-1)*amm4 + 3.
     &    D0/2.D0*epin2*xn*amm2 + 3.D0/2.D0*epin2*xn*amm4 + 1.D0/2.D0*
     &    cf*qcs*amm5*fq2m2 + 5.D0/2.D0*cf*amm2 + 5.D0/2.D0*cf*amm4 - 1.
     &    D0/2.D0*xn**(-1)*pisqo6*amm2 - 1.D0/2.D0*xn**(-1)*pisqo6*amm4
     &     - xn**(-1)*amm2*ffcs2 - xn**(-1)*amm2*fL6m2 - xn**(-1)*amm2
     &     + 1.D0/2.D0*xn*pisqo6*amm2 + 1.D0/2.D0*xn*pisqo6*amm4 -
     &    ddilog(
     & msq**(-1)*tcs)*xn**(-1)*amm2 - ddilog(msq**(-1)*tcs)*xn**(-1)*
     &    amm4 + ddilog(msq**(-1)*tcg)*xn*amm2 + ddilog(msq**(-1)*tcg)*
     &    xn*amm4 - lnrat( - taucs,msq)*taucs*tcs**(-1)*mQ**2*xn**(-1)*
     &    amm1 - 2.D0*lnrat( - taucs,msq)*taucs**2*tcs**(-1)*xn**(-1)*
     &    amm1 + lnrat( - taucs,msq)*tcs**(-1)*mQ**2*xn**(-1)*amm2 -
     &    lnrat( - taucs,msq)*tcs*xn**(-1)*amm3 + lnrat( - taucs,msq)*
     &    epin*xn**(-1)*amm2 + lnrat( - taucs,msq)*epin*xn**(-1)*amm4
     &     - lnrat( - taucs,msq)**2*xn**(-1)*amm2 - lnrat( - taucs,msq)
     &    **2*xn**(-1)*amm4 - 2.D0*lnrat( - taucg,msq)*taucg*tcg**(-1)*
     &    mQ**2*cf*amm1 - lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ**2*
     &    xn**(-1)*amm1 + lnrat( - taucg,msq)*taucg**2*tcg**(-2)*mQ**2*
     &    cf*amm3 - 2.D0*lnrat( - taucg,msq)*tcg**(-1)*mQ**2*cf*qcs*
     &    amm3 - 2.D0*lnrat( - taucg,msq)*tcg**(-1)*mQ**4*cf*amm1 -
     &    lnrat( - taucg,msq)*tcg**(-1)*mQ**4*xn**(-1)*amm1 - lnrat( -
     &    taucg,msq)*epin*xn*amm2 - lnrat( - taucg,msq)*epin*xn*amm4 +
     &    lnrat(
     &  - taucg,msq)**2*xn*amm2 + lnrat( - taucg,msq)**2*xn*amm4 -
     &    lnrat( - taugs,msq)*taucs*taugs**(-1)*cf*amm2 + lnrat( -
     &    taugs,msq)*taucs*cf*amm5 + lnrat( - taugs,msq)*taucs*xn**(-1)
     &    *amm3 + lnrat( - taugs,msq)*taucs**2*taugs**(-1)*cf*amm3 + 1.D
     &    0/2.D0*lnrat( - taugs,msq)*taucs**2*mQ**(-2)*cf*amm5 + 1.D0/2.
     &    D0*lnrat( - taugs,msq)*taugs*cf*amm5 - lnrat( - taugs,msq)*
     &    mQ**2*cf*amm3 - lnrat( - taugs,msq)*mQ**2*cf*amm5 - lnrat( -
     &    taugs,msq)*epin*xn*amm2 - lnrat( - taugs,msq)*epin*xn*amm4 +
     &    1.D0/2.D0*lnrat( - taugs,msq)*xn**(-1)*qcg*amm5 - lnrat( -
     &    taugs,msq)*xn*amm2 + 1.D0/2.D0*lnrat( - taugs,msq)**2*xn*amm2
     &     + 1.D0/2.D0*lnrat( - taugs,msq)**2*xn*amm4 - 3.D0*lnrat( -
     &    qsqhat,msq)*taucs*taucg*taugs**(-1)*qsq**(-1)*qsqhat*cf*amm3
     &     - 2.D0*lnrat( - qsqhat,msq)*taucs*taugs**(-1)*qsq**(-1)*
     &    mQ**2*qsqhat*cf*amm3 + 3.D0*lnrat( - qsqhat,msq)*taucs*
     &    qsq**(-1)*qsqhat*cf*amm3 + 2.D0*lnrat( - qsqhat,msq)*taucs*
     &    qsq**(-1)*qsqhat*cf*amm5
      smm = smm + 3.D0/2.D0*lnrat( - qsqhat,msq)*taucs*qsq**(-1)*qsqhat
     & *xn*amm1 - lnrat( - qsqhat,msq)*taucs**2*taugs**(-1)*qsq**(-1)*
     &    qsqhat*cf*amm3 - 1.D0/2.D0*lnrat( - qsqhat,msq)*taucs**2*
     &    taugs**(-1)*qsq**(-1)*qsqhat*cf*amm5 - 1.D0/2.D0*lnrat( -
     &    qsqhat,msq)*taucs**2*qsq**(-1)*mQ**(-2)*qsqhat*cf*amm5 +
     &    lnrat( - qsqhat,msq)*taucg**(-1)*taugs*qsq**(-1)*mQ**2*qsqhat
     &    *cf*amm3 - lnrat( - qsqhat,msq)*taugs**(-1)*qsq**(-1)*qsqhat*
     &    cf*qcg*amm2 - lnrat( - qsqhat,msq)*taugs**(-1)*qsq**(-1)*
     &    qsqhat*cf*cgs1*amm2 - lnrat( - qsqhat,msq)*taugs*qsq**(-1)*
     &    qsqhat*cf*amm3 - 3.D0/2.D0*lnrat( - qsqhat,msq)*taugs*
     &    qsq**(-1)*qsqhat*cf*amm5 - 1.D0/2.D0*lnrat( - qsqhat,msq)*
     &    taugs*qsq**(-1)*qsqhat*xn*amm1 + 2.D0*lnrat( - qsqhat,msq)*
     &    qsq**(-1)*mQ**2*qsqhat*cf*amm3 - 3.D0*lnrat( - qsqhat,msq)*
     &    qsq**(-1)*mQ**2*qsqhat*cf*amm5 + 3.D0/2.D0*lnrat( - qsqhat,
     &    msq)*qsq**(-1)*mQ**2*qsqhat*xn**(-1)*amm1 - 1.D0/2.D0*lnrat(
     &     - qsqhat,msq)*qsq**(-1)*mQ**2*qsqhat*xn*amm1
      smm = smm + lnrat( - qsqhat,msq)*qsq**(-1)*qsqhat*xn**(-1)*qcg*
     & amm1 + 1.D0/2.D0*lnrat( - qsqhat,msq)*qsq**(-1)*qsqhat*xn**(-1)*
     &    qcg*amm3 + lnrat( - qsqhat,msq)*qsq**(-1)*qsqhat*xn**(-1)*
     &    amm2 - 2.D0*C0fb2m(tcg,msq)*taucs**(-1)*taucg**(-1)*taugs*
     &    mQ**4*xn**(-1)*amm2 + C0fb2m(tcg,msq)*taucs**(-1)*taucg*mQ**2
     &    *xn**(-1)*cgs1*amm3 + C0fb2m(tcg,msq)*taucs**(-1)*mQ**2*
     &    xn**(-1)*qcg*amm2 - C0fb2m(tcg,msq)*taucs**(-1)*mQ**4*
     &    xn**(-1)*amm2 + 2.D0*C0fa2m(tcs,qsq,msq)*taucs**(-1)*
     &    taucg**(-2)*mQ**4*xn**(-1)*qcs**2*amm2 - C0fa2m(tcs,qsq,msq)*
     &    taucs**(-1)*taucg**(-1)*mQ**2*xn**(-1)*qcs**2*amm2 - C0fa2m(
     &    tcs,qsq,msq)*taucs**(-1)*taucg**(-1)*mQ**4*xn**(-1)*qcs*amm2
     &     - C0fa2m(tcs,qsq,msq)*taucs**(-1)*mQ**2*xn**(-1)*qcs*cgs1*
     &    amm3 + C0fa2m(tcs,qsq,msq)*taucs**(-1)*mQ**2*xn**(-1)*qcs*
     &    amm2 - 2.D0*C0fa2m(tcs,qsq,msq)*taucs*mQ**2*xn**(-1)*amm1 +
     &    C0fa2m(tcs,qsq,msq)*taucg**(-1)*taugs*mQ**4*xn**(-1)*amm1 - 2.
     &    D0*C0fa2m(tcs,qsq,msq)*taucg**(-1)*mQ**2*xn**(-1)*qcs*amm2
      smm = smm - C0fa2m(tcs,qsq,msq)*mQ**2*xn**(-1)*amm2 - 18.D0*
     &    ffc002*taucs*taugs**(-1)*cf*qcs*amm5 + 6.D0*ffc002*taucs**2*
     &    taugs**(-1)*mQ**(-2)*cf*qgs*amm5 + 6.D0*ffc002*taucg*
     &    taugs**(-1)*cf*qcg*amm5 + 12.D0*ffc002*mQ**2*cf*amm5 + 24.D0*
     &    ffc00*taucs*taucg*taugs**(-1)*cf*amm3 + 4.D0*ffc00*taucs*
     &    taugs**(-2)*cf*qgs*amm2 + 4.D0*ffc00*taucs*taugs**(-1)*tcs*cf
     &    *amm1 + 16.D0*ffc00*taucs*taugs**(-1)*mQ**2*cf*amm3 + 12.D0*
     &    ffc00*taucs*taugs**(-1)*cf*qcs*amm1 + 8.D0*ffc00*taucs*
     &    taugs**(-1)*cf*amm2 + 12.D0*ffc00*taucs*cf*amm3 - 4.D0*ffc00*
     &    taucs**2*taugs**(-2)*qsqhat*cf*amm3 + 8.D0*ffc00*taucs**2*
     &    taugs**(-1)*cf*amm3 + 12.D0*ffc00*taucg*taugs**(-1)*mQ**2*cf*
     &    amm1 - 8.D0*ffc00*taucg*taugs**(-1)*mQ**2*cf*amm3 - 4.D0*
     &    ffc00*taucg*taugs**(-1)*cf*amm2 - 16.D0*ffc00*taugs**(-1)*
     &    mQ**2*cf*qcs*amm1 - 8.D0*ffc00*taugs**(-1)*mQ**2*cf*amm2 + 2.D
     &    0*ffc00*taugs**(-1)*xn**(-1)*qgs*qcg*amm5 - 4.D0*ffc00*tcg*
     &    xn**(-1)*amm5
      smm = smm - 28.D0*ffc00*mQ**2*cf*amm3 + fI3me*taucs*mQ**2*
     &    xn**(-1)*amm3 - 3.D0*fI3me*taucs*mQ**2*xn*amm3 - fI3me*taucs*
     &    xn*amm2 - fI3me*taucg**(-1)*taugs*mQ**4*xn*amm1 + 2.D0*fI3me*
     &    taucg**(-1)*taugs*mQ**4*xn*amm3 - 2.D0*fI3me*taucg**(-1)*
     &    mQ**2*xn*cgs1*amm2 - fI3me*taucg*mQ**2*xn**(-1)*amm3 - fI3me*
     &    taugs*mQ**2*xn**(-1)*amm1 - 2.D0*fI3me*taugs*mQ**2*xn**(-1)*
     &    amm3 + 2.D0*fI3me*mQ**2*cf*amm2 - fI3me*mQ**2*xn*amm2 + fI3me
     &    *xn*cgs1*amm2


      lomm =-(za(ig,ie)*zb(ig,is)+za(ie,ic)*zb(is,ic))*zb(is,in)
     & /zb(ig,is)/zb(ig,ic)

c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      smm=smm + lomm * xn/six

c--- include finite counterterm due to FDH scheme
c--- gw = gw * ( 1 - 2 * cf * alphas(MSbar) / (2*pi))
      smm=smm - lomm * cf * two

      if (scheme == 'tH-V') then
        smm = smm - lomm * (cf/two + xn/six)
      endif

      virt_mm=smm
      return
      end

      function virt_pp(mQ,ig,is,ie,in,ic,p,musq)
      implicit none
      include 'types.f'
      complex(dp):: virt_pp

c***********************************************************************
c     Author: F. Tramontano                                            *
c     December, 2008                                                   *
c***********************************************************************
c---- One-loop ++ helicity amplitude for W+c production
c---- hg=+1  with s(mu) as gauge vector
c---- hQ=+1  with Q(mu)=Q1(mu)+mQ^2/2/Q.g*g(mu)
c for nwz=+1
c     f(-p1)+f(-p2)--> W^+(n(p3)+e^+(p4)) + Qbar(p5)
c For nwz=-1
c     f(-p1)+f(-p2)--> W^-(e^-(p3)+nbar(p4))+ Q(p5)
c----
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'epinv.f'
      include 'epinv2.f'
      include 'b0.f'
      include 'scheme.f'
      include 'zprods_decl.f'
      integer:: is,ig,ic,ie,in,nu,jpart
      real(dp):: p(mxpart,4),q(mxpart,4)
      real(dp):: dot,taucg,taugs,taucs,msq,mQ
      real(dp):: qsq,tcg,tcs,qsqhat,qcg,qcs,cgs1
      real(dp):: ddilog,epin,epin2,xlog
      complex(dp):: spp,lopp
      complex(dp):: app1,app2,app3,app4
      complex(dp):: ffcg1,ffcg3,ffcs2,fun1,fun2
      complex(dp):: fL6m1,fL6m2,fL6m3,L6m1,L6m2,L6m3
      complex(dp):: lnrat,C0fa2m,C0fb2m,I3me
      real(dp), intent(in) :: musq
c      scheme='dred'
      msq=mQ**2
      xlog=log(musq/msq)
      epin=epinv+xlog
      epin2=epinv2*epinv+epinv*xlog+half*xlog**2
      taugs=+two*dot(p,is,ig)
      taucs=+two*dot(p,ic,is)
      taucg=+two*dot(p,ic,ig)
      qsqhat=taucs+taucg+taugs
      qsq=qsqhat+msq
      tcg=taucg+msq
      tcs=taucs+msq
      qcg=qsqhat-taucg
      qcs=qsqhat-taucs
      cgs1=taucs-msq*taugs/taucg
      ffcg1=fun1(tcg,qsq,msq)
      ffcg3=fun2(qsq,tcg,msq)
      ffcs2=fun2(tcs,qsq,msq)
      fL6m1=L6m1(taucg,taucs,taugs,msq,qsq)
      fL6m2=L6m2(taucg,taucs,taugs,msq,qsq)
      fL6m3=L6m3(taucg,taucs,taugs,msq,qsq)
      do nu=1,4
      do jpart=1,5
      q(jpart,nu)=p(jpart,nu)
      if (jpart==ic) q(ic,nu)=p(ic,nu)-msq/taucg*p(ig,nu)
      enddo
      enddo
      call spinoru(5,q,za,zb)
      app1=za(ig,ie)*za(is,ic)*zb(is,in)/za(ig,is)/za(ig,ic)**2
      app2=za(ig,ie)*zb(ig,in)/za(ig,is)/za(ig,ic)/za(is,ic)/
     &     zb(is,ic)
      app3=za(ig,ie)*zb(is,in)/za(ig,is)/za(ig,ic)**2/zb(is,ic)
      app4=za(is,ie)*zb(is,in)/za(ig,is)/za(ig,ic)/za(is,ic)/
     &     zb(is,ic)
      spp=  + taucs*taucg*taugs**(-1)*mQ*xn**(-1)*app3*fL6m2 + 1.D0/2.D0
     &    *taucs*mQ*xn**(-1)*app3*fL6m1 + 1.D0/2.D0*taucs*mQ*xn*app3*
     &    fL6m3 + 2.D0*taucg**(-2)*taugs**2*qsq*mQ*cf*app4*ffcg3 +
     &    taucg**(-1)*taugs*tcs*mQ*xn**(-1)*app2*ffcs2 + 2.D0*
     &    taucg**(-1)*taugs*qsq*mQ*cf*app3*ffcg3 - taucg**(-1)*taugs*
     &    qsq*mQ*xn**(-1)*app4*ffcg3 + taucg**(-1)*taugs**2*mQ*cf*app4*
     &    ffcg1 - 3.D0*taucg**(-1)*taugs**2*mQ*cf*app4*ffcg3 -
     &    taucg**(-1)*mQ*qsqhat*cf*cgs1*app2 - taucg**(-1)*mQ*cf*qcs**2
     &    *app4 + taucg**(-1)*mQ*cf*qcg*cgs1*app2 - 1.D0/2.D0*
     &    taucg**(-1)*mQ**3*xn**(-1)*qcs*app4*fL6m2 - 1.D0/2.D0*
     &    taucg**(-1)*mQ**3*xn**(-1)*qcg*app2*fL6m2 + taucg*tcg**(-1)*
     &    qsq*mQ*cf*app4 - taucg*tcg**(-1)*mQ*cf*qcg*app3 - taugs**(-1)
     &    *mQ*xn**(-1)*qcg*cgs1*app2*fL6m2 + taugs*mQ*cf*app3*ffcg3 +
     &    taugs*mQ*cf*app3 - 2.D0*taugs*mQ*cf*app4*ffcg3 - 1.D0/2.D0*
     &    tcs*mQ*xn**(-1)*app3*fL6m2 + xlog*mQ*b0*cgs1*app2 - xlog*mQ*
     &    b0*app1
      spp = spp - 3.D0/2.D0*mQ*epin*cf*cgs1*app2 + 3.D0/2.D0*mQ*epin*cf
     &    *app1 - mQ*epin*b0*cgs1*app2 + mQ*epin*b0*app1 + 1.D0/2.D0*mQ
     &    *epin*xn**(-1)*cgs1*app2 - 1.D0/2.D0*mQ*epin*xn**(-1)*app1 -
     &    1.D0/2.D0*mQ*epin*xn*cgs1*app2 + 1.D0/2.D0*mQ*epin*xn*app1 +
     &    1.D0/2.D0*mQ*epin2*xn**(-1)*cgs1*app2 - 1.D0/2.D0*mQ*epin2*
     &    xn**(-1)*app1 - 3.D0/2.D0*mQ*epin2*xn*cgs1*app2 + 3.D0/2.D0*
     &    mQ*epin2*xn*app1 - 3.D0/2.D0*mQ*cf*cgs1*app2 - mQ*cf*cgs1*
     &    app3 + 5.D0/2.D0*mQ*cf*app1 + 1.D0/2.D0*mQ*xn**(-1)*cgs1*
     &    pisqo6*app2 + 1.D0/2.D0*mQ*xn**(-1)*cgs1*app2*fL6m1 + 1.D0/2.D
     &    0*mQ*xn**(-1)*cgs1*app2*fL6m2 + mQ*xn**(-1)*cgs1*app2 + mQ*
     &    xn**(-1)*cgs1*app3*fL6m2 - 1.D0/2.D0*mQ*xn**(-1)*pisqo6*app1
     &     - 1.D0/2.D0*mQ*xn*cgs1*pisqo6*app2 - 1.D0/2.D0*mQ*xn*cgs1*
     &    app2*fL6m3 - mQ*xn*cgs1*app3*fL6m3 + 1.D0/2.D0*mQ*xn*pisqo6*
     &    app1 + ddilog(msq**(-1)*tcs)*mQ*xn**(-1)*cgs1*app2 - ddilog(
     &    msq**(-1)*tcs)*mQ*xn**(-1)*app1 - ddilog(msq**(-1)*tcg)*mQ*xn
     &    *cgs1*app2
      spp = spp + ddilog(msq**(-1)*tcg)*mQ*xn*app1 - I3me(msq,taugs,qsq
     &    )*taucs**(-1)*taugs*mQ**3*xn**(-1)*cgs1*app2 - I3me(msq,taugs
     &    ,qsq)*taucs*taucg**(-1)*taugs*qsq**(-1)*mQ**5*cf*app4 - 1.D0/
     &    2.D0*I3me(msq,taugs,qsq)*taucs*taucg**(-1)*taugs*qsq**(-1)*
     &    mQ**5*xn**(-1)*app4 + I3me(msq,taugs,qsq)*taucs*taucg**(-1)*
     &    taugs*mQ**3*cf*app4 + 1.D0/2.D0*I3me(msq,taugs,qsq)*taucs*
     &    taucg**(-1)*taugs*mQ**3*xn**(-1)*app4 - I3me(msq,taugs,qsq)*
     &    taucs*taucg**(-1)*taugs**2*qsq**(-1)*mQ**3*cf*app4 - 1.D0/2.D0
     &    *I3me(msq,taugs,qsq)*taucs*taucg**(-1)*taugs**2*qsq**(-1)*
     &    mQ**3*xn**(-1)*app4 + I3me(msq,taugs,qsq)*taucs*taucg*taugs*
     &    qsq**(-1)*mQ*cf*app4 + 1.D0/2.D0*I3me(msq,taugs,qsq)*taucs*
     &    taucg*taugs*qsq**(-1)*mQ*xn**(-1)*app4 + I3me(msq,taugs,qsq)*
     &    taucs*taucg*qsq**(-1)*mQ*xn*cgs1*app2 + I3me(msq,taugs,qsq)*
     &    taucs*taucg*qsq**(-1)*mQ**3*cf*app4 + 1.D0/2.D0*I3me(msq,
     &    taugs,qsq)*taucs*taucg*qsq**(-1)*mQ**3*xn**(-1)*app4 - I3me(
     &    msq,taugs,qsq)*taucs*taucg*mQ*cf*app4
      spp = spp - 1.D0/2.D0*I3me(msq,taugs,qsq)*taucs*taucg*mQ*xn**(-1)
     & *app4 + I3me(msq,taugs,qsq)*taucs*taucg**2*qsq**(-1)*mQ*cf*app4
     &     + 1.D0/2.D0*I3me(msq,taugs,qsq)*taucs*taucg**2*qsq**(-1)*mQ*
     &    xn**(-1)*app4 - 2.D0*I3me(msq,taugs,qsq)*taucs*taugs*
     &    qsq**(-1)*mQ**3*cf*app4 - I3me(msq,taugs,qsq)*taucs*taugs*
     &    qsq**(-1)*mQ**3*xn**(-1)*app4 - I3me(msq,taugs,qsq)*taucs**2*
     &    taucg**(-1)*taugs*qsq**(-1)*mQ**3*cf*app4 - 1.D0/2.D0*I3me(
     &    msq,taugs,qsq)*taucs**2*taucg**(-1)*taugs*qsq**(-1)*mQ**3*
     &    xn**(-1)*app4 + 2.D0*I3me(msq,taugs,qsq)*taucs**2*taucg*
     &    qsq**(-1)*mQ*cf*app4 + I3me(msq,taugs,qsq)*taucs**2*taucg*
     &    qsq**(-1)*mQ*xn**(-1)*app4 + I3me(msq,taugs,qsq)*taucs**2*
     &    taugs*qsq**(-1)*mQ*cf*app4 + 1.D0/2.D0*I3me(msq,taugs,qsq)*
     &    taucs**2*taugs*qsq**(-1)*mQ*xn**(-1)*app4 + I3me(msq,taugs,
     &    qsq)*taucs**2*qsq**(-1)*mQ**3*cf*app4 + 1.D0/2.D0*I3me(msq,
     &    taugs,qsq)*taucs**2*qsq**(-1)*mQ**3*xn**(-1)*app4 - I3me(msq,
     &    taugs,qsq)*taucs**2*mQ*cf*app4
      spp = spp - 1.D0/2.D0*I3me(msq,taugs,qsq)*taucs**2*mQ*xn**(-1)*
     & app4 + I3me(msq,taugs,qsq)*taucs**3*qsq**(-1)*mQ*cf*app4 + 1.D0/
     &    2.D0*I3me(msq,taugs,qsq)*taucs**3*qsq**(-1)*mQ*xn**(-1)*app4
     &     + 4.D0*I3me(msq,taugs,qsq)*taucg**(-1)*taugs*mQ**3*cf*cgs1*
     &    app3 + 2.D0*I3me(msq,taugs,qsq)*taucg**(-1)*taugs*mQ**3*
     &    xn**(-1)*cgs1*app3 + I3me(msq,taugs,qsq)*taucg**(-1)*taugs*
     &    mQ**3*xn*cgs1*app2 + I3me(msq,taugs,qsq)*taucg*taugs*
     &    qsq**(-1)*mQ*xn*cgs1*app2 - I3me(msq,taugs,qsq)*taucg*taugs*
     &    qsq**(-1)*mQ**3*cf*app4 - 1.D0/2.D0*I3me(msq,taugs,qsq)*taucg
     &    *taugs*qsq**(-1)*mQ**3*xn**(-1)*app4 - I3me(msq,taugs,qsq)*
     &    taucg*qsq**(-1)*mQ*qsqhat*xn*cgs1*app2 + I3me(msq,taugs,qsq)*
     &    taucg**2*qsq**(-1)*mQ*xn*cgs1*app2 - 1.D0/2.D0*I3me(msq,taugs
     &    ,qsq)*taugs*qsq**(-1)*mQ*qsqhat*xn*cgs1*app2 - 1.D0/2.D0*
     &    I3me(msq,taugs,qsq)*taugs*qsq**(-1)*mQ**3*xn*cgs1*app2 -
     &    I3me(msq,taugs,qsq)*taugs*qsq**(-1)*mQ**5*cf*app4 - 1.D0/2.D0
     &    *I3me(msq,taugs,qsq)*taugs*qsq**(-1)*mQ**5*xn**(-1)*app4
      spp = spp + 1.D0/2.D0*I3me(msq,taugs,qsq)*taugs*mQ*xn*cgs1*app2
     &     + I3me(msq,taugs,qsq)*taugs*mQ**3*cf*app4 + 1.D0/2.D0*I3me(
     &    msq,taugs,qsq)*taugs*mQ**3*xn**(-1)*app4 - I3me(msq,taugs,qsq
     &    )*taugs**2*qsq**(-1)*mQ**3*cf*app4 - 1.D0/2.D0*I3me(msq,taugs
     &    ,qsq)*taugs**2*qsq**(-1)*mQ**3*xn**(-1)*app4 + 1.D0/2.D0*
     &    I3me(msq,taugs,qsq)*qsq**(-1)*mQ*qsqhat**2*xn*cgs1*app2 + 1.D0
     &    /2.D0*I3me(msq,taugs,qsq)*qsq**(-1)*mQ**3*qsqhat*xn*cgs1*app2
     &     - 1.D0/2.D0*I3me(msq,taugs,qsq)*mQ*qsqhat*xn*cgs1*app2 -
     &    lnrat( - taucs,msq)*taucs**2*taucg**(-1)*mQ*xn**(-1)*app2 -
     &    lnrat( - taucs,msq)*taucg*tcs**(-1)*mQ*xn**(-1)*cgs1*app3 +
     &    lnrat( - taucs,msq)*taugs*tcs**(-1)*mQ*xn**(-1)*cgs1*app2 +
     &    lnrat( - taucs,msq)*tcs**(-1)*qsq*mQ*xn**(-1)*cgs1*app4 -
     &    lnrat( - taucs,msq)*mQ*epin*xn**(-1)*cgs1*app2 + lnrat( -
     &    taucs,msq)*mQ*epin*xn**(-1)*app1 + lnrat( - taucs,msq)**2*mQ*
     &    xn**(-1)*cgs1*app2 - lnrat( - taucs,msq)**2*mQ*xn**(-1)*app1
     &     - 2.D0*lnrat(
     &  - taucg,msq)*taucg**(-1)*tcg**(-1)*qsq*mQ*cf*qcs**2*app4 - 2.D0
     &    *lnrat( - taucg,msq)*taucg**(-1)*mQ*cf*qcs**2*app4 + 1.D0/2.D0
     &    *lnrat( - taucg,msq)*taucg*tcg**(-1)*qsq*mQ*xn**(-1)*app2 - 2.
     &    D0*lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ*cf*qcs*app3 - 1.D0/
     &    2.D0*lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ*xn**(-1)*cgs1*
     &    app2 - lnrat( - taucg,msq)*taucg*tcg**(-1)*mQ*xn**(-1)*cgs1*
     &    app3 + lnrat( - taucg,msq)*taucg**2*tcg**(-2)*qsq*mQ*cf*app3
     &     - lnrat( - taucg,msq)*taucg**2*tcg**(-2)*qsq*mQ*cf*app4 - 2.D
     &    0*lnrat( - taucg,msq)*tcg**(-1)*qsq*mQ*cf*qcs*app3 + 4.D0*
     &    lnrat( - taucg,msq)*tcg**(-1)*qsq*mQ*cf*qcs*app4 + lnrat( -
     &    taucg,msq)*tcg**(-1)*qsq*mQ*xn**(-1)*qcs*app4 + 2.D0*lnrat(
     &     - taucg,msq)*tcg**(-1)*mQ*cf*qcs**2*app4 + lnrat( - taucg,
     &    msq)*mQ*epin*xn*cgs1*app2 - lnrat( - taucg,msq)*mQ*epin*xn*
     &    app1 + 2.D0*lnrat( - taucg,msq)*mQ*cf*cgs1*app2 - 1.D0/2.D0*
     &    lnrat( - taucg,msq)*mQ*xn**(-1)*qcs*app2 - lnrat( - taucg,msq
     &    )*mQ*xn**(-1)*cgs1*app2
      spp = spp - lnrat( - taucg,msq)**2*mQ*xn*cgs1*app2 + lnrat( -
     &    taucg,msq)**2*mQ*xn*app1 + lnrat( - taugs,msq)*mQ*epin*xn*
     &    cgs1*app2 - lnrat( - taugs,msq)*mQ*epin*xn*app1 - 1.D0/2.D0*
     &    lnrat( - taugs,msq)**2*mQ*xn*cgs1*app2 + 1.D0/2.D0*lnrat( -
     &    taugs,msq)**2*mQ*xn*app1 + lnrat( - qsqhat,msq)*taucs*
     &    taucg**(-1)*tcs*qsq**(-1)*mQ*qsqhat*xn**(-1)*app2 + lnrat( -
     &    qsqhat,msq)*taucs*taucg**(-1)*mQ*cf*qcg*app4 + 1.D0/2.D0*
     &    lnrat( - qsqhat,msq)*taucs*taucg**(-1)*mQ*xn**(-1)*qcg*app4
     &     - 1.D0/2.D0*lnrat( - qsqhat,msq)*taucs*qsq**(-1)*mQ*qsqhat*
     &    xn**(-1)*app2 + lnrat( - qsqhat,msq)*taucs*mQ*cf*app4 + 1.D0/
     &    2.D0*lnrat( - qsqhat,msq)*taucs*mQ*xn**(-1)*app4 + lnrat( -
     &    qsqhat,msq)*taucg**(-2)*qsq**(-1)*mQ**3*cf*qcs**2*qcg*app4 +
     &    lnrat( - qsqhat,msq)*taucg**(-2)*mQ*cf*qcs**2*qcg*app4 + 1.D0/
     &    2.D0*lnrat( - qsqhat,msq)*taucg**(-1)*taugs*qsq**(-1)*mQ**3*
     &    qsqhat*xn**(-1)*app2 + lnrat( - qsqhat,msq)*taucg**(-1)*
     &    qsq**(-1)*mQ*qsqhat*cf*qcs*qcg*app3
      spp = spp - lnrat( - qsqhat,msq)*taucg**(-1)*qsq**(-1)*mQ**3*
     & qsqhat*cf*qcs*app4 - 1.D0/2.D0*lnrat( - qsqhat,msq)*taucg**(-1)*
     &    qsq**(-1)*mQ**3*qsqhat*xn**(-1)*qcs*app4 + lnrat( - qsqhat,
     &    msq)*taucg**(-1)*qsq**(-1)*mQ**3*cf*qcs*qcg*app4 + lnrat( -
     &    qsqhat,msq)*taucg**(-1)*qsq**(-1)*mQ**3*cf*qcs**2*app4 + 1.D0/
     &    2.D0*lnrat( - qsqhat,msq)*taucg**(-1)*qsq**(-1)*mQ**3*
     &    xn**(-1)*qcs*qcg*app4 + lnrat( - qsqhat,msq)*taucg**(-1)*mQ*
     &    qsqhat*cf*qcs*app3 - lnrat( - qsqhat,msq)*taucg**(-1)*mQ*
     &    qsqhat*cf*qcs*app4 - 1.D0/2.D0*lnrat( - qsqhat,msq)*
     &    taucg**(-1)*mQ*qsqhat*xn**(-1)*qcs*app4 + lnrat( - qsqhat,msq
     &    )*taucg**(-1)*mQ*cf*qcs**2*app4 - lnrat( - qsqhat,msq)*
     &    taucg**(-1)*mQ*cf*qcg**2*app4 - 1.D0/2.D0*lnrat( - qsqhat,msq
     &    )*taucg**(-1)*mQ*xn**(-1)*qcg**2*app4 - lnrat( - qsqhat,msq)*
     &    taucg*mQ*cf*app4 - 1.D0/2.D0*lnrat( - qsqhat,msq)*taucg*mQ*
     &    xn**(-1)*app4 + lnrat( - qsqhat,msq)*qsq**(-1)*mQ*qsqhat*cf*
     &    qcs*app3
      spp = spp + 3.D0/2.D0*lnrat( - qsqhat,msq)*qsq**(-1)*mQ*qsqhat*
     & xn**(-1)*cgs1*app2 + lnrat( - qsqhat,msq)*qsq**(-1)*mQ**3*cf*qcs
     &    *app4 + 1.D0/2.D0*lnrat( - qsqhat,msq)*qsq**(-1)*mQ**3*
     &    xn**(-1)*qcs*app4 + lnrat( - qsqhat,msq)*mQ*qsqhat*cf*app3 -
     &    2.D0*lnrat( - qsqhat,msq)*mQ*cf*qcg*app4 - lnrat( - qsqhat,
     &    msq)*mQ*xn**(-1)*qcg*app4 - C0fb2m(tcg,msq)*taucs**(-1)*taugs
     &    *mQ**3*xn**(-1)*cgs1*app3 - C0fb2m(tcg,msq)*taucs**(-1)*mQ**3
     &    *xn**(-1)*qcs*cgs1*app3 - C0fb2m(tcg,msq)*taucs**(-1)*mQ**3*
     &    xn**(-1)*qcs*cgs1*app4 + 2.D0*C0fb2m(tcg,msq)*mQ**3*xn**(-1)*
     &    cgs1*app2 + 2.D0*C0fa2m(tcs,qsq,msq)*taucs**(-1)*taucg**(-1)*
     &    taugs*mQ**3*xn**(-1)*qcs*cgs1*app3 + C0fa2m(tcs,qsq,msq)*
     &    taucs**(-1)*taucg**(-1)*mQ**3*xn**(-1)*qcs**2*cgs1*app4 +
     &    C0fa2m(tcs,qsq,msq)*taucs**(-1)*mQ**3*xn**(-1)*qcs*cgs1*app3
     &     - C0fa2m(tcs,qsq,msq)*taucg**(-1)*taugs*mQ**3*xn**(-1)*cgs1*
     &    app2 - 2.D0*C0fa2m(tcs,qsq,msq)*mQ**3*xn**(-1)*cgs1*app2


      lopp=-mQ*za(ig,ie)/za(ig,is)/za(ig,ic)**2*
     & (za(is,ic)*zb(is,in)+za(ig,ic)*zb(ig,in))

c--- include finite counterterm to go from DR to MSbar scheme
c--- alphas(DR) = alphas(MSbar) * (1+ (Nc / 6) * alphas(MSbar) / (2*pi))
      spp=spp + lopp * xn/six

c--- include finite counterterm due to FDH scheme
c--- gw = gw * ( 1 - 2*cf * alphas(MSbar) / (2*pi))
      spp=spp - lopp * cf * two

      if (scheme == 'tH-V') then
        spp = spp - lopp * (cf/two + xn/six)
      endif

      virt_pp=spp
      return
      end

      end module
