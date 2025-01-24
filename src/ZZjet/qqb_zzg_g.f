!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_zzg_g(p,msq)
c---  Author: R.K. Ellis
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  + n(p3)+ebar(p4)+e(p5)+nubar(p6)+ g(p7)+ g(p8)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'nf.f'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'ewcharge.f'
      include 'pchoice.f'
      include 'srdiags.f'
      include 'cplx.h'
      include 'interference.f'

      integer polg1,polg2,hq,nq,h34,h56,ii,nmax
      real(dp):: p(mxpart,4),msq(-nf:nf,-nf:nf),
     & fac,q34,q56,s3456,
     & ampsq_dd,ampsq_du,ampsq_ds,ampsq_ud,ampsq_uu,ampsq_uc,
     & ampsq_dd_uu,ampsq_dd_ss,ampsq_uu_dd,ampsq_uu_cc,
     & res,resu,resd
      complex(dp):: LR(2,2),v34(2),v56(2)

      complex(dp):: prop34,prop56,prop3456,
     & coupqe(2,2,2,2),coupsr3456(2,2,2,2),coupsr5634(2,2,2,2),
     & AB(2,2,2,2,2),BA(2,2,2,2,2),
     & q_aAB(2,2,2,2,2),a_qAB(2,2,2,2,2),
     & q_gAB(2,2,2,2,2),g_aAB(2,2,2,2,2),
     & a_gAB(2,2,2,2,2),g_qAB(2,2,2,2,2),
     & q_aBA(2,2,2,2,2),a_qBA(2,2,2,2,2),
     & q_gBA(2,2,2,2,2),g_aBA(2,2,2,2,2),
     & a_gBA(2,2,2,2,2),g_qBA(2,2,2,2,2),
     & g_gAB(2,2,2,2,2),g_gBA(2,2,2,2,2),
     & q_a34AB(2,2,2,2,2),q_a34BA(2,2,2,2,2),
     & q_a56AB(2,2,2,2,2),q_a56BA(2,2,2,2,2),
     & a_q34AB(2,2,2,2,2),a_q34BA(2,2,2,2,2),
     & a_q56AB(2,2,2,2,2),a_q56BA(2,2,2,2,2),
     & g_a34AB(2,2,2,2,2),g_a34BA(2,2,2,2,2),
     & g_a56AB(2,2,2,2,2),g_a56BA(2,2,2,2,2),
     & a_g34AB(2,2,2,2,2),a_g34BA(2,2,2,2,2),
     & a_g56AB(2,2,2,2,2),a_g56BA(2,2,2,2,2),
     & q_g34AB(2,2,2,2,2),q_g34BA(2,2,2,2,2),
     & q_g56AB(2,2,2,2,2),q_g56BA(2,2,2,2,2),
     & g_q34AB(2,2,2,2,2),g_q34BA(2,2,2,2,2),
     & g_q56AB(2,2,2,2,2),g_q56BA(2,2,2,2,2),
     & g_g34AB(2,2,2,2,2),g_g34BA(2,2,2,2,2),
     & g_g56AB(2,2,2,2,2),g_g56BA(2,2,2,2,2)

      integer,parameter::i4(2)=(/4,6/),i6(2)=(/6,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      msq(:,:)=0d0

      fac=V*xn*(4d0*esq**2*gsq)**2
      call spinoru(8,p,za,zb)
c      write(6,*) 'zmass ',zmass
c      write(6,*) 'zwidth',zwidth

c     debug interference not implemented yet.
c      interference=.false.
c      srdiags=.true.
c      write(6,*) 'srdiags',srdiags

      if (interference) then
         nmax=2
      else
         nmax=1
      endif

      do ii=1,nmax

c--   calculate propagators
      s3456=s(3,4)+s(3,5)+s(3,6)+s(4,5)+s(4,6)+s(5,6)
      prop3456=s3456/cplx2(s3456-zmass**2,zmass*zwidth)
      prop34=s(3,i4(ii))/cplx2(s(3,i4(ii))-zmass**2,zmass*zwidth)
      prop56=s(5,i6(ii))/cplx2(s(5,i6(ii))-zmass**2,zmass*zwidth)

c---remember u-ub is the basic process.
      call zzggamps(1,2,3,4,5,6,7,8,q_aAB,q_aBA) !quark-antiquark
      call zzggamps(2,1,3,4,5,6,7,8,a_qAB,a_qBA) !quark-antiquark
      call zzggamps(1,7,3,4,5,6,2,8,q_gAB,q_gBA) !quark-gluon
      call zzggamps(2,7,3,4,5,6,1,8,g_qAB,g_qBA) !gluon-quark
      call zzggamps(7,1,3,4,5,6,2,8,a_gAB,a_gBA) !antiquark-gluon
      call zzggamps(7,2,3,4,5,6,1,8,g_aAB,g_aBA) !gluon-antiquark
      call zzggamps(8,7,3,4,5,6,1,2,g_gAB,g_gBA) !gluon-gluon
      if (srdiags) then
      call a8sr(1,2,3,4,5,6,7,8,q_a56AB,q_a56BA,q_a34AB,q_a34BA)
      call a8sr(2,1,3,4,5,6,7,8,a_q56AB,a_q56BA,a_q34AB,a_q34BA)
      call a8sr(1,7,3,4,5,6,2,8,q_g56AB,q_g56BA,q_g34AB,q_g34BA)
      call a8sr(2,7,3,4,5,6,1,8,g_q56AB,g_q56BA,g_q34AB,g_q34BA)
      call a8sr(7,1,3,4,5,6,2,8,a_g56AB,a_g56BA,a_g34AB,a_g34BA)
      call a8sr(7,2,3,4,5,6,1,8,g_a56AB,g_a56BA,g_a34AB,g_a34BA)
      call a8sr(8,7,3,4,5,6,1,2,g_g56AB,g_g56BA,g_g34AB,g_g34BA)
      endif

      v34(1)=zl1
      v34(2)=zr1
      q34=q1
      v56(1)=zl2
      v56(2)=zr2
      q56=q2
c--- initialize couplings of Z/gamma to quarks and lepton pairs.
      do nq=1,2
      do hq=1,2
      LR(nq,hq)=(2-hq)*zL(nq)+(hq-1)*zR(nq)
      do h34=1,2
      do h56=1,2
      coupqe(nq,hq,h34,h56)=(Q(nq)*q34+LR(nq,hq)*v34(h34)*prop34)
     &                     *(Q(nq)*q56+LR(nq,hq)*v56(h56)*prop56)
      coupsr3456(nq,hq,h34,h56)=(Q(nq)*q34+LR(nq,hq)*v34(h34)*prop3456)
     &                     *(q34*q56+v34(h34)*v56(h56)*prop56)
      coupsr5634(nq,hq,h34,h56)=(Q(nq)*q56+LR(nq,hq)*v56(h56)*prop3456)
     &                     *(q34*q56+v34(h34)*v56(h56)*prop34)
      enddo
      enddo
      enddo
      enddo


c--- first populate the basic set of matrix elements
c--- (for the first generation of quarks)
      do j=-2,2
      do k=-2,2
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) cycle

      if ((j > 0) .and. (k < 0)) then
      do polg1=1,2
      do polg2=1,2
      do hq=1,2
      do h34=1,2
      do h56=1,2
      AB(polg1,polg2,hq,h34,h56)=
     & q_aAB(polg1,polg2,hq,h34,h56)*coupqe(j,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=
     & q_aBA(polg1,polg2,hq,h34,h56)*coupqe(j,hq,h34,h56)
      if (srdiags) then
      AB(polg1,polg2,hq,h34,h56)=AB(polg1,polg2,hq,h34,h56)
     & +q_a56AB(polg1,polg2,hq,h34,h56)*coupsr3456(j,hq,h34,h56)
     & +q_a34AB(polg1,polg2,hq,h34,h56)*coupsr5634(j,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=BA(polg1,polg2,hq,h34,h56)
     & +q_a56BA(polg2,polg1,hq,h34,h56)*coupsr3456(j,hq,h34,h56)
     & +q_a34BA(polg2,polg1,hq,h34,h56)*coupsr5634(j,hq,h34,h56)
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
      call msquarer(AB,BA,res)
c     half for identical gluons in final state
      msq(j,k)=half*fac*aveqq*res

      elseif ((j < 0) .and. (k > 0)) then
      do polg1=1,2
      do polg2=1,2
      do hq=1,2
      do h34=1,2
      do h56=1,2
      AB(polg1,polg2,hq,h34,h56)=
     & a_qAB(polg1,polg2,hq,h34,h56)*coupqe(k,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=
     & a_qBA(polg1,polg2,hq,h34,h56)*coupqe(k,hq,h34,h56)

      if (srdiags) then
      AB(polg1,polg2,hq,h34,h56)=AB(polg1,polg2,hq,h34,h56)
     & +a_q56AB(polg1,polg2,hq,h34,h56)*coupsr3456(k,hq,h34,h56)
     & +a_q34AB(polg1,polg2,hq,h34,h56)*coupsr5634(k,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=BA(polg1,polg2,hq,h34,h56)
     & +a_q56BA(polg2,polg1,hq,h34,h56)*coupsr3456(k,hq,h34,h56)
     & +a_q34BA(polg2,polg1,hq,h34,h56)*coupsr5634(k,hq,h34,h56)
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
      call msquarer(AB,BA,res)
c     half for identical gluons in final state
      msq(j,k)=half*fac*aveqq*res
c      write(6,*) 'j<0,k>0,res',j,k,msq(j,k)

      elseif ((j > 0) .and. (k == 0)) then
      do polg1=1,2
      do polg2=1,2
      do hq=1,2
      do h34=1,2
      do h56=1,2
      AB(polg1,polg2,hq,h34,h56)=
     & q_gAB(polg1,polg2,hq,h34,h56)*coupqe(j,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=
     & q_gBA(polg1,polg2,hq,h34,h56)*coupqe(j,hq,h34,h56)
      if (srdiags) then
      AB(polg1,polg2,hq,h34,h56)=AB(polg1,polg2,hq,h34,h56)
     & +q_g56AB(polg1,polg2,hq,h34,h56)*coupsr3456(j,hq,h34,h56)
     & +q_g34AB(polg1,polg2,hq,h34,h56)*coupsr5634(j,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=BA(polg1,polg2,hq,h34,h56)
     & +q_g56BA(polg2,polg1,hq,h34,h56)*coupsr3456(j,hq,h34,h56)
     & +q_g34BA(polg2,polg1,hq,h34,h56)*coupsr5634(j,hq,h34,h56)
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
      call msquarer(AB,BA,res)
      msq(j,k)=fac*aveqg*res

      elseif ((j == 0) .and. (k > 0)) then
      do polg1=1,2
      do polg2=1,2
      do hq=1,2
      do h34=1,2
      do h56=1,2
      AB(polg1,polg2,hq,h34,h56)=
     & g_qAB(polg1,polg2,hq,h34,h56)*coupqe(k,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=
     & g_qBA(polg1,polg2,hq,h34,h56)*coupqe(k,hq,h34,h56)
      if (srdiags) then
      AB(polg1,polg2,hq,h34,h56)=AB(polg1,polg2,hq,h34,h56)
     & +g_q56AB(polg1,polg2,hq,h34,h56)*coupsr3456(k,hq,h34,h56)
     & +g_q34AB(polg1,polg2,hq,h34,h56)*coupsr5634(k,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=BA(polg1,polg2,hq,h34,h56)
     & +g_q56BA(polg2,polg1,hq,h34,h56)*coupsr3456(k,hq,h34,h56)
     & +g_q34BA(polg2,polg1,hq,h34,h56)*coupsr5634(k,hq,h34,h56)
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
      call msquarer(AB,BA,res)
      msq(j,k)=fac*aveqg*res

      elseif ((j < 0) .and. (k == 0)) then
      do polg1=1,2
      do polg2=1,2
      do hq=1,2
      do h34=1,2
      do h56=1,2
      AB(polg1,polg2,hq,h34,h56)=
     & a_gAB(polg1,polg2,hq,h34,h56)*coupqe(-j,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=
     & a_gBA(polg1,polg2,hq,h34,h56)*coupqe(-j,hq,h34,h56)
      if (srdiags) then
      AB(polg1,polg2,hq,h34,h56)=AB(polg1,polg2,hq,h34,h56)
     & +a_g56AB(polg1,polg2,hq,h34,h56)*coupsr3456(-j,hq,h34,h56)
     & +a_g34AB(polg1,polg2,hq,h34,h56)*coupsr5634(-j,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=BA(polg1,polg2,hq,h34,h56)
     & +a_g56BA(polg2,polg1,hq,h34,h56)*coupsr3456(-j,hq,h34,h56)
     & +a_g34BA(polg2,polg1,hq,h34,h56)*coupsr5634(-j,hq,h34,h56)
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
      call msquarer(AB,BA,res)
      msq(j,k)=fac*aveqg*res

      elseif ((j == 0) .and. (k < 0)) then
      do polg1=1,2
      do polg2=1,2
      do hq=1,2
      do h34=1,2
      do h56=1,2
      AB(polg1,polg2,hq,h34,h56)=
     & g_aAB(polg1,polg2,hq,h34,h56)*coupqe(-k,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=
     & g_aBA(polg1,polg2,hq,h34,h56)*coupqe(-k,hq,h34,h56)
      if (srdiags) then
      AB(polg1,polg2,hq,h34,h56)=AB(polg1,polg2,hq,h34,h56)
     & +g_a56AB(polg1,polg2,hq,h34,h56)*coupsr3456(-k,hq,h34,h56)
     & +g_a34AB(polg1,polg2,hq,h34,h56)*coupsr5634(-k,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=BA(polg1,polg2,hq,h34,h56)
     & +g_a56BA(polg2,polg1,hq,h34,h56)*coupsr3456(-k,hq,h34,h56)
     & +g_a34BA(polg2,polg1,hq,h34,h56)*coupsr5634(-k,hq,h34,h56)
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
      call msquarer(AB,BA,res)
      msq(j,k)=fac*aveqg*res

      elseif ((j == 0) .and. (k == 0)) then
      do polg1=1,2
      do polg2=1,2
      do hq=1,2
      do h34=1,2
      do h56=1,2
      AB(polg1,polg2,hq,h34,h56)=
     & g_gAB(polg1,polg2,hq,h34,h56)*coupqe(1,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=
     & g_gBA(polg1,polg2,hq,h34,h56)*coupqe(1,hq,h34,h56)
      if (srdiags) then
      AB(polg1,polg2,hq,h34,h56)=AB(polg1,polg2,hq,h34,h56)
     & +g_g56AB(polg1,polg2,hq,h34,h56)*coupsr3456(1,hq,h34,h56)
     & +g_g34AB(polg1,polg2,hq,h34,h56)*coupsr5634(1,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=BA(polg1,polg2,hq,h34,h56)
     & +g_g56BA(polg2,polg1,hq,h34,h56)*coupsr3456(1,hq,h34,h56)
     & +g_g34BA(polg2,polg1,hq,h34,h56)*coupsr5634(1,hq,h34,h56)
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
      call msquarer(AB,BA,resd)

      do polg1=1,2
      do polg2=1,2
      do hq=1,2
      do h34=1,2
      do h56=1,2
      AB(polg1,polg2,hq,h34,h56)=
     & g_gAB(polg1,polg2,hq,h34,h56)*coupqe(2,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=
     & g_gBA(polg1,polg2,hq,h34,h56)*coupqe(2,hq,h34,h56)
      if (srdiags) then
      AB(polg1,polg2,hq,h34,h56)=AB(polg1,polg2,hq,h34,h56)
     & +g_g56AB(polg1,polg2,hq,h34,h56)*coupsr3456(2,hq,h34,h56)
     & +g_g34AB(polg1,polg2,hq,h34,h56)*coupsr5634(2,hq,h34,h56)
      BA(polg1,polg2,hq,h34,h56)=BA(polg1,polg2,hq,h34,h56)
     & +g_g56BA(polg2,polg1,hq,h34,h56)*coupsr3456(2,hq,h34,h56)
     & +g_g34BA(polg2,polg1,hq,h34,h56)*coupsr5634(2,hq,h34,h56)
      endif

      enddo
      enddo
      enddo
      enddo
      enddo
      call msquarer(AB,BA,resu)
      msq(j,k)=fac*avegg*(2._dp*resu+3._dp*resd)
      endif

      enddo !endo of k loop
      enddo !endo of j loop

      enddo ! end of interference loop

c---extend to full flavour range
      do j=-nf,nf
      do k=-nf,nf
        if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) cycle
        msq(j,k)=msq(jkswitch(j),jkswitch(k))
      enddo
      enddo

c      return
c      msq(:,:)=0

c---*****************4 quark diagrams********************
      call a8ZZfourqsq(1,2,3,4,5,6,7,8,
     & ampsq_dd,ampsq_du,ampsq_ds,ampsq_ud,ampsq_uu,ampsq_uc,
     & ampsq_dd_uu,ampsq_dd_ss,ampsq_uu_dd,ampsq_uu_cc)
      msq(1,1)=half*ampsq_dd
      msq(1,2)=ampsq_du
      msq(1,3)=ampsq_ds
      msq(2,1)=ampsq_ud
      msq(2,2)=half*ampsq_uu
      msq(2,4)=ampsq_uc

      msq(1,4)=msq(1,2)
      msq(1,5)=msq(1,3)
      msq(2,3)=msq(2,1)
      msq(2,5)=msq(2,1)

      msq(3,1)=msq(1,3)
      msq(3,2)=msq(1,2)
c      msq(3,3)=msq(1,1)
      msq(3,4)=msq(1,2)
      msq(3,5)=msq(1,3)

      msq(4,1)=msq(2,1)
      msq(4,2)=msq(2,4)
      msq(4,3)=msq(2,1)
c      msq(4,4)=msq(2,2)
      msq(4,5)=msq(2,1)

      msq(5,1)=msq(1,3)
      msq(5,2)=msq(1,2)
      msq(5,3)=msq(1,3)
      msq(5,4)=msq(1,2)
c      msq(5,5)=msq(1,1)
      call a8ZZfourqsq(7,2,3,4,5,6,1,8,
     & ampsq_dd,ampsq_du,ampsq_ds,ampsq_ud,ampsq_uu,ampsq_uc,
     & ampsq_dd_uu,ampsq_dd_ss,ampsq_uu_dd,ampsq_uu_cc)
      msq(-1,1)=msq(-1,1)+ampsq_dd+2._dp*ampsq_dd_ss+2._dp*ampsq_uu_dd
      msq(-1,2)=ampsq_du
      msq(-1,3)=ampsq_ds
      msq(-2,1)=ampsq_ud
      msq(-2,2)=msq(-2,2)+ampsq_uu+3._dp*ampsq_dd_uu+ampsq_uu_cc
      msq(-2,4)=ampsq_uc

      msq(-1,4)=msq(-1,2)
      msq(-1,5)=msq(-1,3)
      msq(-2,3)=msq(-2,1)
      msq(-2,5)=msq(-2,1)

      msq(-3,1)=msq(-1,3)
      msq(-3,2)=msq(-1,2)
c      msq(-3,3)=msq(-1,1)
      msq(-3,4)=msq(-1,2)
      msq(-3,5)=msq(-1,3)

      msq(-4,1)=msq(-2,1)
      msq(-4,2)=msq(-2,4)
      msq(-4,3)=msq(-2,1)
c      msq(-4,4)=msq(-2,2)
      msq(-4,5)=msq(-2,1)

      msq(-5,1)=msq(-1,3)
      msq(-5,2)=msq(-1,2)
      msq(-5,3)=msq(-1,3)
      msq(-5,4)=msq(-1,2)
c      msq(-5,5)=msq(-1,1)
      call a8ZZfourqsq(1,8,3,4,5,6,7,2,
     & ampsq_dd,ampsq_du,ampsq_ds,ampsq_ud,ampsq_uu,ampsq_uc,
     & ampsq_dd_uu,ampsq_dd_ss,ampsq_uu_dd,ampsq_uu_cc)
      msq(1,-1)=msq(1,-1)+ampsq_dd+2._dp*ampsq_dd_ss+2._dp*ampsq_dd_uu
      msq(1,-2)=ampsq_du
      msq(1,-3)=ampsq_ds
      msq(2,-1)=ampsq_ud
      msq(2,-2)=msq(2,-2)+ampsq_uu+3._dp*ampsq_uu_dd+ampsq_uu_cc
      msq(2,-4)=ampsq_uc

      msq(1,-4)=msq(1,-2)
      msq(1,-5)=msq(1,-3)
      msq(2,-3)=msq(2,-1)
      msq(2,-5)=msq(2,-1)

      msq(3,-1)=msq(1,-3)
      msq(3,-2)=msq(1,-2)
c      msq(3,-3)=msq(1,-1)
      msq(3,-4)=msq(1,-2)
      msq(3,-5)=msq(1,-3)

      msq(4,-1)=msq(2,-1)
      msq(4,-2)=msq(2,-4)
      msq(4,-3)=msq(2,-1)
c      msq(4,-4)=msq(2,-2)
      msq(4,-5)=msq(2,-1)

      msq(5,-1)=msq(1,-3)
      msq(5,-2)=msq(1,-2)
      msq(5,-3)=msq(1,-3)
      msq(5,-4)=msq(1,-2)
c      msq(5,-5)=msq(1,-1)
      call a8ZZfourqsq(7,8,3,4,5,6,1,2,
     & ampsq_dd,ampsq_du,ampsq_ds,ampsq_ud,ampsq_uu,ampsq_uc,
     & ampsq_dd_uu,ampsq_dd_ss,ampsq_uu_dd,ampsq_uu_cc)
      msq(-1,-1)=half*ampsq_dd
      msq(-1,-2)=ampsq_du
      msq(-1,-3)=ampsq_ds
      msq(-2,-1)=ampsq_ud
      msq(-2,-2)=half*ampsq_uu
      msq(-2,-4)=ampsq_uc

      msq(-1,-4)=msq(-1,-2)
      msq(-1,-5)=msq(-1,-3)
      msq(-2,-3)=msq(-2,-1)
      msq(-2,-5)=msq(-2,-1)

      msq(-3,-1)=msq(-1,-3)
      msq(-3,-2)=msq(-1,-2)
c      msq(-3,-3)=msq(-1,-1)
      msq(-3,-4)=msq(-1,-2)
      msq(-3,-5)=msq(-1,-3)

      msq(-4,-1)=msq(-2,-1)
      msq(-4,-2)=msq(-2,-4)
      msq(-4,-3)=msq(-2,-1)
c      msq(-4,-4)=msq(-2,-2)
      msq(-4,-5)=msq(-2,-1)

      msq(-5,-1)=msq(-1,-3)
      msq(-5,-2)=msq(-1,-2)
      msq(-5,-3)=msq(-1,-3)
      msq(-5,-4)=msq(-1,-2)
c      msq(-5,-5)=msq(-1,-1)

c Extend to other flavours
      do j=3,5
      msq(-j,-j)=msq(-(j-2),-(j-2))
      msq(+j,+j)=msq(+(j-2),+(j-2))
      msq(+j,-j)=msq(+(j-2),-(j-2))
      msq(-j,+j)=msq(-(j-2),+(j-2))
      msq(j,0)=msq((j-2),0)
      msq(0,j)=msq(0,(j-2))
      msq(-j,0)=msq(-(j-2),0)
      msq(0,-j)=msq(0,-(j-2))
      enddo

      return
      end

