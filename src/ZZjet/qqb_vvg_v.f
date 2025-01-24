!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_vvg_v(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  q'(p4)+bar{q'}(p5) + n(p6)+ebar(p7)+ g(p3)
c   for the moment --- radiation only from initial line

c--- Includes both ZZ and WW contributions, taken from qqb_zzg_v.f & qqb_wwg_v.f

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'zcouple.f'
      include 'zcouple_cms.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'pchoice.f'
      include 'scheme.f'
      include 'blha.f'
      integer:: jk,hq,h34,h56,hg,ii
      real(dp):: faczz,fac1,q34,q56,s127,nup,ndn
      complex(dp):: facww
      real(dp):: p(mxpart,4),pww(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: ave
      complex(dp):: v34(2),v56(2),LR(2,nf)
      complex(dp)::
     & aq_0a(2,2,2,2),aq_0sr(2,2,2,2,2),aq_1a(2,2,2,2),aq_1q(2,2,2,2),aq_1ax56(2,2,2,2),aq_1ax34(2,2,2,2),aq_1sr(2,2,2,2,2),
     & qa_0a(2,2,2,2),qa_0sr(2,2,2,2,2),qa_1a(2,2,2,2),qa_1q(2,2,2,2),qa_1ax56(2,2,2,2),qa_1ax34(2,2,2,2),qa_1sr(2,2,2,2,2),
     & ag_0a(2,2,2,2),ag_0sr(2,2,2,2,2),ag_1a(2,2,2,2),ag_1q(2,2,2,2),ag_1ax56(2,2,2,2),ag_1ax34(2,2,2,2),ag_1sr(2,2,2,2,2),
     & qg_0a(2,2,2,2),qg_0sr(2,2,2,2,2),qg_1a(2,2,2,2),qg_1q(2,2,2,2),qg_1ax56(2,2,2,2),qg_1ax34(2,2,2,2),qg_1sr(2,2,2,2,2),
     & gq_0a(2,2,2,2),gq_0sr(2,2,2,2,2),gq_1a(2,2,2,2),gq_1q(2,2,2,2),gq_1ax56(2,2,2,2),gq_1ax34(2,2,2,2),gq_1sr(2,2,2,2,2),
     & ga_0a(2,2,2,2),ga_0sr(2,2,2,2,2),ga_1a(2,2,2,2),ga_1q(2,2,2,2),ga_1ax56(2,2,2,2),ga_1ax34(2,2,2,2),ga_1sr(2,2,2,2,2),
     & aq_1qmLL(2,2,2,2),qa_1qmLL(2,2,2,2),ag_1qmLL(2,2,2,2),qg_1qmLL(2,2,2,2),gq_1qmLL(2,2,2,2),ga_1qmLL(2,2,2,2),
     & aq_1qmLR(2,2,2,2),qa_1qmLR(2,2,2,2),ag_1qmLR(2,2,2,2),qg_1qmLR(2,2,2,2),gq_1qmLR(2,2,2,2),ga_1qmLR(2,2,2,2)
      complex(dp):: amp0,amp1,prop34,prop56,prop127
      complex(dp):: amp0zz(2,2,2,2),amp0ww(2,2,2,2),amp1zz(2,2,2,2),amp1ww(2,2,2,2)
      integer,parameter::i4(2)=(/4,6/),i6(2)=(/6,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      complex(dp):: propwp,propwm,cprop,c(2,2)
      complex(dp)::
     & gub0a(2,2),gub0b(2,2),gub0sr(2,2),gub1a(2,2),gub1b(2,2),gub1ax(2,2),gub1sr(2,2),gub1q(2,2),
     &  gu0a(2,2), gu0b(2,2), gu0sr(2,2), gu1a(2,2), gu1b(2,2), gu1ax(2,2), gu1sr(2,2), gu1q(2,2),
     & gdb0a(2,2),gdb0b(2,2),gdb0sr(2,2),gdb1a(2,2),gdb1b(2,2),gdb1ax(2,2),gdb1sr(2,2),gdb1q(2,2),
     &  gd0a(2,2), gd0b(2,2), gd0sr(2,2), gd1a(2,2), gd1b(2,2), gd1ax(2,2), gd1sr(2,2), gd1q(2,2),
     & ubg0a(2,2),ubg0b(2,2),ubg0sr(2,2),ubg1a(2,2),ubg1b(2,2),ubg1ax(2,2),ubg1sr(2,2),ubg1q(2,2),
     &  ug0a(2,2), ug0b(2,2), ug0sr(2,2), ug1a(2,2), ug1b(2,2), ug1ax(2,2), ug1sr(2,2), ug1q(2,2),
     & dbg0a(2,2),dbg0b(2,2),dbg0sr(2,2),dbg1a(2,2),dbg1b(2,2),dbg1ax(2,2),dbg1sr(2,2),dbg1q(2,2),
     &  dg0a(2,2), dg0b(2,2), dg0sr(2,2), dg1a(2,2), dg1b(2,2), dg1ax(2,2), dg1sr(2,2), dg1q(2,2),
     & ubu0a(2,2),ubu0b(2,2),ubu0sr(2,2),ubu1a(2,2),ubu1b(2,2),ubu1ax(2,2),ubu1sr(2,2),ubu1q(2,2),
     & uub0a(2,2),uub0b(2,2),uub0sr(2,2),uub1a(2,2),uub1b(2,2),uub1ax(2,2),uub1sr(2,2),uub1q(2,2),
     & dbd0a(2,2),dbd0b(2,2),dbd0sr(2,2),dbd1a(2,2),dbd1b(2,2),dbd1ax(2,2),dbd1sr(2,2),dbd1q(2,2),
     & ddb0a(2,2),ddb0b(2,2),ddb0sr(2,2),ddb1a(2,2),ddb1b(2,2),ddb1ax(2,2),ddb1sr(2,2),ddb1q(2,2),
     & amp0a(2,2),amp0b(2,2),amp0sr(2,2),amp1b(2,2),amp1sr(2,2),amp1a(2,2),amp1ax(2,2)
      integer:: minus,mplus,polg,polq,tjk
      parameter(minus=1,mplus=2)
      real(dp), parameter :: mp(nf)=(/-1d0,+1d0,-1d0,+1d0,-1d0/)

      scheme='tH-V'
c--- note that scheme = "tH-V" is also a valid choice in this case
c--- (finite terms relating the two schemes are included in, e.g. AWWjeta.f)

c-set Uncrossed array to zero
c      call writeout(p)
      msq(:,:)=zip

      fac1=two*cf*xn*gsq

c!!! Begin ZZ block
c!!! Begin ZZ block
c!!! Begin ZZ block

      faczz=-four*esq**2

      v34(1)=zl1
      v34(2)=zr1
      q34=q1
      v56(1)=zl2
      v56(2)=zr2
      q56=q2

      LR(1,:)=zL(:)
      LR(2,:)=zR(:)

c--   s returned from sprodx (common block) is 2*dot product
      call spinoru(7,p,za,zb)

      ii=1

c--   calculate propagators
      prop34=s(3,i4(ii))/cplx2(s(3,i4(ii))-zmass**2,zmass*zwidth)
      prop56=s(5,i6(ii))/cplx2(s(5,i6(ii))-zmass**2,zmass*zwidth)
c--- Amplitude returned with arguments (hq,h34,h56,h7)
c---case qbar-q
      if (useblha == 0) then
      call amp_ZZg_v(p,1,2,3,i4(ii),i6(ii),5,7,aq_0a,aq_0sr,aq_1a,aq_1ax56,aq_1ax34,aq_1sr)
      endif
c---case q-qbar
      call amp_ZZg_v(p,2,1,3,i4(ii),i6(ii),5,7,qa_0a,qa_0sr,qa_1a,qa_1ax56,qa_1ax34,qa_1sr)
      call amp_ZZg_v_qloop(p,1,2,3,i4(ii),5,i6(ii),7,aq_1q,aq_1qmLL,aq_1qmLR,qa_1q,qa_1qmLL,qa_1qmLR)

      if (useblha == 0) then
c---case qbar-g
      call amp_ZZg_v(p,1,7,3,i4(ii),i6(ii),5,2,ag_0a,ag_0sr,ag_1a,ag_1ax56,ag_1ax34,ag_1sr)
c---case q-g
      call amp_ZZg_v(p,7,1,3,i4(ii),i6(ii),5,2,qg_0a,qg_0sr,qg_1a,qg_1ax56,qg_1ax34,qg_1sr)
      call amp_ZZg_v_qloop(p,1,7,3,i4(ii),5,i6(ii),2,ag_1q,ag_1qmLL,ag_1qmLR,qg_1q,qg_1qmLL,qg_1qmLR)

c---case g-q
      call amp_ZZg_v(p,7,2,3,i4(ii),i6(ii),5,1,gq_0a,gq_0sr,gq_1a,gq_1ax56,gq_1ax34,gq_1sr)
c---case g-qbar
      call amp_ZZg_v(p,2,7,3,i4(ii),i6(ii),5,1,ga_0a,ga_0sr,ga_1a,ga_1ax56,ga_1ax34,ga_1sr)
      call amp_ZZg_v_qloop(p,7,2,3,i4(ii),5,i6(ii),1,gq_1q,gq_1qmLL,gq_1qmLR,ga_1q,ga_1qmLL,ga_1qmLR)
      endif

c!!! End ZZ block
c!!! End ZZ block
c!!! End ZZ block

c!!! Begin WW block
c!!! Begin WW block
c!!! Begin WW block

c Change the momenta to DKS notation: e-(3) nu~(4) e+(5) nu(6)
      pww(1,:)=p(1,:)
      pww(2,:)=p(2,:)
      pww(3,:)=p(3,:) ! interchanged 3<->5 wrt. qqb_wwg_v.f
      pww(4,:)=p(6,:)
      pww(5,:)=p(4,:)
      pww(6,:)=p(5,:) ! interchanged 3<->5 wrt. qqb_wwg_v.f
      pww(7,:)=p(7,:)

c      fac=gw**4
      facww=esq**2/zxw**2

      call spinoru(7,pww,za,zb)

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7)
      propwp=cplx1(s(3,4))/cplx2(s(3,4)-wmass**2,wmass*wwidth)
      propwm=cplx1(s(5,6))/cplx2(s(5,6)-wmass**2,wmass*wwidth)
      prop127=s127/cplx2(s127-zmass**2,zmass*zwidth)
      cprop=1d0
      cprop=cprop*propwp*propwm

c-- couplings according to hep-ph/9803250 Eqs. 3.4 and 3.6
c-- First argument is left or right, second is d(1) or u(2)
      do j=1,2
        c(minus,j)=+mp(j)*2d0*Q(j)*zxw+mp(j)*zL(j)*zsin2w*prop127
        c(mplus,j)=+mp(j)*2d0*Q(j)*zxw-mp(j)*2d0*Q(j)*zxw*prop127
      enddo

c--- Fill leading order and virtual amplitudes
      if (useblha == 0) then
c--- gluon-quark and gluon-antiquark
      call amp_WWg_v(pww,2,7,3,4,5,6,1,2,gub0a,gub0b,gub0sr,gub1a,gub1b,gub1ax,gub1sr)
      call amp_WWg_v(pww,7,2,3,4,5,6,1,2, gu0a, gu0b, gu0sr, gu1a, gu1b, gu1ax, gu1sr)
      call amp_WWg_v(pww,2,7,6,5,4,3,1,1,gdb0a,gdb0b,gdb0sr,gdb1a,gdb1b,gdb1ax,gdb1sr)
      call amp_WWg_v(pww,7,2,6,5,4,3,1,1, gd0a, gd0b, gd0sr, gd1a, gd1b, gd1ax, gd1sr)
      call amp_wwg_v_qloop(2,7,3,4,5,6,1,gub1q,gu1q,gdb1q,gd1q,gub1ax,gu1ax,gdb1ax,gd1ax)

c--- quark-gluon and antiquark-gluon
      call amp_WWg_v(pww,1,7,3,4,5,6,2,2,ubg0a,ubg0b,ubg0sr,ubg1a,ubg1b,ubg1ax,ubg1sr)
      call amp_WWg_v(pww,7,1,3,4,5,6,2,2, ug0a, ug0b, ug0sr, ug1a, ug1b, ug1ax, ug1sr)
      call amp_WWg_v(pww,1,7,6,5,4,3,2,1,dbg0a,dbg0b,dbg0sr,dbg1a,dbg1b,dbg1ax,dbg1sr)
      call amp_WWg_v(pww,7,1,6,5,4,3,2,1, dg0a, dg0b, dg0sr, dg1a, dg1b, dg1ax, dg1sr)
      call amp_wwg_v_qloop(1,7,3,4,5,6,2,ubg1q,ug1q,dbg1q,dg1q,ubg1ax,ug1ax,dbg1ax,dg1ax)

c--- quark-antiquark and antiquark-quark
      call amp_WWg_v(pww,1,2,3,4,5,6,7,2,ubu0a,ubu0b,ubu0sr,ubu1a,ubu1b,ubu1ax,ubu1sr)
      call amp_WWg_v(pww,1,2,6,5,4,3,7,1,dbd0a,dbd0b,dbd0sr,dbd1a,dbd1b,dbd1ax,dbd1sr)
      endif
      if ((useblha == 0) .or. (blhatype == 2)) then
      call amp_WWg_v(pww,2,1,3,4,5,6,7,2,uub0a,uub0b,uub0sr,uub1a,uub1b,uub1ax,uub1sr)
      endif
      if ((useblha == 0) .or. (blhatype == 1)) then
      call amp_WWg_v(pww,2,1,6,5,4,3,7,1,ddb0a,ddb0b,ddb0sr,ddb1a,ddb1b,ddb1ax,ddb1sr)
      endif
      call amp_wwg_v_qloop(1,2,3,4,5,6,7,ubu1q,uub1q,dbd1q,ddb1q,ubu1ax,uub1ax,dbd1ax,ddb1ax)

c!!! End WW block
c!!! End WW block
c!!! End WW block

c---calculate over a limited flavour range (-2:2)
      do j=-2,2
      do k=-2,2
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) goto 19

c---determine averaging factor for different channels
      if ((j == 0) .or. (k == 0)) then
        jk=j+k
        ave=aveqg
      else
        jk=max(j,k)
        ave=aveqq
      endif

      amp0zz(:,:,:,:)=czip
      amp1zz(:,:,:,:)=czip
      amp0ww(:,:,:,:)=czip
      amp1ww(:,:,:,:)=czip

c numer of quarks circulating in loop
      nup=2._dp
      ndn=3._dp

      if (jk == 0) goto 19

c!! Begin filling of ZZ amplitudes
      do hq=1,2
      do h34=1,2
      do h56=1,2
      do hg=1,2

      amp0=zip
      amp1=zip

c---case qbar-q
      if    ((j < 0).and.(k > 0)) then
      amp0=(prop56*v56(h56)*LR(hq,k)+q56*Q(k))
     &    *(prop34*v34(h34)*LR(hq,k)+q34*Q(k))*aq_0a(hq,h34,h56,hg)
     &    +aq_0sr(k,hq,h34,h56,hg)
      amp1=(prop56*v56(h56)*LR(hq,k)+q56*Q(k))
     &    *(prop34*v34(h34)*LR(hq,k)+q34*Q(k))*aq_1a(hq,h34,h56,hg)
     &    +((prop56*v56(h56)*zL(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zL(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zR(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zR(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))*nup
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))*nup)*aq_1q(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2)))*aq_1qmLL(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2)))*aq_1qmLR(hq,h34,h56,hg)
     &    +(prop34*v34(h34)*LR(hq,k)+q34*Q(k))
     &     *prop56*v56(h56)*aq_1ax56(hq,h34,h56,hg)
     &    +(prop56*v56(h56)*LR(hq,k)+q56*Q(k))
     &     *prop34*v34(h34)*aq_1ax34(hq,h34,h56,hg)
     &    +aq_1sr(k,hq,h34,h56,hg)

cc---case q-qbar
      elseif((j > 0).and.(k < 0)) then
      amp0=(prop56*v56(h56)*LR(hq,j)+q56*Q(j))
     &    *(prop34*v34(h34)*LR(hq,j)+q34*Q(j))*qa_0a(hq,h34,h56,hg)
     &    +qa_0sr(j,hq,h34,h56,hg)
      amp1=(prop56*v56(h56)*LR(hq,j)+q56*Q(j))
     &    *(prop34*v34(h34)*LR(hq,j)+q34*Q(j))*qa_1a(hq,h34,h56,hg)
     &    +((prop56*v56(h56)*zL(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zL(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zR(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zR(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))*nup
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))*nup)*qa_1q(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2)))*qa_1qmLL(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2)))*qa_1qmLR(hq,h34,h56,hg)
     &    +(prop34*v34(h34)*LR(hq,j)+q34*Q(j))
     &     *prop56*v56(h56)*qa_1ax56(hq,h34,h56,hg)
     &    +(prop56*v56(h56)*LR(hq,j)+q56*Q(j))
     &     *prop34*v34(h34)*qa_1ax34(hq,h34,h56,hg)
     &    +qa_1sr(j,hq,h34,h56,hg)

cc---case qbar-g
      elseif((j < 0).and.(k == 0)) then
      amp0=(prop56*v56(h56)*LR(hq,-jk)+q56*Q(-jk))
     &    *(prop34*v34(h34)*LR(hq,-jk)+q34*Q(-jk))*ag_0a(hq,h34,h56,hg)
     &    +ag_0sr(-j,hq,h34,h56,hg)
      amp1=(prop56*v56(h56)*LR(hq,-jk)+q56*Q(-jk))
     &    *(prop34*v34(h34)*LR(hq,-jk)+q34*Q(-jk))*ag_1a(hq,h34,h56,hg)
     &    +((prop56*v56(h56)*zL(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zL(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zR(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zR(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))*nup
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))*nup)*ag_1q(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2)))*ag_1qmLL(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2)))*ag_1qmLR(hq,h34,h56,hg)
     &    +(prop34*v34(h34)*LR(hq,-jk)+q34*Q(-jk))
     &     *prop56*v56(h56)*ag_1ax56(hq,h34,h56,hg)
     &    +(prop56*v56(h56)*LR(hq,-jk)+q56*Q(-jk))
     &     *prop34*v34(h34)*ag_1ax34(hq,h34,h56,hg)
     &    +ag_1sr(-j,hq,h34,h56,hg)

cc---case g-qbar
      elseif((k < 0).and.(j == 0)) then
      amp0=(prop56*v56(h56)*LR(hq,-jk)+q56*Q(-jk))
     &    *(prop34*v34(h34)*LR(hq,-jk)+q34*Q(-jk))*ga_0a(hq,h34,h56,hg)
     &    +ga_0sr(-k,hq,h34,h56,hg)
      amp1=(prop56*v56(h56)*LR(hq,-jk)+q56*Q(-jk))
     &    *(prop34*v34(h34)*LR(hq,-jk)+q34*Q(-jk))*ga_1a(hq,h34,h56,hg)
     &    +((prop56*v56(h56)*zL(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zL(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zR(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zR(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))*nup
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))*nup)*ga_1q(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2)))*ga_1qmLL(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2)))*ga_1qmLR(hq,h34,h56,hg)
     &    +(prop34*v34(h34)*LR(hq,-jk)+q34*Q(-jk))
     &     *prop56*v56(h56)*ga_1ax56(hq,h34,h56,hg)
     &    +(prop56*v56(h56)*LR(hq,-jk)+q56*Q(-jk))
     &     *prop34*v34(h34)*ga_1ax34(hq,h34,h56,hg)
     &    +ga_1sr(-k,hq,h34,h56,hg)

cc---case q-g
      elseif((j > 0).and.(k == 0)) then
      amp0=(prop56*v56(h56)*LR(hq,jk)+q56*Q(jk))
     &    *(prop34*v34(h34)*LR(hq,jk)+q34*Q(jk))*qg_0a(hq,h34,h56,hg)
     &    +qg_0sr(j,hq,h34,h56,hg)
      amp1=(prop56*v56(h56)*LR(hq,jk)+q56*Q(jk))
     &    *(prop34*v34(h34)*LR(hq,jk)+q34*Q(jk))*qg_1a(hq,h34,h56,hg)
     &    +((prop56*v56(h56)*zL(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zL(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zR(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zR(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))*nup
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))*nup)*qg_1q(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2)))*qg_1qmLL(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2)))*qg_1qmLR(hq,h34,h56,hg)
     &    +(prop34*v34(h34)*LR(hq,jk)+q34*Q(jk))
     &     *prop56*v56(h56)*qg_1ax56(hq,h34,h56,hg)
     &    +(prop56*v56(h56)*LR(hq,jk)+q56*Q(jk))
     &     *prop34*v34(h34)*qg_1ax34(hq,h34,h56,hg)
     &    +qg_1sr(j,hq,h34,h56,hg)

      elseif((k > 0).and.(j == 0)) then
cc---case g-q
      amp0=(prop56*v56(h56)*LR(hq,jk)+q56*Q(jk))
     &    *(prop34*v34(h34)*LR(hq,jk)+q34*Q(jk))*gq_0a(hq,h34,h56,hg)
     &    +gq_0sr(k,hq,h34,h56,hg)
      amp1=(prop56*v56(h56)*LR(hq,jk)+q56*Q(jk))
     &    *(prop34*v34(h34)*LR(hq,jk)+q34*Q(jk))*gq_1a(hq,h34,h56,hg)
     &    +((prop56*v56(h56)*zL(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zL(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zR(1)+q56*Q(1))
     &     *(prop34*v34(h34)*zR(1)+q34*Q(1))*ndn
     &     +(prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))*nup
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))*nup)*gq_1q(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2)))*gq_1qmLL(hq,h34,h56,hg)*2._dp
     &    +((prop56*v56(h56)*zL(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zR(2)+q34*Q(2))
     &     +(prop56*v56(h56)*zR(2)+q56*Q(2))
     &     *(prop34*v34(h34)*zL(2)+q34*Q(2)))*gq_1qmLR(hq,h34,h56,hg)
     &    +(prop34*v34(h34)*LR(hq,jk)+q34*Q(jk))
     &     *prop56*v56(h56)*gq_1ax56(hq,h34,h56,hg)
     &    +(prop56*v56(h56)*LR(hq,jk)+q56*Q(jk))
     &     *prop34*v34(h34)*gq_1ax34(hq,h34,h56,hg)
     &    +gq_1sr(k,hq,h34,h56,hg)

      endif

      amp0zz(hq,h34,h56,hg)=amp0*faczz
      amp1zz(hq,h34,h56,hg)=amp1*faczz

      enddo  ! endloop hg
      enddo  ! endloop h56
      enddo  ! endloop h34
      enddo  ! endloop hq

c End filling of ZZ amplitudes

c Begin filling of WW amplitudes

c--- tjk is related to the isospin of the quark in the initial state
c--- and is equal to 2 (u,c) or 1 (d,s,b)
      tjk=2-mod(abs(jk),2)

c--- identify appropriate amplitudes as amp0a, etc.
c--- contribution of closed fermion loops; add to "a" since no further
c--- coupling factors are required, overall factor the same
      if    (j  <  0 .and. tau(jk)  ==  -1d0 .and. k  /=  0) then
        amp0a(:,:)=dbd0a(:,:)
        amp0b(:,:)=dbd0b(:,:)
        amp0sr(:,:)=dbd0sr(:,:)
        amp1a(:,:)=dbd1a(:,:)+dbd1q(:,:)
        amp1b(:,:)=dbd1b(:,:)
        amp1ax(:,:)=dbd1ax(:,:)
        amp1sr(:,:)=dbd1sr(:,:)
      elseif(j  <  0 .and. tau(jk)  ==   1d0 .and. k  /=  0) then
        amp0a(:,:)=ubu0a(:,:)
        amp0b(:,:)=ubu0b(:,:)
        amp0sr(:,:)=ubu0sr(:,:)
        amp1a(:,:)=ubu1a(:,:)+ubu1q(:,:)
        amp1b(:,:)=ubu1b(:,:)
        amp1ax(:,:)=ubu1ax(:,:)
        amp1sr(:,:)=ubu1sr(:,:)
      elseif(j  >  0 .and. tau(jk)  ==  -1d0 .and. k  /=  0) then
        amp0a(:,:)=ddb0a(:,:)
        amp0b(:,:)=ddb0b(:,:)
        amp0sr(:,:)=ddb0sr(:,:)
        amp1a(:,:)=ddb1a(:,:)+ddb1q(:,:)
        amp1b(:,:)=ddb1b(:,:)
        amp1ax(:,:)=ddb1ax(:,:)
        amp1sr(:,:)=ddb1sr(:,:)
      elseif(j  >  0 .and. tau(jk)  ==   1d0 .and. k  /=  0) then
        amp0a(:,:)=uub0a(:,:)
        amp0b(:,:)=uub0b(:,:)
        amp0sr(:,:)=uub0sr(:,:)
        amp1a(:,:)=uub1a(:,:)+uub1q(:,:)
        amp1b(:,:)=uub1b(:,:)
        amp1ax(:,:)=uub1ax(:,:)
        amp1sr(:,:)=uub1sr(:,:)
      elseif(j  ==  0 .and. tau(jk)  ==   1d0 .and. jk  >  0) then
        amp0a(:,:)=gu0a(:,:)
        amp0b(:,:)=gu0b(:,:)
        amp0sr(:,:)=gu0sr(:,:)
        amp1a(:,:)=gu1a(:,:)+gu1q(:,:)
        amp1b(:,:)=gu1b(:,:)
        amp1ax(:,:)=gu1ax(:,:)
        amp1sr(:,:)=gu1sr(:,:)
      elseif(j  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  >  0) then
        amp0a(:,:)=gd0a(:,:)
        amp0b(:,:)=gd0b(:,:)
        amp0sr(:,:)=gd0sr(:,:)
        amp1a(:,:)=gd1a(:,:)+gd1q(:,:)
        amp1b(:,:)=gd1b(:,:)
        amp1ax(:,:)=gd1ax(:,:)
        amp1sr(:,:)=gd1sr(:,:)
      elseif(j  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  <  0) then
        amp0a(:,:)=gub0a(:,:)
        amp0b(:,:)=gub0b(:,:)
        amp0sr(:,:)=gub0sr(:,:)
        amp1a(:,:)=gub1a(:,:)+gub1q(:,:)
        amp1b(:,:)=gub1b(:,:)
        amp1ax(:,:)=gub1ax(:,:)
        amp1sr(:,:)=gub1sr(:,:)
      elseif(j  ==  0 .and. tau(jk)  ==   1d0 .and. jk  <  0) then
        amp0a(:,:)=gdb0a(:,:)
        amp0b(:,:)=gdb0b(:,:)
        amp0sr(:,:)=gdb0sr(:,:)
        amp1a(:,:)=gdb1a(:,:)+gdb1q(:,:)
        amp1b(:,:)=gdb1b(:,:)
        amp1ax(:,:)=gdb1ax(:,:)
        amp1sr(:,:)=gdb1sr(:,:)
      elseif(k  ==  0 .and. tau(jk)  ==   1d0 .and. jk  >  0) then
        amp0a(:,:)=ug0a(:,:)
        amp0b(:,:)=ug0b(:,:)
        amp0sr(:,:)=ug0sr(:,:)
        amp1a(:,:)=ug1a(:,:)+ug1q(:,:)
        amp1b(:,:)=ug1b(:,:)
        amp1ax(:,:)=ug1ax(:,:)
        amp1sr(:,:)=ug1sr(:,:)
      elseif(k  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  >  0) then
        amp0a(:,:)=dg0a(:,:)
        amp0b(:,:)=dg0b(:,:)
        amp0sr(:,:)=dg0sr(:,:)
        amp1a(:,:)=dg1a(:,:)+dg1q(:,:)
        amp1b(:,:)=dg1b(:,:)
        amp1ax(:,:)=dg1ax(:,:)
        amp1sr(:,:)=dg1sr(:,:)
      elseif(k  ==  0 .and. tau(jk)  ==  -1d0 .and. jk  <  0) then
        amp0a(:,:)=ubg0a(:,:)
        amp0b(:,:)=ubg0b(:,:)
        amp0sr(:,:)=ubg0sr(:,:)
        amp1a(:,:)=ubg1a(:,:)+ubg1q(:,:)
        amp1b(:,:)=ubg1b(:,:)
        amp1ax(:,:)=ubg1ax(:,:)
        amp1sr(:,:)=ubg1sr(:,:)
      elseif(k  ==  0 .and. tau(jk)  ==   1d0 .and. jk  <  0) then
        amp0a(:,:)=dbg0a(:,:)
        amp0b(:,:)=dbg0b(:,:)
        amp0sr(:,:)=dbg0sr(:,:)
        amp1a(:,:)=dbg1a(:,:)+dbg1q(:,:)
        amp1b(:,:)=dbg1b(:,:)
        amp1ax(:,:)=dbg1ax(:,:)
        amp1sr(:,:)=dbg1sr(:,:)
      endif

c--- assemble amplitudes with appropriate couplings
      do polq=1,2
      do polg=1,2
        amp0ww(polq,1,1,polg)=facww*(
     &   cprop*(amp0a(polq,polg)+amp0b(polq,polg)*c(polq,tjk))+zxw*amp0sr(polq,polg))
        amp1ww(polq,1,1,polg)=facww*(
     &   cprop*(amp1a(polq,polg)+amp1b(polq,polg)*c(polq,tjk)+amp1ax(polq,polg)*tau(tjk))+zxw*amp1sr(polq,polg))
      enddo
      enddo

c End filling of WW amplitudes

c Now sum and square amplitudes
      do hq=1,2
      do h34=1,2
      do h56=1,2
      do hg=1,2

c--- normal case
c--- (check of LO)
c        msq(j,k)=msq(j,k)+fac1*ave*real(
c     &         (amp0zz(hq,h34,h56,hg)+amp0ww(hq,h34,h56,hg))
c     &   *conjg(amp0zz(hq,h34,h56,hg)+amp0ww(hq,h34,h56,hg)),dp)
c--- (virtual)
        msq(j,k)=msq(j,k)+ason2pi*fac1*ave*real(
     &         (amp1zz(hq,h34,h56,hg)+amp1ww(hq,h34,h56,hg))
     &   *conjg(amp0zz(hq,h34,h56,hg)+amp0ww(hq,h34,h56,hg)),dp)

      enddo  ! endloop hg
      enddo  ! endloop h56
      enddo  ! endloop h34
      enddo  ! endloop hq


   19 continue
      enddo  !endloop j
      enddo  !endloop k


c---extend to full flavour range
      do j=-nf,nf
      do k=-nf,nf
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) goto 20

      msq(j,k)=msq(jkswitch(j),jkswitch(k))

 20   continue
      enddo
      enddo

      return
      end







