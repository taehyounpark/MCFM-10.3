!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_zzg_v(p,msq)
      implicit none
      include 'types.f'
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  q'(p4)+bar{q'}(p5) + n(p6)+ebar(p7)+ g(p3)
c   for the moment --- radiation only from initial line

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
      include 'interference.f'
      include 'pchoice.f'
      include 'scheme.f'
      include 'blha.f'
      integer:: jk,hq,h34,h56,hg,ii,nmax
      real(dp):: fac,fac1,q34,q56,nup,ndn
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
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
      complex(dp):: Uncrossed0(-nf:nf,-nf:nf,2,2,2,2)
      complex(dp):: Uncrossed1(-nf:nf,-nf:nf,2,2,2,2)
      complex(dp):: amp0,amp1,prop34,prop56
      integer,parameter::i4(2)=(/4,6/),i6(2)=(/6,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

c Process with ZZ/WW interference
      if (blhatype > 0) then
        call qqb_vvg_v(p,msq)
        return
      endif

      scheme='tH-V'
c--- note that scheme = "tH-V" is also a valid choice in this case
c--- (finite terms relating the two schemes are included in, e.g. AWWjeta.f)

c-set Uncrossed array to zero
c      call writeout(p)
      msq(:,:)=zip
      Uncrossed0(:,:,:,:,:,:)=czip
      Uncrossed1(:,:,:,:,:,:)=czip

      fac=-four*esq**2
      fac1=two*cf*xn*gsq

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

      if (interference) then
         nmax=2
      else
         nmax=1
      endif

      do ii=1,nmax

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

c---calculate over a limited flavour range (-2:2)
      do j=-2,2
      do k=-2,2
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) goto 19

c---determine averaging factor for different channels
c     vsymfact=symmetry factor
      if ((j == 0) .or. (k == 0)) then
        jk=j+k
        ave=aveqg*vsymfact
      else
        jk=max(j,k)
        ave=aveqq*vsymfact
      endif


c numer of quarks circulating in loop
      nup=2._dp
      ndn=3._dp

      if (jk == 0) goto 19

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

      amp0=amp0*fac
      amp1=amp1*fac

      if (interference .eqv. .false.) then
c--- normal case
c--- (check of LO)
c        msq(j,k)=msq(j,k)+fac1*ave*abs(amp0)**2
c--- (virtual)
        msq(j,k)=msq(j,k)+ason2pi*fac1*ave*real(amp1*conjg(amp0),dp)
      else
cc--- with interference:
cc---    1st pass --> store result
cc---    2nd pass --> fill msq
        if (ii == 1) then
           Uncrossed0(j,k,hq,h34,h56,hg)=amp0
           Uncrossed1(j,k,hq,h34,h56,hg)=amp1
        else
           if (h34 == h56) then
             msq(j,k)=msq(j,k)+ason2pi*fac1*ave*real(
     &        (amp1-Uncrossed1(j,k,hq,h34,h56,hg))*conjg(amp0-Uncrossed0(j,k,hq,h34,h56,hg)),dp)
           else
             msq(j,k)=msq(j,k)+ason2pi*fac1*ave*real(
     &         amp1*conjg(amp0)
     &        +Uncrossed1(j,k,hq,h34,h56,hg)*conjg(Uncrossed0(j,k,hq,h34,h56,hg)),dp)
           endif
        endif
      endif

      enddo  ! endloop hg
      enddo  ! endloop h56
      enddo  ! endloop h34
      enddo  ! endloop hq

   19 continue
      enddo  !endloop j
      enddo  !endloop k
      enddo  !endloop ii


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







