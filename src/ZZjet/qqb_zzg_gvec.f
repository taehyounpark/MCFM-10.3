!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_zzg_gvec(p,n,in,msq)
      implicit none
c----Matrix element for ZZ production
c----averaged over initial colours and spins
c    contracted with the vector n(mu)
c     q(-p1)+qbar(-p2)--> g(p7)+ Z(e-(p3)+e^+(p4))+ Z(mu-(p5)+mu+(p6))
c---
c---ip emitter
c---kp spectator
c---in label of gluon which is contracted with n
      include 'types.f'
      include 'cplx.h'
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
      include 'srdiags.f'
      integer in,hq,h34,h56,jk
      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),n(4),
     & q34,q56,ave,s127,fac,fac1,rescale1,rescale2
      complex(dp):: v34(2),v56(2),LR(2,2)
      complex(dp)::
     & qa12(2,2,2),qa34(2,2,2),qa56(2,2,2),
     & aq12(2,2,2),aq34(2,2,2),aq56(2,2,2),
     & ga12(2,2,2),ga34(2,2,2),ga56(2,2,2),
     & ag12(2,2,2),ag34(2,2,2),ag56(2,2,2),
     & qg12(2,2,2),qg34(2,2,2),qg56(2,2,2),
     & gq12(2,2,2),gq34(2,2,2),gq56(2,2,2),
     & amp,prop34,prop56,prop127
      integer,parameter::jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      msq(:,:)=zip
c     overall factor removed is 2*rt2*e^4*gs*T^C7

c      fac=8._dp*V*esq**4*gsq
      fac=four*esq**2
      fac1=half*V*gsq

      q34=q1
      v34(1)=zl1
      v34(2)=zr1
      q56=q2
      v56(1)=zl2
      v56(2)=zr2

      LR(1,1)=zL(1)
      LR(1,2)=zL(2)
      LR(2,1)=zR(1)
      LR(2,2)=zR(2)

      rescale1=1._dp
      rescale2=1._dp

      call spinoru(7,p,za,zb)

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7)
      prop34=s(3,4)/cplx2(s(3,4)-zmass**2,zmass*zwidth)
      prop56=s(5,6)/cplx2(s(5,6)-zmass**2,zmass*zwidth)
      prop127=s127/cplx2(s127-zmass**2,zmass*zwidth)

      if     (in  ==  1) then
        call   a7ZZn(7,2,3,4,5,6,1,p,n,za,zb,ga12)
        call   a7ZZn(2,7,3,4,5,6,1,p,n,za,zb,gq12)
        if(srdiags) then
        call a7ZZsrn(7,2,3,4,5,6,1,p,n,za,zb,ga34,ga56)
        call a7ZZsrn(2,7,3,4,5,6,1,p,n,za,zb,gq34,gq56)
        endif
      elseif (in  ==  2) then
        call   a7ZZn(1,7,3,4,5,6,2,p,n,za,zb,qg12)
        call   a7ZZn(7,1,3,4,5,6,2,p,n,za,zb,ag12)
        if(srdiags) then
        call a7ZZsrn(1,7,3,4,5,6,2,p,n,za,zb,qg34,qg56)
        call a7ZZsrn(7,1,3,4,5,6,2,p,n,za,zb,ag34,ag56)
        endif
      elseif (in  ==  7) then
        call   a7ZZn(1,2,3,4,5,6,7,p,n,za,zb,qa12)
        call   a7ZZn(2,1,3,4,5,6,7,p,n,za,zb,aq12)
        if(srdiags) then
        call a7ZZsrn(1,2,3,4,5,6,7,p,n,za,zb,qa34,qa56)
        call a7ZZsrn(2,1,3,4,5,6,7,p,n,za,zb,aq34,aq56)
        endif
      else
        write(6,*) 'Invalid value of in=',in,' in qqb_zzg_gvec.f!'
        stop
      endif

c      write(6,*) 'qqb_zzg_gvec:srdiags',srdiags

c---calculate over a limited flavour range (-2:2)
      do j=-2,2
      do k=-2,2
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) goto 19

c---determine averaging factor for different channels
c     vsymfact=symmetry factor
      if ((j == 0) .or. (k == 0)) then
        jk=j+k
        ave=aveqg
      else
        jk=max(j,k)
        ave=aveqq
      endif

      if (jk == 0) goto 19

      do hq=1,2
      do h34=1,2
      do h56=1,2

c---case qbar-q
      if    ((j < 0).and.(k > 0)) then
      amp=(prop56*v56(h56)*lr(hq,k)+q56*q(k))
     & *(prop34*v34(h34)*lr(hq,k)+q34*q(k))*aq12(hq,h34,h56)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,k)+q34*q(k))*aq56(hq,h34,h56)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,k)+q56*q(k))*aq34(hq,h34,h56)
      endif

c---case q-qbar
      elseif((j > 0).and.(k < 0)) then
      amp=(prop56*v56(h56)*lr(hq,j)+q56*q(j))
     & *(prop34*v34(h34)*lr(hq,j)+q34*q(j))*qa12(hq,h34,h56)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,j)+q34*q(j))*qa56(hq,h34,h56)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,j)+q56*q(j))*qa34(hq,h34,h56)
      endif

c---case qbar-g
      elseif((j < 0).and.(k == 0)) then
      amp=(prop56*v56(h56)*lr(hq,-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*lr(hq,-jk)+q34*q(-jk))*ag12(hq,h34,h56)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,-jk)+q34*q(-jk))*ag56(hq,h34,h56)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,-jk)+q56*q(-jk))*ag34(hq,h34,h56)
      endif

c---case g-qbar
      elseif((k < 0).and.(j == 0)) then
      amp=(prop56*v56(h56)*lr(hq,-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*lr(hq,-jk)+q34*q(-jk))*ga12(hq,h34,h56)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,-jk)+q34*q(-jk))*ga56(hq,h34,h56)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,-jk)+q56*q(-jk))*ga34(hq,h34,h56)
      endif

c---case q-g
      elseif((j > 0).and.(k == 0)) then
      amp=(prop56*v56(h56)*lr(hq,+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*lr(hq,+jk)+q34*q(+jk))*qg12(hq,h34,h56)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,+jk)+q34*q(+jk))*qg56(hq,h34,h56)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,+jk)+q56*q(+jk))*qg34(hq,h34,h56)
      endif

      elseif((k > 0).and.(j == 0)) then
c---case g-q
      amp=(prop56*v56(h56)*lr(hq,+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*lr(hq,+jk)+q34*q(+jk))*gq12(hq,h34,h56)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,+jk)+q34*q(+jk))*gq56(hq,h34,h56)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,+jk)+q56*q(+jk))*gq34(hq,h34,h56)
      endif

      endif

      amp=amp*fac

c      msq(j,k)=msq(j,k)+fac1*ave*abs(amp)**2
c      !set-up interference terms
c      if ((interference).and.(ii==1)) then
c      Uncrossed(j,k,hq,h34,h56)=amp
c      elseif (ii==2) then
c      if (h34==h56) then
c      msq(j,k)=msq(j,k)
c     &    -fac1*two*ave*real(conjg(amp)*Uncrossed(j,k,hq,h34,h56))
c      endif
c      endif

      msq(j,k)=msq(j,k)+fac1*ave*abs(amp)**2

      enddo  ! endloop h56
      enddo  ! endloop h34
      enddo  ! endloop hq

   19 continue
      enddo  !endloop j
      enddo  !endloop k


c---extend to full flavour range
      do j=-nf,nf
      do k=-nf,nf
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) cycle
      msq(j,k)=msq(jkswitch(j),jkswitch(k))
      enddo
      enddo

      return
      end

