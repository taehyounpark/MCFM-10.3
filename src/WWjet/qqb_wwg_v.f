!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wwg_v(pin,msq)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'nf.f'
      include 'masses.f'
      include 'scheme.f'
      include 'sprods_com.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'qcdcouple.f'
      include 'zcouple_cms.f'
      include 'blha.f'
      integer j,k,jk,tjk,minus,mplus,polq,polg
      real(dp):: pin(mxpart,4),p(mxpart,4),msq(-nf:nf,-nf:nf),
     & fac,fac1,s127,ave
      complex(dp):: propwp,propwm,prop127,cprop,c(2,2)
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
     & amp0a(2,2),amp0b(2,2),amp0sr(2,2),amp1a(2,2),amp1b(2,2),amp1ax(2,2),amp1sr(2,2),
     & A0(2,2),A1(2,2)
      parameter(minus=1,mplus=2)
      real(dp), parameter :: mp(nf)=(/-1d0,+1d0,-1d0,+1d0,-1d0/)

      scheme='tH-V'
c--- note that scheme = "tH-V" is also a valid choice in this case
c--- (finite terms relating the two schemes are included in, e.g. AWWjeta.f)

      msq(:,:)=0d0

c Change the momenta to DKS notation: e-(3) nu~(4) e+(5) nu(6)
      p(1,:)=pin(1,:)
      p(2,:)=pin(2,:)
      p(3,:)=pin(5,:)
      p(4,:)=pin(6,:)
      p(5,:)=pin(4,:)
      p(6,:)=pin(3,:)
      p(7,:)=pin(7,:)

c      fac=gw**4
      fac=esq**2/abs(zxw)**2
      fac1=two*gsq*cf

      call spinoru(7,p,za,zb)

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
      call amp_WWg_v(p,2,7,3,4,5,6,1,2,gub0a,gub0b,gub0sr,gub1a,gub1b,gub1ax,gub1sr)
      call amp_WWg_v(p,7,2,3,4,5,6,1,2, gu0a, gu0b, gu0sr, gu1a, gu1b, gu1ax, gu1sr)
      call amp_WWg_v(p,2,7,6,5,4,3,1,1,gdb0a,gdb0b,gdb0sr,gdb1a,gdb1b,gdb1ax,gdb1sr)
      call amp_WWg_v(p,7,2,6,5,4,3,1,1, gd0a, gd0b, gd0sr, gd1a, gd1b, gd1ax, gd1sr)
      call amp_wwg_v_qloop(2,7,3,4,5,6,1,gub1q,gu1q,gdb1q,gd1q,gub1ax,gu1ax,gdb1ax,gd1ax)

c--- quark-gluon and antiquark-gluon
      call amp_WWg_v(p,1,7,3,4,5,6,2,2,ubg0a,ubg0b,ubg0sr,ubg1a,ubg1b,ubg1ax,ubg1sr)
      call amp_WWg_v(p,7,1,3,4,5,6,2,2, ug0a, ug0b, ug0sr, ug1a, ug1b, ug1ax, ug1sr)
      call amp_WWg_v(p,1,7,6,5,4,3,2,1,dbg0a,dbg0b,dbg0sr,dbg1a,dbg1b,dbg1ax,dbg1sr)
      call amp_WWg_v(p,7,1,6,5,4,3,2,1, dg0a, dg0b, dg0sr, dg1a, dg1b, dg1ax, dg1sr)
      call amp_wwg_v_qloop(1,7,3,4,5,6,2,ubg1q,ug1q,dbg1q,dg1q,ubg1ax,ug1ax,dbg1ax,dg1ax)

c--- quark-antiquark and antiquark-quark
      call amp_WWg_v(p,1,2,3,4,5,6,7,2,ubu0a,ubu0b,ubu0sr,ubu1a,ubu1b,ubu1ax,ubu1sr)
      call amp_WWg_v(p,1,2,6,5,4,3,7,1,dbd0a,dbd0b,dbd0sr,dbd1a,dbd1b,dbd1ax,dbd1sr)
      endif
      if ((useblha == 0) .or. (blhatype == 2)) then
      call amp_WWg_v(p,2,1,3,4,5,6,7,2,uub0a,uub0b,uub0sr,uub1a,uub1b,uub1ax,uub1sr)
      endif
      if ((useblha == 0) .or. (blhatype == 1)) then
      call amp_WWg_v(p,2,1,6,5,4,3,7,1,ddb0a,ddb0b,ddb0sr,ddb1a,ddb1b,ddb1ax,ddb1sr)
      endif
      call amp_wwg_v_qloop(1,2,3,4,5,6,7,ubu1q,uub1q,dbd1q,ddb1q,ubu1ax,uub1ax,dbd1ax,ddb1ax)

c--- fill matrix elements (no bottom-initiated contribution)
      do j=-4,4
      do k=-4,4

c-- skip gluon-gluon case
      if((j  ==  0).and.(k  ==  0)) goto 19

c-- skip non-diagonal quark flavors except gluon
      if((j  /=  0 .and. k  /=  0) .and. (j  /=  -k)) goto 19

      if (j  ==  0 .or. k  ==  0) then
c--- appropriate values for quark-gluon scattering
        jk=j+k
        ave=xn*aveqg
      else
c--- appropriate values for quark-antiquark scattering
        jk=max(j,k)
        ave=xn*aveqq
      endif

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
        A0(polq,polg)=fac*(
     &   cprop*(amp0a(polq,polg)+amp0b(polq,polg)*c(polq,tjk))+zxw*amp0sr(polq,polg))
        A1(polq,polg)=fac*(
     &   cprop*(amp1a(polq,polg)+amp1b(polq,polg)*c(polq,tjk)+amp1ax(polq,polg)*tau(tjk))+zxw*amp1sr(polq,polg))
      enddo
      enddo


c--- fill matrix elements with appropriate averaging and overall factors
c--- (check of LO)
c      msq(j,k)=fac1*ave*real(
c     & +A0(1,1)*conjg(A0(1,1))+A0(1,2)*conjg(A0(1,2))
c     & +A0(2,1)*conjg(A0(2,1))+A0(2,2)*conjg(A0(2,2)),dp)
c--- (virtual)
      msq(j,k)=ason2pi*fac1*ave*real(
     & +A1(1,1)*conjg(A0(1,1))+A1(1,2)*conjg(A0(1,2))
     & +A1(2,1)*conjg(A0(2,1))+A1(2,2)*conjg(A0(2,2)),dp)

   19 continue

      enddo
      enddo

      return
      end

