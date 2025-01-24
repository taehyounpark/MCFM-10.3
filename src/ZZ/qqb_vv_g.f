!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_vv_g(p,msq)
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
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'srdiags.f'
      include 'pchoice.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      include 'zcouple_cms.f'
      integer:: jk,hq,h34,h56,hg,ii,nmax
      real(dp):: fac,fac1,q34,q56,s127
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf),qdks(mxpart,4)
      real(dp):: ave,rescale1,rescale2
      complex(dp):: aq12(2,2,2,2),aq34(2,2,2,2),aq56(2,2,2,2)
      complex(dp):: qa12(2,2,2,2),qa34(2,2,2,2),qa56(2,2,2,2)
      complex(dp):: gq12(2,2,2,2),gq34(2,2,2,2),gq56(2,2,2,2)
      complex(dp):: ag12(2,2,2,2),ag34(2,2,2,2),ag56(2,2,2,2)
      complex(dp):: ga12(2,2,2,2),ga34(2,2,2,2),ga56(2,2,2,2)
      complex(dp):: qg12(2,2,2,2),qg34(2,2,2,2),qg56(2,2,2,2)
      complex(dp):: amp,prop34,prop56,prop127,v34(2),v56(2),wwamp(5)
      complex(dp):: azz(2,2,2,2),aww(2,2,2,2)
      complex(dp):: ct(2,2),cs_z(2,2),cs_g(2,2),cgamz(2,2),cz(2,2),facww
      complex(dp)::u_ub(5,2,2),d_db(5,2,2),ub_u(5,2,2),db_d(5,2,2),
     &             u_g(5,2,2), d_g(5,2,2), g_ub(5,2,2),g_db(5,2,2),
     &             ub_g(5,2,2),db_g(5,2,2),g_u(5,2,2),g_d(5,2,2),prop12
      integer:: tjk,minus,mplus,jtype,polg,polq
      real(dp), parameter :: mp(nf)=(/-1._dp,+1._dp,-1._dp,+1._dp,-1._dp/)
      parameter(minus=1,mplus=2)
      integer,parameter::i4(2)=(/4,6/),i6(2)=(/6,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      fac=-four*abs(zesq)**2
      fac1=two*cf*gsq

c Fix couplings by removing factors of 1/sqrt(3) (and restore at end)
      zl2=zln
      zr2=zrn
      rescale1=one
      rescale2=one

      v34(1)=zl1
      v34(2)=zr1
      q34=q1
      v56(1)=zl2
      v56(2)=zr2
      q56=q2

c--set msq=0 to initalize
      msq(:,:)=zip

c--   s returned from sprodx (common block) is 2*dot product
      call spinoru(7,p,za,zb)

      nmax=1

      do ii=1,nmax

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7)
      prop127=s127/cplx2(s127-zmass**2,zmass*zwidth)
      prop34=s(3,i4(ii))/cplx2(s(3,i4(ii))-zmass**2,zmass*zwidth)
      prop56=s(5,i6(ii))/cplx2(s(5,i6(ii))-zmass**2,zmass*zwidth)

c--- Amplitude returned with arguments (hq,h34,h56,h7)
c---case qbar-q
      call zzgamp(1,2,3,i4(ii),5,i6(ii),7,za,zb,aq12,aq34,aq56)
c---case q-qbar
      call zzgamp(2,1,3,i4(ii),5,i6(ii),7,za,zb,qa12,qa34,qa56)
c---case qbar-g
      call zzgamp(1,7,3,i4(ii),5,i6(ii),2,za,zb,ag12,ag34,ag56)
c---case q-g
      call zzgamp(7,1,3,i4(ii),5,i6(ii),2,za,zb,qg12,qg34,qg56)
c---case g-q
      call zzgamp(7,2,3,i4(ii),5,i6(ii),1,za,zb,gq12,gq34,gq56)
c---case g-qbar
      call zzgamp(2,7,3,i4(ii),5,i6(ii),1,za,zb,ga12,ga34,ga56)

      enddo

c!!!!!!!!!!!!!!!!!! WW BLOCK
c!!!!!!!!!!!!!!!!!! WW BLOCK
c!!!!!!!!!!!!!!!!!! WW BLOCK

      facww=-abs(zesq)**2/zxw**2

c----Change the momenta to DKS notation
c   swapped possibility if we want to swap momenta for hadronic case
c   We have --- f(p1) + f'(p2)-->mu^-(p3)+nubar(p4)+e^+(p6)+nu(p5)+g(p7)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)+g(p7)
c----
c   or normal configuration
c   We have --- f(p1) + f'(p2)-->mu^-(p5)+nubar(p6)+e^+(p4)+nu(p3)+g(p7)
c   DKS have--- ubar(q1)+u(q2)-->mu^-(q3)+nubar(q4)+e^+(q5)+nu(q6)+g(p7)

c----all other cases
      do j=1,4
      qdks(1,j)=p(1,j)
      qdks(2,j)=p(2,j)
      qdks(3,j)=p(3,j) ! interchanged 3<->5 wrt. qqb_ww_g.f
      qdks(4,j)=p(6,j)
      qdks(5,j)=p(4,j)
      qdks(6,j)=p(5,j) ! interchanged 3<->5 wrt. qqb_ww_g.f
      qdks(7,j)=p(7,j)
      enddo

c--   s returned from sprodx (common block) is 2*dot product
      call spinoru(7,qdks,za,zb)

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7)
      prop12=s127/cplx2(s127-zmass**2,zmass*zwidth)

c-- couplings according to 3.4 and 3.6
      do j=1,2
      ct(minus,j)=1._dp
      ct(mplus,j)=0._dp
      cs_z(minus,j)=+mp(j)*zl(j)*zsin2w*prop12
      cs_z(mplus,j)=-mp(j)*2._dp*Q(j)*zxw*prop12
      cs_g(minus,j)=+mp(j)*2._dp*Q(j)*zxw
      cs_g(mplus,j)=+mp(j)*2._dp*Q(j)*zxw
      cz(minus,j)=two*zxw*zln*zL(j)*prop12
      cz(mplus,j)=two*zxw*zln*zR(j)*prop12
      cgamz(minus,j)=two*zxw*(-Q(j)+zle*zL(j)*prop12)
      cgamz(mplus,j)=two*zxw*(-Q(j)+zle*zR(j)*prop12)
      enddo
c--
c      l(j)=(tau(j)-two*Q(j)*xw)/sin2w ; r(j)=(-two*Q(j)*xw)/sin2w
c      le=(-1._dp-two*(-1._dp)*xw)/sin2w ; re=(-two*(-1._dp)*xw)/sin2w
c      ln=(+1._dp-two*(+0._dp)*xw)/sin2w ; rn=0._dp
c---

      xdelg1_z=delg1_z
      xdelg1_g=delg1_g
      xdelk_z=delk_z
      xdelk_g=delk_g
      xlambda_z=lambda_z
      xlambda_g=lambda_g

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

c!!!!!!!!!!!!!!!!!! END OF WW BLOCK
c!!!!!!!!!!!!!!!!!! END OF WW BLOCK
c!!!!!!!!!!!!!!!!!! END OF WW BLOCK

      aww(:,:,:,:)=czip
      azz(:,:,:,:)=czip

c---calculate over a limited flavour range (-2:2)
      do j=-2,2
      do k=-2,2
c-- skip gluon-gluon case
      if((j == 0).and.(k == 0)) goto 19
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) goto 19
      jk=max(j,k)
      ave=xn*aveqq
      if (j == 0 .or. k == 0) then
          jk=j+k
          ave=xn*aveqg
      endif

      do hq=1,2
      do h34=1,2
      do h56=1,2
      do hg=1,2

      amp=zip

c---case qbar-q
      if    ((j < 0).and.(k > 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*zL(k)+q56*q(k))
     & *(prop34*v34(h34)*zL(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zL(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zL(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*zR(k)+q56*q(k))
     & *(prop34*v34(h34)*zR(k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zR(k)+q34*q(k))*aq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zR(k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif
      endif

c---case q-qbar
      elseif((j > 0).and.(k < 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*zL(j)+q56*q(j))
     & *(prop34*v34(h34)*zL(j)+q34*q(j))*qa12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zL(j)+q34*q(j))*qa56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zL(j)+q56*q(j))*qa34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*zR(j)+q56*q(j))
     & *(prop34*v34(h34)*zR(j)+q34*q(j))*qa12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zR(j)+q34*q(j))*qa56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zR(j)+q56*q(j))*qa34(hq,h34,h56,hg)
      endif
      endif

c---case qbar-g
      elseif((j < 0).and.(k == 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*zL(-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*zL(-jk)+q34*q(-jk))*ag12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zL(-jk)+q34*q(-jk))*ag56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zL(-jk)+q56*q(-jk))*ag34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*zR(-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*zR(-jk)+q34*q(-jk))*ag12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zR(-jk)+q34*q(-jk))*ag56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zR(-jk)+q56*q(-jk))*ag34(hq,h34,h56,hg)
      endif
      endif

c---case g-qbar
      elseif((k < 0).and.(j == 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*zL(-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*zL(-jk)+q34*q(-jk))*ga12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zL(-jk)+q34*q(-jk))*ga56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zL(-jk)+q56*q(-jk))*ga34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*zR(-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*zR(-jk)+q34*q(-jk))*ga12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zR(-jk)+q34*q(-jk))*ga56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zR(-jk)+q56*q(-jk))*ga34(hq,h34,h56,hg)
      endif
      endif

c---case q-g
      elseif((j > 0).and.(k == 0)) then
      if (hq == 1) then
      amp=(prop56*v56(h56)*zL(+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*zL(+jk)+q34*q(+jk))*qg12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zL(+jk)+q34*q(+jk))*qg56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zL(+jk)+q56*q(+jk))*qg34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*zR(+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*zR(+jk)+q34*q(+jk))*qg12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zR(+jk)+q34*q(+jk))*qg56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zR(+jk)+q56*q(+jk))*qg34(hq,h34,h56,hg)
      endif
      endif

      elseif((k > 0).and.(j == 0)) then
c---case g-q
      if (hq == 1) then
      amp=(prop56*v56(h56)*zL(+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*zL(+jk)+q34*q(+jk))*gq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zL(+jk)+q34*q(+jk))*gq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zL(+jk)+q56*q(+jk))*gq34(hq,h34,h56,hg)
      endif
      elseif (hq == 2) then
      amp=(prop56*v56(h56)*zR(+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*zR(+jk)+q34*q(+jk))*gq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*zR(+jk)+q34*q(+jk))*gq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*zR(+jk)+q56*q(+jk))*gq34(hq,h34,h56,hg)
      endif
      endif

      endif

      azz(hq,h34,h56,hg)=amp*fac

      enddo  ! endloop hg
      enddo  ! endloop h56
      enddo  ! endloop h34
      enddo  ! endloop hq

c!!!!!!!!!!!!!!!!!! WW BLOCK
c!!!!!!!!!!!!!!!!!! WW BLOCK
c!!!!!!!!!!!!!!!!!! WW BLOCK


      do polg=1,2
      do polq=1,2
c---sum is over diagram type t,s(Z),e,n,s(photon)
      do jtype=1,5
          if    (j < 0 .and. tau(jk) == -1._dp .and. k  /=  0) then
            wwamp(jtype)=db_d(jtype,polg,polq)
          elseif(j < 0 .and. tau(jk) ==  1._dp .and. k  /=  0) then
            wwamp(jtype)=ub_u(jtype,polg,polq)
          elseif(j > 0 .and. tau(jk) == -1._dp .and. k  /=  0) then
            wwamp(jtype)=d_db(jtype,polg,polq)
          elseif(j > 0 .and. tau(jk) ==  1._dp .and. k  /=  0) then
            wwamp(jtype)=u_ub(jtype,polg,polq)
          elseif(j == 0 .and. tau(jk) ==  1._dp .and. jk > 0) then
            wwamp(jtype)=g_u(jtype,polg,polq)
          elseif(j == 0 .and. tau(jk) == -1._dp .and. jk > 0) then
            wwamp(jtype)=g_d(jtype,polg,polq)
          elseif(j == 0 .and. tau(jk) == -1._dp .and. jk < 0) then
            wwamp(jtype)=g_ub(jtype,polg,polq)
          elseif(j == 0 .and. tau(jk) ==  1._dp .and. jk < 0) then
            wwamp(jtype)=g_db(jtype,polg,polq)
          elseif(k == 0 .and. tau(jk) ==  1._dp .and. jk > 0) then
            wwamp(jtype)=u_g(jtype,polg,polq)
          elseif(k == 0 .and. tau(jk) == -1._dp .and. jk > 0) then
            wwamp(jtype)=d_g(jtype,polg,polq)
          elseif(k == 0 .and. tau(jk) == -1._dp .and. jk < 0) then
            wwamp(jtype)=ub_g(jtype,polg,polq)
          elseif(k == 0 .and. tau(jk) ==  1._dp .and. jk < 0) then
            wwamp(jtype)=db_g(jtype,polg,polq)
          endif
      enddo

c---tjk is equal to 2 (u,c) or 1 (d,s,b)
      tjk=2-mod(abs(jk),2)

      aww(polq,1,1,polg)=facww
     & *(ct(polq,tjk)*wwamp(1)+cs_z(polq,tjk)*wwamp(2)+cs_g(polq,tjk)*wwamp(5)
     &  +cz(polq,tjk)*wwamp(3)+cgamz(polq,tjk)*wwamp(4))

      enddo
      enddo

c!!!!!!!!!!!!!!!!!! END OF WW BLOCK
c!!!!!!!!!!!!!!!!!! END OF WW BLOCK
c!!!!!!!!!!!!!!!!!! END OF WW BLOCK

      do hq=1,2
      do h34=1,2
      do h56=1,2
      do hg=1,2
c      msq(j,k)=msq(j,k)+fac1*ave*abs(aww(hq,h34,h56,hg))**2
      msq(j,k)=msq(j,k)+fac1*ave*abs(azz(hq,h34,h56,hg)+aww(hq,h34,h56,hg))**2
      enddo
      enddo
      enddo
      enddo

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

c      msq(:,:)=msq(:,:)*3d0 ! artificially inflate by factor of three

      return
      end




