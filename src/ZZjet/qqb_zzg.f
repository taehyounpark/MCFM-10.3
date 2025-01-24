!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_zzg(p,msq)
c Note: LO routine for ZZ+jet process
c Differs from real radiation in ZZ process in the following ways:
c       (1) Does not include the possibility of anomalous couplings
c       (2) Is computed in the complex mass scheme
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
      include 'srdiags.f'
      include 'interference.f'
      include 'pchoice.f'
      integer:: jk,hq,h34,h56,hg,ii,nmax
      real(dp):: fac,fac1,q34,q56,s127
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: ave,rescale1,rescale2,oprat
      complex(dp):: v34(2),v56(2),LR(2,2)
      complex(dp):: aq12(2,2,2,2),aq34(2,2,2,2),aq56(2,2,2,2)
      complex(dp):: qa12(2,2,2,2),qa34(2,2,2,2),qa56(2,2,2,2)
      complex(dp):: gq12(2,2,2,2),gq34(2,2,2,2),gq56(2,2,2,2)
      complex(dp):: ag12(2,2,2,2),ag34(2,2,2,2),ag56(2,2,2,2)
      complex(dp):: ga12(2,2,2,2),ga34(2,2,2,2),ga56(2,2,2,2)
      complex(dp):: qg12(2,2,2,2),qg34(2,2,2,2),qg56(2,2,2,2)
      complex(dp):: Uncrossed(-nf:nf,-nf:nf,2,2,2,2)
      complex(dp):: amp,prop34,prop56,prop127
      integer,parameter::i4(2)=(/4,6/),i6(2)=(/6,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

c-set Uncrossed array to zero
c      call writeout(p)
      msq(:,:)=zip
      Uncrossed(:,:,:,:,:,:)=czip

c      write(6,*) 'g,gsq',sqrt(gsq),gsq
c      write(6,*) 'esq',esq
c      write(6,*) esq/(4.d0*pi),(0.3082739d0)**2/(4d0*pi)
c      write(6,*) (4.d0*pi)/esq
c      write(6,*) 'sqrt(esq)',sqrt(esq)
c      write(6,*) 'gau',sqrt(esq)*2d0/3d0
c      write(6,*) 'gad',-sqrt(esq)/3d0
c      write(6,*) 'gzl,gzr',-sqrt(esq)*l1,-sqrt(esq)*r1
c      write(6,*) 'gzul,gzur',-sqrt(esq)*L(2),-sqrt(esq)*R(2)
c      write(6,*) 'gzdl,gzdr',-sqrt(esq)*L(1),-sqrt(esq)*R(1)
c      write(6,*) 'qqb_zzg: srdiags',srdiags
c      write(6,*)
c      pause
      fac=-four*esq**2
      fac1=two*cf*xn*gsq

      v34(1)=zl1
      v34(2)=zr1
      q34=q1
      v56(1)=zl2
      v56(2)=zr2
      q56=q2

      LR(1,1)=zL(1)
      LR(1,2)=zL(2)
      LR(2,1)=zR(1)
      LR(2,2)=zR(2)

c      write(6,*) 'fac',fac
c      write(6,*) 'q34,v34',q34,v34
c      write(6,*) 'q56,v56',q56,v56

c----setup factor to avoid summing over too many neutrinos
c----if coupling enters twice
c      if (q34 == zip) then
c      rescale1=one/sqrt(three)
c      else
c      rescale1=one
c      endif
c      if (q56 == zip) then
c      rescale2=one/sqrt(three)
c      else
c      rescale2=one
c      endif

      rescale1=one
      rescale2=one

c--   s returned from sprodx (common block) is 2*dot product
      call spinoru(7,p,za,zb)

      if (interference) then
         nmax=2
      else
         nmax=1
      endif

      do ii=1,nmax

c--   calculate propagators
      s127=s(1,2)+s(1,7)+s(2,7)
      prop127=s127/cplx2(s127-zmass**2,zmass*zwidth)
      prop34=s(3,i4(ii))/cplx2(s(3,i4(ii))-zmass**2,zmass*zwidth)
      prop56=s(5,i6(ii))/cplx2(s(5,i6(ii))-zmass**2,zmass*zwidth)
c      write(6,*) 's127',s127
c      write(6,*) 'prop34',prop34
c      write(6,*) 'prop56',prop56
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

      if (jk == 0) goto 19

      do hq=1,2
      do h34=1,2
      do h56=1,2
      do hg=1,2

      amp=zip

c---case qbar-q
      if    ((j < 0).and.(k > 0)) then
      amp=(prop56*v56(h56)*lr(hq,k)+q56*q(k))
     & *(prop34*v34(h34)*lr(hq,k)+q34*q(k))*aq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,k)+q34*q(k))*aq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,k)+q56*q(k))*aq34(hq,h34,h56,hg)
      endif

c---case q-qbar
      elseif((j > 0).and.(k < 0)) then
      amp=(prop56*v56(h56)*lr(hq,j)+q56*q(j))
     & *(prop34*v34(h34)*lr(hq,j)+q34*q(j))*qa12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,j)+q34*q(j))*qa56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,j)+q56*q(j))*qa34(hq,h34,h56,hg)
      endif

c---case qbar-g
      elseif((j < 0).and.(k == 0)) then
      amp=(prop56*v56(h56)*lr(hq,-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*lr(hq,-jk)+q34*q(-jk))*ag12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,-jk)+q34*q(-jk))*ag56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,-jk)+q56*q(-jk))*ag34(hq,h34,h56,hg)
      endif

c---case g-qbar
      elseif((k < 0).and.(j == 0)) then
      amp=(prop56*v56(h56)*lr(hq,-jk)+q56*q(-jk))
     & *(prop34*v34(h34)*lr(hq,-jk)+q34*q(-jk))*ga12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,-jk)+q34*q(-jk))*ga56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,-jk)+q56*q(-jk))*ga34(hq,h34,h56,hg)
      endif

c---case q-g
      elseif((j > 0).and.(k == 0)) then
      amp=(prop56*v56(h56)*lr(hq,+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*lr(hq,+jk)+q34*q(+jk))*qg12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,+jk)+q34*q(+jk))*qg56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,+jk)+q56*q(+jk))*qg34(hq,h34,h56,hg)
      endif

      elseif((k > 0).and.(j == 0)) then
c---case g-q
      amp=(prop56*v56(h56)*lr(hq,+jk)+q56*q(+jk))
     & *(prop34*v34(h34)*lr(hq,+jk)+q34*q(+jk))*gq12(hq,h34,h56,hg)
      if (srdiags) then
      amp=amp
     & +(prop56*v34(h34)*v56(h56)+q34*q56)*rescale1
     & *(prop127*v34(h34)*lr(hq,+jk)+q34*q(+jk))*gq56(hq,h34,h56,hg)
     & +(prop34*v34(h34)*v56(h56)+q34*q56)*rescale2
     & *(prop127*v56(h56)*lr(hq,+jk)+q56*q(+jk))*gq34(hq,h34,h56,hg)
      endif

      endif

      amp=amp*fac

c      msq(j,k)=msq(j,k)+fac1*ave*abs(amp)**2
c      !set-up interference terms
c      if ((interference).and.(ii==1)) then
c      Uncrossed(j,k,hq,h34,h56,hg)=amp
c      elseif (ii==2) then
c      if (h34==h56) then
c      msq(j,k)=msq(j,k)
c     &    -fac1*two*ave*real(conjg(amp)*Uncrossed(j,k,hq,h34,h56,hg))
c      endif
c      endif


      if (interference .eqv. .false.) then
c--- normal case
        msq(j,k)=msq(j,k)+fac1*ave*abs(amp)**2
      else
c--- with interference:
c---    1st pass --> store result
c---    2nd pass --> fill msq
        if (ii == 1) then
          Uncrossed(j,k,hq,h34,h56,hg)=amp
        else
          if (h34 == h56) then
            oprat=one
     &           -two*real(conjg(amp)*Uncrossed(j,k,hq,h34,h56,hg))
     &            /(abs(amp)**2+abs(Uncrossed(j,k,hq,h34,h56,hg))**2)
          else
            oprat=one
          endif
          if (bw34_56) then
            msq(j,k)=msq(j,k)
     &        +fac1*ave*two*abs(Uncrossed(j,k,hq,h34,h56,hg))**2*oprat
          else
            msq(j,k)=msq(j,k)+fac1*ave*two*abs(amp)**2*oprat
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
      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) cycle
      msq(j,k)=msq(jkswitch(j),jkswitch(k))
      enddo
      enddo

      return
      end




