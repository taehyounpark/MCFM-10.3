!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_zz_v(p,msqv)
      implicit none
      include 'types.f'
c----Author R.K.Ellis December 1998
c----modified by JMC to include supplementary diagrams February 1999
c--- calculate the virtual matrix element squared
c----and subtraction terms for ZZ production
c----modified by RKE (following suggestion of GZ)
c----to correct supplementary diagrams April 2011
c----NB: we also include virtual photons
c    averaged over initial colours and spins
c    u(-p1)+dbar(-p2)-->\mu^-(p4)+\mu^+(p5)+e^-(p6)+e^+(p7)
c    Notation to allow room for p3 --- gluon emission.
c----No statistical factor of 1/2 included.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'masses.f'
      include 'zcouple.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'scheme.f'
      include 'ewcharge.f'
      include 'srdiags.f'
      include 'interference.f'
      include 'zcouple_cms.f'
      include 'blha.f'

      real(dp):: msqv(-nf:nf,-nf:nf),
     & p(mxpart,4),qdks(mxpart,4),virt,
     & fac,facnlo,ave,rescale1,rescale2

      complex(dp):: v2(2),v1(2)
      complex(dp):: qqb(2,2,2),qbq(2,2,2),lqqb(2,2,2),lqbq(2,2,2)
      complex(dp):: qqb1(2,2,2),qbq1(2,2,2),qqb2(2,2,2),qbq2(2,2,2)
      complex(dp):: a6trees,a6loops
      complex(dp):: aqqb,aqbq,bqqb,bqbq,Vpole,Vpole12,suppl
      complex(dp):: prop12,prop34,prop56

      integer:: j,k,polq,pol1,pol2,ii,nmax
      integer,parameter::i4(2)=(/4,5/),i5(2)=(/5,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      complex(dp):: aqqb_SAVE(-nf:nf,-nf:nf,2,2,2)
      complex(dp):: bqqb_SAVE(-nf:nf,-nf:nf,2,2,2)
      complex(dp):: aqbq_SAVE(-nf:nf,-nf:nf,2,2,2)
      complex(dp):: bqbq_SAVE(-nf:nf,-nf:nf,2,2,2)

c Process with ZZ/WW interference
      if (blhatype == 1) then
        call qqb_vv_v(p,msqv)
        return
      endif

      scheme='dred'

      fac=-four*abs(zesq)**2
      ave=aveqq*xn*vsymfact
      facnlo=ason2pi*cf

      v1(1)=zl1
      v1(2)=zr1
      v2(1)=zl2
      v2(2)=zr2

c----setup factor to avoid summing over too many neutrinos
c----if coupling enters twice
c      if (q1 == zip) then
c      rescale1=one/sqrt(three)
c      else
      rescale1=one
c      endif
c      if (q2 == zip) then
c      rescale2=one/sqrt(three)
c      else
      rescale2=one
c      endif

c--set msq=0 to initalize
      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=zip
      aqqb_SAVE(j,k,:,:,:)=zip
      aqbq_SAVE(j,k,:,:,:)=zip
      bqqb_SAVE(j,k,:,:,:)=zip
      bqbq_SAVE(j,k,:,:,:)=zip

      enddo
      enddo

      if (interference) then
         nmax=2
      else
         nmax=1
      endif


c----Change the momenta to DKS notation
c   We have --- q(-p1)+qbar(-p2)-->l(p3)+lbar(p4) + l'(p5)+lbar'(p6)
c   DKS have--- q(q2) +qbar(q1) -->mu^-(q3)+mu^+(q4)+e^-(q6)+e^+(q5)
         do j=1,4
            qdks(1,j)=p(1,j)
            qdks(2,j)=p(2,j)
            qdks(3,j)=p(3,j)
            qdks(4,j)=p(4,j)
            qdks(5,j)=p(6,j)
            qdks(6,j)=p(5,j)
         enddo

      call spinoru(6,qdks,za,zb)


      do ii=1,nmax


c--   s returned from sprod (common block) is 2*dot product
c--   calculate propagators
c----No Baur-Zeppenfeld
      prop12=s(1,2)/cplx2(s(1,2)-zmass**2,zmass*zwidth)
      prop34=s(3,i4(ii))/cplx2(s(3,i4(ii))-zmass**2,zmass*zwidth)
      prop56=s(i5(ii),6)/cplx2(s(i5(ii),6)-zmass**2,zmass*zwidth)

c-- here the labels correspond to the polarizations of the
c-- quark, lepton 4 and lepton 6 respectively

      qbq(1,1,1)=A6trees(1,2,6,i5(ii),i4(ii),3,za,zb)
      qbq(1,1,2)=A6trees(1,2,6,i5(ii),3,i4(ii),za,zb)
      qbq(1,2,1)=A6trees(1,2,i5(ii),6,i4(ii),3,za,zb)
      qbq(1,2,2)=A6trees(1,2,i5(ii),6,3,i4(ii),za,zb)

      qqb(1,1,1)=A6trees(2,1,6,i5(ii),i4(ii),3,za,zb)
      qqb(1,1,2)=A6trees(2,1,6,i5(ii),3,i4(ii),za,zb)
      qqb(1,2,1)=A6trees(2,1,i5(ii),6,i4(ii),3,za,zb)
      qqb(1,2,2)=A6trees(2,1,i5(ii),6,3,i4(ii),za,zb)

      lqbq(1,1,1)=A6loops(1,2,6,i5(ii),i4(ii),3,za,zb)
      lqbq(1,1,2)=A6loops(1,2,6,i5(ii),3,i4(ii),za,zb)
      lqbq(1,2,1)=A6loops(1,2,i5(ii),6,i4(ii),3,za,zb)
      lqbq(1,2,2)=A6loops(1,2,i5(ii),6,3,i4(ii),za,zb)

      lqqb(1,1,1)=A6loops(2,1,6,i5(ii),i4(ii),3,za,zb)
      lqqb(1,1,2)=A6loops(2,1,6,i5(ii),3,i4(ii),za,zb)
      lqqb(1,2,1)=A6loops(2,1,i5(ii),6,i4(ii),3,za,zb)
      lqqb(1,2,2)=A6loops(2,1,i5(ii),6,3,i4(ii),za,zb)

      if (srdiags) then
c---for supplementary diagrams.
      qbq1(1,1,1)=+A6trees(3,i4(ii),1,2,i5(ii),6,za,zb)
      qbq2(1,1,1)=+A6trees(6,i5(ii),1,2,i4(ii),3,za,zb)
      qbq1(1,1,2)=-A6trees(i4(ii),3,1,2,i5(ii),6,za,zb)
      qbq2(1,1,2)=+A6trees(6,i5(ii),1,2,3,i4(ii),za,zb)
      qbq1(1,2,1)=+A6trees(3,i4(ii),1,2,6,i5(ii),za,zb)
      qbq2(1,2,1)=-A6trees(i5(ii),6,1,2,i4(ii),3,za,zb)
      qbq1(1,2,2)=-A6trees(i4(ii),3,1,2,6,i5(ii),za,zb)
      qbq2(1,2,2)=-A6trees(i5(ii),6,1,2,3,i4(ii),za,zb)

      qqb1(1,1,1)=-A6trees(3,i4(ii),2,1,i5(ii),6,za,zb)
      qqb2(1,1,1)=-A6trees(6,i5(ii),2,1,i4(ii),3,za,zb)
      qqb1(1,1,2)=+A6trees(i4(ii),3,2,1,i5(ii),6,za,zb)
      qqb2(1,1,2)=-A6trees(6,i5(ii),2,1,3,i4(ii),za,zb)
      qqb1(1,2,1)=-A6trees(3,i4(ii),2,1,6,i5(ii),za,zb)
      qqb2(1,2,1)=+A6trees(i5(ii),6,2,1,i4(ii),3,za,zb)
      qqb1(1,2,2)=+A6trees(i4(ii),3,2,1,6,i5(ii),za,zb)
      qqb2(1,2,2)=+A6trees(i5(ii),6,2,1,3,i4(ii),za,zb)
c---loop diagrams just tree*Vpole since they're all triangle-type
      Vpole12=Vpole(s(1,2))
      endif

      do j=1,2
      do k=1,2
      qbq(2,j,k)=-qqb(1,j,k)
      qqb(2,j,k)=-qbq(1,j,k)
      lqbq(2,j,k)=-lqqb(1,j,k)
      lqqb(2,j,k)=-lqbq(1,j,k)
      qbq1(2,j,k)=-qqb1(1,j,k)
      qqb1(2,j,k)=-qbq1(1,j,k)
      qbq2(2,j,k)=-qqb2(1,j,k)
      qqb2(2,j,k)=-qbq2(1,j,k)
      enddo
      enddo

c for interference calculation, only generate (34) resonance, so instead of
c |A1+A2|^2, accummulate |A1+A2|^2+A1^2-A2^2

c----calculate over limited flavor range -2:2
      do j=-2,2
      k=-j
      virt=zip
      if (j==0) go to 20

      if ((j > 0).and.(k < 0)) then
      do polq=1,2
      do pol1=1,2
      do pol2=1,2
      if     (polq == 1) then
       aqqb=(prop56*v2(pol1)*zl(j)+q2*q(j))
     &     *(prop34*v1(pol2)*zl(j)+q1*q(j))* qqb(polq,pol1,pol2)
       bqqb=(prop56*v2(pol1)*zl(j)+q2*q(j))
     &     *(prop34*v1(pol2)*zl(j)+q1*q(j))*lqqb(polq,pol1,pol2)
         if (srdiags) then
         suppl=-(
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zl(j)+q1*q(j))*qqb1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zl(j)+q2*q(j))*qqb2(polq,pol1,pol2))

         aqqb=aqqb+suppl
         bqqb=bqqb+suppl*Vpole12
         endif
      elseif (polq == 2) then
       aqqb=(prop56*v2(pol1)*zr(j)+q2*q(j))
     &     *(prop34*v1(pol2)*zr(j)+q1*q(j))* qqb(polq,pol1,pol2)
       bqqb=(prop56*v2(pol1)*zr(j)+q2*q(j))
     &     *(prop34*v1(pol2)*zr(j)+q1*q(j))*lqqb(polq,pol1,pol2)
         if (srdiags) then
         suppl=-(
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zr(j)+q1*q(j))*qqb1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zr(j)+q2*q(j))*qqb2(polq,pol1,pol2))

         aqqb=aqqb+suppl
         bqqb=bqqb+suppl*Vpole12
         endif
      endif

      aqqb=FAC*aqqb
      bqqb=FAC*bqqb

c      virt=virt+facnlo*ave*two*real(conjg(aqqb)*bqqb)
c      !interference terms
c      if ((ii==1).and.(interference)) then
c         aqqb_SAVE(j,k,polq,pol1,pol2)=aqqb
c         bqqb_SAVE(j,k,polq,pol1,pol2)=bqqb
c      elseif (ii==2) then
c         if (pol1==pol2) then
c      virt=virt-two*facnlo*ave*real(
c     &           conjg(aqqb)*bqqb_SAVE(j,k,polq,pol1,pol2)
c     &          +conjg(bqqb)*aqqb_SAVE(j,k,polq,pol1,pol2))
c         endif
c      endif


      if (interference .eqv. .false.) then
c--- normal case
        virt=virt+facnlo*ave*two*real(conjg(aqqb)*bqqb)
      else
c--- with interference:
c---    1st pass --> store result
c---    2nd pass --> fill msq
        if (ii == 1) then
          aqqb_SAVE(j,k,polq,pol1,pol2)=aqqb
          bqqb_SAVE(j,k,polq,pol1,pol2)=bqqb
        else
          if (useblha == 1) then
             virt=virt+facnlo*ave*two*real(
     &          conjg(aqqb_SAVE(j,k,polq,pol1,pol2))*bqqb_SAVE(j,k,polq,pol1,pol2)
     &         +conjg(aqqb)*bqqb)
             if (pol1 == pol2) then
             virt=virt-facnlo*ave*two*real(
     &          conjg(aqqb_SAVE(j,k,polq,pol1,pol2))*bqqb
     &         +conjg(aqqb)*bqqb_SAVE(j,k,polq,pol1,pol2))
             endif
          else
            if (pol1 == pol2) then
              virt=virt+facnlo*ave*two*real(
     &         conjg(aqqb_SAVE(j,k,polq,pol1,pol2)-aqqb)*(bqqb_SAVE(j,k,polq,pol1,pol2)-bqqb)
     &        +conjg(aqqb_SAVE(j,k,polq,pol1,pol2))*bqqb_SAVE(j,k,polq,pol1,pol2)
     &        -conjg(aqqb)*bqqb)
            else
              virt=virt+facnlo*ave*two*real(
     &         2._dp*conjg(aqqb_SAVE(j,k,polq,pol1,pol2))*bqqb_SAVE(j,k,polq,pol1,pol2))
            endif
          endif
        endif
      endif


      enddo
      enddo
      enddo

      elseif ((j < 0).and.(k > 0)) then

      do polq=1,2
      do pol1=1,2
      do pol2=1,2
      if     (polq == 1) then
       aqbq=(prop56*v2(pol1)*zl(k)+q2*q(k))
     &     *(prop34*v1(pol2)*zl(k)+q1*q(k))* qbq(polq,pol1,pol2)
       bqbq=(prop56*v2(pol1)*zl(k)+q2*q(k))
     &     *(prop34*v1(pol2)*zl(k)+q1*q(k))*lqbq(polq,pol1,pol2)
         if (srdiags) then
         suppl=
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zl(k)+q1*q(k))*qbq1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zl(k)+q2*q(k))*qbq2(polq,pol1,pol2)
         aqbq=aqbq+suppl
         bqbq=bqbq+suppl*Vpole12
         endif
      elseif (polq == 2) then
       aqbq=(prop56*v2(pol1)*zr(k)+q2*q(k))
     &     *(prop34*v1(pol2)*zr(k)+q1*q(k))* qbq(polq,pol1,pol2)
       bqbq=(prop56*v2(pol1)*zr(k)+q2*q(k))
     &     *(prop34*v1(pol2)*zr(k)+q1*q(k))*lqbq(polq,pol1,pol2)
         if (srdiags) then
         suppl=
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zr(k)+q1*q(k))*qbq1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zr(k)+q2*q(k))*qbq2(polq,pol1,pol2)
         aqbq=aqbq+suppl
         bqbq=bqbq+suppl*Vpole12
         endif
      endif

      aqbq=FAC*aqbq
      bqbq=FAC*bqbq

      if (interference .eqv. .false.) then
c--- normal case
        virt=virt+facnlo*ave*two*real(conjg(aqbq)*bqbq)
      else
c--- with interference:
c---    1st pass --> store result
c---    2nd pass --> fill msq
        if (ii == 1) then
          aqbq_SAVE(j,k,polq,pol1,pol2)=aqbq
          bqbq_SAVE(j,k,polq,pol1,pol2)=bqbq
        else
          if (pol1 == pol2) then
            virt=virt+facnlo*ave*two*real(
     &       conjg(aqbq_SAVE(j,k,polq,pol1,pol2)-aqbq)*(bqbq_SAVE(j,k,polq,pol1,pol2)-bqbq)
     &      +conjg(aqbq_SAVE(j,k,polq,pol1,pol2))*bqbq_SAVE(j,k,polq,pol1,pol2)
     &      -conjg(aqbq)*bqbq)
          else
            virt=virt+facnlo*ave*two*real(
     &       2._dp*conjg(aqbq_SAVE(j,k,polq,pol1,pol2))*bqbq_SAVE(j,k,polq,pol1,pol2))
          endif
        endif
      endif

      enddo
      enddo
      enddo

      endif

      msqv(j,k)=msqv(j,k)+virt

 20   continue
      enddo
      enddo !ii

c----extend to full flavor range
      do j=-nf,nf
      k=-j

      if (j==0) go to 21

      msqv(j,k)=msqv(jkswitch(j),jkswitch(k))

 21      continue
      enddo

      return
      end
