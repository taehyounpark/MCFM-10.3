!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_zz(p,msq)
      implicit none
      include 'types.f'
c----Author R.K.Ellis December 1998
c----modified by JMC to include supplementary diagrams February 1999
c----modified by RKE (following suggestion of GZ)
c----to correct supplementary diagrams April 2011
c----Matrix element for ZZ production
c----NB: we also include virtual photons
c    in the notation of DKS
c    averaged over initial colours and spins
c    u(-p1)+dbar(-p2)-->e^-(p3) + e^+(p4)   + \mu^-(p5)+ \mu^+(p6)
c    q(-p1)+qbar(-p2)-->l'(p3)  + lbar'(p4) + l(p5)    + lbar(p6)
c    with Z-leptons couplings l1 for (5,6) and l2 for(3,4)
c          and lepton charges q2 for (5,6) and q1 for (3,4)
c----No statistical factor of 1/2 included.
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'masses.f'
      include 'zcouple_cms.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'ewcharge.f'
      include 'srdiags.f'
      include 'interference.f'

      real(dp):: msq(-nf:nf,-nf:nf),p(mxpart,4),qdks(mxpart,4),
     & ave,rescale1,rescale2,q1,q2

      complex(dp):: qqb(2,2,2),qbq(2,2,2),q_qb,qb_q,v1(2),v2(2)
      complex(dp):: qqb1(2,2,2),qbq1(2,2,2),qqb2(2,2,2),qbq2(2,2,2)
      complex(dp):: a6trees,prop12,prop34,prop56
      complex(dp):: Uncrossed(-nf:nf,-nf:nf,2,2,2)

      real(dp):: FAC
      integer:: j,k,polq,pol1,pol2
      integer:: ii,nmax
      integer,parameter::i4(2)=(/4,5/),i5(2)=(/5,4/),
     & jkswitch(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

c     vsymfact=symmetry factor
      fac=-four*abs(zesq)**2
      ave=aveqq*xn*vsymfact
      q1=rq1
      q2=rq2
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
      msq(j,k)=zip
      Uncrossed(j,k,:,:,:)=zip
      enddo
      enddo

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

      if (interference) then
         nmax=2
      else
         nmax=1
      endif

      do ii=1,nmax

c--   s returned from sprod (common block) is 2*dot product
c--   calculate propagators
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
      endif

      qbq(2,:,:)=-qqb(1,:,:)
      qqb(2,:,:)=-qbq(1,:,:)
      qbq1(2,:,:)=-qqb1(1,:,:)
      qqb1(2,:,:)=-qbq1(1,:,:)
      qbq2(2,:,:)=-qqb2(1,:,:)
      qqb2(2,:,:)=-qbq2(1,:,:)

c for interference calculation, only generate (34) resonance, so instead of
c |A1+A2|^2, accummulate |A1+A2|^2+A1^2-A2^2

      do j=-2,2
      k=-j

      if (j==0) go to 20

      if ((j > 0).and.(k < 0)) then
      do polq=1,2
      do pol1=1,2
      do pol2=1,2
      if     (polq == 1) then
       q_qb=(prop56*v2(pol1)*zl(j)+q2*q(j))
     &     *(prop34*v1(pol2)*zl(j)+q1*q(j))*qqb(polq,pol1,pol2)
        if (srdiags) then
          q_qb=q_qb-(
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zl(j)+q1*q(j))*qqb1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zl(j)+q2*q(j))*qqb2(polq,pol1,pol2))

        endif
      elseif (polq == 2) then
       q_qb=(prop56*v2(pol1)*zr(j)+q2*q(j))
     &     *(prop34*v1(pol2)*zr(j)+q1*q(j))*qqb(polq,pol1,pol2)
        if (srdiags) then
         q_qb=q_qb-(
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zr(j)+q1*q(j))*qqb1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zr(j)+q2*q(j))*qqb2(polq,pol1,pol2))

        endif
      endif
      q_qb=FAC*q_qb

      if (interference .eqv. .false.) then
c--- normal case
        msq(j,k)=msq(j,k)+ave*abs(q_qb)**2
      else
c--- with interference:
c---    1st pass --> store result
c---    2nd pass --> fill msq
        if (ii == 1) then
          Uncrossed(j,k,polq,pol1,pol2)=q_qb
        else
          if (pol1 == pol2) then
            msq(j,k)=msq(j,k)+ave*(abs(Uncrossed(j,k,polq,pol1,pol2)-q_qb)**2
     &                            +abs(Uncrossed(j,k,polq,pol1,pol2))**2-abs(q_qb)**2)
          else
            msq(j,k)=msq(j,k)+ave*(2._dp*abs(Uncrossed(j,k,polq,pol1,pol2))**2)
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
       qb_q=(prop56*v2(pol1)*zl(k)+q2*q(k))
     &     *(prop34*v1(pol2)*zl(k)+q1*q(k))*qbq(polq,pol1,pol2)
        if (srdiags) then
         qb_q=qb_q
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zl(k)+q1*q(k))*qbq1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zl(k)+q2*q(k))*qbq2(polq,pol1,pol2)

        endif
      elseif (polq == 2) then
       qb_q=(prop56*v2(pol1)*zr(k)+q2*q(k))
     &     *(prop34*v1(pol2)*zr(k)+q1*q(k))*qbq(polq,pol1,pol2)
        if (srdiags) then
         qb_q=qb_q
     &       +(prop56*v1(pol2)*v2(pol1)+q1*q2)*rescale1
     &       *(prop12*v1(pol2)*zr(k)+q1*q(k))*qbq1(polq,pol1,pol2)
     &       +(prop34*v1(pol2)*v2(pol1)+q1*q2)*rescale2
     &       *(prop12*v2(pol1)*zr(k)+q2*q(k))*qbq2(polq,pol1,pol2)

        endif
      endif
      qb_q=FAC*qb_q

      if (interference .eqv. .false.) then
c--- normal case
        msq(j,k)=msq(j,k)+ave*abs(qb_q)**2
      else
c--- with interference:
c---    1st pass --> store result
c---    2nd pass --> fill msq
        if (ii == 1) then
          Uncrossed(j,k,polq,pol1,pol2)=qb_q
        else
          if (pol1 == pol2) then
            msq(j,k)=msq(j,k)+ave*(abs(Uncrossed(j,k,polq,pol1,pol2)-qb_q)**2
     &                            +abs(Uncrossed(j,k,polq,pol1,pol2))**2-abs(qb_q)**2)
          else
            msq(j,k)=msq(j,k)+ave*(2._dp*abs(Uncrossed(j,k,polq,pol1,pol2))**2)
          endif
        endif
      endif

      enddo !end pol2 loop
      enddo !end pol1 loop
      enddo !end polq loop

      endif

 20   continue
      enddo !end j-loop
      enddo   ! end ii loop

c----extend matrix element to full flavor range
      do j=-nf,nf
      k=-j
      if (j==0) go to 21
      msq(j,k)=msq(jkswitch(j),jkswitch(k))
 21   continue
      enddo

      return
      end
