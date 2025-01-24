!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_wz_g(P,msq)
      implicit none
      include 'types.f'
c---Author John Campbell Fri Feb 19 11:06:08 CST 1999
c---Modified to include supplementary diagrams by JC on Feb 24
c---Matrix element squared averaged over initial colors and spins
c---  averaged(summed) over initial(final) colours and spins
c For nwz=+1
c     u(-p1)+dbar(-p2)-->mu^-(p5)+mu^+(p6)+n(p3)+e^+(p4)+g(p7)
c For nwz=-1
c     d(-p1)+ubar(-p2)-->mu^-(p5)+mu^+(p6)+e^-(p3)+nbar(p4)+g(p7)
c---
c   for the moment --- radiation only from initial line

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'masses.f'
      include 'ckm.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      include 'zerowidth.f'
      include 'ewcharge.f'
      include 'anomcoup.f'
      include 'xanomcoup.f'
      include 'nwz.f'
      include 'plabel.f'
      include 'pchoice.f'
      include 'zcouple_cms.f'
      include 'ewinput.f'
      include 'interference.f'
      integer:: polg,polz,minus,mplus,jp,kp,ii,nmax,i3(2),i4(2),i5(2),i6(2)
      real(dp):: FAC,FACM,FAC1
      complex(dp):: prop12,prop34,prop56
      real(dp):: P(mxpart,4),qdks(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: ave,s127,wwflag,xfac
      complex(dp)::
     &  qu_qb(10,2,2),qu_gg(10,2,2),gg_qb(10,2,2),
     &  qb_qu(10,2,2),qb_gg(10,2,2),gg_qu(10,2,2),
     &  props,propw,propz,cprop,A(2,2),
     &  A_save(-nf:nf,-nf:nf,2,2)
      complex(dp):: v2(2),cl1,cl2,en1,en2
      complex(dp):: ZgLR(nf,2),c1(2),c2(2),cotw
      parameter(minus=1,mplus=2)
      FAC=-2._dp*abs(zesq/zxw)*abs(zesq)
      FAC1=two*gsq*cf
      if ((nwz==1) .or. (nwz == -1)) then
      FACM=nwz*FAC
      else
      write(6,*) 'nwz  /=  +1 or -1 in qqb_wz_g.f'
      stop
      endif
      if     (nwz==-1) then
        cl1=1._dp
        cl2=0._dp
        en1=zle
        en2=zln
      elseif (nwz==+1) then
        cl1=0._dp
        cl2=1._dp
        en1=zln
        en2=zle
      endif
      v2(1)=zl1
      v2(2)=zr1

c--- wwflag=1 for most cases, indicating presence of diagram with 2 W's
      wwflag=1._dp
c--- but for Z -> bbbar this diagram contains |V_tb|**2 which we take 0
      if (plabel(5) == 'bq') then
        wwflag=0._dp
      endif

c-- if Z -> neutrinos, we need to switch c1 and c2
      if (plabel(5) == 'nl') then
        cl1=1._dp-cl1
        cl2=1._dp-cl2
      endif

      do jp=-nf,nf
      do kp=-nf,nf
      msq(jp,kp)=0._dp
      enddo
      enddo

c----Change the momenta to DKS notation
c   We have --- d(-p1)+ubar(-p2)-->nu(p3)+e^+(p4)+mu^-(p5)+mu^+(p6)+g(p7)
c   DKS have--- u( q2)+dbar( q1)-->nu(q3)+e^+(q4)+mu^-(q6)+mu^+(q5)+g(p7)

      do jp=1,4
      qdks(1,jp)=p(1,jp)
      qdks(2,jp)=p(2,jp)
      qdks(3,jp)=p(3,jp)
      qdks(4,jp)=p(4,jp)
      qdks(5,jp)=p(6,jp)
      qdks(6,jp)=p(5,jp)
      qdks(7,jp)=p(7,jp)
      enddo

      call spinoru(7,qdks,za,zb)

      i3(1)=3
      i4(1)=4
      i5(1)=5
      i6(1)=6

      if (interference) then
         nmax=2
         if (nwz == +1) then
           i3(2)=3
           i4(2)=5
           i5(2)=4
           i6(2)=6
         else
           i3(2)=6
           i4(2)=4
           i5(2)=5
           i6(2)=3
         endif
      else
         nmax=1
      endif

      do ii=1,nmax

c--   s returned from sprodx (common block) is 2*dot product

c--   calculate propagators
      cotw=sqrt((one-zxw)/zxw)
      s127=s(1,2)+s(1,7)+s(2,7)
      if     ((zerowidth  .eqv. .true.) .or. (ewscheme > 3)) then
      prop12=s127/cplx2(s127-wmass**2,wmass*wwidth)
      prop34=s(i3(ii),i4(ii))/cplx2(s(i3(ii),i4(ii))-wmass**2,wmass*wwidth)
      prop56=s(i5(ii),i6(ii))/cplx2(s(i5(ii),i6(ii))-zmass**2,zmass*zwidth)
      cprop=cone
      elseif (zerowidth .neqv. .true.) then
      prop12=cplx1(s127/(s127-wmass**2))
      prop34=cplx1(s(i3(ii),i4(ii))/(s(i3(ii),i4(ii))-wmass**2))
      prop56=cplx1(s(i5(ii),i6(ii))/(s(i5(ii),i6(ii))-zmass**2))
      props=(s127-wmass**2)/cplx2(s127-wmass**2,wmass*wwidth)
      propw=(s(i3(ii),i4(ii))-wmass**2)/cplx2(s(i3(ii),i4(ii))-wmass**2,wmass*wwidth)
      propz=(s(i5(ii),i6(ii))-zmass**2)/cplx2(s(i5(ii),i6(ii))-zmass**2,zmass*zwidth)
      cprop=props*propw*propz
      endif

c--- DEBUG to compare with Madgraph
c      prop12=s127/cplx2(s127-wmass**2,wmass*wwidth)
c      prop34=s(3,4)/cplx2(s(3,4)-wmass**2,wmass*wwidth)
c      prop56=s(5,6)/cplx2(s(5,6)-zmass**2,zmass*zwidth)
c      cprop=cone
c--- DEBUG to compare with Madgraph

c--- apply a dipole form factor to anomalous couplings (only if tevscale > 0)
      if (tevscale > 0._dp) then
        xfac=1._dp/(1._dp+s127/(tevscale*1d3)**2)**2
      else
        xfac=1._dp
      endif
      xdelg1_z=xfac*delg1_z
      xdelg1_g=xfac*delg1_g
      xdelk_z=xfac*delk_z
      xdelk_g=xfac*delk_g
      xlambda_z=xfac*lambda_z
      xlambda_g=xfac*lambda_g

c---case dbar-u
      call wzamps(1,2,i3(ii),i4(ii),i5(ii),i6(ii),7,za,zb,qb_qu)
c---case u-dbar
      call wzamps(2,1,i3(ii),i4(ii),i5(ii),i6(ii),7,za,zb,qu_qb)

c---case g-u
      call wzamps(7,2,i3(ii),i4(ii),i5(ii),i6(ii),1,za,zb,gg_qu)
c---case u-g
      call wzamps(7,1,i3(ii),i4(ii),i5(ii),i6(ii),2,za,zb,qu_gg)

c---case dbar-g
      call wzamps(1,7,i3(ii),i4(ii),i5(ii),i6(ii),2,za,zb,qb_gg)
c---case g-dbar
      call wzamps(2,7,i3(ii),i4(ii),i5(ii),i6(ii),1,za,zb,gg_qb)

c---set up left/right handed couplings for both Z and gamma*
c---note that the second label corresponds to the helicity
c---of the LEPTON coupling v2, NOT the quarks (all L)
      do j=1,nf
        ZgLR(j,minus)=zL(j)*v2(1)*prop56+Q(j)*q1
        ZgLR(j,mplus)=zL(j)*v2(2)*prop56+Q(j)*q1
      enddo

      do polz=1,2
      if(nwz==1) then
        c1(polz)=ZgLR(2,polz)
        c2(polz)=ZgLR(1,polz)
      else
        c1(polz)=ZgLR(1,polz)
        c2(polz)=ZgLR(2,polz)
      endif
      enddo

      do j=-nf,nf
      if (((j==+1).or.(j==+3).or.(j==+5)
     & .or.(j==-2).or.(j==-4)) .and. (nwz == +1))
     & go to 20
      if (((j==-1).or.(j==-3).or.(j==-5)
     & .or.(j==+2).or.(j==+4)) .and. (nwz == -1))
     & go to 20
      do k=-nf,nf
      if (((k==+1).or.(k==+3).or.(k==+5)
     & .or.(k==-2).or.(k==-4)) .and. (nwz == +1))
     & go to 19
      if (((k==-1).or.(k==-3).or.(k==-5)
     & .or.(k==+2).or.(k==+4)) .and. (nwz == -1))
     & go to 19

      if     ((j > 0) .and. (k < 0)) then

c---case u-db
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((ZgLR(+j,polz)*qu_qb(2,polg,polz)
     &                  +ZgLR(-k,polz)*qu_qb(3,polg,polz))*FAC
     &                  +(cotw*v2(polz)*prop56*qu_qb(1,polg,polz)
     &                                     +q1*qu_qb(10,polg,polz))
     &                    *prop12*FACM)*prop34
     &            +FAC*((en1*v2(polz)*prop56+q1*(-1._dp)*cl1)
     &                   *prop12*qu_qb(5,polg,polz)
     &                 +(en2*v2(polz)*prop56+q1*(-1._dp)*cl2)
     &                   *prop12*qu_qb(4,polg,polz)
     &          +wwflag*0.5_dp*prop34*prop12/zxw*qu_qb(6,polg,polz)*cl1
     &          +wwflag*0.5_dp*prop34*prop12/zxw*qu_qb(7,polg,polz)*cl2)

c          A(polg,polz)=((L(+j)*qu_qb(2,polg,polz)
c     &                  +L(-k)*qu_qb(3,polg,polz))*FAC
c     &                  +cotw*prop12*qu_qb(1,polg,polz)*FACM)
c     &                 *prop34*prop56*v2(polz)
          enddo
          enddo
          ave=xn*aveqq*Vsq(j,k)

      elseif ((j < 0) .and. (k > 0)) then


c---case db-u
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((ZgLR(+k,polz)*qb_qu(2,polg,polz)
     &                  +ZgLR(-j,polz)*qb_qu(3,polg,polz))*FAC
     &                  +(cotw*v2(polz)*prop56*qb_qu(1,polg,polz)
     &                                     +q1*qb_qu(10,polg,polz))
     &                    *prop12*FACM)*prop34
     &            +FAC*((en1*v2(polz)*prop56+q1*(-1._dp)*cl1)
     &                   *prop12*qb_qu(5,polg,polz)
     &                 +(en2*v2(polz)*prop56+q1*(-1._dp)*cl2)
     &                   *prop12*qb_qu(4,polg,polz)
     &          +wwflag*0.5_dp*prop34*prop12/zxw*qb_qu(6,polg,polz)*cl1
     &          +wwflag*0.5_dp*prop34*prop12/zxw*qb_qu(7,polg,polz)*cl2)

c          A(polg,polz)=((L(+k)*qb_qu(2,polg,polz)
c     &                  +L(-j)*qb_qu(3,polg,polz))*FAC
c     &                  +cotw*prop12*qb_qu(1,polg,polz)*FACM)
c     &                 *prop34*prop56*v2(polz)
          enddo
          enddo
          ave=xn*aveqq*Vsq(j,k)

      elseif ((j > 0) .and. (k == 0)) then
c---case u-g
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((c1(polz)*qu_gg(2,polg,polz)
     &                  +c2(polz)*qu_gg(3,polg,polz))*FAC
     &                  +(cotw*v2(polz)*prop56*qu_gg(1,polg,polz)
     &                                     +q1*qu_gg(10,polg,polz))
     &                    *prop12*FACM)*prop34
     &            +FAC*((en1*v2(polz)*prop56+q1*(-1._dp)*cl1)
     &                   *prop12*qu_gg(5,polg,polz)
     &                 +(en2*v2(polz)*prop56+q1*(-1._dp)*cl2)
     &                   *prop12*qu_gg(4,polg,polz)
     &          +wwflag*0.5_dp*prop34*prop12/zxw*qu_gg(6,polg,polz)*cl1
     &          +wwflag*0.5_dp*prop34*prop12/zxw*qu_gg(7,polg,polz)*cl2)

c          A(polg,polz)=((c1*qu_gg(2,polg,polz)
c     &                  +c2*qu_gg(3,polg,polz))*FAC
c     &                  +cotw*prop12*qu_gg(1,polg,polz)*FACM)
c     &                 *prop34*prop56*v2(polz)
          enddo
          enddo
          ave=xn*aveqg*Vsum(j)
      elseif ((j < 0) .and. (k == 0)) then
c---case db-g
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((c1(polz)*qb_gg(2,polg,polz)
     &                  +c2(polz)*qb_gg(3,polg,polz))*FAC
     &                  +(cotw*v2(polz)*prop56*qb_gg(1,polg,polz)
     &                                     +q1*qb_gg(10,polg,polz))
     &                    *prop12*FACM)*prop34
     &            +FAC*((en1*v2(polz)*prop56+q1*(-1._dp)*cl1)
     &                   *prop12*qb_gg(5,polg,polz)
     &                 +(en2*v2(polz)*prop56+q1*(-1._dp)*cl2)
     &                   *prop12*qb_gg(4,polg,polz)
     &          +wwflag*0.5_dp*prop34*prop12/zxw*qb_gg(6,polg,polz)*cl1
     &          +wwflag*0.5_dp*prop34*prop12/zxw*qb_gg(7,polg,polz)*cl2)

c          A(polg,polz)=((c1*qb_gg(2,polg,polz)
c     &                  +c2*qb_gg(3,polg,polz))*FAC
c     &                  +cotw*prop12*qb_gg(1,polg,polz)*FACM)
c     &                 *prop34*prop56*v2(polz)
          enddo
          enddo

          ave=xn*aveqg*Vsum(j)
      elseif ((j == 0) .and. (k > 0)) then
c---case g-u
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((c1(polz)*gg_qu(2,polg,polz)
     &                  +c2(polz)*gg_qu(3,polg,polz))*FAC
     &                  +(cotw*v2(polz)*prop56*gg_qu(1,polg,polz)
     &                                     +q1*gg_qu(10,polg,polz))
     &                    *prop12*FACM)*prop34
     &            +FAC*((en1*v2(polz)*prop56+q1*(-1._dp)*cl1)
     &                   *prop12*gg_qu(5,polg,polz)
     &                 +(en2*v2(polz)*prop56+q1*(-1._dp)*cl2)
     &                   *prop12*gg_qu(4,polg,polz)
     &          +wwflag*0.5_dp*prop34*prop12/zxw*gg_qu(6,polg,polz)*cl1
     &          +wwflag*0.5_dp*prop34*prop12/zxw*gg_qu(7,polg,polz)*cl2)

c          A(polg,polz)=((c1(polz)*gg_qu(2,polg,polz)
c     &                  +c2(polz)*gg_qu(3,polg,polz))*FAC
c     &                  +cotw*prop12*gg_qu(1,polg,polz)*FACM)
c     &                 *prop34*prop56*v2(polz)
          enddo
          enddo
          ave=xn*aveqg*Vsum(k)
      elseif ((j == 0) .and. (k < 0)) then
c---case g-db
          do polg=1,2
          do polz=1,2
          A(polg,polz)=((c1(polz)*gg_qb(2,polg,polz)
     &                  +c2(polz)*gg_qb(3,polg,polz))*FAC
     &                  +(cotw*v2(polz)*prop56*gg_qb(1,polg,polz)
     &                                     +q1*gg_qb(10,polg,polz))
     &                    *prop12*FACM)*prop34
     &            +FAC*((en1*v2(polz)*prop56+q1*(-1._dp)*cl1)
     &                   *prop12*gg_qb(5,polg,polz)
     &                 +(en2*v2(polz)*prop56+q1*(-1._dp)*cl2)
     &                   *prop12*gg_qb(4,polg,polz)
     &          +wwflag*0.5_dp*prop34*prop12/zxw*gg_qb(6,polg,polz)*cl1
     &          +wwflag*0.5_dp*prop34*prop12/zxw*gg_qb(7,polg,polz)*cl2)

c          A(polg,polz)=((c1*gg_qb(2,polg,polz)
c     &                  +c2*gg_qb(3,polg,polz))*FAC
c     &                  +cotw*prop12*gg_qb(1,polg,polz)*FACM)
c     &                 *prop34*prop56*v2(polz)

          enddo
          enddo

          ave=xn*aveqg*Vsum(k)

      else
          ave=0._dp
      endif

      if (ave>0._dp) then

        if (interference .eqv. .false.) then
c--- normal case
        msq(j,k)=FAC1*ave*vsymfact*abs(cprop)**2
     &            *(abs(A(mplus,minus))**2+abs(A(minus,minus))**2
     &             +abs(A(mplus,mplus))**2+abs(A(minus,mplus))**2)
        else
c--- with interference:
c---    1st pass --> store result
c---    2nd pass --> fill msq
          if (ii == 1) then
            A_save(j,k,:,:)=A(:,:)
          else
            msq(j,k)=FAC1*ave*vsymfact*abs(cprop)**2
     &            *(abs(A_save(j,k,mplus,minus)-A(mplus,minus))**2
     &             +abs(A_save(j,k,mplus,minus))**2-abs(A(mplus,minus))**2
     &             +abs(A_save(j,k,minus,minus)-A(minus,minus))**2
     &             +abs(A_save(j,k,minus,minus))**2-abs(A(minus,minus))**2
     &             +two*abs(A_save(j,k,mplus,mplus))**2
     &             +two*abs(A_save(j,k,minus,mplus))**2)
          endif
        endif

      endif

 19   continue
      enddo
 20   continue
      enddo

      enddo   ! end ii loop

      return
      end







