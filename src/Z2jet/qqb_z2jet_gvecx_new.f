!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_z2jet_gvecx_new(p,n,in,msqv,msqv_cs,ppmsqvx)
      implicit none
      include 'types.f'

c-----Nov 10 99 --- checked that gives the right matrix element
c     when summed over polarizations.

c----Matrix element for Z+2jet production
c----averaged over initial colours and spins
c    line 6 contracted with the vector n(mu)
c     q(-p1)+qbar(-p2)--> g(p5)+ g(p6)+Z(f(p3)+af(p4))

c---- Note that this routine is the counterpart of qqb_z2jet_gvec
c---- and should be used when we are calculating the 4Q contribution.
c---- This routine is necessary because it returns the contribution
c---- for the correct ordering of particles in the final state - for
c---- example, gg -> qqb as well as gg -> qbq

c--- Note: this differs from the "old" qqb_z2jet_gvecx.f by the fact
c---  that the array msqvx with 11^4 = 14641 elements has been
c---  replaced by the array ppmsqvx with 80 elements

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'ewcouple.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'pp.f'
      include 'zcouple_cms.f'

c ip is the label of the emitter
c in is the label of the contracted line
      integer:: j,k,pq,pl,nquark,in,icol
      real(dp):: fac,n(4),ggbit1,ggbit2
      complex(dp):: prop,zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),
     & qqbZgg_cs(0:2,2,2),qbqZgg_cs(0:2,2,2),qgZqg_cs(0:2,2,2),
     & gqZqg_cs(0:2,2,2),qbgZqbg_cs(0:2,2,2),gqbZqbg_cs(0:2,2,2),
     & ggZqbq_cs(0:2,2,2),ggZqqb_cs(0:2,2,2),
     & qqbZgg(2,2),qbqZgg(2,2),qgZqg(2,2),gqZqg(2,2),
     & qbgZqbg(2,2),gqbZqbg(2,2),ggZqbq(2,2),
     & qgZgq(2,2),qgZgq_cs(0:2,2,2),qbgZgqb(2,2),qbgZgqb_cs(0:2,2,2),
     & gqZgq(2,2),gqZgq_cs(0:2,2,2),gqbZgqb(2,2),gqbZgqb_cs(0:2,2,2)
      real(dp):: msqv_cs(0:2,-nf:nf,-nf:nf)
      real(dp):: ppmsqvx(ppmax)
      integer,parameter::
     & jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/),
     & kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      include 'cplx.h'

      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      do icol=0,2
        msqv_cs(icol,j,k)=0._dp
      enddo
      enddo
      enddo

      do icol=0,2
      do pq=1,2
      do pl=1,2
      qqbZgg_cs(icol,pq,pl) =0._dp
      gqZqg_cs(icol,pq,pl)  =0._dp
      qgZqg_cs(icol,pq,pl)  =0._dp
      qbqZgg_cs(icol,pq,pl) =0._dp
      gqbZqbg_cs(icol,pq,pl)=0._dp
      qbgZqbg_cs(icol,pq,pl)=0._dp
      ggZqbq_cs(icol,pq,pl) =0._dp
      ggZqqb_cs(icol,pq,pl) =0._dp
      enddo
      enddo
      enddo

      call spinoru(6,p,za,zb)
c   zab=<i-|k|j-> zba=<i+|k|j+> where k is an arbitrary 4-vector
c---Conventions of Bern, Dixon, Kosower, Weinzierl,
c---ie, za(i,j)*zb(j,i)=s(i,j)
      call spinork(6,p,zab,zba,n)

      if (in == 1) then
cargument 1-4 represent (1) incoming quark line
c                       (2) incoming anti-quark line
c                       (3) outgoing gluon line
c                       (4) outgoing gluon line contracted with n
           call z2jetsqn(5,6,2,1,p,n,za,zb,zab,zba,ggZqbq)
           call storezcsv(ggZqbq_cs)
           call z2jetsqn(2,5,6,1,p,n,za,zb,zab,zba,gqZqg)
           call storezcsv(gqZqg_cs)
           call z2jetsqn(5,2,6,1,p,n,za,zb,zab,zba,gqbZqbg)
           call storezcsv(gqbZqbg_cs)
c--- extra term for gg -> qqb
           call z2jetsqn(6,5,2,1,p,n,za,zb,zab,zba,ggZqbq)
           call storezcsv(ggZqqb_cs)
c--- extra term for gq -> gq (individual color structures not needed)
           call z2jetsqn(2,6,5,1,p,n,za,zb,zab,zba,gqZgq)
c--- extra term for gqb -> gqb (individual color structures not needed)
           call z2jetsqn(6,2,5,1,p,n,za,zb,zab,zba,gqbZgqb)
      elseif (in == 2)  then
           call z2jetsqn(5,6,1,2,p,n,za,zb,zab,zba,ggZqbq)
           call storezcsv_swapcol(ggZqbq_cs)
           call z2jetsqn(5,1,6,2,p,n,za,zb,zab,zba,qbgZqbg)
           call storezcsv(qbgZqbg_cs)
           call z2jetsqn(1,5,6,2,p,n,za,zb,zab,zba,qgZqg)
           call storezcsv(qgZqg_cs)
c--- extra term for gg -> qqb
           call z2jetsqn(6,5,1,2,p,n,za,zb,zab,zba,ggZqbq)
           call storezcsv(ggZqqb_cs)
c--- extra term for qg -> gq (individual color structures not needed)
           call z2jetsqn(1,6,5,2,p,n,za,zb,zab,zba,qgZgq)
c--- extra term for qbg -> gqb (individual color structures not needed)
           call z2jetsqn(6,1,5,2,p,n,za,zb,zab,zba,qbgZgqb)
      elseif (in == 5) then
          call z2jetsqn(1,2,6,5,p,n,za,zb,zab,zba,qqbZgg)
          call storezcsv(qqbZgg_cs)
          call z2jetsqn(2,1,6,5,p,n,za,zb,zab,zba,qbqZgg)
          call storezcsv(qbqZgg_cs)
c The variable names have been changed in the following calls to reflect the fact
c that the gluon is actually in position 5, but we must copy arrays at the end in
c order to ensure that msqv_cs is always filled
          call z2jetsqn(1,6,2,5,p,n,za,zb,zab,zba,qgZgq)
          call storezcsv_swapcol(qgZgq_cs)
          call z2jetsqn(2,6,1,5,p,n,za,zb,zab,zba,gqZgq)
          call storezcsv_swapcol(gqZgq_cs)
          call z2jetsqn(6,1,2,5,p,n,za,zb,zab,zba,qbgZgqb)
          call storezcsv_swapcol(qbgZgqb_cs)
          call z2jetsqn(6,2,1,5,p,n,za,zb,zab,zba,gqbZgqb)
          call storezcsv_swapcol(gqbZgqb_cs)
          qgZqg_cs(:,:,:)=qgZgq_cs(:,:,:)
          gqZqg_cs(:,:,:)=gqZgq_cs(:,:,:)
          qbgZqbg_cs(:,:,:)=qbgZgqb_cs(:,:,:)
          gqbZqbg_cs(:,:,:)=gqbZgqb_cs(:,:,:)
      elseif (in == 6) then
          call z2jetsqn(1,2,5,6,p,n,za,zb,zab,zba,qqbZgg)
          call storezcsv_swapcol(qqbZgg_cs)
          call z2jetsqn(2,1,5,6,p,n,za,zb,zab,zba,qbqZgg)
          call storezcsv_swapcol(qbqZgg_cs)
          call z2jetsqn(1,5,2,6,p,n,za,zb,zab,zba,qgZqg)
          call storezcsv_swapcol(qgZqg_cs)
          call z2jetsqn(2,5,1,6,p,n,za,zb,zab,zba,gqZqg)
          call storezcsv_swapcol(gqZqg_cs)
          call z2jetsqn(5,1,2,6,p,n,za,zb,zab,zba,qbgZqbg)
          call storezcsv_swapcol(qbgZqbg_cs)
          call z2jetsqn(5,2,1,6,p,n,za,zb,zab,zba,gqbZqbg)
          call storezcsv_swapcol(gqbZqbg_cs)
      endif

      prop=s(3,4)/cplx2((s(3,4)-zmass**2),zmass*zwidth)
      fac=v*xn/four*(esq*gsq)**2*two

      do icol=0,2
      do pq=1,2
      do pl=1,2
      qqbZgg_cs(icol,pq,pl) =half*aveqq*fac*qqbZgg_cs(icol,pq,pl)
      qbqZgg_cs(icol,pq,pl) =half*aveqq*fac*qbqZgg_cs(icol,pq,pl)
      gqZqg_cs(icol,pq,pl)  =aveqg*fac*gqZqg_cs(icol,pq,pl)
      qgZqg_cs(icol,pq,pl)  =aveqg*fac*qgZqg_cs(icol,pq,pl)
      gqbZqbg_cs(icol,pq,pl)=aveqg*fac*gqbZqbg_cs(icol,pq,pl)
      qbgZqbg_cs(icol,pq,pl)=aveqg*fac*qbgZqbg_cs(icol,pq,pl)
      ggZqbq_cs(icol,pq,pl) =avegg*fac*ggZqbq_cs(icol,pq,pl)
      ggZqqb_cs(icol,pq,pl) =avegg*fac*ggZqqb_cs(icol,pq,pl)
      enddo
      enddo
      enddo

c Additional swapped matrix elements that are not broken into color structures
      qgZgq(:,:)=aveqg*fac*qgZgq(:,:)
      gqZgq(:,:)=aveqg*fac*gqZgq(:,:)
      qbgZgqb(:,:)=aveqg*fac*qbgZgqb(:,:)
      gqbZgqb(:,:)=aveqg*fac*gqbZgqb(:,:)

      ppmsqvx(:)=zip

      do j=-nf,nf
      do k=-nf,nf
      if( (j  /=  0) .and. (k  /=  0) .and. (j  /=  -k)) goto 19

      if     ((j == 0) .and. (k == 0)) then

          do icol=0,2
            msqv_cs(icol,j,k)=0._dp
            do nquark=1,2
            ggbit1=
     &      +abs(Q(nquark)*q1+zL(nquark)*zl1*prop)**2*ggZqbq_cs(icol,1,1)
     &      +abs(Q(nquark)*q1+zR(nquark)*zr1*prop)**2*ggZqbq_cs(icol,2,2)
     &      +abs(Q(nquark)*q1+zL(nquark)*zr1*prop)**2*ggZqbq_cs(icol,1,2)
     &      +abs(Q(nquark)*q1+zR(nquark)*zl1*prop)**2*ggZqbq_cs(icol,2,1)
            ggbit2=
     &      +abs(Q(nquark)*q1+zL(nquark)*zl1*prop)**2*ggZqqb_cs(icol,1,1)
     &      +abs(Q(nquark)*q1+zR(nquark)*zr1*prop)**2*ggZqqb_cs(icol,2,2)
     &      +abs(Q(nquark)*q1+zL(nquark)*zr1*prop)**2*ggZqqb_cs(icol,1,2)
     &      +abs(Q(nquark)*q1+zR(nquark)*zl1*prop)**2*ggZqqb_cs(icol,2,1)
            if (icol == 0) then
              ppmsqvx(pp(j,k,-nquark,nquark))=ggbit1
              ppmsqvx(pp(j,k,nquark,-nquark))=ggbit2
            else
              ppmsqvx(pp(j,k,-nquark,nquark))=ppmsqvx(pp(j,k,-nquark,nquark))+ggbit1
              ppmsqvx(pp(j,k,nquark,-nquark))=ppmsqvx(pp(j,k,nquark,-nquark))+ggbit2
            endif
            if (nquark == 1) then
              msqv_cs(icol,j,k)=msqv_cs(icol,j,k)+3._dp*ggbit1
            else
              msqv_cs(icol,j,k)=msqv_cs(icol,j,k)+2._dp*ggbit1
            endif
            enddo
          enddo
      elseif ((j > 0) .and. (k < 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(j)*q1+zL(j)*zl1*prop)**2*qqbZgg_cs(icol,1,1)
     &       +abs(Q(j)*q1+zR(j)*zr1*prop)**2*qqbZgg_cs(icol,2,2)
     &       +abs(Q(j)*q1+zL(j)*zr1*prop)**2*qqbZgg_cs(icol,1,2)
     &       +abs(Q(j)*q1+zR(j)*zl1*prop)**2*qqbZgg_cs(icol,2,1)
          enddo
      elseif ((j < 0) .and. (k > 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(k)*q1+zL(k)*zl1*prop)**2*qbqZgg_cs(icol,1,1)
     &       +abs(Q(k)*q1+zR(k)*zr1*prop)**2*qbqZgg_cs(icol,2,2)
     &       +abs(Q(k)*q1+zL(k)*zr1*prop)**2*qbqZgg_cs(icol,1,2)
     &       +abs(Q(k)*q1+zR(k)*zl1*prop)**2*qbqZgg_cs(icol,2,1)
          enddo
      elseif ((j > 0) .and. (k == 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(j)*q1+zL(j)*zl1*prop)**2*qgZqg_cs(icol,1,1)
     &       +abs(Q(j)*q1+zR(j)*zr1*prop)**2*qgZqg_cs(icol,2,2)
     &       +abs(Q(j)*q1+zL(j)*zr1*prop)**2*qgZqg_cs(icol,1,2)
     &       +abs(Q(j)*q1+zR(j)*zl1*prop)**2*qgZqg_cs(icol,2,1)
          enddo
          if (j <= 4) then
          ppmsqvx(pp(jj(j),k,jj(j),k))=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
          ppmsqvx(pp(jj(j),k,k,jj(j)))=
     &       +abs(Q(j)*q1+zL(j)*zl1*prop)**2*qgZgq(1,1)
     &       +abs(Q(j)*q1+zR(j)*zr1*prop)**2*qgZgq(2,2)
     &       +abs(Q(j)*q1+zL(j)*zr1*prop)**2*qgZgq(1,2)
     &       +abs(Q(j)*q1+zR(j)*zl1*prop)**2*qgZgq(2,1)
          endif

      elseif ((j < 0) .and. (k == 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(-j)*q1+zL(-j)*zl1*prop)**2*qbgZqbg_cs(icol,1,1)
     &       +abs(Q(-j)*q1+zR(-j)*zr1*prop)**2*qbgZqbg_cs(icol,2,2)
     &       +abs(Q(-j)*q1+zL(-j)*zr1*prop)**2*qbgZqbg_cs(icol,1,2)
     &       +abs(Q(-j)*q1+zR(-j)*zl1*prop)**2*qbgZqbg_cs(icol,2,1)
          enddo
          if (j >= -4) then
          ppmsqvx(pp(jj(j),k,jj(j),k))=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
          ppmsqvx(pp(jj(j),k,k,jj(j)))=
     &       +abs(Q(-j)*q1+zL(-j)*zl1*prop)**2*qbgZgqb(1,1)
     &       +abs(Q(-j)*q1+zR(-j)*zr1*prop)**2*qbgZgqb(2,2)
     &       +abs(Q(-j)*q1+zL(-j)*zr1*prop)**2*qbgZgqb(1,2)
     &       +abs(Q(-j)*q1+zR(-j)*zl1*prop)**2*qbgZgqb(2,1)
          endif

      elseif ((j == 0) .and. (k > 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(k)*q1+zL(k)*zl1*prop)**2*gqZqg_cs(icol,1,1)
     &       +abs(Q(k)*q1+zR(k)*zr1*prop)**2*gqZqg_cs(icol,2,2)
     &       +abs(Q(k)*q1+zL(k)*zr1*prop)**2*gqZqg_cs(icol,1,2)
     &       +abs(Q(k)*q1+zR(k)*zl1*prop)**2*gqZqg_cs(icol,2,1)
          enddo
          if (k <= 4) then
          ppmsqvx(pp(j,kk(k),kk(k),j))=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
          ppmsqvx(pp(j,kk(k),j,kk(k)))=
     &       +abs(Q(k)*q1+zL(k)*zl1*prop)**2*gqZgq(1,1)
     &       +abs(Q(k)*q1+zR(k)*zr1*prop)**2*gqZgq(2,2)
     &       +abs(Q(k)*q1+zL(k)*zr1*prop)**2*gqZgq(1,2)
     &       +abs(Q(k)*q1+zR(k)*zl1*prop)**2*gqZgq(2,1)
          endif

      elseif ((j == 0) .and. (k < 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(-k)*q1+zL(-k)*zl1*prop)**2*gqbZqbg_cs(icol,1,1)
     &       +abs(Q(-k)*q1+zR(-k)*zr1*prop)**2*gqbZqbg_cs(icol,2,2)
     &       +abs(Q(-k)*q1+zL(-k)*zr1*prop)**2*gqbZqbg_cs(icol,1,2)
     &       +abs(Q(-k)*q1+zR(-k)*zl1*prop)**2*gqbZqbg_cs(icol,2,1)
          enddo
          if (k >= -4) then
          ppmsqvx(pp(j,kk(k),kk(k),j))=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)
          ppmsqvx(pp(j,kk(k),j,kk(k)))=
     &       +abs(Q(-k)*q1+zL(-k)*zl1*prop)**2*gqbZgqb(1,1)
     &       +abs(Q(-k)*q1+zR(-k)*zr1*prop)**2*gqbZgqb(2,2)
     &       +abs(Q(-k)*q1+zL(-k)*zr1*prop)**2*gqbZgqb(1,2)
     &       +abs(Q(-k)*q1+zR(-k)*zl1*prop)**2*gqbZgqb(2,1)
          endif
      endif
      msqv(j,k)=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)

   19 continue
      enddo
      enddo

      return
      end
