!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      subroutine qqb_z2jet_gvec(p,n,in,msqv)
      implicit none
      include 'types.f'

c-----Nov 10 99 --- checked that gives the right matrix element
c     when summed over polarizations.

c----Matrix element for Z+2jet production
c----averaged over initial colours and spins
c    line 6 contracted with the vector n(mu)
c     q(-p1)+qbar(-p2)--> g(p5)+ g(p6)+Z(f(p3)+af(p4))

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'masses.f'
      include 'qcdcouple.f'
      include 'zcouple.f'
      include 'ewcharge.f'
      include 'zprods_decl.f'
      include 'sprods_com.f'
      include 'msqv_cs.f'
      include 'flags.f'
      include 'zcouple_cms.f'

c ip is the label of the emitter
c in is the label of the contracted line
      integer:: j,k,pq,pl,nquark,in,icol
      real(dp):: fac,n(4)
      complex(dp):: prop,zab(mxpart,mxpart),zba(mxpart,mxpart)
      real(dp):: msqv(-nf:nf,-nf:nf),p(mxpart,4),ggtemp(0:2),
     & qqbZgg_cs(0:2,2,2),qbqZgg_cs(0:2,2,2),qgZqg_cs(0:2,2,2),
     & gqZqg_cs(0:2,2,2),qbgZqbg_cs(0:2,2,2),gqbZqbg_cs(0:2,2,2),
     & ggZqbq_cs(0:2,2,2),
     & qqbZgg(2,2),qbqZgg(2,2),qgZqg(2,2),gqZqg(2,2),
     & qbgZqbg(2,2),gqbZqbg(2,2),ggZqbq(2,2)
      include 'cplx.h'


      do j=-nf,nf
      do k=-nf,nf
      msqv(j,k)=0._dp
      do icol=0,2
        msqv_cs(icol,j,k)=0._dp
      enddo
      enddo
      enddo

      call spinoru(6,p,za,zb)
c   zab=<i-|k|j-> zba=<i+|k|j+> where k is an arbitrary 4-vector
c---Conventions of Bern, Dixon, Kosower, Weinzierl,
c---ie, za(i,j)*zb(j,i)=s(i,j)
      call spinork(6,p,zab,zba,n)

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
      enddo
      enddo
      enddo

c--- note the interchange of color structures in contributions for calls with q <-> qb
c--- (these call storezcsv_swapcol instead of storezcsv)
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
      elseif (in == 2)  then
           call z2jetsqn(5,6,1,2,p,n,za,zb,zab,zba,ggZqbq)
           call storezcsv_swapcol(ggZqbq_cs)
           call z2jetsqn(5,1,6,2,p,n,za,zb,zab,zba,qbgZqbg)
           call storezcsv(qbgZqbg_cs)
           call z2jetsqn(1,5,6,2,p,n,za,zb,zab,zba,qgZqg)
           call storezcsv(qgZqg_cs)
      elseif (in == 5) then
          call z2jetsqn(1,2,6,5,p,n,za,zb,zab,zba,qqbZgg)
          call storezcsv(qqbZgg_cs)
          call z2jetsqn(2,1,6,5,p,n,za,zb,zab,zba,qbqZgg)
          call storezcsv(qbqZgg_cs)
c The rest of the in=5 cases should not be used since they correspond to
c a gluon in position 5 (abnormal)
          call z2jetsqn(1,6,2,5,p,n,za,zb,zab,zba,qgZqg)
          call storezcsv(qgZqg_cs)
          call z2jetsqn(2,6,1,5,p,n,za,zb,zab,zba,gqZqg)
          call storezcsv(gqZqg_cs)
          call z2jetsqn(6,1,2,5,p,n,za,zb,zab,zba,qbgZqbg)
          call storezcsv(qbgZqbg_cs)
          call z2jetsqn(6,2,1,5,p,n,za,zb,zab,zba,gqbZqbg)
          call storezcsv(gqbZqbg_cs)
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
      fac=v*xn/four*(abs(zesq)*gsq)**2*two

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
      enddo
      enddo
      enddo

c DEBUG: only keep one color ordering
c      do icol=0,2,2
c      qqbZgg_cs(icol,:,:) =zip
c      qbqZgg_cs(icol,:,:) =zip
c      gqZqg_cs(icol,:,:)  =zip
c      qgZqg_cs(icol,:,:)  =zip
c      gqbZqbg_cs(icol,:,:)=zip
c      qbgZqbg_cs(icol,:,:)=zip
c      ggZqbq_cs(icol,:,:) =zip
c      enddo

      if (Gflag) then
      do j=-nf,nf
      do k=-nf,nf
      if( (j  /=  0) .and. (k  /=  0) .and. (j  /=  -k)) goto 19

      if     ((j == 0) .and. (k == 0)) then

           do icol=0,2
           ggtemp(icol)=0._dp
             do nquark=1,nf
             ggtemp(icol)=ggtemp(icol)
     &       +abs(Q(nquark)*q1+zL(nquark)*zl1*prop)**2*ggZqbq_cs(icol,1,1)
     &       +abs(Q(nquark)*q1+zR(nquark)*zr1*prop)**2*ggZqbq_cs(icol,2,2)
     &       +abs(Q(nquark)*q1+zL(nquark)*zr1*prop)**2*ggZqbq_cs(icol,1,2)
     &       +abs(Q(nquark)*q1+zR(nquark)*zl1*prop)**2*ggZqbq_cs(icol,2,1)
             enddo
           msqv_cs(icol,j,k)=ggtemp(icol)
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
      elseif ((j < 0) .and. (k == 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(-j)*q1+zL(-j)*zl1*prop)**2*qbgZqbg_cs(icol,1,1)
     &       +abs(Q(-j)*q1+zR(-j)*zr1*prop)**2*qbgZqbg_cs(icol,2,2)
     &       +abs(Q(-j)*q1+zL(-j)*zr1*prop)**2*qbgZqbg_cs(icol,1,2)
     &       +abs(Q(-j)*q1+zR(-j)*zl1*prop)**2*qbgZqbg_cs(icol,2,1)
          enddo
      elseif ((j == 0) .and. (k > 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(k)*q1+zL(k)*zl1*prop)**2*gqZqg_cs(icol,1,1)
     &       +abs(Q(k)*q1+zR(k)*zr1*prop)**2*gqZqg_cs(icol,2,2)
     &       +abs(Q(k)*q1+zL(k)*zr1*prop)**2*gqZqg_cs(icol,1,2)
     &       +abs(Q(k)*q1+zR(k)*zl1*prop)**2*gqZqg_cs(icol,2,1)
          enddo
      elseif ((j == 0) .and. (k < 0)) then
          do icol=0,2
             msqv_cs(icol,j,k)=
     &       +abs(Q(-k)*q1+zL(-k)*zl1*prop)**2*gqbZqbg_cs(icol,1,1)
     &       +abs(Q(-k)*q1+zR(-k)*zr1*prop)**2*gqbZqbg_cs(icol,2,2)
     &       +abs(Q(-k)*q1+zL(-k)*zr1*prop)**2*gqbZqbg_cs(icol,1,2)
     &       +abs(Q(-k)*q1+zR(-k)*zl1*prop)**2*gqbZqbg_cs(icol,2,1)
          enddo
      endif
      msqv(j,k)=msqv_cs(0,j,k)+msqv_cs(1,j,k)+msqv_cs(2,j,k)

   19 continue
      enddo
      enddo
      endif

      return
      end


      subroutine storezcsv(p1p2)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the Z2jet_gvec matrix elements into elements of p1p2

      include 'mmsqv_cs.f'
      integer:: icol,pq,pl
      real(dp):: p1p2(0:2,2,2)

      do pq=1,2
      do pl=1,2
      do icol=0,2
        p1p2(icol,pq,pl)=mmsqv_cs(icol,pq,pl)
      enddo
      enddo
      enddo

      return
      end



      subroutine storezcsv_swapcol(p1p2)
      implicit none
      include 'types.f'
c-- this routine transfers the information on the colour structure
c-- for the Z2jet_gvec matrix elements into elements of p1p2
c-- but with color orderings 1 and 2 swapped

      include 'mmsqv_cs.f'
      integer:: icol,pq,pl
      real(dp):: p1p2(0:2,2,2)

      do pq=1,2
      do pl=1,2
      p1p2(0,pq,pl)=mmsqv_cs(0,pq,pl)
      do icol=1,2
        p1p2(icol,pq,pl)=mmsqv_cs(3-icol,pq,pl)
      enddo
      enddo
      enddo

      return
      end

