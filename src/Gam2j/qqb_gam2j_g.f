!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine qqb_gam2j_g(p,msq)
      implicit none
      include 'types.f'
c***********************************************************************
c     Author: J.M.Campbell                                             *
c     May, 2016.                                                       *
c***********************************************************************
c---Matrix element squared averaged over initial colors and spins
c     q(-p1)+qbar(-p2) -->  gamma + f(p3)+f(p4)+g(p5)

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'zprods_com.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'flags.f'
      include 'lc.f'
      include 'first.f'
      include 'mpicommon.f'
c      include 'selectpdfs.f'
      integer:: j,k,nquark
      real(dp):: P(mxpart,4),msq(-nf:nf,-nf:nf)
      real(dp):: msq_gg,msq_qqb,msq_qbq,msq_qg,msq_gq,msq_gqb,msq_qbg
      real(dp):: fac
      real(dp):: msqi_qq(2),msqi_qbqb(2),
     &                 msqi_qqb(2),msqi_qbq(2),
     &                 msqi_qqbs(2),msqi_qbqs(2),
     &                 msqi_qg(2),msqi_qbg(2),
     &                 msqi_gqb(2),msqi_gq(2),ggtemp
      real(dp):: msqn_qq(2,2),msqn_qbqb(2,2),
     &                 msqn_qqb(2,2),msqn_qbq(2,2),
     &                 msqn_qqbs(2,2),msqn_qbqs(2,2),
     &                 msqn_qg(2,2),msqn_qbg(2,2),
     &                 msqn_gqb(2,2),msqn_gq(2,2)
      real(dp):: ampsq_1gam3g
      logical:: qb1,qb2,ab1,ab2,gb1,gb2
      integer,parameter::jj(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)
      integer,parameter::kk(-nf:nf)=(/-1,-2,-1,-2,-1,0,1,2,1,2,1/)

      if (first) then
        first=.false.
!$omp master
        if (rank == 0) then
        if ((Gflag) .or. (QandGflag)) then
          write(*,*) 'Using QQGG+G (REAL) matrix elements'
c          write(*,*) '[LC is     N   ]'
c          write(*,*) '[SLC is   1/N  ]'
c          write(*,*) '[SSLC is 1/N**3]'
        endif
        if ((Qflag) .or. (QandGflag)) then
          write(*,*) 'Using QQBQQB+G (REAL) matrix elements'
c          write(*,*) '[LC is   1 ]'
c          write(*,*) '[SLC is 1/N]'
        endif
        if     (colourchoice == 1) then
          write(*,*) 'Leading colour only in REAL'
        elseif (colourchoice == 2) then
          write(*,*) 'Sub-leading colour only in REAL'
        elseif (colourchoice == 3) then
          write(*,*) 'Sub-sub-leading colour only in REAL'
        elseif (colourchoice == 0) then
          write(*,*) 'Total of all colour structures in REAL'
        else
          write(*,*) 'Bad colourchoice'
          stop
        endif
        endif
!$omp end master
      endif

      msq(:,:)=zip

c---call spinor routine and load common block twopij
      call spinoru(6,p,za,zb)

      qb1=.true.
      qb2=.true.
      ab1=.true.
      ab2=.true.
      gb1=.true.
      gb2=.true.
c      if (selectpdfs(1)==1 .or. selectpdfs(1)==2 .or. selectpdfs(1)==99) qb1=.true.
c      if (selectpdfs(2)==1 .or. selectpdfs(2)==2 .or. selectpdfs(2)==99) qb2=.true.
c      if (selectpdfs(1)==-1 .or. selectpdfs(1)==2 .or.  selectpdfs(1)==99) ab1=.true.
c      if (selectpdfs(2)==-1 .or. selectpdfs(2)==2 .or.  selectpdfs(2)==99) ab2=.true.
c      if (selectpdfs(1)==0 .or. selectpdfs(1)==99) gb1=.true.
c      if (selectpdfs(2)==0 .or. selectpdfs(2)==99) gb2=.true.

      if (Gflag) then
c***********************************************************************
c     Calculate contributions from the QQGGG matrix elements            *
c***********************************************************************

      if (qb1 .and. ab2) then
        msq_qqb=ampsq_1gam3g(4,5,6,3,1,2,za,zb) ! q  qb -> a g  g  g
      else
        msq_qqb=zip
      endif
      if (ab1 .and. qb2) then
        msq_qbq=ampsq_1gam3g(4,5,6,3,2,1,za,zb) ! qb q  -> a g  g  g
      else
        msq_qbq=zip
      endif
      if (qb1 .and. gb2) then
        msq_qg= ampsq_1gam3g(2,5,6,3,1,4,za,zb) ! q  g  -> a q  g  g
      else
        msq_qg=zip
      endif
      if (ab1 .and. gb2) then
        msq_qbg=ampsq_1gam3g(2,5,6,3,4,1,za,zb) ! qb g  -> a qb g  g
      else
        msq_qbg=zip
      endif

c      msq_gg= ampsq_1gam3g(1,2,6,3,5,4,za,zb) ! g  g  -> a q  qb g
c      msq_gq= ampsq_1gam3g(1,4,6,3,2,5,za,zb) ! g  q  -> a g  q  g
c      msq_gqb=ampsq_1gam3g(1,4,6,3,5,2,za,zb) ! g  qb -> a g  qb g

c--- The calls above have been cross-checked with a MG calculation of the matrix elements
c--- These permutations are to match qqb_z2jet_g.f, for simplicity in the subtraction routine
      if (gb1 .and. gb2) then
        msq_gg= ampsq_1gam3g(1,2,6,3,4,5,za,zb) ! g  g  -> a qb q  g
      else
        msq_gg=zip
      endif
      if (gb1 .and. qb2) then
        msq_gq= ampsq_1gam3g(1,5,6,3,2,4,za,zb) ! g  q  -> a q  g  g
      else
        msq_gq=zip
      endif
      if (gb1 .and. ab2) then
        msq_gqb=ampsq_1gam3g(1,5,6,3,4,2,za,zb) ! g  qb -> a qb g  g
      else
        msq_gqb=zip
      endif

      fac=4._dp*esq*gsq**3*xn**3*CF

      do j=-nf,nf
      do k=-nf,nf

      if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) cycle

c-- Note: appropriate statistical factors (for gluons in final state)
c--       are included below (/six and /two)
      if     ((j == 0) .and. (k == 0)) then
          ggtemp=0._dp
          do nquark=1,nf
            ggtemp=ggtemp+Q(nquark)**2*msq_gg
          enddo
          msq(j,k)=avegg*fac*ggtemp
      elseif ((j > 0) .and. (k < 0)) then
          msq(j,k)=aveqq*fac*Q(j)**2*msq_qqb/six
      elseif ((j < 0) .and. (k > 0)) then
          msq(j,k)=aveqq*fac*Q(k)**2*msq_qbq/six
      elseif ((j > 0) .and. (k == 0)) then
          msq(j,k)=aveqg*fac*Q(j)**2*msq_qg/two
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=aveqg*fac*Q(-j)**2*msq_qbg/two
      elseif ((j == 0) .and. (k > 0)) then
          msq(j,k)=aveqg*fac*Q(k)**2*msq_gq/two
      elseif ((j == 0) .and. (k < 0)) then
          msq(j,k)=aveqg*fac*Q(-k)**2*msq_gqb/two
      endif
      enddo
      enddo

      endif

      if (Qflag) then
      fac=4._dp*esq*gsq**3*xn**2*Cf
c--- extra factor of 2**3=8 to compensate for Ta normalization

c--- note: the following two arrays end up being overall 1<->2 symmetric ! CHECK THIS
      if (qb1 .and. ab2) then
        call ampsq_1gam1g2q(1,4,5,2,6,3,za,zb,msqn_qqbs,msqi_qqbs) ! q  rb -> a q  rb g
        call ampsq_1gam1g2q(1,2,5,4,6,3,za,zb,msqn_qqb,msqi_qqb)   ! q  qb -> a r  rb g
      else
        msqn_qqbs=zip; msqi_qqbs=zip
        msqn_qqb=zip; msqi_qqb=zip
      endif
c      call ampsq_1gam1g2q(4,1,2,5,6,3,za,zb,msqn_qbqs,msqi_qbqs) ! qb r  -> a qb r  g
c      call ampsq_1gam1g2q(2,1,4,5,6,3,za,zb,msqn_qbq,msqi_qbq)   ! qb q  -> a rb r  g

      if (ab1 .and. ab2) then
        call ampsq_1gam1g2q(4,1,5,2,6,3,za,zb,msqn_qbqb,msqi_qbqb) ! qb rb -> a qb rb g
      else
        msqn_qbqb=zip; msqi_qbqb=zip
      endif
      if (qb1 .and. qb2) then
        call ampsq_1gam1g2q(1,4,2,5,6,3,za,zb,msqn_qq,msqi_qq)     ! q  r  -> a q  r  g
      else
        msqn_qq=zip; msqi_qq=zip
      endif

c      call ampsq_1gam1g2q(1,6,4,5,2,3,za,zb,msqn_qg,msqi_qg)     ! q  g  -> a rb r  q
c      call ampsq_1gam1g2q(6,2,4,5,1,3,za,zb,msqn_gqb,msqi_gqb)   ! g  qb -> a rb r  qb
c      call ampsq_1gam1g2q(2,6,4,5,1,3,za,zb,msqn_gq,msqi_gq)     ! g  q  -> a rb r  q
c      call ampsq_1gam1g2q(6,1,4,5,2,3,za,zb,msqn_qbg,msqi_qbg)   ! qb g  -> a rb r  qb

c--- The calls above have been cross-checked with a MG calculation of the matrix elements
c--- These permutations are to match qqb_z2jet_g.f, for simplicity in the subtraction routine
      if (ab1 .and. qb2) then
        call ampsq_1gam1g2q(5,1,2,4,6,3,za,zb,msqn_qbqs,msqi_qbqs) ! qb r  -> a r  qb g
        call ampsq_1gam1g2q(2,1,5,4,6,3,za,zb,msqn_qbq,msqi_qbq)   ! qb q  -> a r  rb g
      else
        msqn_qbqs=zip; msqi_qbqs=zip
        msqn_qbq=zip; msqi_qbq=zip
      endif

      if (gb1 .and. ab2) then
        call ampsq_1gam1g2q(6,2,5,4,1,3,za,zb,msqn_gqb,msqi_gqb)   ! g  qb -> a r  rb qb
      else
        msqn_gqb=zip; msqi_gqb=zip
      endif
      if (gb1 .and. qb2) then
        call ampsq_1gam1g2q(2,6,5,4,1,3,za,zb,msqn_gq,msqi_gq)     ! g  q  -> a r  rb q
      else
        msqn_gq=zip; msqi_gq=zip
      endif

      if (qb1 .and. gb2) then
        call ampsq_1gam1g2q(1,6,5,4,2,3,za,zb,msqn_qg,msqi_qg)     ! q  g  -> a r  rb q
      else
        msqn_qg=zip; msqi_qg=zip
      endif
      if (ab1 .and. gb2) then
        call ampsq_1gam1g2q(6,1,5,4,2,3,za,zb,msqn_qbg,msqi_qbg)   ! qb g  -> a r  rb qb
      else
        msqn_qbg=zip; msqi_qbg=zip
      endif

      do j=-nf,nf
      do k=-nf,nf

      if ((j > 0) .and. (k < 0)) then
c-qqb
           if (k==-j) then
            msq(j,k)=msq(j,k)+fac*aveqq*(msqi_qqbs(jj(j))
     &      +real(1+jj(j),dp)*msqn_qqb(jj(j),1)
     &      +real(3-jj(j),dp)*msqn_qqb(jj(j),2))
           else
            msq(j,k)=msq(j,k)+fac*aveqq*msqn_qqbs(jj(j),-kk(k))
           endif
      elseif ((j < 0) .and. (k > 0)) then
c-qbq
          if (j ==-k) then
            msq(j,k)=msq(j,k)+fac*aveqq*(msqi_qbqs(kk(k))
     &      +real(1+kk(k),dp)*msqn_qbq(kk(k),1)
     &      +real(3-kk(k),dp)*msqn_qbq(kk(k),2))
          else
           msq(j,k)=msq(j,k)+fac*aveqq*msqn_qbqs(-jj(j),kk(k))
          endif
      elseif ((j > 0) .and. (k == 0)) then
c-qg
          msq(j,k)=msq(j,k)
     &     +fac*aveqg*(half*msqi_qg(jj(j))
     &      +real(1+jj(j),dp)*msqn_qg(jj(j),1)
     &      +real(3-jj(j),dp)*msqn_qg(jj(j),2)
     & )
c-qbg
      elseif ((j < 0) .and. (k == 0)) then
          msq(j,k)=msq(j,k)
     &     +fac*aveqg*(half*msqi_qbg(-jj(j))
     &      +real(1-jj(j),dp)*msqn_qbg(-jj(j),1)
     &      +real(3+jj(j),dp)*msqn_qbg(-jj(j),2)
     & )
      elseif ((j == 0) .and. (k > 0)) then
c-gq
          msq(j,k)=msq(j,k)
     &      +fac*aveqg*(half*msqi_gq(kk(k))
     &      +real(1+kk(k),dp)*msqn_gq(kk(k),1)
     &      +real(3-kk(k),dp)*msqn_gq(kk(k),2)
     & )
      elseif ((j == 0) .and. (k < 0)) then
c-gqb
          msq(j,k)=msq(j,k)
     &      +fac*aveqg*(half*msqi_gqb(-kk(k))
     &      +real(1-kk(k),dp)*msqn_gqb(-kk(k),1)
     &      +real(3+kk(k),dp)*msqn_gqb(-kk(k),2)
     & )
      elseif ((j > 0) .and. (k > 0)) then
c-qq
          if (j==k) then
          msq(j,k)=msq(j,k)+half*fac*aveqq*msqi_qq(jj(j))
          else
          msq(j,k)=msq(j,k)+fac*aveqq*msqn_qq(jj(j),kk(k))
         endif
      elseif ((j < 0) .and. (k < 0)) then
c-qbqb
          if (j==k) then
          msq(j,k)=msq(j,k)+half*fac*aveqq*msqi_qbqb(-jj(j))
          else
          msq(j,k)=msq(j,k)+fac*aveqq*msqn_qbqb(-jj(j),-kk(k))
          endif
      endif

      enddo
      enddo

      endif

      return
      end
