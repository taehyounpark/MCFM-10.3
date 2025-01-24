!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

c---  CW Virtual routine for Gamma + 2 j
c---- June 16
      subroutine qqb_gam2j_v(p,msqv)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'ewcharge.f'
      include 'zprods_com.f'
      include 'scheme.f'
      include 'epinv.f'
      include 'flags.f'
      include 'msq_cs.f'
      include 'lc.f'
      include 'first.f'
      include 'ppmax.f'
      include 'mpicommon.f'
      include 'blha.f'
      real(dp):: p(mxpart,4)
      real(dp):: msq0(-nf:nf,-nf:nf),msqv(-nf:nf,-nf:nf)
      real(dp):: mqq(0:2,-nf:nf,-nf:nf),msqx_cs(0:2,-nf:nf,-nf:nf)
      real(dp):: subuv(0:2)
      real(dp):: ppmsqx(0:2,ppmax)
      integer i,j,k,cs
      common/mqq/mqq
      integer:: rvcolourchoice
      common/rvcolourchoice/rvcolourchoice
!$omp threadprivate(/mqq/,/rvcolourchoice/)
      real(dp):: Qsum_all,Qsum_lin,fac,ga4q_noid,ga4q_id
      real(dp) :: qqb_gagg,qbq_gagg,qg_gaqg,qbg_gaqbg,
     & gq_gaqg,gqb_gaqbg,gg_gaqbq

      complex(dp):: qR_v1_a(2,2,2),qR_v2_a(2,2,2)
     &     ,qR_v1_b(2,2,2),qR_v2_b(2,2,2)
      complex(dp):: qR_a(2,2,2),qR_b(2,2,2)
      complex(dp):: qq_v1_a(2,2,2),qq_v2_a(2,2,2)
     &     ,qq_v1_b(2,2,2),qq_v2_b(2,2,2)
      complex(dp):: qq_a(2,2,2),qq_b(2,2,2)

      complex(dp):: qbR_v1_a(2,2,2),qbR_v2_a(2,2,2)
     &     ,qbR_v1_b(2,2,2),qbR_v2_b(2,2,2)
      complex(dp):: qbR_a(2,2,2),qbR_b(2,2,2)
      complex(dp):: qbq_v1_a(2,2,2),qbq_v2_a(2,2,2)
     &     ,qbq_v1_b(2,2,2),qbq_v2_b(2,2,2)
      complex(dp):: qbq_a(2,2,2),qbq_b(2,2,2)

      complex(dp):: qRb_v1_a(2,2,2),qRb_v2_a(2,2,2)
     &     ,qRb_v1_b(2,2,2),qRb_v2_b(2,2,2)
      complex(dp):: qRb_a(2,2,2),qRb_b(2,2,2)
      complex(dp):: qqb_v1_a(2,2,2),qqb_v2_a(2,2,2)
     &     ,qqb_v1_b(2,2,2),qqb_v2_b(2,2,2)
      complex(dp):: qqb_a(2,2,2),qqb_b(2,2,2)

      complex(dp):: qbRb_v1_a(2,2,2),qbRb_v2_a(2,2,2)
     &     ,qbRb_v1_b(2,2,2),qbRb_v2_b(2,2,2)
      complex(dp):: qbRb_a(2,2,2),qbRb_b(2,2,2)
      complex(dp):: qbqb_v1_a(2,2,2),qbqb_v2_a(2,2,2)
     &     ,qbqb_v1_b(2,2,2),qbqb_v2_b(2,2,2)
      complex(dp):: qbqb_a(2,2,2),qbqb_b(2,2,2)
      real(dp) :: qqb_gagg_nf,qg_gaqg_nf,gq_gaqg_nf,gg_gaqbq_nf
     &,qbq_gagg_nf,qbg_gaqbg_nf,gqb_gaqbg_nf
      logical:: qb1,qb2,ab1,ab2,gb1,gb2

      scheme = 'dred'

      if (Qflag .and. Gflag) then
        write(6,*) 'Both Qflag and Gflag cannot be true'
        write(6,*) 'They are set in file options.DAT'
        write(6,*) 'Failed in qqb_z2jet_v.f'
        stop
      endif
      if (first) then
        first=.false.
!$omp master
        if (rank == 0) then
        if ((Gflag) .or. (QandGflag)) then
          write(*,*) 'Using QQGG (VIRTUAL) matrix elements'
c          write(*,*) '[LC is     N   ]'
c          write(*,*) '[SLC is   1/N  ]'
c          write(*,*) '[SSLC is 1/N**3]'
        endif
        if ((Qflag) .or. (QandGflag)) then
          write(*,*) 'Using QQBQQB (VIRTUAL) matrix elements'
c          write(*,*) '[LC is   1 ]'
c          write(*,*) '[SLC is 1/N]'
        endif
        if     (rvcolourchoice == 1) then
          write(*,*) 'Leading colour only in VIRTUAL'
        elseif (rvcolourchoice == 2) then
          write(*,*) 'Sub-leading colour only in VIRTUAL'
        elseif (rvcolourchoice == 3) then
          write(*,*) 'Sub-sub-leading colour only in VIRTUAL'
        elseif (rvcolourchoice == 0) then
          write(*,*) 'Total of all colour structures in VIRTUAL'
        else
          write(*,*) 'Bad colourchoice'
          stop
        endif
        endif
!$omp end master
      endif

      Qsum_all=zip
      Qsum_lin=zip
      do j=1,nf
         Qsum_all=Qsum_all+Q(j)**2
         Qsum_lin=Qsum_lin+Q(j)
      enddo

      do j=-nf,nf
         do k=-nf,nf
            msqv(j,k)=0._dp
         enddo
      enddo

      call spinoru(5,p,za,zb)
      call qqb_gam2jx_new(p,msq0,mqq,ppmsqx,msqx_cs)

      qb1=.true.
      qb2=.true.
      ab1=.true.
      ab2=.true.
      gb1=.true.
      gb2=.true.

      if (Gflag) then
c----UV counterterm contains the finite renormalization to arrive
c----at MS bar scheme.
      if     (colourchoice == 1) then
        subuv(1)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
        subuv(2)=subuv(1)
      elseif (colourchoice == 2) then
        subuv(0)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
      elseif (colourchoice == 3) then
c--- all zero already
      elseif (colourchoice == 0) then
        subuv(1)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
        subuv(2)=subuv(1)
        subuv(0)=subuv(1)
      endif

c--- transfer lowest order matrix elements
c--- NB: this breaks the routine if Qflag = Gflag = .true.

      do cs=0,2
        do j=-nf,nf
        do k=-nf,nf
        msq_cs(cs,j,k)=msqx_cs(cs,j,k)
        enddo
        enddo
      enddo

      if ((qb1 .and. ab2) .or. (ab1 .and. qb2)) then
         if ((useblha == 0).or.((blhafl(1) /= 21).and.(blhafl(2) /= 21))) then
        call qqbgg_ga_v(2,1,4,5,3,za,zb,qqb_gagg,qqb_gagg_nf)
         endif
      else
        qqb_gagg=zip
        qqb_gagg_nf=zip
      endif
      if ((qb1 .and. gb2) .or. (ab1 .and. gb2)) then
         if ((useblha == 0).or.((blhafl(1) /= 21).and.(blhafl(2) == 21))) then
        call qqbgg_ga_v(4,1,2,5,3,za,zb,qg_gaqg,qg_gaqg_nf)
         endif
      else
        qg_gaqg=zip
        qg_gaqg_nf=zip
      endif
      if ((gb1 .and. qb2) .or. (gb1 .and. ab2)) then
         if ((useblha == 0).or.((blhafl(1) == 21).and.(blhafl(2) /= 21))) then
        call qqbgg_ga_v(4,2,1,5,3,za,zb,gq_gaqg,gq_gaqg_nf)
         endif
      else
        gq_gaqg=zip
        gq_gaqg_nf=zip
      endif
      if (gb1 .and. gb2) then
         if ((useblha == 0).or.((blhafl(1) == 21).and.(blhafl(2) == 21))) then
        call qqbgg_ga_v(5,4,1,2,3,za,zb,gg_gaqbq,gg_gaqbq_nf)
         endif
      else
        gg_gaqbq=zip
        gg_gaqbq_nf=zip
      endif

      qbq_gagg=qqb_gagg
      qbg_gaqbg=qg_gaqg
      gqb_gaqbg=gq_gaqg
      qbq_gagg_nf=qqb_gagg_nf
      qbg_gaqbg_nf=qg_gaqg_nf
      gqb_gaqbg_nf=gq_gaqg_nf

      fac=2._dp*esq*gsq**2*V*ason2pi

      qqb_gagg=fac*half*aveqq*qqb_gagg
      qbq_gagg=fac*half*aveqq*qbq_gagg
      gq_gaqg=aveqg*fac*gq_gaqg
      qg_gaqg=aveqg*fac*qg_gaqg

      gqb_gaqbg=aveqg*fac*gqb_gaqbg
      qbg_gaqbg=aveqg*fac*qbg_gaqbg
      gg_gaqbq=avegg*fac*gg_gaqbq

      qqb_gagg_nf=fac*half*aveqq*qqb_gagg_nf
      qbq_gagg_nf=fac*half*aveqq*qbq_gagg_nf
      gq_gaqg_nf=aveqg*fac*gq_gaqg_nf
      qg_gaqg_nf=aveqg*fac*qg_gaqg_nf

      gqb_gaqbg_nf=aveqg*fac*gqb_gaqbg_nf
      qbg_gaqbg_nf=aveqg*fac*qbg_gaqbg_nf
      gg_gaqbq_nf=avegg*fac*gg_gaqbq_nf

      endif

c***********************************************************************
c     Endpoint contributions from QQQQ matrix elements                 *
c***********************************************************************
      if (Qflag) then
        subuv(1)=2._dp*xn*(epinv*(11._dp-2._dp*real(nf,dp)/xn)-1._dp)/6._dp
        subuv(2)=subuv(1)
        subuv(0)=subuv(1)

c--- transfer lowest order matrix elements
c--- NB: this breaks the routine if Qflag = Gflag = .true.

      do cs=0,2
        do j=-nf,nf
        do k=-nf,nf
        msq_cs(cs,j,k)=mqq(cs,j,k)
        enddo
        enddo
      enddo

      endif

      if (Qflag) then
        if (useblha == 0) then
        if (qb1 .and. ab2) then
c---- q Qb -> q Qb
          call ga_qqbQQb(1,4,5,2,za,zb,qRb_a,qRb_b)
          call ga_qqbQQb_v(1,4,5,2,3,za,zb,qRb_v1_a,qRb_v2_a
     &      ,qRb_v1_b,qRb_v2_b)
c     q qb -> QQb
          call ga_qqbQQb(1,2,5,4,za,zb,qqb_a,qqb_b)
          call ga_qqbQQb_v(1,2,5,4,3,za,zb,qqb_v1_a,qqb_v2_a
     &     ,qqb_v1_b,qqb_v2_b)
        else
          qRb_a=zip; qRb_b=zip
          qRb_v1_a=zip; qRb_v2_a=zip; qRb_v1_b=zip; qRb_v2_b=zip
          qqb_a=zip; qqb_b=zip
          qqb_v1_a=zip; qqb_v2_a=zip; qqb_v1_b=zip; qqb_v2_b=zip
        endif

        if (ab1 .and. qb2) then
c---- qb Q -> qb Q
          call ga_qqbQQb(5,1,2,4,za,zb,qbR_a,qbR_b)
          call ga_qqbQQb_v(5,1,2,4,3,za,zb,qbR_v1_a,qbR_v2_a
     &      ,qbR_v1_b,qbR_v2_b)
c     qb q -> QQb
          call ga_qqbQQb(2,1,5,4,za,zb,qbq_a,qbq_b)
          call ga_qqbQQb_v(2,1,5,4,3,za,zb,qbq_v1_a,qbq_v2_a
     &     ,qbq_v1_b,qbq_v2_b)
        else
          qbR_a=zip; qbR_b=zip
          qbR_v1_a=zip; qbR_v2_a=zip; qbR_v1_b=zip; qbR_v2_b=zip
          qbq_a=zip; qbq_b=zip
          qbq_v1_a=zip; qbq_v2_a=zip; qbq_v1_b=zip; qbq_v2_b=zip
        endif
        endif

c---- q Q -> q Q
        if (qb1 .and. qb2) then
          call ga_qqbQQb(1,4,2,5,za,zb,qR_a,qR_b)
          call ga_qqbQQb_v(1,4,2,5,3,za,zb,qR_v1_a,qR_v2_a,qR_v1_b,qR_v2_b)
c     qq->qq
          if ((useblha == 0).or.(blhatype == 2)) then
          call ga_qqbQQb(1,5,2,4,za,zb,qq_a,qq_b)
          call ga_qqbQQb_v(1,5,2,4,3,za,zb,qq_v1_a,qq_v2_a,qq_v1_b,qq_v2_b)
          endif
        else
          qR_a=zip; qR_b=zip
          qR_v1_a=zip; qR_v2_a=zip; qR_v1_b=zip; qR_v2_b=zip
          qq_a=zip; qq_b=zip
          qq_v1_a=zip; qq_v2_a=zip; qq_v1_b=zip; qq_v2_b=zip
        endif

        if (useblha == 0) then
c---- qb Qb -> qb Qb
        if (ab1 .and. ab2) then
          call ga_qqbQQb(4,1,5,2,za,zb,qbRb_a,qbRb_b)
          call ga_qqbQQb_v(4,1,5,2,3,za,zb,qbRb_v1_a,qbRb_v2_a
     &      ,qbRb_v1_b,qbRb_v2_b)
c     qbqb->qbqb
          call ga_qqbQQb(5,1,4,2,za,zb,qbqb_a,qbqb_b)
          call ga_qqbQQb_v(5,1,4,2,3,za,zb,qbqb_v1_a,qbqb_v2_a,qbqb_v1_b
     &     ,qbqb_v2_b)
        else
          qbRb_a=zip; qbRb_b=zip
          qbRb_v1_a=zip; qbRb_v2_a=zip; qbRb_v1_b=zip; qbRb_v2_b=zip
          qbqb_a=zip; qbqb_b=zip
          qbqb_v1_a=zip; qbqb_v2_a=zip; qbqb_v1_b=zip; qbqb_v2_b=zip
        endif
        endif

      endif

      if(Gflag) then

c---- assemble gg pieces
      do j=-nf,nf
      do k=-nf,nf

         if( j  /=  0 .and. k  /=  0 .and. j  /=  -k) goto 19
         if     ((j == 0) .and. (k == 0)) then
            msqv(j,k)=Qsum_all*gg_gaqbq+Qsum_lin*gg_gaqbq_nf
         elseif ((j > 0) .and. (k < 0)) then
            msqv(j,k)=qqb_gagg*Q(j)**2+Q(j)*qqb_gagg_nf
         elseif ((j < 0) .and. (k > 0)) then
            msqv(j,k)=qbq_gagg*Q(-j)**2+qbq_gagg_nf*Q(-j)
         elseif ((j > 0) .and. (k == 0)) then
            msqv(j,k)=qg_gaqg*Q(j)**2+qg_gaqg_nf*Q(j)
         elseif ((j < 0) .and. (k == 0)) then
            msqv(j,k)=qbg_gaqbg*Q(-j)**2+qbg_gaqbg_nf*Q(-j)
         elseif ((j == 0) .and. (k > 0)) then
            msqv(j,k)=gq_gaqg*Q(k)**2+gq_gaqg_nf*Q(k)
         elseif ((j == 0) .and. (k < 0)) then
            msqv(j,k)=gqb_gaqbg*Q(-k)**2+gqb_gaqbg_nf*Q(-k)
         endif

 19      continue
      enddo
      enddo

      endif


      if(Qflag) then

      do j=-nf,nf
      do k=-nf,nf

         if ((j > 0) .and. (k > 0)) then
c---- QQ case
            if (j  /=  k) then
               msqv(j,k) =msqv(j,k) + ga4q_noid(j,k,qR_a,qR_b,qR_v1_a,qR_v1_b
     &              ,qR_v2_a,qR_v2_b)
            elseif(j==k) then
c-------identical quarks
               msqv(j,k) =msqv(j,k) + ga4q_id(j,k,qR_a,qR_b,qq_a,qq_b,
     &  qR_v1_a,qR_v1_b,qR_v2_a,qR_v2_b,qq_v1_a,qq_v1_b,qq_v2_a,qq_v2_b)
            endif

         elseif((j < 0) .and. (k<0)) then
            if(j /= k) then
            msqv(j,k) =msqv(j,k) + ga4q_noid(j,k,qbRb_a,qbRb_b,qbRb_v1_a
     &             ,qbRb_v1_b ,qbRb_v2_a,qbRb_v2_b)
         elseif(j==k) then
c-------identical quarks
               msqv(j,k) =msqv(j,k) + ga4q_id(j,k,qbRb_a,qbRb_b,qbqb_a,qbqb_b,
     &  qbRb_v1_a,qbRb_v1_b,qbRb_v2_a,qbRb_v2_b,qbqb_v1_a,qbqb_v1_b,qbqb_v2_a,qbqb_v2_b)
            endif
         elseif((j>0).and.(k<0)) then
            if(j /= -k) then
             msqv(j,k) =msqv(j,k) + ga4q_noid(j,k,qRb_a,qRb_b,qRb_v1_a
     &             ,qRb_v1_b ,qRb_v2_a,qRb_v2_b)
            elseif(j==-k) then
              do i = 1,nf
                if(i /= j) then
             msqv(j,k) =msqv(j,k) + ga4q_noid(j,i,qqb_a,qqb_b,qqb_v1_a
     &                  ,qqb_v1_b ,qqb_v2_a,qqb_v2_b)
                else
              msqv(j,k)= msqv(j,k) +2._dp*ga4q_id(i,i,qRb_a,qRb_b,qqb_a,qqb_b,
     &  qRb_v1_a,qRb_v1_b,qRb_v2_a,qRb_v2_b,qqb_v1_a,qqb_v1_b,qqb_v2_a,qqb_v2_b)
                endif
              enddo
            endif
       elseif((j<0).and.(k>0)) then
          if(j /= -k) then
             msqv(j,k) =msqv(j,k) + ga4q_noid(j,k,qbR_a,qbR_b,qbR_v1_a
     &            ,qbR_v1_b ,qbR_v2_a,qbR_v2_b)
          elseif(j==-k) then
            do i=1,nf
              if(i /= k) then
             msqv(j,k) =msqv(j,k) + ga4q_noid(k,i,qbq_a,qbq_b,qbq_v1_a
     &                  ,qbq_v1_b ,qbq_v2_a,qbq_v2_b)
              else
              msqv(j,k)= msqv(j,k) + 2._dp*ga4q_id(j,j,qbR_a,qbR_b,qbq_a,qbq_b,
     &  qbR_v1_a,qbR_v1_b,qbR_v2_a,qbR_v2_b,qbq_v1_a,qbq_v1_b,qbq_v2_a,qbq_v2_b)
              endif
            enddo
          endif

        endif

      enddo
      enddo

      endif

c***********************************************************************
c     UV contributions are included here                               *
c     This is the correction to put the answer in UV renormalized      *
c     dred scheme with msbar coupling                                  *
c***********************************************************************
      do j=-nf,nf
      do k=-nf,nf

      do cs=0,2
         msqv(j,k)=msqv(j,k)-ason2pi*subuv(cs)*msq_cs(cs,j,k)
      enddo

      enddo
      enddo
      return
      end

      function ga4q_noid(j,k,lo_a,lo_b,v_del1_a,v_del1_b
     &     ,v_del2_a,v_del2_b)
      implicit none
      include 'types.f'
      real(dp) ::  ga4q_noid
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'nf.f'
      include 'ewcharge.f'
      complex(dp) :: lo_a(2,2,2),lo_b(2,2,2)
      complex(dp) :: v_del1_a(2,2,2),v_del2_a(2,2,2)
      complex(dp) :: v_del1_b(2,2,2),v_del2_b(2,2,2)
      complex(dp) :: lo_tot(2,2,2),v_del1_tot(2,2,2)
      integer :: j,k,a,b
      integer :: h1,h2,h3
      real(dp) :: fac
      a=abs(j)
      b=abs(k)

      fac=V*2._dp*gsq**2*esq*aveqq*ason4pi

c-----LO amplitude
      lo_tot(:,:,:) = +lo_a(:,:,:)*Q(a)+lo_b(:,:,:)*Q(b)

c-----total del_1 piece
      v_del1_tot(:,:,:) = +v_del1_a(:,:,:)*Q(a)+v_del1_b(:,:,:)*Q(b)
c----- and del_2 piece
c      v_del2_tot(:,:,:) = v_del2_a(:,:,:)*Q(a)+v_del2_b(:,:,:)*Q(b)

c---- note that v_del2 doesnt actually interfer with the LO stucture so we just get the following
      ga4q_noid=zip

      do h1=1,2
         do h2=1,2
            do h3=1,2
               ga4q_noid = ga4q_noid
     &            +real(conjg(lo_tot(h1,h2,h3))*v_del1_tot(h1,h2,h3),dp)
     &            +real(conjg(v_del1_tot(h1,h2,h3))*lo_tot(h1,h2,h3),dp)

            enddo
         enddo
      enddo
      ga4q_noid=ga4q_noid*fac
      return
      end



      function ga4q_id(j,k,lo_a,lo_b,lo2_a,lo2_b,v_del1_a,v_del1_b
     &     ,v_del2_a,v_del2_b,v2_del1_a,v2_del1_b
     &     ,v2_del2_a,v2_del2_b)
      implicit none
      include 'types.f'
      real(dp) ::  ga4q_id
      include 'constants.f'
      include 'qcdcouple.f'
      include 'ewcouple.f'
      include 'nf.f'
      include 'ewcharge.f'
      complex(dp) :: lo_a(2,2,2),lo_b(2,2,2)
      complex(dp) :: lo2_a(2,2,2),lo2_b(2,2,2)
      complex(dp) :: v_del1_a(2,2,2),v_del2_a(2,2,2)
      complex(dp) :: v_del1_b(2,2,2),v_del2_b(2,2,2)
      complex(dp) :: v2_del1_a(2,2,2),v2_del2_a(2,2,2)
      complex(dp) :: v2_del1_b(2,2,2),v2_del2_b(2,2,2)
      complex(dp) :: lo_tot(2,2,2),v1_del1_tot(2,2,2),v1_del2_tot(2,2,2)
      complex(dp) :: lo2_tot(2,2,2),v2_del1_tot(2,2,2),v2_del2_tot(2,2,2)
      integer :: j,k,a,b
      integer :: h1,h2,h3
      real(dp) :: fac
      a=abs(j)
      b=abs(k)

c----- overall color factor is still V, but with half for identical quarks
      fac=V*2._dp*gsq**2*esq*aveqq*ason4pi*half

c-----LO amplitude
      lo_tot(:,:,:) = lo_a(:,:,:)*Q(a)+lo_b(:,:,:)*Q(b)
      lo2_tot(:,:,:) = (lo2_a(:,:,:)*Q(a)+lo2_b(:,:,:)*Q(b))

c-----v1, total del_1 piece
      v1_del1_tot(:,:,:) = v_del1_a(:,:,:)*Q(a)+v_del1_b(:,:,:)*Q(b)
c----- and del_2 piece
      v1_del2_tot(:,:,:) = (v_del2_a(:,:,:)*Q(a)+v_del2_b(:,:,:)*Q(b))


c-----v2, total del_1 piece
      v2_del1_tot(:,:,:) = (v2_del1_a(:,:,:)*Q(a)+v2_del1_b(:,:,:)*Q(b))
c----- v2, and del_2 piece
      v2_del2_tot(:,:,:) = (v2_del2_a(:,:,:)*Q(a)+v2_del2_b(:,:,:)*Q(b))

c---  answer is
c(-1 + Nc^2)  (ALOC2 (v1del2 + v2de1l) + ALOC (v1de1l + v2del2))

      ga4q_id=zip

c==== two copies of non-identical bit
      do h1=1,2
         do h2=1,2
            do h3=1,2
               ga4q_id = ga4q_id
     &           +real(conjg(lo_tot(h1,h2,h3))*v1_del1_tot(h1,h2,h3),dp)
     &           +real(conjg(v1_del1_tot(h1,h2,h3))*lo_tot(h1,h2,h3),dp)
     &          +real(conjg(lo2_tot(h1,h2,h3))*v2_del1_tot(h1,h2,h3),dp)
     &          +real(conjg(v2_del1_tot(h1,h2,h3))*lo2_tot(h1,h2,h3),dp)

            enddo
         enddo
      enddo
c===== 1/xn intf
      do h1=1,2
         h2=h1
c     do h2=1,2
            do h3=1,2
               ga4q_id = ga4q_id
     &    +(real(conjg(lo_tot(h1,h2,h3))*v2_del2_tot(h1,h2,h3),dp)
     &    +real(conjg(v2_del2_tot(h1,h2,h3))*lo_tot(h1,h2,h3),dp))*(-one)
     &    + (real(conjg(lo2_tot(h1,h2,h3))*v1_del2_tot(h1,h2,h3),dp)
     &     +real(conjg(v1_del2_tot(h1,h2,h3))*lo2_tot(h1,h2,h3),dp))*(-one)
c            enddo
         enddo
      enddo

      ga4q_id=ga4q_id*fac
 !     write(6,*) epinv2,ga4q_noid
      return
      end
