!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine gen_realps(r,p,pswt,*)
          use Multichannel
      implicit none
      include 'types.f'
c--- calls appropriate subroutine to choose a
c--- real-emission phase space point for the process
c--- specified by the variable 'case'
c---
c---    input: r of random numbers, r
c---    output: momenta p and phase-space weight pswt
c---    note: common block 'npart' is also filled here
c---
c---    alternate return if generated point should be discarded

      include 'constants.f'
      include 'nf.f'
      include 'mxpart.f'
      include 'cplx.h'
      include 'breit.f'
      include 'dm_params.f'
      include 'frag.f'
      include 'ipsgen.f'
      include 'limits.f'
      include 'masses.f'
      include 'npart.f'
      include 'kpart.f'
      include 'phasemin.f'
      include 'kprocess.f'
      include 'energy.f'
      include 'taucut.f'
      include 'vegas_common.f'
      integer:: ii,nj
      real(dp):: r(mxdim),p(mxpart,4),pvec(4),pswt,
     & ptmp,m3,m4,m5,wt34,wt345,wt346,wt3456,wtprop,zaprop,
     & s34,s345,s346,s3456,dot,wtips(4),s126,wt126

      real(dp) :: plo(mxpart,4), pone(mxpart,4)
      real(dp) :: pswtdip

      real(dp), allocatable :: dipconfig(:,:)

c--- statement function
      wtprop(s34,wmass,wwidth)=(s34-wmass**2)**2+(wmass*wwidth)**2
      zaprop(s34,wmass,wwidth)=s34*((s34-wmass**2)**2+(wmass*wwidth)**2)

      wtips = 0._dp

c--- processes that use "gen3"
      if     ( (kcase==kW_only)
     &    .or. (kcase==kZ_only)
     &    .or. (kcase==kWln_ew)
     &    .or. (kcase==kggfus0)
     &    .or. (kcase==kHigaga)
     &   .or.  (kcase==kWcsbar)
     &   .or.  (kcase==kWcs_ms) ) then
        npart=3
c        if (new_pspace) then
c          call gen3a(r,p,pswt,*999)
c        else
          call gen3(r,p,pswt,*999)
c        endif

c--- processes that use "gen3jet"
      elseif ((kcase==kgamgam) .or. (kcase==kgg2gam)) then
        !npart=3
c        call gen3(r,p,pswt,*999)
        !call gen3jetgaga(r,p,pswt,*999)
c         npart = 3
c         call gen_photons_jets_res(r,2,1,p,pswt,*999)

          npart = 2
          call gen_photons_jets_res(r,2,0,plo,pswt,*999)
          allocate(dipconfig(2,3))
          dipconfig(1,:) = [1,5,2]
          dipconfig(2,:) = [2,5,1]

          npart=npart+1
          if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),
     &                  r(ndim+2),plo,p,pswtdip,dipconfig_in=dipconfig)) then
              goto 999
          endif
          pswt=pswt*pswtdip

c--- processes that use "gen3m"
      elseif ( (kcase==ktt_tot)
     &    .or. (kcase==ktt_mix)
     &    .or. (kcase==kbb_tot)
     &    .or. (kcase==kcc_tot) ) then
        m3=mass2
        m4=mass2
        m5=0._dp
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)

      elseif ((kcase==ktwojet) .or. (kcase==ktwo_ew)
     &    .or.(kcase==kdirgam) .or. (kcase==khflgam)) then
        npart=3
        if(frag) then
c---       this phase space does a better job for the photons
           call gen_photons_jets(r,1,2,p,pswt,*999)
        else
c---       this phase space does a better job for the jets
           call gen3jet(r,p,pswt,*999)
        endif


c--- processes that use "gen4"
      elseif ( (kcase==kW_cjet)
     &   .or.  (kcase==kWbfrmc)
     &   .or.  (kcase==kW_tndk)
     &   .or.  (kcase==kepem3j)   ) then
        npart=4
        call gen4(r,p,pswt,*999)

c--- processes that use "gen4"
      elseif ((kcase==kqg_tbq) .or. (kcase==kqq_tbg)) then
        npart=4
        call gen4(r,p,pswt,*999)

c--- processes that use "gen4mdk"
      elseif (kcase==k4ftwdk) then
        npart=6
        call gen4mdk(r,p,pswt,*999)

      elseif (kcase==kHi_Zga) then
        npart=4
        call gen_HZgamj(r,p,pswt,*999)

c--- processes that use "gen4mdkrad"
      elseif (kcase==kdk_4ft) then
        npart=6
        call gen4mdkrad(r,p,pswt,*999)

c--- processes that use "gen4mdk"
      elseif ( (kcase==kZ_tdkj)
     &     .or.(kcase==kH_tdkj)) then
        npart=7
        call gen5mdk(r,p,pswt,*999)

c--- processes that use "gen4_3M"
      elseif ( (kcase==ktottth) ) then
        m3=mt
        m4=mt
        m5=hmass
        npart=4
        call gen4_3m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen5"
      elseif ( (kcase==kWbbmas)
     & .or. (kcase==kWttmas)
     & .or. (kcase==kWWqqdk)
     & .or. (kcase==kZHbbdk)
     & .or. (kcase==kWHbbdk)) then
        npart=5
        call gen5(r,p,pswt,*999)
      elseif ((kcase==kZ_tjet) .or. (kcase==kH_tjet)) then
        npart=5
        taumin=(mt/sqrts)**2
        call gen5(r,p,pswt,*999)

      elseif (kcase==kHWWdkW) then
        npart=5
        call gen5h(r,p,pswt,*999)

c--- processes that use "gen6"
      elseif (
     &      (kcase==kW_twdk)
     & .or. (kcase==kWtdkay)
     & .or. (kcase==kHWWjet)
     & .or. (kcase==kHZZjet)
     & .or. (kcase==kWH1jet)
     & .or. (kcase==kZH1jet)
     & ) then
        npart=6
        if ((usescet) .and. (abovecut)) then
c          if (tauboost) then
c            call gen6(r,p,pswt,*999)
c          else
c no need to call different routine for tauboost since uses t0=1.e-15
            call genVHjjtaucut(r,p,pswt,*999)
c          endif
        else
          call gen6(r,p,pswt,*999)
        endif

c--- processes that use "gen7"
      elseif (
     &      (kcase==kqq_HWW)
     & .or. (kcase==kqq_HZZ)
     & .or. (kcase==kWH__WW)
     & .or. (kcase==kWH__ZZ)
     & .or. (kcase==kZH__WW)
     & .or. (kcase==kZH__ZZ)
     & .or. (kcase==kHWW2jt)
     & .or. (kcase==kHZZ2jt)
     & .or. (kcase==ktt_ldk)
     & .or. (kcase==ktt_hdk)
     & .or. (kcase==ktt_udk)
     & .or. (kcase==ktthWdk)
     & ) then
        npart=7
        call gen7(r,p,pswt,*999)

c--- processes that use "gen7m"
      elseif ( (kcase==ktt_bbl)
     &    .or. (kcase==ktt_bbh)
     &    .or. (kcase==ktt_bbu)
     & ) then
        m3=mt
        m4=mt
        m5=0._dp
        npart=7
      call gen7m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen9"
      elseif ( (kcase==kqq_ttw)) then
        npart=9
      call gen9_rap(r,p,pswt,*999)

      elseif ( (kcase==kttwldk)) then
        npart=9
      call gen9dk_rap(r,p,pswt,*999)

c--- processes that use "gen_njets" with an argument of "2"
c     elseif (kcase == kggfus1) then
c         npart = 3
c         call gen_njets(r,1,pone,pswt,*999)

c         allocate(dipconfig(6,3))
c         dipconfig(1,:) = (/ 1,6,2 /)
c         dipconfig(2,:) = (/ 2,6,1 /)
c         dipconfig(3,:) = (/ 1,6,5 /)
c         dipconfig(4,:) = (/ 2,6,5 /)
c         dipconfig(5,:) = (/ 5,6,1 /)
c         dipconfig(6,:) = (/ 5,6,2 /)
c         npart=npart+1
c         if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),
c    &                  r(ndim+2),pone,p,pswtdip,dipconfig_in=dipconfig)) then
c             goto 999
c         endif
c         pswt=pswt*pswtdip


      elseif (kcase == kW_1jet .or. kcase == kZ_1jet) then
          npart = 2
          call gen2(r,plo,pswt,*999)
          allocate(dipconfig(2,3))
          dipconfig(1,:) = (/ 1,5,2 /)
          dipconfig(2,:) = (/ 2,5,1 /)
          npart=npart+1
          if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),
     &                  r(ndim+1),plo,pone,pswtdip,dipconfig_in=dipconfig)) then
              goto 999
          endif
          pswt=pswt*pswtdip

          deallocate(dipconfig)
          allocate(dipconfig(6,3))
          dipconfig(1,:) = (/ 1,6,2 /)
          dipconfig(2,:) = (/ 2,6,1 /)
          dipconfig(3,:) = (/ 1,6,5 /)
          dipconfig(4,:) = (/ 2,6,5 /)
          dipconfig(5,:) = (/ 5,6,1 /)
          dipconfig(6,:) = (/ 5,6,2 /)
          npart=npart+1
          if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),
     &                  r(ndim+2),pone,p,pswtdip,dipconfig_in=dipconfig)) then
              goto 999
          endif
          pswt=pswt*pswtdip

      elseif ( (kcase==kW_1jet)
     &    .or. (kcase==kWcjet0)
     &    .or. (kcase==kZ_1jet)
     &    .or. (kcase==kH_1jet)
     &    .or. (kcase==kggfus1)
     &    .or. (kcase==khjetma)
     &    .or. (kcase==kHgagaj)
     &    .or. (kcase==kgQ__ZQ) ) then
        npart=4
        if (origkpart == kresummed) then
            call gen_njets(r,2,p,pswt,*999)
        else
c          if (new_pspace) then
c            call gen4a(r,p,pswt,*999)
c          else
            if ((usescet) .and. (abovecut)) then
c               if (tauboost) then
c                 call gen_njets(r,2,p,pswt,*999)
c               else
c no need to call different routine for tauboost since uses t0=1.e-15
                call gen4taucut(r,p,pswt,*999)
c               endif
            else
              call gen_njets(r,2,p,pswt,*999)
            endif
c          endif
        endif

c--- processes that use "gen_njets" with an argument of "3"
      elseif ( (kcase==kWbbbar)
     &    .or. (kcase==kW_2jet)
     &    .or. (kcase==kZ_2jet)
     &    .or. (kcase==kZbbbar)
     &    .or. (kcase==kW_bjet)
     &    .or. (kcase==kZ_bjet)
     &    .or. (kcase==kqq_Hqq)
     &    .or. (kcase==kqq_Hgg)
     &    .or. (kcase==kggfus2)
     &    .or. (kcase==kgagajj)) then
        npart=5
        call gen_njets(r,3,p,pswt,*999)

      elseif ((kcase==kWga_ew) .or. (kcase==kWgajew)) then
        if (kcase == kWgajew) then
          nj=1
        else
          nj=0
        endif
        npart=4+nj
        if  (ipsgen == 1) then
            call gen_Vphotons_jets(r,2,nj,p,pswt,*999) ! (34) BW
            if (nj == 1) then
              pvec(:)=p(6,:)
              p(6,:)=p(7,:)
              p(7,:)=pvec(:)
            endif
        elseif  (ipsgen == 2) then
            call gen_Vphotons_jets_dkrad2(r,2,nj,p,pswt,*999) ! (3456) BW
            if (nj == 1) then
              pvec(:)=p(6,:)
              p(6,:)=p(5,:)
              p(5,:)=pvec(:)
            endif
        elseif  (ipsgen == 3) then
           call gen_Vphotons_jets_dkrad(r,2,nj,p,pswt,*999) ! (345) BW
        elseif  (ipsgen == 4) then
           call gen_Vphotons_jets_dkrad(r,2,nj,p,pswt,*999) ! (346) BW
           do ii=1,4
              ptmp=p(5,ii)
              p(5,ii)=p(6+nj,ii)
              p(6+nj,ii)=ptmp
           enddo
        else
           write(6,*) 'Parameter ipsgen should be 1 or 2 or 3 or 4'
           write(6,*) 'ipsgen = ',ipsgen
           stop
        endif
        s34=2._dp*dot(p,3,4)
        s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
        s346=s34+2._dp*dot(p,3,6+nj)+2._dp*dot(p,4,6+nj)
        s3456=s345+s346-s34+2._dp*dot(p,5,6+nj)
        wt34=wtprop(s34,wmass,wwidth)
        wt345=wtprop(s345,wmass,wwidth)
        wt346=wtprop(s346,wmass,wwidth)
        wt3456=wtprop(s3456,wmass,wwidth)
        wtips(1)=wt345*wt346*wt3456
        wtips(2)=wt34*wt345*wt346
        wtips(3)=wt34*wt346*wt3456
        wtips(4)=wt34*wt345*wt3456
        pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2)+wtips(3)+wtips(4))

c--- processes that use "gen_Vphotons_jets"
      elseif ( (kcase==kWgamma)
     &   .or.  (kcase==kZgamma)   ) then
        npart=4
        if (ipsgen == 1) then
          call gen_Vphotons_jets(r,1,1,p,pswt,*999)
        elseif (ipsgen == 2) then
         call gen_Vphotons_jets_dkrad(r,1,1,p,pswt,*999)
        else
          write(6,*) 'Invalid value of ipsgen in gen_lops.f!'
          stop
        endif
        if ((kcase == kZgamma)) then
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=zaprop(s34,zmass,zwidth)
          wt345=zaprop(s345,zmass,zwidth)
          wtips(1)=1._dp/wt34
          wtips(2)=1._dp/wt345
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
        endif

c--- processes that use "gen_photons_jets"
      ! gen_photons_jets_res is gen_photons_jets but only with the ptgam
      ! cuts and no further tweaks.
      elseif ((kcase == kgmgmjt)) then
          !npart=4
          !call gen_photons_jets_res(r,2,2,p,pswt,*999)

          npart = 2
          call gen_photons_jets_res(r,2,0,plo,pswt,*999)
          allocate(dipconfig(2,3))
          dipconfig(1,:) = (/ 1,5,2 /)
          dipconfig(2,:) = (/ 2,5,1 /)
          npart=npart+1
          if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),
     &                  r(ndim+1),plo,pone,pswtdip,dipconfig_in=dipconfig)) then
              goto 999
          endif
          pswt=pswt*pswtdip

          deallocate(dipconfig)
          allocate(dipconfig(6,3))
          dipconfig(1,:) = (/ 1,6,2 /)
          dipconfig(2,:) = (/ 2,6,1 /)
          dipconfig(3,:) = (/ 1,6,5 /)
          dipconfig(4,:) = (/ 2,6,5 /)
          dipconfig(5,:) = (/ 5,6,1 /)
          dipconfig(6,:) = (/ 5,6,2 /)
          npart=npart+1
          if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),
     &                  r(ndim+2),pone,p,pswtdip,dipconfig_in=dipconfig)) then
              goto 999
          endif
          pswt=pswt*pswtdip

c!      elseif (kcase==kgmgmjt) then
c!        npart=4
c!        if ((usescet) .and. (abovecut)) then
c!!          if (tauboost) then
c!!            call gen_photons_jets(r,2,2,p,pswt,*999)
c!!          else
c!! no need to call different routine for tauboost since uses t0=1.e-12
c!            call gen4taucut(r,p,pswt,*999)
c!!          endif
c!        else
c!          call gen_photons_jets(r,2,2,p,pswt,*999)
c!        endif

      elseif (kcase==ktrigam) then
        npart=4
        call gen_photons_jets(r,3,1,p,pswt,*999)

c--- special treatment for Z+gamma+gamma
      elseif ((kcase==kW_2gam) .or. (kcase==kZ_2gam)) then
        npart=5
        if  (ipsgen == 1) then
            call gen_Vphotons_jets(r,2,1,p,pswt,*999) !AA+AB
        elseif  (ipsgen == 2) then
            call gen_Vphotons_jets_dkrad2(r,2,1,p,pswt,*999) !BB+BC
        elseif  (ipsgen == 3) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999) !CC+AC+CD
        elseif  (ipsgen == 4) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999) !D.e+_dpA.e+_dpBD
           do ii=1,4
              ptmp=p(5,ii)
              p(5,ii)=p(6,ii)
              p(6,ii)=ptmp
           enddo
        else
           write(6,*) 'Parameter ipsgen should be 1 or 2 or 3 or 4'
           write(6,*) 'ipsgen = ',ipsgen
           stop
        endif
        if (kcase==kW_2gam) then
c          if (vetow_2gam(p)) goto 999 ! partition PS according to ipsgen
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          s346=s34+2._dp*dot(p,3,6)+2._dp*dot(p,4,6)
          s3456=s345+s346-s34+2._dp*dot(p,5,6)
          wt34=wtprop(s34,wmass,wwidth)
          wt345=wtprop(s345,wmass,wwidth)
          wt346=wtprop(s346,wmass,wwidth)
          wt3456=wtprop(s3456,wmass,wwidth)
          wtips(1)=wt345*wt346*wt3456
          wtips(2)=wt34*wt345*wt346
          wtips(3)=wt34*wt346*wt3456
          wtips(4)=wt34*wt345*wt3456
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2)+wtips(3)+wtips(4))
        endif
      elseif((kcase==kdm_jet).or.(kcase==kdm_gam)) then
         m3=xmass
         m4=xmass
         npart=4
         call gen4m(r,p,m3,m4,0._dp,0._dp,pswt,*999)
      elseif (kcase==kZgajet) then
        npart=5
        if (ipsgen == 1) then
          call gen_Vphotons_jets(r,1,2,p,pswt,*999)
        else
          call gen_Vphotons_jets_dkrad(r,1,2,p,pswt,*999)
        endif

        s34=2._dp*dot(p,3,4)
        s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
        wt34=zaprop(s34,zmass,zwidth)
        wt345=zaprop(s345,zmass,zwidth)
        wtips(1)=wt345
        wtips(2)=wt34
        pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2))

c--- processes that use "gen_stop" with an argument of "1"
      elseif ( (kcase==kttdkay)
     &    .or. (kcase==ktdecay) ) then
        npart=5
        call gen_stop(r,1,p,pswt,*999)
      elseif (kcase == ktopanom) then
        npart=5

        select case(ipsgen)
        case (1)
            ! radiation in production
            call gen_stop(r,2,p,pswt,*999)
        case (2)
            ! radiation in decay
            call gen_stop(r,1,p,pswt,*999)
c       case (3)
c           ! in between
c           call pstopReal(r,p,pswt,*999)
        case default
            write(6,*) 'Abort in gen_realps'
            stop
        end select

c       call gen5(r,p,pswt,*999)

        s345=2._dp*dot(p,3,4)+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
        s126=2._dp*dot(p,1,2)+2._dp*dot(p,1,6)+2._dp*dot(p,2,6)

        wt345 = 1._dp/((s345-mt**2)**2 + (mt*twidth)**2)
        wt126 = 1._dp/((s126-mt**2)**2 + (mt*twidth)**2)

        wtips = 0._dp
        wtips(1) = wt345
        wtips(2) = wt126/dot(p,5,7)
c       wtips(3) = wt126*wt345

        pswt = pswt * wtips(ipsgen)/sum(wtips)

c--- processes that use "gen_stop" with an argument of "2"
      elseif (kcase==kbq_tpq) then
        npart=5
        call gen_stop(r,2,p,pswt,*999)

c--- processes that use "gen_stop" with an argument of "2"
      elseif (kcase==kt_bbar) then
        npart=5
        call gen_stop(r,2,p,pswt,*999)

c--- DEFAULT: processes that use "gen5"
      else
        npart=5
c        if (new_pspace) then
c          call gen5a(r,p,pswt,*999)
c        else
          call gen5(r,p,pswt,*999)
c        endif
      endif


      return

c--- alternate return
  999 continue
      return 1

      end

