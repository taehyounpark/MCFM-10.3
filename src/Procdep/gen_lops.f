!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine gen_lops(r,p,pswt,*)
          use Multichannel
          use VVconfig_m
          use SCET, only: useQT
      implicit none
      include 'types.f'
c--- calls appropriate subroutine to choose a
c--- leading-order phase space point for the process
c--- specified by the variable 'case'
c---
c---    input: vector of random numbers, r
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
      include 'ipsgen.f'
      include 'limits.f'
      include 'masses.f'
c      include 'mxdim.f'
      include 'npart.f'
      include 'phasemin.f'
      include 'kprocess.f'
      include 'zerowidth.f'
      include 'energy.f'
      include 'new_pspace.f'
      include 'vegas_common.f'
      include 'taucut.f'
      include 'kpart.f'
      include 'hdecaymode.f'
      integer:: ii
      real(dp):: r(mxdim),p(mxpart,4),pswt,
     & ptmp,m3,m4,m5,wt34,wt345,wt346,wt3456,wtprop,zaprop,
     & s34,s345,s346,s3456,dot,wtips(4),plo(mxpart,4),pswtdip,
     & rmass,rwidth

      real(dp), allocatable :: dipconfig(:,:)

c--- statement function
      wtprop(s34,rmass,rwidth)=(s34-rmass**2)**2+(rmass*rwidth)**2
      zaprop(s34,rmass,rwidth)=s34*((s34-rmass**2)**2+(rmass*rwidth)**2)

      pswt=zip

c--- processes that use "gen2"
      if     ( (kcase==kW_only)
     &    .or. (kcase==kZ_only)
     &    .or. (kcase==kWln_ew)
     &    .or. (kcase==kWln_aq)
     &    .or. (kcase==kHigaga)
     &    .or. (kcase==kWcsbar)
     &    .or. (kcase==kWcs_ms)
     &    .or. (kcase==kgg2lep)
     &    .or. (kcase==kvlchk2) ) then
        if (kcase==kvlchk2) then
          wsqmin=0._dp
          wsqmax=sqrts**2
        endif
        npart=2
c        if (new_pspace) then
c          call gen2a(r,p,pswt,*999)
c        else
          call gen2(r,p,pswt,*999)
c        endif

c--- processes that use "gen2jet"
      elseif ((kcase==ktwojet)
     &   .or. (kcase==kdirgam)
     &   .or. (kcase==khflgam)
     &   .or. (kcase==kgamgam)
     &   .or. (kcase==kgg2gam)
     &   .or. (kcase==ktwo_ew)) then
        npart=2
        call gen_photons_jets_res(r,2,0,p,pswt,*999)
        !call gen2jet(r,p,pswt,*999)

c--- processes that use "gen2m"
      elseif ( (kcase==kggfus0)
     &    .or. (kcase==kggscl0)
     &    .or. (kcase==ktt_tot)
     &    .or. (kcase==ktt_mix)
     &    .or. (kcase==kbb_tot)
     &    .or. (kcase==kcc_tot) ) then
        npart=2
        call gen2m(r,p,pswt,*999)

c--- processes that use "gen3"
      elseif ( (kcase==kW_cjet)
     &   .or.  (kcase==kWbfrmc)
     &   .or.  (kcase==kW_tndk)
     &   .or.  (kcase==kvlchwn)
     &   .or.  (kcase==kepem3j)
     &   .or.  (kcase==kvlchk3) ) then
        if (kcase==kvlchk3) then
          wsqmin=0._dp
          wsqmax=sqrts**2
        endif
        npart=3
        call gen3(r,p,pswt,*999)

      elseif((kcase==kdm_jet).or.(kcase==kdm_gam)) then
         m3=xmass
         m4=m3
         m5=0._dp
         npart=3
         call gen3m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen3h"
      elseif (kcase==kHi_Zga) then
        npart=3
        call gen3h(r,p,pswt,*999)

      elseif (kcase==kHi_Zaj) then
        npart=4
        call gen_HZgamj(r,p,pswt,*999)

c---  processes that use "gen3jet"
      elseif ((kcase == kZgamma) .and. (origkpart == kresummed)) then
        npart=3
        call gen_Vphotons_jets(r,1,0,p,pswt,*999)
        write(6,*) 'Abort in gen_lops'
        stop
        !call gen3(r,p,pswt,*999)
      elseif ( (kcase==kWgamma)
     &   .or.  (kcase==kZgamma)
     &   .or.  (kcase==kWgaj_a)
     &   .or.  (kcase==kW_frag)
     &   .or.  (kcase==kZ_frag)) then
        npart=3
        if (ipsgen == 1) then
          call gen_Vphotons_jets(r,1,0,p,pswt,*999)
        elseif (ipsgen > 1) then
          call gen_Vphotons_jets_dkrad(r,1,0,p,pswt,*999)
c          call genVphoton(r,p,pswt,*999)
        else
          write(6,*) 'Invalid value of ipsgen in gen_lops.f!'
          stop
        endif
        if ((kcase == kZgamma) .and. ((decayChannel() == decayElAntiEl))) then
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=zaprop(s34,zmass,zwidth)
          wt345=zaprop(s345,zmass,zwidth)
          wtips(1)=wt345
          wtips(2)=wt34
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
        endif
        if (((kcase == kWgamma).or.(kcase == kWgaj_a).or.(kcase == kWga_ew))) then
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=wtprop(s34,wmass,wwidth)
          wt345=wtprop(s345,wmass,wwidth)
          wtips(1)=wt345
          wtips(2)=wt34
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
        endif

c--- processes that use "gen3jetgaga"
      elseif ((kcase == kgmgmjt .or. kcase==kgg2gamjt)) then
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

c!!      elseif (kcase==kgmgmjt) then
c!!        npart=3
c!!c        call gen3jetgaga(r,p,pswt,*999)
c!!c        call gen_photons_jets(r,2,1,p,pswt,*999)
c!!        if     ((usescet) .and. (abovecut) .and. (kpart /= kreal)) then
c!!          if (tauboost) then
c!!            call gen3(r,p,pswt,*999) ! gen3taucut tailored to non-boosted case
c!!          else
c!!            call gen3taucut(r,p,pswt,*999)
c!!          endif
c!!        else
c          call gen3(r,p,pswt,*999)
c!!        endif

      elseif (kcase==kgmgmjj) then
        npart=4
c        call gen3jetgaga(r,p,pswt,*999)
        call gen_photons_jets(r,2,2,p,pswt,*999)

c--- processes that use "gen3jetgaga"
      elseif (kcase==ktrigam) then
        npart=3
c        call gen3jetgaga(r,p,pswt,*999)
c        call gen_photons_jets(r,3,0,p,pswt,*999)
        call gen3(r,p,pswt,*999)

      elseif(kcase==kfourga) then
        npart=4
        call gen_photons_jets(r,4,0,p,pswt,*999)

c--- special treatment for W+gamma+jet and Z+gamma+jet
      elseif ( (kcase==kWgajet) .or. (kcase==kZgajet) ) then
        npart=4
        if     (ipsgen == 1) then
          if     (((usescet) .and. (abovecut)) .or. (origkpart==kresummed)) then
            call genVgataucut(r,p,pswt,*999)
          else
            call gen4(r,p,pswt,*999)
c            call gen_Vphotons_jets(r,1,1,p,pswt,*999) ! MCFM-8.0 call
          endif
        elseif (ipsgen == 2) then
          if     (((usescet) .and. (abovecut)) .or. (origkpart==kresummed)) then
            call genVgataucut_dkrad(r,p,pswt,*999)
          else
            call gen4(r,p,pswt,*999)
c            call gen_Vphotons_jets_dkrad(r,1,1,p,pswt,*999)
          endif
        else
          write(6,*) 'Parameter ipsgen should be 1 or 2'
          write(6,*) 'ipsgen = ',ipsgen
          stop
        endif
        if ((kcase == kZgajet) .and. (decayChannel() == decayElAntiEl)) then
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=zaprop(s34,zmass,zwidth)
          wt345=zaprop(s345,zmass,zwidth)
          wtips(1)=wt345
          wtips(2)=wt34
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
        endif
        if ((kcase == kWgajet) .or. (kcase==kWgajew)) then
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=wtprop(s34,wmass,wwidth)
          wt345=wtprop(s345,wmass,wwidth)
          wtips(1)=wt345
          wtips(2)=wt34
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
        endif
c      endif

c--- special treatment for W+gamma (EW corrections)
      elseif (kcase==kWga_ew) then
        npart=3
        if  ((ipsgen == 1) .or. (ipsgen >= 3)) then
            call gen_Vphotons_jets(r,1,0,p,pswt,*999) ! (34) BW
        elseif  (ipsgen == 2) then
           call gen_Vphotons_jets_dkrad(r,1,0,p,pswt,*999) ! (345) BW
        endif

        if (maxipsgen == 2) then
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=wtprop(s34,wmass,wwidth)
          wt345=wtprop(s345,wmass,wwidth)
          wtips(1)=wt345
          wtips(2)=wt34
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
        endif

c--- special treatment for W+gamma+jet (EW corrections)
      elseif ((kcase==kWgajew) .or. (kcase==kWgajja))then
        npart=4
        if  ((ipsgen == 1) .or. (ipsgen >= 3)) then
            call gen_Vphotons_jets(r,1,1,p,pswt,*999) ! (34) BW
        elseif  (ipsgen == 2) then
           call gen_Vphotons_jets_dkrad(r,1,1,p,pswt,*999) ! (345) BW
        endif

        if (maxipsgen == 2) then
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=wtprop(s34,wmass,wwidth)
          wt345=wtprop(s345,wmass,wwidth)
          wtips(1)=wt345
          wtips(2)=wt34
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
        endif

c--- special treatment for Z+gamma+gamma
      elseif ( (kcase==kZ_2gam) .or. (kcase==kW_2gam)) then
        npart=4
        if  (ipsgen == 1) then
            call gen_Vphotons_jets(r,2,0,p,pswt,*999) !AA+AB
        elseif  (ipsgen == 2) then
            call gen_Vphotons_jets_dkrad2(r,2,0,p,pswt,*999) !BB+BC
        elseif  (ipsgen == 3) then
           call gen_Vphotons_jets_dkrad(r,2,0,p,pswt,*999) !CC+AC+CD
        elseif  (ipsgen == 4) then
           call gen_Vphotons_jets_dkrad(r,2,0,p,pswt,*999) !D.e+_dpA.e+_dpBD
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

c--- special treatment for Z+gamma+gamma+jet
      elseif ( (kcase==kZ2gajt) ) then
        npart=5
        if (ipsgen == 1) then
           call gen_Vphotons_jets(r,2,1,p,pswt,*999)  !AA+AB
        elseif (ipsgen == 2) then
           call gen_Vphotons_jets_dkrad2(r,2,1,p,pswt,*999)  !BB+BC
        elseif (ipsgen == 3) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999)  !CC+AC+CD
        elseif (ipsgen == 4) then
           call gen_Vphotons_jets_dkrad(r,2,1,p,pswt,*999)  !D.e+_dpA.e+_dpBD
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

c--- special treatment for Z+gamma+jet+jet
      elseif ( (kcase==kZga2jt) .or. (kcase==kWga2jt) ) then
        npart=5
c        if (ipsgen == 1) then
c          call gen_Vphotons_jets(r,1,2,p,pswt,*999)
c        elseif (ipsgen == 2) then
c          call gen_Vphotons_jets_dkrad(r,1,2,p,pswt,*999)
c        else
c          write(6,*) 'Parameter ipsgen should be 1 or 2'
c          write(6,*) 'ipsgen = ',ipsgen
c          stop
c        endif
        call gen5(r,p,pswt,*999)
        if (maxipsgen == 2) then
          s34=2._dp*dot(p,3,4)
          s345=s34+2._dp*dot(p,3,5)+2._dp*dot(p,4,5)
          wt34=zaprop(s34,zmass,zwidth)
          wt345=zaprop(s345,zmass,zwidth)
          wtips(1)=wt345
          wtips(2)=wt34
          pswt=pswt*wtips(ipsgen)/(wtips(1)+wtips(2))
        endif

      elseif ((kcase==kthrjet) .or. (kcase==kgamjet)) then
c        call gen3jet(r,p,pswt,*999)
        call gen_photons_jets(r,1,2,p,pswt,*999)
        npart=3

c--- processes that use "gen3m"
      elseif ( (kcase==ktottth) ) then
        m3=mt
        m4=mt
        m5=hmass
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)

      elseif ( (kcase==ktotttz) ) then
        m3=mt
        m4=mt
        m5=zmass
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen3m"
      elseif (kcase==ktt_glu) then
        m3=mass2
        m4=mass2
        m5=0._dp
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen3m"
      elseif ((kcase==kqg_tbq) .or. (kcase==kqq_tbg)) then
        m3=mt
        m4=mb
        m5=0._dp
        npart=3
        call gen3m(r,p,m3,m4,m5,pswt,*999)

c--- processes that use gen3mdk
      elseif ( (kcase==k4ftwdk) .or. (kcase==kdk_4ft) ) then
        m3=mt
        m4=mb
        m5=0._dp
        npart=5
        call gen3mdk(r,p,m3,m4,m5,pswt,*999)

c--- processes that use "gen3m_rap"
      elseif ( (kcase==kvlchm3) ) then
        taumin=(2._dp*mt/sqrts)**2
        m3=mt
        m4=mt
        npart=3
        call gen3m_rap(r,p,m3,m4,pswt,*999)

c--- processes that use "gen4"
      elseif (kcase==kqgtbqq) then
        npart=4
        call gen4(r,p,pswt,*999)

      elseif (kcase==kH_tjet) then
        npart=4
        taumin=((mt+hmass)/sqrts)**2
        call gen4(r,p,pswt,*999)

      elseif (kcase==kZ_tjet) then
        npart=4
        if (zerowidth) then
          taumin=((mt+zmass)/sqrts)**2
        else
          taumin=((mt+sqrt(wsqmin))/sqrts)**2
        endif
        call gen4(r,p,pswt,*999)

      elseif((kcase==kdm_gaj).or.(kcase==kdm2jet)) then
         m3=xmass
         m4=xmass
         npart=4
         call gen4m(r,p,m3,m4,0._dp,0._dp,pswt,*999)

c gg -> VV contributions in NNLO VV calculations
      elseif ( ((kcase == kWWqqbr) .or. (kcase == kZZlept))
     &    .and. (kpart == knnlo) .and. (ipsgen == 2)) then
        npart=4
        call gen4handc(r,p,pswt,*999)

      elseif ( (kcase==kggWW4l) ) then
        npart=4
        if (ipsgen == 1) then
          call gen4(r,p,pswt,*999)
        else
          call gen4handc(r,p,pswt,*999)
        endif
        s3456=2._dp*dot(p,1,2)
        wtips(1)=(s3456-hmass**2)**2/s3456**2
        wtips(2)=one-wtips(1)
        pswt=pswt*wtips(ipsgen)

c--- processes that use "gen4handc"
      elseif ( (kcase==kHZZint)
     &    .or. (kcase==kHZZHpi)
     &    .or. (kcase==kggZZ4l)
     &    .or. (kcase==kHWWint)
     &    .or. (kcase==kHWWHpi)
     &    .or. (kcase==kggWW4l)
     &    .or. (kcase==kHZZ_tb) ) then
!     ) then
        npart=4
c        call gen4_intf(r,p,pswt,*999)
        call gen4handc(r,p,pswt,*999)

c--- processes that use "gen4h"
      elseif ( (kcase==kHWW_4l)
     &    .or. (kcase==kHWWdkW)
     &    .or. (kcase==kHWW2lq)
     &    .or. (kcase==kHWW_tb)
     &    .or. (kcase==kHZZ_4l)
     &    .or. (kcase==kHmZZ4l) ) then
!     &    .or. (kcase==kHZZ_tb) ) then
        npart=4
        call gen4h(r,p,pswt,*999)

c--- processes that use "gen4hvv"
      elseif ( (kcase==kHVV_tb) ) then
        npart=4
        call gen4hvv(r,p,pswt,*999)

c--- processes that use "gen4handcvv"
      elseif ( (kcase==kggVV4l) ) then
        npart=4
        call gen4handcvv(r,p,pswt,*999)

c--- processes that use "gen4vv"
      elseif ( (kcase==kggVVbx) .or. (kcase==kVVlept) ) then
        npart=4
        call gen4vv(r,p,pswt,*999)

c--- processes that use "gen4mdk"
      elseif ((kcase==k4ftjet)
     &   .or. (kcase==kZ_tdkj)
     &   .or. (kcase==kH_tdkj)) then
        npart=6
        call gen4mdk(r,p,pswt,*999)

c--- processes that use "gen5mdk"
      elseif  (kcase==kZtdk2j) then
        npart=7
        call gen5mdk(r,p,pswt,*999)

c--- processes that use "gen5"
      elseif ( (kcase==kW_twdk)
     & .or.    (kcase==kHWWjet)
     & .or.    (kcase==kHZZjet)
     & .or.    (kcase==kWtdkay)
     & .or.    (kcase==kZt2jet)
     & .or.    (kcase==kWbbjem)
     & .or.    (kcase==kvlchwt)
     & .or.    (kcase==kvlchk5)
     & ) then
        npart=5
        call gen5(r,p,pswt,*999)

c--- processes that use "gen5"
      elseif (kcase==kHZZqgI)  then
        npart=5
        call gen5(r,p,pswt,*999)

c--- processes that use "gen6"
      elseif ( (kcase==ktt_bbl)
     &    .or. (kcase==ktt_bbh)
     &    .or. (kcase==ktt_bbu)
     &    .or. (kcase==ktt_ldk)
     &    .or. (kcase==ktt_hdk)
     &    .or. (kcase==ktt_udk)
     &    .or. (kcase==ktthWdk)
     &    .or. (kcase==kWtbwdk)
     &    .or. (kcase==kttZbbl)
     &    .or. (kcase==ktautau)
     &    .or. (kcase==kvlchk6)
     &    .or. (kcase==kvlchwg)
     &    .or. (kcase==kvlchwh)
     &    .or. (kcase==kqq_HWW)
     &    .or. (kcase==kqq_HZZ)
     &    .or. (kcase==kqqWWqq)
     &    .or. (kcase==kqqWWss)
     &    .or. (kcase==kqqWZqq)
     &    .or. (kcase==kqqZZqq)
     &    .or. (kcase==kqqVVqq)
     &    .or. (kcase==kHWW2jt)
     &    .or. (kcase==kHZZ2jt)
     &    .or. (kcase==kWH__ZZ)
     &    .or. (kcase==kWH__WW)
     &    .or. (kcase==kZH__WW)
     &    .or. (kcase==kZH__ZZ)
     &    .or. (kcase==kWpWp2j)
     &    .or. (kcase==kWpmZjj)
     &    .or. (kcase==kWpmZbj)
     &    .or. (kcase==kWpmZbb)
     &    .or. (kcase==kWW2jet)) then
        npart=6
        call gen6(r,p,pswt,*999)

c--- processes that use "gen7"
      elseif ( (kcase==kqq_ttg) ) then
        m3=mt
        m4=mt
        m5=0._dp
        npart=7
        call gen7m(r,p,m3,m4,m5,pswt,*999)

      elseif ( (kcase==kWpWp3j)
     &     .or.(kcase==kHWW3jt)
     &     .or.(kcase==kHZZ3jt) )  then
        npart=7
      call gen7(r,p,pswt,*999)

c--- processes that use "gen8"
      elseif ( (kcase==kqq_tth)
     &    .or. (kcase==kvlchk8) ) then
        npart=8
        call gen8(r,p,pswt,*999)

c--- processes that use "gen8"
      elseif ( (kcase==kqq_ttz)
     &    .or. (kcase==kqqtthz) ) then
        npart=8
c        call gen8(r,p,pswt,*999)
        call genttvdk(r,p,pswt,*999)

      elseif ( (kcase==kqq_ttw)  .or. (kcase==kttwldk) ) then
        npart=8
        taumin=(2._dp*mt/sqrts)**2
        call gen8_rap(r,p,pswt,*999)

      elseif (kcase==ktth_ww) then
        npart=10
        call gen10(r,p,pswt,*999)
c--- processes that use "gen_njets" with an argument of "1"
c     elseif (kcase == kggfus1) then
c         npart = 2
c         call gen2m(r,plo,pswt,*999)
c         allocate(dipconfig(2,3))
c         dipconfig(1,:) = (/ 1,5,2 /)
c         dipconfig(2,:) = (/ 2,5,1 /)
c         npart=npart+1
c         if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),
c    &                  r(ndim+1),plo,p,pswtdip,dipconfig_in=dipconfig)) then
c             goto 999
c         endif
c         pswt=pswt*pswtdip

      ! this seems less efficient than gen_njets, at least for Z+1jet
c     elseif (kcase == kW_1jet .or. kcase == kZ_1jet) then
c         npart = 2
c         call gen2(r,plo,pswt,*999)
c         allocate(dipconfig(2,3))
c         dipconfig(1,:) = (/ 1,5,2 /)
c         dipconfig(2,:) = (/ 2,5,1 /)
c         npart=npart+1
c         if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),
c    &                  r(ndim+1),plo,p,pswtdip,dipconfig_in=dipconfig)) then
c             goto 999
c         endif
c         pswt=pswt*pswtdip

      elseif ( (kcase==kW_1jet)
     &    .or. (kcase==kWcjet0)
     &    .or. (kcase==kZ_1jet)
     &    .or. (kcase==kH_1jet)
     &    .or. (kcase==khttjet)
     &    .or. (kcase==kattjet)
     &    .or. (kcase==kggfus1)
     &    .or. (kcase==khjetma)
     &    .or. (kcase==kHgagaj)
     &    .or. (kcase==kgQ__ZQ) ) then
        npart=3
        if (origkpart == kresummed .or. (usescet .and. useQT)) then
            call gen_njets(r,1,p,pswt,*999)
        else
            if ((usescet) .and. (abovecut)) then
              if ((tauboost .eqv. .false.) .and. (kpart .ne. kreal)) then
                call gen3taucut(r,p,pswt,*999)
              else
                call gen3(r,p,pswt,*999)
              endif
            else
              call gen_njets(r,1,p,pswt,*999)
            endif
        endif

c This code will also work, but is less efficient at benchmark points
c        if     ((usescet) .and. (abovecut)) then
c          call gen3taucut(r,p,pswt,*999)
c        else
c          call gen_njets(r,1,p,pswt,*999)
c        endif

c--- processes that use "gen_njets" with an argument of "2"
      elseif ( (kcase==kWbbbar)
     &    .or. (kcase==kW_2jet)
     &    .or. (kcase==kZ_2jet)
     &    .or. (kcase==kZbbbar)
     &    .or. (kcase==kqq_Hqq)
     &    .or. (kcase==kqq_Hgg)
     &    .or. (kcase==kggfus2)
     &    .or. (kcase==kgagajj)
     &    .or. (kcase==kh2jmas)
     &    .or. (kcase==kW_bjet)
     &    .or. (kcase==kWcjetg)
     &    .or. (kcase==kZ_bjet) ) then
        npart=4
        if ((usescet) .and. (abovecut)) then
c          call genVjtaucut(r,p,pswt,*999)
          call gen_njets(r,2,p,pswt,*999)
        else
          call gen_njets(r,2,p,pswt,*999)
        endif

c--- processes that use "gen_njets" with an argument of "3"
      elseif ( (kcase==kW_3jet)
     &    .or. (kcase==kWbbjet)
     &    .or. (kcase==kZ_3jet)
     &    .or. (kcase==kZbbjet)
     &    .or. (kcase==kWb2jet)
     &    .or. (kcase==kqqHqqg)
     &    .or. (kcase==kggfus3)
     &    .or. (kcase==kZbjetg)) then
        npart=5
        call gen_njets(r,3,p,pswt,*999)
c--- processes that use "gen_stop" with an argument of "1" (number of extra jets)
      elseif (kcase == ktopanom) then
        npart = 4
        !call gen_stop(r,1,p,pswt,*999) ! 718
        call pstop(r,p,pswt,*999)
      elseif ( (kcase==kbq_tpq) .or. (kcase==ktopanom)
     &    .or. (kcase==kttdkay)
     &    .or. (kcase==kt_bbar)
     &    .or. (kcase==ktdecay) ) then
        npart=4
        call gen_stop(r,1,p,pswt,*999)
        !call pstop(r,p,pswt,*999)

c--- WH processes with 2-body Higgs decays
      elseif ( (kcase==kWHbbar) .or. (kcase==kZHbbar)
     &    .or. (kcase==kWHbbdk) .or. (kcase==kZHbbdk)
     &    .or. (kcase==kWHgaga) .or. (kcase==kZHgaga) ) then
        npart=4
        call genVH(r,p,pswt,*999)

c--- WH+jet processes with 2-body Higgs decays
      elseif ( (kcase==kWH1jet)  .or. (kcase==kZH1jet) ) then
        npart=5
        if ((usescet).and.(abovecut).and.(tauboost.eqv. .false.).and.(kpart  /=  kreal)) then
          call genVHjtaucut(r,p,pswt,*999)
        else
          if (hdecaymode == 'wpwm') then
            npart=7
            call gen7(r,p,pswt,*999)
          else
            call gen5(r,p,pswt,*999)
          endif
        endif

c--- VV+jet processes decays
      elseif ((kcase==kWW_jet) .or. (kcase==kWZ_jet) .or.(kcase==kZZ_jet)) then
        npart=5
        if     (((usescet) .and. (abovecut)) .or. (origkpart==kresummed))  then
          call genVVjtaucut(r,p,pswt,*999)
        else
          call gen5(r,p,pswt,*999)
        endif

c--- DEFAULT: processes that use "gen4"
      else
        if ((kcase==kvlchk4) .or. (kcase==kvlchkm)) then
          wsqmin=0._dp
          wsqmax=sqrts**2
        endif
        npart=4
        call gen4(r,p,pswt,*999)
      endif

c--- ensure no entry beyond npart+2
      p(npart+3,:)=0._dp

      return

c--- alternate return
  999 continue
      pswt=zip
      return 1

      end

