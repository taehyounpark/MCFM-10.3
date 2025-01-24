!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine masscuts(p,*)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'mxpart.f'
      include 'limits.f'
      include 'kprocess.f'
      include 'cutoff.f'
      include 'interference.f'
      include 'nqcdjets.f'
      include 'nproc.f'
      include 'masses.f'
      include 'first.f'
      include 'mpicommon.f'
      include 'nwz.f'
      include 'breit.f'
      logical:: VBSprocess
      real(dp):: p(mxpart,4),s34,s56,s36,s45,s3456,s78,pt34sq,sum1,sum2,P1,P2
      real(dp):: pt,ptlZ1,ptlZ2,ptlW
      save VBSprocess
!$omp threadprivate(VBSprocess)

      if (first) then
      first=.false.
c--- do not allow a cut on m56 for W/Z+gamma processes, tau pairs, or DM
      if ( (kcase==kWgamma) .or. (kcase==kZgamma)
     & .or.(kcase==ktautau) .or. (kcase==kdm_jet)
     & .or.(kcase==kdm_gam) ) then
        bbsqmin=0._dp
        bbsqmax=81d8
      endif

c--- do not allow a cut on m34 for direct photon process, or tau pairs
      if ((kcase==kdirgam) .or. (kcase==kgamjet)
     &.or.(kcase==kgam_2j) .or. (kcase==ktautau)) then
        wsqmin=0._dp
        wsqmax=81.e8_dp
      endif

c---- do not allow cuts on m34 if doing gamma gamma (will be done elsewhere)
      if((kcase==kgamgam) .or. (kcase==kgg2gam) .or. (kcase==kgmgmjt)) then
         return
      endif

c---- do not allow any of these cuts if doing gam+2jets
      if(kcase==kgam_2j) then
         return
      endif

c---- check to see whether this is a VBS process
      if ( ((nproc >= 220) .and. (nproc <= 229) .and. (nproc /= 221))
     & .or. (nproc == 2201) .or. (nproc == 2221)
     & .or. (nproc == 2231) .or. (nproc == 2241)
     & .or. (nproc == 2251) .or. (nproc == 2281)
     & .or. (nproc == 2291) ) then
        VBSprocess=.true.
      else
        VBSprocess=.false.
      endif

!$omp master
      if (rank == 0) then
      write(6,*)
      write(6,*) '****************** Basic mass cuts *****************'
      write(6,*) '*                                                  *'
      write(6,99) sqrt(wsqmin),'m34',sqrt(wsqmax)
      if ((nqcdjets < 2) .or. (VBSprocess)) then
      write(6,99) sqrt(bbsqmin),'m56',sqrt(bbsqmax)
      else
      write(6,98) sqrt(bbsqmin),'m(jet1,jet2)',sqrt(bbsqmax)
      endif
      if (interference) then
      write(6,99) sqrt(wsqmin),'m36',sqrt(wsqmax)
      write(6,99) sqrt(wsqmin),'m45',sqrt(wsqmax)
      endif
      write(6,97) m3456min,'m3456',m3456max
      if (VBSprocess) then
       write(6,*)'*               m(jet1,jet2) > 100 GeV             *'
      endif
      if (pt34min > 0._dp) then
       write(6,99) pt34min,'pt34',pt34max
      endif
      write(6,*) '****************************************************'
      endif
!$omp end master
      endif

c for ZZ interference case, only apply cut for leptons closest to Z
      if ((kcase == kZZlept) .and. (interference)) then
        s34=+(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &      -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
        s56=+(p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &      -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2
        s45=+(p(4,4)+p(5,4))**2-(p(4,1)+p(5,1))**2
     &      -(p(4,2)+p(5,2))**2-(p(4,3)+p(5,3))**2
        s36=+(p(3,4)+p(6,4))**2-(p(3,1)+p(6,1))**2
     &      -(p(3,2)+p(6,2))**2-(p(3,3)+p(6,3))**2
        if (closestZ) then
! Find Z-pair assignment with one pair closest to Z mass
          sum1=min((s34-zmass**2)**2,(s56-zmass**2)**2)
          sum2=min((s45-zmass**2)**2,(s36-zmass**2)**2)
        else
! Assign according to sum of differences from Z mass, c.f. 1711.06631
          sum1=(s34-zmass**2)**2+(s56-zmass**2)**2
          sum2=(s45-zmass**2)**2+(s36-zmass**2)**2
        endif
        if (sum1 < sum2) then
          if ((s34 < bbsqmin).or.(s34>bbsqmax)) return 1
          if ((s56 < bbsqmin).or.(s56>bbsqmax)) return 1
        else
          if ((s45 < bbsqmin).or.(s45>bbsqmax)) return 1
          if ((s36 < bbsqmin).or.(s36>bbsqmax)) return 1
        endif
        if (s34 < moppmin**2) return 1
        if (s56 < moppmin**2) return 1
        if (s45 < moppmin**2) return 1
        if (s36 < moppmin**2) return 1
        goto 44
      endif

       if ((kcase == kWZbbar) .and. (interference)) then
        s34=+(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &      -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
        s56=+(p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &      -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2
        s45=+(p(4,4)+p(5,4))**2-(p(4,1)+p(5,1))**2
     &      -(p(4,2)+p(5,2))**2-(p(4,3)+p(5,3))**2
        s36=+(p(3,4)+p(6,4))**2-(p(3,1)+p(6,1))**2
     &      -(p(3,2)+p(6,2))**2-(p(3,3)+p(6,3))**2
        if (closestZ) then
! WZ interference analyis per 2110.11231
          if (nwz == +1) then
            if (abs(s56-zmass**2) < abs(s45-zmass**2)) then
              if ((s56 < bbsqmin).or.(s56>bbsqmax)) return 1
              ptlZ1=max(pt(5,p),pt(6,p))
              ptlZ2=min(pt(5,p),pt(6,p))
              ptlW=pt(4,p)
            else
              if ((s45 < bbsqmin).or.(s45>bbsqmax)) return 1
              ptlZ1=max(pt(5,p),pt(4,p))
              ptlZ2=min(pt(5,p),pt(4,p))
              ptlW=pt(6,p)
            endif
            if (s56 < moppmin**2) return 1
            if (s45 < moppmin**2) return 1
          else
            if (abs(s56-zmass**2) < abs(s36-zmass**2)) then
              if ((s56 < bbsqmin).or.(s56>bbsqmax)) return 1
              ptlZ1=max(pt(5,p),pt(6,p))
              ptlZ2=min(pt(5,p),pt(6,p))
              ptlW=pt(3,p)
            else
              if ((s36 < bbsqmin).or.(s36>bbsqmax)) return 1
              ptlZ1=max(pt(3,p),pt(6,p))
              ptlZ2=min(pt(3,p),pt(6,p))
              ptlW=pt(5,p)
            endif
            if (s56 < moppmin**2) return 1
            if (s36 < moppmin**2) return 1
          endif
          if (ptlZ1 < 25._dp) return 1
          if (ptlZ2 < 10._dp) return 1
          if (ptlW < 25._dp) return 1
        else
c--- for WZ interference case, make cut based on ATLAS resonant-shape procedure
          P1=abs(s34-wmass**2+im*wmass*wwidth)**2*abs(s56-zmass**2+im*zmass*zwidth)**2
          if (nwz == +1) then
            P2=abs(s36-wmass**2+im*wmass*wwidth)**2*abs(s45-zmass**2+im*zmass*zwidth)**2
            if (P2 < P1) then
              if ((s45 < bbsqmin).or.(s45>bbsqmax)) return 1
            else
              if ((s56 < bbsqmin).or.(s56>bbsqmax)) return 1
            endif
          else
            P2=abs(s45-wmass**2+im*wmass*wwidth)**2*abs(s36-zmass**2+im*zmass*zwidth)**2
            if (P2 < P1) then
              if ((s36 < bbsqmin).or.(s36>bbsqmax)) return 1
            else
              if ((s56 < bbsqmin).or.(s56>bbsqmax)) return 1
            endif
          endif
        endif
        goto 44
      endif

c--- only apply cuts on s34 if vectors 3 and 4 are defined
      if ((abs(p(3,4)) > 1.e-8_dp) .and. (abs(p(4,4)) > 1.e-8_dp)) then
        s34=+(p(3,4)+p(4,4))**2-(p(3,1)+p(4,1))**2
     &      -(p(3,2)+p(4,2))**2-(p(3,3)+p(4,3))**2
c--- do not accept s34<cutoff either
        if ((s34 < max(wsqmin,cutoff)).or.(s34 > wsqmax)) return 1
c--- optional cut on pt34
        pt34sq=(p(3,1)+p(4,1))**2+(p(3,2)+p(4,2))**2
        if (pt34sq < pt34min**2) return 1
        if (pt34sq > pt34max**2) return 1
      endif

c---- do not allow any further cuts if doing W,Z,H+2jets
      if((kcase==kW_2jet) .or. (kcase==kZ_2jet) .or. (kcase==kggfus2) .or. (kcase==kgagajj)) then
         return
      endif

c--- only apply cuts on s56 if vectors 5 and 6 are defined and generated with BW
      if ((abs(p(5,4)) > 1.e-8_dp) .and. (abs(p(6,4)) > 1.e-8_dp)
     &    .and. (n2 == 1)) then
        s56=+(p(5,4)+p(6,4))**2-(p(5,1)+p(6,1))**2
     &      -(p(5,2)+p(6,2))**2-(p(5,3)+p(6,3))**2
c--- do not accept s56<cutoff either
c        if ((s56 < max(bbsqmin,cutoff)).or.(s56>bbsqmax)) return 1
        if ((s56 < bbsqmin).or.(s56>bbsqmax)) return 1
      endif

c      if (interference) then
c      s45=+(p(4,4)+p(5,4))**2-(p(4,1)+p(5,1))**2
c     &    -(p(4,2)+p(5,2))**2-(p(4,3)+p(5,3))**2
c      s36=+(p(3,4)+p(6,4))**2-(p(3,1)+p(6,1))**2
c     &    -(p(3,2)+p(6,2))**2-(p(3,3)+p(6,3))**2
c        if (wsqmin  /=  bbsqmin) then
c        write(6,*) 'masscuts: min. cuts must be equal for interference'
c        stop
c        endif
c        if ((s45 < bbsqmin) .or. (s45 > bbsqmax)) return 1
c        if ((s36 < bbsqmin) .or. (s36 > bbsqmax)) return 1
c      endif

   44 continue

c--- only apply cuts on s3456 if vectors 3, 4, 5 and 6 are defined
      if ((abs(p(3,4)) > 1.e-8_dp) .and. (abs(p(4,4)) > 1.e-8_dp)
     &    .and. (abs(p(5,4)) > 1.e-8_dp) .and. (abs(p(6,4)) > 1.e-8_dp))
     &     then

      s3456=+(p(3,4)+p(4,4)+p(5,4)+p(6,4))**2
     &      -(p(3,1)+p(4,1)+p(5,1)+p(6,1))**2
     &      -(p(3,2)+p(4,2)+p(5,2)+p(6,2))**2
     &      -(p(3,3)+p(4,3)+p(5,3)+p(6,3))**2

      if (s3456 < m3456min**2) return 1
      if (s3456 > m3456max**2) return 1

      endif

      if (VBSprocess) then
        s78=+(p(7,4)+p(8,4))**2-(p(7,1)+p(8,1))**2
     &      -(p(7,2)+p(8,2))**2-(p(7,3)+p(8,3))**2
        if (s78 < 1d4) return 1
      endif

   96 format(' *                       ',a5,' > ',f8.2,'           *')
   97 format(' *          ',f8.2,'  <  ',a5,' < ',f8.2,'           *')
   98 format(' *      ',f8.2,'  <   ',a12,'  < ',f8.2,'      *')
   99 format(' *          ',f8.2,'  <   ',a3,'  < ',f8.2,'           *')

      return
      end

