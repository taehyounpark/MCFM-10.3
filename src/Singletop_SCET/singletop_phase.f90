!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
 
module singletop_phase
      implicit none

      public :: gen_singletop

      private

      private :: gen_stop_new, gen_stop_decay

      contains

      function gen_singletop(r,p,wt)
          use types
          use Multichannel
          use singletop2_nnlo_vars
          implicit none
          include 'mxpart.f'
          include 'nproc.f'
          include 'kpart.f'
          include 'taucut.f'! usescet
          include 'npart.f'
          include 'vegas_common.f'! ndim
          include 'ipsgen.f'

          logical :: gen_singletop
          real(dp), intent(in) :: r(mxdim)
          real(dp), intent(out) :: p(mxpart,4)
          real(dp), intent(out) :: wt

          real(dp) :: wtdip, plo(mxpart,4), pone(mxpart,4)
          real(dp), allocatable :: dipconfig(:,:)

          real(dp) :: beamsign

          gen_singletop = .false.

#define ptgen_sqrt 1
#define ptgen_linear 2
#define ptgen_sqrt_invsq 3
    
          if (nproc == 1610) then
                  if (currentContrib == 4) then
                      if (kpart == kreal) then
                          if (ipsgen == 4) then
                              npart = 6
                              call gen_stop_new(r,3,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

!                             allocate(dipconfig(1,3))
!                             dipconfig(1,:) = [1,6,2]

!                             npart=npart+1
!                             if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
!                                           plo,p,wtdip, dipconfig_in=dipconfig)) then
!                                 gen_singletop = .false.
!                                 return
!                             endif
!                             wt=wt*wtdip

                              gen_singletop = .true.
                              return
                          elseif (ipsgen == 5) then
                              npart = 5
                              call gen_stop_new(r,2,ptgen_sqrt,p,wt)

                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              ! does this interfere with slicing?
!                             allocate(dipconfig(6,3))
!                             dipconfig(1,:) = [6,7,1]
!                             dipconfig(2,:) = [1,7,6]
!                             dipconfig(3,:) = [1,6,7]
!                             dipconfig(4,:) = [6,7,2]
!                             dipconfig(5,:) = [2,7,6]
!                             dipconfig(6,:) = [2,6,7]

!                             npart=npart+1
!                             if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
!                                           plo,p,wtdip, dipconfig_in=dipconfig)) then
!                                 gen_singletop = .false.
!                                 return
!                             endif
!                             wt=wt*wtdip

                              gen_singletop = .true.
                              return
                          endif
                      elseif (kpart == kvirt) then
                          if (ipsgen == 4) then
                              npart = 5
                              call gen_stop_new(r,2,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              gen_singletop = .true.
                              return
                          elseif (ipsgen == 5) then
                              npart = 4
                              call gen_stop_new(r,1,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              gen_singletop = .true.
                              return
                          endif
                      endif
                  elseif (currentContrib == 5) then
                      if (kpart == kvirt) then
                          if (ipsgen == 7) then ! VV
                              npart = 4
                              call gen_stop_new(r,1,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              gen_singletop = .true.
                              return
                          elseif (ipsgen == 6) then ! VR
                              npart = 5
                              call gen_stop_decay(r,1,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              gen_singletop = .true.
                              return
                          endif
                      elseif (kpart == kreal) then
                          if (ipsgen == 6) then ! RR
! LABORDAY ADDED
#define DEBUGTHIS 0
#if (DEBUGTHIS == 1)
                              npart = 5
                              call gen_stop_decay(r,1,ptgen_sqrt,plo,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              allocate(dipconfig(3,3))
                              ! adjust channel here
                              dipconfig(1,:) = [1,6,8]
                              dipconfig(2,:) = [6,8,1]
                              dipconfig(3,:) = [1,8,6]

                              npart=npart+1
                              if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                            plo,p,wtdip, dipconfig_in=dipconfig)) then
                                  gen_singletop = .false.
                                  return
                              endif
                              wt=wt*wtdip
#else
                              npart = 6
                              call gen_stop_decay(r,2,ptgen_sqrt,p,wt)

                              p(9,:) = p(8,:)
                              p(8,:) = p(7,:)
                              p(7,:) = p(9,:)
                              p(9,:) = 0._dp

                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif
#endif

                              gen_singletop = .true.
                              return
                          elseif (ipsgen == 7) then ! RV
                              npart = 5
                              call gen_stop_new(r,2,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

!                             allocate(dipconfig(1,3))
!                             dipconfig(1,:) = [6,7,1]

!                             npart=npart+1
!                             if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
!                                           plo,p,wtdip, dipconfig_in=dipconfig)) then
!                                 gen_singletop = .false.
!                                 return
!                             endif
!                             wt=wt*wtdip

                              gen_singletop = .true.
                              return
                          endif
                      endif
                  elseif (currentContrib == 6) then
                      if (kpart == kvirt) then
                          if (ipsgen == 8) then ! VR
                              npart = 5
                              call gen_stop_decay(r,1,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              gen_singletop = .true.
                              return
                          elseif (ipsgen == 9) then ! VV
                              npart = 4
                              call gen_stop_new(r,1,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              gen_singletop = .true.
                              return
                          endif
                      elseif (kpart == kreal) then
                          if (ipsgen == 8) then ! RR
#define DEBUGTHIS 0
#if (DEBUGTHIS == 1)
                              npart = 5
                              call gen_stop_decay(r,1,ptgen_sqrt,plo,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              allocate(dipconfig(1,3))
                              ! adjust channel here
                              dipconfig(1,:) = [2,8,1]

                              npart=npart+1
                              if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                            plo,p,wtdip, dipconfig_in=dipconfig)) then
                                  gen_singletop = .false.
                                  return
                              endif
                              wt=wt*wtdip
#else
                              npart = 6
                              call gen_stop_decay(r,2,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

                              p(9,:) = p(8,:)
                              p(8,:) = p(7,:)
                              p(7,:) = p(9,:)
                              p(9,:) = 0._dp
#endif

                              gen_singletop = .true.
                              return
                          elseif (ipsgen == 9) then ! RV
                              npart = 5
                              call gen_stop_new(r,2,ptgen_sqrt,p,wt)
                              if (wt == 0._dp) then
                                  gen_singletop = .false.
                                  return
                              endif

!                             allocate(dipconfig(1,3))
!                             dipconfig(1,:) = [6,7,1]

!                             npart=npart+1
!                             if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
!                                           plo,p,wtdip, dipconfig_in=dipconfig)) then
!                                 gen_singletop = .false.
!                                 return
!                             endif
!                             wt=wt*wtdip


                              gen_singletop = .true.
                              return
                          endif
                      endif
                  elseif ((kpart == klord .or. kpart == kvirt .or. kpart == ksnlo &
                        .or. kpart == ksnloV .or. kpart == knnlo) .and. &
                      (currentContrib == 1 .or. currentContrib == 2)) then
                      npart = 4
                      call gen_stop_new(r,1,ptgen_sqrt,p,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif

                      gen_singletop = .true.
                      return
                  elseif ((kpart == klord .or. kpart == kvirt .or. kpart == ksnlo &
                          .or. kpart == ksnloV .or. kpart == knnlo ) .and. &
                          (currentContrib == 3)) then
                      npart = 4
                      call gen_stop_new(r,1,ptgen_sqrt,p,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif

                      gen_singletop = .true.
                      return
                  elseif ((kpart == kreal) .and. (currentContrib == 2)) then
                      if (maxbeams == 2) then
                          npart = 5
                          call gen_stop_new(r,2,ptgen_sqrt,p,wt)
                          if (wt == 0._dp) then
                              gen_singletop = .false.
                              return
                          endif
                      else
                          ! to tweak
                          error stop "to tweak"
                          npart = 4
                          call gen_stop_new(r,1,ptgen_sqrt,plo,wt)
                          if (wt == 0._dp) then
                              gen_singletop = .false.
                              return
                          endif

                          allocate(dipconfig(2,3))
                          if (beams_enabled(1) == 1) then
                              dipconfig(1,:) = [1,7,2]
                          else
                              dipconfig(1,:) = [2,7,1]
                          endif

                          npart=npart+1
                          if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                        plo,p,wtdip, dipconfig_in=dipconfig)) then
                              gen_singletop = .false.
                              return
                          endif
                          wt=wt*wtdip
                      endif

                      gen_singletop = .true.
                      return
                  elseif ((kpart == kreal) .and. (currentContrib == 1)) then
                      if (maxbeams == 2) then
                          npart = 5
                          call gen_stop_new(r,2,ptgen_sqrt,p,wt)
                          if (wt == 0._dp) then
                              gen_singletop = .false.
                              return
                          endif
                      else
                          npart = 4
                          call gen_stop_new(r,1,ptgen_sqrt,plo,wt)
                          if (wt == 0._dp) then
                              gen_singletop = .false.
                              return
                          endif

                          allocate(dipconfig(2,3))
                          if (beams_enabled(1) == 1) then
                              dipconfig(1,:) = [6,7,1]
                          else
                              dipconfig(1,:) = [6,7,2]
                          endif

                          npart=npart+1
                          if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                        plo,p,wtdip, dipconfig_in=dipconfig)) then
                              gen_singletop = .false.
                              return
                          endif
                          wt=wt*wtdip
                      endif

                      gen_singletop = .true.
                      return
                  elseif ((kpart == kreal) .and. (currentContrib == 3)) then
                      npart = 5
                      call gen_stop_decay(r,1,ptgen_sqrt,p,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif

                      gen_singletop = .true.
                      return
                  else
                      error stop "to do"
                  endif

 !            endif

          elseif (nproc == 1650) then
              if ((kpart == klord) .and. (origkpart == ksnlo) .and. (currentContrib == 2)) then
                  if (maxbeams == 2) then
                      npart = 5
                      call gen_stop_new(r,2,ptgen_sqrt,p,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif
                  else
                      npart = 4
                      call gen_stop_new(r,1,ptgen_sqrt,plo,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif
 
                      allocate(dipconfig(1,3))
                      if (beams_enabled(1) == 1) then
                          beamsign = 1
                          dipconfig(1,:) = [1,7,2]
                      else
                          beamsign = -1
                          dipconfig(1,:) = [2,7,1]
                      endif
 
                      npart=npart+1
                      if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                    plo,p,wtdip, dipconfig_in=dipconfig)) then
                          gen_singletop = .false.
                          return
                      endif
                      wt=wt*wtdip
                  endif

                  gen_singletop = .true.
                  return

              elseif ((kpart == klord) .and. (origkpart == ksnlo) .and. (currentContrib == 1)) then
                  if (maxbeams == 2) then
                      npart = 5
                      call gen_stop_new(r,2,ptgen_sqrt,p,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif
                  else
                      npart = 4
                      call gen_stop_new(r,1,ptgen_sqrt,plo,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif
 
                      allocate(dipconfig(2,3))
                      if (beams_enabled(1) == 1) then
                          beamsign = 1
                          dipconfig(1,:) = [6,7,1]
                          dipconfig(2,:) = [1,7,6]
                      else
                          beamsign = -1
                          dipconfig(1,:) = [6,7,2]
                          dipconfig(2,:) = [2,7,6]
                      endif
 
                      npart=npart+1
                      if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                    plo,p,wtdip, dipconfig_in=dipconfig)) then
                          gen_singletop = .false.
                          return
                      endif
                      wt=wt*wtdip

                      if (taups > 0) then
#define Energy(i) (p(i,4))
#define pzcomp(i) (beamsign*p(i,3))
                  tauphase: block
                      real(dp) :: taubeam, taujet, p67(4)
                      real(dp) :: wtjet, wtbeam

                      p67(:) = p(6,:) + p(7,:)
                      taubeam = min(Energy(7)-pzcomp(7), Energy(6)-pzcomp(6))
                      taujet = Energy(6)+Energy(7) - sqrt(p67(1)**2+p67(2)**2+p67(3)**2)

                      wtjet = 1._dp/taujet
                      wtbeam = 1._dp/taubeam

                      if (taups == 1) then
                          wt = wt * wtbeam/(wtbeam+wtjet)
                      elseif (taups == 2) then
                          wt = wt * wtjet/(wtbeam+wtjet)
                      else
                          continue
                      endif
                  end block tauphase
#undef Energy
#undef pzcomp
                      endif
                  endif

                  gen_singletop = .true.
                  return
              elseif ((kpart == klord .or. kpart == kvirt) .and. (currentContrib == 1)) then
                  if (maxbeams == 2) then
                      npart = 5
                      call gen_stop_new(r,2,ptgen_sqrt,p,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif
                  else
                      npart = 4
                      call gen_stop_new(r,1,ptgen_sqrt,plo,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif

                      if (beams_enabled(1) == 1) then
                          allocate(dipconfig(2,3))
                          dipconfig(1,:) = [6,7,1]
                          dipconfig(2,:) = [1,7,6]
                      else
                          allocate(dipconfig(2,3))
                          dipconfig(1,:) = [6,7,2]
                          dipconfig(2,:) = [2,7,6]
                      endif

                      npart=npart+1
                      if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                    plo,p,wtdip, dipconfig_in=dipconfig)) then
                          gen_singletop = .false.
                          return
                      endif
                      wt=wt*wtdip
                  endif

                  gen_singletop = .true.
                  return
              elseif ((kpart == klord .or. kpart == kvirt) .and. (currentContrib == 2)) then
                  ! this seems to be better than attaching a dipole splitting
                  npart = 5
                  call gen_stop_new(r,2,ptgen_sqrt,p,wt)
                  if (wt == 0._dp) then
                      gen_singletop = .false.
                      return
                  endif

                  gen_singletop = .true.
                  return
              elseif (kpart == kreal .and. origkpart == knnlo .and. currentContrib == 1) then
#define REALPS 2
#if (REALPS == 1)
                  npart = 6
                  call gen_stop_new(r,3,ptgen_sqrt_invsq,p,wt)
                  if (wt == 0._dp) then
                      gen_singletop = .false.
                      return
                  endif
#else
                  npart = 4
                  call gen_stop_new(r,1,ptgen_sqrt,plo,wt)
                  if (wt == 0._dp) then
                      gen_singletop = .false.
                      return
                  endif

                  if (taups > 0) then
                      if (beams_enabled(1) == 1) then
                          beamsign = 1d0
                          allocate(dipconfig(1,3))
                          dipconfig(1,:) = [6,7,1]
                      else
                          beamsign = -1d0
                          allocate(dipconfig(1,3))
                          dipconfig(1,:) = [6,7,2]
                      endif
                  else
                      allocate(dipconfig(1,3))
                      dipconfig(1,:) = [6,7,1]
                      !dipconfig(2,:) = [6,7,2]
                  endif

                  npart=npart+1
                  if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),r(ndim+1), &
                                plo,pone,wtdip, dipconfig_in=dipconfig)) then
                      gen_singletop = .false.
                      return
                  endif
                  wt=wt*wtdip

                  deallocate(dipconfig)
                  if (taups > 0) then
                      allocate(dipconfig(1,3))
                      if (iand(quarkChannel, partons_enabled) > 0) then
                          ! best for quark channel
                          if (beams_enabled(1) == 1) then
                              dipconfig(1,:) = [1,8,7]
                          else
                              dipconfig(1,:) = [2,8,7]
                          endif
                      else
                          ! best for gluon channel
                          dipconfig(1,:) = [7,8,6]
                      endif
                  else
                      allocate(dipconfig(3,3))
                      dipconfig(1,:) = [1,8,7]
                      dipconfig(2,:) = [2,8,7]
                      dipconfig(3,:) = [7,8,6]
                  endif

                  npart=npart+1
                  if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                pone,p,wtdip, dipconfig_in=dipconfig)) then
                      gen_singletop = .false.
                      return
                  endif
                  wt=wt*wtdip
#endif
#undef REALPS
                  if (taups > 0) then
#define E(i) (p(i,4))
#define pz(i) (beamsign*p(i,3))
                  tauphasennlo: block
                      real(dp) :: tauallbeam,taubeamjet,taualljet
                      real(dp) :: p67(4), p68(4), p78(4), p678(4)
                      real(dp) :: wtallbeam,wtbeamjet,wtalljet

                      p67(:) = p(6,:) + p(7,:)
                      p68(:) = p(6,:) + p(8,:)
                      p78(:) = p(7,:) + p(8,:)
                      p678(:) = p(6,:) + p(7,:) + p(8,:)

                      tauallbeam = min(E(6)+E(7)-pz(6)-pz(7), E(6)+E(8)-pz(6)-pz(8), E(7)+E(8)-pz(7)-pz(8))
                      taubeamjet = min(E(6)+E(7)+E(8) - sqrt(p67(1)**2+p67(2)**2+p67(3)**2) - pz(8), &
                          E(6)+E(7)+E(8) - sqrt(p68(1)**2+p68(2)**2+p68(3)**2) - pz(7), &
                          E(6)+E(7)+E(8) - sqrt(p78(1)**2+p78(2)**2+p78(3)**2) - pz(6))
                      taualljet = E(6)+E(7)+E(8) - sqrt(p678(1)**2+p678(2)**2+p678(3)**2)

                      wtallbeam = 1._dp/tauallbeam
                      wtbeamjet = 1._dp/taubeamjet
                      wtalljet = 1._dp/taualljet

                      if (taups == 1) then
                          wt = wt * wtallbeam/(wtalljet+wtallbeam+wtbeamjet)
                      elseif (taups == 2) then
                          wt = wt * wtbeamjet/(wtalljet+wtallbeam+wtbeamjet)
                      elseif (taups == 3) then
                          wt = wt * wtalljet/(wtalljet+wtallbeam+wtbeamjet)
                      else
                          continue
                          !error stop "wrong taups in phase space generation"
                      endif

                  end block tauphasennlo
#undef E
#undef pz
                  endif

                  gen_singletop = .true.
                  return

              elseif (kpart == kreal .and. currentContrib == 2) then
#define REALPS 1
#if (REALPS == 2)
                      npart = 6
                      call gen_stop_new(r,3,ptgen_sqrt,p,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif
#else
                      npart = 4
                      call gen_stop_new(r,1,ptgen_sqrt,plo,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif

                      allocate(dipconfig(1,3))
                      if (beams_enabled(1) == 1) then
                          dipconfig(1,:) = [1,7,2]
                      else
                          dipconfig(1,:) = [2,7,1]
                      endif

                      npart=npart+1
                      if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),r(ndim+1), &
                                    plo,pone,wtdip, dipconfig_in=dipconfig)) then
                          gen_singletop = .false.
                          return
                      endif
                      wt=wt*wtdip

                      deallocate(dipconfig)
                      allocate(dipconfig(2,3))
                      if (beams_enabled(1) == 1) then
                          dipconfig(1,:) = [1,8,7]
                          dipconfig(2,:) = [7,8,1]
                      else
                          dipconfig(1,:) = [2,8,7]
                          dipconfig(2,:) = [7,8,2]
                      endif

                      npart=npart+1
                      if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                    pone,p,wtdip, dipconfig_in=dipconfig)) then
                          gen_singletop = .false.
                          return
                      endif
                      wt=wt*wtdip
#endif
#undef REALPS

                  gen_singletop = .true.
                  return
              elseif (kpart == kreal .and. currentContrib == 1) then
                  if (maxbeams == 2) then
                      npart = 6
                      call gen_stop_new(r,3,ptgen_sqrt,p,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif
                  else
                      npart = 4
                      call gen_stop_new(r,1,ptgen_sqrt,plo,wt)
                      if (wt == 0._dp) then
                          gen_singletop = .false.
                          return
                      endif

                      allocate(dipconfig(2,3))
                      if (beams_enabled(1) == 1) then
                          dipconfig(1,:) = [6,7,1]
                          dipconfig(2,:) = [1,7,6]
                      else
                          dipconfig(1,:) = [6,7,2]
                          dipconfig(2,:) = [2,7,6]
                      endif

                      npart=npart+1
                      if (.not. multichan(r(ndim-5),r(ndim-4),r(ndim-3),r(ndim+1), &
                                    plo,pone,wtdip, dipconfig_in=dipconfig)) then
                          gen_singletop = .false.
                          return
                      endif
                      wt=wt*wtdip

                      deallocate(dipconfig)
                      allocate(dipconfig(2,3))
                      if (beams_enabled(1) == 1) then
                          dipconfig(1,:) = [1,8,7]
                          dipconfig(2,:) = [7,8,6]
                      else
                          dipconfig(1,:) = [2,8,7]
                          dipconfig(2,:) = [7,8,6]
                      endif

                      npart=npart+1
                      if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                    pone,p,wtdip, dipconfig_in=dipconfig)) then
                          gen_singletop = .false.
                          return
                      endif
                      wt=wt*wtdip
                  endif

                  gen_singletop = .true.
                  return
              elseif (((kpart == klord) .or. (kpart == kvirt)) .and. (currentContrib == 3)) then
                  npart = 5
                  call gen_stop_decay(r,1,ptgen_sqrt,p,wt)
                  if (wt == 0._dp) then
                      gen_singletop = .false.
                      return
                  endif

                  gen_singletop = .true.
                  return
              elseif ((kpart == kreal) .and. (currentContrib == 3)) then
#define REALPS 1
#if (REALPS == 1)
                  npart = 5
                  call gen_stop_decay(r,1,ptgen_sqrt,plo,wt)
                  if (wt == 0._dp) then
                      gen_singletop = .false.
                      return
                  endif

                  npart=npart+1
                  allocate(dipconfig(4,3))
                  dipconfig(1,:) = [7,8,5]
                  dipconfig(2,:) = [8,7,5]
                  dipconfig(3,:) = [5,8,7]
                  dipconfig(4,:) = [5,7,8]

                  if (.not. multichan(r(ndim-2),r(ndim-1),r(ndim),r(ndim+2), &
                                plo,p,wtdip, dipconfig_in=dipconfig)) then
                      gen_singletop = .false.
                      return
                  endif
                  wt=wt*wtdip
#else
                 error stop "todo: start from gen_stop_new and radiate twice in decay"
                 ! this would require multichan for massive spectators
#endif


                  gen_singletop = .true.
                  return
              else
                  write (*,*) "KPART", kpart, origkpart, currentContrib, kreal, knnlo
                  error stop "to do"
              endif

          else
              error stop "gen_singletop is only for nprocbelow 1610 and 1650"
          endif

      end function

    subroutine gen_stop_new(r,njets,ptgen,p,wt)
    implicit none
    include 'types.f'
          
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    include 'masses.f'
    include 'mxdim.f'
    include 'limits.f'
    include 'zerowidth.f'
    include 'kprocess.f'
    include 'reset.f'
    include 'kpart.f'
    include 'x1x2.f'
    include 'energy.f'
    include 'notag.f'
    include 'taucut.f'
    include 'tiny.f'
    include 'debug.f'
!---- Generate phase space for 2-->2+n process
!---- with (345) being a top and 6,..,5+n the jets
!---- r(mxdim),p1(4),p2(4) are inputs reversed in sign
!---- from physical values
!---- phase space for -p1-p2 --> p3+p4+p5+p6
!---- with all 2 pi's (ie 1/(2*pi)^(3n-4), where n is the number
!----  of final state particles)
!---- This routine has a minimum of 4 final state particles, hence
!---- the twopi**2 correction factor is given by the ratio of
!---- (1/twopi)**(3n-4) present in the phase space and the factor
!---- of [(1/twopi)**2]**(n-1) from the number of branchings
!---- For the specific case 'ttdkay' where one of the jets is
!---- associated with the top quark decay, we must add an extra
!---- factor of (1/twopi) since the number of jets generated is
!---- larger than the value of 'njets' passed
    real(dp), intent(in) :: r(mxdim)
    integer, intent(in) :: njets
    integer, intent(in) :: ptgen
    real(dp), intent(out) :: p(mxpart,4)
    real(dp), intent(out) :: wt

    real(dp):: psumjet(4),pcm(4),Q(4)
    real(dp):: wt0
    real(dp):: pt,etamax,etamin
    real(dp):: y,sinhy,coshy,phi,mv2,wtbw,mjets
    real(dp):: ybar,ptsumjet2,ycm,sumpst,q0st,rshat,dely,xjac
    real(dp):: plstar,estar,plstarsq,y5starmax,y5starmin
    real(dp):: bm(4),wp(4),nn(4),ep(4),wtwp,wtepnn
    integer:: j,nu,ijet
    parameter(wt0=1._dp/twopi**2)


          
    p(1:8,:) = 0._dp
    pcm = 0._dp
    psumjet = 0._dp
    wt=2._dp*pi
                
    do ijet=1,njets
    !--- generate the pt of jet number ijet
    !--- rapidity limited by E=pT*coshy
        wt=wt/16._dp/pi**3

    !       ! some other ways to generate in pT
    !       ! directly in pT
    !       pt = sqrts*r(ijet)/2._dp
    !       xjac = r(ijet)*sqrts**2/4._dp
    !       ! 1/pt
    !       pt = sqrts*r(ijet)/(sqrts+2*r(ijet)-sqrts*r(ijet))
    !       xjac = r(ijet)*sqrts**3 / (2*r(ijet) - sqrts*(r(ijet)-1))**3
    !       ! -log(pT)
    !       pt = exp(log(sqrts/2._dp) - (1-r(ijet))/r(ijet))
    !       xjac = exp(2._dp*log(sqrts/2._dp) - 2._dp*(1-r(ijet))/r(ijet))/r(ijet)**2
        if (ptgen == ptgen_sqrt) then
            pt = r(ijet)**2*sqrts/2._dp
            xjac = pt*sqrts/2._dp*2._dp*r(ijet)
        elseif (ptgen == ptgen_linear) then
            pt = sqrts*r(ijet)/2._dp
            xjac = r(ijet)*sqrts**2/4._dp
        elseif (ptgen == ptgen_sqrt_invsq) then
            if (ijet == 1) then
                pt = r(ijet)**2*sqrts/2._dp
                xjac = pt*sqrts/2._dp*2._dp*r(ijet)
            else
                pt = 1._dp / sqrt(4._dp/sqrts**2 + (1-r(ijet))/r(ijet))
                xjac = sqrts**4 / 2._dp / (sqrts**2*(r(ijet)-1._dp) - 4*r(ijet))**2
            endif
        else
            pt = r(ijet)**2*sqrts/2._dp
            xjac = pt*sqrts/2._dp*2._dp*r(ijet)
        endif

    ! for nnlo double real currentContrib=1
!       if ((kcase == kbq_tpq_jet) .AND. (kpart == kreal) .AND. usescet .AND. (currentContrib == 1)) then
!           if (ijet == 1) then
!           ! in 1/sqrt(pT)
!               pt = r(ijet)**2*sqrts/2._dp
!               xjac = pt*sqrts/2._dp*2._dp*r(ijet)
!           else
!           ! 1/pt^2
!               pt = 1._dp / sqrt(4._dp/sqrts**2 + (1-r(ijet))/r(ijet))
!               xjac = sqrts**4 / 2._dp / (sqrts**2*(r(ijet)-1._dp) - 4*r(ijet))**2
!           endif
!       ! for nlo real, currentContrib == 2
!       elseif ((kcase == kbq_tpq) .AND. (kpart == kreal) .AND. ( .NOT. usescet) .AND. (currentContrib == 2)) then
!           error stop "todo"
!       else
!           error stop "please explicitly setup gen_stop for this case"
!       ! all genpt(r(ijet),ptbreak,.true.,pt,xjac)
!       endif

        wt=wt*xjac

        etamax=sqrts/2._dp/pt
        if (etamax**2 <= 1._dp) then
                if (debug) then
                    write(6,*) 'etamax**2 <= 1._dp in gen_stop.f',etamax**2
                endif
            wt=0._dp
            return
        endif

        etamax=log(etamax+sqrt(etamax**2-1._dp))
        etamax=min(etamax,100d0)
        y=etamax*(2._dp*r(njets+ijet)-1._dp)
        wt=wt*2._dp*etamax
                
        sinhy=sinh(y)
        coshy=sqrt(1._dp+sinhy**2)
                
        p(5+ijet,4)=pt*coshy
                
        phi=2._dp*pi*r(2*njets+ijet)
        wt=wt*2._dp*pi
                
        p(5+ijet,1)=pt*cos(phi)
        p(5+ijet,2)=pt*sin(phi)
        p(5+ijet,3)=pt*sinhy
                
        do nu=1,4
            psumjet(nu)=psumjet(nu)+p(5+ijet,nu)
        enddo
    enddo
          
    ! top generated zero width
      mv2 = mt**2
      wtbw = pi*mt*twidth

    wt=wt*wtbw
!--- for one jet, mjets must be exactly zero
    if (njets == 1) then
        mjets=0._dp
    else
    !--- invariant mass of jets
        mjets=psumjet(4)**2-psumjet(1)**2-psumjet(2)**2-psumjet(3)**2
    !--- check that argument of upcoming sqrt is not negative
        if (mjets < 0._dp) then
            wt=0._dp
            return
        endif
        mjets=sqrt(mjets)
    endif
          
    if (psumjet(4)-psumjet(3) == 0._dp) then
        wt=0._dp
        return
    endif
    ybar=(psumjet(4)+psumjet(3))/(psumjet(4)-psumjet(3))
!--- check that argument of upcoming log is not negative or infinite
    if (ybar <= 0._dp) then
        wt=0._dp
        return
    endif
    ybar=0.5_dp*log(ybar)

    ptsumjet2=psumjet(1)**2+psumjet(2)**2
    plstarsq=((sqrts**2-mv2-mjets**2)**2 &
    -4._dp*(mjets**2*mv2+ptsumjet2*sqrts**2))/(4._dp*sqrts**2)
!--- check that argument of upcoming sqrt is not negative
    if (plstarsq < 0._dp) then
        wt=0._dp
        return
    endif
    plstar=sqrt(plstarsq)
    Estar=plstarsq+ptsumjet2+mjets**2
!--- check that argument of upcoming sqrt is not negative
    if (Estar < 0._dp) then
        wt=0._dp
        return
    endif
    Estar=sqrt(Estar)
    if (Estar-plstar == 0._dp) then
        wt=0._dp
        return
    endif
    y5starmax=(Estar+plstar)/(Estar-plstar)
!--- check that argument of upcoming log is not negative or infinite
    if (y5starmax <= 0._dp) then
        wt=0._dp
        return
    endif
    y5starmax=0.5_dp*log(y5starmax)
    y5starmin=-y5starmax

    etamax=ybar-y5starmin
    etamin=ybar-y5starmax
    dely=etamax-etamin
    ycm=etamin+r(3*njets+1)*dely
    sinhy=sinh(ycm)
    coshy=sqrt(1._dp+sinhy**2)
          
!--- now make the initial state momenta
    sumpst=ptsumjet2+(psumjet(3)*coshy-psumjet(4)*sinhy)**2
    q0st=mv2+sumpst
!--- check that argument of upcoming sqrt is not negative
    if (q0st < 0._dp) then
        wt=0._dp
        return
    endif
    q0st=sqrt(q0st)
    rshat=mjets**2+sumpst
!--- check that argument of upcoming sqrt is not negative
    if (rshat < 0._dp) then
        wt=0._dp
        return
    endif
    rshat=q0st+sqrt(rshat)
    pcm(4)=rshat*coshy
    pcm(3)=rshat*sinhy
                
    xx(1)=(pcm(4)+pcm(3))/sqrts
    xx(2)=(pcm(4)-pcm(3))/sqrts
          
    if   ((xx(1) > 1._dp) .OR. (xx(2) > 1._dp)) then
        wt=0._dp
        return
    endif
          
    wt=wt*dely
    do j=1,4
        Q(j)=pcm(j)-psumjet(j)
    enddo
          
    p(1,4)=-xx(1)*sqrts/2._dp
    p(1,3)=p(1,4)
    p(2,4)=-xx(2)*sqrts/2._dp
    p(2,3)=-p(2,4)
          
    wt=wt*rshat/(sqrts**2*q0st)

    ! decay of the top-quark piece here

        call phi1_2m(mb,r(3*njets+2),r(3*njets+3),r(3*njets+4),zip, &
        Q,bm,wp,wtwp,*999)
        call phi3m0(r(3*njets+5),r(3*njets+6),wp,nn,ep,wtepnn,*999)
        wt=wt0*wt*wtwp*wtepnn

    do nu=1,4
        p(3,nu)=nn(nu)
        p(4,nu)=ep(nu)
        p(5,nu)=bm(nu)
    enddo

    return
          
999 wt=0._dp
    return
          
    end subroutine gen_stop_new

    subroutine gen_stop_decay(r,njets,ptgen,p,wt)
    implicit none
    include 'types.f'
          
    include 'constants.f'
    include 'nf.f'
    include 'mxpart.f'
    include 'cplx.h'
    include 'masses.f'
    include 'mxdim.f'
    include 'limits.f'
    include 'zerowidth.f'
    include 'kprocess.f'
    include 'reset.f'
    include 'kpart.f'
    include 'x1x2.f'
    include 'energy.f'
    include 'notag.f'
    include 'taucut.f'
    include 'tiny.f'
    include 'debug.f'
!---- Generate phase space for 2-->2+n process
!---- with (345) being a top and 6,..,5+n the jets
!---- r(mxdim),p1(4),p2(4) are inputs reversed in sign
!---- from physical values
!---- phase space for -p1-p2 --> p3+p4+p5+p6
!---- with all 2 pi's (ie 1/(2*pi)^(3n-4), where n is the number
!----  of final state particles)
!---- This routine has a minimum of 4 final state particles, hence
!---- the twopi**2 correction factor is given by the ratio of
!---- (1/twopi)**(3n-4) present in the phase space and the factor
!---- of [(1/twopi)**2]**(n-1) from the number of branchings
!---- For the specific case 'ttdkay' where one of the jets is
!---- associated with the top quark decay, we must add an extra
!---- factor of (1/twopi) since the number of jets generated is
!---- larger than the value of 'njets' passed
    real(dp), intent(in) :: r(mxdim)
    integer, intent(in) :: njets
    integer, intent(in) :: ptgen
    real(dp), intent(out) :: p(mxpart,4)
    real(dp), intent(out) :: wt

    real(dp):: psumjet(4),pcm(4),Q(4)
    real(dp):: wt0,wtbg
    real(dp):: pt,etamax,etamin
    real(dp):: y,sinhy,coshy,phi,mv2,wtbw,mjets
    real(dp):: ybar,ptsumjet2,ycm,sumpst,q0st,rshat,dely,xjac
    real(dp):: plstar,estar,plstarsq,y5starmax,y5starmin
    real(dp):: bm(4),wp(4),nn(4),ep(4),pbg(4),g(4),wtwp,wtepnn
    integer:: j,nu,ijet,in
    parameter(wt0=1._dp/twopi**2)

    integer, parameter :: ndecayjets = 1

          
    p(1:8,:) = 0._dp
    pcm = 0._dp
    psumjet = 0._dp
    wt=2._dp*pi
                
    do ijet=1,njets
    !--- generate the pt of jet number ijet
    !--- rapidity limited by E=pT*coshy
        wt=wt/16._dp/pi**3

    !       ! some other ways to generate in pT
    !       ! directly in pT
    !       pt = sqrts*r(ijet)/2._dp
    !       xjac = r(ijet)*sqrts**2/4._dp
    !       ! 1/pt
    !       pt = sqrts*r(ijet)/(sqrts+2*r(ijet)-sqrts*r(ijet))
    !       xjac = r(ijet)*sqrts**3 / (2*r(ijet) - sqrts*(r(ijet)-1))**3
    !       ! -log(pT)
    !       pt = exp(log(sqrts/2._dp) - (1-r(ijet))/r(ijet))
    !       xjac = exp(2._dp*log(sqrts/2._dp) - 2._dp*(1-r(ijet))/r(ijet))/r(ijet)**2
        if (ptgen == ptgen_sqrt) then
            pt = r(ijet)**2*sqrts/2._dp
            xjac = pt*sqrts/2._dp*2._dp*r(ijet)
        elseif (ptgen == ptgen_linear) then
            pt = sqrts*r(ijet)/2._dp
            xjac = r(ijet)*sqrts**2/4._dp
        elseif (ptgen == ptgen_sqrt_invsq) then
            if (ijet == 1) then
                pt = r(ijet)**2*sqrts/2._dp
                xjac = pt*sqrts/2._dp*2._dp*r(ijet)
            else
                pt = 1._dp / sqrt(4._dp/sqrts**2 + (1-r(ijet))/r(ijet))
                xjac = sqrts**4 / 2._dp / (sqrts**2*(r(ijet)-1._dp) - 4*r(ijet))**2
            endif
        else
            pt = r(ijet)**2*sqrts/2._dp
            xjac = pt*sqrts/2._dp*2._dp*r(ijet)
        endif

    ! for nnlo double real currentContrib=1
!       if ((kcase == kbq_tpq_jet) .AND. (kpart == kreal) .AND. usescet .AND. (currentContrib == 1)) then
!           if (ijet == 1) then
!           ! in 1/sqrt(pT)
!               pt = r(ijet)**2*sqrts/2._dp
!               xjac = pt*sqrts/2._dp*2._dp*r(ijet)
!           else
!           ! 1/pt^2
!               pt = 1._dp / sqrt(4._dp/sqrts**2 + (1-r(ijet))/r(ijet))
!               xjac = sqrts**4 / 2._dp / (sqrts**2*(r(ijet)-1._dp) - 4*r(ijet))**2
!           endif
!       ! for nlo real, currentContrib == 2
!       elseif ((kcase == kbq_tpq) .AND. (kpart == kreal) .AND. ( .NOT. usescet) .AND. (currentContrib == 2)) then
!           error stop "todo"
!       else
!           error stop "please explicitly setup gen_stop for this case"
!       ! all genpt(r(ijet),ptbreak,.true.,pt,xjac)
!       endif

        wt=wt*xjac

        etamax=sqrts/2._dp/pt
        if (etamax**2 <= 1._dp) then
                if (debug) then
                    write(6,*) 'etamax**2 <= 1._dp in gen_stop.f',etamax**2
                endif
            wt=0._dp
            return
        endif

        etamax=log(etamax+sqrt(etamax**2-1._dp))
        etamax=min(etamax,100d0)
        y=etamax*(2._dp*r(njets+ijet)-1._dp)
        wt=wt*2._dp*etamax
                
        sinhy=sinh(y)
        coshy=sqrt(1._dp+sinhy**2)
                
        p(5+ijet,4)=pt*coshy
                
        phi=2._dp*pi*r(2*njets+ijet)
        wt=wt*2._dp*pi
                
        p(5+ijet,1)=pt*cos(phi)
        p(5+ijet,2)=pt*sin(phi)
        p(5+ijet,3)=pt*sinhy
                
        do nu=1,4
            psumjet(nu)=psumjet(nu)+p(5+ijet,nu)
        enddo
    enddo
          
    ! top generated zero width
      mv2 = mt**2
      wtbw = pi*mt*twidth

    wt=wt*wtbw
!--- for one jet, mjets must be exactly zero
    if (njets == 1) then
        mjets=0._dp
    else
    !--- invariant mass of jets
        mjets=psumjet(4)**2-psumjet(1)**2-psumjet(2)**2-psumjet(3)**2
    !--- check that argument of upcoming sqrt is not negative
        if (mjets < 0._dp) then
            wt=0._dp
            return
        endif
        mjets=sqrt(mjets)
    endif
          
    if (psumjet(4)-psumjet(3) == 0._dp) then
        wt=0._dp
        return
    endif
    ybar=(psumjet(4)+psumjet(3))/(psumjet(4)-psumjet(3))
!--- check that argument of upcoming log is not negative or infinite
    if (ybar <= 0._dp) then
        wt=0._dp
        return
    endif
    ybar=0.5_dp*log(ybar)

    ptsumjet2=psumjet(1)**2+psumjet(2)**2
    plstarsq=((sqrts**2-mv2-mjets**2)**2 &
    -4._dp*(mjets**2*mv2+ptsumjet2*sqrts**2))/(4._dp*sqrts**2)
!--- check that argument of upcoming sqrt is not negative
    if (plstarsq < 0._dp) then
        wt=0._dp
        return
    endif
    plstar=sqrt(plstarsq)
    Estar=plstarsq+ptsumjet2+mjets**2
!--- check that argument of upcoming sqrt is not negative
    if (Estar < 0._dp) then
        wt=0._dp
        return
    endif
    Estar=sqrt(Estar)
    if (Estar-plstar == 0._dp) then
        wt=0._dp
        return
    endif
    y5starmax=(Estar+plstar)/(Estar-plstar)
!--- check that argument of upcoming log is not negative or infinite
    if (y5starmax <= 0._dp) then
        wt=0._dp
        return
    endif
    y5starmax=0.5_dp*log(y5starmax)
    y5starmin=-y5starmax

    etamax=ybar-y5starmin
    etamin=ybar-y5starmax
    dely=etamax-etamin
    ycm=etamin+r(3*njets+1)*dely
    sinhy=sinh(ycm)
    coshy=sqrt(1._dp+sinhy**2)
          
!--- now make the initial state momenta
    sumpst=ptsumjet2+(psumjet(3)*coshy-psumjet(4)*sinhy)**2
    q0st=mv2+sumpst
!--- check that argument of upcoming sqrt is not negative
    if (q0st < 0._dp) then
        wt=0._dp
        return
    endif
    q0st=sqrt(q0st)
    rshat=mjets**2+sumpst
!--- check that argument of upcoming sqrt is not negative
    if (rshat < 0._dp) then
        wt=0._dp
        return
    endif
    rshat=q0st+sqrt(rshat)
    pcm(4)=rshat*coshy
    pcm(3)=rshat*sinhy
                
    xx(1)=(pcm(4)+pcm(3))/sqrts
    xx(2)=(pcm(4)-pcm(3))/sqrts
          
    if   ((xx(1) > 1._dp) .OR. (xx(2) > 1._dp)) then
        wt=0._dp
        return
    endif
          
    wt=wt*dely
    do j=1,4
        Q(j)=pcm(j)-psumjet(j)
    enddo
          
    p(1,4)=-xx(1)*sqrts/2._dp
    p(1,3)=p(1,4)
    p(2,4)=-xx(2)*sqrts/2._dp
    p(2,3)=-p(2,4)
          
    wt=wt*rshat/(sqrts**2*q0st)

    ! decay of the top-quark piece here
    ! Q = top

    in=3*njets+2
    call phi1_2(r(in),r(in+1),r(in+2),r(in+3),Q,pbg,wp,wtwp,*999)
    in=in+4
      call phi3m0(r(in),r(in+1),pbg,bm,g,wtbg,*999)
      call phi3m0(r(in+2),r(in+3),wp,nn,ep,wtepnn,*999)
        wt=wt0*wt*wtwp*wtbg*wtepnn/twopi

    do nu=1,4
      p(5+njets+1,nu)=g(nu)
    enddo

    do nu=1,4
        p(3,nu)=nn(nu)
        p(4,nu)=ep(nu)
        p(5,nu)=bm(nu)
    enddo

    return
          
999 wt=0._dp
    return
          
    end subroutine gen_stop_decay

end module
