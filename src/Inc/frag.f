
!--- include file for controlling fragmentation/Frixione variables

!--- the following are set by the input file
      logical:: frag ! whether fragmentation process should be included
      character(8):: fragset         ! label for fragmentation set
      real(dp):: frag_scale ! fragmentation scale
      real(dp):: cone_ang   ! cone size for Frixione isolation
      real(dp):: epsilon_h  ! energy fraction for isolation
      real(dp):: frag_scalestart ! frag. scale value in input file
      real(dp):: n_pow   ! exponent for smooth-cone (Frixione) isolation

!-- the following is used when computing fragmentation processes
      real(dp):: z_frag ! energy fraction carried by photon
      logical:: rescale ! Indicates if p_part->1/z*p_gamma or not
!      real(dp):: p_phys(mxpart,4) ! Physical momenta with jets
!                                        ! rescaled by factor of z to photons

      common/fraginputs/frag_scale,cone_ang,epsilon_h,frag_scalestart,frag,fragset
      common/n_pow_common/n_pow
!$omp threadprivate(/n_pow_common/)
      common/fragvars/z_frag,rescale
!===== logical:: variable to specific fragintmore
      logical:: fragint_mode
      common/fragint_mode/fragint_mode
!$omp threadprivate(/fragvars/)

      logical :: fixed_coneenergy
      common/fixed_coneenergy/fixed_coneenergy

      logical :: photiso_hybrid
      real(dp) :: R_inner
      common/photiso_hybrid/R_inner,photiso_hybrid
