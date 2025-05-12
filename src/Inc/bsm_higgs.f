      real(dp) :: t1, t2, t3, t4, t5, t6
      real(dp) :: w1, w2, w3, w4, w5
      real(dp) :: cx, cx_sm, mx
      real(dp) :: c6_cfg(3), c6_init, c6_step, c6, c6_sm
      real(dp) :: ct_cfg(3), ct_init, ct_step, ct, ct_sm
      real(dp) :: cg_cfg(3), cg_init, cg_step, cg, cg_sm
      integer  :: c6_nval, ct_nval, cg_nval
      character(len=255) :: bsm_higgs_scenario

      common /bsm_higgs/ t1, t2, t3, t4, t5, t6
      common /bsm_higgs/ w1, w2, w3, w4, w5
      common /bsm_higgs/ cx, cx_sm, mx
      common /bsm_higgs/ c6_cfg, c6_init, c6_step, c6, c6_sm
      common /bsm_higgs/ ct_cfg, ct_init, ct_step, ct, ct_sm
      common /bsm_higgs/ cg_cfg, cg_init, cg_step, cg, cg_sm
      common /bsm_higgs/ c6_nval, ct_nval, cg_nval
      common /bsm_higgs/ bsm_higgs_scenario
