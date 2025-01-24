      complex(dp):: coupfac(8,2,2,2,2) ! args are (diagtype,qtype,h12,h56,h78)
      common/coupfac_vv_common/coupfac
!$omp threadprivate(/coupfac_vv_common/)
