
!--- Common block that identifies the heavy flavour being used
!--- in the current process
      integer:: flav
      common/heavyflav/flav
!$omp threadprivate(/heavyflav/)
