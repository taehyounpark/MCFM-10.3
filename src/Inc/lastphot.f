
!--- integer:: index of last particle that is a photon,
!--- used in dipolesfrag to shuffle fragmented photon to last position
      integer:: lastphot
      common/lastphot/lastphot
!$omp threadprivate(/lastphot/)
