      real(dp) :: delg1_z,delg1_g,lambda_g,lambda_z
      real(dp) :: h1Z,h2Z,h3Z,h4Z,h1gam,h2gam,h3gam,h4gam
      real(dp) :: delk_g,delk_z,tevscale
!      real(dp) :: h1tZ,h2tZ,h3tZ,h4tZ,h1tgam,h2tgam,h3tgam,h4tgam
      logical :: anomtgc
      real(dp) :: hitZ(4), hitgam(4)
      common/anomcoup1/delg1_z,delg1_g,lambda_g,lambda_z,delk_g,delk_z
      common/anomcoup2/h1Z,h2Z,h3Z,h4Z,h1gam,h2gam,h3gam,h4gam,tevscale,hitZ,hitgam,anomtgc
!$omp threadprivate(/anomcoup1/,/anomcoup2/)

