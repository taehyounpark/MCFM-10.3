
      include 'maxd.f'
      integer:: ndmax
!----ndmax=The maximum number of dipoles for the problem at hand
      real(dp):: ptilde(0:maxd,mxpart,4)
      real(dp):: ptildejet(mxpart,4,0:maxd)
      common/ptildes/ptilde,ptildejet,ndmax
!$omp threadprivate(/ptildes/)
