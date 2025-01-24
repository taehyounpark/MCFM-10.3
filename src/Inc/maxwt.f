
! --- Common block for keeping track of weights, used
! --- if unweighting is selected :
      real(dp):: wtmax,newwt
      logical:: evtgen
      logical, parameter :: unweight = .false.
      logical:: skipnt
      integer:: nevtrequested
      common/maxwt/wtmax,newwt,nevtrequested,evtgen,skipnt

