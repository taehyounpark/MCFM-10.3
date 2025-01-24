
      integer, parameter :: klord=1, kvirt=2, kreal=3, ktota=4, kfrag=5, ktodk=6
      integer, parameter :: ksnlo=7, knnlo=8
      integer, parameter :: kresummed=9
      integer, parameter :: kn3lo = 10
      integer, parameter :: ksnloR=1, ksnloV=2
      integer, parameter :: knnloVV=1, knnloRV=2, knnloRR=3
      integer kpart,ksnlopart,knnlopart
      integer:: origKpart
      logical:: coeffonly
      common/kpart/kpart
      common/origKpart/origKpart
      common/ksnlopart/ksnlopart
      common/knnlopart/knnlopart
      common/coeffonly/coeffonly

      integer, parameter :: kresexp = 1, kresonly = 2, kresabove = 3
      integer, parameter :: kresmatchcorr = 4
      integer :: krespart
      common/krespart/krespart
      integer :: kresorder
      common/kresorder/kresorder

