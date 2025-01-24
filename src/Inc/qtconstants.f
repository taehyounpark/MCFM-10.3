      real(dp)::beta0,beta1,beta2,Gamma0(0:1),Gamma1(0:1),Gamma2(0:1),
     & tgamman0(0:1),tgamman1(0:1),tgamman2(0:1),
     & gammaB0(0:1),gammaB1(0:1),gammaB2(0:1),
     & gammaS0(0:1),gammaS1(0:1),gammaS2(0:1),
     & tgammaB0(0:1),tgammaB1(0:1),tgammaB2(0:1),
     & tgammaS0(0:1),tgammaS1(0:1),tgammaS2(0:1),
     & tildes0(0:1),tildes1(0:1),tildes2(0:1),tildeS3(0:1),S2(0:1)
      common/qtconstants/beta0,beta1,beta2,Gamma0,Gamma1,Gamma2,
     & tgamman0,tgamman1,tgamman2,
     & gammaB0,gammaB1,gammaB2,
     & gammaS0,gammaS1,gammaS2,
     & tgammaB0,tgammaB1,tgammaB2,
     & tgammaS0,tgammaS1,tgammaS2,
     & tildes0,tildes1,tildes2,tildes3,S2
!$omp threadprivate(/qtconstants/)
