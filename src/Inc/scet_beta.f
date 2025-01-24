      include 'nflav.f'

      real(dp):: be0,be1,Ga1,Gaq1,Gag1,gBq1,gBg0,gBg1,gams1

      be0 = 11/three*CA-4/three*TR*nflav
      be1= 34/three*CA**2-(20/three*CA+4*CF)*TR*nflav

      Ga1 = 4/three*((four-pisq)*CA+5*be0)
      Gaq1 = Ga1*CF
      gBq1 = CF*(CA*(146/nine-80*zeta3)+CF*(three-4*pisq+48*zeta3)+be0*(121/nine+2*pisq/three))

      Gag1 = Ga1*CA
      gBg0 = 2*be0
      gBg1 = CA*(CA*(182/nine-32*zeta3)+be0*(94/nine-2/three*pisq))+2*be1
      gams1=CA*(-64/nine+28*zeta3)+be0*(-56/nine+2*zeta2)
