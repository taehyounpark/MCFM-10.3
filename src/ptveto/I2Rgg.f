!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      I2R(gg,delt)= + R10 * (  - 149.D0/3870720.D0*CA**2 )
      I2R(gg,delt) = I2R(gg,delt) + R8 * (  - 30571386962947.D0/
     &    327168056033280000.D0*CA*nfl + 2692634479.D0/35407798272000.D0
     &    *CA*ln2*nfl - 4001.D0/59875200.D0*CA*ln2**2*nfl - 
     &    269462059416217.D0/163584028016640000.D0*CA**2 + 1.D0/518400.D
     &    0*CA**2*zeta3 + 48471314251.D0/70815596544000.D0*CA**2*ln2 - 
     &    28529.D0/119750400.D0*CA**2*ln2**2 )
      I2R(gg,delt) = I2R(gg,delt) + R6 * ( 651413779.D0/3840721920000.D0
     &    *CA*nfl - 1272767.D0/3048192000.D0*CA*ln2*nfl + 11.D0/12600.D0
     &    *CA*ln2**2*nfl - 193271596459.D0/3840721920000.D0*CA**2 - 
     &    25763593.D0/3048192000.D0*CA**2*ln2 + 23.D0/4200.D0*CA**2*
     &    ln2**2 )
      I2R(gg,delt) = I2R(gg,delt) + R4 * ( 13532213.D0/5334336000.D0*CA
     &    *nfl - 57311.D0/6350400.D0*CA*ln2*nfl - 53.D0/3780.D0*CA*
     &    ln2**2*nfl + 9559878361.D0/5334336000.D0*CA**2 + 1.D0/72.D0*
     &    CA**2*zeta3 - 23695409.D0/12700800.D0*CA**2*ln2 - 587.D0/7560.
     &    D0*CA**2*ln2**2 )
      I2R(gg,delt) = I2R(gg,delt) + R2 * (  - 13033.D0/81000.D0*CA*nfl
     &     + 881.D0/5400.D0*CA*ln2*nfl + 7.D0/45.D0*CA*ln2**2*nfl + 
     &    368533.D0/81000.D0*CA**2 + 3*CA**2*zeta3 - 23381.D0/5400.D0*
     &    CA**2*ln2 - 52.D0/45.D0*CA**2*ln2**2 )
      I2R(gg,delt) = I2R(gg,delt) + lnR * ( 326.D0/27.D0*CA*nfl - 116.D0
     &    /9.D0*CA*ln2*nfl - 16.D0/3.D0*CA*ln2**2*nfl - 1622.D0/27.D0*
     &    CA**2 + 8*CA**2*zeta3 + 548.D0/9.D0*CA**2*ln2 + 88.D0/3.D0*
     &    CA**2*ln2**2 )
      I2R(gg,delt) = I2R(gg,delt) - 325.D0/9.D0*CA*nfl + 29.D0/54.D0*CA
     &    *pisq*nfl + 640.D0/27.D0*CA*ln2*nfl + 4.D0/9.D0*CA*ln2*pisq*
     &    nfl + 20*CA*ln2**2*nfl + 8.D0/9.D0*CA*ln2**2*ln64*nfl + 16.D0/
     &    9.D0*CA*ln2**3*nfl + 1621.D0/9.D0*CA**2 + 128*CA**2*Li4xhalf
     &     - 137.D0/54.D0*CA**2*pisq - 47.D0/30.D0*CA**2*pisq**2 + 32*
     &    CA**2*zeta3*lnmuonpt - 3232.D0/27.D0*CA**2*ln2 - 22.D0/9.D0*
     &    CA**2*ln2*pisq + 112*CA**2*ln2*zeta3 - 92*CA**2*ln2**2 - 16.D0
     &    /3.D0*CA**2*ln2**2*pisq - 352.D0/9.D0*CA**2*ln2**3 + 16.D0/3.D
     &    0*CA**2*ln2**4

      I2R(gg,plus)= + R8 * (  - 5664846191.D0/35407798272000.D0*CA*nfl
     &     + 4001.D0/29937600.D0*CA*ln2*nfl )
      I2R(gg,plus) = I2R(gg,plus) + R6 * ( 7001023.D0/3048192000.D0*CA*
     &    nfl - 11.D0/6300.D0*CA*ln2*nfl + 74801417.D0/3048192000.D0*
     &    CA**2 - 23.D0/2100.D0*CA**2*ln2 )
      I2R(gg,plus) = I2R(gg,plus) + R4 * (  - 168401.D0/6350400.D0*CA*
     &    nfl + 53.D0/1890.D0*CA*ln2*nfl + 16018321.D0/12700800.D0*
     &    CA**2 - 1.D0/108.D0*CA**2*pisq + 587.D0/3780.D0*CA**2*ln2 )
      I2R(gg,plus) = I2R(gg,plus) + R2 * ( 3071.D0/5400.D0*CA*nfl - 14.D
     &    0/45.D0*CA*ln2*nfl + 1429.D0/5400.D0*CA**2 - 2*CA**2*pisq + 
     &    104.D0/45.D0*CA**2*ln2 )
      I2R(gg,plus) = I2R(gg,plus) + picube * (  - 4*CA*nfl*BCGGCATF1 - 
     &    8*CA**2*BCGGCACA1 - 4*CF*nfl*BCGGCFTF1 )
      I2R(gg,plus) = I2R(gg,plus) + lnR * (  - 92.D0/9.D0*CA*nfl + 32.D0
     &    /3.D0*CA*ln2*nfl + 524.D0/9.D0*CA**2 - 16.D0/3.D0*CA**2*pisq
     &     - 176.D0/3.D0*CA**2*ln2 )
      I2R(gg,plus) = I2R(gg,plus) - 26.D0/9.D0*CA*nfl + 92.D0/9.D0*CA*
     &    lnpi*nfl + 116.D0/9.D0*CA*ln2*nfl - 32.D0/3.D0*CA*ln2*lnpi*
     &    nfl - 32.D0/3.D0*CA*ln2**2*nfl + 26.D0/9.D0*CA**2 - 524.D0/9.D
     &    0*CA**2*lnpi + 16.D0/3.D0*CA**2*pisq*lnpi - 548.D0/9.D0*CA**2
     &    *ln2 + 176.D0/3.D0*CA**2*ln2*lnpi + 16.D0/3.D0*CA**2*ln2*pisq
     &     + 176.D0/3.D0*CA**2*ln2**2

      if (x < xmin) then
      I2R(gg,rglr)= + R8 * ( 384581.D0/8360755200.D0*CA*x**(-1)*nfl + 
     &    2500865593.D0/2212987392000.D0*CA*nfl - 1010854877.D0/
     &    292626432000.D0*CA*omx*nfl + 49532206637.D0/8851949568000.D0*
     &    CA*omx**2*nfl - 1098589549.D0/147532492800.D0*CA*omx**3*nfl
     &     + 71623561.D0/9835499520.D0*CA*omx**4*nfl - 77879803.D0/
     &    24588748800.D0*CA*omx**5*nfl + 1.D0/388800.D0*CA*ln2*x**(-1)*
     &    nfl - 7099.D0/7484400.D0*CA*ln2*nfl + 1511.D0/544320.D0*CA*
     &    ln2*omx*nfl - 150719.D0/29937600.D0*CA*ln2*omx**2*nfl + 61.D0/
     &    9240.D0*CA*ln2*omx**3*nfl - 17.D0/2772.D0*CA*ln2*omx**4*nfl
     &     + 19.D0/6930.D0*CA*ln2*omx**5*nfl - 1.D0/6144.D0*CF*x**(-1)*
     &    nfl + 1.D0/6144.D0*CF*nfl - 41653.D0/19353600.D0*CF*omx*nfl
     &     + 7549723.D0/406425600.D0*CF*omx**2*nfl - 8339.D0/302400.D0*
     &    CF*omx**3*nfl + 4283.D0/1587600.D0*CF*omx**4*nfl + 91859.D0/
     &    10160640.D0*CF*omx**5*nfl + 1.D0/540.D0*CF*ln2*omx*nfl - 13.D0
     &    /945.D0*CF*ln2*omx**2*nfl + 1.D0/45.D0*CF*ln2*omx**3*nfl - 1.D
     &    0/315.D0*CF*ln2*omx**4*nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + R8 * (  - 1.D0/126.D0*CF*ln2*omx**5
     &    *nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + R6 * ( 5.D0/13824.D0*CA*x**(-1)*nfl
     &     - 36358423.D0/3048192000.D0*CA*nfl + 8015921.D0/304819200.D0
     &    *CA*omx*nfl - 578261.D0/21772800.D0*CA*omx**2*nfl + 1838579.D0
     &    /76204800.D0*CA*omx**3*nfl - 249667.D0/19051200.D0*CA*omx**4*
     &    nfl + 61.D0/6300.D0*CA*ln2*nfl - 41.D0/1890.D0*CA*ln2*omx*nfl
     &     + 7.D0/270.D0*CA*ln2*omx**2*nfl - 41.D0/1890.D0*CA*ln2*
     &    omx**3*nfl + 11.D0/945.D0*CA*ln2*omx**4*nfl - 160.D0/3.D0*
     &    CA**2*omx**(-7)*Li2mx + 296.D0/9.D0*CA**2*omx**(-7)*lnx + 160.
     &    D0/3.D0*CA**2*omx**(-7)*Li2omx - 184.D0/9.D0*CA**2*omx**(-6)
     &     + 80.D0/3.D0*CA**2*omx**(-6)*lnopxon2 + 880.D0/3.D0*CA**2*
     &    omx**(-6)*Li2mx - 1628.D0/9.D0*CA**2*omx**(-6)*lnx - 880.D0/3.
     &    D0*CA**2*omx**(-6)*Li2omx + 1040.D0/9.D0*CA**2*omx**(-5) - 
     &    376.D0/3.D0*CA**2*omx**(-5)*lnopxon2 - 2140.D0/3.D0*CA**2*
     &    omx**(-5)*Li2mx + 3977.D0/9.D0*CA**2*omx**(-5)*lnx + 2140.D0/
     &    3.D0*CA**2*omx**(-5)*Li2omx - 7577.D0/27.D0*CA**2*omx**(-4)
     &     + 2306.D0/9.D0*CA**2*omx**(-4)*lnopxon2 )
      I2R(gg,rglr) = I2R(gg,rglr) + R6 * ( 1008*CA**2*omx**(-4)*Li2mx
     &     - 1892.D0/3.D0*CA**2*omx**(-4)*lnx - 1008*CA**2*omx**(-4)*
     &    Li2omx + 20767.D0/54.D0*CA**2*omx**(-3) - 2693.D0/9.D0*CA**2*
     &    omx**(-3)*lnopxon2 - 2704.D0/3.D0*CA**2*omx**(-3)*Li2mx + 
     &    1736.D0/3.D0*CA**2*omx**(-3)*lnx + 2704.D0/3.D0*CA**2*
     &    omx**(-3)*Li2omx - 17161.D0/54.D0*CA**2*omx**(-2) + 1928.D0/9.
     &    D0*CA**2*omx**(-2)*lnopxon2 + 1528.D0/3.D0*CA**2*omx**(-2)*
     &    Li2mx - 1070.D0/3.D0*CA**2*omx**(-2)*lnx - 1528.D0/3.D0*CA**2
     &    *omx**(-2)*Li2omx + 853.D0/6.D0*CA**2*omx**(-1) - 847.D0/9.D0
     &    *CA**2*omx**(-1)*lnopxon2 - 164*CA**2*omx**(-1)*Li2mx + 1240.D
     &    0/9.D0*CA**2*omx**(-1)*lnx + 164*CA**2*omx**(-1)*Li2omx + 83.D
     &    0/4608.D0*CA**2*x**(-1) - 61.D0/48.D0*CA**2*opx**(-4) + 37.D0/
     &    48.D0*CA**2*opx**(-3) - 223.D0/192.D0*CA**2*opx**(-2) - 7.D0/
     &    6.D0*CA**2*opx**(-1) - 83980501697.D0/3048192000.D0*CA**2 + 
     &    100.D0/9.D0*CA**2*lnopxon2 + 88.D0/3.D0*CA**2*Li2mx - 64.D0/9.
     &    D0*CA**2*lnx )
      I2R(gg,rglr) = I2R(gg,rglr) + R6 * (  - 88.D0/3.D0*CA**2*Li2omx
     &     + 6688831123.D0/304819200.D0*CA**2*omx + 155.D0/9.D0*CA**2*
     &    omx*lnopxon2 - 44.D0/3.D0*CA**2*omx*Li2mx - 203.D0/9.D0*CA**2
     &    *omx*lnx + 44.D0/3.D0*CA**2*omx*Li2omx - 2542644463.D0/
     &    152409600.D0*CA**2*omx**2 - 20.D0/3.D0*CA**2*omx**2*lnopxon2
     &     + 20.D0/3.D0*CA**2*omx**2*Li2mx + 20.D0/3.D0*CA**2*omx**2*
     &    lnx - 20.D0/3.D0*CA**2*omx**2*Li2omx + 143265127.D0/76204800.D
     &    0*CA**2*omx**3 + 539159.D0/3810240.D0*CA**2*omx**4 - 40.D0/9.D
     &    0*CA**2*pisq*omx**(-7) + 220.D0/9.D0*CA**2*pisq*omx**(-6) - 
     &    535.D0/9.D0*CA**2*pisq*omx**(-5) + 84*CA**2*pisq*omx**(-4) - 
     &    676.D0/9.D0*CA**2*pisq*omx**(-3) + 382.D0/9.D0*CA**2*pisq*
     &    omx**(-2) - 41.D0/3.D0*CA**2*pisq*omx**(-1) + 22.D0/9.D0*
     &    CA**2*pisq - 11.D0/9.D0*CA**2*pisq*omx + 5.D0/9.D0*CA**2*pisq
     &    *omx**2 + 160.D0/3.D0*CA**2*ln2*omx**(-6) - 800.D0/3.D0*CA**2
     &    *ln2*omx**(-5) + 5260.D0/9.D0*CA**2*ln2*omx**(-4) - 6622.D0/9.
     &    D0*CA**2*ln2*omx**(-3) )
      I2R(gg,rglr) = I2R(gg,rglr) + R6 * ( 5152.D0/9.D0*CA**2*ln2*
     &    omx**(-2) - 2395.D0/9.D0*CA**2*ln2*omx**(-1) + 372229.D0/6300.
     &    D0*CA**2*ln2 - 8804.D0/945.D0*CA**2*ln2*omx + 9943.D0/945.D0*
     &    CA**2*ln2*omx**2 - 2024.D0/945.D0*CA**2*ln2*omx**3 - 25.D0/
     &    189.D0*CA**2*ln2*omx**4 + 473.D0/26460.D0*CF*omx*nfl - 18163.D
     &    0/151200.D0*CF*omx**2*nfl + 7631.D0/105840.D0*CF*omx**3*nfl
     &     + 14863.D0/529200.D0*CF*omx**4*nfl - 1.D0/63.D0*CF*ln2*omx*
     &    nfl + 4.D0/45.D0*CF*ln2*omx**2*nfl - 4.D0/63.D0*CF*ln2*omx**3
     &    *nfl - 8.D0/315.D0*CF*ln2*omx**4*nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + R4 * (  - 179.D0/3240.D0*CA*x**(-1)
     &    *nfl + 935581.D0/6350400.D0*CA*nfl - 1182347.D0/6350400.D0*CA
     &    *omx*nfl + 1737949.D0/12700800.D0*CA*omx**2*nfl - 113.D0/2205.
     &    D0*CA*omx**3*nfl + 1.D0/54.D0*CA*ln2*x**(-1)*nfl - 223.D0/
     &    1890.D0*CA*ln2*nfl + 311.D0/1890.D0*CA*ln2*omx*nfl - 103.D0/
     &    945.D0*CA*ln2*omx**2*nfl + 1.D0/21.D0*CA*ln2*omx**3*nfl - 24*
     &    CA**2*omx**(-5)*Li2mx + 14*CA**2*omx**(-5)*lnx + 24*CA**2*
     &    omx**(-5)*Li2omx - 10*CA**2*omx**(-4) + 12*CA**2*omx**(-4)*
     &    lnopxon2 + 96*CA**2*omx**(-4)*Li2mx - 56*CA**2*omx**(-4)*lnx
     &     - 96*CA**2*omx**(-4)*Li2omx + 41*CA**2*omx**(-3) - 38*CA**2*
     &    omx**(-3)*lnopxon2 - 482.D0/3.D0*CA**2*omx**(-3)*Li2mx + 289.D
     &    0/3.D0*CA**2*omx**(-3)*lnx + 482.D0/3.D0*CA**2*omx**(-3)*
     &    Li2omx - 389.D0/6.D0*CA**2*omx**(-2) + 145.D0/3.D0*CA**2*
     &    omx**(-2)*lnopxon2 + 436.D0/3.D0*CA**2*omx**(-2)*Li2mx - 93*
     &    CA**2*omx**(-2)*lnx - 436.D0/3.D0*CA**2*omx**(-2)*Li2omx + 
     &    601.D0/12.D0*CA**2*omx**(-1) )
      I2R(gg,rglr) = I2R(gg,rglr) + R4 * (  - 104.D0/3.D0*CA**2*
     &    omx**(-1)*lnopxon2 - 58*CA**2*omx**(-1)*Li2mx + 125.D0/3.D0*
     &    CA**2*omx**(-1)*lnx + 58*CA**2*omx**(-1)*Li2omx + 42959.D0/
     &    51840.D0*CA**2*x**(-1) - 49.D0/8.D0*CA**2*opx**(-2) - 95.D0/8.
     &    D0*CA**2*opx**(-1) - 77744551.D0/6350400.D0*CA**2 + 9*CA**2*
     &    lnopxon2 + 52.D0/3.D0*CA**2*Li2mx + 19.D0/3.D0*CA**2*lnx - 52.
     &    D0/3.D0*CA**2*Li2omx + 415571551.D0/12700800.D0*CA**2*omx + 4.
     &    D0/3.D0*CA**2*omx*lnopxon2 - 74.D0/3.D0*CA**2*omx*Li2mx - 
     &    CA**2*omx*lnx + 74.D0/3.D0*CA**2*omx*Li2omx - 67030661.D0/
     &    6350400.D0*CA**2*omx**2 + 25.D0/3.D0*CA**2*omx**2*lnopxon2 + 
     &    8*CA**2*omx**2*Li2mx - 25.D0/3.D0*CA**2*omx**2*lnx - 8*CA**2*
     &    omx**2*Li2omx - 117731.D0/117600.D0*CA**2*omx**3 - 2*CA**2*
     &    pisq*omx**(-5) + 8*CA**2*pisq*omx**(-4) - 241.D0/18.D0*CA**2*
     &    pisq*omx**(-3) + 109.D0/9.D0*CA**2*pisq*omx**(-2) - 29.D0/6.D0
     &    *CA**2*pisq*omx**(-1) - 1.D0/108.D0*CA**2*pisq*x**(-1) + 79.D0
     &    /54.D0*CA**2*pisq )
      I2R(gg,rglr) = I2R(gg,rglr) + R4 * (  - 223.D0/108.D0*CA**2*pisq*
     &    omx + 73.D0/108.D0*CA**2*pisq*omx**2 + 24*CA**2*ln2*omx**(-4)
     &     - 84*CA**2*ln2*omx**(-3) + 362.D0/3.D0*CA**2*ln2*omx**(-2)
     &     - 91*CA**2*ln2*omx**(-1) - 11.D0/108.D0*CA**2*ln2*x**(-1) + 
     &    40543.D0/1890.D0*CA**2*ln2 - 55243.D0/3780.D0*CA**2*ln2*omx
     &     + 69901.D0/3780.D0*CA**2*ln2*omx**2 + 38.D0/35.D0*CA**2*ln2*
     &    omx**3 + 1.D0/8.D0*CF*x**(-1)*nfl - 1.D0/8.D0*CF*nfl - 661.D0/
     &    3600.D0*CF*omx*nfl + 1939.D0/3600.D0*CF*omx**2*nfl + 2.D0/15.D
     &    0*CF*ln2*omx*nfl - 8.D0/15.D0*CF*ln2*omx**2*nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + R2 * ( 55.D0/36.D0*CA*x**(-1)*nfl
     &     - 479.D0/216.D0*CA*nfl + 2419.D0/1800.D0*CA*omx*nfl - 371.D0/
     &    225.D0*CA*omx**2*nfl - 4.D0/3.D0*CA*ln2*x**(-1)*nfl + 16.D0/9.
     &    D0*CA*ln2*nfl - 16.D0/15.D0*CA*ln2*omx*nfl + 22.D0/15.D0*CA*
     &    ln2*omx**2*nfl - 16*CA**2*omx**(-3)*Li2mx + 8*CA**2*omx**(-3)
     &    *lnx + 16*CA**2*omx**(-3)*Li2omx - 8*CA**2*omx**(-2) + 8*
     &    CA**2*omx**(-2)*lnopxon2 + 40*CA**2*omx**(-2)*Li2mx - 24*
     &    CA**2*omx**(-2)*lnx - 40*CA**2*omx**(-2)*Li2omx + 16*CA**2*
     &    omx**(-1) - 16*CA**2*omx**(-1)*lnopxon2 - 24*CA**2*omx**(-1)*
     &    Li2mx + 40*CA**2*omx**(-1)*lnx + 24*CA**2*omx**(-1)*Li2omx + 
     &    289.D0/12.D0*CA**2*x**(-1) - 8*CA**2*opx**(-1)*Li2mx + 8*
     &    CA**2*opx**(-1)*Li2omx + 419.D0/216.D0*CA**2 + 24*CA**2*
     &    lnopxon2 - 24*CA**2*Li2mx - 56*CA**2*lnx + 24*CA**2*Li2omx - 
     &    78119.D0/1800.D0*CA**2*omx - 32*CA**2*omx*lnopxon2 + 32*CA**2
     &    *omx*Li2mx + 56*CA**2*omx*lnx - 32*CA**2*omx*Li2omx + 989.D0/
     &    150.D0*CA**2*omx**2 )
      I2R(gg,rglr) = I2R(gg,rglr) + R2 * ( 24*CA**2*omx**2*lnopxon2 - 
     &    24*CA**2*omx**2*Li2mx - 24*CA**2*omx**2*lnx + 24*CA**2*omx**2
     &    *Li2omx - 4.D0/3.D0*CA**2*pisq*omx**(-3) + 10.D0/3.D0*CA**2*
     &    pisq*omx**(-2) - 2*CA**2*pisq*omx**(-1) - 2*CA**2*pisq*
     &    x**(-1) - 2.D0/3.D0*CA**2*pisq*opx**(-1) + 2*CA**2*pisq + 2.D0
     &    /3.D0*CA**2*pisq*omx + 16*CA**2*ln2*omx**(-2) - 32*CA**2*ln2*
     &    omx**(-1) - 22*CA**2*ln2*x**(-1) + 158.D0/9.D0*CA**2*ln2 + 
     &    206.D0/15.D0*CA**2*ln2*omx + 76.D0/5.D0*CA**2*ln2*omx**2 - 46.
     &    D0/9.D0*CF*x**(-1)*nfl + 46.D0/9.D0*CF*nfl + 23.D0/9.D0*CF*
     &    omx*nfl + 23.D0/9.D0*CF*omx**2*nfl + 16.D0/3.D0*CF*ln2*
     &    x**(-1)*nfl - 16.D0/3.D0*CF*ln2*nfl - 8.D0/3.D0*CF*ln2*omx*
     &    nfl - 8.D0/3.D0*CF*ln2*omx**2*nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + picube * (  - 4*CA*omx**(-1)*nfl*
     &    BCGGCATFx + 4*CA*omx**(-1)*nfl*BCGGCATF1 + 4*CA*nfl*BCGGCATFx
     &     - 8*CA**2*omx**(-1)*BCGGCACAx + 8*CA**2*omx**(-1)*BCGGCACA1
     &     + 8*CA**2*BCGGCACAx - 4*CF*omx**(-1)*nfl*BCGGCFTFx + 4*CF*
     &    omx**(-1)*nfl*BCGGCFTF1 + 4*CF*nfl*BCGGCFTFx )
      I2R(gg,rglr) = I2R(gg,rglr) + lnR * (  - 92.D0/9.D0*CA*x**(-1)*
     &    nfl + 184.D0/9.D0*CA*nfl - 92.D0/9.D0*CA*omx*nfl + 92.D0/9.D0
     &    *CA*omx**2*nfl + 32.D0/3.D0*CA*ln2*x**(-1)*nfl - 64.D0/3.D0*
     &    CA*ln2*nfl + 32.D0/3.D0*CA*ln2*omx*nfl - 32.D0/3.D0*CA*ln2*
     &    omx**2*nfl + 524.D0/9.D0*CA**2*x**(-1) - 1048.D0/9.D0*CA**2
     &     + 524.D0/9.D0*CA**2*omx - 524.D0/9.D0*CA**2*omx**2 - 16.D0/3.
     &    D0*CA**2*pisq*x**(-1) + 32.D0/3.D0*CA**2*pisq - 16.D0/3.D0*
     &    CA**2*pisq*omx + 16.D0/3.D0*CA**2*pisq*omx**2 - 176.D0/3.D0*
     &    CA**2*ln2*x**(-1) + 352.D0/3.D0*CA**2*ln2 - 176.D0/3.D0*CA**2
     &    *ln2*omx + 176.D0/3.D0*CA**2*ln2*omx**2 )
      I2R(gg,rglr) = I2R(gg,rglr) + 26.D0/9.D0*CA*nfl + 92.D0/9.D0*CA*
     &    lnpi*x**(-1)*nfl - 184.D0/9.D0*CA*lnpi*nfl + 92.D0/9.D0*CA*
     &    lnpi*omx*nfl - 92.D0/9.D0*CA*lnpi*omx**2*nfl + 92.D0/9.D0*CA*
     &    ln2*x**(-1)*nfl - 208.D0/9.D0*CA*ln2*nfl + 92.D0/9.D0*CA*ln2*
     &    omx*nfl - 92.D0/9.D0*CA*ln2*omx**2*nfl - 32.D0/3.D0*CA*ln2*
     &    lnpi*x**(-1)*nfl + 64.D0/3.D0*CA*ln2*lnpi*nfl - 32.D0/3.D0*CA
     &    *ln2*lnpi*omx*nfl + 32.D0/3.D0*CA*ln2*lnpi*omx**2*nfl - 32.D0/
     &    3.D0*CA*ln2**2*x**(-1)*nfl + 64.D0/3.D0*CA*ln2**2*nfl - 32.D0/
     &    3.D0*CA*ln2**2*omx*nfl + 32.D0/3.D0*CA*ln2**2*omx**2*nfl - 26.
     &    D0/9.D0*CA**2 - 524.D0/9.D0*CA**2*lnpi*x**(-1) + 1048.D0/9.D0
     &    *CA**2*lnpi - 524.D0/9.D0*CA**2*lnpi*omx + 524.D0/9.D0*CA**2*
     &    lnpi*omx**2 + 16.D0/3.D0*CA**2*pisq*lnpi*x**(-1) - 32.D0/3.D0
     &    *CA**2*pisq*lnpi + 16.D0/3.D0*CA**2*pisq*lnpi*omx - 16.D0/3.D0
     &    *CA**2*pisq*lnpi*omx**2 - 524.D0/9.D0*CA**2*ln2*x**(-1) + 
     &    1072.D0/9.D0*CA**2*ln2 - 524.D0/9.D0*CA**2*ln2*omx + 524.D0/9.
     &    D0*CA**2*ln2*omx**2
      I2R(gg,rglr) = I2R(gg,rglr) + 176.D0/3.D0*CA**2*ln2*lnpi*x**(-1)
     &     - 352.D0/3.D0*CA**2*ln2*lnpi + 176.D0/3.D0*CA**2*ln2*lnpi*
     &    omx - 176.D0/3.D0*CA**2*ln2*lnpi*omx**2 + 16.D0/3.D0*CA**2*
     &    ln2*pisq*x**(-1) - 32.D0/3.D0*CA**2*ln2*pisq + 16.D0/3.D0*
     &    CA**2*ln2*pisq*omx - 16.D0/3.D0*CA**2*ln2*pisq*omx**2 + 176.D0
     &    /3.D0*CA**2*ln2**2*x**(-1) - 352.D0/3.D0*CA**2*ln2**2 + 176.D0
     &    /3.D0*CA**2*ln2**2*omx - 176.D0/3.D0*CA**2*ln2**2*omx**2
      else
      I2R(gg,rglr)= + R8 * ( 41642550023.D0/35407798272000.D0*CA*nfl - 
     &    498697271.D0/146313216000.D0*CA*omx*nfl + 7398426929.D0/
     &    1311399936000.D0*CA*omx**2*nfl - 10481311649.D0/1416311930880.
     &    D0*CA*omx**3*nfl + 51894704027.D0/7081559654400.D0*CA*omx**4*
     &    nfl - 22103643157.D0/7081559654400.D0*CA*omx**5*nfl + 384581.D
     &    0/8360755200.D0*CA*omx**6*nfl + 384581.D0/8360755200.D0*CA*
     &    omx**7*nfl + 384581.D0/8360755200.D0*CA*omx**8*nfl + 384581.D0
     &    /8360755200.D0*CA*omx**9*nfl - 28319.D0/29937600.D0*CA*ln2*
     &    nfl + 3781.D0/1360800.D0*CA*ln2*omx*nfl - 8369.D0/1663200.D0*
     &    CA*ln2*omx**2*nfl + 197717.D0/29937600.D0*CA*ln2*omx**3*nfl
     &     - 183523.D0/29937600.D0*CA*ln2*omx**4*nfl + 82157.D0/
     &    29937600.D0*CA*ln2*omx**5*nfl + 1.D0/388800.D0*CA*ln2*omx**6*
     &    nfl + 1.D0/388800.D0*CA*ln2*omx**7*nfl + 1.D0/388800.D0*CA*
     &    ln2*omx**8*nfl + 1.D0/388800.D0*CA*ln2*omx**9*nfl - 44803.D0/
     &    19353600.D0*CF*omx*nfl + 7483573.D0/406425600.D0*CF*omx**2*
     &    nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + R8 * (  - 268423.D0/9676800.D0*CF*
     &    omx**3*nfl + 515149.D0/203212800.D0*CF*omx**4*nfl + 360821.D0/
     &    40642560.D0*CF*omx**5*nfl - 1.D0/6144.D0*CF*omx**6*nfl - 1.D0/
     &    6144.D0*CF*omx**7*nfl - 1.D0/6144.D0*CF*omx**8*nfl - 1.D0/
     &    6144.D0*CF*omx**9*nfl + 1.D0/540.D0*CF*ln2*omx*nfl - 13.D0/
     &    945.D0*CF*ln2*omx**2*nfl + 1.D0/45.D0*CF*ln2*omx**3*nfl - 1.D0
     &    /315.D0*CF*ln2*omx**4*nfl - 1.D0/126.D0*CF*ln2*omx**5*nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + R6 * (  - 35255923.D0/3048192000.D0
     &    *CA*nfl + 8126171.D0/304819200.D0*CA*omx*nfl - 285193.D0/
     &    10886400.D0*CA*omx**2*nfl + 3732283.D0/152409600.D0*CA*omx**3
     &    *nfl - 1942211.D0/152409600.D0*CA*omx**4*nfl + 5.D0/13824.D0*
     &    CA*omx**5*nfl + 5.D0/13824.D0*CA*omx**6*nfl + 5.D0/13824.D0*
     &    CA*omx**7*nfl + 5.D0/13824.D0*CA*omx**8*nfl + 5.D0/13824.D0*
     &    CA*omx**9*nfl + 61.D0/6300.D0*CA*ln2*nfl - 41.D0/1890.D0*CA*
     &    ln2*omx*nfl + 7.D0/270.D0*CA*ln2*omx**2*nfl - 41.D0/1890.D0*
     &    CA*ln2*omx**3*nfl + 11.D0/945.D0*CA*ln2*omx**4*nfl + 
     &    452322403.D0/3048192000.D0*CA**2 - 152091251.D0/304819200.D0*
     &    CA**2*omx + 62895233.D0/304819200.D0*CA**2*omx**2 + 1903711.D0
     &    /5080320.D0*CA**2*omx**3 - 3473878679.D0/18441561600.D0*CA**2
     &    *omx**4 - 216313487.D0/2235340800.D0*CA**2*omx**5 - 
     &    1915704559037.D0/99731965132800.D0*CA**2*omx**6 + 
     &    353098375787.D0/24932991283200.D0*CA**2*omx**7 + 
     &    2018612965331.D0/66487976755200.D0*CA**2*omx**8 )
      I2R(gg,rglr) = I2R(gg,rglr) + R6 * ( 1091332337089.D0/
     &    28494847180800.D0*CA**2*omx**9 + 10021157.D0/4354560.D0*CA**2
     &    *omx**10 - 1162447.D0/870912.D0*CA**2*omx**11 - 571.D0/6300.D0
     &    *CA**2*ln2 + 388.D0/945.D0*CA**2*ln2*omx - 229.D0/945.D0*
     &    CA**2*ln2*omx**2 - 37.D0/126.D0*CA**2*ln2*omx**3 + 2059.D0/
     &    10395.D0*CA**2*ln2*omx**4 + 23.D0/315.D0*CA**2*ln2*omx**5 + 
     &    752.D0/135135.D0*CA**2*ln2*omx**6 - 4021.D0/270270.D0*CA**2*
     &    ln2*omx**7 - 311.D0/15015.D0*CA**2*ln2*omx**8 - 829.D0/38610.D
     &    0*CA**2*ln2*omx**9 - 43.D0/54.D0*CA**2*ln2*omx**10 + 20.D0/27.
     &    D0*CA**2*ln2*omx**11 + 473.D0/26460.D0*CF*omx*nfl - 18163.D0/
     &    151200.D0*CF*omx**2*nfl + 7631.D0/105840.D0*CF*omx**3*nfl + 
     &    14863.D0/529200.D0*CF*omx**4*nfl - 1.D0/63.D0*CF*ln2*omx*nfl
     &     + 4.D0/45.D0*CF*ln2*omx**2*nfl - 4.D0/63.D0*CF*ln2*omx**3*
     &    nfl - 8.D0/315.D0*CF*ln2*omx**4*nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + R4 * ( 584741.D0/6350400.D0*CA*nfl
     &     - 1533187.D0/6350400.D0*CA*omx*nfl + 115141.D0/1411200.D0*CA
     &    *omx**2*nfl - 16907.D0/158760.D0*CA*omx**3*nfl - 179.D0/3240.D
     &    0*CA*omx**4*nfl - 179.D0/3240.D0*CA*omx**5*nfl - 179.D0/3240.D
     &    0*CA*omx**6*nfl - 179.D0/3240.D0*CA*omx**7*nfl - 179.D0/3240.D
     &    0*CA*omx**8*nfl - 179.D0/3240.D0*CA*omx**9*nfl - 94.D0/945.D0
     &    *CA*ln2*nfl + 173.D0/945.D0*CA*ln2*omx*nfl - 19.D0/210.D0*CA*
     &    ln2*omx**2*nfl + 25.D0/378.D0*CA*ln2*omx**3*nfl + 1.D0/54.D0*
     &    CA*ln2*omx**4*nfl + 1.D0/54.D0*CA*ln2*omx**5*nfl + 1.D0/54.D0
     &    *CA*ln2*omx**6*nfl + 1.D0/54.D0*CA*ln2*omx**7*nfl + 1.D0/54.D0
     &    *CA*ln2*omx**8*nfl + 1.D0/54.D0*CA*ln2*omx**9*nfl - 42933799.D
     &    0/12700800.D0*CA**2 + 28860487.D0/6350400.D0*CA**2*omx + 
     &    59173.D0/78400.D0*CA**2*omx**2 - 549883.D0/907200.D0*CA**2*
     &    omx**3 + 4890559.D0/3175200.D0*CA**2*omx**4 + 162470933.D0/
     &    101606400.D0*CA**2*omx**5 + 11161480831.D0/7376624640.D0*
     &    CA**2*omx**6 )
      I2R(gg,rglr) = I2R(gg,rglr) + R4 * ( 25996676701.D0/18441561600.D0
     &    *CA**2*omx**7 + 7289530947001.D0/5540664729600.D0*CA**2*
     &    omx**8 + 61654701651919.D0/49865982566400.D0*CA**2*omx**9 + 
     &    44143619.D0/17418240.D0*CA**2*omx**10 + 304201.D0/1451520.D0*
     &    CA**2*omx**11 + 1.D0/108.D0*CA**2*pisq - 1.D0/54.D0*CA**2*
     &    pisq*omx - 1.D0/108.D0*CA**2*pisq*omx**3 - 1.D0/108.D0*CA**2*
     &    pisq*omx**4 - 1.D0/108.D0*CA**2*pisq*omx**5 - 1.D0/108.D0*
     &    CA**2*pisq*omx**6 - 1.D0/108.D0*CA**2*pisq*omx**7 - 1.D0/108.D
     &    0*CA**2*pisq*omx**8 - 1.D0/108.D0*CA**2*pisq*omx**9 + 6277.D0/
     &    3780.D0*CA**2*ln2 - 5491.D0/1890.D0*CA**2*ln2*omx - 6.D0/35.D0
     &    *CA**2*ln2*omx**2 + 503.D0/540.D0*CA**2*ln2*omx**3 - 3733.D0/
     &    3780.D0*CA**2*ln2*omx**4 - 1619.D0/1890.D0*CA**2*ln2*omx**5
     &     - 5813.D0/8316.D0*CA**2*ln2*omx**6 - 12017.D0/20790.D0*CA**2
     &    *ln2*omx**7 - 264317.D0/540540.D0*CA**2*ln2*omx**8 - 114281.D0
     &    /270270.D0*CA**2*ln2*omx**9 - 47.D0/27.D0*CA**2*ln2*omx**10
     &     + 8.D0/9.D0*CA**2*ln2*omx**11 )
      I2R(gg,rglr) = I2R(gg,rglr) + R4 * (  - 211.D0/3600.D0*CF*omx*nfl
     &     + 2389.D0/3600.D0*CF*omx**2*nfl + 1.D0/8.D0*CF*omx**3*nfl + 
     &    1.D0/8.D0*CF*omx**4*nfl + 1.D0/8.D0*CF*omx**5*nfl + 1.D0/8.D0
     &    *CF*omx**6*nfl + 1.D0/8.D0*CF*omx**7*nfl + 1.D0/8.D0*CF*
     &    omx**8*nfl + 1.D0/8.D0*CF*omx**9*nfl + 2.D0/15.D0*CF*ln2*omx*
     &    nfl - 8.D0/15.D0*CF*ln2*omx**2*nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + R2 * (  - 149.D0/216.D0*CA*nfl + 
     &    1723.D0/600.D0*CA*omx*nfl - 109.D0/900.D0*CA*omx**2*nfl + 55.D
     &    0/36.D0*CA*omx**3*nfl + 55.D0/36.D0*CA*omx**4*nfl + 55.D0/36.D
     &    0*CA*omx**5*nfl + 55.D0/36.D0*CA*omx**6*nfl + 55.D0/36.D0*CA*
     &    omx**7*nfl + 55.D0/36.D0*CA*omx**8*nfl + 55.D0/36.D0*CA*
     &    omx**9*nfl + 4.D0/9.D0*CA*ln2*nfl - 12.D0/5.D0*CA*ln2*omx*nfl
     &     + 2.D0/15.D0*CA*ln2*omx**2*nfl - 4.D0/3.D0*CA*ln2*omx**3*nfl
     &     - 4.D0/3.D0*CA*ln2*omx**4*nfl - 4.D0/3.D0*CA*ln2*omx**5*nfl
     &     - 4.D0/3.D0*CA*ln2*omx**6*nfl - 4.D0/3.D0*CA*ln2*omx**7*nfl
     &     - 4.D0/3.D0*CA*ln2*omx**8*nfl - 4.D0/3.D0*CA*ln2*omx**9*nfl
     &     + 2477.D0/216.D0*CA**2 + 77831.D0/1800.D0*CA**2*omx - 887.D0/
     &    120.D0*CA**2*omx**2 + 517.D0/12.D0*CA**2*omx**3 + 3050357.D0/
     &    88200.D0*CA**2*omx**4 + 1378987.D0/44100.D0*CA**2*omx**5 + 
     &    149505029.D0/5080320.D0*CA**2*omx**6 + 51289423.D0/1814400.D0
     &    *CA**2*omx**7 + 56319150883.D0/2049062400.D0*CA**2*omx**8 + 
     &    165554049803.D0/6147187200.D0*CA**2*omx**9 )
      I2R(gg,rglr) = I2R(gg,rglr) + R2 * (  - 527393.D0/145152.D0*CA**2
     &    *omx**10 + 1162447.D0/241920.D0*CA**2*omx**11 + 2*CA**2*pisq
     &     - 4*CA**2*pisq*omx - 2*CA**2*pisq*omx**3 - 2*CA**2*pisq*
     &    omx**4 - 2*CA**2*pisq*omx**5 - 2*CA**2*pisq*omx**6 - 2*CA**2*
     &    pisq*omx**7 - 2*CA**2*pisq*omx**8 - 2*CA**2*pisq*omx**9 - 124.
     &    D0/9.D0*CA**2*ln2 - 584.D0/15.D0*CA**2*ln2*omx + 8*CA**2*ln2*
     &    omx**2 - 42*CA**2*ln2*omx**3 - 3404.D0/105.D0*CA**2*ln2*
     &    omx**4 - 3068.D0/105.D0*CA**2*ln2*omx**5 - 1735.D0/63.D0*
     &    CA**2*ln2*omx**6 - 1192.D0/45.D0*CA**2*ln2*omx**7 - 9921.D0/
     &    385.D0*CA**2*ln2*omx**8 - 87478.D0/3465.D0*CA**2*ln2*omx**9
     &     + 5.D0/9.D0*CA**2*ln2*omx**10 - 8.D0/3.D0*CA**2*ln2*omx**11
     &     - 23.D0/9.D0*CF*omx*nfl - 23.D0/9.D0*CF*omx**2*nfl - 46.D0/9.
     &    D0*CF*omx**3*nfl - 46.D0/9.D0*CF*omx**4*nfl - 46.D0/9.D0*CF*
     &    omx**5*nfl - 46.D0/9.D0*CF*omx**6*nfl - 46.D0/9.D0*CF*omx**7*
     &    nfl - 46.D0/9.D0*CF*omx**8*nfl - 46.D0/9.D0*CF*omx**9*nfl + 8.
     &    D0/3.D0*CF*ln2*omx*nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + R2 * ( 8.D0/3.D0*CF*ln2*omx**2*nfl
     &     + 16.D0/3.D0*CF*ln2*omx**3*nfl + 16.D0/3.D0*CF*ln2*omx**4*
     &    nfl + 16.D0/3.D0*CF*ln2*omx**5*nfl + 16.D0/3.D0*CF*ln2*omx**6
     &    *nfl + 16.D0/3.D0*CF*ln2*omx**7*nfl + 16.D0/3.D0*CF*ln2*
     &    omx**8*nfl + 16.D0/3.D0*CF*ln2*omx**9*nfl )
      I2R(gg,rglr) = I2R(gg,rglr) + picube * (  - 4*CA*omx**(-1)*nfl*
     &    BCGGCATFx + 4*CA*omx**(-1)*nfl*BCGGCATF1 + 4*CA*nfl*BCGGCATFx
     &     - 8*CA**2*omx**(-1)*BCGGCACAx + 8*CA**2*omx**(-1)*BCGGCACA1
     &     + 8*CA**2*BCGGCACAx - 4*CF*omx**(-1)*nfl*BCGGCFTFx + 4*CF*
     &    omx**(-1)*nfl*BCGGCFTF1 + 4*CF*nfl*BCGGCFTFx )
      I2R(gg,rglr) = I2R(gg,rglr) + lnR * ( 92.D0/9.D0*CA*nfl - 184.D0/
     &    9.D0*CA*omx*nfl - 92.D0/9.D0*CA*omx**3*nfl - 92.D0/9.D0*CA*
     &    omx**4*nfl - 92.D0/9.D0*CA*omx**5*nfl - 92.D0/9.D0*CA*omx**6*
     &    nfl - 92.D0/9.D0*CA*omx**7*nfl - 92.D0/9.D0*CA*omx**8*nfl - 
     &    92.D0/9.D0*CA*omx**9*nfl - 32.D0/3.D0*CA*ln2*nfl + 64.D0/3.D0
     &    *CA*ln2*omx*nfl + 32.D0/3.D0*CA*ln2*omx**3*nfl + 32.D0/3.D0*
     &    CA*ln2*omx**4*nfl + 32.D0/3.D0*CA*ln2*omx**5*nfl + 32.D0/3.D0
     &    *CA*ln2*omx**6*nfl + 32.D0/3.D0*CA*ln2*omx**7*nfl + 32.D0/3.D0
     &    *CA*ln2*omx**8*nfl + 32.D0/3.D0*CA*ln2*omx**9*nfl - 524.D0/9.D
     &    0*CA**2 + 1048.D0/9.D0*CA**2*omx + 524.D0/9.D0*CA**2*omx**3
     &     + 524.D0/9.D0*CA**2*omx**4 + 524.D0/9.D0*CA**2*omx**5 + 524.D
     &    0/9.D0*CA**2*omx**6 + 524.D0/9.D0*CA**2*omx**7 + 524.D0/9.D0*
     &    CA**2*omx**8 + 524.D0/9.D0*CA**2*omx**9 + 16.D0/3.D0*CA**2*
     &    pisq - 32.D0/3.D0*CA**2*pisq*omx - 16.D0/3.D0*CA**2*pisq*
     &    omx**3 - 16.D0/3.D0*CA**2*pisq*omx**4 - 16.D0/3.D0*CA**2*pisq
     &    *omx**5 )
      I2R(gg,rglr) = I2R(gg,rglr) + lnR * (  - 16.D0/3.D0*CA**2*pisq*
     &    omx**6 - 16.D0/3.D0*CA**2*pisq*omx**7 - 16.D0/3.D0*CA**2*pisq
     &    *omx**8 - 16.D0/3.D0*CA**2*pisq*omx**9 + 176.D0/3.D0*CA**2*
     &    ln2 - 352.D0/3.D0*CA**2*ln2*omx - 176.D0/3.D0*CA**2*ln2*
     &    omx**3 - 176.D0/3.D0*CA**2*ln2*omx**4 - 176.D0/3.D0*CA**2*ln2
     &    *omx**5 - 176.D0/3.D0*CA**2*ln2*omx**6 - 176.D0/3.D0*CA**2*
     &    ln2*omx**7 - 176.D0/3.D0*CA**2*ln2*omx**8 - 176.D0/3.D0*CA**2
     &    *ln2*omx**9 )
      I2R(gg,rglr) = I2R(gg,rglr) + 26.D0/9.D0*CA*nfl - 92.D0/9.D0*CA*
     &    lnpi*nfl + 184.D0/9.D0*CA*lnpi*omx*nfl + 92.D0/9.D0*CA*lnpi*
     &    omx**3*nfl + 92.D0/9.D0*CA*lnpi*omx**4*nfl + 92.D0/9.D0*CA*
     &    lnpi*omx**5*nfl + 92.D0/9.D0*CA*lnpi*omx**6*nfl + 92.D0/9.D0*
     &    CA*lnpi*omx**7*nfl + 92.D0/9.D0*CA*lnpi*omx**8*nfl + 92.D0/9.D
     &    0*CA*lnpi*omx**9*nfl - 116.D0/9.D0*CA*ln2*nfl + 184.D0/9.D0*
     &    CA*ln2*omx*nfl + 92.D0/9.D0*CA*ln2*omx**3*nfl + 92.D0/9.D0*CA
     &    *ln2*omx**4*nfl + 92.D0/9.D0*CA*ln2*omx**5*nfl + 92.D0/9.D0*
     &    CA*ln2*omx**6*nfl + 92.D0/9.D0*CA*ln2*omx**7*nfl + 92.D0/9.D0
     &    *CA*ln2*omx**8*nfl + 92.D0/9.D0*CA*ln2*omx**9*nfl + 32.D0/3.D0
     &    *CA*ln2*lnpi*nfl - 64.D0/3.D0*CA*ln2*lnpi*omx*nfl - 32.D0/3.D0
     &    *CA*ln2*lnpi*omx**3*nfl - 32.D0/3.D0*CA*ln2*lnpi*omx**4*nfl
     &     - 32.D0/3.D0*CA*ln2*lnpi*omx**5*nfl - 32.D0/3.D0*CA*ln2*lnpi
     &    *omx**6*nfl - 32.D0/3.D0*CA*ln2*lnpi*omx**7*nfl - 32.D0/3.D0*
     &    CA*ln2*lnpi*omx**8*nfl - 32.D0/3.D0*CA*ln2*lnpi*omx**9*nfl + 
     &    32.D0/3.D0*CA*ln2**2*nfl
      I2R(gg,rglr) = I2R(gg,rglr) - 64.D0/3.D0*CA*ln2**2*omx*nfl - 32.D0
     &    /3.D0*CA*ln2**2*omx**3*nfl - 32.D0/3.D0*CA*ln2**2*omx**4*nfl
     &     - 32.D0/3.D0*CA*ln2**2*omx**5*nfl - 32.D0/3.D0*CA*ln2**2*
     &    omx**6*nfl - 32.D0/3.D0*CA*ln2**2*omx**7*nfl - 32.D0/3.D0*CA*
     &    ln2**2*omx**8*nfl - 32.D0/3.D0*CA*ln2**2*omx**9*nfl - 26.D0/9.
     &    D0*CA**2 + 524.D0/9.D0*CA**2*lnpi - 1048.D0/9.D0*CA**2*lnpi*
     &    omx - 524.D0/9.D0*CA**2*lnpi*omx**3 - 524.D0/9.D0*CA**2*lnpi*
     &    omx**4 - 524.D0/9.D0*CA**2*lnpi*omx**5 - 524.D0/9.D0*CA**2*
     &    lnpi*omx**6 - 524.D0/9.D0*CA**2*lnpi*omx**7 - 524.D0/9.D0*
     &    CA**2*lnpi*omx**8 - 524.D0/9.D0*CA**2*lnpi*omx**9 - 16.D0/3.D0
     &    *CA**2*pisq*lnpi + 32.D0/3.D0*CA**2*pisq*lnpi*omx + 16.D0/3.D0
     &    *CA**2*pisq*lnpi*omx**3 + 16.D0/3.D0*CA**2*pisq*lnpi*omx**4
     &     + 16.D0/3.D0*CA**2*pisq*lnpi*omx**5 + 16.D0/3.D0*CA**2*pisq*
     &    lnpi*omx**6 + 16.D0/3.D0*CA**2*pisq*lnpi*omx**7 + 16.D0/3.D0*
     &    CA**2*pisq*lnpi*omx**8 + 16.D0/3.D0*CA**2*pisq*lnpi*omx**9 + 
     &    548.D0/9.D0*CA**2*ln2
      I2R(gg,rglr) = I2R(gg,rglr) - 1048.D0/9.D0*CA**2*ln2*omx - 524.D0/
     &    9.D0*CA**2*ln2*omx**3 - 524.D0/9.D0*CA**2*ln2*omx**4 - 524.D0/
     &    9.D0*CA**2*ln2*omx**5 - 524.D0/9.D0*CA**2*ln2*omx**6 - 524.D0/
     &    9.D0*CA**2*ln2*omx**7 - 524.D0/9.D0*CA**2*ln2*omx**8 - 524.D0/
     &    9.D0*CA**2*ln2*omx**9 - 176.D0/3.D0*CA**2*ln2*lnpi + 352.D0/3.
     &    D0*CA**2*ln2*lnpi*omx + 176.D0/3.D0*CA**2*ln2*lnpi*omx**3 + 
     &    176.D0/3.D0*CA**2*ln2*lnpi*omx**4 + 176.D0/3.D0*CA**2*ln2*
     &    lnpi*omx**5 + 176.D0/3.D0*CA**2*ln2*lnpi*omx**6 + 176.D0/3.D0
     &    *CA**2*ln2*lnpi*omx**7 + 176.D0/3.D0*CA**2*ln2*lnpi*omx**8 + 
     &    176.D0/3.D0*CA**2*ln2*lnpi*omx**9 - 16.D0/3.D0*CA**2*ln2*pisq
     &     + 32.D0/3.D0*CA**2*ln2*pisq*omx + 16.D0/3.D0*CA**2*ln2*pisq*
     &    omx**3 + 16.D0/3.D0*CA**2*ln2*pisq*omx**4 + 16.D0/3.D0*CA**2*
     &    ln2*pisq*omx**5 + 16.D0/3.D0*CA**2*ln2*pisq*omx**6 + 16.D0/3.D
     &    0*CA**2*ln2*pisq*omx**7 + 16.D0/3.D0*CA**2*ln2*pisq*omx**8 + 
     &    16.D0/3.D0*CA**2*ln2*pisq*omx**9 - 176.D0/3.D0*CA**2*ln2**2
     &     + 352.D0/3.D0*CA**2*ln2**2*omx
      I2R(gg,rglr) = I2R(gg,rglr) + 176.D0/3.D0*CA**2*ln2**2*omx**3 + 
     &    176.D0/3.D0*CA**2*ln2**2*omx**4 + 176.D0/3.D0*CA**2*ln2**2*
     &    omx**5 + 176.D0/3.D0*CA**2*ln2**2*omx**6 + 176.D0/3.D0*CA**2*
     &    ln2**2*omx**7 + 176.D0/3.D0*CA**2*ln2**2*omx**8 + 176.D0/3.D0
     &    *CA**2*ln2**2*omx**9
      endif
