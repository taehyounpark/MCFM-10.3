!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      if (x < xmin) then
      I2R(qg,rglr)= + R8 * (  - 35*CA*omx**(-8)*Li2mx + 533.D0/24.D0*CA
     &    *omx**(-8)*lnx + 35*CA*omx**(-8)*Li2omx - 307.D0/24.D0*CA*
     &    omx**(-7) + 35.D0/2.D0*CA*omx**(-7)*lnopxon2 + 220*CA*
     &    omx**(-7)*Li2mx - 422.D0/3.D0*CA*omx**(-7)*lnx - 220*CA*
     &    omx**(-7)*Li2omx + 1307.D0/16.D0*CA*omx**(-6) - 385.D0/4.D0*
     &    CA*omx**(-6)*lnopxon2 - 1255.D0/2.D0*CA*omx**(-6)*Li2mx + 
     &    6483.D0/16.D0*CA*omx**(-6)*lnx + 1255.D0/2.D0*CA*omx**(-6)*
     &    Li2omx - 8365.D0/36.D0*CA*omx**(-5) + 2863.D0/12.D0*CA*
     &    omx**(-5)*lnopxon2 + 6419.D0/6.D0*CA*omx**(-5)*Li2mx - 100859.
     &    D0/144.D0*CA*omx**(-5)*lnx - 6419.D0/6.D0*CA*omx**(-5)*Li2omx
     &     + 37169.D0/96.D0*CA*omx**(-4) - 1401.D0/4.D0*CA*omx**(-4)*
     &    lnopxon2 - 865739.D0/720.D0*CA*omx**(-4)*Li2mx + 6948889.D0/
     &    8640.D0*CA*omx**(-4)*lnx + 865739.D0/720.D0*CA*omx**(-4)*
     &    Li2omx - 7159903.D0/17280.D0*CA*omx**(-3) + 481079.D0/1440.D0
     &    *CA*omx**(-3)*lnopxon2 + 41441.D0/45.D0*CA*omx**(-3)*Li2mx - 
     &    459431.D0/720.D0*CA*omx**(-3)*lnx )
      I2R(qg,rglr) = I2R(qg,rglr) + R8 * (  - 41441.D0/45.D0*CA*
     &    omx**(-3)*Li2omx + 669827.D0/2304.D0*CA*omx**(-2) - 612473.D0/
     &    2880.D0*CA*omx**(-2)*lnopxon2 - 2055437.D0/4320.D0*CA*
     &    omx**(-2)*Li2mx + 1013009.D0/2880.D0*CA*omx**(-2)*lnx + 
     &    2055437.D0/4320.D0*CA*omx**(-2)*Li2omx - 2641049.D0/20736.D0*
     &    CA*omx**(-1) + 255829.D0/2880.D0*CA*omx**(-1)*lnopxon2 + 
     &    224587.D0/1440.D0*CA*omx**(-1)*Li2mx - 562813.D0/4320.D0*CA*
     &    omx**(-1)*lnx - 224587.D0/1440.D0*CA*omx**(-1)*Li2omx + 115.D0
     &    /512.D0*CA*opx**(-6) - 1373.D0/3072.D0*CA*opx**(-5) + 499.D0/
     &    1536.D0*CA*opx**(-4) - 1817.D0/18432.D0*CA*opx**(-3) - 931.D0/
     &    55296.D0*CA*opx**(-2) - 87673.D0/368640.D0*CA*opx**(-1) + 
     &    24727244027.D0/812851200.D0*CA - 7439.D0/360.D0*CA*lnopxon2
     &     - 7679.D0/270.D0*CA*Li2mx + 114697.D0/4320.D0*CA*lnx + 7679.D
     &    0/270.D0*CA*Li2omx - 149513417.D0/30105600.D0*CA*omx - 881.D0/
     &    1440.D0*CA*omx*lnopxon2 + 6839.D0/2160.D0*CA*omx*Li2mx + 1021.
     &    D0/1080.D0*CA*omx*lnx )
      I2R(qg,rglr) = I2R(qg,rglr) + R8 * (  - 6839.D0/2160.D0*CA*omx*
     &    Li2omx + 15259508777.D0/6502809600.D0*CA*omx**2 + 13.D0/9.D0*
     &    CA*omx**2*lnopxon2 - 13.D0/18.D0*CA*omx**2*Li2mx - 13.D0/9.D0
     &    *CA*omx**2*lnx + 13.D0/18.D0*CA*omx**2*Li2omx - 206603081.D0/
     &    1625702400.D0*CA*omx**3 - 26084869.D0/203212800.D0*CA*omx**4
     &     - 8029621.D0/406425600.D0*CA*omx**5 - 35.D0/12.D0*CA*pisq*
     &    omx**(-8) + 55.D0/3.D0*CA*pisq*omx**(-7) - 1255.D0/24.D0*CA*
     &    pisq*omx**(-6) + 6419.D0/72.D0*CA*pisq*omx**(-5) - 865739.D0/
     &    8640.D0*CA*pisq*omx**(-4) + 41441.D0/540.D0*CA*pisq*omx**(-3)
     &     - 2055437.D0/51840.D0*CA*pisq*omx**(-2) + 224587.D0/17280.D0
     &    *CA*pisq*omx**(-1) - 7679.D0/3240.D0*CA*pisq + 6839.D0/25920.D
     &    0*CA*pisq*omx - 13.D0/216.D0*CA*pisq*omx**2 + 35*CA*ln2*
     &    omx**(-7) - 405.D0/2.D0*CA*ln2*omx**(-6) + 3175.D0/6.D0*CA*
     &    ln2*omx**(-5) - 2462.D0/3.D0*CA*ln2*omx**(-4) + 596639.D0/720.
     &    D0*CA*ln2*omx**(-3) - 89217.D0/160.D0*CA*ln2*omx**(-2) + 
     &    349229.D0/1440.D0*CA*ln2*omx**(-1) )
      I2R(qg,rglr) = I2R(qg,rglr) + R8 * (  - 460171.D0/7560.D0*CA*ln2
     &     + 106501.D0/15120.D0*CA*ln2*omx - 16823.D0/15120.D0*CA*ln2*
     &    omx**2 + 3209.D0/15120.D0*CA*ln2*omx**3 + 37.D0/280.D0*CA*ln2
     &    *omx**4 + 23.D0/1260.D0*CA*ln2*omx**5 - 8698309.D0/
     &    12386304000.D0*CF + 165291043.D0/43352064000.D0*CF*omx - 
     &    41619673.D0/14450688000.D0*CF*omx**2 - 147337.D0/12902400.D0*
     &    CF*omx**3 + 4718647.D0/203212800.D0*CF*omx**4 - 1217729.D0/
     &    101606400.D0*CF*omx**5 - 1.D0/3110400.D0*CF*pisq + 1.D0/
     &    1555200.D0*CF*pisq*omx - 1.D0/1555200.D0*CF*pisq*omx**2 + 167.
     &    D0/345600.D0*CF*ln2 - 1243.D0/403200.D0*CF*ln2*omx + 4421.D0/
     &    1209600.D0*CF*ln2*omx**2 + 7.D0/1080.D0*CF*ln2*omx**3 - 11.D0/
     &    630.D0*CF*ln2*omx**4 + 13.D0/1260.D0*CF*ln2*omx**5 )
      I2R(qg,rglr) = I2R(qg,rglr) + R6 * (  - 40.D0/3.D0*CA*omx**(-6)*
     &    Li2mx + 74.D0/9.D0*CA*omx**(-6)*lnx + 40.D0/3.D0*CA*omx**(-6)
     &    *Li2omx - 46.D0/9.D0*CA*omx**(-5) + 20.D0/3.D0*CA*omx**(-5)*
     &    lnopxon2 + 188.D0/3.D0*CA*omx**(-5)*Li2mx - 355.D0/9.D0*CA*
     &    omx**(-5)*lnx - 188.D0/3.D0*CA*omx**(-5)*Li2omx + 24*CA*
     &    omx**(-4) - 26*CA*omx**(-4)*lnopxon2 - 391.D0/3.D0*CA*
     &    omx**(-4)*Li2mx + 3041.D0/36.D0*CA*omx**(-4)*lnx + 391.D0/3.D0
     &    *CA*omx**(-4)*Li2omx - 5177.D0/108.D0*CA*omx**(-3) + 805.D0/
     &    18.D0*CA*omx**(-3)*lnopxon2 + 155*CA*omx**(-3)*Li2mx - 3805.D0
     &    /36.D0*CA*omx**(-3)*lnx - 155*CA*omx**(-3)*Li2omx + 11153.D0/
     &    216.D0*CA*omx**(-2) - 1567.D0/36.D0*CA*omx**(-2)*lnopxon2 - 
     &    331.D0/3.D0*CA*omx**(-2)*Li2mx + 347.D0/4.D0*CA*omx**(-2)*lnx
     &     + 331.D0/3.D0*CA*omx**(-2)*Li2omx - 2783.D0/108.D0*CA*
     &    omx**(-1) + 899.D0/36.D0*CA*omx**(-1)*lnopxon2 + 43*CA*
     &    omx**(-1)*Li2mx - 1505.D0/36.D0*CA*omx**(-1)*lnx - 43*CA*
     &    omx**(-1)*Li2omx )
      I2R(qg,rglr) = I2R(qg,rglr) + R6 * ( 53.D0/96.D0*CA*opx**(-4) - 
     &    199.D0/576.D0*CA*opx**(-3) + 45.D0/128.D0*CA*opx**(-2) + 1133.
     &    D0/2304.D0*CA*opx**(-1) + 66619.D0/11520.D0*CA - 14.D0/9.D0*
     &    CA*lnopxon2 - 32.D0/3.D0*CA*Li2mx - 4.D0/9.D0*CA*lnx + 32.D0/
     &    3.D0*CA*Li2omx - 3439173.D0/313600.D0*CA*omx - 155.D0/18.D0*
     &    CA*omx*lnopxon2 + 22.D0/3.D0*CA*omx*Li2mx + 203.D0/18.D0*CA*
     &    omx*lnx - 22.D0/3.D0*CA*omx*Li2omx + 142427221.D0/16934400.D0
     &    *CA*omx**2 + 10.D0/3.D0*CA*omx**2*lnopxon2 - 10.D0/3.D0*CA*
     &    omx**2*Li2mx - 10.D0/3.D0*CA*omx**2*lnx + 10.D0/3.D0*CA*
     &    omx**2*Li2omx - 9479033.D0/8467200.D0*CA*omx**3 - 38833.D0/
     &    1058400.D0*CA*omx**4 - 10.D0/9.D0*CA*pisq*omx**(-6) + 47.D0/9.
     &    D0*CA*pisq*omx**(-5) - 391.D0/36.D0*CA*pisq*omx**(-4) + 155.D0
     &    /12.D0*CA*pisq*omx**(-3) - 331.D0/36.D0*CA*pisq*omx**(-2) + 
     &    43.D0/12.D0*CA*pisq*omx**(-1) - 8.D0/9.D0*CA*pisq + 11.D0/18.D
     &    0*CA*pisq*omx - 5.D0/18.D0*CA*pisq*omx**2 + 40.D0/3.D0*CA*ln2
     &    *omx**(-5) )
      I2R(qg,rglr) = I2R(qg,rglr) + R6 * (  - 56*CA*ln2*omx**(-4) + 931.
     &    D0/9.D0*CA*ln2*omx**(-3) - 1933.D0/18.D0*CA*ln2*omx**(-2) + 
     &    1139.D0/18.D0*CA*ln2*omx**(-1) - 155.D0/9.D0*CA*ln2 + 1471.D0/
     &    315.D0*CA*ln2*omx - 6667.D0/1260.D0*CA*ln2*omx**2 + 1537.D0/
     &    1260.D0*CA*ln2*omx**3 + 23.D0/630.D0*CA*ln2*omx**4 + 574541.D0
     &    /33868800.D0*CF - 120377.D0/3386880.D0*CF*omx + 6619.D0/
     &    1209600.D0*CF*omx**2 + 23813.D0/423360.D0*CF*omx**3 - 21989.D0
     &    /529200.D0*CF*omx**4 - 11.D0/1260.D0*CF*ln2 + 1.D0/42.D0*CF*
     &    ln2*omx - 1.D0/180.D0*CF*ln2*omx**2 - 5.D0/126.D0*CF*ln2*
     &    omx**3 + 23.D0/630.D0*CF*ln2*omx**4 )
      I2R(qg,rglr) = I2R(qg,rglr) + R4 * (  - 6*CA*omx**(-4)*Li2mx + 7.D
     &    0/2.D0*CA*omx**(-4)*lnx + 6*CA*omx**(-4)*Li2omx - 5.D0/2.D0*
     &    CA*omx**(-3) + 3*CA*omx**(-3)*lnopxon2 + 18*CA*omx**(-3)*
     &    Li2mx - 23.D0/2.D0*CA*omx**(-3)*lnx - 18*CA*omx**(-3)*Li2omx
     &     + 27.D0/4.D0*CA*omx**(-2) - 13.D0/2.D0*CA*omx**(-2)*lnopxon2
     &     - 133.D0/6.D0*CA*omx**(-2)*Li2mx + 199.D0/12.D0*CA*omx**(-2)
     &    *lnx + 133.D0/6.D0*CA*omx**(-2)*Li2omx - 143.D0/24.D0*CA*
     &    omx**(-1) + 79.D0/12.D0*CA*omx**(-1)*lnopxon2 + 21.D0/2.D0*CA
     &    *omx**(-1)*Li2mx - 23.D0/6.D0*CA*omx**(-1)*lnx - 21.D0/2.D0*
     &    CA*omx**(-1)*Li2omx + 37.D0/16.D0*CA*opx**(-2) + 167.D0/32.D0
     &    *CA*opx**(-1) + 1477.D0/288.D0*CA + 7.D0/4.D0*CA*lnopxon2 - 
     &    26.D0/3.D0*CA*Li2mx - 113.D0/12.D0*CA*lnx + 26.D0/3.D0*CA*
     &    Li2omx - 121661.D0/7200.D0*CA*omx - 2.D0/3.D0*CA*omx*lnopxon2
     &     + 37.D0/3.D0*CA*omx*Li2mx + 1.D0/2.D0*CA*omx*lnx - 37.D0/3.D0
     &    *CA*omx*Li2omx + 42599.D0/9600.D0*CA*omx**2 - 25.D0/6.D0*CA*
     &    omx**2*lnopxon2 )
      I2R(qg,rglr) = I2R(qg,rglr) + R4 * (  - 4*CA*omx**2*Li2mx + 25.D0/
     &    6.D0*CA*omx**2*lnx + 4*CA*omx**2*Li2omx + 8659.D0/14400.D0*CA
     &    *omx**3 - 1.D0/2.D0*CA*pisq*omx**(-4) + 3.D0/2.D0*CA*pisq*
     &    omx**(-3) - 133.D0/72.D0*CA*pisq*omx**(-2) + 7.D0/8.D0*CA*
     &    pisq*omx**(-1) - 13.D0/18.D0*CA*pisq + 37.D0/36.D0*CA*pisq*
     &    omx - 1.D0/3.D0*CA*pisq*omx**2 + 6*CA*ln2*omx**(-3) - 15*CA*
     &    ln2*omx**(-2) + 91.D0/6.D0*CA*ln2*omx**(-1) - 25.D0/6.D0*CA*
     &    ln2 + 42.D0/5.D0*CA*ln2*omx - 87.D0/10.D0*CA*ln2*omx**2 - 19.D
     &    0/30.D0*CA*ln2*omx**3 - 12841.D0/115200.D0*CF + 631.D0/3840.D0
     &    *CF*omx - 751.D0/19200.D0*CF*omx**2 - 541.D0/7200.D0*CF*
     &    omx**3 - 1.D0/432.D0*CF*pisq + 1.D0/216.D0*CF*pisq*omx - 1.D0/
     &    216.D0*CF*pisq*omx**2 + 71.D0/240.D0*CF*ln2 - 11.D0/24.D0*CF*
     &    ln2*omx + 11.D0/40.D0*CF*ln2*omx**2 + 1.D0/15.D0*CF*ln2*
     &    omx**3 )
      I2R(qg,rglr) = I2R(qg,rglr) + R2 * (  - 4*CA*omx**(-2)*Li2mx + 2*
     &    CA*omx**(-2)*lnx + 4*CA*omx**(-2)*Li2omx - 2*CA*omx**(-1) + 2
     &    *CA*omx**(-1)*lnopxon2 + 2*CA*omx**(-1)*Li2mx - 8*CA*
     &    omx**(-1)*lnx - 2*CA*omx**(-1)*Li2omx - 2*CA - 6*CA*lnopxon2
     &     + 6*CA*Li2mx + 22*CA*lnx - 6*CA*Li2omx + 157.D0/18.D0*CA*omx
     &     + 16*CA*omx*lnopxon2 - 16*CA*omx*Li2mx - 28*CA*omx*lnx + 16*
     &    CA*omx*Li2omx - 65.D0/36.D0*CA*omx**2 - 12*CA*omx**2*lnopxon2
     &     + 12*CA*omx**2*Li2mx + 12*CA*omx**2*lnx - 12*CA*omx**2*
     &    Li2omx - 1.D0/3.D0*CA*pisq*omx**(-2) + 1.D0/6.D0*CA*pisq*
     &    omx**(-1) + 1.D0/6.D0*CA*pisq - 2.D0/3.D0*CA*pisq*omx + 1.D0/
     &    3.D0*CA*pisq*omx**2 + 4*CA*ln2*omx**(-1) - 4*CA*ln2 + 22.D0/3.
     &    D0*CA*ln2*omx - 31.D0/3.D0*CA*ln2*omx**2 - 91.D0/24.D0*CF + 
     &    181.D0/36.D0*CF*omx - 35.D0/18.D0*CF*omx**2 - 1.D0/6.D0*CF*
     &    pisq + 1.D0/3.D0*CF*pisq*omx - 1.D0/3.D0*CF*pisq*omx**2 + 5.D0
     &    /18.D0*CF*ln4096 - 5.D0/9.D0*CF*ln4096*omx + 5.D0/18.D0*CF*
     &    ln4096*omx**2 )
      I2R(qg,rglr) = I2R(qg,rglr) + R2 * ( 7.D0/6.D0*CF*ln2 + 1.D0/3.D0
     &    *CF*ln2*omx )
      I2R(qg,rglr) = I2R(qg,rglr) + picube * ( 4*CA*omx**(-1)*BCQGCATF1
     &     - 4*CA*omx**(-1)*BCQGCATFx + 4*CA*BCQGCATFx + 4*CF*omx**(-1)
     &    *BCQGCFTF1 - 4*CF*omx**(-1)*BCQGCFTFx + 4*CF*BCQGCFTFx )
      I2R(qg,rglr) = I2R(qg,rglr) + lnR * ( 12*CF - 24*CF*omx + 24*CF*
     &    omx**2 - 4.D0/3.D0*CF*pisq + 8.D0/3.D0*CF*pisq*omx - 8.D0/3.D0
     &    *CF*pisq*omx**2 - 4.D0/3.D0*CF*ln512 + 8.D0/3.D0*CF*ln512*omx
     &     - 8.D0/3.D0*CF*ln512*omx**2 )
      I2R(qg,rglr) = I2R(qg,rglr) - 12*CF*lnpi + 24*CF*lnpi*omx - 24*CF
     &    *lnpi*omx**2 + 4.D0/3.D0*CF*pisq*lnpi - 8.D0/3.D0*CF*pisq*
     &    lnpi*omx + 8.D0/3.D0*CF*pisq*lnpi*omx**2 - 12*CF*ln2 + 24*CF*
     &    ln2*omx - 24*CF*ln2*omx**2 + 12*CF*ln2*lnpi - 24*CF*ln2*lnpi*
     &    omx + 24*CF*ln2*lnpi*omx**2 + 4.D0/3.D0*CF*ln2*pisq - 8.D0/3.D
     &    0*CF*ln2*pisq*omx + 8.D0/3.D0*CF*ln2*pisq*omx**2 + 12*CF*
     &    ln2**2 - 24*CF*ln2**2*omx + 24*CF*ln2**2*omx**2
      else
      I2R(qg,rglr)= + R8 * ( 26143.D0/54190080.D0*CA - 605800943.D0/
     &    97542144000.D0*CA*omx + 80811613.D0/10838016000.D0*CA*omx**2
     &     + 2561307120019.D0/47210397696000.D0*CA*omx**3 - 10843261241.
     &    D0/134120448000.D0*CA*omx**4 - 4096184161799.D0/
     &    6838763323392000.D0*CA*omx**5 + 930429342095299.D0/
     &    95742686527488000.D0*CA*omx**6 + 418738589462629.D0/
     &    71807014895616000.D0*CA*omx**7 + 188725409924693.D0/
     &    47871343263744000.D0*CA*omx**8 + 3782765592904494199.D0/
     &    1328142547509313536000.D0*CA*omx**9 - 1681468283.D0/
     &    12541132800.D0*CA*omx**10 + 23481991.D0/104509440.D0*CA*
     &    omx**11 - 1.D0/3024.D0*CA*ln2 + 4147.D0/907200.D0*CA*ln2*omx
     &     - 4847.D0/604800.D0*CA*ln2*omx**2 - 773453.D0/19958400.D0*CA
     &    *ln2*omx**3 + 14141.D0/201600.D0*CA*ln2*omx**4 - 2851.D0/
     &    579150.D0*CA*ln2*omx**5 - 1819177.D0/172972800.D0*CA*ln2*
     &    omx**6 - 8336059.D0/1556755200.D0*CA*ln2*omx**7 - 190433.D0/
     &    64864800.D0*CA*ln2*omx**8 )
      I2R(qg,rglr) = I2R(qg,rglr) + R8 * (  - 22338079.D0/13232419200.D0
     &    *CA*ln2*omx**9 + 1271.D0/4860.D0*CA*ln2*omx**10 - 13.D0/162.D0
     &    *CA*ln2*omx**11 - 8698309.D0/12386304000.D0*CF + 165291043.D0/
     &    43352064000.D0*CF*omx - 41619673.D0/14450688000.D0*CF*omx**2
     &     - 147337.D0/12902400.D0*CF*omx**3 + 4718647.D0/203212800.D0*
     &    CF*omx**4 - 1217729.D0/101606400.D0*CF*omx**5 - 1.D0/3110400.D
     &    0*CF*pisq + 1.D0/1555200.D0*CF*pisq*omx - 1.D0/1555200.D0*CF*
     &    pisq*omx**2 + 167.D0/345600.D0*CF*ln2 - 1243.D0/403200.D0*CF*
     &    ln2*omx + 4421.D0/1209600.D0*CF*ln2*omx**2 + 7.D0/1080.D0*CF*
     &    ln2*omx**3 - 11.D0/630.D0*CF*ln2*omx**4 + 13.D0/1260.D0*CF*
     &    ln2*omx**5 )
      I2R(qg,rglr) = I2R(qg,rglr) + R6 * (  - 443.D0/43200.D0*CA + 
     &    236311.D0/3386880.D0*CA*omx + 23243.D0/2419200.D0*CA*omx**2
     &     - 395028323.D0/1219276800.D0*CA*omx**3 + 198348163.D0/
     &    1219276800.D0*CA*omx**4 + 11386886827.D0/147532492800.D0*CA*
     &    omx**5 + 380302079.D0/10928332800.D0*CA*omx**6 + 
     &    1245384067583.D0/79785572106240.D0*CA*omx**7 + 744560139341.D0
     &    /132975953510400.D0*CA*omx**8 + 8512031713.D0/113979388723200.
     &    D0*CA*omx**9 - 10021157.D0/8709120.D0*CA*omx**10 + 1162447.D0/
     &    1741824.D0*CA*omx**11 + 1.D0/180.D0*CA*ln2 - 13.D0/252.D0*CA*
     &    ln2*omx + 1.D0/45.D0*CA*ln2*omx**2 + 467.D0/1890.D0*CA*ln2*
     &    omx**3 - 1213.D0/7560.D0*CA*ln2*omx**4 - 4847.D0/83160.D0*CA*
     &    ln2*omx**5 - 57.D0/3080.D0*CA*ln2*omx**6 - 955.D0/216216.D0*
     &    CA*ln2*omx**7 + 199.D0/180180.D0*CA*ln2*omx**8 + 64.D0/19305.D
     &    0*CA*ln2*omx**9 + 43.D0/108.D0*CA*ln2*omx**10 - 10.D0/27.D0*
     &    CA*ln2*omx**11 + 574541.D0/33868800.D0*CF - 120377.D0/3386880.
     &    D0*CF*omx )
      I2R(qg,rglr) = I2R(qg,rglr) + R6 * ( 6619.D0/1209600.D0*CF*omx**2
     &     + 23813.D0/423360.D0*CF*omx**3 - 21989.D0/529200.D0*CF*
     &    omx**4 - 11.D0/1260.D0*CF*ln2 + 1.D0/42.D0*CF*ln2*omx - 1.D0/
     &    180.D0*CF*ln2*omx**2 - 5.D0/126.D0*CF*ln2*omx**3 + 23.D0/630.D
     &    0*CF*ln2*omx**4 )
      I2R(qg,rglr) = I2R(qg,rglr) + R4 * ( 121.D0/288.D0*CA - 93011.D0/
     &    86400.D0*CA*omx - 883.D0/9600.D0*CA*omx**2 + 21627547.D0/
     &    16934400.D0*CA*omx**3 - 81503.D0/2822400.D0*CA*omx**4 - 
     &    334499.D0/2419200.D0*CA*omx**5 - 2876939.D0/19353600.D0*CA*
     &    omx**6 - 497200189.D0/3688312320.D0*CA*omx**7 - 948680749.D0/
     &    8196249600.D0*CA*omx**8 - 19416008963471.D0/199463930265600.D0
     &    *CA*omx**9 - 44143619.D0/34836480.D0*CA*omx**10 - 304201.D0/
     &    2903040.D0*CA*omx**11 - 1.D0/4.D0*CA*ln2 + 161.D0/180.D0*CA*
     &    ln2*omx - 17.D0/120.D0*CA*ln2*omx**2 - 2573.D0/2520.D0*CA*ln2
     &    *omx**3 + 6.D0/35.D0*CA*ln2*omx**4 + 11.D0/60.D0*CA*ln2*
     &    omx**5 + 109.D0/720.D0*CA*ln2*omx**6 + 4043.D0/33264.D0*CA*
     &    ln2*omx**7 + 227.D0/2310.D0*CA*ln2*omx**8 + 87119.D0/1081080.D
     &    0*CA*ln2*omx**9 + 47.D0/54.D0*CA*ln2*omx**10 - 4.D0/9.D0*CA*
     &    ln2*omx**11 - 12841.D0/115200.D0*CF + 631.D0/3840.D0*CF*omx
     &     - 751.D0/19200.D0*CF*omx**2 - 541.D0/7200.D0*CF*omx**3 - 1.D0
     &    /432.D0*CF*pisq )
      I2R(qg,rglr) = I2R(qg,rglr) + R4 * ( 1.D0/216.D0*CF*pisq*omx - 1.D
     &    0/216.D0*CF*pisq*omx**2 + 71.D0/240.D0*CF*ln2 - 11.D0/24.D0*
     &    CF*ln2*omx + 11.D0/40.D0*CF*ln2*omx**2 + 1.D0/15.D0*CF*ln2*
     &    omx**3 )
      I2R(qg,rglr) = I2R(qg,rglr) + R2 * ( 4*CA - 155.D0/12.D0*CA*omx
     &     + 839.D0/36.D0*CA*omx**2 - 4147.D0/800.D0*CA*omx**3 - 4877.D0
     &    /2400.D0*CA*omx**4 - 97373.D0/88200.D0*CA*omx**5 - 80491.D0/
     &    117600.D0*CA*omx**6 - 46556807.D0/101606400.D0*CA*omx**7 - 
     &    32894273.D0/101606400.D0*CA*omx**8 - 1954611499.D0/8196249600.
     &    D0*CA*omx**9 + 527393.D0/290304.D0*CA*omx**10 - 1162447.D0/
     &    483840.D0*CA*omx**11 - 1.D0/3.D0*CA*pisq + 2.D0/3.D0*CA*pisq*
     &    omx - 2.D0/3.D0*CA*pisq*omx**2 - 4*CA*ln2 + 13*CA*ln2*omx - 
     &    71.D0/3.D0*CA*ln2*omx**2 + 57.D0/10.D0*CA*ln2*omx**3 + 19.D0/
     &    10.D0*CA*ln2*omx**4 + 101.D0/105.D0*CA*ln2*omx**5 + 41.D0/70.D
     &    0*CA*ln2*omx**6 + 499.D0/1260.D0*CA*ln2*omx**7 + 361.D0/1260.D
     &    0*CA*ln2*omx**8 + 251.D0/1155.D0*CA*ln2*omx**9 - 5.D0/18.D0*
     &    CA*ln2*omx**10 + 4.D0/3.D0*CA*ln2*omx**11 - 91.D0/24.D0*CF + 
     &    181.D0/36.D0*CF*omx - 35.D0/18.D0*CF*omx**2 - 1.D0/6.D0*CF*
     &    pisq + 1.D0/3.D0*CF*pisq*omx - 1.D0/3.D0*CF*pisq*omx**2 + 5.D0
     &    /18.D0*CF*ln4096 )
      I2R(qg,rglr) = I2R(qg,rglr) + R2 * (  - 5.D0/9.D0*CF*ln4096*omx
     &     + 5.D0/18.D0*CF*ln4096*omx**2 + 7.D0/6.D0*CF*ln2 + 1.D0/3.D0
     &    *CF*ln2*omx )
      I2R(qg,rglr) = I2R(qg,rglr) + picube * ( 4*CA*omx**(-1)*BCQGCATF1
     &     - 4*CA*omx**(-1)*BCQGCATFx + 4*CA*BCQGCATFx + 4*CF*omx**(-1)
     &    *BCQGCFTF1 - 4*CF*omx**(-1)*BCQGCFTFx + 4*CF*BCQGCFTFx )
      I2R(qg,rglr) = I2R(qg,rglr) + lnR * ( 12*CF - 24*CF*omx + 24*CF*
     &    omx**2 - 4.D0/3.D0*CF*pisq + 8.D0/3.D0*CF*pisq*omx - 8.D0/3.D0
     &    *CF*pisq*omx**2 - 4.D0/3.D0*CF*ln512 + 8.D0/3.D0*CF*ln512*omx
     &     - 8.D0/3.D0*CF*ln512*omx**2 )
      I2R(qg,rglr) = I2R(qg,rglr) - 12*CF*lnpi + 24*CF*lnpi*omx - 24*CF
     &    *lnpi*omx**2 + 4.D0/3.D0*CF*pisq*lnpi - 8.D0/3.D0*CF*pisq*
     &    lnpi*omx + 8.D0/3.D0*CF*pisq*lnpi*omx**2 - 12*CF*ln2 + 24*CF*
     &    ln2*omx - 24*CF*ln2*omx**2 + 12*CF*ln2*lnpi - 24*CF*ln2*lnpi*
     &    omx + 24*CF*ln2*lnpi*omx**2 + 4.D0/3.D0*CF*ln2*pisq - 8.D0/3.D
     &    0*CF*ln2*pisq*omx + 8.D0/3.D0*CF*ln2*pisq*omx**2 + 12*CF*
     &    ln2**2 - 24*CF*ln2**2*omx + 24*CF*ln2**2*omx**2
      endif
