!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      if (x < xmin) then
      I2R(gq,rglr)= + R8 * (  - 70*CF*CA*omx**(-8)*Li2mx + 533.D0/12.D0
     &    *CF*CA*omx**(-8)*lnx + 70*CF*CA*omx**(-8)*Li2omx - 307.D0/12.D
     &    0*CF*CA*omx**(-7) + 35*CF*CA*omx**(-7)*lnopxon2 + 470*CF*CA*
     &    omx**(-7)*Li2mx - 3553.D0/12.D0*CF*CA*omx**(-7)*lnx - 470*CF*
     &    CA*omx**(-7)*Li2omx + 1429.D0/8.D0*CF*CA*omx**(-6) - 415.D0/2.
     &    D0*CF*CA*omx**(-6)*lnopxon2 - 1435*CF*CA*omx**(-6)*Li2mx + 
     &    7191.D0/8.D0*CF*CA*omx**(-6)*lnx + 1435*CF*CA*omx**(-6)*
     &    Li2omx - 40039.D0/72.D0*CF*CA*omx**(-5) + 1661.D0/3.D0*CF*CA*
     &    omx**(-5)*lnopxon2 + 7916.D0/3.D0*CF*CA*omx**(-5)*Li2mx - 
     &    59197.D0/36.D0*CF*CA*omx**(-5)*lnx - 7916.D0/3.D0*CF*CA*
     &    omx**(-5)*Li2omx + 49205.D0/48.D0*CF*CA*omx**(-4) - 880*CF*CA
     &    *omx**(-4)*lnopxon2 - 1169219.D0/360.D0*CF*CA*omx**(-4)*Li2mx
     &     + 8701189.D0/4320.D0*CF*CA*omx**(-4)*lnx + 1169219.D0/360.D0
     &    *CF*CA*omx**(-4)*Li2omx - 10747843.D0/8640.D0*CF*CA*omx**(-3)
     &     + 663599.D0/720.D0*CF*CA*omx**(-3)*lnopxon2 + 1007519.D0/360.
     &    D0*CF*CA*omx**(-3)*Li2mx )
      I2R(gq,rglr) = I2R(gq,rglr) + R8 * (  - 498319.D0/288.D0*CF*CA*
     &    omx**(-3)*lnx - 1007519.D0/360.D0*CF*CA*omx**(-3)*Li2omx + 
     &    1985723.D0/1920.D0*CF*CA*omx**(-2) - 949199.D0/1440.D0*CF*CA*
     &    omx**(-2)*lnopxon2 - 3685193.D0/2160.D0*CF*CA*omx**(-2)*Li2mx
     &     + 304043.D0/288.D0*CF*CA*omx**(-2)*lnx + 3685193.D0/2160.D0*
     &    CF*CA*omx**(-2)*Li2omx - 7647787.D0/12960.D0*CF*CA*omx**(-1)
     &     + 231599.D0/720.D0*CF*CA*omx**(-1)*lnopxon2 + 86533.D0/120.D0
     &    *CF*CA*omx**(-1)*Li2mx - 43573.D0/96.D0*CF*CA*omx**(-1)*lnx
     &     - 86533.D0/120.D0*CF*CA*omx**(-1)*Li2omx - 11441.D0/33177600.
     &    D0*CF*CA*x**(-1) - 115.D0/256.D0*CF*CA*opx**(-6) + 1271.D0/
     &    1536.D0*CF*CA*opx**(-5) - 875.D0/1536.D0*CF*CA*opx**(-4) + 
     &    1081.D0/9216.D0*CF*CA*opx**(-3) + 847.D0/6912.D0*CF*CA*
     &    opx**(-2) + 185911.D0/552960.D0*CF*CA*opx**(-1) + 
     &    355211709841.D0/1625702400.D0*CF*CA - 155999.D0/1440.D0*CF*CA
     &    *lnopxon2 - 422999.D0/2160.D0*CF*CA*Li2mx + 150359.D0/1080.D0
     &    *CF*CA*lnx )
      I2R(gq,rglr) = I2R(gq,rglr) + R8 * ( 422999.D0/2160.D0*CF*CA*
     &    Li2omx - 1942098398873.D0/48771072000.D0*CF*CA*omx + 254.D0/9.
     &    D0*CF*CA*omx*lnopxon2 + 235.D0/9.D0*CF*CA*omx*Li2mx - 33*CF*
     &    CA*omx*lnx - 235.D0/9.D0*CF*CA*omx*Li2omx - 1445700449.D0/
     &    812851200.D0*CF*CA*omx**2 - 47.D0/9.D0*CF*CA*omx**2*lnopxon2
     &     + 2.D0/9.D0*CF*CA*omx**2*Li2mx + 47.D0/9.D0*CF*CA*omx**2*lnx
     &     - 2.D0/9.D0*CF*CA*omx**2*Li2omx + 97290457.D0/406425600.D0*
     &    CF*CA*omx**3 + 42147523.D0/203212800.D0*CF*CA*omx**4 + 346247.
     &    D0/4838400.D0*CF*CA*omx**5 - 35.D0/6.D0*CF*CA*pisq*omx**(-8)
     &     + 235.D0/6.D0*CF*CA*pisq*omx**(-7) - 1435.D0/12.D0*CF*CA*
     &    pisq*omx**(-6) + 1979.D0/9.D0*CF*CA*pisq*omx**(-5) - 1169219.D
     &    0/4320.D0*CF*CA*pisq*omx**(-4) + 1007519.D0/4320.D0*CF*CA*
     &    pisq*omx**(-3) - 3685193.D0/25920.D0*CF*CA*pisq*omx**(-2) + 
     &    86533.D0/1440.D0*CF*CA*pisq*omx**(-1) - 422999.D0/25920.D0*CF
     &    *CA*pisq + 235.D0/108.D0*CF*CA*pisq*omx + 1.D0/54.D0*CF*CA*
     &    pisq*omx**2 )
      I2R(gq,rglr) = I2R(gq,rglr) + R8 * ( 70*CF*CA*ln2*omx**(-7) - 435
     &    *CF*CA*ln2*omx**(-6) + 3670.D0/3.D0*CF*CA*ln2*omx**(-5) - 
     &    6181.D0/3.D0*CF*CA*ln2*omx**(-4) + 829199.D0/360.D0*CF*CA*ln2
     &    *omx**(-3) - 142231.D0/80.D0*CF*CA*ln2*omx**(-2) + 341419.D0/
     &    360.D0*CF*CA*ln2*omx**(-1) - 5079467.D0/15120.D0*CF*CA*ln2 + 
     &    30902947.D0/453600.D0*CF*CA*ln2*omx - 24869.D0/7560.D0*CF*CA*
     &    ln2*omx**2 - 829.D0/1890.D0*CF*CA*ln2*omx**3 - 149.D0/630.D0*
     &    CF*CA*ln2*omx**4 - 1.D0/15.D0*CF*CA*ln2*omx**5 + 4987.D0/
     &    265420800.D0*CF**2*x**(-1) + 65864681.D0/130056192000.D0*
     &    CF**2 - 1018205861.D0/390168576000.D0*CF**2*omx - 1429121.D0/
     &    812851200.D0*CF**2*omx**2 + 35692463.D0/812851200.D0*CF**2*
     &    omx**3 - 3535121.D0/40642560.D0*CF**2*omx**4 + 3174257.D0/
     &    67737600.D0*CF**2*omx**5 - 1.D0/777600.D0*CF**2*pisq*x**(-1)
     &     + 1.D0/1555200.D0*CF**2*pisq + 1.D0/1555200.D0*CF**2*pisq*
     &    omx - 1.D0/86400.D0*CF**2*ln2*x**(-1) - 139.D0/403200.D0*
     &    CF**2*ln2 )
      I2R(gq,rglr) = I2R(gq,rglr) + R8 * ( 1433.D0/725760.D0*CF**2*ln2*
     &    omx + 1.D0/1890.D0*CF**2*ln2*omx**2 - 29.D0/945.D0*CF**2*ln2*
     &    omx**3 + 17.D0/252.D0*CF**2*ln2*omx**4 - 17.D0/420.D0*CF**2*
     &    ln2*omx**5 )
      I2R(gq,rglr) = I2R(gq,rglr) + R6 * (  - 80.D0/3.D0*CF*CA*
     &    omx**(-6)*Li2mx + 148.D0/9.D0*CF*CA*omx**(-6)*lnx + 80.D0/3.D0
     &    *CF*CA*omx**(-6)*Li2omx - 92.D0/9.D0*CF*CA*omx**(-5) + 40.D0/
     &    3.D0*CF*CA*omx**(-5)*lnopxon2 + 424.D0/3.D0*CF*CA*omx**(-5)*
     &    Li2mx - 770.D0/9.D0*CF*CA*omx**(-5)*lnx - 424.D0/3.D0*CF*CA*
     &    omx**(-5)*Li2omx + 172.D0/3.D0*CF*CA*omx**(-4) - 60*CF*CA*
     &    omx**(-4)*lnopxon2 - 998.D0/3.D0*CF*CA*omx**(-4)*Li2mx + 3581.
     &    D0/18.D0*CF*CA*omx**(-4)*lnx + 998.D0/3.D0*CF*CA*omx**(-4)*
     &    Li2omx - 7409.D0/54.D0*CF*CA*omx**(-3) + 1057.D0/9.D0*CF*CA*
     &    omx**(-3)*lnopxon2 + 1384.D0/3.D0*CF*CA*omx**(-3)*Li2mx - 817.
     &    D0/3.D0*CF*CA*omx**(-3)*lnx - 1384.D0/3.D0*CF*CA*omx**(-3)*
     &    Li2omx + 20363.D0/108.D0*CF*CA*omx**(-2) - 2401.D0/18.D0*CF*
     &    CA*omx**(-2)*lnopxon2 - 1240.D0/3.D0*CF*CA*omx**(-2)*Li2mx + 
     &    2123.D0/9.D0*CF*CA*omx**(-2)*lnx + 1240.D0/3.D0*CF*CA*
     &    omx**(-2)*Li2omx - 18103.D0/108.D0*CF*CA*omx**(-1) + 841.D0/9.
     &    D0*CF*CA*omx**(-1)*lnopxon2 )
      I2R(gq,rglr) = I2R(gq,rglr) + R6 * ( 240*CF*CA*omx**(-1)*Li2mx - 
     &    1297.D0/9.D0*CF*CA*omx**(-1)*lnx - 240*CF*CA*omx**(-1)*Li2omx
     &     + 7.D0/384.D0*CF*CA*x**(-1) - 53.D0/48.D0*CF*CA*opx**(-4) + 
     &    133.D0/288.D0*CF*CA*opx**(-3) - 419.D0/576.D0*CF*CA*opx**(-2)
     &     - 1667.D0/1152.D0*CF*CA*opx**(-1) + 1365713.D0/17280.D0*CF*
     &    CA - 985.D0/18.D0*CF*CA*lnopxon2 - 214.D0/3.D0*CF*CA*Li2mx + 
     &    1513.D0/18.D0*CF*CA*lnx + 214.D0/3.D0*CF*CA*Li2omx + 
     &    101055709.D0/8467200.D0*CF*CA*omx + 100.D0/3.D0*CF*CA*omx*
     &    lnopxon2 - 20.D0/3.D0*CF*CA*omx*Li2mx - 128.D0/3.D0*CF*CA*omx
     &    *lnx + 20.D0/3.D0*CF*CA*omx*Li2omx - 89921651.D0/4233600.D0*
     &    CF*CA*omx**2 - 28.D0/3.D0*CF*CA*omx**2*lnopxon2 + 8*CF*CA*
     &    omx**2*Li2mx + 28.D0/3.D0*CF*CA*omx**2*lnx - 8*CF*CA*omx**2*
     &    Li2omx + 1216207.D0/529200.D0*CF*CA*omx**3 + 8987.D0/29400.D0
     &    *CF*CA*omx**4 - 20.D0/9.D0*CF*CA*pisq*omx**(-6) + 106.D0/9.D0
     &    *CF*CA*pisq*omx**(-5) - 499.D0/18.D0*CF*CA*pisq*omx**(-4) + 
     &    346.D0/9.D0*CF*CA*pisq*omx**(-3) )
      I2R(gq,rglr) = I2R(gq,rglr) + R6 * (  - 310.D0/9.D0*CF*CA*pisq*
     &    omx**(-2) + 20*CF*CA*pisq*omx**(-1) - 107.D0/18.D0*CF*CA*pisq
     &     - 5.D0/9.D0*CF*CA*pisq*omx + 2.D0/3.D0*CF*CA*pisq*omx**2 + 
     &    80.D0/3.D0*CF*CA*ln2*omx**(-5) - 128*CF*CA*ln2*omx**(-4) + 
     &    2438.D0/9.D0*CF*CA*ln2*omx**(-3) - 3019.D0/9.D0*CF*CA*ln2*
     &    omx**(-2) + 2372.D0/9.D0*CF*CA*ln2*omx**(-1) - 1150.D0/9.D0*
     &    CF*CA*ln2 + 12947.D0/630.D0*CF*CA*ln2*omx + 7729.D0/630.D0*CF
     &    *CA*ln2*omx**2 - 857.D0/315.D0*CF*CA*ln2*omx**3 - 32.D0/105.D0
     &    *CF*CA*ln2*omx**4 + 1.D0/6912.D0*CF**2*x**(-1) - 122533.D0/
     &    16934400.D0*CF**2 + 128677.D0/5644800.D0*CF**2*omx + 26413.D0/
     &    705600.D0*CF**2*omx**2 - 288313.D0/1058400.D0*CF**2*omx**3 + 
     &    58841.D0/264600.D0*CF**2*omx**4 + 1.D0/210.D0*CF**2*ln2 - 11.D
     &    0/630.D0*CF**2*ln2*omx - 2.D0/105.D0*CF**2*ln2*omx**2 + 64.D0/
     &    315.D0*CF**2*ln2*omx**3 - 62.D0/315.D0*CF**2*ln2*omx**4 )
      I2R(gq,rglr) = I2R(gq,rglr) + R4 * (  - 12*CF*CA*omx**(-4)*Li2mx
     &     + 7*CF*CA*omx**(-4)*lnx + 12*CF*CA*omx**(-4)*Li2omx - 5*CF*
     &    CA*omx**(-3) + 6*CF*CA*omx**(-3)*lnopxon2 + 48*CF*CA*
     &    omx**(-3)*Li2mx - 26*CF*CA*omx**(-3)*lnx - 48*CF*CA*omx**(-3)
     &    *Li2omx + 45.D0/2.D0*CF*CA*omx**(-2) - 19*CF*CA*omx**(-2)*
     &    lnopxon2 - 241.D0/3.D0*CF*CA*omx**(-2)*Li2mx + 253.D0/6.D0*CF
     &    *CA*omx**(-2)*lnx + 241.D0/3.D0*CF*CA*omx**(-2)*Li2omx - 449.D
     &    0/12.D0*CF*CA*omx**(-1) + 133.D0/6.D0*CF*CA*omx**(-1)*
     &    lnopxon2 + 80*CF*CA*omx**(-1)*Li2mx - 349.D0/6.D0*CF*CA*
     &    omx**(-1)*lnx - 80*CF*CA*omx**(-1)*Li2omx + 55.D0/288.D0*CF*
     &    CA*x**(-1) - 37.D0/8.D0*CF*CA*opx**(-2) - 197.D0/16.D0*CF*CA*
     &    opx**(-1) + 789.D0/32.D0*CF*CA - 27*CF*CA*lnopxon2 - 37.D0/3.D
     &    0*CF*CA*Li2mx + 403.D0/6.D0*CF*CA*lnx + 37.D0/3.D0*CF*CA*
     &    Li2omx + 2141503.D0/43200.D0*CF*CA*omx + 145.D0/6.D0*CF*CA*
     &    omx*lnopxon2 - 40*CF*CA*omx*Li2mx - 193.D0/6.D0*CF*CA*omx*lnx
     &     + 40*CF*CA*omx*Li2omx )
      I2R(gq,rglr) = I2R(gq,rglr) + R4 * (  - 214141.D0/7200.D0*CF*CA*
     &    omx**2 + 16*CF*CA*omx**2*Li2mx - 16*CF*CA*omx**2*Li2omx + 
     &    4541.D0/3600.D0*CF*CA*omx**3 - CF*CA*pisq*omx**(-4) + 4*CF*CA
     &    *pisq*omx**(-3) - 241.D0/36.D0*CF*CA*pisq*omx**(-2) + 20.D0/3.
     &    D0*CF*CA*pisq*omx**(-1) - 37.D0/36.D0*CF*CA*pisq - 10.D0/3.D0
     &    *CF*CA*pisq*omx + 4.D0/3.D0*CF*CA*pisq*omx**2 + 12*CF*CA*ln2*
     &    omx**(-3) - 42*CF*CA*ln2*omx**(-2) + 181.D0/3.D0*CF*CA*ln2*
     &    omx**(-1) + 2.D0/3.D0*CF*CA*ln2*x**(-1) - 53*CF*CA*ln2 - 913.D
     &    0/90.D0*CF*CA*ln2*omx + 451.D0/15.D0*CF*CA*ln2*omx**2 - 22.D0/
     &    15.D0*CF*CA*ln2*omx**3 + 815.D0/1152.D0*CF**2*x**(-1) - 1599.D
     &    0/6400.D0*CF**2 - 100837.D0/172800.D0*CF**2*omx - 3859.D0/
     &    7200.D0*CF**2*omx**2 + 8341.D0/7200.D0*CF**2*omx**3 - 1.D0/
     &    108.D0*CF**2*pisq*x**(-1) + 1.D0/216.D0*CF**2*pisq + 1.D0/216.
     &    D0*CF**2*pisq*omx - 3.D0/4.D0*CF**2*ln2*x**(-1) + 41.D0/120.D0
     &    *CF**2*ln2 + 187.D0/360.D0*CF**2*ln2*omx + 4.D0/15.D0*CF**2*
     &    ln2*omx**2 )
      I2R(gq,rglr) = I2R(gq,rglr) + R4 * (  - 16.D0/15.D0*CF**2*ln2*
     &    omx**3 )
      I2R(gq,rglr) = I2R(gq,rglr) + R2 * (  - 8*CF*CA*omx**(-2)*Li2mx
     &     + 4*CF*CA*omx**(-2)*lnx + 8*CF*CA*omx**(-2)*Li2omx - 4*CF*CA
     &    *omx**(-1) + 4*CF*CA*omx**(-1)*lnopxon2 + 28*CF*CA*omx**(-1)*
     &    Li2mx - 4*CF*CA*omx**(-1)*lnx - 28*CF*CA*omx**(-1)*Li2omx - 
     &    11.D0/3.D0*CF*CA*x**(-1) + 59.D0/3.D0*CF*CA - 24*CF*CA*Li2mx
     &     - 12*CF*CA*lnx + 24*CF*CA*Li2omx - 163.D0/18.D0*CF*CA*omx - 
     &    12*CF*CA*omx*lnopxon2 - 20*CF*CA*omx*Li2mx + 28*CF*CA*omx*lnx
     &     + 20*CF*CA*omx*Li2omx - 8*CF*CA*omx**2 + 16*CF*CA*omx**2*
     &    lnopxon2 - 16*CF*CA*omx**2*lnx - 2.D0/3.D0*CF*CA*pisq*
     &    omx**(-2) + 7.D0/3.D0*CF*CA*pisq*omx**(-1) - 4.D0/3.D0*CF*CA*
     &    pisq*x**(-1) - 4.D0/3.D0*CF*CA*pisq - CF*CA*pisq*omx + 8*CF*
     &    CA*ln2*omx**(-1) + 4*CF*CA*ln2*x**(-1) - 20*CF*CA*ln2 - 10.D0/
     &    3.D0*CF*CA*ln2*omx + 24*CF*CA*ln2*omx**2 + 145.D0/6.D0*CF**2*
     &    x**(-1) - 63.D0/4.D0*CF**2 - 331.D0/36.D0*CF**2*omx + 8*CF**2
     &    *omx**2 - 2.D0/3.D0*CF**2*pisq*x**(-1) + 1.D0/3.D0*CF**2*pisq
     &     + 1.D0/3.D0*CF**2*pisq*omx )
      I2R(gq,rglr) = I2R(gq,rglr) + R2 * (  - 22*CF**2*ln2*x**(-1) + 15
     &    *CF**2*ln2 + 25.D0/3.D0*CF**2*ln2*omx - 8*CF**2*ln2*omx**2 )
      I2R(gq,rglr) = I2R(gq,rglr) + picube * ( 8*CF*CA*omx**(-1)*
     &    BCGQCFCA1 - 8*CF*CA*omx**(-1)*BCGQCFCAx + 8*CF*CA*BCGQCFCAx
     &     + 8*CF**2*omx**(-1)*BCGQCFCF1 - 8*CF**2*omx**(-1)*BCGQCFCFx
     &     + 8*CF**2*BCGQCFCFx )
      I2R(gq,rglr) = I2R(gq,rglr) + lnR * ( 48*CF**2*x**(-1) - 24*CF**2
     &     - 24*CF**2*omx - 16.D0/3.D0*CF**2*pisq*x**(-1) + 8.D0/3.D0*
     &    CF**2*pisq + 8.D0/3.D0*CF**2*pisq*omx - 16.D0/3.D0*CF**2*
     &    ln512*x**(-1) + 8.D0/3.D0*CF**2*ln512 + 8.D0/3.D0*CF**2*ln512
     &    *omx )
      I2R(gq,rglr) = I2R(gq,rglr) - 48*CF**2*lnpi*x**(-1) + 24*CF**2*
     &    lnpi + 24*CF**2*lnpi*omx + 16.D0/3.D0*CF**2*pisq*lnpi*x**(-1)
     &     - 8.D0/3.D0*CF**2*pisq*lnpi - 8.D0/3.D0*CF**2*pisq*lnpi*omx
     &     - 48*CF**2*ln2*x**(-1) + 24*CF**2*ln2 + 24*CF**2*ln2*omx + 
     &    48*CF**2*ln2*lnpi*x**(-1) - 24*CF**2*ln2*lnpi - 24*CF**2*ln2*
     &    lnpi*omx + 16.D0/3.D0*CF**2*ln2*pisq*x**(-1) - 8.D0/3.D0*
     &    CF**2*ln2*pisq - 8.D0/3.D0*CF**2*ln2*pisq*omx + 48*CF**2*
     &    ln2**2*x**(-1) - 24*CF**2*ln2**2 - 24*CF**2*ln2**2*omx
      else
      I2R(gq,rglr)= + R8 * (  - 26143.D0/27095040.D0*CF*CA + 50912017.D0
     &    /24385536000.D0*CF*CA*omx + 684541541.D0/16257024000.D0*CF*CA
     &    *omx**2 - 3313243403669.D0/23605198848000.D0*CF*CA*omx**3 + 
     &    819311169179.D0/7868399616000.D0*CF*CA*omx**4 + 
     &    752999051400907.D0/23935671631872000.D0*CF*CA*omx**5 - 
     &    949461192225683.D0/47871343263744000.D0*CF*CA*omx**6 - 
     &    1690977175572821.D0/143614029791232000.D0*CF*CA*omx**7 - 
     &    7874077199017.D0/997319651328000.D0*CF*CA*omx**8 - 
     &    3733510562518919849.D0/664071273754656768000.D0*CF*CA*omx**9
     &     + 8562227.D0/13063680.D0*CF*CA*omx**10 - 15649297.D0/
     &    26127360.D0*CF*CA*omx**11 + 1.D0/1512.D0*CF*CA*ln2 - 743.D0/
     &    226800.D0*CF*CA*ln2*omx - 2531.D0/100800.D0*CF*CA*ln2*omx**2
     &     + 551249.D0/4989600.D0*CF*CA*ln2*omx**3 - 499577.D0/4989600.D
     &    0*CF*CA*ln2*omx**4 - 2499313.D0/129729600.D0*CF*CA*ln2*omx**5
     &     + 1717439.D0/86486400.D0*CF*CA*ln2*omx**6 + 3609157.D0/
     &    389188800.D0*CF*CA*ln2*omx**7 )
      I2R(gq,rglr) = I2R(gq,rglr) + R8 * ( 132151.D0/28828800.D0*CF*CA*
     &    ln2*omx**8 + 7584757.D0/3308104800.D0*CF*CA*ln2*omx**9 + 949.D
     &    0/324.D0*CF*CA*ln2*omx**10 + 2.D0/81.D0*CF*CA*ln2*omx**11 + 
     &    22769437.D0/43352064000.D0*CF**2 - 1010874971.D0/390168576000.
     &    D0*CF**2*omx - 22621573.D0/13005619200.D0*CF**2*omx**2 + 
     &    63480419.D0/1445068800.D0*CF**2*omx**3 - 376998119.D0/
     &    4335206400.D0*CF**2*omx**4 + 609701707.D0/13005619200.D0*
     &    CF**2*omx**5 + 4987.D0/265420800.D0*CF**2*omx**6 + 4987.D0/
     &    265420800.D0*CF**2*omx**7 + 4987.D0/265420800.D0*CF**2*omx**8
     &     + 4987.D0/265420800.D0*CF**2*omx**9 - 1.D0/1555200.D0*CF**2*
     &    pisq - 1.D0/1555200.D0*CF**2*pisq*omx - 1.D0/777600.D0*CF**2*
     &    pisq*omx**2 - 1.D0/777600.D0*CF**2*pisq*omx**3 - 1.D0/777600.D
     &    0*CF**2*pisq*omx**4 - 1.D0/777600.D0*CF**2*pisq*omx**5 - 1.D0/
     &    777600.D0*CF**2*pisq*omx**6 - 1.D0/777600.D0*CF**2*pisq*
     &    omx**7 - 1.D0/777600.D0*CF**2*pisq*omx**8 - 1.D0/777600.D0*
     &    CF**2*pisq*omx**9 )
      I2R(gq,rglr) = I2R(gq,rglr) + R8 * (  - 431.D0/1209600.D0*CF**2*
     &    ln2 + 7123.D0/3628800.D0*CF**2*ln2*omx + 313.D0/604800.D0*
     &    CF**2*ln2*omx**2 - 2063.D0/67200.D0*CF**2*ln2*omx**3 + 40793.D
     &    0/604800.D0*CF**2*ln2*omx**4 - 24487.D0/604800.D0*CF**2*ln2*
     &    omx**5 - 1.D0/86400.D0*CF**2*ln2*omx**6 - 1.D0/86400.D0*CF**2
     &    *ln2*omx**7 - 1.D0/86400.D0*CF**2*ln2*omx**8 - 1.D0/86400.D0*
     &    CF**2*ln2*omx**9 )
      I2R(gq,rglr) = I2R(gq,rglr) + R6 * ( 443.D0/21600.D0*CF*CA + 
     &    27679.D0/470400.D0*CF*CA*omx - 1098187.D0/2822400.D0*CF*CA*
     &    omx**2 + 42587707.D0/87091200.D0*CF*CA*omx**3 - 3662467.D0/
     &    43545600.D0*CF*CA*omx**4 - 1058294137.D0/9220780800.D0*CF*CA*
     &    omx**5 - 3772841929.D0/147532492800.D0*CF*CA*omx**6 + 
     &    2355631566227.D0/199463930265600.D0*CF*CA*omx**7 + 726822119.D
     &    0/24625176576.D0*CF*CA*omx**8 + 606720251617.D0/
     &    15957114421248.D0*CF*CA*omx**9 + 16279.D0/4860.D0*CF*CA*
     &    omx**10 - 1269757.D0/725760.D0*CF*CA*omx**11 - 1.D0/90.D0*CF*
     &    CA*ln2 - 1.D0/70.D0*CF*CA*ln2*omx + 59.D0/210.D0*CF*CA*ln2*
     &    omx**2 - 58.D0/135.D0*CF*CA*ln2*omx**3 + 59.D0/540.D0*CF*CA*
     &    ln2*omx**4 + 1024.D0/10395.D0*CF*CA*ln2*omx**5 + 163.D0/10395.
     &    D0*CF*CA*ln2*omx**6 - 1366.D0/135135.D0*CF*CA*ln2*omx**7 - 
     &    655.D0/36036.D0*CF*CA*ln2*omx**8 - 1081.D0/54054.D0*CF*CA*ln2
     &    *omx**9 + 7.D0/27.D0*CF*CA*ln2*omx**10 + 8.D0/9.D0*CF*CA*ln2*
     &    omx**11 )
      I2R(gq,rglr) = I2R(gq,rglr) + R6 * (  - 120083.D0/16934400.D0*
     &    CF**2 + 388481.D0/16934400.D0*CF**2*omx + 318181.D0/8467200.D0
     &    *CF**2*omx**2 - 2305279.D0/8467200.D0*CF**2*omx**3 + 1884137.D
     &    0/8467200.D0*CF**2*omx**4 + 1.D0/6912.D0*CF**2*omx**5 + 1.D0/
     &    6912.D0*CF**2*omx**6 + 1.D0/6912.D0*CF**2*omx**7 + 1.D0/6912.D
     &    0*CF**2*omx**8 + 1.D0/6912.D0*CF**2*omx**9 + 1.D0/210.D0*
     &    CF**2*ln2 - 11.D0/630.D0*CF**2*ln2*omx - 2.D0/105.D0*CF**2*
     &    ln2*omx**2 + 64.D0/315.D0*CF**2*ln2*omx**3 - 62.D0/315.D0*
     &    CF**2*ln2*omx**4 )
      I2R(gq,rglr) = I2R(gq,rglr) + R4 * ( 23.D0/144.D0*CF*CA - 20963.D0
     &    /21600.D0*CF*CA*omx + 6899.D0/2880.D0*CF*CA*omx**2 - 4472597.D
     &    0/8467200.D0*CF*CA*omx**3 + 197879.D0/338688.D0*CF*CA*omx**4
     &     + 7224187.D0/8467200.D0*CF*CA*omx**5 + 55978019.D0/67737600.D
     &    0*CF*CA*omx**6 + 10987637953.D0/14753249280.D0*CF*CA*omx**7
     &     + 2030391919.D0/3073593600.D0*CF*CA*omx**8 + 11702351550749.D
     &    0/19946393026560.D0*CF*CA*omx**9 + 15987107.D0/2903040.D0*CF*
     &    CA*omx**10 - 518587.D0/362880.D0*CF*CA*omx**11 + 1.D0/2.D0*CF
     &    *CA*ln2 + 91.D0/90.D0*CF*CA*ln2*omx - 5.D0/4.D0*CF*CA*ln2*
     &    omx**2 + 412.D0/315.D0*CF*CA*ln2*omx**3 + 55.D0/252.D0*CF*CA*
     &    ln2*omx**4 + 2.D0/35.D0*CF*CA*ln2*omx**5 + 331.D0/2520.D0*CF*
     &    CA*ln2*omx**6 + 460.D0/2079.D0*CF*CA*ln2*omx**7 + 8221.D0/
     &    27720.D0*CF*CA*ln2*omx**8 + 9619.D0/27027.D0*CF*CA*ln2*omx**9
     &     - 22.D0/9.D0*CF*CA*ln2*omx**10 + 16.D0/9.D0*CF*CA*ln2*
     &    omx**11 + 26359.D0/57600.D0*CF**2 + 21413.D0/172800.D0*CF**2*
     &    omx )
      I2R(gq,rglr) = I2R(gq,rglr) + R4 * ( 4939.D0/28800.D0*CF**2*
     &    omx**2 + 5971.D0/3200.D0*CF**2*omx**3 + 815.D0/1152.D0*CF**2*
     &    omx**4 + 815.D0/1152.D0*CF**2*omx**5 + 815.D0/1152.D0*CF**2*
     &    omx**6 + 815.D0/1152.D0*CF**2*omx**7 + 815.D0/1152.D0*CF**2*
     &    omx**8 + 815.D0/1152.D0*CF**2*omx**9 - 1.D0/216.D0*CF**2*pisq
     &     - 1.D0/216.D0*CF**2*pisq*omx - 1.D0/108.D0*CF**2*pisq*omx**2
     &     - 1.D0/108.D0*CF**2*pisq*omx**3 - 1.D0/108.D0*CF**2*pisq*
     &    omx**4 - 1.D0/108.D0*CF**2*pisq*omx**5 - 1.D0/108.D0*CF**2*
     &    pisq*omx**6 - 1.D0/108.D0*CF**2*pisq*omx**7 - 1.D0/108.D0*
     &    CF**2*pisq*omx**8 - 1.D0/108.D0*CF**2*pisq*omx**9 - 49.D0/120.
     &    D0*CF**2*ln2 - 83.D0/360.D0*CF**2*ln2*omx - 29.D0/60.D0*CF**2
     &    *ln2*omx**2 - 109.D0/60.D0*CF**2*ln2*omx**3 - 3.D0/4.D0*CF**2
     &    *ln2*omx**4 - 3.D0/4.D0*CF**2*ln2*omx**5 - 3.D0/4.D0*CF**2*
     &    ln2*omx**6 - 3.D0/4.D0*CF**2*ln2*omx**7 - 3.D0/4.D0*CF**2*ln2
     &    *omx**8 - 3.D0/4.D0*CF**2*ln2*omx**9 )
      I2R(gq,rglr) = I2R(gq,rglr) + R2 * (  - 8*CF*CA + 12*CF*CA*omx - 
     &    47.D0/18.D0*CF*CA*omx**2 + 3453.D0/400.D0*CF*CA*omx**3 + 1409.
     &    D0/300.D0*CF*CA*omx**4 + 432679.D0/176400.D0*CF*CA*omx**5 + 
     &    355.D0/336.D0*CF*CA*omx**6 + 1321541.D0/10160640.D0*CF*CA*
     &    omx**7 - 13121809.D0/25401600.D0*CF*CA*omx**8 - 577838059.D0/
     &    585446400.D0*CF*CA*omx**9 + 194515.D0/290304.D0*CF*CA*omx**10
     &     + 511.D0/288.D0*CF*CA*omx**11 - 2.D0/3.D0*CF*CA*pisq - 2.D0/
     &    3.D0*CF*CA*pisq*omx - 4.D0/3.D0*CF*CA*pisq*omx**2 - 4.D0/3.D0
     &    *CF*CA*pisq*omx**3 - 4.D0/3.D0*CF*CA*pisq*omx**4 - 4.D0/3.D0*
     &    CF*CA*pisq*omx**5 - 4.D0/3.D0*CF*CA*pisq*omx**6 - 4.D0/3.D0*
     &    CF*CA*pisq*omx**7 - 4.D0/3.D0*CF*CA*pisq*omx**8 - 4.D0/3.D0*
     &    CF*CA*pisq*omx**9 + 8*CF*CA*ln2 - 12*CF*CA*ln2*omx + 10.D0/3.D
     &    0*CF*CA*ln2*omx**2 - 43.D0/5.D0*CF*CA*ln2*omx**3 - 22.D0/5.D0
     &    *CF*CA*ln2*omx**4 - 239.D0/105.D0*CF*CA*ln2*omx**5 - CF*CA*
     &    ln2*omx**6 - 19.D0/126.D0*CF*CA*ln2*omx**7 + 143.D0/315.D0*CF
     &    *CA*ln2*omx**8 )
      I2R(gq,rglr) = I2R(gq,rglr) + R2 * ( 299.D0/330.D0*CF*CA*ln2*
     &    omx**9 - 20.D0/9.D0*CF*CA*ln2*omx**10 + 101.D0/12.D0*CF**2 + 
     &    539.D0/36.D0*CF**2*omx + 193.D0/6.D0*CF**2*omx**2 + 145.D0/6.D
     &    0*CF**2*omx**3 + 145.D0/6.D0*CF**2*omx**4 + 145.D0/6.D0*CF**2
     &    *omx**5 + 145.D0/6.D0*CF**2*omx**6 + 145.D0/6.D0*CF**2*omx**7
     &     + 145.D0/6.D0*CF**2*omx**8 + 145.D0/6.D0*CF**2*omx**9 - 1.D0/
     &    3.D0*CF**2*pisq - 1.D0/3.D0*CF**2*pisq*omx - 2.D0/3.D0*CF**2*
     &    pisq*omx**2 - 2.D0/3.D0*CF**2*pisq*omx**3 - 2.D0/3.D0*CF**2*
     &    pisq*omx**4 - 2.D0/3.D0*CF**2*pisq*omx**5 - 2.D0/3.D0*CF**2*
     &    pisq*omx**6 - 2.D0/3.D0*CF**2*pisq*omx**7 - 2.D0/3.D0*CF**2*
     &    pisq*omx**8 - 2.D0/3.D0*CF**2*pisq*omx**9 - 7*CF**2*ln2 - 41.D
     &    0/3.D0*CF**2*ln2*omx - 30*CF**2*ln2*omx**2 - 22*CF**2*ln2*
     &    omx**3 - 22*CF**2*ln2*omx**4 - 22*CF**2*ln2*omx**5 - 22*CF**2
     &    *ln2*omx**6 - 22*CF**2*ln2*omx**7 - 22*CF**2*ln2*omx**8 - 22*
     &    CF**2*ln2*omx**9 )
      I2R(gq,rglr) = I2R(gq,rglr) + picube * ( 8*CF*CA*omx**(-1)*
     &    BCGQCFCA1 - 8*CF*CA*omx**(-1)*BCGQCFCAx + 8*CF*CA*BCGQCFCAx
     &     + 8*CF**2*omx**(-1)*BCGQCFCF1 - 8*CF**2*omx**(-1)*BCGQCFCFx
     &     + 8*CF**2*BCGQCFCFx )
      I2R(gq,rglr) = I2R(gq,rglr) + lnR * ( 24*CF**2 + 24*CF**2*omx + 
     &    48*CF**2*omx**2 + 48*CF**2*omx**3 + 48*CF**2*omx**4 + 48*
     &    CF**2*omx**5 + 48*CF**2*omx**6 + 48*CF**2*omx**7 + 48*CF**2*
     &    omx**8 + 48*CF**2*omx**9 - 8.D0/3.D0*CF**2*pisq - 8.D0/3.D0*
     &    CF**2*pisq*omx - 16.D0/3.D0*CF**2*pisq*omx**2 - 16.D0/3.D0*
     &    CF**2*pisq*omx**3 - 16.D0/3.D0*CF**2*pisq*omx**4 - 16.D0/3.D0
     &    *CF**2*pisq*omx**5 - 16.D0/3.D0*CF**2*pisq*omx**6 - 16.D0/3.D0
     &    *CF**2*pisq*omx**7 - 16.D0/3.D0*CF**2*pisq*omx**8 - 16.D0/3.D0
     &    *CF**2*pisq*omx**9 - 8.D0/3.D0*CF**2*ln512 - 8.D0/3.D0*CF**2*
     &    ln512*omx - 16.D0/3.D0*CF**2*ln512*omx**2 - 16.D0/3.D0*CF**2*
     &    ln512*omx**3 - 16.D0/3.D0*CF**2*ln512*omx**4 - 16.D0/3.D0*
     &    CF**2*ln512*omx**5 - 16.D0/3.D0*CF**2*ln512*omx**6 - 16.D0/3.D
     &    0*CF**2*ln512*omx**7 - 16.D0/3.D0*CF**2*ln512*omx**8 - 16.D0/
     &    3.D0*CF**2*ln512*omx**9 )
      I2R(gq,rglr) = I2R(gq,rglr) - 24*CF**2*lnpi - 24*CF**2*lnpi*omx
     &     - 48*CF**2*lnpi*omx**2 - 48*CF**2*lnpi*omx**3 - 48*CF**2*
     &    lnpi*omx**4 - 48*CF**2*lnpi*omx**5 - 48*CF**2*lnpi*omx**6 - 
     &    48*CF**2*lnpi*omx**7 - 48*CF**2*lnpi*omx**8 - 48*CF**2*lnpi*
     &    omx**9 + 8.D0/3.D0*CF**2*pisq*lnpi + 8.D0/3.D0*CF**2*pisq*
     &    lnpi*omx + 16.D0/3.D0*CF**2*pisq*lnpi*omx**2 + 16.D0/3.D0*
     &    CF**2*pisq*lnpi*omx**3 + 16.D0/3.D0*CF**2*pisq*lnpi*omx**4 + 
     &    16.D0/3.D0*CF**2*pisq*lnpi*omx**5 + 16.D0/3.D0*CF**2*pisq*
     &    lnpi*omx**6 + 16.D0/3.D0*CF**2*pisq*lnpi*omx**7 + 16.D0/3.D0*
     &    CF**2*pisq*lnpi*omx**8 + 16.D0/3.D0*CF**2*pisq*lnpi*omx**9 - 
     &    24*CF**2*ln2 - 24*CF**2*ln2*omx - 48*CF**2*ln2*omx**2 - 48*
     &    CF**2*ln2*omx**3 - 48*CF**2*ln2*omx**4 - 48*CF**2*ln2*omx**5
     &     - 48*CF**2*ln2*omx**6 - 48*CF**2*ln2*omx**7 - 48*CF**2*ln2*
     &    omx**8 - 48*CF**2*ln2*omx**9 + 24*CF**2*ln2*lnpi + 24*CF**2*
     &    ln2*lnpi*omx + 48*CF**2*ln2*lnpi*omx**2 + 48*CF**2*ln2*lnpi*
     &    omx**3
      I2R(gq,rglr) = I2R(gq,rglr) + 48*CF**2*ln2*lnpi*omx**4 + 48*CF**2
     &    *ln2*lnpi*omx**5 + 48*CF**2*ln2*lnpi*omx**6 + 48*CF**2*ln2*
     &    lnpi*omx**7 + 48*CF**2*ln2*lnpi*omx**8 + 48*CF**2*ln2*lnpi*
     &    omx**9 + 8.D0/3.D0*CF**2*ln2*pisq + 8.D0/3.D0*CF**2*ln2*pisq*
     &    omx + 16.D0/3.D0*CF**2*ln2*pisq*omx**2 + 16.D0/3.D0*CF**2*ln2
     &    *pisq*omx**3 + 16.D0/3.D0*CF**2*ln2*pisq*omx**4 + 16.D0/3.D0*
     &    CF**2*ln2*pisq*omx**5 + 16.D0/3.D0*CF**2*ln2*pisq*omx**6 + 16.
     &    D0/3.D0*CF**2*ln2*pisq*omx**7 + 16.D0/3.D0*CF**2*ln2*pisq*
     &    omx**8 + 16.D0/3.D0*CF**2*ln2*pisq*omx**9 + 24*CF**2*ln2**2
     &     + 24*CF**2*ln2**2*omx + 48*CF**2*ln2**2*omx**2 + 48*CF**2*
     &    ln2**2*omx**3 + 48*CF**2*ln2**2*omx**4 + 48*CF**2*ln2**2*
     &    omx**5 + 48*CF**2*ln2**2*omx**6 + 48*CF**2*ln2**2*omx**7 + 48
     &    *CF**2*ln2**2*omx**8 + 48*CF**2*ln2**2*omx**9
      endif
