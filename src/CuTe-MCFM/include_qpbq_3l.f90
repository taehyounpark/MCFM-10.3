!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
! The 3-loop beam function is taken from the paper:
! 'Unpolarized Quark and Gluon TMD PDFs and FFs at N3LO'
! Ming-xing Luo, Tong-Zhi Yang, Hua Xing Zhu, Yu Jiao Zhu
! https://arxiv.org/abs/2012.03256
!
! This include file is obtained from the ancillary file
! TMDPDFN.m, and is obtained from the expression TMDpdfN["qq"]
! therein after multiplication by the factor exp(-h_i),
! where h_i is given by the solution of Eq. (5) in
! arXiv:2009.11437, and by identifying Lh = -Lp.
! This removes the leading logarithmic behavior, passing
! from the functions "B" to "\bar B" (c.f. Eq. (4))
!
! The labelling of 'type_in' is as follows:
!       1       regular term
!       2       coefficient of delta(1-z)
!       3       coefficient of (1/(1-z))_+
!       4       coefficient of (log(1-z)/(1-z))_+
!       5       coefficient of (log^2(1-z)/(1-z))_+
select case (powAs)
    case (0)
        select case (powLperp)
            case (0)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (1)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (2)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (3)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case default
                Ibar_select = 0._dp
        end select
    case (1)
        select case (powLperp)
            case (0)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (1)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (2)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (3)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case default
                Ibar_select = 0._dp
        end select
    case (2)
        select case (powLperp)
            case (0)
                select case (type_in)
                    case (1)
                        Ibar_select = z**3*(-3.6040303008016474_dp) + (-3.2270183101427925_dp) + z**5*(-0.0013164390324173112_dp) + &
                        z**4*(0.05637832379670145_dp) + z*(0.5603571563984158_dp) + (2.6451727005890224_dp)/z + &
                        z**2*(3.5704568978264084_dp) + (z*(-5.333333333333333_dp) + z**3*(1.215623396613648_dp) + &
                        (2.6666666666666665_dp) + z**2*(5.654771516227999_dp))*Log(z) + (z**2*(-1.773058773058773_dp) + &
                        (-0.6666666666666666_dp) + z*(-0.6666666666666666_dp) + z**3*(-0.14346862340441793_dp))*Log(z)**2 + &
                        (z**3*(-0.0170296134428453_dp) + z**2*(0.00030671710459053266_dp) + (0.4444444444444444_dp) + &
                        z*(0.4444444444444444_dp))*Log(z)**3

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (1)
                select case (type_in)
                    case (1)
                        Ibar_select = z*(-16._dp) + (-7.703703703703703_dp)/z + (10.666666666666666_dp) + &
                        z**2*(13.037037037037036_dp) + (z*(-8._dp) + z**2*(-7.111111111111111_dp) + (-2.6666666666666665_dp))*Log(z) &
                        + ((2.6666666666666665_dp) + z*(2.6666666666666665_dp))*Log(z)**2

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (2)
                select case (type_in)
                    case (1)
                        Ibar_select = z**2*(-1.7777777777777777_dp) + z*(-1.3333333333333333_dp) + (1.3333333333333333_dp) + &
                        (1.7777777777777777_dp)/z + ((2.6666666666666665_dp) + z*(2.6666666666666665_dp))*Log(z)

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (3)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case default
                Ibar_select = 0._dp
        end select
    case (3)
        select case (powLperp)
            case (0)
                select case (type_in)
                    case (1)
                        Ibar_select = (z**3*(-3.3748641023051906e8_dp) + z**6*(-39310.43325413669_dp) + &
                        z**2*(-3397.2559799966602_dp) + (-472.18580115215207_dp) + Nf*(z*(-59.09868138406361_dp) + &
                        z**3*(-24.98539332433859_dp) + z**5*(-0.9090821349290608_dp) + z**6*(0.10213731794968792_dp) + &
                        (5.88282387997262_dp) + z**4*(18.534746870913505_dp) + z**2*(60.47344877255524_dp)) + &
                        z*(494.1845796659956_dp) + z**5*(3.2415703570137727e6_dp) + z**4*(3.3428752556419367e8_dp))/z + &
                        (z**3*(-2.0326153123472175e8_dp) + z**2*(-1.3739440190163258e8_dp) + (-420.5501349733193_dp) + &
                        z*(-257.507993015641_dp) + (-78.98467140066984_dp)/z + Nf*((-30.966480714928217_dp) + &
                        z*(-12.003517751965253_dp) + z**3*(-10.540989729225023_dp) + z**2*(5.1478051820853_dp)))*Log(z) + &
                        (z**2*(-2.3899558551923998e7_dp) + z*(-52.57566145124676_dp) + Nf*((-10.549308734586743_dp) + &
                        z*(-8.771530956808965_dp) + z**2*(0.1997362426567558_dp) + z**3*(3.5376079174563064_dp)) + &
                        (89.39183109189284_dp) + z**3*(5.3702064049717955e7_dp))*Log(z)**2 + (z**3*(-9.716068704083942e6_dp) + &
                        z**2*(-2.225269152974023e6_dp) + (-46.79702185625944_dp) + z*(-3.7846761772470945_dp) + &
                        Nf*((-1.6131687242798354_dp) + z*(-1.6131687242798354_dp) + z**3*(-0.3148828372387587_dp) + &
                        z**2*(0.7031832520791511_dp)))*Log(z)**3 + (z**2*(-111179.04088538211_dp) + Nf*((-0.24691358024691357_dp) + &
                        z*(-0.24691358024691357_dp) + z**2*(-0.002601203056413591_dp) + z**3*(0.11712315648343229_dp)) + &
                        z*(2.345679012345679_dp) + (6.419753086419753_dp) + z**3*(887693.1712435037_dp))*Log(z)**4 + &
                        (z**3*(-85139.58559616876_dp) + z**2*(-2394.7528571136036_dp) + (-0.4444444444444444_dp) + &
                        z*(1.2740740740740741_dp))*Log(z)**5 + ((-2792.4680875169915_dp) + z**2*(-1078.8475886795727_dp) + &
                        z**3*(-875.2059658293723_dp) + Nf*(z**2*(-7.70122877162304_dp) + (-5.544759836563321_dp) + &
                        z**3*(3.185230352303523_dp) + z*(10.060758255882838_dp)) + z*(4746.521642025937_dp))*Log(z*(-1._dp) + &
                        (1._dp)) + (z*(-12517.727530004448_dp) + z**3*(-4430.35070098663_dp) + Nf*(z*(-1.1444946128399929_dp) + &
                        z**3*(0.030580661136560605_dp) + (0.3406240567413631_dp) + z**2*(0.773289894962069_dp)) + &
                        (4046.7576591163565_dp) + z**2*(12901.320571874721_dp))*Log(z*(-1._dp) + (1._dp))**2 + &
                        (z*(-42.897419068067926_dp) + z**3*(-13.159651273668437_dp) + Nf*(z**2*(-0.6827091430256285_dp) + &
                        (-0.32694334958049026_dp) + z**3*(0.22726561277195179_dp) + z*(0.7823868798341669_dp)) + &
                        (15.140488835471352_dp) + z**2*(40.91658150626502_dp))*Log(z*(-1._dp) + (1._dp))**3 + &
                        (z**2*(-7.504455857694885_dp) + (-2.2065175424656394_dp) + z**3*(2.525512367491166_dp) + &
                        z*(7.185461032669358_dp))*Log(z*(-1._dp) + (1._dp))**4

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (1)
                select case (type_in)
                    case (1)
                        Ibar_select = (z**3*(-4.1134098173339283e6_dp) + z**5*(-136507.60461517842_dp) + z*(-471.45491946456224_dp) &
                        + z**2*(-134.5277776315147_dp) + Nf*(z*(-47.84250324323074_dp) + z**4*(-4.589616073204694_dp) + &
                        (-1.1565265637483262_dp) + z**5*(-1.0078269216569447_dp) + z**3*(0.17401837928153718_dp) + &
                        z**6*(0.2108997733320193_dp) + z**2*(54.211554625724375_dp)) + (483.9673467888446_dp) + &
                        z**6*(3110.7370438220105_dp) + z**4*(4.2469287000109935e6_dp))/z + (z**3*(-2.338617730777145e6_dp) + &
                        z**2*(-1.6447847163410652e6_dp) + z*(-99.2512121237033_dp) + Nf*((-29.870148425982446_dp) + &
                        z*(-16.833111388945408_dp) + z**3*(4.0626629422718805_dp) + z**2*(12.756509161041466_dp)) + &
                        (98.11244277743864_dp)/z + (449.6164685898825_dp))*Log(z) + (z**2*(-267474.5257407032_dp) + &
                        (-158.3322686966518_dp) + Nf*((-7.703703703703703_dp) + z*(-7.703703703703703_dp) + &
                        z**3*(-0.415594042725251_dp) + z**2*(2.7071215217733062_dp)) + z*(7.941645856015171_dp) + &
                        z**3*(741485.6163247644_dp))*Log(z)**2 + (z**3*(-90685.81400758731_dp) + z**2*(-21036.33587666063_dp) + &
                        z*(-5.62962962962963_dp) + Nf*((-1.1851851851851851_dp) + z*(-1.1851851851851851_dp) + &
                        z**2*(0.0169144727867958_dp) + z**3*(0.3615537167475894_dp)) + (31.111111111111114_dp))*Log(z)**3 + &
                        (z**2*(-678.7419680866079_dp) + (-2.9629629629629632_dp) + z*(8.296296296296296_dp) + &
                        z**3*(14371.52674381814_dp))*Log(z)**4 + (z**2*(-1419.9417460688342_dp) + (-537.5383226650899_dp) + &
                        Nf*(z*(-4.559282020884582_dp) + z**3*(-0.17257441490570324_dp) + (0.7118723215079577_dp) + &
                        z**2*(4.019984114282327_dp)) + z**3*(428.56402183496346_dp) + z*(1528.9160468989608_dp))*Log(z*(-1._dp) + &
                        (1._dp)) + (z*(-74.74156853609831_dp) + z**3*(-21.393301000303122_dp) + Nf*(z**2*(-3.202431556813283_dp) + &
                        (-1.9599822537344245_dp) + z**3*(1.0656690959838735_dp) + z*(4.096744714563834_dp)) + &
                        (31.340800434564265_dp) + z**2*(64.79406910183718_dp))*Log(z*(-1._dp) + (1._dp))**2 + &
                        (z**2*(-10.102286050707095_dp) + (-1.9908893803281575_dp) + z**3*(3.3149575944487277_dp) + &
                        z*(8.778217836586524_dp))*Log(z*(-1._dp) + (1._dp))**3

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (2)
                select case (type_in)
                    case (1)
                        Ibar_select = (z**3*(-1879.0460834615833_dp) + (-220.95532436198914_dp) + z**2*(-158.81336518685896_dp) + &
                        z**6*(-20.021939275220372_dp) + Nf*(z*(-15.922728604826144_dp) + z**5*(-0.7777828865694975_dp) + &
                        z**4*(-0.47407834101382496_dp) + z**6*(0.1798342541436464_dp) + (1.1851851851851851_dp) + &
                        z**3*(5.623854999025531_dp) + z**2*(10.185715290166272_dp)) + z*(263.26740393874053_dp) + &
                        z**5*(307.82203849281643_dp) + z**4*(1707.7472716829411_dp))/z + (z**3*(-1482.192546009433_dp) + &
                        z**2*(-807.3621590952911_dp) + (-145.55928700679456_dp) + z*(-127.39461783961922_dp) + &
                        (-46.222222222222236_dp)/z + Nf*((-7.703703703703703_dp) + z*(-7.703703703703703_dp) + &
                        z**3*(-0.15299579924828655_dp) + z**2*(1.909043965771614_dp)))*Log(z) + (z**2*(-127.36102354904506_dp) + &
                        z*(-14.222222222222221_dp) + Nf*((-1.7777777777777777_dp) + z*(-1.7777777777777777_dp) + &
                        z**2*(-0.04630079326533915_dp) + z**3*(0.5050636719141682_dp)) + (38.222222222222214_dp) + &
                        z**3*(214.53622606133652_dp))*Log(z)**2 + (z**3*(-84.36428149024246_dp) + (-5.925925925925927_dp) + &
                        z**2*(-5.406404346394467_dp) + z*(15.407407407407407_dp))*Log(z)**3 + (z*(-177.77963268425535_dp) + &
                        z**3*(-30.466535499847353_dp) + Nf*(z**2*(-3.511929408639212_dp) + (-2.929127721470389_dp) + &
                        z**3*(1.1802897324733002_dp) + z*(5.260767397636301_dp)) + (87.4343263699818_dp) + &
                        z**2*(120.81184181412094_dp))*Log(z*(-1._dp) + (1._dp))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (3)
                select case (type_in)
                    case (1)
                        Ibar_select = (-50.62385710459351_dp) + z**4*(-4.407436322826193_dp) + z**3*(-2.6864439571778664_dp) + &
                        z**2*(-0.5265500760136049_dp) + Nf*((-0.7901234567901234_dp)/z + (-0.5925925925925926_dp) + &
                        z*(0.5925925925925926_dp) + z**2*(0.7901234567901234_dp)) + z**5*(1.019060773480663_dp) + &
                        z*(18.114115066319957_dp) + (39.11111111111111_dp)/z + (z**2*(-8.935170641153789_dp) + &
                        Nf*((-1.1851851851851851_dp) + z*(-1.1851851851851851_dp)) + z**3*(-0.8669762220634759_dp) + &
                        (7.111111111111111_dp)/z + z*(9.481481481481485_dp) + (19.555555555555557_dp))*Log(z) + &
                        ((-2.9629629629629632_dp) + z**2*(-0.2623712088231026_dp) + z**3*(2.862027514358221_dp) + &
                        z*(13.037037037037036_dp))*Log(z)**2 + (z**2*(-19.900933497371888_dp) + (-16.598390495119958_dp) + &
                        z**3*(6.688308538163001_dp) + z*(29.811015454328846_dp))*Log(z*(-1._dp) + (1._dp))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case default
                Ibar_select = 0._dp
        end select
    case (4)
        select case (powLperp)
            case (2)
                select case (type_in)
                    case (1)
                        Ibar_select = z**2*(-120491.36783272114_dp) + (-54978.43935530921_dp) + z*(-19520.45146856377_dp) + &
                        (nf*(-386.14301560237607_dp))/z + nf**2*z**2*(-10.534979423868313_dp) + (nf**2*(0.3664031069582032_dp))/z + &
                        nf**2*z*(1.597477690395161_dp) + nf*(38.20799772468365_dp) + nf**2*(51.461231332483514_dp) + &
                        nf*z**2*(1197.1006285558358_dp) + nf*z*(1316.7400472679383_dp) + (22375.471060762604_dp)/z + &
                        ((-2837.889696814507_dp)/z + z**2*(-484.7810524461406_dp) + nf*((-23.39461783961922_dp)/z + &
                        z**2*(-23.39461783961922_dp) + (263.18945069571623_dp) + z*(263.18945069571623_dp)) + &
                        (10791.001219116066_dp) + z*(13144.109863484431_dp))*Hr1(-1) + z*(-18313.654706113237_dp)*Hr1(0) + &
                        nf*(-2162.287418183806_dp)*Hr1(0) + nf*z**2*(-1159.7330744041965_dp)*Hr1(0) + &
                        nf*z*(-107.9955611740408_dp)*Hr1(0) + (nf*(-19.34231284582478_dp)*Hr1(0))/z + &
                        nf**2*z**2*(2.3703703703703702_dp)*Hr1(0) + nf**2*z*(26.314592870426893_dp)*Hr1(0) + &
                        nf**2*(39.35162990746393_dp)*Hr1(0) + ((5458.699640174869_dp)*Hr1(0))/z + (68218.40561505272_dp)*Hr1(0) + &
                        z**2*(89022.34445378877_dp)*Hr1(0) + z**2*(-12247.672701091973_dp)*Hr1(1) + &
                        nf*z*(-704.936318030071_dp)*Hr1(1) + (nf*(-169.28655538377376_dp)*Hr1(1))/z + &
                        nf**2*(-24.88888888888889_dp)*Hr1(1) + nf**2*z**2*(-2.3703703703703702_dp)*Hr1(1) + &
                        (nf**2*(2.3703703703703702_dp)*Hr1(1))/z + nf**2*z*(24.88888888888889_dp)*Hr1(1) + &
                        nf*z**2*(234.2742097047614_dp)*Hr1(1) + nf*(639.9486637090833_dp)*Hr1(1) + &
                        ((2588.7836896544595_dp)*Hr1(1))/z + z*(2852.977187539238_dp)*Hr1(1) + (6805.911823898275_dp)*Hr1(1) + &
                        (-1957.3496925814748_dp)*Hr2(-1,-1) + z*(-1957.3496925814748_dp)*Hr2(-1,-1) + &
                        ((4356.597722133534_dp)*Hr2(-1,-1))/z + z**2*(4356.597722133534_dp)*Hr2(-1,-1) + &
                        z**2*(-22232.42778950975_dp)*Hr2(-1,0) + z*(-17393.18096586072_dp)*Hr2(-1,0) + &
                        ((-16578.10680185543_dp)*Hr2(-1,0))/z + (-11738.859978206401_dp)*Hr2(-1,0) + &
                        (nf*(79.53909465020575_dp)*Hr2(-1,0))/z + nf*z**2*(100.8724279835391_dp)*Hr2(-1,0) + &
                        nf*(289.18518518518516_dp)*Hr2(-1,0) + nf*z*(310.51851851851853_dp)*Hr2(-1,0) + &
                        z**2*(-7403.096845248393_dp)*Hr2(0,-1) + ((-2105.51560556573_dp)*Hr2(0,-1))/z + &
                        nf*z*(-198.85425163676337_dp)*Hr2(0,-1) + nf*(198.85425163676337_dp)*Hr2(0,-1) + &
                        z*(769.2745768945342_dp)*Hr2(0,-1) + (5052.086162197381_dp)*Hr2(0,-1) + (-12332.391147006294_dp)*Hr2(0,0) + &
                        z**2*(-4288.541551294237_dp)*Hr2(0,0) + nf*z*(-422.3378623916079_dp)*Hr2(0,0) + &
                        nf*(-422.29380634086874_dp)*Hr2(0,0) + nf**2*z**2*(-4.7407407407407405_dp)*Hr2(0,0) + &
                        nf**2*(15.407407407407407_dp)*Hr2(0,0) + nf**2*z*(15.407407407407407_dp)*Hr2(0,0) + &
                        ((588.6746566646319_dp)*Hr2(0,0))/z + nf*z**2*(785.2510288065844_dp)*Hr2(0,0) + &
                        z*(4551.873018643857_dp)*Hr2(0,0) + ((-2153.499958396644_dp)*Hr2(0,1))/z + &
                        z*(-1727.7841905755781_dp)*Hr2(0,1) + nf*z**2*(-218.8641975308642_dp)*Hr2(0,1) + &
                        (nf*(-16.59259259259259_dp)*Hr2(0,1))/z + nf**2*(-15.407407407407407_dp)*Hr2(0,1) + &
                        nf**2*z*(-15.407407407407407_dp)*Hr2(0,1) + nf**2*z**2*(4.7407407407407405_dp)*Hr2(0,1) + &
                        nf*(42.563994084833155_dp)*Hr2(0,1) + nf*z*(452.4405372947097_dp)*Hr2(0,1) + (14076.44507375338_dp)*Hr2(0,1) &
                        + z**2*(17159.872264625257_dp)*Hr2(0,1) + z*(-17376.03077628846_dp)*Hr2(1,0) + &
                        ((-13913.798714631412_dp)*Hr2(1,0))/z + nf*(-529.3827160493827_dp)*Hr2(1,0) + &
                        nf*z**2*(-312.2304526748971_dp)*Hr2(1,0) + nf**2*z*(-7.111111111111111_dp)*Hr2(1,0) + &
                        (nf**2*(-4.7407407407407405_dp)*Hr2(1,0))/z + nf**2*z**2*(4.7407407407407405_dp)*Hr2(1,0) + &
                        nf**2*(7.111111111111111_dp)*Hr2(1,0) + (nf*(324.082304526749_dp)*Hr2(1,0))/z + &
                        nf*z*(517.5308641975308_dp)*Hr2(1,0) + (13278.25299851068_dp)*Hr2(1,0) + &
                        z**2*(18011.57649240919_dp)*Hr2(1,0) + ((-4.7407407407407405_dp)*(z**2*(-640.4882719846497_dp) + &
                        (-9.191563279726404_dp) + nf**2*((-1._dp) + z*(-0.75_dp) + z**2*(0.75_dp) + z**3*(1._dp)) + &
                        nf*(z**3*(-58.166666666666664_dp) + z*(-47.083333333333336_dp) + (42.666666666666664_dp) + &
                        z**2*(62.583333333333336_dp)) + z*(265.67577198464966_dp) + z**3*(384.00406327972644_dp))*Hr2(1,1))/z + &
                        (nf*(-47.407407407407405_dp)*Hr3(-1,-1,0))/z + nf*z**2*(-47.407407407407405_dp)*Hr3(-1,-1,0) + &
                        nf*(35.55555555555556_dp)*Hr3(-1,-1,0) + nf*z*(35.55555555555556_dp)*Hr3(-1,-1,0) + &
                        ((2189.9588477366256_dp)*Hr3(-1,-1,0))/z + z**2*(3386.9958847736625_dp)*Hr3(-1,-1,0) + &
                        (5835.061728395061_dp)*Hr3(-1,-1,0) + z*(7032.098765432099_dp)*Hr3(-1,-1,0) + &
                        z*(-792.6913580246913_dp)*Hr3(-1,0,0) + (nf*(-94.81481481481481_dp)*Hr3(-1,0,0))/z + &
                        nf*z**2*(-94.81481481481481_dp)*Hr3(-1,0,0) + nf*(-3.5555555555555554_dp)*Hr3(-1,0,0) + &
                        nf*z*(-3.5555555555555554_dp)*Hr3(-1,0,0) + (221.82716049382717_dp)*Hr3(-1,0,0) + &
                        z**2*(4059.3909465020574_dp)*Hr3(-1,0,0) + ((5073.909465020576_dp)*Hr3(-1,0,0))/z + &
                        z*(-4211.753086419753_dp)*Hr3(-1,0,1) + (-3379.753086419753_dp)*Hr3(-1,0,1) + &
                        nf*(-142.22222222222223_dp)*Hr3(-1,0,1) + nf*z*(-142.22222222222223_dp)*Hr3(-1,0,1) + &
                        (nf*(-9.481481481481481_dp)*Hr3(-1,0,1))/z + nf*z**2*(-9.481481481481481_dp)*Hr3(-1,0,1) + &
                        z**2*(279.7037037037037_dp)*Hr3(-1,0,1) + ((1111.7037037037037_dp)*Hr3(-1,0,1))/z + &
                        (-5661.497517187851_dp)*Hr3(0,-1,-1) + z*(5661.497517187851_dp)*Hr3(0,-1,-1) + &
                        z*(-16426.4937025299_dp)*Hr3(0,-1,0) + ((-2090.6666666666665_dp)*Hr3(0,-1,0))/z + &
                        nf*z**2*(-121.67901234567901_dp)*Hr3(0,-1,0) + nf*z*(-67.55555555555556_dp)*Hr3(0,-1,0) + &
                        (nf*(-1.5802469135802468_dp)*Hr3(0,-1,0))/z + nf*(138.66666666666666_dp)*Hr3(0,-1,0) + &
                        z**2*(8073.744855967078_dp)*Hr3(0,-1,0) + (13269.950492653357_dp)*Hr3(0,-1,0) + &
                        z*(-5501.634295283787_dp)*Hr3(0,0,-1) + (10909.690119209095_dp)*Hr3(0,0,-1) + &
                        z**2*(-10975.604938271605_dp)*Hr3(0,0,0) + nf*(-855.5061728395061_dp)*Hr3(0,0,0) + &
                        nf*z*(-405.1358024691358_dp)*Hr3(0,0,0) + nf**2*(7.111111111111111_dp)*Hr3(0,0,0) + &
                        nf**2*z*(7.111111111111111_dp)*Hr3(0,0,0) + nf*z**2*(83.75308641975309_dp)*Hr3(0,0,0) + &
                        z*(3979.456209638797_dp)*Hr3(0,0,0) + (7176.850004696378_dp)*Hr3(0,0,0) + &
                        z**2*(-8290.502057613168_dp)*Hr3(0,0,1) + ((-1109.3333333333333_dp)*Hr3(0,0,1))/z + &
                        (-732.1737179986515_dp)*Hr3(0,0,1) + nf*(-95.20987654320987_dp)*Hr3(0,0,1) + &
                        nf*z*(-25.28395061728395_dp)*Hr3(0,0,1) + nf**2*(-7.111111111111111_dp)*Hr3(0,0,1) + &
                        nf**2*z*(-7.111111111111111_dp)*Hr3(0,0,1) + nf*z**2*(61.629629629629626_dp)*Hr3(0,0,1) + &
                        z*(508.8991722523715_dp)*Hr3(0,0,1) + z*(-7808.208924129093_dp)*Hr3(0,1,0) + &
                        (-6198.134850055019_dp)*Hr3(0,1,0) + ((-3320.8888888888887_dp)*Hr3(0,1,0))/z + &
                        z**2*(-2144.7901234567903_dp)*Hr3(0,1,0) + nf*z*(-57.28395061728395_dp)*Hr3(0,1,0) + &
                        nf*(-39.50617283950617_dp)*Hr3(0,1,0) + nf*z**2*(-26.864197530864196_dp)*Hr3(0,1,0) + &
                        (nf*(83.75308641975309_dp)*Hr3(0,1,0))/z + ((-2042.0740740740741_dp)*Hr3(0,1,1))/z + &
                        z*(-1500.530776348283_dp)*Hr3(0,1,1) + nf*z**2*(-124.8395061728395_dp)*Hr3(0,1,1) + &
                        nf*z*(-118.91358024691358_dp)*Hr3(0,1,1) + nf*(-13.432098765432098_dp)*Hr3(0,1,1) + &
                        nf**2*(7.111111111111111_dp)*Hr3(0,1,1) + nf**2*z*(7.111111111111111_dp)*Hr3(0,1,1) + &
                        (nf*(56.888888888888886_dp)*Hr3(0,1,1))/z + (165.4445322936923_dp)*Hr3(0,1,1) + &
                        z**2*(3033.0205761316874_dp)*Hr3(0,1,1) + (-2051.753086419753_dp)*Hr3(1,0,0) + &
                        z**2*(-449.3168724279835_dp)*Hr3(1,0,0) + nf*(-114.96296296296296_dp)*Hr3(1,0,0) + &
                        (nf*(-44.24691358024691_dp)*Hr3(1,0,0))/z + nf*z**2*(44.24691358024691_dp)*Hr3(1,0,0) + &
                        nf*z*(114.96296296296296_dp)*Hr3(1,0,0) + ((1250.5020576131687_dp)*Hr3(1,0,0))/z + &
                        z*(1250.567901234568_dp)*Hr3(1,0,0) + ((-3921.514403292181_dp)*Hr3(1,0,1))/z + &
                        z*(-2910.222222222222_dp)*Hr3(1,0,1) + nf*z**2*(-93.23456790123457_dp)*Hr3(1,0,1) + &
                        nf*(-1.1851851851851851_dp)*Hr3(1,0,1) + nf*z*(1.1851851851851851_dp)*Hr3(1,0,1) + &
                        (nf*(93.23456790123457_dp)*Hr3(1,0,1))/z + (2080.5925925925926_dp)*Hr3(1,0,1) + &
                        z**2*(4751.14403292181_dp)*Hr3(1,0,1) + ((-4703.736625514403_dp)*Hr3(1,1,0))/z + &
                        z*(-3619.358024691358_dp)*Hr3(1,1,0) + nf*z**2*(-109.03703703703704_dp)*Hr3(1,1,0) + &
                        nf*z*(-10.666666666666666_dp)*Hr3(1,1,0) + nf*(10.666666666666666_dp)*Hr3(1,1,0) + &
                        (nf*(109.03703703703704_dp)*Hr3(1,1,0))/z + (2789.7283950617284_dp)*Hr3(1,1,0) + &
                        z**2*(5533.366255144033_dp)*Hr3(1,1,0) + z*(-783.604938271605_dp)*Hr3(1,1,1) + &
                        z**2*(-247.7037037037037_dp)*Hr3(1,1,1) + ((-214.5185185185185_dp)*Hr3(1,1,1))/z + &
                        (nf*(-61.629629629629626_dp)*Hr3(1,1,1))/z + nf*(-46.22222222222222_dp)*Hr3(1,1,1) + &
                        nf*z*(46.22222222222222_dp)*Hr3(1,1,1) + nf*z**2*(61.629629629629626_dp)*Hr3(1,1,1) + &
                        (1245.8271604938273_dp)*Hr3(1,1,1) + (-393.48148148148147_dp)*Hr4(-1,-1,-1,0) + &
                        z*(-393.48148148148147_dp)*Hr4(-1,-1,-1,0) + ((524.641975308642_dp)*Hr4(-1,-1,-1,0))/z + &
                        z**2*(524.641975308642_dp)*Hr4(-1,-1,-1,0) + ((-3136.7901234567903_dp)*Hr4(-1,-1,0,0))/z + &
                        z**2*(-3136.7901234567903_dp)*Hr4(-1,-1,0,0) + (1954.3703703703704_dp)*Hr4(-1,-1,0,0) + &
                        z*(1954.3703703703704_dp)*Hr4(-1,-1,0,0) + ((-2386.1728395061727_dp)*Hr4(-1,-1,0,1))/z + &
                        z**2*(-2386.1728395061727_dp)*Hr4(-1,-1,0,1) + (993.1851851851852_dp)*Hr4(-1,-1,0,1) + &
                        z*(993.1851851851852_dp)*Hr4(-1,-1,0,1) + ((-1238.9135802469136_dp)*Hr4(-1,0,-1,0))/z + &
                        z**2*(-1238.9135802469136_dp)*Hr4(-1,0,-1,0) + (929.1851851851852_dp)*Hr4(-1,0,-1,0) + &
                        z*(929.1851851851852_dp)*Hr4(-1,0,-1,0) + (-1460.148148148148_dp)*Hr4(-1,0,0,0) + &
                        z*(-1460.148148148148_dp)*Hr4(-1,0,0,0) + ((2477.8271604938273_dp)*Hr4(-1,0,0,0))/z + &
                        z**2*(2477.8271604938273_dp)*Hr4(-1,0,0,0) + (-861.6296296296297_dp)*Hr4(-1,0,0,1) + &
                        z*(-861.6296296296297_dp)*Hr4(-1,0,0,1) + ((2741.7283950617284_dp)*Hr4(-1,0,0,1))/z + &
                        z**2*(2741.7283950617284_dp)*Hr4(-1,0,0,1) + (-18.962962962962962_dp)*Hr4(-1,0,1,0) + &
                        z*(-18.962962962962962_dp)*Hr4(-1,0,1,0) + ((1087.20987654321_dp)*Hr4(-1,0,1,0))/z + &
                        z**2*(1087.20987654321_dp)*Hr4(-1,0,1,0) + ((707.9506172839506_dp)*Hr4(-1,0,1,1))/z + &
                        z**2*(707.9506172839506_dp)*Hr4(-1,0,1,1) + (1061.9259259259259_dp)*Hr4(-1,0,1,1) + &
                        z*(1061.9259259259259_dp)*Hr4(-1,0,1,1) + z**2*(-2408.296296296296_dp)*Hr4(0,-1,-1,0) + &
                        ((-512._dp)*Hr4(0,-1,-1,0))/z + nf*z*(-71.11111111111111_dp)*Hr4(0,-1,-1,0) + &
                        nf*(71.11111111111111_dp)*Hr4(0,-1,-1,0) + (440.8888888888889_dp)*Hr4(0,-1,-1,0) + &
                        z*(1749.3333333333333_dp)*Hr4(0,-1,-1,0) + (-3491.5555555555557_dp)*Hr4(0,-1,0,0) + &
                        nf*z*(-92.44444444444444_dp)*Hr4(0,-1,0,0) + z*(-86.51851851851852_dp)*Hr4(0,-1,0,0) + &
                        nf*(92.44444444444444_dp)*Hr4(0,-1,0,0) + ((1792._dp)*Hr4(0,-1,0,0))/z + &
                        z**2*(3509.7283950617284_dp)*Hr4(0,-1,0,0) + z*(-1476.7407407407406_dp)*Hr4(0,-1,0,1) + &
                        (-967.1111111111111_dp)*Hr4(0,-1,0,1) + nf*(-85.33333333333333_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(85.33333333333333_dp)*Hr4(0,-1,0,1) + ((1024._dp)*Hr4(0,-1,0,1))/z + &
                        z**2*(3296.395061728395_dp)*Hr4(0,-1,0,1) + (-4618.666666666667_dp)*Hr4(0,0,-1,0) + &
                        z*(-378.0740740740741_dp)*Hr4(0,0,-1,0) + nf*z*(-142.22222222222223_dp)*Hr4(0,0,-1,0) + &
                        nf*(213.33333333333334_dp)*Hr4(0,0,-1,0) + z**2*(821.7283950617284_dp)*Hr4(0,0,-1,0) + &
                        (-1611.2592592592594_dp)*Hr4(0,0,0,0) + z**2*(-1112.4938271604938_dp)*Hr4(0,0,0,0) + &
                        nf*z*(-219.25925925925927_dp)*Hr4(0,0,0,0) + nf*(-84.14814814814815_dp)*Hr4(0,0,0,0) + &
                        z*(1557.9259259259259_dp)*Hr4(0,0,0,0) + z**2*(-2607.4074074074074_dp)*Hr4(0,0,0,1) + &
                        z*(-2105.4814814814813_dp)*Hr4(0,0,0,1) + nf*(-218.07407407407408_dp)*Hr4(0,0,0,1) + &
                        nf*z*(-47.407407407407405_dp)*Hr4(0,0,0,1) + (5313.777777777777_dp)*Hr4(0,0,0,1) + &
                        z*(-5007.407407407408_dp)*Hr4(0,0,1,0) + z**2*(-1670.320987654321_dp)*Hr4(0,0,1,0) + &
                        nf*(-16.59259259259259_dp)*Hr4(0,0,1,0) + nf*z*(125.62962962962963_dp)*Hr4(0,0,1,0) + &
                        ((256._dp)*Hr4(0,0,1,0))/z + (1086.8148148148148_dp)*Hr4(0,0,1,0) + z*(-3309.037037037037_dp)*Hr4(0,0,1,1) + &
                        (-3070.814814814815_dp)*Hr4(0,0,1,1) + z**2*(-1703.5061728395062_dp)*Hr4(0,0,1,1) + &
                        nf*z*(182.5185185185185_dp)*Hr4(0,0,1,1) + nf*(246.5185185185185_dp)*Hr4(0,0,1,1) + &
                        ((384._dp)*Hr4(0,0,1,1))/z + z**2*(-1286.320987654321_dp)*Hr4(0,1,0,0) + &
                        z*(-865.7777777777778_dp)*Hr4(0,1,0,0) + nf*(-120.88888888888889_dp)*Hr4(0,1,0,0) + &
                        nf*z*(-120.88888888888889_dp)*Hr4(0,1,0,0) + ((625.7777777777778_dp)*Hr4(0,1,0,0))/z + &
                        (837.3333333333334_dp)*Hr4(0,1,0,0) + z*(-2935.703703703704_dp)*Hr4(0,1,0,1) + &
                        (-2212.740740740741_dp)*Hr4(0,1,0,1) + z**2*(-982.9135802469136_dp)*Hr4(0,1,0,1) + &
                        ((-199.11111111111111_dp)*Hr4(0,1,0,1))/z + nf*(92.44444444444444_dp)*Hr4(0,1,0,1) + &
                        nf*z*(92.44444444444444_dp)*Hr4(0,1,0,1) + z*(-3172.740740740741_dp)*Hr4(0,1,1,0) + &
                        (-2603.8518518518517_dp)*Hr4(0,1,1,0) + z**2*(-919.7037037037037_dp)*Hr4(0,1,1,0) + &
                        ((-341.3333333333333_dp)*Hr4(0,1,1,0))/z + nf*(116.14814814814815_dp)*Hr4(0,1,1,0) + &
                        nf*z*(116.14814814814815_dp)*Hr4(0,1,1,0) + ((-625.7777777777778_dp)*Hr4(0,1,1,1))/z + &
                        nf*(-92.44444444444444_dp)*Hr4(0,1,1,1) + nf*z*(-92.44444444444444_dp)*Hr4(0,1,1,1) + &
                        z**2*(339.75308641975306_dp)*Hr4(0,1,1,1) + (600.8888888888889_dp)*Hr4(0,1,1,1) + &
                        z*(1556.148148148148_dp)*Hr4(0,1,1,1) + ((-771.1604938271605_dp)*Hr4(1,0,-1,0))/z + &
                        (-578.3703703703703_dp)*Hr4(1,0,-1,0) + z*(578.3703703703703_dp)*Hr4(1,0,-1,0) + &
                        z**2*(771.1604938271605_dp)*Hr4(1,0,-1,0) + z**2*(-1365.3333333333333_dp)*Hr4(1,0,0,0) + &
                        z*(-1194.6666666666667_dp)*Hr4(1,0,0,0) + (1194.6666666666667_dp)*Hr4(1,0,0,0) + &
                        ((1365.3333333333333_dp)*Hr4(1,0,0,0))/z + z*(-839.1111111111111_dp)*Hr4(1,0,0,1) + &
                        z**2*(-739.5555555555555_dp)*Hr4(1,0,0,1) + ((739.5555555555555_dp)*Hr4(1,0,0,1))/z + &
                        (839.1111111111111_dp)*Hr4(1,0,0,1) + z*(-355.55555555555554_dp)*Hr4(1,0,1,0) + &
                        z**2*(-170.66666666666666_dp)*Hr4(1,0,1,0) + ((170.66666666666666_dp)*Hr4(1,0,1,0))/z + &
                        (355.55555555555554_dp)*Hr4(1,0,1,0) + ((-598.9135802469136_dp)*Hr4(1,0,1,1))/z + &
                        (-221.62962962962962_dp)*Hr4(1,0,1,1) + z*(221.62962962962962_dp)*Hr4(1,0,1,1) + &
                        z**2*(598.9135802469136_dp)*Hr4(1,0,1,1) + z*(-585.4814814814815_dp)*Hr4(1,1,0,0) + &
                        z**2*(-401.38271604938274_dp)*Hr4(1,1,0,0) + ((401.38271604938274_dp)*Hr4(1,1,0,0))/z + &
                        (585.4814814814815_dp)*Hr4(1,1,0,0) + ((-772.7407407407408_dp)*Hr4(1,1,0,1))/z + (-352._dp)*Hr4(1,1,0,1) + &
                        z*(352._dp)*Hr4(1,1,0,1) + z**2*(772.7407407407408_dp)*Hr4(1,1,0,1) + &
                        ((-978.1728395061729_dp)*Hr4(1,1,1,0))/z + (-506.0740740740741_dp)*Hr4(1,1,1,0) + &
                        z*(506.0740740740741_dp)*Hr4(1,1,1,0) + z**2*(978.1728395061729_dp)*Hr4(1,1,1,0) + &
                        ((-410.8641975308642_dp)*Hr4(1,1,1,1))/z + (-308.14814814814815_dp)*Hr4(1,1,1,1) + &
                        z*(308.14814814814815_dp)*Hr4(1,1,1,1) + z**2*(410.8641975308642_dp)*Hr4(1,1,1,1) + &
                        (-786.9629629629629_dp)*Hr5(0,-1,-1,-1,0) + z*(786.9629629629629_dp)*Hr5(0,-1,-1,-1,0) + &
                        z*(-4439.7037037037035_dp)*Hr5(0,-1,-1,0,0) + (4439.7037037037035_dp)*Hr5(0,-1,-1,0,0) + &
                        z*(-3048.296296296296_dp)*Hr5(0,-1,-1,0,1) + (3048.296296296296_dp)*Hr5(0,-1,-1,0,1) + &
                        z*(-1858.3703703703704_dp)*Hr5(0,-1,0,-1,0) + (1858.3703703703704_dp)*Hr5(0,-1,0,-1,0) + &
                        (-3451.259259259259_dp)*Hr5(0,-1,0,0,0) + z*(3451.259259259259_dp)*Hr5(0,-1,0,0,0) + &
                        (-3316.1481481481483_dp)*Hr5(0,-1,0,0,1) + z*(3316.1481481481483_dp)*Hr5(0,-1,0,0,1) + &
                        (-1099.851851851852_dp)*Hr5(0,-1,0,1,0) + z*(1099.851851851852_dp)*Hr5(0,-1,0,1,0) + &
                        z*(-1294.2222222222222_dp)*Hr5(0,0,-1,-1,0) + (4285.62962962963_dp)*Hr5(0,0,-1,-1,0) + &
                        (-3344.5925925925926_dp)*Hr5(0,0,-1,0,0) + z*(5560.888888888889_dp)*Hr5(0,0,-1,0,0) + &
                        (-4489.481481481482_dp)*Hr5(0,0,-1,0,1) + z*(2697.4814814814813_dp)*Hr5(0,0,-1,0,1) + &
                        z*(2258.962962962963_dp)*Hr5(0,0,0,-1,0) + (2576.5925925925926_dp)*Hr5(0,0,0,-1,0) + &
                        nf*(-156.44444444444446_dp)*Hr5(0,0,0,0,0) + nf*z*(-156.44444444444446_dp)*Hr5(0,0,0,0,0) + &
                        (1533.6296296296296_dp)*Hr5(0,0,0,0,0) + z*(4899.555555555556_dp)*Hr5(0,0,0,0,0) + &
                        (-1278.8148148148148_dp)*Hr5(0,0,0,0,1) + z*(4621.037037037037_dp)*Hr5(0,0,0,0,1) + &
                        (284.44444444444446_dp)*Hr5(0,0,0,1,0) + z*(3323.259259259259_dp)*Hr5(0,0,0,1,0) + &
                        z*(2256.5925925925926_dp)*Hr5(0,0,0,1,1) + (2832.5925925925926_dp)*Hr5(0,0,0,1,1) + &
                        (382.81481481481484_dp)*Hr5(0,0,1,0,0) + z*(3307.8518518518517_dp)*Hr5(0,0,1,0,0) + &
                        z*(343.7037037037037_dp)*Hr5(0,0,1,0,1) + (1474.3703703703704_dp)*Hr5(0,0,1,0,1) + &
                        z*(-130.37037037037038_dp)*Hr5(0,0,1,1,0) + (1640.2962962962963_dp)*Hr5(0,0,1,1,0) + &
                        z*(-1069.037037037037_dp)*Hr5(0,0,1,1,1) + (-749.0370370370371_dp)*Hr5(0,0,1,1,1) + &
                        (-1156.7407407407406_dp)*Hr5(0,1,0,-1,0) + z*(-1156.7407407407406_dp)*Hr5(0,1,0,-1,0) + &
                        (2161.777777777778_dp)*Hr5(0,1,0,0,0) + z*(2161.777777777778_dp)*Hr5(0,1,0,0,0) + &
                        (1298.962962962963_dp)*Hr5(0,1,0,0,1) + z*(1298.962962962963_dp)*Hr5(0,1,0,0,1) + &
                        (407.7037037037037_dp)*Hr5(0,1,0,1,0) + z*(407.7037037037037_dp)*Hr5(0,1,0,1,0) + &
                        (-746.6666666666666_dp)*Hr5(0,1,0,1,1) + z*(-746.6666666666666_dp)*Hr5(0,1,0,1,1) + &
                        (791.7037037037037_dp)*Hr5(0,1,1,0,0) + z*(791.7037037037037_dp)*Hr5(0,1,1,0,0) + &
                        (-1007.4074074074074_dp)*Hr5(0,1,1,0,1) + z*(-1007.4074074074074_dp)*Hr5(0,1,1,0,1) + &
                        (-1315.5555555555557_dp)*Hr5(0,1,1,1,0) + z*(-1315.5555555555557_dp)*Hr5(0,1,1,1,0) + &
                        (-616.2962962962963_dp)*Hr5(0,1,1,1,1) + z*(-616.2962962962963_dp)*Hr5(0,1,1,1,1)

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (3)
                select case (type_in)
                    case (1)
                        Ibar_select = z*(-7434.326224233369_dp) + (-5155.828741530615_dp)/z + nf*(-495.2047887829962_dp) + &
                        nf*z**2*(-319.11655903545085_dp) + nf**2*z*(-4.397193323026427_dp) + (nf**2*(-0.7901234567901234_dp))/z + &
                        nf**2*z**2*(0.7901234567901234_dp) + nf**2*(12.195399269566165_dp) + (nf*(133.0382998603344_dp))/z + &
                        nf*z*(591.4711093186862_dp) + (13120.599418865742_dp) + z**2*(16787.908742689935_dp) + &
                        ((-858.6691214467646_dp)/z + z**2*(-858.6691214467646_dp) + (644.0018410850735_dp) + &
                        z*(644.0018410850735_dp))*Hr1(-1) + (-4241.7244158687845_dp)*Hr1(0) + z**2*(-3702.6734113000393_dp)*Hr1(0) + &
                        z*(-2475.1117665035963_dp)*Hr1(0) + ((-1561.7472233738827_dp)*Hr1(0))/z + nf*(-68.37954716676218_dp)*Hr1(0) &
                        + nf*z*(-58.600941725110154_dp)*Hr1(0) + nf**2*z**2*(-1.5802469135802468_dp)*Hr1(0) + &
                        nf**2*(5.135802469135802_dp)*Hr1(0) + nf**2*z*(5.135802469135802_dp)*Hr1(0) + &
                        (nf*(25.371742112482853_dp)*Hr1(0))/z + nf*z**2*(203.36899862825788_dp)*Hr1(0) + &
                        z*(-5313.752610926626_dp)*Hr1(1) + ((-2554.700222468202_dp)*Hr1(1))/z + nf*(-125.36625514403292_dp)*Hr1(1) + &
                        nf*z**2*(-92.3127572016461_dp)*Hr1(1) + (nf**2*(-1.5802469135802468_dp)*Hr1(1))/z + &
                        nf**2*(-1.1851851851851851_dp)*Hr1(1) + nf**2*z*(1.1851851851851851_dp)*Hr1(1) + &
                        nf**2*z**2*(1.5802469135802468_dp)*Hr1(1) + (nf*(78.88065843621399_dp)*Hr1(1))/z + &
                        nf*z*(138.79835390946502_dp)*Hr1(1) + (3359.908989527449_dp)*Hr1(1) + z**2*(4508.543843867379_dp)*Hr1(1) + &
                        (nf*(-33.18518518518518_dp)*Hr2(-1,0))/z + nf*z**2*(-33.18518518518518_dp)*Hr2(-1,0) + &
                        nf*(24.88888888888889_dp)*Hr2(-1,0) + nf*z*(24.88888888888889_dp)*Hr2(-1,0) + &
                        (1180.1152263374486_dp)*Hr2(-1,0) + z*(1180.1152263374486_dp)*Hr2(-1,0) + &
                        ((1479.9012345679012_dp)*Hr2(-1,0))/z + z**2*(1479.9012345679012_dp)*Hr2(-1,0) + &
                        z*(-1288.003682170147_dp)*Hr2(0,-1) + (1288.003682170147_dp)*Hr2(0,-1) + &
                        z**2*(-2090.1399176954733_dp)*Hr2(0,0) + z*(-206.90435107440533_dp)*Hr2(0,0) + &
                        nf*(-203.85185185185185_dp)*Hr2(0,0) + ((-184.88888888888889_dp)*Hr2(0,0))/z + &
                        nf*z*(-165.1358024691358_dp)*Hr2(0,0) + nf**2*(2.3703703703703702_dp)*Hr2(0,0) + &
                        nf**2*z*(2.3703703703703702_dp)*Hr2(0,0) + nf*z**2*(25.28395061728395_dp)*Hr2(0,0) + &
                        (2142.9126400932455_dp)*Hr2(0,0) + z**2*(-2136.230452674897_dp)*Hr2(0,1) + (-1517.7458069533318_dp)*Hr2(0,1) &
                        + z*(-704.7746135376939_dp)*Hr2(0,1) + nf*(-18.962962962962962_dp)*Hr2(0,1) + &
                        nf*z*(-9.876543209876543_dp)*Hr2(0,1) + (nf*(-4.7407407407407405_dp)*Hr2(0,1))/z + &
                        nf**2*(-2.3703703703703702_dp)*Hr2(0,1) + nf**2*z*(-2.3703703703703702_dp)*Hr2(0,1) + &
                        nf*z**2*(27.91769547325103_dp)*Hr2(0,1) + ((98.76543209876543_dp)*Hr2(0,1))/z + &
                        (-978.238683127572_dp)*Hr2(1,0) + z**2*(-700.312757201646_dp)*Hr2(1,0) + &
                        (nf*(-8.954732510288066_dp)*Hr2(1,0))/z + nf*(-6.716049382716049_dp)*Hr2(1,0) + &
                        nf*z*(6.716049382716049_dp)*Hr2(1,0) + nf*z**2*(8.954732510288066_dp)*Hr2(1,0) + &
                        z*(768.0658436213992_dp)*Hr2(1,0) + ((910.4855967078189_dp)*Hr2(1,0))/z + ((-638.6831275720165_dp)/z + &
                        z*(-274.7654320987654_dp) + (-145.58024691358025_dp) + nf*(z**2*(-38.452674897119344_dp) + &
                        z*(-28.839506172839506_dp) + (28.839506172839506_dp) + (38.452674897119344_dp)/z) + &
                        z**2*(1059.0288065843622_dp))*Hr2(1,1) + ((-399.275720164609_dp)*Hr3(-1,-1,0))/z + &
                        z**2*(-399.275720164609_dp)*Hr3(-1,-1,0) + (299.4567901234568_dp)*Hr3(-1,-1,0) + &
                        z*(299.4567901234568_dp)*Hr3(-1,-1,0) + (-270.61728395061726_dp)*Hr3(-1,0,0) + &
                        z*(-270.61728395061726_dp)*Hr3(-1,0,0) + ((360.8230452674897_dp)*Hr3(-1,0,0))/z + &
                        z**2*(360.8230452674897_dp)*Hr3(-1,0,0) + (-241.77777777777777_dp)*Hr3(-1,0,1) + &
                        z*(-241.77777777777777_dp)*Hr3(-1,0,1) + ((322.3703703703704_dp)*Hr3(-1,0,1))/z + &
                        z**2*(322.3703703703704_dp)*Hr3(-1,0,1) + (-811.8518518518518_dp)*Hr3(0,-1,0) + &
                        nf*z*(-49.77777777777778_dp)*Hr3(0,-1,0) + nf*(49.77777777777778_dp)*Hr3(0,-1,0) + &
                        z**2*(143.27572016460906_dp)*Hr3(0,-1,0) + ((256._dp)*Hr3(0,-1,0))/z + z*(512.395061728395_dp)*Hr3(0,-1,0) + &
                        (-459.85185185185185_dp)*Hr3(0,0,0) + z**2*(-168.559670781893_dp)*Hr3(0,0,0) + &
                        nf*z*(-56.888888888888886_dp)*Hr3(0,0,0) + nf*(-4.7407407407407405_dp)*Hr3(0,0,0) + &
                        z*(730.074074074074_dp)*Hr3(0,0,0) + z**2*(-383.4732510288066_dp)*Hr3(0,0,1) + &
                        z*(-262.71604938271605_dp)*Hr3(0,0,1) + nf*(-70.32098765432099_dp)*Hr3(0,0,1) + &
                        nf*z*(-27.65432098765432_dp)*Hr3(0,0,1) + ((170.66666666666666_dp)*Hr3(0,0,1))/z + &
                        (1263.0123456790122_dp)*Hr3(0,0,1) + z**2*(-405.59670781893004_dp)*Hr3(0,1,0) + &
                        z*(-182.91358024691357_dp)*Hr3(0,1,0) + nf*(-13.432098765432098_dp)*Hr3(0,1,0) + &
                        nf*z*(-13.432098765432098_dp)*Hr3(0,1,0) + ((369.77777777777777_dp)*Hr3(0,1,0))/z + &
                        (398.61728395061726_dp)*Hr3(0,1,0) + z*(-1232.9876543209878_dp)*Hr3(0,1,1) + &
                        (-597.7283950617284_dp)*Hr3(0,1,1) + z**2*(-524.641975308642_dp)*Hr3(0,1,1) + &
                        nf*(57.67901234567901_dp)*Hr3(0,1,1) + nf*z*(57.67901234567901_dp)*Hr3(0,1,1) + &
                        ((322.3703703703704_dp)*Hr3(0,1,1))/z + z**2*(-313.4156378600823_dp)*Hr3(1,0,0) + &
                        z*(-235.06172839506172_dp)*Hr3(1,0,0) + (235.06172839506172_dp)*Hr3(1,0,0) + &
                        ((313.4156378600823_dp)*Hr3(1,0,0))/z + z**2*(-351.86831275720164_dp)*Hr3(1,0,1) + &
                        z*(-263.9012345679012_dp)*Hr3(1,0,1) + (263.9012345679012_dp)*Hr3(1,0,1) + &
                        ((351.86831275720164_dp)*Hr3(1,0,1))/z + z**2*(-351.86831275720164_dp)*Hr3(1,1,0) + &
                        z*(-263.9012345679012_dp)*Hr3(1,1,0) + (263.9012345679012_dp)*Hr3(1,1,0) + &
                        ((351.86831275720164_dp)*Hr3(1,1,0))/z + z*(-598.9135802469136_dp)*Hr4(0,-1,-1,0) + &
                        (598.9135802469136_dp)*Hr4(0,-1,-1,0) + (-541.2345679012345_dp)*Hr4(0,-1,0,0) + &
                        z*(541.2345679012345_dp)*Hr4(0,-1,0,0) + (-483.55555555555554_dp)*Hr4(0,-1,0,1) + &
                        z*(483.55555555555554_dp)*Hr4(0,-1,0,1) + (276.5432098765432_dp)*Hr4(0,0,-1,0) + &
                        z*(875.4567901234568_dp)*Hr4(0,0,-1,0) + nf*(-28.444444444444443_dp)*Hr4(0,0,0,0) + &
                        nf*z*(-28.444444444444443_dp)*Hr4(0,0,0,0) + (270.22222222222223_dp)*Hr4(0,0,0,0) + &
                        z*(1058.7654320987654_dp)*Hr4(0,0,0,0) + (-505.679012345679_dp)*Hr4(0,0,0,1) + &
                        z*(1400.0987654320988_dp)*Hr4(0,0,0,1) + (219.65432098765433_dp)*Hr4(0,0,1,0) + &
                        z*(1158.320987654321_dp)*Hr4(0,0,1,0) + (850.1728395061729_dp)*Hr4(0,0,1,1) + &
                        z*(850.1728395061729_dp)*Hr4(0,0,1,1) + (470.12345679012344_dp)*Hr4(0,1,0,0) + &
                        z*(470.12345679012344_dp)*Hr4(0,1,0,0) + (527.8024691358024_dp)*Hr4(0,1,0,1) + &
                        z*(527.8024691358024_dp)*Hr4(0,1,0,1) + (527.8024691358024_dp)*Hr4(0,1,1,0) + &
                        z*(527.8024691358024_dp)*Hr4(0,1,1,0)

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (4)
                select case (type_in)
                    case (1)
                        Ibar_select = (-932.2479503960094_dp) + z**2*(-690.3278556100735_dp) + (nf*(-34.28257887517147_dp))/z + &
                        nf**2*z**2*(-0.3950617283950617_dp) + nf**2*z*(-0.2962962962962963_dp) + nf**2*(0.2962962962962963_dp) + &
                        (nf**2*(0.3950617283950617_dp))/z + nf*z*(24.266335291614528_dp) + nf*(25.12230236980383_dp) + &
                        nf*z**2*(34.28257887517147_dp) + z*(57.05093411498649_dp) + (661.7344338060944_dp)/z + &
                        z**2*(-195.55555555555554_dp)*Hr1(0) + nf*(-32._dp)*Hr1(0) + nf*z*(-24.493827160493826_dp)*Hr1(0) + &
                        (nf*(-5.794238683127572_dp)*Hr1(0))/z + nf**2*(0.5925925925925926_dp)*Hr1(0) + &
                        nf**2*z*(0.5925925925925926_dp)*Hr1(0) + nf*z**2*(4.2139917695473255_dp)*Hr1(0) + &
                        z*(126.66524250855474_dp)*Hr1(0) + ((195.55555555555554_dp)*Hr1(0))/z + (338.9096730971788_dp)*Hr1(0) + &
                        z**2*(-391.1111111111111_dp)*Hr1(1) + (-243.9506172839506_dp)*Hr1(1) + &
                        (nf*(-10.008230452674898_dp)*Hr1(1))/z + nf*(-7.506172839506172_dp)*Hr1(1) + &
                        nf*z*(7.506172839506172_dp)*Hr1(1) + nf*z**2*(10.008230452674898_dp)*Hr1(1) + &
                        z*(243.9506172839506_dp)*Hr1(1) + ((391.1111111111111_dp)*Hr1(1))/z + (-63.407407407407405_dp)*Hr2(0,0) + &
                        nf*z*(-15.802469135802468_dp)*Hr2(0,0) + z**2*(-12.641975308641975_dp)*Hr2(0,0) + &
                        nf*(0.7901234567901234_dp)*Hr2(0,0) + ((21.333333333333332_dp)*Hr2(0,0))/z + &
                        z*(271.4074074074074_dp)*Hr2(0,0) + z**2*(-44.24691358024691_dp)*Hr2(0,1) + &
                        nf*(-15.012345679012345_dp)*Hr2(0,1) + nf*z*(-15.012345679012345_dp)*Hr2(0,1) + &
                        ((61.629629629629626_dp)*Hr2(0,1))/z + z*(168.2962962962963_dp)*Hr2(0,1) + (247.7037037037037_dp)*Hr2(0,1) + &
                        z**2*(-52.93827160493827_dp)*Hr2(1,0) + z*(-39.7037037037037_dp)*Hr2(1,0) + (39.7037037037037_dp)*Hr2(1,0) + &
                        ((52.93827160493827_dp)*Hr2(1,0))/z + (z**2*(-105.87654320987654_dp) + z*(-79.4074074074074_dp) + &
                        (79.4074074074074_dp) + (105.87654320987654_dp)/z)*Hr2(1,1) + (nf*(-2.3703703703703702_dp) + &
                        nf*z*(-2.3703703703703702_dp) + (22.51851851851852_dp) + z*(161.1851851851852_dp))*Hr3(0,0,0) + &
                        (-59.25925925925926_dp)*Hr3(0,0,1) + z*(218.07407407407408_dp)*Hr3(0,0,1) + (79.4074074074074_dp)*Hr3(0,1,0) &
                        + z*(79.4074074074074_dp)*Hr3(0,1,0) + (158.8148148148148_dp)*Hr3(0,1,1) + &
                        z*(158.8148148148148_dp)*Hr3(0,1,1)

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (5)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (6)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case default
                Ibar_select = 0._dp
        end select
    case (5)
        select case (powLperp)
            case (2)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (3)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (4)
                select case (type_in)
                    case (1)
                        Ibar_select = (-107237.23285764939_dp)/z + z*(-96657.4524285029_dp) + nf*z**2*(-19010.08237291654_dp) + &
                        nf*(-15375.893641597377_dp) + nf**2*z*(-618.6479412412657_dp) + (nf**2*(-90.31270205667106_dp))/z + &
                        nf**3*(-8.130266179710777_dp) + nf**3*z**2*(-0.5267489711934157_dp) + (nf**3*(0.5267489711934157_dp))/z + &
                        nf**3*z*(2.9314622153509506_dp) + nf**2*z**2*(218.97174416874213_dp) + nf**2*(686.3844640805536_dp) + &
                        (nf*(5930.765024741579_dp))/z + nf*z*(15493.327598901551_dp) + (307249.7397778502_dp) + &
                        z**2*(428328.0097519824_dp) + ((-43386.49376656064_dp)/z + z**2*(-43386.49376656064_dp) + &
                        (-20951.490131705166_dp) + z*(-20951.490131705166_dp) + nf*((-689.9246094369187_dp) + &
                        z*(-689.9246094369187_dp) + (919.8994792492248_dp)/z + z**2*(919.8994792492248_dp)))*Hr1(-1) + &
                        (-107217.65784948491_dp)*Hr1(0) + z*(-59095.240716199405_dp)*Hr1(0) + ((-37135.86714933466_dp)*Hr1(0))/z + &
                        z**2*(-31316.683156337313_dp)*Hr1(0) + nf**2*z**2*(-165.17969821673526_dp)*Hr1(0) + &
                        (nf**2*(-16.095107453132144_dp)*Hr1(0))/z + nf**3*(-3.4238683127572016_dp)*Hr1(0) + &
                        nf**3*z*(-3.4238683127572016_dp)*Hr1(0) + nf**3*z**2*(1.0534979423868314_dp)*Hr1(0) + &
                        nf**2*z*(204.681567187674_dp)*Hr1(0) + nf**2*(246.1607583188965_dp)*Hr1(0) + &
                        (nf*(1574.4864368705591_dp)*Hr1(0))/z + nf*z*(2767.8560486661927_dp)*Hr1(0) + &
                        nf*(4840.747927481836_dp)*Hr1(0) + nf*z**2*(6780.612470334416_dp)*Hr1(0) + z*(-166845.52762410985_dp)*Hr1(1) &
                        + ((-84150.20390796897_dp)*Hr1(1))/z + nf*(-6244.310981130706_dp)*Hr1(1) + &
                        nf*z**2*(-4521.370452228104_dp)*Hr1(1) + (nf**2*(-79.2318244170096_dp)*Hr1(1))/z + &
                        nf**2*z*(-65.7119341563786_dp)*Hr1(1) + nf**3*z**2*(-1.0534979423868314_dp)*Hr1(1) + &
                        nf**3*z*(-0.7901234567901234_dp)*Hr1(1) + nf**3*(0.7901234567901234_dp)*Hr1(1) + &
                        (nf**3*(1.0534979423868314_dp)*Hr1(1))/z + nf**2*(60.181069958847736_dp)*Hr1(1) + &
                        nf**2*z**2*(84.76268861454047_dp)*Hr1(1) + (nf*(2947.664005040175_dp)*Hr1(1))/z + &
                        nf*z*(7818.017428318634_dp)*Hr1(1) + (111976.93730964694_dp)*Hr1(1) + z**2*(139018.79422243187_dp)*Hr1(1) + &
                        (-6381.5318662516875_dp)*Hr2(-1,-1) + z*(-6381.5318662516875_dp)*Hr2(-1,-1) + &
                        ((8508.709155002249_dp)*Hr2(-1,-1))/z + z**2*(8508.709155002249_dp)*Hr2(-1,-1) + &
                        (nf*(-1898.9300411522634_dp)*Hr2(-1,0))/z + nf*z**2*(-1898.9300411522634_dp)*Hr2(-1,0) + &
                        nf*(-231.02331961591221_dp)*Hr2(-1,0) + nf*z*(-231.02331961591221_dp)*Hr2(-1,0) + &
                        nf**2*(-17.77777777777778_dp)*Hr2(-1,0) + nf**2*z*(-17.77777777777778_dp)*Hr2(-1,0) + &
                        (nf**2*(23.703703703703702_dp)*Hr2(-1,0))/z + nf**2*z**2*(23.703703703703702_dp)*Hr2(-1,0) + &
                        ((32947.38054981388_dp)*Hr2(-1,0))/z + z**2*(32947.38054981388_dp)*Hr2(-1,0) + &
                        (37232.5756987507_dp)*Hr2(-1,0) + z*(37232.5756987507_dp)*Hr2(-1,0) + z*(-26387.642145092956_dp)*Hr2(0,-1) + &
                        ((-5874.648479726604_dp)*Hr2(0,-1))/z + z**2*(-2634.0606752756453_dp)*Hr2(0,-1) + &
                        nf*(-1379.8492188738373_dp)*Hr2(0,-1) + nf*z*(1379.8492188738373_dp)*Hr2(0,-1) + &
                        (32769.17401134464_dp)*Hr2(0,-1) + z**2*(-61870.87779773989_dp)*Hr2(0,0) + z*(-11220.7378394324_dp)*Hr2(0,0) &
                        + ((-7227.46389234387_dp)*Hr2(0,0))/z + nf*(-4048.326875186134_dp)*Hr2(0,0) + &
                        nf*z*(-2487.6901131257973_dp)*Hr2(0,0) + nf**2*z**2*(-16.15363511659808_dp)*Hr2(0,0) + &
                        nf**3*(-1.5802469135802468_dp)*Hr2(0,0) + nf**3*z*(-1.5802469135802468_dp)*Hr2(0,0) + &
                        nf**2*z*(152.92710038413804_dp)*Hr2(0,0) + (nf*(168.8230452674897_dp)*Hr2(0,0))/z + &
                        nf**2*(246.29335552817096_dp)*Hr2(0,0) + nf*z**2*(2268.3859167809787_dp)*Hr2(0,0) + &
                        (13558.756999334992_dp)*Hr2(0,0) + z**2*(-48264.867526517206_dp)*Hr2(0,1) + &
                        z*(-40107.344922655466_dp)*Hr2(0,1) + (-36030.79214784408_dp)*Hr2(0,1) + ((-10847.4903471162_dp)*Hr2(0,1))/z &
                        + nf*(-1005.1432242641541_dp)*Hr2(0,1) + (nf*(-29.93689986282579_dp)*Hr2(0,1))/z + &
                        nf**2*(-19.665294924554185_dp)*Hr2(0,1) + nf**2*z**2*(-19.13854595336077_dp)*Hr2(0,1) + &
                        nf**2*z*(-18.348422496570645_dp)*Hr2(0,1) + (nf**2*(-1.2290809327846366_dp)*Hr2(0,1))/z + &
                        nf**3*(1.5802469135802468_dp)*Hr2(0,1) + nf**3*z*(1.5802469135802468_dp)*Hr2(0,1) + &
                        nf*z*(261.60025996081157_dp)*Hr2(0,1) + nf*z**2*(2378.534979423868_dp)*Hr2(0,1) + &
                        (-14190.406911222282_dp)*Hr2(1,0) + z**2*(-7778.877662037416_dp)*Hr2(1,0) + &
                        (nf*(-883.387288523091_dp)*Hr2(1,0))/z + nf*(-349.71742112482855_dp)*Hr2(1,0) + &
                        nf**2*z**2*(-3.6872427983539096_dp)*Hr2(1,0) + nf**2*z*(-2.765432098765432_dp)*Hr2(1,0) + &
                        nf**2*(2.765432098765432_dp)*Hr2(1,0) + (nf**2*(3.6872427983539096_dp)*Hr2(1,0))/z + &
                        nf*z*(472.71330589849106_dp)*Hr2(1,0) + nf*z**2*(760.3914037494285_dp)*Hr2(1,0) + &
                        z*(3569.8307795350383_dp)*Hr2(1,0) + ((18399.453793724657_dp)*Hr2(1,0))/z + &
                        ((36.69684499314129_dp)*(z**2*(-1378.6831536252023_dp) + (-740.7924777820269_dp) + nf**2*((-1._dp) + &
                        z*(-0.75_dp) + z**2*(0.75_dp) + z**3*(1._dp)) + nf*(z**3*(-49.698564593301434_dp) + &
                        z*(-30.211722488038276_dp) + z**2*(36.915071770334926_dp) + (42.995215311004785_dp)) + &
                        z*(799.8554024290299_dp) + z**3*(1319.620228978199_dp))*Hr2(1,1))/z + &
                        ((-16857.54732510288_dp)*Hr3(-1,-1,0))/z + z**2*(-16857.54732510288_dp)*Hr3(-1,-1,0) + &
                        (-12013.958847736625_dp)*Hr3(-1,-1,0) + z*(-12013.958847736625_dp)*Hr3(-1,-1,0) + &
                        nf*(-317.36625514403295_dp)*Hr3(-1,-1,0) + nf*z*(-317.36625514403295_dp)*Hr3(-1,-1,0) + &
                        (nf*(423.1550068587106_dp)*Hr3(-1,-1,0))/z + nf*z**2*(423.1550068587106_dp)*Hr3(-1,-1,0) + &
                        (nf*(-385.40466392318245_dp)*Hr3(-1,0,0))/z + nf*z**2*(-385.40466392318245_dp)*Hr3(-1,0,0) + &
                        nf*(289.05349794238685_dp)*Hr3(-1,0,0) + nf*z*(289.05349794238685_dp)*Hr3(-1,0,0) + &
                        (10916.938271604939_dp)*Hr3(-1,0,0) + z*(10916.938271604939_dp)*Hr3(-1,0,0) + &
                        ((15342.35390946502_dp)*Hr3(-1,0,0))/z + z**2*(15342.35390946502_dp)*Hr3(-1,0,0) + &
                        (nf*(-347.65432098765433_dp)*Hr3(-1,0,1))/z + nf*z**2*(-347.65432098765433_dp)*Hr3(-1,0,1) + &
                        nf*(260.74074074074076_dp)*Hr3(-1,0,1) + nf*z*(260.74074074074076_dp)*Hr3(-1,0,1) + &
                        (9819.91769547325_dp)*Hr3(-1,0,1) + z*(9819.91769547325_dp)*Hr3(-1,0,1) + &
                        ((13827.160493827161_dp)*Hr3(-1,0,1))/z + z**2*(13827.160493827161_dp)*Hr3(-1,0,1) + &
                        (-12763.063732503375_dp)*Hr3(0,-1,-1) + z*(12763.063732503375_dp)*Hr3(0,-1,-1) + &
                        (-4826.63461072905_dp)*Hr3(0,-1,0) + nf*z*(-1442.238683127572_dp)*Hr3(0,-1,0) + &
                        (nf*(-306.5679012345679_dp)*Hr3(0,-1,0))/z + nf*z**2*(-116.58710562414267_dp)*Hr3(0,-1,0) + &
                        nf**2*(-35.55555555555556_dp)*Hr3(0,-1,0) + nf**2*z*(35.55555555555556_dp)*Hr3(0,-1,0) + &
                        nf*(1759.6049382716049_dp)*Hr3(0,-1,0) + z**2*(7020.510288065843_dp)*Hr3(0,-1,0) + &
                        ((9837.037037037036_dp)*Hr3(0,-1,0))/z + z*(16840.593458465675_dp)*Hr3(0,-1,0) + &
                        z*(-19599.490945636546_dp)*Hr3(0,0,-1) + (-6836.427213133172_dp)*Hr3(0,0,-1) + &
                        (-19253.7677596292_dp)*Hr3(0,0,0) + z**2*(-9823.45861911294_dp)*Hr3(0,0,0) + &
                        nf*z*(-2418.406368571906_dp)*Hr3(0,0,0) + ((-554.6666666666666_dp)*Hr3(0,0,0))/z + &
                        nf**2*(30.551440329218106_dp)*Hr3(0,0,0) + nf**2*z*(54.25514403292181_dp)*Hr3(0,0,0) + &
                        nf*z**2*(137.6570644718793_dp)*Hr3(0,0,0) + nf*(1116.7670350503977_dp)*Hr3(0,0,0) + &
                        z*(9550.786750661508_dp)*Hr3(0,0,0) + z**2*(-20674.194787379973_dp)*Hr3(0,0,1) + &
                        nf*(-2723.1604938271603_dp)*Hr3(0,0,1) + nf*z*(-1618.1728395061727_dp)*Hr3(0,0,1) + &
                        (nf*(-153.28395061728395_dp)*Hr3(0,0,1))/z + nf**2*z*(22.650205761316872_dp)*Hr3(0,0,1) + &
                        nf**2*(45.56378600823045_dp)*Hr3(0,0,1) + nf*z**2*(293.5747599451303_dp)*Hr3(0,0,1) + &
                        z*(1691.9187584707015_dp)*Hr3(0,0,1) + ((3496.296296296296_dp)*Hr3(0,0,1))/z + &
                        (15234.588458676459_dp)*Hr3(0,0,1) + z**2*(-18438.935528120714_dp)*Hr3(0,1,0) + &
                        nf*(-1212.8395061728395_dp)*Hr3(0,1,0) + (-1125.0662229932595_dp)*Hr3(0,1,0) + &
                        nf*z*(-801.5802469135803_dp)*Hr3(0,1,0) + (nf*(-320.9657064471879_dp)*Hr3(0,1,0))/z + &
                        nf**2*(5.530864197530864_dp)*Hr3(0,1,0) + nf**2*z*(5.530864197530864_dp)*Hr3(0,1,0) + &
                        nf*z**2*(320.7901234567901_dp)*Hr3(0,1,0) + ((11371.456790123457_dp)*Hr3(0,1,0))/z + &
                        z*(14702.489332562294_dp)*Hr3(0,1,0) + z**2*(-22751.517146776405_dp)*Hr3(0,1,1) + &
                        (-20207.046019804086_dp)*Hr3(0,1,1) + z*(-8489.91021733495_dp)*Hr3(0,1,1) + &
                        (nf*(-126.41975308641975_dp)*Hr3(0,1,1))/z + nf**2*(-55.04526748971193_dp)*Hr3(0,1,1) + &
                        nf**2*z*(-55.04526748971193_dp)*Hr3(0,1,1) + nf*z**2*(346.6008230452675_dp)*Hr3(0,1,1) + &
                        nf*(618.6666666666666_dp)*Hr3(0,1,1) + nf*z*(833.3168724279835_dp)*Hr3(0,1,1) + &
                        ((5742.617283950617_dp)*Hr3(0,1,1))/z + z**2*(-12259.672610882488_dp)*Hr3(1,0,0) + &
                        (-9564.04938271605_dp)*Hr3(1,0,0) + (nf*(-285.49794238683126_dp)*Hr3(1,0,0))/z + &
                        nf*(-214.12345679012347_dp)*Hr3(1,0,0) + nf*z*(214.12345679012347_dp)*Hr3(1,0,0) + &
                        nf*z**2*(285.49794238683126_dp)*Hr3(1,0,0) + z*(8899.818930041152_dp)*Hr3(1,0,0) + &
                        ((12923.903063557385_dp)*Hr3(1,0,0))/z + (-12045.234567901234_dp)*Hr3(1,0,1) + &
                        z**2*(-10197.947873799725_dp)*Hr3(1,0,1) + (nf*(-159.25377229080934_dp)*Hr3(1,0,1))/z + &
                        nf*(-119.440329218107_dp)*Hr3(1,0,1) + nf*z*(119.440329218107_dp)*Hr3(1,0,1) + &
                        nf*z**2*(159.25377229080934_dp)*Hr3(1,0,1) + z*(10052.543209876543_dp)*Hr3(1,0,1) + &
                        ((12190.639231824416_dp)*Hr3(1,0,1))/z + (-12045.234567901234_dp)*Hr3(1,1,0) + &
                        z**2*(-10197.947873799725_dp)*Hr3(1,1,0) + (nf*(-159.25377229080934_dp)*Hr3(1,1,0))/z + &
                        nf*(-119.440329218107_dp)*Hr3(1,1,0) + nf*z*(119.440329218107_dp)*Hr3(1,1,0) + &
                        nf*z**2*(159.25377229080934_dp)*Hr3(1,1,0) + z*(10052.543209876543_dp)*Hr3(1,1,0) + &
                        ((12190.639231824416_dp)*Hr3(1,1,0))/z + ((-6745.371742112483_dp)*Hr3(1,1,1))/z + &
                        (-4152.493827160494_dp)*Hr3(1,1,1) + nf*z**2*(-491.9835390946502_dp)*Hr3(1,1,1) + &
                        nf*z*(-368.98765432098764_dp)*Hr3(1,1,1) + z*(167.11111111111111_dp)*Hr3(1,1,1) + &
                        nf*(368.98765432098764_dp)*Hr3(1,1,1) + (nf*(491.9835390946502_dp)*Hr3(1,1,1))/z + &
                        z**2*(10730.754458161866_dp)*Hr3(1,1,1) + (-2926.617283950617_dp)*Hr4(-1,-1,-1,0) + &
                        z*(-2926.617283950617_dp)*Hr4(-1,-1,-1,0) + ((3902.1563786008232_dp)*Hr4(-1,-1,-1,0))/z + &
                        z**2*(3902.1563786008232_dp)*Hr4(-1,-1,-1,0) + ((-3561.8765432098767_dp)*Hr4(-1,-1,0,0))/z + &
                        z**2*(-3561.8765432098767_dp)*Hr4(-1,-1,0,0) + (2671.4074074074074_dp)*Hr4(-1,-1,0,0) + &
                        z*(2671.4074074074074_dp)*Hr4(-1,-1,0,0) + ((-3221.59670781893_dp)*Hr4(-1,-1,0,1))/z + &
                        z**2*(-3221.59670781893_dp)*Hr4(-1,-1,0,1) + (2416.1975308641977_dp)*Hr4(-1,-1,0,1) + &
                        z*(2416.1975308641977_dp)*Hr4(-1,-1,0,1) + ((-1951.0781893004116_dp)*Hr4(-1,0,-1,0))/z + &
                        z**2*(-1951.0781893004116_dp)*Hr4(-1,0,-1,0) + (1463.3086419753085_dp)*Hr4(-1,0,-1,0) + &
                        z*(1463.3086419753085_dp)*Hr4(-1,0,-1,0) + (-1208.0987654320988_dp)*Hr4(-1,0,0,0) + &
                        z*(-1208.0987654320988_dp)*Hr4(-1,0,0,0) + ((1610.798353909465_dp)*Hr4(-1,0,0,0))/z + &
                        z**2*(1610.798353909465_dp)*Hr4(-1,0,0,0) + (-2160.9876543209875_dp)*Hr4(-1,0,0,1) + &
                        z*(-2160.9876543209875_dp)*Hr4(-1,0,0,1) + ((2881.3168724279835_dp)*Hr4(-1,0,0,1))/z + &
                        z**2*(2881.3168724279835_dp)*Hr4(-1,0,0,1) + (-952.8888888888889_dp)*Hr4(-1,0,1,0) + &
                        z*(-952.8888888888889_dp)*Hr4(-1,0,1,0) + ((1270.5185185185185_dp)*Hr4(-1,0,1,0))/z + &
                        z**2*(1270.5185185185185_dp)*Hr4(-1,0,1,0) + (-1905.7777777777778_dp)*Hr4(-1,0,1,1) + &
                        z*(-1905.7777777777778_dp)*Hr4(-1,0,1,1) + ((2541.037037037037_dp)*Hr4(-1,0,1,1))/z + &
                        z**2*(2541.037037037037_dp)*Hr4(-1,0,1,1) + z*(-7480.098765432099_dp)*Hr4(0,-1,-1,0) + &
                        ((-2705.382716049383_dp)*Hr4(0,-1,-1,0))/z + z**2*(-1196.7736625514403_dp)*Hr4(0,-1,-1,0) + &
                        nf*(-634.7325102880659_dp)*Hr4(0,-1,-1,0) + nf*z*(634.7325102880659_dp)*Hr4(0,-1,-1,0) + &
                        (10406.716049382716_dp)*Hr4(0,-1,-1,0) + (-9472.395061728395_dp)*Hr4(0,-1,0,0) + &
                        nf*z*(-578.1069958847737_dp)*Hr4(0,-1,0,0) + nf*(578.1069958847737_dp)*Hr4(0,-1,0,0) + &
                        z**2*(1099.851851851852_dp)*Hr4(0,-1,0,0) + ((2462.0246913580245_dp)*Hr4(0,-1,0,0))/z + &
                        z*(6800.9876543209875_dp)*Hr4(0,-1,0,0) + (-8538.074074074075_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(-521.4814814814815_dp)*Hr4(0,-1,0,1) + nf*(521.4814814814815_dp)*Hr4(0,-1,0,1) + &
                        z**2*(1002.9300411522634_dp)*Hr4(0,-1,0,1) + ((2218.6666666666665_dp)*Hr4(0,-1,0,1))/z + &
                        z*(6121.876543209876_dp)*Hr4(0,-1,0,1) + nf*z*(-779.5884773662551_dp)*Hr4(0,0,-1,0) + &
                        nf*(-144.8559670781893_dp)*Hr4(0,0,-1,0) + z**2*(269.69547325102883_dp)*Hr4(0,0,-1,0) + &
                        ((1024._dp)*Hr4(0,0,-1,0))/z + (4830.024691358025_dp)*Hr4(0,0,-1,0) + &
                        z*(13773.432098765432_dp)*Hr4(0,0,-1,0) + nf*z*(-1407.4732510288065_dp)*Hr4(0,0,0,0) + &
                        nf*(-711.1111111111111_dp)*Hr4(0,0,0,0) + z**2*(-382.06858710562415_dp)*Hr4(0,0,0,0) + &
                        nf**2*(31.604938271604937_dp)*Hr4(0,0,0,0) + nf**2*z*(31.604938271604937_dp)*Hr4(0,0,0,0) + &
                        (5021.497942386832_dp)*Hr4(0,0,0,0) + z*(16696.88888888889_dp)*Hr4(0,0,0,0) + &
                        (-9244.04938271605_dp)*Hr4(0,0,0,1) + z**2*(-1297.9094650205761_dp)*Hr4(0,0,0,1) + &
                        nf*z*(-997.1358024691358_dp)*Hr4(0,0,0,1) + nf**2*(-6.320987654320987_dp)*Hr4(0,0,0,1) + &
                        nf**2*z*(-6.320987654320987_dp)*Hr4(0,0,0,1) + nf*(277.3333333333333_dp)*Hr4(0,0,0,1) + &
                        ((768._dp)*Hr4(0,0,0,1))/z + z*(16070.32098765432_dp)*Hr4(0,0,0,1) + &
                        z**2*(-1942.6502057613168_dp)*Hr4(0,0,1,0) + nf*z*(-863.0781893004115_dp)*Hr4(0,0,1,0) + &
                        nf*(-437.201646090535_dp)*Hr4(0,0,1,0) + ((2360.8888888888887_dp)*Hr4(0,0,1,0))/z + &
                        (5542.716049382716_dp)*Hr4(0,0,1,0) + z*(11531.456790123457_dp)*Hr4(0,0,1,0) + &
                        z**2*(-3177.349794238683_dp)*Hr4(0,0,1,1) + z*(-1516.2469135802469_dp)*Hr4(0,0,1,1) + &
                        nf*(-1097.7448559670781_dp)*Hr4(0,0,1,1) + nf*z*(-287.0781893004115_dp)*Hr4(0,0,1,1) + &
                        ((2958.222222222222_dp)*Hr4(0,0,1,1))/z + (18444.64197530864_dp)*Hr4(0,0,1,1) + &
                        z**2*(-2083.467764060357_dp)*Hr4(0,1,0,0) + nf*(-428.24691358024694_dp)*Hr4(0,1,0,0) + &
                        nf*z*(-428.24691358024694_dp)*Hr4(0,1,0,0) + ((2812.8395061728397_dp)*Hr4(0,1,0,0))/z + &
                        z*(4173.432098765432_dp)*Hr4(0,1,0,0) + (7845.6625514403295_dp)*Hr4(0,1,0,0) + &
                        z**2*(-3810.5020576131687_dp)*Hr4(0,1,0,1) + nf*(-238.880658436214_dp)*Hr4(0,1,0,1) + &
                        nf*z*(-238.880658436214_dp)*Hr4(0,1,0,1) + z*(177.3827160493827_dp)*Hr4(0,1,0,1) + &
                        ((4326.716049382716_dp)*Hr4(0,1,0,1))/z + (6280.2962962962965_dp)*Hr4(0,1,0,1) + &
                        z**2*(-3810.5020576131687_dp)*Hr4(0,1,1,0) + nf*(-238.880658436214_dp)*Hr4(0,1,1,0) + &
                        nf*z*(-238.880658436214_dp)*Hr4(0,1,1,0) + z*(177.3827160493827_dp)*Hr4(0,1,1,0) + &
                        ((4326.716049382716_dp)*Hr4(0,1,1,0))/z + (6280.2962962962965_dp)*Hr4(0,1,1,0) + &
                        z*(-14025.481481481482_dp)*Hr4(0,1,1,1) + (-7499.061728395061_dp)*Hr4(0,1,1,1) + &
                        z**2*(-4890.3374485596705_dp)*Hr4(0,1,1,1) + nf*(737.9753086419753_dp)*Hr4(0,1,1,1) + &
                        nf*z*(737.9753086419753_dp)*Hr4(0,1,1,1) + ((3811.5555555555557_dp)*Hr4(0,1,1,1))/z + &
                        z**2*(-1552.8559670781892_dp)*Hr4(1,0,0,0) + z*(-1164.641975308642_dp)*Hr4(1,0,0,0) + &
                        (1164.641975308642_dp)*Hr4(1,0,0,0) + ((1552.8559670781892_dp)*Hr4(1,0,0,0))/z + &
                        z**2*(-3445.991769547325_dp)*Hr4(1,0,0,1) + z*(-2584.4938271604938_dp)*Hr4(1,0,0,1) + &
                        (2584.4938271604938_dp)*Hr4(1,0,0,1) + ((3445.991769547325_dp)*Hr4(1,0,0,1))/z + &
                        z**2*(-3786.2716049382716_dp)*Hr4(1,0,1,0) + z*(-2839.703703703704_dp)*Hr4(1,0,1,0) + &
                        (2839.703703703704_dp)*Hr4(1,0,1,0) + ((3786.2716049382716_dp)*Hr4(1,0,1,0))/z + &
                        z**2*(-3786.2716049382716_dp)*Hr4(1,0,1,1) + z*(-2839.703703703704_dp)*Hr4(1,0,1,1) + &
                        (2839.703703703704_dp)*Hr4(1,0,1,1) + ((3786.2716049382716_dp)*Hr4(1,0,1,1))/z + &
                        z**2*(-3445.991769547325_dp)*Hr4(1,1,0,0) + z*(-2584.4938271604938_dp)*Hr4(1,1,0,0) + &
                        (2584.4938271604938_dp)*Hr4(1,1,0,0) + ((3445.991769547325_dp)*Hr4(1,1,0,0))/z + &
                        z**2*(-3786.2716049382716_dp)*Hr4(1,1,0,1) + z*(-2839.703703703704_dp)*Hr4(1,1,0,1) + &
                        (2839.703703703704_dp)*Hr4(1,1,0,1) + ((3786.2716049382716_dp)*Hr4(1,1,0,1))/z + &
                        z**2*(-3786.2716049382716_dp)*Hr4(1,1,1,0) + z*(-2839.703703703704_dp)*Hr4(1,1,1,0) + &
                        (2839.703703703704_dp)*Hr4(1,1,1,0) + ((3786.2716049382716_dp)*Hr4(1,1,1,0))/z + &
                        (-5853.234567901234_dp)*Hr5(0,-1,-1,-1,0) + z*(5853.234567901234_dp)*Hr5(0,-1,-1,-1,0) + &
                        z*(-5342.814814814815_dp)*Hr5(0,-1,-1,0,0) + (5342.814814814815_dp)*Hr5(0,-1,-1,0,0) + &
                        z*(-4832.395061728395_dp)*Hr5(0,-1,-1,0,1) + (4832.395061728395_dp)*Hr5(0,-1,-1,0,1) + &
                        z*(-2926.617283950617_dp)*Hr5(0,-1,0,-1,0) + (2926.617283950617_dp)*Hr5(0,-1,0,-1,0) + &
                        (-2416.1975308641977_dp)*Hr5(0,-1,0,0,0) + z*(2416.1975308641977_dp)*Hr5(0,-1,0,0,0) + &
                        (-4321.975308641975_dp)*Hr5(0,-1,0,0,1) + z*(4321.975308641975_dp)*Hr5(0,-1,0,0,1) + &
                        (-1905.7777777777778_dp)*Hr5(0,-1,0,1,0) + z*(1905.7777777777778_dp)*Hr5(0,-1,0,1,0) + &
                        (-3811.5555555555557_dp)*Hr5(0,-1,0,1,1) + z*(3811.5555555555557_dp)*Hr5(0,-1,0,1,1) + &
                        z*(-9013.728395061727_dp)*Hr5(0,0,-1,-1,0) + (-3160.4938271604938_dp)*Hr5(0,0,-1,-1,0) + &
                        (2868.1481481481483_dp)*Hr5(0,0,-1,0,0) + z*(8210.962962962964_dp)*Hr5(0,0,-1,0,0) + &
                        (2575.8024691358023_dp)*Hr5(0,0,-1,0,1) + z*(7408.197530864198_dp)*Hr5(0,0,-1,0,1) + &
                        (-1144.0987654320988_dp)*Hr5(0,0,0,-1,0) + nf*z*(-113.77777777777777_dp)*Hr5(0,0,0,-1,0) + &
                        nf*(113.77777777777777_dp)*Hr5(0,0,0,-1,0) + z*(7231.20987654321_dp)*Hr5(0,0,0,-1,0) + &
                        (-1013.4650205761317_dp)*Hr5(0,0,0,0,0) + nf*z*(-417.18518518518516_dp)*Hr5(0,0,0,0,0) + &
                        nf*(151.7037037037037_dp)*Hr5(0,0,0,0,0) + z*(8309.991769547325_dp)*Hr5(0,0,0,0,0) + &
                        nf*(-297.08641975308643_dp)*Hr5(0,0,0,0,1) + nf*z*(-297.08641975308643_dp)*Hr5(0,0,0,0,1) + &
                        (3367.5061728395062_dp)*Hr5(0,0,0,0,1) + z*(15258.864197530864_dp)*Hr5(0,0,0,0,1) + &
                        (-2277.135802469136_dp)*Hr5(0,0,0,1,0) + nf*(-82.17283950617283_dp)*Hr5(0,0,0,1,0) + &
                        nf*z*(-82.17283950617283_dp)*Hr5(0,0,0,1,0) + z*(15116.641975308641_dp)*Hr5(0,0,0,1,0) + &
                        (-6479.0123456790125_dp)*Hr5(0,0,0,1,1) + z*(17414.320987654322_dp)*Hr5(0,0,0,1,1) + &
                        (-803.8189300411523_dp)*Hr5(0,0,1,0,0) + z*(8995.292181069959_dp)*Hr5(0,0,1,0,0) + &
                        (1926.320987654321_dp)*Hr5(0,0,1,0,1) + z*(12820.543209876543_dp)*Hr5(0,0,1,0,1) + &
                        (1926.320987654321_dp)*Hr5(0,0,1,1,0) + z*(12820.543209876543_dp)*Hr5(0,0,1,1,0) + &
                        (9067.456790123457_dp)*Hr5(0,0,1,1,1) + z*(9067.456790123457_dp)*Hr5(0,0,1,1,1) + &
                        (2329.283950617284_dp)*Hr5(0,1,0,0,0) + z*(2329.283950617284_dp)*Hr5(0,1,0,0,0) + &
                        (5168.9876543209875_dp)*Hr5(0,1,0,0,1) + z*(5168.9876543209875_dp)*Hr5(0,1,0,0,1) + &
                        (5679.407407407408_dp)*Hr5(0,1,0,1,0) + z*(5679.407407407408_dp)*Hr5(0,1,0,1,0) + &
                        (5679.407407407408_dp)*Hr5(0,1,0,1,1) + z*(5679.407407407408_dp)*Hr5(0,1,0,1,1) + &
                        (5168.9876543209875_dp)*Hr5(0,1,1,0,0) + z*(5168.9876543209875_dp)*Hr5(0,1,1,0,0) + &
                        (5679.407407407408_dp)*Hr5(0,1,1,0,1) + z*(5679.407407407408_dp)*Hr5(0,1,1,0,1) + &
                        (5679.407407407408_dp)*Hr5(0,1,1,1,0) + z*(5679.407407407408_dp)*Hr5(0,1,1,1,0)

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (5)
                select case (type_in)
                    case (1)
                        Ibar_select = (-17848.732953161365_dp) + z**2*(-10847.762647276682_dp) + (nf*(-909.4256397571146_dp))/z + &
                        nf**2*z*(-25.625255982016185_dp) + nf**2*z**2*(-24.97960676726109_dp) + nf**2*(-6.434035131536076_dp) + &
                        (nf**3*(-0.21069958847736625_dp))/z + nf**3*(-0.1580246913580247_dp) + nf**3*z*(0.1580246913580247_dp) + &
                        nf**3*z**2*(0.21069958847736625_dp) + (nf**2*(24.97960676726109_dp))/z + z*(68.1354758928428_dp) + &
                        nf*z*(251.23800496872934_dp) + nf*z**2*(952.7490061267798_dp) + nf*(1238.275162930763_dp) + &
                        (10166.453905323033_dp)/z + z**2*(-4610.9553510738615_dp)*Hr1(0) + nf*(-533.4027576349986_dp)*Hr1(0) + &
                        (nf*(-249.73168724279836_dp)*Hr1(0))/z + nf*z*(-214.89560452615643_dp)*Hr1(0) + &
                        nf**2*z**2*(-2.528395061728395_dp)*Hr1(0) + nf**3*(-0.3160493827160494_dp)*Hr1(0) + &
                        nf**3*z*(-0.3160493827160494_dp)*Hr1(0) + (nf**2*(3.9681755829903977_dp)*Hr1(0))/z + &
                        nf**2*z*(24.046090534979424_dp)*Hr1(0) + nf**2*(28.918518518518518_dp)*Hr1(0) + &
                        nf*z**2*(215.73882030178328_dp)*Hr1(0) + z*(1035.7540733617961_dp)*Hr1(0) + (2943.451757073187_dp)*Hr1(0) + &
                        ((3750.3420169119977_dp)*Hr1(0))/z + z**2*(-8111.754777696588_dp)*Hr1(1) + (-7269.978155410688_dp)*Hr1(1) + &
                        (nf*(-465.4705075445816_dp)*Hr1(1))/z + nf*z*(-63.60493827160494_dp)*Hr1(1) + &
                        nf**2*z**2*(-6.496570644718793_dp)*Hr1(1) + nf**2*z*(-4.872427983539095_dp)*Hr1(1) + &
                        nf**2*(4.872427983539095_dp)*Hr1(1) + (nf**2*(6.496570644718793_dp)*Hr1(1))/z + &
                        nf*(63.60493827160494_dp)*Hr1(1) + nf*z**2*(465.4705075445816_dp)*Hr1(1) + z*(7269.978155410688_dp)*Hr1(1) + &
                        ((8111.754777696588_dp)*Hr1(1))/z + (-1494.045094878264_dp)*Hr2(0,0) + z**2*(-538.8641975308642_dp)*Hr2(0,0) &
                        + nf*z*(-419.49359596692744_dp)*Hr2(0,0) + (nf*(-24.9679012345679_dp)*Hr2(0,0))/z + &
                        nf**2*(1.0008230452674898_dp)*Hr2(0,0) + nf**2*z*(8.7440329218107_dp)*Hr2(0,0) + &
                        nf*z**2*(11.799176954732511_dp)*Hr2(0,0) + nf*(102.22492255159105_dp)*Hr2(0,0) + ((704._dp)*Hr2(0,0))/z + &
                        z*(3257.2920145289104_dp)*Hr2(0,0) + z**2*(-1755.6543209876543_dp)*Hr2(0,1) + &
                        nf*(-421.1358024691358_dp)*Hr2(0,1) + nf*z*(-338.17283950617286_dp)*Hr2(0,1) + &
                        (nf*(-68.47736625514403_dp)*Hr2(0,1))/z + nf**2*(9.74485596707819_dp)*Hr2(0,1) + &
                        nf**2*z*(9.74485596707819_dp)*Hr2(0,1) + nf*z**2*(42.13991769547325_dp)*Hr2(0,1) + &
                        ((2085.925925925926_dp)*Hr2(0,1))/z + (2427.154800289737_dp)*Hr2(0,1) + z*(4447.368792059285_dp)*Hr2(0,1) + &
                        z**2*(-1920.79012345679_dp)*Hr2(1,0) + (-1010.1069958847737_dp)*Hr2(1,0) + &
                        (nf*(-55.30864197530864_dp)*Hr2(1,0))/z + nf*(-41.48148148148148_dp)*Hr2(1,0) + &
                        nf*z*(41.48148148148148_dp)*Hr2(1,0) + nf*z**2*(55.30864197530864_dp)*Hr2(1,0) + &
                        z*(1010.1069958847737_dp)*Hr2(1,0) + ((1920.79012345679_dp)*Hr2(1,0))/z + (z**2*(-3841.58024691358_dp) + &
                        (-2020.2139917695474_dp) + nf*((-110.61728395061728_dp)/z + (-82.96296296296296_dp) + &
                        z*(82.96296296296296_dp) + z**2*(110.61728395061728_dp)) + z*(2020.2139917695474_dp) + &
                        (3841.58024691358_dp)/z)*Hr2(1,1) + (nf*(z*(-174.45925925925926_dp) + (-58.31111111111111_dp)) + &
                        z**2*(-17.979698216735255_dp) + nf**2*((2.528395061728395_dp) + z*(2.528395061728395_dp)) + (51.2_dp)/z + &
                        (387.02880658436214_dp) + z*(2615.8353909465022_dp))*Hr3(0,0,0) + (-968.6913580246913_dp)*Hr3(0,0,1) + &
                        nf*z*(-199.11111111111111_dp)*Hr3(0,0,1) + z**2*(-84.2798353909465_dp)*Hr3(0,0,1) + &
                        nf*(33.18518518518518_dp)*Hr3(0,0,1) + ((199.11111111111111_dp)*Hr3(0,0,1))/z + &
                        z*(3380.1481481481483_dp)*Hr3(0,0,1) + z**2*(-160.13168724279836_dp)*Hr3(0,1,0) + &
                        nf*(-82.96296296296296_dp)*Hr3(0,1,0) + nf*z*(-82.96296296296296_dp)*Hr3(0,1,0) + &
                        ((274.962962962963_dp)*Hr3(0,1,0))/z + z*(1042.567901234568_dp)*Hr3(0,1,0) + &
                        (1368.888888888889_dp)*Hr3(0,1,0) + z**2*(-320.2633744855967_dp)*Hr3(0,1,1) + &
                        nf*(-165.92592592592592_dp)*Hr3(0,1,1) + nf*z*(-165.92592592592592_dp)*Hr3(0,1,1) + &
                        ((549.925925925926_dp)*Hr3(0,1,1))/z + z*(2085.135802469136_dp)*Hr3(0,1,1) + &
                        (2737.777777777778_dp)*Hr3(0,1,1) + z**2*(-145.0315500685871_dp)*Hr3(1,0,0) + &
                        z*(-108.77366255144032_dp)*Hr3(1,0,0) + (108.77366255144032_dp)*Hr3(1,0,0) + &
                        ((145.0315500685871_dp)*Hr3(1,0,0))/z + z**2*(-435.0946502057613_dp)*Hr3(1,0,1) + &
                        z*(-326.320987654321_dp)*Hr3(1,0,1) + (326.320987654321_dp)*Hr3(1,0,1) + &
                        ((435.0946502057613_dp)*Hr3(1,0,1))/z + z**2*(-435.0946502057613_dp)*Hr3(1,1,0) + &
                        z*(-326.320987654321_dp)*Hr3(1,1,0) + (326.320987654321_dp)*Hr3(1,1,0) + &
                        ((435.0946502057613_dp)*Hr3(1,1,0))/z + z**2*(-870.1893004115226_dp)*Hr3(1,1,1) + &
                        z*(-652.641975308642_dp)*Hr3(1,1,1) + (652.641975308642_dp)*Hr3(1,1,1) + &
                        ((870.1893004115226_dp)*Hr3(1,1,1))/z + (-50.67325102880658_dp)*Hr4(0,0,0,0) + &
                        nf*z*(-26.548148148148147_dp)*Hr4(0,0,0,0) + nf*(7.5851851851851855_dp)*Hr4(0,0,0,0) + &
                        z*(777.0600823045268_dp)*Hr4(0,0,0,0) + nf*(-18.962962962962962_dp)*Hr4(0,0,0,1) + &
                        nf*z*(-18.962962962962962_dp)*Hr4(0,0,0,1) + (216.49382716049382_dp)*Hr4(0,0,0,1) + &
                        z*(1453.8271604938273_dp)*Hr4(0,0,0,1) + (-292.34567901234567_dp)*Hr4(0,0,1,0) + &
                        z*(944.9876543209876_dp)*Hr4(0,0,1,0) + (-584.6913580246913_dp)*Hr4(0,0,1,1) + &
                        z*(1889.9753086419753_dp)*Hr4(0,0,1,1) + (217.54732510288065_dp)*Hr4(0,1,0,0) + &
                        z*(217.54732510288065_dp)*Hr4(0,1,0,0) + (652.641975308642_dp)*Hr4(0,1,0,1) + &
                        z*(652.641975308642_dp)*Hr4(0,1,0,1) + (652.641975308642_dp)*Hr4(0,1,1,0) + &
                        z*(652.641975308642_dp)*Hr4(0,1,1,0) + (1305.283950617284_dp)*Hr4(0,1,1,1) + &
                        z*(1305.283950617284_dp)*Hr4(0,1,1,1)

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (6)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case default
                Ibar_select = 0._dp
        end select
    case (6)
        select case (powLperp)
            case (2)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (3)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (4)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (5)
                select case (type_in)
                    case (1)
                        Ibar_select = (0._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case (6)
                select case (type_in)
                    case (1)
                        Ibar_select = (-276342.93048001593_dp) + z**2*(-150595.7294935212_dp) + (nf*(-19416.34072939799_dp))/z + &
                        z*(-7733.380949321348_dp) + nf**2*(-1010.7765588511406_dp) + nf**2*z**2*(-936.3362909285732_dp) + &
                        nf**2*z*(-546.1520765774793_dp) + (nf**3*(-17.225341665396535_dp))/z + nf**3*(-0.5995716142376857_dp) + &
                        nf**4*z**2*(-0.11705532693187014_dp) + nf**4*z*(-0.0877914951989026_dp) + nf**4*(0.0877914951989026_dp) + &
                        (nf**4*(0.11705532693187014_dp))/z + nf**3*z**2*(17.225341665396535_dp) + nf**3*z*(20.778632927752863_dp) + &
                        (nf**2*(894.0526853517799_dp))/z + nf*z*(1994.013293373795_dp) + nf*z**2*(21022.417179890606_dp) + &
                        nf*(31542.89309669743_dp) + (148537.61617422395_dp)/z + z**2*(-93646.69234958959_dp)*Hr1(0) + &
                        nf*(-7842.566922527337_dp)*Hr1(0) + (nf*(-6684.803817022826_dp)*Hr1(0))/z + &
                        nf*z*(-1737.1800398309024_dp)*Hr1(0) + nf**2*z**2*(-179.62139917695472_dp)*Hr1(0) + &
                        nf**3*(-22.28148148148148_dp)*Hr1(0) + nf**3*z*(-19.214631915866484_dp)*Hr1(0) + &
                        (nf**3*(-2.5908245694253926_dp)*Hr1(0))/z + nf**4*(0.1755829903978052_dp)*Hr1(0) + &
                        nf**4*z*(0.1755829903978052_dp)*Hr1(0) + nf**3*z**2*(1.4983081847279378_dp)*Hr1(0) + &
                        (nf**2*(231.7266270385612_dp)*Hr1(0))/z + nf**2*z*(332.04124367589344_dp)*Hr1(0) + &
                        nf**2*(647.621870370288_dp)*Hr1(0) + z*(1681.3043684622971_dp)*Hr1(0) + &
                        nf*z**2*(7280.491409874511_dp)*Hr1(0) + (18837.682745143822_dp)*Hr1(0) + ((61963.749240709236_dp)*Hr1(0))/z &
                        + z**2*(-144109.62722332927_dp)*Hr1(1) + (-139814.9133117352_dp)*Hr1(1) + &
                        (nf*(-13651.672564039956_dp)*Hr1(1))/z + nf*z*(-8947.919681527264_dp)*Hr1(1) + &
                        nf**2*z**2*(-411.3480262155159_dp)*Hr1(1) + nf**2*z*(-65.69144947416552_dp)*Hr1(1) + &
                        (nf**3*(-4.08913275415333_dp)*Hr1(1))/z + nf**3*(-3.066849565614998_dp)*Hr1(1) + &
                        nf**3*z*(3.066849565614998_dp)*Hr1(1) + nf**3*z**2*(4.08913275415333_dp)*Hr1(1) + &
                        nf**2*(65.69144947416552_dp)*Hr1(1) + (nf**2*(411.3480262155159_dp)*Hr1(1))/z + &
                        nf*(8947.919681527264_dp)*Hr1(1) + nf*z**2*(13651.672564039956_dp)*Hr1(1) + z*(139814.91331173517_dp)*Hr1(1) &
                        + ((144109.62722332927_dp)*Hr1(1))/z + (-28551.191250098786_dp)*Hr2(0,0) + &
                        z**2*(-16020.392161490992_dp)*Hr2(0,0) + nf*z*(-6704.01805496849_dp)*Hr2(0,0) + &
                        (nf*(-1161.1654320987655_dp)*Hr2(0,0))/z + nf**2*(-74.59927600935261_dp)*Hr2(0,0) + &
                        nf**2*z**2*(-8.67770156988264_dp)*Hr2(0,0) + nf**3*z*(-4.82267946959305_dp)*Hr2(0,0) + &
                        nf**3*(-1.3110196616369456_dp)*Hr2(0,0) + (nf**2*(21.530376467001982_dp)*Hr2(0,0))/z + &
                        nf**2*z*(396.33186253660074_dp)*Hr2(0,0) + nf*z**2*(742.3921963115379_dp)*Hr2(0,0) + &
                        nf*(3004.482012537933_dp)*Hr2(0,0) + ((15266.820563156514_dp)*Hr2(0,0))/z + &
                        z*(33818.017799142515_dp)*Hr2(0,0) + z**2*(-47763.54896655856_dp)*Hr2(0,1) + &
                        nf*z*(-7507.628644871091_dp)*Hr2(0,1) + nf*(-6432.299884011465_dp)*Hr2(0,1) + &
                        (nf*(-3296.1492455418384_dp)*Hr2(0,1))/z + nf**2*z**2*(-31.277183356195703_dp)*Hr2(0,1) + &
                        nf**3*(-6.133699131229996_dp)*Hr2(0,1) + nf**3*z*(-6.133699131229996_dp)*Hr2(0,1) + &
                        (nf**2*(56.982533150434385_dp)*Hr2(0,1))/z + nf**2*z*(374.58875171467764_dp)*Hr2(0,1) + &
                        nf**2*(440.7835390946502_dp)*Hr2(0,1) + nf*z**2*(2458.602773967383_dp)*Hr2(0,1) + &
                        (9718.71601512541_dp)*Hr2(0,1) + ((45351.788206984995_dp)*Hr2(0,1))/z + z*(88488.53046464309_dp)*Hr2(0,1) + &
                        z**2*(-47577.66108802618_dp)*Hr2(1,0) + (-38619.912848818036_dp)*Hr2(1,0) + &
                        (nf*(-2877.3760097546105_dp)*Hr2(1,0))/z + nf*z*(-537.6643804298126_dp)*Hr2(1,0) + &
                        nf**2*z**2*(-44.12985825331504_dp)*Hr2(1,0) + nf**2*z*(-33.09739368998628_dp)*Hr2(1,0) + &
                        nf**2*(33.09739368998628_dp)*Hr2(1,0) + (nf**2*(44.12985825331504_dp)*Hr2(1,0))/z + &
                        nf*(537.6643804298126_dp)*Hr2(1,0) + nf*z**2*(2877.3760097546105_dp)*Hr2(1,0) + &
                        z*(38619.912848818036_dp)*Hr2(1,0) + ((47577.66108802618_dp)*Hr2(1,0))/z + &
                        ((-88.25971650663008_dp)*((-1055.015140078642_dp) + z**2*(-892.477537513963_dp) + nf**2*((-1._dp) + &
                        z*(-0.75_dp) + z**2*(0.75_dp) + z**3*(1._dp)) + nf*(z**3*(-65.2024756852343_dp) + z*(-12.18368700265252_dp) &
                        + z**2*(12.18368700265252_dp) + (65.2024756852343_dp)) + z*(892.477537513963_dp) + &
                        z**3*(1055.015140078642_dp))*Hr2(1,1))/z + (z**2*(-1059.1868312757201_dp) + nf**3*((-2.0367626886145405_dp) &
                        + z*(-2.0367626886145405_dp)) + nf*(z*(-4637.203574531355_dp) + (-1315.0771034189374_dp) + &
                        (-77.74814814814815_dp)/z + z**2*(21.975186709343088_dp)) + nf**2*((83.89355281207133_dp) + &
                        z*(159.85075445816187_dp)) + (1971.2_dp)/z + (6553.184664420731_dp) + z*(36256.66780606679_dp))*Hr3(0,0,0) + &
                        (-12617.096793453718_dp)*Hr3(0,0,1) + nf*z*(-6109.866666666667_dp)*Hr3(0,0,1) + &
                        z**2*(-4589.037037037037_dp)*Hr3(0,0,1) + (nf*(-292.240329218107_dp)*Hr3(0,0,1))/z + &
                        nf**2*(-9.76241426611797_dp)*Hr3(0,0,1) + nf*z**2*(102.25953360768176_dp)*Hr3(0,0,1) + &
                        nf**2*z*(142.1519890260631_dp)*Hr3(0,0,1) + nf*(1361.3124828532236_dp)*Hr3(0,0,1) + &
                        ((7770.074074074074_dp)*Hr3(0,0,1))/z + z*(57029.789796567406_dp)*Hr3(0,0,1) + &
                        z**2*(-7856.9876543209875_dp)*Hr3(0,1,0) + nf*(-2593.7119341563784_dp)*Hr3(0,1,0) + &
                        nf*z*(-2154.8422496570643_dp)*Hr3(0,1,0) + (nf*(-387.57018747142206_dp)*Hr3(0,1,0))/z + &
                        nf**2*(66.19478737997257_dp)*Hr3(0,1,0) + nf**2*z*(66.19478737997257_dp)*Hr3(0,1,0) + &
                        nf*z**2*(197.5893918609968_dp)*Hr3(0,1,0) + ((11038.024691358025_dp)*Hr3(0,1,0))/z + &
                        (19399.998719373518_dp)*Hr3(0,1,0) + z*(28072.671787503376_dp)*Hr3(0,1,0) + &
                        z**2*(-15713.975308641975_dp)*Hr3(0,1,1) + nf*(-5187.423868312757_dp)*Hr3(0,1,1) + &
                        nf*z*(-4309.6844993141285_dp)*Hr3(0,1,1) + (nf*(-775.1403749428441_dp)*Hr3(0,1,1))/z + &
                        nf**2*(132.38957475994513_dp)*Hr3(0,1,1) + nf**2*z*(132.38957475994513_dp)*Hr3(0,1,1) + &
                        nf*z**2*(395.1787837219936_dp)*Hr3(0,1,1) + ((22076.04938271605_dp)*Hr3(0,1,1))/z + &
                        (35740.01993498384_dp)*Hr3(0,1,1) + z*(53085.36607124355_dp)*Hr3(0,1,1) + &
                        z**2*(-6298.3374485596705_dp)*Hr3(1,0,0) + (-2890.8910227099527_dp)*Hr3(1,0,0) + &
                        (nf*(-195.0531931108063_dp)*Hr3(1,0,0))/z + nf*(-146.2898948331047_dp)*Hr3(1,0,0) + &
                        nf*z*(146.2898948331047_dp)*Hr3(1,0,0) + nf*z**2*(195.0531931108063_dp)*Hr3(1,0,0) + &
                        z*(2890.8910227099527_dp)*Hr3(1,0,0) + ((6298.3374485596705_dp)*Hr3(1,0,0))/z + &
                        z**2*(-18895.012345679013_dp)*Hr3(1,0,1) + (-8672.673068129858_dp)*Hr3(1,0,1) + &
                        (nf*(-585.1595793324188_dp)*Hr3(1,0,1))/z + nf*(-438.86968449931413_dp)*Hr3(1,0,1) + &
                        nf*z*(438.86968449931413_dp)*Hr3(1,0,1) + nf*z**2*(585.1595793324188_dp)*Hr3(1,0,1) + &
                        z*(8672.673068129858_dp)*Hr3(1,0,1) + ((18895.012345679013_dp)*Hr3(1,0,1))/z + &
                        z**2*(-18895.012345679013_dp)*Hr3(1,1,0) + (-8672.673068129858_dp)*Hr3(1,1,0) + &
                        (nf*(-585.1595793324188_dp)*Hr3(1,1,0))/z + nf*(-438.86968449931413_dp)*Hr3(1,1,0) + &
                        nf*z*(438.86968449931413_dp)*Hr3(1,1,0) + nf*z**2*(585.1595793324188_dp)*Hr3(1,1,0) + &
                        z*(8672.673068129858_dp)*Hr3(1,1,0) + ((18895.012345679013_dp)*Hr3(1,1,0))/z + &
                        z**2*(-37790.02469135803_dp)*Hr3(1,1,1) + (-17345.346136259715_dp)*Hr3(1,1,1) + &
                        (nf*(-1170.3191586648377_dp)*Hr3(1,1,1))/z + nf*(-877.7393689986283_dp)*Hr3(1,1,1) + &
                        nf*z*(877.7393689986283_dp)*Hr3(1,1,1) + nf*z**2*(1170.3191586648377_dp)*Hr3(1,1,1) + &
                        z*(17345.346136259715_dp)*Hr3(1,1,1) + ((37790.02469135803_dp)*Hr3(1,1,1))/z + &
                        nf*z*(-1449.1683584819386_dp)*Hr4(0,0,0,0) + (-1137.3388203017832_dp)*Hr4(0,0,0,0) + &
                        z**2*(-19.97744246303917_dp)*Hr4(0,0,0,0) + nf**2*(-8.427983539094651_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(34.55473251028807_dp)*Hr4(0,0,0,0) + ((102.4_dp)*Hr4(0,0,0,0))/z + &
                        nf*(219.46703246456332_dp)*Hr4(0,0,0,0) + z*(16995.965249199817_dp)*Hr4(0,0,0,0) + &
                        nf*z*(-2092.1064471879286_dp)*Hr4(0,0,0,1) + nf*(-659.8759945130315_dp)*Hr4(0,0,0,1) + &
                        z**2*(-115.86916628562719_dp)*Hr4(0,0,0,1) + nf**2*(26.126748971193415_dp)*Hr4(0,0,0,1) + &
                        nf**2*z*(26.126748971193415_dp)*Hr4(0,0,0,1) + ((500.6222222222222_dp)*Hr4(0,0,0,1))/z + &
                        (4871.403749428441_dp)*Hr4(0,0,0,1) + z*(31208.354823959762_dp)*Hr4(0,0,0,1) + &
                        (-6392.186556927298_dp)*Hr4(0,0,1,0) + nf*z*(-1154.9849108367628_dp)*Hr4(0,0,1,0) + &
                        z**2*(-288.424325560128_dp)*Hr4(0,0,1,0) + nf*(277.24554183813444_dp)*Hr4(0,0,1,0) + &
                        ((948.1481481481482_dp)*Hr4(0,0,1,0))/z + z*(19479.703703703704_dp)*Hr4(0,0,1,0) + &
                        (-12784.373113854595_dp)*Hr4(0,0,1,1) + nf*z*(-2309.9698216735255_dp)*Hr4(0,0,1,1) + &
                        z**2*(-576.848651120256_dp)*Hr4(0,0,1,1) + nf*(554.4910836762689_dp)*Hr4(0,0,1,1) + &
                        ((1896.2962962962963_dp)*Hr4(0,0,1,1))/z + z*(38959.40740740741_dp)*Hr4(0,0,1,1) + &
                        z**2*(-400.1731443377534_dp)*Hr4(0,1,0,0) + nf*(-292.5797896662094_dp)*Hr4(0,1,0,0) + &
                        nf*z*(-292.5797896662094_dp)*Hr4(0,1,0,0) + ((839.9890260631001_dp)*Hr4(0,1,0,0))/z + &
                        z*(3897.444901691815_dp)*Hr4(0,1,0,0) + (4827.566529492456_dp)*Hr4(0,1,0,0) + &
                        z**2*(-1200.5194330132601_dp)*Hr4(0,1,0,1) + nf*(-877.7393689986283_dp)*Hr4(0,1,0,1) + &
                        nf*z*(-877.7393689986283_dp)*Hr4(0,1,0,1) + ((2519.9670781893005_dp)*Hr4(0,1,0,1))/z + &
                        z*(11692.334705075446_dp)*Hr4(0,1,0,1) + (14482.699588477366_dp)*Hr4(0,1,0,1) + &
                        z**2*(-1200.5194330132601_dp)*Hr4(0,1,1,0) + nf*(-877.7393689986283_dp)*Hr4(0,1,1,0) + &
                        nf*z*(-877.7393689986283_dp)*Hr4(0,1,1,0) + ((2519.9670781893005_dp)*Hr4(0,1,1,0))/z + &
                        z*(11692.334705075446_dp)*Hr4(0,1,1,0) + (14482.699588477366_dp)*Hr4(0,1,1,0) + &
                        z**2*(-2401.0388660265203_dp)*Hr4(0,1,1,1) + nf*(-1755.4787379972565_dp)*Hr4(0,1,1,1) + &
                        nf*z*(-1755.4787379972565_dp)*Hr4(0,1,1,1) + ((5039.934156378601_dp)*Hr4(0,1,1,1))/z + &
                        z*(23384.66941015089_dp)*Hr4(0,1,1,1) + (28965.399176954732_dp)*Hr4(0,1,1,1) + &
                        z**2*(-310.04054260021337_dp)*Hr4(1,0,0,0) + z*(-232.53040695016003_dp)*Hr4(1,0,0,0) + &
                        (232.53040695016003_dp)*Hr4(1,0,0,0) + ((310.04054260021337_dp)*Hr4(1,0,0,0))/z + &
                        z**2*(-1240.1621704008535_dp)*Hr4(1,0,0,1) + z*(-930.1216278006401_dp)*Hr4(1,0,0,1) + &
                        (930.1216278006401_dp)*Hr4(1,0,0,1) + ((1240.1621704008535_dp)*Hr4(1,0,0,1))/z + &
                        z**2*(-1860.2432556012802_dp)*Hr4(1,0,1,0) + z*(-1395.1824417009602_dp)*Hr4(1,0,1,0) + &
                        (1395.1824417009602_dp)*Hr4(1,0,1,0) + ((1860.2432556012802_dp)*Hr4(1,0,1,0))/z + &
                        z**2*(-3720.4865112025605_dp)*Hr4(1,0,1,1) + z*(-2790.3648834019205_dp)*Hr4(1,0,1,1) + &
                        (2790.3648834019205_dp)*Hr4(1,0,1,1) + ((3720.4865112025605_dp)*Hr4(1,0,1,1))/z + &
                        z**2*(-1240.1621704008535_dp)*Hr4(1,1,0,0) + z*(-930.1216278006401_dp)*Hr4(1,1,0,0) + &
                        (930.1216278006401_dp)*Hr4(1,1,0,0) + ((1240.1621704008535_dp)*Hr4(1,1,0,0))/z + &
                        z**2*(-3720.4865112025605_dp)*Hr4(1,1,0,1) + z*(-2790.3648834019205_dp)*Hr4(1,1,0,1) + &
                        (2790.3648834019205_dp)*Hr4(1,1,0,1) + ((3720.4865112025605_dp)*Hr4(1,1,0,1))/z + &
                        z**2*(-3720.4865112025605_dp)*Hr4(1,1,1,0) + z*(-2790.3648834019205_dp)*Hr4(1,1,1,0) + &
                        (2790.3648834019205_dp)*Hr4(1,1,1,0) + ((3720.4865112025605_dp)*Hr4(1,1,1,0))/z + &
                        z**2*(-7440.973022405121_dp)*Hr4(1,1,1,1) + z*(-5580.729766803841_dp)*Hr4(1,1,1,1) + &
                        (5580.729766803841_dp)*Hr4(1,1,1,1) + ((7440.973022405121_dp)*Hr4(1,1,1,1))/z + &
                        nf*z*(-170.24526748971192_dp)*Hr5(0,0,0,0,0) + nf*(-22.334156378600824_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*(0.842798353909465_dp)*Hr5(0,0,0,0,0) + nf**2*z*(0.842798353909465_dp)*Hr5(0,0,0,0,0) + &
                        (103.2193872885231_dp)*Hr5(0,0,0,0,0) + z*(3110.113214449017_dp)*Hr5(0,0,0,0,0) + &
                        (-519.3510745313215_dp)*Hr5(0,0,0,0,1) + nf*z*(-225.86995884773663_dp)*Hr5(0,0,0,0,1) + &
                        nf*(69.9522633744856_dp)*Hr5(0,0,0,0,1) + z*(7384.411888431641_dp)*Hr5(0,0,0,0,1) + &
                        nf*(-77.95884773662551_dp)*Hr5(0,0,0,1,0) + nf*z*(-77.95884773662551_dp)*Hr5(0,0,0,1,0) + &
                        (1062.628257887517_dp)*Hr5(0,0,0,1,0) + z*(6732.554183813443_dp)*Hr5(0,0,0,1,0) + &
                        nf*(-155.91769547325103_dp)*Hr5(0,0,0,1,1) + nf*z*(-155.91769547325103_dp)*Hr5(0,0,0,1,1) + &
                        (2125.256515775034_dp)*Hr5(0,0,0,1,1) + z*(13465.108367626886_dp)*Hr5(0,0,0,1,1) + &
                        (-959.8536808413352_dp)*Hr5(0,0,1,0,0) + z*(2820.0969364426155_dp)*Hr5(0,0,1,0,0) + &
                        (-2879.5610425240056_dp)*Hr5(0,0,1,0,1) + z*(8460.290809327846_dp)*Hr5(0,0,1,0,1) + &
                        (-2879.5610425240056_dp)*Hr5(0,0,1,1,0) + z*(8460.290809327846_dp)*Hr5(0,0,1,1,0) + &
                        (-5759.122085048011_dp)*Hr5(0,0,1,1,1) + z*(16920.58161865569_dp)*Hr5(0,0,1,1,1) + &
                        (465.06081390032006_dp)*Hr5(0,1,0,0,0) + z*(465.06081390032006_dp)*Hr5(0,1,0,0,0) + &
                        (1860.2432556012802_dp)*Hr5(0,1,0,0,1) + z*(1860.2432556012802_dp)*Hr5(0,1,0,0,1) + &
                        (2790.3648834019205_dp)*Hr5(0,1,0,1,0) + z*(2790.3648834019205_dp)*Hr5(0,1,0,1,0) + &
                        (5580.729766803841_dp)*Hr5(0,1,0,1,1) + z*(5580.729766803841_dp)*Hr5(0,1,0,1,1) + &
                        (1860.2432556012802_dp)*Hr5(0,1,1,0,0) + z*(1860.2432556012802_dp)*Hr5(0,1,1,0,0) + &
                        (5580.729766803841_dp)*Hr5(0,1,1,0,1) + z*(5580.729766803841_dp)*Hr5(0,1,1,0,1) + &
                        (5580.729766803841_dp)*Hr5(0,1,1,1,0) + z*(5580.729766803841_dp)*Hr5(0,1,1,1,0) + &
                        (11161.459533607682_dp)*Hr5(0,1,1,1,1) + z*(11161.459533607682_dp)*Hr5(0,1,1,1,1)

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case default
                        Ibar_select = 0._dp
                end select
            case default
                Ibar_select = 0._dp
        end select
    case default
        Ibar_select = 0._dp
end select
