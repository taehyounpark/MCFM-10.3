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
                        Ibar_select = z*(2.6666666666666665_dp)
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
                        Ibar_select = (-5.333333333333333_dp)/z + z*(-2.6666666666666665_dp) + (5.333333333333334_dp)
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
                        Ibar_select = (z**3*(-86.2200328407225_dp) + z**5*(-72.93175074183976_dp) + (-44.34762556464744_dp) + &
                        z*(-43.76700883311665_dp) + z**7*(-1.2441960447119518_dp) + z**6*(7.84991423670669_dp) + &
                        Nf*(z**3*(-25.001344086021508_dp) + z*(-16.98765432098765_dp) + z**5*(-0.7109533468559838_dp) + &
                        z**6*(-0.06066734074823053_dp) + z**7*(-0.0201765447667087_dp) + z**4*(10.77611140867702_dp) + &
                        (11.061728395061728_dp) + z**2*(23.51085568326948_dp)) + z**2*(104.68828446318862_dp) + &
                        z**4*(124.82050788764909_dp))/z + ((-66.66666666666667_dp) + z*(-39.111111111111114_dp) + &
                        z**2*(2.688618925831202_dp) + z**3*(67.72398685651697_dp))*Log(z) + (z**2*(0.1952662721893491_dp) + &
                        z*(9.333333333333334_dp) + z**3*(11.018819503849445_dp) + (12.444444444444445_dp))*Log(z)**2 + &
                        ((-4.148148148148148_dp) + z*(-3.259259259259259_dp) + z**2*(0.001443001443001443_dp) + &
                        z**3*(2.7276657060518734_dp))*Log(z)**3 + (z*(-174.33771036747885_dp) + z**3*(-101.27627627627628_dp) + &
                        (-4.802616287732057_dp) + Nf*(z*(-28.970198612185694_dp) + z**3*(-3.6768968456947997_dp) + &
                        (16.202206438801_dp) + z**2*(17.630074204264677_dp)) + z**2*(259.9721584870427_dp))*Log(z*(-1._dp) + &
                        (1._dp)) + (z**2*(-106.35646527385659_dp) + (-79.07451533973273_dp) + Nf*(z*(-5.521815574158203_dp) + &
                        z**3*(-0.31933701657458563_dp) + z**2*(2.7954688674965356_dp) + (3.934572612125142_dp)) + &
                        z**3*(26.974308300395258_dp) + z*(153.56778342430516_dp))*Log(z*(-1._dp) + (1._dp))**2 + &
                        (z*(-19.089851828036377_dp) + z**3*(-8.67537122375832_dp) + (2.985018079916806_dp) + &
                        z**2*(23.298723490396412_dp))*Log(z*(-1._dp) + (1._dp))**3

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
                        Ibar_select = (z**2*(-334.49068211099967_dp) + z**4*(-178.48953068592058_dp) + z*(-16.03917053805621_dp) + &
                        z**5*(-7.902693310165074_dp) + z**7*(0.21675531914893617_dp) + z**6*(2.719227674979887_dp) + &
                        (10.318945069571619_dp) + Nf*(z**3*(-24.312534209085932_dp) + z*(-18.962962962962962_dp) + &
                        z**5*(-0.8577745025792188_dp) + z**7*(-0.08584686774941995_dp) + z**6*(-0.030690537084398978_dp) + &
                        z**4*(10.164082687338501_dp) + (11.851851851851851_dp) + z**2*(29.937580437580444_dp)) + &
                        z**3*(430.70112079701113_dp))/z + (z**2*(0.3026389668725435_dp) + z**3*(10.917067307692308_dp) + &
                        (49.77777777777778_dp) + z*(62.22222222222222_dp))*Log(z) + ((-24.88888888888889_dp) + &
                        z*(-19.555555555555557_dp) + z**2*(0.03954802259887006_dp) + z**3*(1.697350069735007_dp))*Log(z)**2 + &
                        (z**2*(-395.10632248462554_dp) + (-325.42735156739957_dp) + Nf*(z*(-31.278707683663434_dp) + &
                        z**3*(-4.588987217305801_dp) + (18.67819356651215_dp) + z**2*(20.74505689001264_dp)) + &
                        z**3*(88.17281879194631_dp) + z*(563.0275219267455_dp))*Log(z*(-1._dp) + (1._dp)) + ((-28.05023086242759_dp) &
                        + z**3*(-8.163377192982455_dp) + z**2*(6.0543013013151_dp) + z*(21.27041786520606_dp))*Log(z*(-1._dp) + &
                        (1._dp))**2

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
                        Ibar_select = z**2*(-151.83575418994414_dp) + (-82.66666666666669_dp)/z + z**4*(-9.073199527744983_dp) + &
                        z**6*(-1.1011826544021024_dp) + z**5*(1.3794506612410986_dp) + (24.8888888888889_dp) + &
                        z**3*(70.8823244552058_dp) + z*(152.8594733461785_dp) + (z*(-35.55555555555556_dp) + (-32._dp)/z + &
                        (-24.88888888888889_dp) + z**2*(0.15075707702435814_dp) + z**3*(1.5782652043868395_dp))*Log(z) + &
                        (z*(-209.44223210403788_dp) + z**3*(-31.727268369841344_dp) + (123.52414853376492_dp) + &
                        z**2*(140.7564630512254_dp))*Log(z*(-1._dp) + (1._dp))

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
                        Ibar_select = (z**5*(-94511.77090261626_dp) + z*(-7169.740093699399_dp) + z**7*(-342.0538838532933_dp) + &
                        Nf**2*(z**2*(-31.173181882262902_dp) + (-27.279733843564706_dp) + z**3*(-26.247694334650856_dp) + &
                        z**7*(0.02021772939346812_dp) + z**6*(0.1076417419884963_dp) + z**5*(0.8288920056100982_dp) + &
                        z**4*(21.88598781549173_dp) + z*(47.82294346360314_dp)) + z**6*(3406.906761190654_dp) + (9340.6870238563_dp) &
                        + z**2*(20999.973574346088_dp) + z**3*(30566.664811476672_dp) + Nf*(z**5*(-37315.55468746373_dp) + &
                        z**2*(-29146.324885331553_dp) + z**7*(-147.50397072375537_dp) + z*(-65.96810489118803_dp) + &
                        (244.58197514129375_dp) + z**6*(1148.2378807116193_dp) + z**3*(31138.04110829745_dp) + &
                        z**4*(34514.804302991746_dp)) + z**4*(35210.52952516608_dp))/z + (z*(478.2703489846124_dp) + &
                        (2860.5342887466486_dp)/z + (10576.670708024894_dp) + z**2*(20646.511902525468_dp) + &
                        Nf*(z**2*(-36449.28966223132_dp) + (-206.82543760377865_dp) + (66.85499822000071_dp)/z + &
                        z*(333.0650805773116_dp) + z**3*(43827.543301759135_dp)) + z**3*(151301.8124360537_dp))*Log(z) + &
                        ((-1170.7587186262226_dp) + z**2*(595.7540407497181_dp) + (923.1797016265684_dp)/z + &
                        z*(1376.8557494240706_dp) + Nf*(z**2*(-11432.693452380952_dp) + (-104.95370601500125_dp) + &
                        z*(-8.308237388561654_dp) + z**3*(1416.7368821292775_dp)) + z**3*(6192.450254632087_dp))*Log(z)**2 + &
                        (z**2*(-244.7558051372226_dp) + z*(-153.21994034339727_dp) + (911.287682068523_dp) + &
                        Nf*(z**2*(-1575.987044534413_dp) + (-39.46227709190672_dp) + z*(25.29492455418381_dp) + &
                        z**3*(8704.259040590407_dp)) + z**3*(16390.839774993303_dp))*Log(z)**3 + ((-58.76543209876543_dp) + &
                        z**2*(-25.033492822966508_dp) + z*(3.720164609053498_dp) + Nf*(z**2*(-105.47199533255542_dp) + &
                        (-2.831275720164609_dp) + z*(1.1193415637860082_dp) + z**3*(10.557038834951456_dp)) + &
                        z**3*(710.5834378920954_dp))*Log(z)**4 + (z*(-13.948971193415638_dp) + z**2*(-0.7566742944317315_dp) + &
                        (10.587654320987655_dp) + Nf*(z**2*(-2.85130890052356_dp) + (-0.9481481481481482_dp) + &
                        z*(0.4740740740740741_dp) + z**3*(208.85885486018643_dp)) + z**3*(304.8141745894555_dp))*Log(z)**5 + &
                        (z*(-83287.57934917185_dp) + z**3*(-31559.24995948433_dp) + Nf**2*(z*(-62.77291538634297_dp) + &
                        z**3*(-27.913618677042802_dp) + (3.009895268230336_dp) + z**2*(79.38034249885915_dp)) + &
                        (22452.174841868346_dp) + Nf*(z**2*(-64093.50016612267_dp) + (-28875.27347334142_dp) + &
                        z**3*(17191.63440541754_dp) + z*(76164.02824730723_dp)) + z**2*(90454.92088371639_dp))*Log(z*(-1._dp) + &
                        (1._dp)) + ((-6113.997999590717_dp) + z**2*(-4353.898241035487_dp) + Nf**2*(z**2*(-39.295505999708794_dp) + &
                        (-28.17249743907884_dp) + z**3*(10.49977588525325_dp) + z*(52.622548541188706_dp)) + &
                        z**3*(276.44338194242954_dp) + z*(9371.97012076143_dp) + Nf*(z*(-21862.470882213318_dp) + &
                        z**3*(-8352.664173860585_dp) + (6974.898594789018_dp) + z**2*(23363.031642269452_dp)))*Log(z*(-1._dp) + &
                        (1._dp))**2 + (z*(-8481.811550085245_dp) + z**3*(-3621.1800041306424_dp) + Nf**2*((-2.5620097169717466_dp) + &
                        z**3*(-1.2904073587385019_dp) + z*(0.8706491122420283_dp) + z**2*(1.9941136424805657_dp)) + &
                        Nf*(z**2*(-1671.4604442749549_dp) + (-523.4906387395719_dp) + z**3*(535.6529411764706_dp) + &
                        z*(1687.6547947873016_dp)) + (2013.5752176738144_dp) + z**2*(9886.596158383645_dp))*Log(z*(-1._dp) + &
                        (1._dp))**3 + (z**3*(-121.32318710832587_dp) + (-81.98127527815578_dp) + z*(-76.90796165489496_dp) + &
                        z**2*(249.6362923541338_dp) + Nf*(z*(-597.0222306164422_dp) + z**3*(-197.15028901734104_dp) + &
                        (204.2569584538716_dp) + z**2*(592.7962196161257_dp)))*Log(z*(-1._dp) + (1._dp))**4 + &
                        (z*(-255.22204900946286_dp) + z**3*(-89.44044502617801_dp) + (79.18709828793872_dp) + &
                        z**2*(263.006259945233_dp))*Log(z*(-1._dp) + (1._dp))**5

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
                        Ibar_select = (z**3*(-66939.26451567911_dp) + z**5*(-43640.85076677764_dp) + z**2*(-40820.90029940498_dp) + &
                        (-8329.890560190655_dp) + z**7*(-302.56850852711375_dp) + Nf**2*(z**2*(-77.97445672763347_dp) + &
                        z**4*(-35.920374470295265_dp) + (-12.378600823045268_dp) + z**7*(0.06725397943196539_dp) + &
                        z**6*(0.20222611593983328_dp) + z**5*(2.369844489519946_dp) + z*(32.13168724279836_dp) + &
                        z**3*(83.3378094939442_dp)) + z**6*(2159.5786990167953_dp) + Nf*(z**3*(-5014.506794109375_dp) + &
                        z**5*(-573.9266693274021_dp) + z*(-490.2772938483839_dp) + z**7*(-3.733700223649448_dp) + &
                        z**6*(45.16999768678995_dp) + (518.7498168642495_dp) + z**2*(582.477231789651_dp) + &
                        z**4*(5377.9879395868875_dp)) + z*(11757.52097459584_dp) + z**4*(142844.75300591116_dp))/z + &
                        (z**3*(-47650.75189886759_dp) + z**2*(-44262.56434834571_dp) + (-4425.496954200217_dp) + &
                        (-1721.4369617902225_dp)/z + Nf*(z**3*(-2557.2758045167684_dp) + z**2*(-1953.7330542758182_dp) + &
                        (-198.90742127519263_dp) + (62.84347343625167_dp)/z + z*(124.40761548401011_dp)) + &
                        z*(3168.551435857338_dp))*Log(z) + (z**2*(-6535.2817158676135_dp) + z*(-695.7032048642793_dp) + &
                        Nf*(z**2*(-303.81571620501387_dp) + (-88.09876543209877_dp) + z*(47.30864197530864_dp) + &
                        z**3*(1199.3961589988419_dp)) + (3229.449605856966_dp) + z**3*(37861.97093586552_dp))*Log(z)**2 + &
                        (z**3*(-670.403748053927_dp) + z**2*(-454.6027935456744_dp) + (-286.2222222222223_dp) + &
                        Nf*(z**3*(-83.65199321396027_dp) + z**2*(-21.332259023296757_dp) + (-11.45679012345679_dp) + &
                        z*(3.3580246913580245_dp)) + z*(180.2469135802469_dp))*Log(z)**3 + (z*(-89.23456790123457_dp) + &
                        z**2*(-12.508541392904073_dp) + Nf*((-5.925925925925926_dp) + z**2*(-0.6308139534883721_dp) + &
                        z*(2.962962962962963_dp) + z**3*(29.153427638737757_dp)) + (66.17283950617283_dp) + &
                        z**3*(1014.0617760617761_dp))*Log(z)**4 + (z**2*(-59144.057807899546_dp) + (-37843.51155524647_dp) + &
                        Nf**2*(z**2*(-58.76691672456219_dp) + (-54.007357506349955_dp) + z**3*(12.256322818982666_dp) + &
                        z*(89.45622301686777_dp)) + Nf*(z**2*(-908.5945922619609_dp) + z*(274.65703639520024_dp) + &
                        z**3*(426.1626588310546_dp) + (518.0412908934857_dp)) + z**3*(13525.511008216847_dp) + &
                        z*(80976.31519378009_dp))*Log(z*(-1._dp) + (1._dp)) + (z*(-22516.296518935298_dp) + &
                        z**3*(-10122.749948133172_dp) + Nf**2*((-13.115238157128097_dp) + z**2*(-9.31822567503274_dp) + &
                        z**3*(1.0644567219152854_dp) + z*(18.40604414728259_dp)) + Nf*(z*(-842.8843485514335_dp) + &
                        z**3*(-123.98394647948501_dp) + (490.70946029523355_dp) + z**2*(560.3069828838331_dp)) + &
                        (5134.980188237667_dp) + z**2*(26925.88322640089_dp))*Log(z*(-1._dp) + (1._dp))**2 + &
                        ((-246.6140084175432_dp) + z**3*(-65.88836403461971_dp) + z**2*(-28.91135442237791_dp) + &
                        Nf*(z**2*(-164.32072592121926_dp) + (-48.06436653669776_dp) + z**3*(55.9553401860879_dp) + &
                        z*(158.60259177800194_dp)) + z*(259.33965280046664_dp))*Log(z*(-1._dp) + (1._dp))**3 + &
                        (z*(-946.8394094774993_dp) + z**3*(-313.005082592122_dp) + (321.3616078871331_dp) + &
                        z**2*(941.4458471454511_dp))*Log(z*(-1._dp) + (1._dp))**4

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
                        Ibar_select = (z**5*(-2828.0740380454235_dp) + z*(-1921.906781696768_dp) + z**2*(-1704.9094560837823_dp) + &
                        z**7*(-59.52085659115983_dp) + Nf**2*(z**2*(-9.386600886600887_dp) + (-3.950617283950617_dp) + &
                        z**4*(-3.3880275624461667_dp) + z**6*(0.010230179028132991_dp) + z**7*(0.028615622583139988_dp) + &
                        z**5*(0.2859248341930729_dp) + z*(6.320987654320989_dp) + z**3*(8.104178371929093_dp)) + &
                        z**6*(214.83343351494656_dp) + Nf*(z*(-229.36726055725282_dp) + z**4*(-178.08205258067804_dp) + &
                        z**2*(-118.18652630559994_dp) + z**6*(-7.277252911524441_dp) + z**7*(4.054541893131286_dp) + &
                        z**5*(114.61486608653492_dp) + z**3*(187.97051798236524_dp) + (308.80249919264116_dp)) + &
                        z**4*(993.2566921070128_dp) + (1738.3991733727978_dp) + z**3*(3011.831480913418_dp))/z + &
                        (z*(18.837203522249037_dp) + (61.91367041742972_dp)/z + Nf*(z**3*(-119.74668719745489_dp) + (-64._dp) + &
                        z**2*(-53.043003643547664_dp) + (111.80246913580247_dp)/z + z*(164.74074074074073_dp)) + &
                        z**2*(358.55947574193874_dp) + z**3*(2470.547336337676_dp) + (3284.4219009299295_dp))*Log(z) + &
                        ((-398.2222222222225_dp) + Nf*(z*(-8.296296296296296_dp) + z**2*(-4.8845898593719355_dp) + &
                        (-4.740740740740743_dp) + z**3*(2.286801370196382_dp)) + z**2*(31.869564066314016_dp) + &
                        z**3*(304.94555354111714_dp) + z*(312.88888888888886_dp))*Log(z)**2 + (z*(-153.28395061728395_dp) + &
                        z**2*(1.9324055666003976_dp) + Nf*((-9.481481481481481_dp) + z**3*(-6.832326283987914_dp) + &
                        z**2*(-0.23712183156173347_dp) + z*(4.7407407407407405_dp)) + (105.87654320987656_dp) + &
                        z**3*(134.59642401021708_dp))*Log(z)**3 + (z**2*(-5382.39168676986_dp) + (-990.4770491448967_dp) + &
                        Nf**2*(z**2*(-6.915018963337549_dp) + (-6.22606452217072_dp) + z**3*(1.5296624057686001_dp) + &
                        z*(10.42623589455448_dp)) + Nf*(z**2*(-160.9406367492202_dp) + (-142.3551668762234_dp) + &
                        z**3*(40.700142343905874_dp) + z*(199.7808464667229_dp)) + z**3*(1543.4946737675941_dp) + &
                        z*(5460.799670155016_dp))*Log(z*(-1._dp) + (1._dp)) + (z*(-5983.700427542028_dp) + &
                        z**3*(-1373.5855449003866_dp) + Nf*(z**2*(-206.6029605030009_dp) + (-162.01908920860333_dp) + &
                        z**3*(50.14378749905067_dp) + z*(290.62641036070187_dp)) + (2939.872256135636_dp) + &
                        z**2*(4901.858160751223_dp))*Log(z*(-1._dp) + (1._dp))**2 + (z*(-384.98661958398407_dp) + &
                        z**3*(-2.737920937042457_dp) + z**2*(158.0816726790372_dp) + (306.6799048790264_dp))*Log(z*(-1._dp) + &
                        (1._dp))**3

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
                        Ibar_select = ((-1194.651740861711_dp) + z**3*(-548.7415523954633_dp) + z**4*(-365.93375082544543_dp) + &
                        z**6*(-13.901324599354476_dp) + z**7*(17.992705441428722_dp) + Nf*(z**4*(-71.06210158935419_dp) + &
                        z*(-39.50617283950617_dp) + z**2*(-24.010340118847466_dp) + z**6*(-7.181876292953364_dp) + &
                        z**7*(1.6204092319884238_dp) + z**5*(32.278259138366096_dp) + (34.69958847736624_dp) + &
                        z**3*(70.79186246049774_dp)) + z**5*(292.83877689474707_dp) + z**2*(905.1280202572963_dp) + &
                        z*(1112.1290850662403_dp))/z + ((-682.666666666667_dp)/z + (-223.99999999999997_dp) + &
                        z**3*(-47.29348722280993_dp) + z**2*(-5.9592490345760245_dp) + Nf*(z**3*(-3.706019918386917_dp) + &
                        z*(1.5802469135802417_dp) + (3.9506172839506166_dp) + z**2*(7.784015356314111_dp) + &
                        (13.432098765432102_dp)/z) + z*(26.46913580246906_dp))*Log(z) + (z*(-143.80246913580245_dp) + (-64._dp)/z + &
                        z**3*(-5.101605685706764_dp) + z**2*(0.383617193836172_dp) + Nf*(z**3*(-11.12646793134598_dp) + &
                        (-4.7407407407407405_dp) + z**2*(0.7233727810650887_dp) + z*(2.3703703703703702_dp)) + &
                        (52.93827160493826_dp))*Log(z)**2 + (z*(-1677.4861530471046_dp) + z**3*(-130.7999645636753_dp) + &
                        Nf*(z**2*(-45.82480713284435_dp) + (-37.63522883033307_dp) + z**3*(10.8108385339717_dp) + &
                        z*(65.93314804648969_dp)) + z**2*(894.5973443567198_dp) + (984.2072917725787_dp))*Log(z*(-1._dp) + (1._dp)) &
                        + (z**2*(-1021.2985497692816_dp) + (-690.9097544699171_dp) + z**3*(270.28081740276855_dp) + &
                        z*(1336.8410670833434_dp))*Log(z*(-1._dp) + (1._dp))**2

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
                        Ibar_select = ((-448840.96646975284_dp) + nf*z*(-31331.967578396743_dp) + &
                        nf**2*z**2*(-1134.8301644655623_dp) + nf**2*(-1007.3157828200781_dp) + nf**3*z*(-13.16872427983539_dp) + &
                        nf**3*z**2*(6.584362139917695_dp) + nf**3*(13.16872427983539_dp) + nf**2*z**3*(44.41241760643143_dp) + &
                        nf*z**3*(1172.801094028936_dp) + nf**2*z*(2150.309137727712_dp) + nf*z**2*(8511.8313666295_dp) + &
                        nf*(21930.47976848368_dp) + z**2*(167534.88287690943_dp) + z**3*(814262.7426688131_dp) + &
                        z*(2.3386935150672398e6_dp) + (z**2*(-82030.65262636205_dp) + z**3*(-18028.22868862438_dp) + &
                        nf**2*((-116.97308919809608_dp) + z*(-116.97308919809608_dp) + z**2*(-58.48654459904804_dp)) + &
                        nf*(z**2*(-3778.844691389914_dp) + z*(-1645.9994248071587_dp) + z**3*(402.0408399104933_dp) + &
                        (599.0174204688925_dp)) + z*(77085.34938778623_dp) + (132336.624120582_dp))*Hr1(-1) + &
                        z*(-2.3586903824615926e6_dp)*Hr1(0) + z**3*(-536025.5910362407_dp)*Hr1(0) + (-131107.12646633538_dp)*Hr1(0) &
                        + nf*z**2*(-3610.2342588963497_dp)*Hr1(0) + nf**2*(-247.67891370372882_dp)*Hr1(0) + &
                        nf**2*z**3*(69.70644718792867_dp)*Hr1(0) + nf**2*z**2*(131.331918081109_dp)*Hr1(0) + &
                        nf**2*z*(897.884940834345_dp)*Hr1(0) + nf*z**3*(5817.959843441199_dp)*Hr1(0) + &
                        nf*(6547.055322143833_dp)*Hr1(0) + nf*z*(24595.3479647672_dp)*Hr1(0) + z**2*(58353.4796403572_dp)*Hr1(0) + &
                        (-56708.400259556714_dp)*Hr1(1) + z**2*(-54666.26939681153_dp)*Hr1(1) + &
                        nf*z**2*(-8070.9281772369695_dp)*Hr1(1) + nf*(-6863.2890545543_dp)*Hr1(1) + &
                        nf**2*z*(-909.2218970132891_dp)*Hr1(1) + nf**3*(-15.802469135802466_dp)*Hr1(1) + &
                        nf**3*z**2*(-7.901234567901233_dp)*Hr1(1) + nf**3*z*(15.802469135802466_dp)*Hr1(1) + &
                        nf**2*z**3*(32.83401920438958_dp)*Hr1(1) + nf**2*z**2*(496.0430472720766_dp)*Hr1(1) + &
                        nf**2*(643.1039271916155_dp)*Hr1(1) + nf*z**3*(985.3966903522478_dp)*Hr1(1) + &
                        nf*z*(10806.880214259429_dp)*Hr1(1) + z*(38978.27081080509_dp)*Hr1(1) + z**3*(118714.30600772313_dp)*Hr1(1) &
                        + z*(-476566.67010474874_dp)*Hr2(-1,-1) + (-454184.5193371887_dp)*Hr2(-1,-1) + &
                        z**2*(-129294.3090425439_dp)*Hr2(-1,-1) + z**3*(-30740.527841259653_dp)*Hr2(-1,-1) + &
                        nf*(-389.91029732698695_dp)*Hr2(-1,-1) + nf*z*(-389.91029732698695_dp)*Hr2(-1,-1) + &
                        nf*z**2*(-194.95514866349347_dp)*Hr2(-1,-1) + nf*z**3*(-440.88888888888886_dp)*Hr2(-1,0) + &
                        nf**2*(-129.05349794238683_dp)*Hr2(-1,0) + nf**2*z*(-71.11111111111113_dp)*Hr2(-1,0) + &
                        nf**2*z**2*(-14.22222222222222_dp)*Hr2(-1,0) + nf**2*z**3*(-1.0534979423868314_dp)*Hr2(-1,0) + &
                        nf*z**2*(4053.281990852368_dp)*Hr2(-1,0) + nf*(5398.679208042184_dp)*Hr2(-1,0) + &
                        nf*z*(8105.5104837623485_dp)*Hr2(-1,0) + z**3*(151636.0816955407_dp)*Hr2(-1,0) + &
                        z**2*(191686.2273276768_dp)*Hr2(-1,0) + (743810.8599029998_dp)*Hr2(-1,0) + &
                        z*(758477.3230696958_dp)*Hr2(-1,0) + nf*z*(-6269.7575810179505_dp)*Hr2(0,-1) + &
                        nf*(-1881.9670350982574_dp)*Hr2(0,-1) + nf*z**2*(514.0317419760777_dp)*Hr2(0,-1) + &
                        z**3*(60092.97502403524_dp)*Hr2(0,-1) + z**2*(117885.54486860191_dp)*Hr2(0,-1) + &
                        (270064.6684269045_dp)*Hr2(0,-1) + z*(483386.71228956647_dp)*Hr2(0,-1) + &
                        z**2*(-81526.44374603553_dp)*Hr2(0,0) + z**3*(-52124.24008158415_dp)*Hr2(0,0) + &
                        nf*z*(-19590.11051056837_dp)*Hr2(0,0) + (-17464.21450082806_dp)*Hr2(0,0) + &
                        nf*z**3*(-5591.0013717421125_dp)*Hr2(0,0) + nf*z**2*(-294.37583467602235_dp)*Hr2(0,0) + &
                        nf**2*z**2*(-97.74248856027367_dp)*Hr2(0,0) + nf**2*z**3*(2.1069958847736627_dp)*Hr2(0,0) + &
                        nf*(311.6525454325509_dp)*Hr2(0,0) + nf**2*z*(617.2133721822757_dp)*Hr2(0,0) + &
                        z*(292510.8626011058_dp)*Hr2(0,0) + z*(-362584.25421264605_dp)*Hr2(0,1) + &
                        z**3*(-154926.08835327407_dp)*Hr2(0,1) + z**2*(-21259.916146243893_dp)*Hr2(0,1) + &
                        nf*z*(-7465.643472000413_dp)*Hr2(0,1) + nf**2*z*(-183.76954732510285_dp)*Hr2(0,1) + &
                        nf*z**3*(-97.49245541838134_dp)*Hr2(0,1) + nf**2*z**3*(-43.7201646090535_dp)*Hr2(0,1) + &
                        nf**2*(27.65432098765432_dp)*Hr2(0,1) + nf**2*z**2*(39.93415637860082_dp)*Hr2(0,1) + &
                        nf*z**2*(513.8385814544104_dp)*Hr2(0,1) + nf*(5026.557619418598_dp)*Hr2(0,1) + &
                        (87048.94765474467_dp)*Hr2(0,1) + z*(-535428.4165642223_dp)*Hr2(1,0) + &
                        z**3*(-129872.55767382082_dp)*Hr2(1,0) + nf*(-9406.302869873709_dp)*Hr2(1,0) + &
                        nf*z**2*(-3173.2200220424793_dp)*Hr2(1,0) + nf**2*z*(-410.3374485596708_dp)*Hr2(1,0) + &
                        nf**2*z**3*(-36.8724279835391_dp)*Hr2(1,0) + nf**2*z**2*(180.67489711934155_dp)*Hr2(1,0) + &
                        nf**2*(306.56790123456784_dp)*Hr2(1,0) + nf*z**3*(1504.3072702331958_dp)*Hr2(1,0) + &
                        nf*z*(9769.25485889977_dp)*Hr2(1,0) + z**2*(114240.16172832364_dp)*Hr2(1,0) + &
                        (572333.0508149313_dp)*Hr2(1,0) + (z*(-58708.28894734273_dp) + z**2*(-15896.12046597938_dp) + &
                        nf**3*(z*(-9.481481481481481_dp) + z**2*(4.7407407407407405_dp) + (9.481481481481481_dp)) + &
                        nf**2*((-169.4814814814815_dp) + z**2*(8.032921810699587_dp) + z**3*(27.91769547325103_dp) + &
                        z*(185.4156378600823_dp)) + z**3*(858.1874784310539_dp) + nf*(z*(-6795.969157487582_dp) + &
                        z**3*(-2026.9739368998628_dp) + z**2*(365.4907515832969_dp) + (4592.6221067331235_dp)) + &
                        (122628.56250812128_dp))*Hr2(1,1) + z**2*(114277.5093425755_dp)*Hr3(-1,-1,-1) + &
                        (228555.018685151_dp)*Hr3(-1,-1,-1) + z*(228555.018685151_dp)*Hr3(-1,-1,-1) + &
                        z*(-310795.88848718145_dp)*Hr3(-1,-1,0) + (-276253.2712032309_dp)*Hr3(-1,-1,0) + &
                        z**2*(-149860.36399667716_dp)*Hr3(-1,-1,0) + z**3*(-34637.432098765436_dp)*Hr3(-1,-1,0) + &
                        nf*z**2*(-2719.341563786008_dp)*Hr3(-1,-1,0) + nf*z*(-2350.090534979424_dp)*Hr3(-1,-1,0) + &
                        nf*(-1260.2469135802467_dp)*Hr3(-1,-1,0) + nf**2*(-28.44444444444444_dp)*Hr3(-1,-1,0) + &
                        nf**2*z*(-28.44444444444444_dp)*Hr3(-1,-1,0) + nf**2*z**2*(-14.22222222222222_dp)*Hr3(-1,-1,0) + &
                        nf*z**3*(307.6213991769547_dp)*Hr3(-1,-1,0) + (-173200.7534746298_dp)*Hr3(-1,0,-1) + &
                        z*(-173200.7534746298_dp)*Hr3(-1,0,-1) + z**2*(-86600.3767373149_dp)*Hr3(-1,0,-1) + &
                        (-151208.35005038403_dp)*Hr3(-1,0,0) + z*(-124810.52288989021_dp)*Hr3(-1,0,0) + &
                        z**3*(-15776.395061728397_dp)*Hr3(-1,0,0) + nf**2*z**2*(7.11111111111111_dp)*Hr3(-1,0,0) + &
                        nf**2*(14.22222222222222_dp)*Hr3(-1,0,0) + nf**2*z*(14.22222222222222_dp)*Hr3(-1,0,0) + &
                        nf*z**3*(530.9629629629629_dp)*Hr3(-1,0,0) + nf*z**2*(5269.991769547324_dp)*Hr3(-1,0,0) + &
                        nf*(7097.415637860082_dp)*Hr3(-1,0,0) + nf*z*(10219.588477366255_dp)*Hr3(-1,0,0) + &
                        z**2*(41878.78793777095_dp)*Hr3(-1,0,0) + (-43902.28314532855_dp)*Hr3(-1,0,1) + &
                        z*(-7309.690552735955_dp)*Hr3(-1,0,1) + nf*(-183.04526748971193_dp)*Hr3(-1,0,1) + &
                        nf*z**3*(-90.60082304526749_dp)*Hr3(-1,0,1) + nf**2*z**2*(28.44444444444444_dp)*Hr3(-1,0,1) + &
                        nf**2*(56.88888888888888_dp)*Hr3(-1,0,1) + nf**2*z*(56.88888888888888_dp)*Hr3(-1,0,1) + &
                        nf*z*(636.8395061728395_dp)*Hr3(-1,0,1) + nf*z**2*(1343.20987654321_dp)*Hr3(-1,0,1) + &
                        z**3*(6143.999999999999_dp)*Hr3(-1,0,1) + z**2*(33675.72262486659_dp)*Hr3(-1,0,1) + &
                        (-165867.84048290024_dp)*Hr3(0,-1,-1) + z**2*(-124724.50590895659_dp)*Hr3(0,-1,-1) + &
                        z*(-72674.080617786_dp)*Hr3(0,-1,-1) + z*(-268514.7571514078_dp)*Hr3(0,-1,0) + &
                        z**3*(-38848.79012345678_dp)*Hr3(0,-1,0) + nf**2*(-88.49382716049382_dp)*Hr3(0,-1,0) + &
                        nf**2*z**2*(-47.407407407407405_dp)*Hr3(0,-1,0) + nf**2*z*(-18.962962962962962_dp)*Hr3(0,-1,0) + &
                        nf*z**3*(811.1934156378601_dp)*Hr3(0,-1,0) + nf*z**2*(3988.4115226337453_dp)*Hr3(0,-1,0) + &
                        nf*(6772.148148148148_dp)*Hr3(0,-1,0) + nf*z*(7662.748971193415_dp)*Hr3(0,-1,0) + &
                        (165841.31025686016_dp)*Hr3(0,-1,0) + z**2*(167993.52269709698_dp)*Hr3(0,-1,0) + &
                        z*(-171650.2101922595_dp)*Hr3(0,0,-1) + nf*z*(-148.16591298425504_dp)*Hr3(0,0,-1) + &
                        nf*z**2*(-74.08295649212752_dp)*Hr3(0,0,-1) + (80851.79925372402_dp)*Hr3(0,0,-1) + &
                        z**2*(85280.53038086304_dp)*Hr3(0,0,-1) + z**2*(-12160.59992620608_dp)*Hr3(0,0,0) + &
                        nf*z**3*(-1070.3539094650205_dp)*Hr3(0,0,0) + nf**2*z**2*(-230.71604938271605_dp)*Hr3(0,0,0) + &
                        nf**2*z*(423.5061728395062_dp)*Hr3(0,0,0) + nf*z**2*(4763.5361792565_dp)*Hr3(0,0,0) + &
                        nf*z*(6485.806718636753_dp)*Hr3(0,0,0) + z**3*(91226.07407407407_dp)*Hr3(0,0,0) + &
                        z*(207503.07488398952_dp)*Hr3(0,0,0) + nf*z**3*(-314.4691358024691_dp)*Hr3(0,0,1) + &
                        nf**2*z*(-69.53086419753086_dp)*Hr3(0,0,1) + nf**2*z**2*(-12.641975308641973_dp)*Hr3(0,0,1) + &
                        nf*z*(13.195015494467752_dp)*Hr3(0,0,1) + nf*z**2*(2199.723479907087_dp)*Hr3(0,0,1) + &
                        nf*(3109.925925925926_dp)*Hr3(0,0,1) + z**2*(13241.997498988538_dp)*Hr3(0,0,1) + &
                        z**3*(50767.8024691358_dp)*Hr3(0,0,1) + (55236.911862555935_dp)*Hr3(0,0,1) + &
                        z*(322220.608246542_dp)*Hr3(0,0,1) + nf*z*(-2334.0246913580245_dp)*Hr3(0,1,0) + &
                        nf*(-2110.683127572016_dp)*Hr3(0,1,0) + nf**2*z*(-41.08641975308642_dp)*Hr3(0,1,0) + &
                        nf**2*z**2*(34.76543209876543_dp)*Hr3(0,1,0) + nf**2*(139.06172839506172_dp)*Hr3(0,1,0) + &
                        nf*z**3*(183.83539094650206_dp)*Hr3(0,1,0) + nf*z**2*(2316.641975308642_dp)*Hr3(0,1,0) + &
                        z**3*(6251.456790123457_dp)*Hr3(0,1,0) + z**2*(74261.05415926456_dp)*Hr3(0,1,0) + &
                        (143231.94737401695_dp)*Hr3(0,1,0) + z*(226758.8904541099_dp)*Hr3(0,1,0) + &
                        z**3*(-37969.38271604938_dp)*Hr3(0,1,1) + z*(-25058.89409556614_dp)*Hr3(0,1,1) + &
                        nf*(-243.0946502057613_dp)*Hr3(0,1,1) + nf**2*z*(-194.37037037037038_dp)*Hr3(0,1,1) + &
                        nf**2*z**2*(4.7407407407407405_dp)*Hr3(0,1,1) + nf**2*(58.469135802469125_dp)*Hr3(0,1,1) + &
                        nf*z**3*(1041.3827160493827_dp)*Hr3(0,1,1) + nf*z**2*(2341.794238683128_dp)*Hr3(0,1,1) + &
                        nf*z*(4330.139917695473_dp)*Hr3(0,1,1) + z**2*(26472.82408541913_dp)*Hr3(0,1,1) + &
                        (67246.48059549867_dp)*Hr3(0,1,1) + (-14418.882795151978_dp)*Hr3(1,0,-1) + &
                        z**2*(-7209.441397575989_dp)*Hr3(1,0,-1) + z*(14418.882795151978_dp)*Hr3(1,0,-1) + &
                        (-110732.34028698462_dp)*Hr3(1,0,0) + z**2*(-10282.466439788619_dp)*Hr3(1,0,0) + &
                        nf*z*(-9258.732510288064_dp)*Hr3(1,0,0) + z**3*(-4472.888888888889_dp)*Hr3(1,0,0) + &
                        nf*z**3*(-279.1769547325103_dp)*Hr3(1,0,0) + nf**2*(-45.82716049382716_dp)*Hr3(1,0,0) + &
                        nf**2*z**2*(-22.91358024691358_dp)*Hr3(1,0,0) + nf**2*z*(45.82716049382716_dp)*Hr3(1,0,0) + &
                        nf*z**2*(4494.551440329218_dp)*Hr3(1,0,0) + nf*(7626.798353909465_dp)*Hr3(1,0,0) + &
                        z*(116645.52547216983_dp)*Hr3(1,0,0) + z*(-161075.14870203787_dp)*Hr3(1,0,1) + &
                        z**3*(-48170.66666666666_dp)*Hr3(1,0,1) + nf*(-2376.8230452674898_dp)*Hr3(1,0,1) + &
                        nf**2*z*(-45.82716049382716_dp)*Hr3(1,0,1) + nf**2*z**2*(22.91358024691358_dp)*Hr3(1,0,1) + &
                        nf**2*(45.82716049382716_dp)*Hr3(1,0,1) + nf*z**3*(696.8888888888889_dp)*Hr3(1,0,1) + &
                        nf*z**2*(1342.5514403292182_dp)*Hr3(1,0,1) + nf*z*(2489.4156378600824_dp)*Hr3(1,0,1) + &
                        z**2*(32344.191634969546_dp)*Hr3(1,0,1) + (159334.70425759343_dp)*Hr3(1,0,1) + &
                        z*(-184064.60010121838_dp)*Hr3(1,1,0) + z**3*(-48853.333333333336_dp)*Hr3(1,1,0) + &
                        nf*(-3106.8971193415637_dp)*Hr3(1,1,0) + nf**2*z*(-58.469135802469125_dp)*Hr3(1,1,0) + &
                        nf**2*z**2*(29.234567901234563_dp)*Hr3(1,1,0) + nf**2*(58.469135802469125_dp)*Hr3(1,1,0) + &
                        nf*z**3*(734.8148148148148_dp)*Hr3(1,1,0) + nf*z**2*(1130.798353909465_dp)*Hr3(1,1,0) + &
                        nf*z*(3241.61316872428_dp)*Hr3(1,1,0) + z**2*(40398.71980369561_dp)*Hr3(1,1,0) + &
                        (182605.43960739122_dp)*Hr3(1,1,0) + z*(-25774.13463222268_dp)*Hr3(1,1,1) + &
                        z**3*(-16872.296296296296_dp)*Hr3(1,1,1) + nf*z*(-6466.3703703703695_dp)*Hr3(1,1,1) + &
                        nf**2*(-216.4938271604938_dp)*Hr3(1,1,1) + nf*z**3*(-182.78189300411523_dp)*Hr3(1,1,1) + &
                        nf**2*z**2*(-108.2469135802469_dp)*Hr3(1,1,1) + nf**2*z*(216.4938271604938_dp)*Hr3(1,1,1) + &
                        z**2*(1003.80805685208_dp)*Hr3(1,1,1) + nf*z**2*(3363.16049382716_dp)*Hr3(1,1,1) + &
                        nf*(7236.609053497942_dp)*Hr3(1,1,1) + (14475.171669259715_dp)*Hr3(1,1,1) + &
                        z*(-62072.098765432085_dp)*Hr4(-1,-1,-1,0) + (-47962.074074074066_dp)*Hr4(-1,-1,-1,0) + &
                        z**2*(-8359.506172839507_dp)*Hr4(-1,-1,-1,0) + z**3*(-4949.333333333333_dp)*Hr4(-1,-1,-1,0) + &
                        nf*(-916.5432098765432_dp)*Hr4(-1,-1,-1,0) + nf*z*(-916.5432098765432_dp)*Hr4(-1,-1,-1,0) + &
                        nf*z**2*(-458.2716049382716_dp)*Hr4(-1,-1,-1,0) + nf*(-1223.111111111111_dp)*Hr4(-1,-1,0,0) + &
                        nf*z*(-1223.111111111111_dp)*Hr4(-1,-1,0,0) + nf*z**2*(-611.5555555555555_dp)*Hr4(-1,-1,0,0) + &
                        z**2*(17904.197530864196_dp)*Hr4(-1,-1,0,0) + z**3*(24841.48148148148_dp)*Hr4(-1,-1,0,0) + &
                        z*(222137.67901234564_dp)*Hr4(-1,-1,0,0) + (233197.03703703708_dp)*Hr4(-1,-1,0,0) + &
                        nf*(-221.23456790123456_dp)*Hr4(-1,-1,0,1) + nf*z*(-221.23456790123456_dp)*Hr4(-1,-1,0,1) + &
                        nf*z**2*(-110.61728395061728_dp)*Hr4(-1,-1,0,1) + z**3*(16213.333333333334_dp)*Hr4(-1,-1,0,1) + &
                        z**2*(21537.185185185186_dp)*Hr4(-1,-1,0,1) + (146360.88888888885_dp)*Hr4(-1,-1,0,1) + &
                        z*(152912.59259259255_dp)*Hr4(-1,-1,0,1) + nf*(-767.9999999999999_dp)*Hr4(-1,0,-1,0) + &
                        nf*z*(-767.9999999999999_dp)*Hr4(-1,0,-1,0) + nf*z**2*(-383.99999999999994_dp)*Hr4(-1,0,-1,0) + &
                        z**2*(8455.901234567902_dp)*Hr4(-1,0,-1,0) + z**3*(10799.407407407407_dp)*Hr4(-1,0,-1,0) + &
                        z*(103028.93827160494_dp)*Hr4(-1,0,-1,0) + (108532.14814814813_dp)*Hr4(-1,0,-1,0) + &
                        (-184391.11111111112_dp)*Hr4(-1,0,0,0) + z*(-169404.04938271604_dp)*Hr4(-1,0,0,0) + &
                        z**3*(-20034.37037037037_dp)*Hr4(-1,0,0,0) + z**2*(-11860.938271604939_dp)*Hr4(-1,0,0,0) + &
                        nf*z**2*(643.9506172839505_dp)*Hr4(-1,0,0,0) + nf*(1287.901234567901_dp)*Hr4(-1,0,0,0) + &
                        nf*z*(1287.901234567901_dp)*Hr4(-1,0,0,0) + z*(-188609.58024691354_dp)*Hr4(-1,0,0,1) + &
                        (-184820.14814814815_dp)*Hr4(-1,0,0,1) + z**2*(-26940.839506172837_dp)*Hr4(-1,0,0,1) + &
                        z**3*(-20906.666666666668_dp)*Hr4(-1,0,0,1) + nf*z**2*(311.30864197530866_dp)*Hr4(-1,0,0,1) + &
                        nf*(622.6172839506173_dp)*Hr4(-1,0,0,1) + nf*z*(622.6172839506173_dp)*Hr4(-1,0,0,1) + &
                        z*(-87340.24691358025_dp)*Hr4(-1,0,1,0) + (-81524.14814814815_dp)*Hr4(-1,0,1,0) + &
                        z**2*(-17805.432098765432_dp)*Hr4(-1,0,1,0) + z**3*(-9130.666666666666_dp)*Hr4(-1,0,1,0) + &
                        nf*z**2*(199.11111111111111_dp)*Hr4(-1,0,1,0) + nf*(398.22222222222223_dp)*Hr4(-1,0,1,0) + &
                        nf*z*(398.22222222222223_dp)*Hr4(-1,0,1,0) + z*(-66373.53086419753_dp)*Hr4(-1,0,1,1) + &
                        (-36683.851851851854_dp)*Hr4(-1,0,1,1) + z**2*(-27684.345679012345_dp)*Hr4(-1,0,1,1) + &
                        z**3*(-6143.999999999999_dp)*Hr4(-1,0,1,1) + nf*(-1194.6666666666667_dp)*Hr4(-1,0,1,1) + &
                        nf*z*(-1194.6666666666667_dp)*Hr4(-1,0,1,1) + nf*z**2*(-597.3333333333334_dp)*Hr4(-1,0,1,1) + &
                        nf*z*(-2755.9506172839506_dp)*Hr4(0,-1,-1,0) + nf*(-670.0246913580247_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**2*(460.64197530864203_dp)*Hr4(0,-1,-1,0) + z**2*(1647.4074074074074_dp)*Hr4(0,-1,-1,0) + &
                        z**3*(22888.296296296296_dp)*Hr4(0,-1,-1,0) + (37333.333333333336_dp)*Hr4(0,-1,-1,0) + &
                        z*(168105.8765432099_dp)*Hr4(0,-1,-1,0) + (-158293.33333333334_dp)*Hr4(0,-1,0,0) + &
                        z*(-147697.77777777778_dp)*Hr4(0,-1,0,0) + z**3*(-32104.296296296292_dp)*Hr4(0,-1,0,0) + &
                        z**2*(-23371.851851851847_dp)*Hr4(0,-1,0,0) + nf*z**2*(57.28395061728395_dp)*Hr4(0,-1,0,0) + &
                        nf*z*(1460.148148148148_dp)*Hr4(0,-1,0,0) + nf*(2929.777777777778_dp)*Hr4(0,-1,0,0) + &
                        z*(-161584.19753086418_dp)*Hr4(0,-1,0,1) + (-74538.66666666667_dp)*Hr4(0,-1,0,1) + &
                        z**3*(-25087.999999999996_dp)*Hr4(0,-1,0,1) + z**2*(-20520.296296296296_dp)*Hr4(0,-1,0,1) + &
                        nf*z**2*(-82.17283950617283_dp)*Hr4(0,-1,0,1) + nf*(809.0864197530863_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(2433.5802469135797_dp)*Hr4(0,-1,0,1) + (-20352._dp)*Hr4(0,0,-1,0) + &
                        z**3*(-14885.925925925925_dp)*Hr4(0,0,-1,0) + z**2*(-10736.197530864198_dp)*Hr4(0,0,-1,0) + &
                        nf*z*(-1294.2222222222222_dp)*Hr4(0,0,-1,0) + nf*z**2*(692.5432098765432_dp)*Hr4(0,0,-1,0) + &
                        nf*(1308.4444444444443_dp)*Hr4(0,0,-1,0) + z*(31936.000000000004_dp)*Hr4(0,0,-1,0) + &
                        nf*z*(-3136.7901234567903_dp)*Hr4(0,0,0,0) + nf**2*z**2*(-56.88888888888888_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(113.77777777777776_dp)*Hr4(0,0,0,0) + nf*z**2*(2085.5308641975307_dp)*Hr4(0,0,0,0) + &
                        z**2*(4423.901234567901_dp)*Hr4(0,0,0,0) + z**3*(18052.74074074074_dp)*Hr4(0,0,0,0) + &
                        z*(95419.25925925926_dp)*Hr4(0,0,0,0) + nf*z**2*(-1245.827160493827_dp)*Hr4(0,0,0,1) + &
                        nf**2*z*(-151.7037037037037_dp)*Hr4(0,0,0,1) + nf**2*z**2*(75.85185185185185_dp)*Hr4(0,0,0,1) + &
                        nf*z*(1837.037037037037_dp)*Hr4(0,0,0,1) + z**3*(34417.77777777778_dp)*Hr4(0,0,0,1) + &
                        z*(46674.567901234564_dp)*Hr4(0,0,0,1) + z**2*(64026.864197530864_dp)*Hr4(0,0,0,1) + &
                        (-8832._dp)*Hr4(0,0,1,0) + nf*z**2*(-2405.135802469136_dp)*Hr4(0,0,1,0) + &
                        nf*(-767.9999999999999_dp)*Hr4(0,0,1,0) + nf*z*(-679.5061728395061_dp)*Hr4(0,0,1,0) + &
                        nf**2*z*(-37.925925925925924_dp)*Hr4(0,0,1,0) + nf**2*z**2*(18.962962962962962_dp)*Hr4(0,0,1,0) + &
                        z**3*(28890.07407407407_dp)*Hr4(0,0,1,0) + z**2*(81443.16049382715_dp)*Hr4(0,0,1,0) + &
                        z*(102607.4074074074_dp)*Hr4(0,0,1,0) + (-7104._dp)*Hr4(0,0,1,1) + nf*z*(-4612.740740740741_dp)*Hr4(0,0,1,1) &
                        + nf*z**2*(-2304._dp)*Hr4(0,0,1,1) + nf*(-1265.7777777777776_dp)*Hr4(0,0,1,1) + &
                        nf**2*z**2*(-9.481481481481481_dp)*Hr4(0,0,1,1) + nf**2*z*(18.962962962962962_dp)*Hr4(0,0,1,1) + &
                        z**3*(28946.96296296296_dp)*Hr4(0,0,1,1) + z**2*(85160.29629629628_dp)*Hr4(0,0,1,1) + &
                        z*(204654.6172839506_dp)*Hr4(0,0,1,1) + (-82474.66666666669_dp)*Hr4(0,1,0,0) + &
                        nf*z**2*(-174.41975308641975_dp)*Hr4(0,1,0,0) + nf*z*(607.604938271605_dp)*Hr4(0,1,0,0) + &
                        nf*(1065.0864197530864_dp)*Hr4(0,1,0,0) + z**3*(20186.074074074073_dp)*Hr4(0,1,0,0) + &
                        z**2*(22996.938271604937_dp)*Hr4(0,1,0,0) + z*(71173.92592592593_dp)*Hr4(0,1,0,0) + &
                        nf*(-3043.555555555555_dp)*Hr4(0,1,0,1) + nf*z**2*(-2287.4074074074074_dp)*Hr4(0,1,0,1) + &
                        nf*z*(37.925925925925924_dp)*Hr4(0,1,0,1) + z**3*(20167.11111111111_dp)*Hr4(0,1,0,1) + &
                        (21994.666666666668_dp)*Hr4(0,1,0,1) + z**2*(68897.58024691358_dp)*Hr4(0,1,0,1) + &
                        z*(118219.85185185185_dp)*Hr4(0,1,0,1) + nf*(-3157.333333333333_dp)*Hr4(0,1,1,0) + &
                        nf*z**2*(-2445.432098765432_dp)*Hr4(0,1,1,0) + nf*z*(12.641975308641973_dp)*Hr4(0,1,1,0) + &
                        z**3*(19787.85185185185_dp)*Hr4(0,1,1,0) + (24042.666666666668_dp)*Hr4(0,1,1,0) + &
                        z**2*(71678.81481481482_dp)*Hr4(0,1,1,0) + z*(116462.6172839506_dp)*Hr4(0,1,1,0) + &
                        z*(-20157.62962962963_dp)*Hr4(0,1,1,1) + nf*(-3760.9876543209875_dp)*Hr4(0,1,1,1) + &
                        nf*z**2*(-1582.6172839506173_dp)*Hr4(0,1,1,1) + nf*z*(5142.123456790123_dp)*Hr4(0,1,1,1) + &
                        z**3*(8391.111111111113_dp)*Hr4(0,1,1,1) + z**2*(34881.18518518518_dp)*Hr4(0,1,1,1) + &
                        (49642.66666666666_dp)*Hr4(0,1,1,1) + z*(-22079.20987654321_dp)*Hr4(1,0,-1,0) + &
                        z**2*(-6012.049382716049_dp)*Hr4(1,0,-1,0) + z**3*(-4579.555555555555_dp)*Hr4(1,0,-1,0) + &
                        nf*z*(-265.48148148148147_dp)*Hr4(1,0,-1,0) + nf*z**2*(132.74074074074073_dp)*Hr4(1,0,-1,0) + &
                        nf*(265.48148148148147_dp)*Hr4(1,0,-1,0) + (30580.14814814815_dp)*Hr4(1,0,-1,0) + &
                        (-116683.85185185185_dp)*Hr4(1,0,0,0) + nf*z*(-1240.4938271604938_dp)*Hr4(1,0,0,0) + &
                        nf*z**2*(620.2469135802469_dp)*Hr4(1,0,0,0) + nf*(1240.4938271604938_dp)*Hr4(1,0,0,0) + &
                        z**2*(1366.9135802469134_dp)*Hr4(1,0,0,0) + z**3*(11557.925925925925_dp)*Hr4(1,0,0,0) + &
                        z*(98033.38271604938_dp)*Hr4(1,0,0,0) + (-76984.88888888889_dp)*Hr4(1,0,0,1) + &
                        nf*(-524.641975308642_dp)*Hr4(1,0,0,1) + nf*z**2*(-262.320987654321_dp)*Hr4(1,0,0,1) + &
                        nf*z*(524.641975308642_dp)*Hr4(1,0,0,1) + z**3*(8609.185185185182_dp)*Hr4(1,0,0,1) + &
                        z**2*(14700.64197530864_dp)*Hr4(1,0,0,1) + z*(64019.35802469135_dp)*Hr4(1,0,0,1) + &
                        (-30584.888888888887_dp)*Hr4(1,0,1,0) + nf*(-2784.395061728395_dp)*Hr4(1,0,1,0) + &
                        nf*z**2*(-1392.1975308641975_dp)*Hr4(1,0,1,0) + nf*z*(2784.395061728395_dp)*Hr4(1,0,1,0) + &
                        z**3*(6855.111111111111_dp)*Hr4(1,0,1,0) + z*(20807.11111111111_dp)*Hr4(1,0,1,0) + &
                        z**2*(30364.44444444444_dp)*Hr4(1,0,1,0) + z*(-48923.654320987655_dp)*Hr4(1,0,1,1) + &
                        nf*(-4144.987654320987_dp)*Hr4(1,0,1,1) + nf*z**2*(-2072.4938271604933_dp)*Hr4(1,0,1,1) + &
                        z**3*(-1270.5185185185182_dp)*Hr4(1,0,1,1) + nf*z*(4144.987654320987_dp)*Hr4(1,0,1,1) + &
                        z**2*(36288.39506172839_dp)*Hr4(1,0,1,1) + (51804.44444444444_dp)*Hr4(1,0,1,1) + &
                        (-82192.5925925926_dp)*Hr4(1,1,0,0) + nf*(-398.22222222222223_dp)*Hr4(1,1,0,0) + &
                        nf*z**2*(-199.11111111111111_dp)*Hr4(1,1,0,0) + nf*z*(398.22222222222223_dp)*Hr4(1,1,0,0) + &
                        z**3*(9092.74074074074_dp)*Hr4(1,1,0,0) + z**2*(14774.123456790123_dp)*Hr4(1,1,0,0) + &
                        z*(66617.28395061729_dp)*Hr4(1,1,0,0) + z*(-50396.444444444445_dp)*Hr4(1,1,0,1) + &
                        nf*(-4246.123456790124_dp)*Hr4(1,1,0,1) + nf*z**2*(-2123.061728395062_dp)*Hr4(1,1,0,1) + &
                        z**3*(-1080.8888888888887_dp)*Hr4(1,1,0,1) + nf*z*(4246.123456790124_dp)*Hr4(1,1,0,1) + &
                        z**2*(37562.07407407407_dp)*Hr4(1,1,0,1) + (52250.07407407407_dp)*Hr4(1,1,0,1) + &
                        z*(-58790.71604938272_dp)*Hr4(1,1,1,0) + nf*(-4536.888888888889_dp)*Hr4(1,1,1,0) + &
                        nf*z**2*(-2268.4444444444443_dp)*Hr4(1,1,1,0) + z**3*(-1460.148148148148_dp)*Hr4(1,1,1,0) + &
                        nf*z*(4536.888888888889_dp)*Hr4(1,1,1,0) + z**2*(40068.345679012345_dp)*Hr4(1,1,1,0) + &
                        (60233.481481481474_dp)*Hr4(1,1,1,0) + z*(-31554.37037037037_dp)*Hr4(1,1,1,1) + &
                        z**3*(-1991.111111111111_dp)*Hr4(1,1,1,1) + nf*(-1921.5802469135801_dp)*Hr4(1,1,1,1) + &
                        nf*z**2*(-960.7901234567901_dp)*Hr4(1,1,1,1) + nf*z*(1921.5802469135801_dp)*Hr4(1,1,1,1) + &
                        z**2*(16227.555555555553_dp)*Hr4(1,1,1,1) + (35986.96296296296_dp)*Hr4(1,1,1,1) + &
                        z**2*(33005.03703703703_dp)*Hr5(-1,-1,-1,-1,0) + (66010.07407407406_dp)*Hr5(-1,-1,-1,-1,0) + &
                        z*(66010.07407407406_dp)*Hr5(-1,-1,-1,-1,0) + (-135885.4320987654_dp)*Hr5(-1,-1,-1,0,0) + &
                        z*(-135885.4320987654_dp)*Hr5(-1,-1,-1,0,0) + z**2*(-67942.7160493827_dp)*Hr5(-1,-1,-1,0,0) + &
                        (-105939.75308641975_dp)*Hr5(-1,-1,-1,0,1) + z*(-105939.75308641975_dp)*Hr5(-1,-1,-1,0,1) + &
                        z**2*(-52969.87654320987_dp)*Hr5(-1,-1,-1,0,1) + (-61999.4074074074_dp)*Hr5(-1,-1,0,-1,0) + &
                        z*(-61999.4074074074_dp)*Hr5(-1,-1,0,-1,0) + z**2*(-30999.7037037037_dp)*Hr5(-1,-1,0,-1,0) + &
                        z**2*(46987.06172839506_dp)*Hr5(-1,-1,0,0,0) + (93974.12345679013_dp)*Hr5(-1,-1,0,0,0) + &
                        z*(93974.12345679013_dp)*Hr5(-1,-1,0,0,0) + z**2*(60918.518518518526_dp)*Hr5(-1,-1,0,0,1) + &
                        (121837.03703703705_dp)*Hr5(-1,-1,0,0,1) + z*(121837.03703703705_dp)*Hr5(-1,-1,0,0,1) + &
                        z**2*(27810.76543209876_dp)*Hr5(-1,-1,0,1,0) + (55621.53086419752_dp)*Hr5(-1,-1,0,1,0) + &
                        z*(55621.53086419752_dp)*Hr5(-1,-1,0,1,0) + z**2*(34297.67901234568_dp)*Hr5(-1,-1,0,1,1) + &
                        (68595.35802469136_dp)*Hr5(-1,-1,0,1,1) + z*(68595.35802469136_dp)*Hr5(-1,-1,0,1,1) + &
                        (-61999.4074074074_dp)*Hr5(-1,0,-1,-1,0) + z*(-61999.4074074074_dp)*Hr5(-1,0,-1,-1,0) + &
                        z**2*(-30999.7037037037_dp)*Hr5(-1,0,-1,-1,0) + z**2*(46550.91358024691_dp)*Hr5(-1,0,-1,0,0) + &
                        (93101.82716049382_dp)*Hr5(-1,0,-1,0,0) + z*(93101.82716049382_dp)*Hr5(-1,0,-1,0,0) + &
                        z**2*(37146.864197530864_dp)*Hr5(-1,0,-1,0,1) + (74293.72839506173_dp)*Hr5(-1,0,-1,0,1) + &
                        z*(74293.72839506173_dp)*Hr5(-1,0,-1,0,1) + z**2*(16681.086419753086_dp)*Hr5(-1,0,0,-1,0) + &
                        (33362.17283950617_dp)*Hr5(-1,0,0,-1,0) + z*(33362.17283950617_dp)*Hr5(-1,0,0,-1,0) + &
                        (-33387.456790123455_dp)*Hr5(-1,0,0,0,0) + z*(-33387.456790123455_dp)*Hr5(-1,0,0,0,0) + &
                        z**2*(-16693.728395061727_dp)*Hr5(-1,0,0,0,0) + (-66253.43209876545_dp)*Hr5(-1,0,0,0,1) + &
                        z*(-66253.43209876545_dp)*Hr5(-1,0,0,0,1) + z**2*(-33126.716049382725_dp)*Hr5(-1,0,0,0,1) + &
                        (-47024.98765432099_dp)*Hr5(-1,0,0,1,0) + z*(-47024.98765432099_dp)*Hr5(-1,0,0,1,0) + &
                        z**2*(-23512.493827160495_dp)*Hr5(-1,0,0,1,0) + (-59072.79012345679_dp)*Hr5(-1,0,0,1,1) + &
                        z*(-59072.79012345679_dp)*Hr5(-1,0,0,1,1) + z**2*(-29536.395061728395_dp)*Hr5(-1,0,0,1,1) + &
                        (-13311.999999999998_dp)*Hr5(-1,0,1,0,0) + z*(-13311.999999999998_dp)*Hr5(-1,0,1,0,0) + &
                        z**2*(-6655.999999999999_dp)*Hr5(-1,0,1,0,0) + (-19200._dp)*Hr5(-1,0,1,0,1) + z*(-19200._dp)*Hr5(-1,0,1,0,1) &
                        + z**2*(-9600._dp)*Hr5(-1,0,1,0,1) + (-19200._dp)*Hr5(-1,0,1,1,0) + z*(-19200._dp)*Hr5(-1,0,1,1,0) + &
                        z**2*(-9600._dp)*Hr5(-1,0,1,1,0) + (-11093.333333333334_dp)*Hr5(-1,0,1,1,1) + &
                        z*(-11093.333333333334_dp)*Hr5(-1,0,1,1,1) + z**2*(-5546.666666666667_dp)*Hr5(-1,0,1,1,1) + &
                        (-52053.33333333333_dp)*Hr5(0,-1,-1,-1,0) + z*(-47634.96296296296_dp)*Hr5(0,-1,-1,-1,0) + &
                        z**2*(-32366.617283950614_dp)*Hr5(0,-1,-1,-1,0) + z*(-5216.395061728395_dp)*Hr5(0,-1,-1,0,0) + &
                        z**2*(90566.32098765433_dp)*Hr5(0,-1,-1,0,0) + (106211.55555555556_dp)*Hr5(0,-1,-1,0,0) + &
                        z*(20363.061728395063_dp)*Hr5(0,-1,-1,0,1) + z**2*(59640.09876543209_dp)*Hr5(0,-1,-1,0,1) + &
                        (74808.88888888889_dp)*Hr5(0,-1,-1,0,1) + z*(853.3333333333334_dp)*Hr5(0,-1,0,-1,0) + &
                        z**2*(42162.567901234564_dp)*Hr5(0,-1,0,-1,0) + (51000.88888888888_dp)*Hr5(0,-1,0,-1,0) + &
                        (-73187.55555555555_dp)*Hr5(0,-1,0,0,0) + z**2*(-66908.44444444444_dp)*Hr5(0,-1,0,0,0) + &
                        z*(22991.012345679013_dp)*Hr5(0,-1,0,0,0) + (-83370.66666666667_dp)*Hr5(0,-1,0,0,1) + &
                        z**2*(-70147.95061728395_dp)*Hr5(0,-1,0,0,1) + z*(-8917.333333333334_dp)*Hr5(0,-1,0,0,1) + &
                        (-34958.22222222222_dp)*Hr5(0,-1,0,1,0) + z**2*(-28760.493827160495_dp)*Hr5(0,-1,0,1,0) + &
                        z*(-8624.987654320988_dp)*Hr5(0,-1,0,1,0) + z*(-50454.12345679013_dp)*Hr5(0,-1,0,1,1) + &
                        (-33564.44444444445_dp)*Hr5(0,-1,0,1,1) + z**2*(-19595.061728395063_dp)*Hr5(0,-1,0,1,1) + &
                        z*(-73758.02469135802_dp)*Hr5(0,0,-1,-1,0) + nf*z*(-180.14814814814815_dp)*Hr5(0,0,-1,-1,0) + &
                        nf*z**2*(-90.07407407407408_dp)*Hr5(0,0,-1,-1,0) + z**2*(28423.111111111106_dp)*Hr5(0,0,-1,-1,0) + &
                        (29696._dp)*Hr5(0,0,-1,-1,0) + z**2*(-72926.024691358_dp)*Hr5(0,0,-1,0,0) + (-50432._dp)*Hr5(0,0,-1,0,0) + &
                        nf*z**2*(727.7037037037037_dp)*Hr5(0,0,-1,0,0) + nf*z*(1455.4074074074074_dp)*Hr5(0,0,-1,0,0) + &
                        z*(54989.432098765436_dp)*Hr5(0,0,-1,0,0) + z**2*(-37632.79012345679_dp)*Hr5(0,0,-1,0,1) + &
                        (-34304._dp)*Hr5(0,0,-1,0,1) + z*(67471.8024691358_dp)*Hr5(0,0,-1,0,1) + &
                        z*(-36897.185185185175_dp)*Hr5(0,0,0,-1,0) + z**2*(-25446.716049382718_dp)*Hr5(0,0,0,-1,0) + &
                        (-4608._dp)*Hr5(0,0,0,-1,0) + nf*z**2*(789.3333333333333_dp)*Hr5(0,0,0,-1,0) + &
                        nf*z*(3038.814814814815_dp)*Hr5(0,0,0,-1,0) + z**2*(-49943.7037037037_dp)*Hr5(0,0,0,0,0) + &
                        z*(-20685.432098765432_dp)*Hr5(0,0,0,0,0) + nf*z**2*(2785.185185185185_dp)*Hr5(0,0,0,0,0) + &
                        nf*z*(2920.296296296296_dp)*Hr5(0,0,0,0,0) + z**2*(-63375.8024691358_dp)*Hr5(0,0,0,0,1) + &
                        nf*z*(-3055.4074074074074_dp)*Hr5(0,0,0,0,1) + nf*z**2*(1527.7037037037037_dp)*Hr5(0,0,0,0,1) + &
                        z*(39367.90123456789_dp)*Hr5(0,0,0,0,1) + z**2*(-59235.95061728395_dp)*Hr5(0,0,0,1,0) + &
                        (-18432._dp)*Hr5(0,0,0,1,0) + nf*z*(-1668.7407407407404_dp)*Hr5(0,0,0,1,0) + &
                        nf*z**2*(663.7037037037036_dp)*Hr5(0,0,0,1,0) + z*(732.4444444444445_dp)*Hr5(0,0,0,1,0) + &
                        z**2*(-55900.44444444444_dp)*Hr5(0,0,0,1,1) + z*(-42544.98765432099_dp)*Hr5(0,0,0,1,1) + &
                        (-25344._dp)*Hr5(0,0,0,1,1) + nf*z*(-208.59259259259255_dp)*Hr5(0,0,0,1,1) + &
                        nf*z**2*(104.29629629629628_dp)*Hr5(0,0,0,1,1) + z**2*(-59749.92592592592_dp)*Hr5(0,0,1,0,0) + &
                        (-38656._dp)*Hr5(0,0,1,0,0) + z*(-11212.641975308641_dp)*Hr5(0,0,1,0,0) + &
                        nf*z*(-1083.2592592592591_dp)*Hr5(0,0,1,0,0) + nf*z**2*(541.6296296296296_dp)*Hr5(0,0,1,0,0) + &
                        (-43263.99999999999_dp)*Hr5(0,0,1,0,1) + z**2*(-42666.666666666664_dp)*Hr5(0,0,1,0,1) + &
                        z*(-22215.11111111111_dp)*Hr5(0,0,1,0,1) + nf*z**2*(-85.33333333333333_dp)*Hr5(0,0,1,0,1) + &
                        nf*z*(170.66666666666666_dp)*Hr5(0,0,1,0,1) + (-43263.99999999999_dp)*Hr5(0,0,1,1,0) + &
                        z**2*(-42034.567901234564_dp)*Hr5(0,0,1,1,0) + z*(-21772.64197530864_dp)*Hr5(0,0,1,1,0) + &
                        nf*z**2*(-85.33333333333333_dp)*Hr5(0,0,1,1,0) + nf*z*(170.66666666666666_dp)*Hr5(0,0,1,1,0) + &
                        (-40192._dp)*Hr5(0,0,1,1,1) + z**2*(-26462.81481481481_dp)*Hr5(0,0,1,1,1) + &
                        nf*z**2*(-47.407407407407405_dp)*Hr5(0,0,1,1,1) + nf*z*(94.81481481481481_dp)*Hr5(0,0,1,1,1) + &
                        z*(29970.96296296296_dp)*Hr5(0,0,1,1,1) + z**2*(13495.308641975309_dp)*Hr5(0,1,0,-1,0) + &
                        (13738.666666666664_dp)*Hr5(0,1,0,-1,0) + z*(14225.382716049382_dp)*Hr5(0,1,0,-1,0) + &
                        (-44287.99999999999_dp)*Hr5(0,1,0,0,0) + z**2*(-43142.71604938272_dp)*Hr5(0,1,0,0,0) + &
                        z*(-24753.777777777777_dp)*Hr5(0,1,0,0,0) + (-51968._dp)*Hr5(0,1,0,0,1) + &
                        z**2*(-43618.370370370365_dp)*Hr5(0,1,0,0,1) + z*(915.753086419753_dp)*Hr5(0,1,0,0,1) + &
                        (-53617.77777777778_dp)*Hr5(0,1,0,1,0) + z**2*(-39793.77777777778_dp)*Hr5(0,1,0,1,0) + &
                        z*(14250.666666666666_dp)*Hr5(0,1,0,1,0) + (-33422.22222222222_dp)*Hr5(0,1,0,1,1) + &
                        z**2*(-16726.913580246914_dp)*Hr5(0,1,0,1,1) + z*(41247.60493827161_dp)*Hr5(0,1,0,1,1) + &
                        (-53418.66666666666_dp)*Hr5(0,1,1,0,0) + z**2*(-44187.25925925926_dp)*Hr5(0,1,1,0,0) + &
                        z*(-2298.4691358024693_dp)*Hr5(0,1,1,0,0) + (-33991.11111111111_dp)*Hr5(0,1,1,0,1) + &
                        z**2*(-16948.14814814815_dp)*Hr5(0,1,1,0,1) + z*(39983.40740740741_dp)*Hr5(0,1,1,0,1) + &
                        (-32853.33333333333_dp)*Hr5(0,1,1,1,0) + z**2*(-15273.086419753085_dp)*Hr5(0,1,1,1,0) + &
                        z*(40046.61728395062_dp)*Hr5(0,1,1,1,0) + (-10752._dp)*Hr5(0,1,1,1,1) + &
                        z**2*(-2863.4074074074074_dp)*Hr5(0,1,1,1,1) + z*(23646.814814814818_dp)*Hr5(0,1,1,1,1) + &
                        (-9016.888888888889_dp)*Hr5(1,0,-1,-1,0) + z**2*(-4508.444444444444_dp)*Hr5(1,0,-1,-1,0) + &
                        z*(9016.888888888889_dp)*Hr5(1,0,-1,-1,0) + z*(-6637.037037037037_dp)*Hr5(1,0,-1,0,0) + &
                        z**2*(3318.5185185185187_dp)*Hr5(1,0,-1,0,0) + (6637.037037037037_dp)*Hr5(1,0,-1,0,0) + &
                        z*(-4257.185185185185_dp)*Hr5(1,0,-1,0,1) + z**2*(2128.5925925925926_dp)*Hr5(1,0,-1,0,1) + &
                        (4257.185185185185_dp)*Hr5(1,0,-1,0,1) + z*(-10891.06172839506_dp)*Hr5(1,0,0,-1,0) + &
                        z**2*(5445.53086419753_dp)*Hr5(1,0,0,-1,0) + (10891.06172839506_dp)*Hr5(1,0,0,-1,0) + &
                        (-21782.12345679012_dp)*Hr5(1,0,0,0,0) + z**2*(-10891.06172839506_dp)*Hr5(1,0,0,0,0) + &
                        z*(21782.12345679012_dp)*Hr5(1,0,0,0,0) + (-42587.654320987655_dp)*Hr5(1,0,0,0,1) + &
                        z**2*(-21293.827160493827_dp)*Hr5(1,0,0,0,1) + z*(42587.654320987655_dp)*Hr5(1,0,0,0,1) + &
                        (-52423.1111111111_dp)*Hr5(1,0,0,1,0) + z**2*(-26211.55555555555_dp)*Hr5(1,0,0,1,0) + &
                        z*(52423.1111111111_dp)*Hr5(1,0,0,1,0) + (-49177.28395061728_dp)*Hr5(1,0,0,1,1) + &
                        z**2*(-24588.64197530864_dp)*Hr5(1,0,0,1,1) + z*(49177.28395061728_dp)*Hr5(1,0,0,1,1) + &
                        (-53886.41975308642_dp)*Hr5(1,0,1,0,0) + z**2*(-26943.20987654321_dp)*Hr5(1,0,1,0,0) + &
                        z*(53886.41975308642_dp)*Hr5(1,0,1,0,0) + (-51525.530864197535_dp)*Hr5(1,0,1,0,1) + &
                        z**2*(-25762.765432098768_dp)*Hr5(1,0,1,0,1) + z*(51525.530864197535_dp)*Hr5(1,0,1,0,1) + &
                        (-49882.074074074066_dp)*Hr5(1,0,1,1,0) + z**2*(-24941.037037037033_dp)*Hr5(1,0,1,1,0) + &
                        z*(49882.074074074066_dp)*Hr5(1,0,1,1,0) + (-23883.851851851847_dp)*Hr5(1,0,1,1,1) + &
                        z**2*(-11941.925925925923_dp)*Hr5(1,0,1,1,1) + z*(23883.851851851847_dp)*Hr5(1,0,1,1,1) + &
                        z*(-12765.234567901232_dp)*Hr5(1,1,0,-1,0) + z**2*(6382.617283950616_dp)*Hr5(1,1,0,-1,0) + &
                        (12765.234567901232_dp)*Hr5(1,1,0,-1,0) + (-48304.98765432098_dp)*Hr5(1,1,0,0,0) + &
                        z**2*(-24152.49382716049_dp)*Hr5(1,1,0,0,0) + z*(48304.98765432098_dp)*Hr5(1,1,0,0,0) + &
                        (-53358.61728395062_dp)*Hr5(1,1,0,0,1) + z**2*(-26679.30864197531_dp)*Hr5(1,1,0,0,1) + &
                        z*(53358.61728395062_dp)*Hr5(1,1,0,0,1) + (-50703.8024691358_dp)*Hr5(1,1,0,1,0) + &
                        z**2*(-25351.9012345679_dp)*Hr5(1,1,0,1,0) + z*(50703.8024691358_dp)*Hr5(1,1,0,1,0) + &
                        (-25527.308641975305_dp)*Hr5(1,1,0,1,1) + z**2*(-12763.654320987653_dp)*Hr5(1,1,0,1,1) + &
                        z*(25527.308641975305_dp)*Hr5(1,1,0,1,1) + (-51282.17283950617_dp)*Hr5(1,1,1,0,0) + &
                        z**2*(-25641.086419753086_dp)*Hr5(1,1,1,0,0) + z*(51282.17283950617_dp)*Hr5(1,1,1,0,0) + &
                        (-24705.58024691358_dp)*Hr5(1,1,1,0,1) + z**2*(-12352.79012345679_dp)*Hr5(1,1,1,0,1) + &
                        z*(24705.58024691358_dp)*Hr5(1,1,1,0,1) + (-21418.666666666668_dp)*Hr5(1,1,1,1,0) + &
                        z**2*(-10709.333333333334_dp)*Hr5(1,1,1,1,0) + z*(21418.666666666668_dp)*Hr5(1,1,1,1,0) + &
                        (-5214.814814814815_dp)*Hr5(1,1,1,1,1) + z**2*(-2607.4074074074074_dp)*Hr5(1,1,1,1,1) + &
                        z*(5214.814814814815_dp)*Hr5(1,1,1,1,1))/z

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
                        Ibar_select = (z*(-447914.31258482905_dp) + z**3*(-112932.47676040782_dp) + nf**2*(-281.3550710804751_dp) + &
                        nf**2*z**2*(-217.95847019424687_dp) + nf**3*z*(-2.633744855967078_dp) + nf**3*z**2*(1.316872427983539_dp) + &
                        nf**3*(2.633744855967078_dp) + nf**2*z**3*(9.803383630544122_dp) + nf**2*z*(350.51468394790066_dp) + &
                        nf*z**2*(764.6242235143792_dp) + nf*z**3*(1622.310412411921_dp) + nf*z*(1812.6505450228306_dp) + &
                        nf*(4221.881914702317_dp) + z**2*(4979.742206068593_dp) + (66792.8649909556_dp) + &
                        (nf*((-995.5709591749068_dp) + z*(-995.5709591749068_dp) + z**2*(-497.7854795874534_dp)) + &
                        z**3*(6997.590136028327_dp) + z**2*(12353.36329893726_dp) + z*(78022.19418703933_dp) + &
                        (88600.89378718419_dp))*Hr1(-1) + z**2*(-8417.452149705636_dp)*Hr1(0) + nf*z*(-5399.828100461422_dp)*Hr1(0) &
                        + nf*z**3*(-1295.070873342478_dp)*Hr1(0) + nf**2*(-86.82578875171467_dp)*Hr1(0) + &
                        nf**2*z**2*(-11.225942880758922_dp)*Hr1(0) + nf**2*z**3*(5.969821673525377_dp)*Hr1(0) + &
                        nf**2*z*(95.1432437862092_dp)*Hr1(0) + nf*z**2*(2033.475680200342_dp)*Hr1(0) + &
                        nf*(3490.4400055566853_dp)*Hr1(0) + (9396.366276418172_dp)*Hr1(0) + z**3*(10058.80702479694_dp)*Hr1(0) + &
                        z*(151071.50260846564_dp)*Hr1(0) + z*(-90931.10675691547_dp)*Hr1(1) + z**3*(-37367.940118225946_dp)*Hr1(1) + &
                        nf*z*(-879.5818077260777_dp)*Hr1(1) + nf**2*z*(-115.97256515775035_dp)*Hr1(1) + &
                        nf**2*z**3*(-11.939643347050755_dp)*Hr1(1) + nf**3*(-1.5802469135802466_dp)*Hr1(1) + &
                        nf**3*z**2*(-0.7901234567901233_dp)*Hr1(1) + nf**3*z*(1.5802469135802466_dp)*Hr1(1) + &
                        nf**2*z**2*(39.286694101508914_dp)*Hr1(1) + nf**2*(67.33607681755831_dp)*Hr1(1) + &
                        nf*z**3*(448.99497027892085_dp)*Hr1(1) + nf*z**2*(654.7044841099523_dp)*Hr1(1) + &
                        nf*(1271.14650822905_dp)*Hr1(1) + z**2*(24296.104470719005_dp)*Hr1(1) + (90069.09865561892_dp)*Hr1(1) + &
                        (-24326.936683894415_dp)*Hr2(-1,-1) + z*(-24326.936683894415_dp)*Hr2(-1,-1) + &
                        z**2*(-12163.468341947208_dp)*Hr2(-1,-1) + z*(-57525.84494223027_dp)*Hr2(-1,0) + &
                        (-55593.071279678836_dp)*Hr2(-1,0) + z**3*(-9623.703703703704_dp)*Hr2(-1,0) + &
                        z**2*(-4483.087080168636_dp)*Hr2(-1,0) + nf**2*(-14.22222222222222_dp)*Hr2(-1,0) + &
                        nf**2*z*(-14.22222222222222_dp)*Hr2(-1,0) + nf**2*z**2*(-7.11111111111111_dp)*Hr2(-1,0) + &
                        nf*z**3*(218.07407407407405_dp)*Hr2(-1,0) + nf*z**2*(908.641975308642_dp)*Hr2(-1,0) + &
                        nf*(2149.9259259259256_dp)*Hr2(-1,0) + nf*z*(2662.716049382716_dp)*Hr2(-1,0) + &
                        z*(-19325.687270180264_dp)*Hr2(0,-1) + (20992.77040808498_dp)*Hr2(0,-1) + &
                        z**2*(21826.311977037338_dp)*Hr2(0,-1) + nf*z**3*(-286.9026063100137_dp)*Hr2(0,0) + &
                        nf**2*z**2*(-28.971193415637856_dp)*Hr2(0,0) + nf**2*z*(73.7448559670782_dp)*Hr2(0,0) + &
                        (247.6546816697189_dp)*Hr2(0,0) + nf*(609.9753086419753_dp)*Hr2(0,0) + &
                        nf*z**2*(1571.3555740379593_dp)*Hr2(0,0) + nf*z*(1689.1135895399632_dp)*Hr2(0,0) + &
                        z**2*(2288.5429946606237_dp)*Hr2(0,0) + z**3*(18724.872427983537_dp)*Hr2(0,0) + &
                        z*(34165.3732118983_dp)*Hr2(0,0) + nf*z*(-2021.3113854595333_dp)*Hr2(0,1) + &
                        nf*z**3*(-245.46502057613168_dp)*Hr2(0,1) + (-102.09480649656687_dp)*Hr2(0,1) + &
                        nf**2*z**2*(-2.8971193415637857_dp)*Hr2(0,1) + nf**2*(21.596707818930042_dp)*Hr2(0,1) + &
                        nf**2*z*(24.757201646090536_dp)*Hr2(0,1) + nf*(1259.5445816186557_dp)*Hr2(0,1) + &
                        nf*z**2*(1597.190672153635_dp)*Hr2(0,1) + z**2*(9482.476979626963_dp)*Hr2(0,1) + &
                        z**3*(15590.189300411523_dp)*Hr2(0,1) + z*(100143.0173928455_dp)*Hr2(0,1) + &
                        (-52690.269872947516_dp)*Hr2(1,0) + z**2*(-6011.505306844132_dp)*Hr2(1,0) + &
                        nf*z*(-1829.048010973937_dp)*Hr2(1,0) + nf*z**3*(-77.6076817558299_dp)*Hr2(1,0) + &
                        nf**2*z*(-15.275720164609051_dp)*Hr2(1,0) + nf**2*z**2*(7.637860082304526_dp)*Hr2(1,0) + &
                        nf**2*(15.275720164609051_dp)*Hr2(1,0) + nf*z**2*(947.9725651577503_dp)*Hr2(1,0) + &
                        nf*(1246.3758573388202_dp)*Hr2(1,0) + z**3*(4265.61316872428_dp)*Hr2(1,0) + &
                        z*(55651.257527268506_dp)*Hr2(1,0) + (z*(-18747.353116088358_dp) + z**3*(-9298.699588477366_dp) + &
                        nf**2*(z*(-63.7366255144033_dp) + z**2*(31.86831275720165_dp) + (63.7366255144033_dp)) + &
                        nf*((-1606.9355281207133_dp) + z**3*(232.47187928669413_dp) + z**2*(399.18792866941016_dp) + &
                        z*(1148.4883401920438_dp)) + z**2*(7091.536640348704_dp) + (20938.10208728177_dp))*Hr2(1,1) + &
                        nf*(-565.7283950617283_dp)*Hr3(-1,-1,0) + nf*z*(-565.7283950617283_dp)*Hr3(-1,-1,0) + &
                        nf*z**2*(-282.86419753086415_dp)*Hr3(-1,-1,0) + z**2*(2207.341563786008_dp)*Hr3(-1,-1,0) + &
                        z**3*(3501.827160493827_dp)*Hr3(-1,-1,0) + z*(31118.748971193414_dp)*Hr3(-1,-1,0) + &
                        (36366.22222222222_dp)*Hr3(-1,-1,0) + (-30488.493827160495_dp)*Hr3(-1,0,0) + &
                        z*(-25961.086419753086_dp)*Hr3(-1,0,0) + z**3*(-3002.469135802469_dp)*Hr3(-1,0,0) + &
                        z**2*(-1539.5555555555554_dp)*Hr3(-1,0,0) + nf*z**2*(222.02469135802465_dp)*Hr3(-1,0,0) + &
                        nf*(444.0493827160493_dp)*Hr3(-1,0,0) + nf*z*(444.0493827160493_dp)*Hr3(-1,0,0) + &
                        (-24610.765432098764_dp)*Hr3(-1,0,1) + z*(-20803.42386831276_dp)*Hr3(-1,0,1) + &
                        z**3*(-2503.1111111111113_dp)*Hr3(-1,0,1) + z**2*(-871.7695473251028_dp)*Hr3(-1,0,1) + &
                        nf*z**2*(161.18518518518516_dp)*Hr3(-1,0,1) + nf*(322.3703703703703_dp)*Hr3(-1,0,1) + &
                        nf*z*(322.3703703703703_dp)*Hr3(-1,0,1) + (-27335.111111111106_dp)*Hr3(0,-1,0) + &
                        z**2*(-5137.909465020576_dp)*Hr3(0,-1,0) + z**3*(-1965.827160493827_dp)*Hr3(0,-1,0) + &
                        nf*z*(-357.1358024691358_dp)*Hr3(0,-1,0) + nf*z**2*(205.4320987654321_dp)*Hr3(0,-1,0) + &
                        nf*(654.2222222222222_dp)*Hr3(0,-1,0) + z*(5761.0534979423865_dp)*Hr3(0,-1,0) + &
                        z**2*(-2784.9218106995886_dp)*Hr3(0,0,0) + nf*z*(-1071.4074074074072_dp)*Hr3(0,0,0) + &
                        nf**2*z**2*(-18.962962962962962_dp)*Hr3(0,0,0) + nf**2*z*(37.925925925925924_dp)*Hr3(0,0,0) + &
                        nf*z**2*(347.6543209876543_dp)*Hr3(0,0,0) + z**3*(3034.074074074074_dp)*Hr3(0,0,0) + &
                        z*(19866.337448559672_dp)*Hr3(0,0,0) + (-18816._dp)*Hr3(0,0,1) + z*(-910.7489711934155_dp)*Hr3(0,0,1) + &
                        nf*z**2*(-244.41152263374482_dp)*Hr3(0,0,1) + nf**2*z*(-25.283950617283946_dp)*Hr3(0,0,1) + &
                        nf**2*z**2*(12.641975308641973_dp)*Hr3(0,0,1) + nf*(407.7037037037037_dp)*Hr3(0,0,1) + &
                        nf*z*(738.5020576131687_dp)*Hr3(0,0,1) + z**2*(5388.641975308642_dp)*Hr3(0,0,1) + &
                        z**3*(6377.876543209875_dp)*Hr3(0,0,1) + (-38656._dp)*Hr3(0,1,0) + nf*z*(-321.31687242798347_dp)*Hr3(0,1,0) &
                        + nf*z**2*(23.176954732510286_dp)*Hr3(0,1,0) + nf*(555.1934156378601_dp)*Hr3(0,1,0) + &
                        z**2*(2185.4814814814813_dp)*Hr3(0,1,0) + z**3*(5720.493827160493_dp)*Hr3(0,1,0) + &
                        z*(23355.522633744855_dp)*Hr3(0,1,0) + (-26609.777777777777_dp)*Hr3(0,1,1) + &
                        nf*z*(-1283.1604938271603_dp)*Hr3(0,1,1) + nf*z**2*(-515.1604938271605_dp)*Hr3(0,1,1) + &
                        nf*(-52.674897119341566_dp)*Hr3(0,1,1) + z**3*(7667.358024691357_dp)*Hr3(0,1,1) + &
                        z**2*(15430.58436213992_dp)*Hr3(0,1,1) + z*(55673.152263374475_dp)*Hr3(0,1,1) + &
                        (-30858.271604938273_dp)*Hr3(1,0,0) + z**2*(-503.440329218107_dp)*Hr3(1,0,0) + &
                        nf*z*(-298.6666666666667_dp)*Hr3(1,0,0) + nf*z**2*(149.33333333333334_dp)*Hr3(1,0,0) + &
                        nf*(298.6666666666667_dp)*Hr3(1,0,0) + z**3*(3128.8888888888882_dp)*Hr3(1,0,0) + &
                        z*(26636.37860082304_dp)*Hr3(1,0,0) + (-39194.86419753086_dp)*Hr3(1,0,1) + &
                        nf*z*(-38.97942386831275_dp)*Hr3(1,0,1) + nf*z**2*(19.489711934156375_dp)*Hr3(1,0,1) + &
                        nf*(38.97942386831275_dp)*Hr3(1,0,1) + z**2*(2261.59670781893_dp)*Hr3(1,0,1) + &
                        z**3*(4323.555555555555_dp)*Hr3(1,0,1) + z*(33394.30452674897_dp)*Hr3(1,0,1) + &
                        (-39194.86419753086_dp)*Hr3(1,1,0) + nf*z*(-38.97942386831275_dp)*Hr3(1,1,0) + &
                        nf*z**2*(19.489711934156375_dp)*Hr3(1,1,0) + nf*(38.97942386831275_dp)*Hr3(1,1,0) + &
                        z**2*(2261.59670781893_dp)*Hr3(1,1,0) + z**3*(4323.555555555555_dp)*Hr3(1,1,0) + &
                        z*(33394.30452674897_dp)*Hr3(1,1,0) + (-7376.592592592591_dp)*Hr3(1,1,1) + &
                        nf*(-1144.0987654320988_dp)*Hr3(1,1,1) + nf*z**2*(-572.0493827160494_dp)*Hr3(1,1,1) + &
                        nf*z*(1144.0987654320988_dp)*Hr3(1,1,1) + z**3*(2085.925925925926_dp)*Hr3(1,1,1) + &
                        z*(4800.79012345679_dp)*Hr3(1,1,1) + z**2*(10298.469135802468_dp)*Hr3(1,1,1) + &
                        (-12266.930041152264_dp)*Hr4(-1,-1,-1,0) + z*(-12266.930041152264_dp)*Hr4(-1,-1,-1,0) + &
                        z**2*(-6133.465020576132_dp)*Hr4(-1,-1,-1,0) + z**2*(5230.617283950617_dp)*Hr4(-1,-1,0,0) + &
                        (10461.234567901234_dp)*Hr4(-1,-1,0,0) + z*(10461.234567901234_dp)*Hr4(-1,-1,0,0) + &
                        z**2*(4327.769547325102_dp)*Hr4(-1,-1,0,1) + (8655.539094650205_dp)*Hr4(-1,-1,0,1) + &
                        z*(8655.539094650205_dp)*Hr4(-1,-1,0,1) + z**2*(3066.732510288066_dp)*Hr4(-1,0,-1,0) + &
                        (6133.465020576132_dp)*Hr4(-1,0,-1,0) + z*(6133.465020576132_dp)*Hr4(-1,0,-1,0) + &
                        (-4327.769547325102_dp)*Hr4(-1,0,0,0) + z*(-4327.769547325102_dp)*Hr4(-1,0,0,0) + &
                        z**2*(-2163.884773662551_dp)*Hr4(-1,0,0,0) + (-6849.843621399176_dp)*Hr4(-1,0,0,1) + &
                        z*(-6849.843621399176_dp)*Hr4(-1,0,0,1) + z**2*(-3424.921810699588_dp)*Hr4(-1,0,0,1) + &
                        (-2522.074074074074_dp)*Hr4(-1,0,1,0) + z*(-2522.074074074074_dp)*Hr4(-1,0,1,0) + &
                        z**2*(-1261.037037037037_dp)*Hr4(-1,0,1,0) + (-5044.148148148148_dp)*Hr4(-1,0,1,1) + &
                        z*(-5044.148148148148_dp)*Hr4(-1,0,1,1) + z**2*(-2522.074074074074_dp)*Hr4(-1,0,1,1) + &
                        z*(-9624.75720164609_dp)*Hr4(0,-1,-1,0) + (10505.481481481482_dp)*Hr4(0,-1,-1,0) + &
                        z**2*(10945.843621399177_dp)*Hr4(0,-1,-1,0) + z**2*(-9370.864197530864_dp)*Hr4(0,-1,0,0) + &
                        (-9007.407407407407_dp)*Hr4(0,-1,0,0) + z*(8280.493827160491_dp)*Hr4(0,-1,0,0) + &
                        z**2*(-7795.88477366255_dp)*Hr4(0,-1,0,1) + (-7509.333333333332_dp)*Hr4(0,-1,0,1) + &
                        z*(6936.230452674897_dp)*Hr4(0,-1,0,1) + z**2*(-9890.23868312757_dp)*Hr4(0,0,-1,0) + &
                        (-4608._dp)*Hr4(0,0,-1,0) + z*(-4022.2551440329216_dp)*Hr4(0,0,-1,0) + &
                        nf*z**2*(170.66666666666666_dp)*Hr4(0,0,-1,0) + nf*z*(341.3333333333333_dp)*Hr4(0,0,-1,0) + &
                        z**2*(-10842.600823045268_dp)*Hr4(0,0,0,0) + z*(-3760.987654320987_dp)*Hr4(0,0,0,0) + &
                        nf*z*(530.9629629629629_dp)*Hr4(0,0,0,0) + nf*z**2*(587.8518518518517_dp)*Hr4(0,0,0,0) + &
                        z**2*(-18967.17695473251_dp)*Hr4(0,0,0,1) + (-4608._dp)*Hr4(0,0,0,1) + &
                        nf*z*(-771.1604938271605_dp)*Hr4(0,0,0,1) + nf*z**2*(385.58024691358025_dp)*Hr4(0,0,0,1) + &
                        z*(10233.67901234568_dp)*Hr4(0,0,0,1) + z**2*(-19638.25514403292_dp)*Hr4(0,0,1,0) + &
                        (-12117.333333333334_dp)*Hr4(0,0,1,0) + z*(-4622.748971193416_dp)*Hr4(0,0,1,0) + &
                        nf*z*(-278.12345679012344_dp)*Hr4(0,0,1,0) + nf*z**2*(139.06172839506172_dp)*Hr4(0,0,1,0) + &
                        z**2*(-20949.860082304527_dp)*Hr4(0,0,1,1) + z*(-16373.465020576132_dp)*Hr4(0,0,1,1) + &
                        (-15872._dp)*Hr4(0,0,1,1) + nf*z*(-63.209876543209866_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(31.604938271604933_dp)*Hr4(0,0,1,1) + (-11908.74074074074_dp)*Hr4(0,1,0,0) + &
                        z**2*(-11274.534979423866_dp)*Hr4(0,1,0,0) + z*(-5610.930041152263_dp)*Hr4(0,1,0,0) + &
                        (-20536.88888888889_dp)*Hr4(0,1,0,1) + z**2*(-17636.60905349794_dp)*Hr4(0,1,0,1) + &
                        z*(-3638.781893004115_dp)*Hr4(0,1,0,1) + (-20536.88888888889_dp)*Hr4(0,1,1,0) + &
                        z**2*(-17636.60905349794_dp)*Hr4(0,1,1,0) + z*(-3638.781893004115_dp)*Hr4(0,1,1,0) + &
                        (-21390.22222222222_dp)*Hr4(0,1,1,1) + z**2*(-14361.283950617284_dp)*Hr4(0,1,1,1) + &
                        z*(9949.234567901232_dp)*Hr4(0,1,1,1) + (-5044.148148148148_dp)*Hr4(1,0,0,0) + &
                        z**2*(-2522.074074074074_dp)*Hr4(1,0,0,0) + z*(5044.148148148148_dp)*Hr4(1,0,0,0) + &
                        (-13295.14403292181_dp)*Hr4(1,0,0,1) + z**2*(-6647.572016460905_dp)*Hr4(1,0,0,1) + &
                        z*(13295.14403292181_dp)*Hr4(1,0,0,1) + (-16501.991769547323_dp)*Hr4(1,0,1,0) + &
                        z**2*(-8250.995884773662_dp)*Hr4(1,0,1,0) + z*(16501.991769547323_dp)*Hr4(1,0,1,0) + &
                        (-20705.448559670782_dp)*Hr4(1,0,1,1) + z**2*(-10352.724279835391_dp)*Hr4(1,0,1,1) + &
                        z*(20705.448559670782_dp)*Hr4(1,0,1,1) + (-13295.14403292181_dp)*Hr4(1,1,0,0) + &
                        z**2*(-6647.572016460905_dp)*Hr4(1,1,0,0) + z*(13295.14403292181_dp)*Hr4(1,1,0,0) + &
                        (-20705.448559670782_dp)*Hr4(1,1,0,1) + z**2*(-10352.724279835391_dp)*Hr4(1,1,0,1) + &
                        z*(20705.448559670782_dp)*Hr4(1,1,0,1) + (-20705.448559670782_dp)*Hr4(1,1,1,0) + &
                        z**2*(-10352.724279835391_dp)*Hr4(1,1,1,0) + z*(20705.448559670782_dp)*Hr4(1,1,1,0) + &
                        (-16813.827160493827_dp)*Hr4(1,1,1,1) + z**2*(-8406.913580246914_dp)*Hr4(1,1,1,1) + &
                        z*(16813.827160493827_dp)*Hr4(1,1,1,1))/z

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
                        Ibar_select = ((-17380.831746559503_dp) + nf*z*(-1729.9279179591092_dp) + nf*z**3*(-210.63374485596708_dp) + &
                        nf**2*(-17.009602194787384_dp) + nf**2*z**2*(-14.024691358024691_dp) + nf*z**2*(-2.4945311899723492_dp) + &
                        nf**2*z**3*(2.1947873799725657_dp) + nf**2*z*(30.02469135802469_dp) + nf*(1184.6374326205594_dp) + &
                        z**2*(2673.2903851206215_dp) + z**3*(3853.3179297139823_dp) + z*(34783.1896098702_dp) + &
                        (-11176.300540915414_dp)*Hr1(0) + nf*z**3*(-39.59396433470508_dp)*Hr1(0) + &
                        nf**2*(-6.584362139917697_dp)*Hr1(0) + nf**2*z*(-0.7901234567901235_dp)*Hr1(0) + &
                        nf**2*z**2*(1.580246913580247_dp)*Hr1(0) + nf*z**2*(238.32501348389837_dp)*Hr1(0) + &
                        nf*z*(366.54338867006334_dp)*Hr1(0) + z**2*(526.2355267446983_dp)*Hr1(0) + nf*(597.9259259259259_dp)*Hr1(0) &
                        + z*(1054.2750234410257_dp)*Hr1(0) + z**3*(1812.148148148148_dp)*Hr1(0) + (-19772.228174545202_dp)*Hr1(1) + &
                        z**2*(-2080.8795193713668_dp)*Hr1(1) + nf*z*(-862.551440329218_dp)*Hr1(1) + &
                        nf*z**3*(-73.5692729766804_dp)*Hr1(1) + nf**2*(-5.530864197530864_dp)*Hr1(1) + &
                        nf**2*z**2*(-2.765432098765432_dp)*Hr1(1) + nf**2*z*(5.530864197530864_dp)*Hr1(1) + &
                        nf*z**2*(244.67489711934158_dp)*Hr1(1) + nf*(752.6803840877916_dp)*Hr1(1) + &
                        z**3*(3006.8148148148152_dp)*Hr1(1) + z*(20177.956569606937_dp)*Hr1(1) + (-3456._dp)*Hr2(0,0) + &
                        z**2*(-2207.341563786008_dp)*Hr2(0,0) + nf*z*(-176.72427983539097_dp)*Hr2(0,0) + &
                        nf**2*z**2*(-3.160493827160494_dp)*Hr2(0,0) + nf**2*z*(6.320987654320988_dp)*Hr2(0,0) + &
                        nf*(101.92592592592592_dp)*Hr2(0,0) + nf*z**2*(117.99176954732512_dp)*Hr2(0,0) + &
                        z**3*(210.17283950617286_dp)*Hr2(0,0) + z*(2121.0864197530864_dp)*Hr2(0,0) + &
                        (-8579.555555555555_dp)*Hr2(0,1) + z**2*(-2389.3333333333335_dp)*Hr2(0,1) + &
                        z*(-1765.1358024691358_dp)*Hr2(0,1) + nf*z**2*(124.31275720164611_dp)*Hr2(0,1) + &
                        nf*z*(125.89300411522635_dp)*Hr2(0,1) + nf*(220.70781893004119_dp)*Hr2(0,1) + &
                        z**3*(605.2345679012346_dp)*Hr2(0,1) + (-6083.16049382716_dp)*Hr2(1,0) + &
                        z**2*(-675.358024691358_dp)*Hr2(1,0) + nf*z*(-122.99588477366255_dp)*Hr2(1,0) + &
                        nf*z**2*(61.49794238683128_dp)*Hr2(1,0) + nf*(122.99588477366255_dp)*Hr2(1,0) + &
                        z**3*(523.0617283950618_dp)*Hr2(1,0) + z*(5379.555555555556_dp)*Hr2(1,0) + ((-12166.32098765432_dp) + &
                        z**2*(-1350.716049382716_dp) + nf*(z*(-245.9917695473251_dp) + z**2*(122.99588477366255_dp) + &
                        (245.9917695473251_dp)) + z**3*(1046.1234567901236_dp) + z*(10759.111111111111_dp))*Hr2(1,1) + &
                        ((-384.00000000000006_dp) + z*((-313.41563786008237_dp) + nf*(44.24691358024691_dp)) + &
                        z**2*((-1727.7366255144036_dp) + nf*(63.20987654320988_dp)))*Hr3(0,0,0) + &
                        z**2*(-2913.9753086419755_dp)*Hr3(0,0,1) + (-1322.6666666666667_dp)*Hr3(0,0,1) + &
                        nf*z*(-82.17283950617283_dp)*Hr3(0,0,1) + nf*z**2*(41.08641975308642_dp)*Hr3(0,0,1) + &
                        z*(1120.3950617283951_dp)*Hr3(0,0,1) + z**2*(-1675.061728395062_dp)*Hr3(0,1,0) + &
                        (-1569.1851851851852_dp)*Hr3(0,1,0) + z*(-1357.4320987654323_dp)*Hr3(0,1,0) + &
                        z**2*(-3350.123456790124_dp)*Hr3(0,1,1) + (-3138.3703703703704_dp)*Hr3(0,1,1) + &
                        z*(-2714.8641975308647_dp)*Hr3(0,1,1) + (-664.2304526748972_dp)*Hr3(1,0,0) + &
                        z**2*(-332.1152263374486_dp)*Hr3(1,0,0) + z*(664.2304526748972_dp)*Hr3(1,0,0) + &
                        (-1992.6913580246915_dp)*Hr3(1,0,1) + z**2*(-996.3456790123457_dp)*Hr3(1,0,1) + &
                        z*(1992.6913580246915_dp)*Hr3(1,0,1) + (-1992.6913580246915_dp)*Hr3(1,1,0) + &
                        z**2*(-996.3456790123457_dp)*Hr3(1,1,0) + z*(1992.6913580246915_dp)*Hr3(1,1,0) + &
                        (-3985.382716049383_dp)*Hr3(1,1,1) + z**2*(-1992.6913580246915_dp)*Hr3(1,1,1) + &
                        z*(3985.382716049383_dp)*Hr3(1,1,1))/z

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
                        Ibar_select = (z*(-1.2332664974292908e7_dp) + z**3*(-2.836421951404494e6_dp) + nf**2*(-8223.910681080284_dp) &
                        + nf**2*z**2*(-5907.744228671624_dp) + nf**2*z**3*(-624.7427780405193_dp) + nf**3*z*(-321.6797175400653_dp) &
                        + nf**3*z**3*(-8.691358024691358_dp) + nf**4*(-1.755829903978052_dp) + nf**4*z**2*(-0.877914951989026_dp) + &
                        nf**4*z*(1.755829903978052_dp) + nf**3*z**2*(184.4774162471891_dp) + nf**3*(235.58958298143708_dp) + &
                        nf*(1459.0217798840895_dp) + nf**2*z*(3385.586185664163_dp) + nf*z**2*(12085.617597418563_dp) + &
                        nf*z**3*(120403.31093805254_dp) + z**2*(170225.68315235176_dp) + nf*z*(356142.0329395592_dp) + &
                        (1.7213019903156676e6_dp) + (nf*(z*(-105404.39522932813_dp) + (-98001.77975481717_dp) + &
                        z**2*(-32753.231513645558_dp) + z**3*(-7712.329406980315_dp)) + nf**2*(z**2*(337.7056408515404_dp) + &
                        (675.4112817030808_dp) + z*(675.4112817030808_dp)) + z**3*(376468.4836307601_dp) + &
                        z**2*(413615.40813680243_dp) + z*(2.9943032989432043e6_dp) + (3.110204629524099e6_dp))*Hr1(-1) + &
                        nf*z*(-278757.90387485403_dp)*Hr1(0) + z**3*(-217938.94887629937_dp)*Hr1(0) + &
                        nf*z**3*(-36545.6092192487_dp)*Hr1(0) + nf**2*(-4635.014059083678_dp)*Hr1(0) + &
                        nf**2*z**2*(-2129.3199862924075_dp)*Hr1(0) + nf**3*z*(-76.41439147710997_dp)*Hr1(0) + &
                        nf**3*z**2*(-13.940952409593153_dp)*Hr1(0) + nf**3*z**3*(-2.633744855967078_dp)*Hr1(0) + &
                        nf**3*(64.61454046639231_dp)*Hr1(0) + nf**2*z**3*(996.7797591830515_dp)*Hr1(0) + &
                        nf**2*z*(3292.6036769854923_dp)*Hr1(0) + nf*z**2*(17586.531293702046_dp)*Hr1(0) + &
                        z**2*(29805.126651826668_dp)*Hr1(0) + nf*(63493.840197599224_dp)*Hr1(0) + (363910.4787446698_dp)*Hr1(0) + &
                        z*(2.262103018823811e6_dp)*Hr1(0) + z*(-3.576360892768053e6_dp)*Hr1(1) + &
                        z**3*(-1.088057496140819e6_dp)*Hr1(1) + nf*(-19333.469893539466_dp)*Hr1(1) + &
                        nf*z**2*(-15814.079496329801_dp)*Hr1(1) + nf**2*z**2*(-1314.073851620496_dp)*Hr1(1) + &
                        nf**2*z**3*(-437.52842554488643_dp)*Hr1(1) + nf**2*(-83.21369927818131_dp)*Hr1(1) + &
                        nf**3*(-71.82807498856882_dp)*Hr1(1) + nf**3*z**2*(-54.85505258344765_dp)*Hr1(1) + &
                        nf**4*z*(-1.0534979423868314_dp)*Hr1(1) + nf**4*z**2*(0.5267489711934157_dp)*Hr1(1) + &
                        nf**4*(1.0534979423868314_dp)*Hr1(1) + nf**3*z**3*(9.305898491083676_dp)*Hr1(1) + &
                        nf**3*z*(123.31778692272519_dp)*Hr1(1) + nf**2*z*(553.7370951019886_dp)*Hr1(1) + &
                        nf*z*(24239.348139012094_dp)*Hr1(1) + nf*z**3*(28079.87651141408_dp)*Hr1(1) + &
                        z**2*(1.067861780891705e6_dp)*Hr1(1) + (3.509591427938627e6_dp)*Hr1(1) + (-1.35711273690269e6_dp)*Hr2(-1,-1) &
                        + z*(-1.2091461114037354e6_dp)*Hr2(-1,-1) + z**2*(-239400.40057195962_dp)*Hr2(-1,-1) + &
                        z**3*(-96039.23856827388_dp)*Hr2(-1,-1) + nf*z**2*(10927.308288199427_dp)*Hr2(-1,-1) + &
                        nf*(21854.616576398854_dp)*Hr2(-1,-1) + nf*z*(21854.616576398854_dp)*Hr2(-1,-1) + &
                        z*(-1.0971988794695265e6_dp)*Hr2(-1,0) + (-998366.2074899839_dp)*Hr2(-1,0) + &
                        z**3*(-197333.69433238884_dp)*Hr2(-1,0) + z**2*(-153129.04400694996_dp)*Hr2(-1,0) + &
                        nf**2*z*(-2339.2921810699595_dp)*Hr2(-1,0) + nf**2*(-1748.9382716049386_dp)*Hr2(-1,0) + &
                        nf**2*z**2*(-945.9094650205761_dp)*Hr2(-1,0) + nf**2*z**3*(-165.92592592592592_dp)*Hr2(-1,0) + &
                        nf**3*z**2*(4.7407407407407405_dp)*Hr2(-1,0) + nf**3*(9.481481481481481_dp)*Hr2(-1,0) + &
                        nf**3*z*(9.481481481481481_dp)*Hr2(-1,0) + nf*z**2*(9977.68287181325_dp)*Hr2(-1,0) + &
                        nf*z**3*(13164.949245541839_dp)*Hr2(-1,0) + nf*(89503.58339337501_dp)*Hr2(-1,0) + &
                        nf*z*(91016.81613228677_dp)*Hr2(-1,0) + z*(-478456.87182410626_dp)*Hr2(0,-1) + &
                        nf*(-23136.988220940944_dp)*Hr2(0,-1) + nf*z**2*(-13025.314510261082_dp)*Hr2(0,-1) + &
                        nf*z*(15945.309403576517_dp)*Hr2(0,-1) + z**3*(53149.105862305296_dp)*Hr2(0,-1) + &
                        z**2*(474864.0359464889_dp)*Hr2(0,-1) + (1.0865153181863117e6_dp)*Hr2(0,-1) + &
                        nf*z**3*(-22871.908550525834_dp)*Hr2(0,0) + z**2*(-18086.811614252394_dp)*Hr2(0,0) + &
                        nf*z*(-17801.584943831687_dp)*Hr2(0,0) + nf**2*z**2*(-1215.2006124001132_dp)*Hr2(0,0) + &
                        nf**2*(-704.8193872885231_dp)*Hr2(0,0) + nf**3*z*(-71.63786008230454_dp)*Hr2(0,0) + &
                        nf**3*z**2*(31.60493827160494_dp)*Hr2(0,0) + nf**2*z**3*(184.36213991769546_dp)*Hr2(0,0) + &
                        nf**2*z*(197.08726226989026_dp)*Hr2(0,0) + nf*(22411.734262750135_dp)*Hr2(0,0) + &
                        nf*z**2*(25768.141882956534_dp)*Hr2(0,0) + (32934.457952627185_dp)*Hr2(0,0) + &
                        z**3*(478964.8029663421_dp)*Hr2(0,0) + z*(1.494241028286621e6_dp)*Hr2(0,0) + &
                        nf*z*(-98104.3927703713_dp)*Hr2(0,1) + nf*z**3*(-20052.952903520803_dp)*Hr2(0,1) + &
                        nf**2*(-687.846364883402_dp)*Hr2(0,1) + nf**2*z**2*(-131.58481938728852_dp)*Hr2(0,1) + &
                        nf**3*(-18.43621399176955_dp)*Hr2(0,1) + nf**3*z*(-13.87105624142661_dp)*Hr2(0,1) + &
                        nf**3*z**2*(9.305898491083676_dp)*Hr2(0,1) + nf**2*z**3*(123.90306355738456_dp)*Hr2(0,1) + &
                        nf**2*z*(2278.774577046182_dp)*Hr2(0,1) + nf*z**2*(15717.412337943231_dp)*Hr2(0,1) + &
                        nf*(50645.18817166078_dp)*Hr2(0,1) + z**2*(218827.93529159785_dp)*Hr2(0,1) + &
                        z**3*(256136.28960854342_dp)*Hr2(0,1) + (561130.1976219945_dp)*Hr2(0,1) + &
                        z*(2.0902394001692485e6_dp)*Hr2(0,1) + (-972710.2355061886_dp)*Hr2(1,0) + &
                        nf*z*(-83285.84497031062_dp)*Hr2(1,0) + z**2*(-54083.7927886745_dp)*Hr2(1,0) + &
                        nf*z**3*(-6879.283036122542_dp)*Hr2(1,0) + nf**2*(-581.8088705989941_dp)*Hr2(1,0) + &
                        nf**2*z**2*(-408.8742569730224_dp)*Hr2(1,0) + nf**3*(-10.710562414266118_dp)*Hr2(1,0) + &
                        nf**3*z**2*(-5.355281207133059_dp)*Hr2(1,0) + nf**3*z*(10.710562414266118_dp)*Hr2(1,0) + &
                        nf**2*z**3*(16.329218106995885_dp)*Hr2(1,0) + nf**2*z*(748.8321902149062_dp)*Hr2(1,0) + &
                        z**3*(1305.9983220446086_dp)*Hr2(1,0) + nf*z**2*(12889.949462750188_dp)*Hr2(1,0) + &
                        nf*(84676.11108828044_dp)*Hr2(1,0) + z*(1.002710249460488e6_dp)*Hr2(1,0) + (z*(-1.0332021039261314e6_dp) + &
                        z**3*(-458243.5408830746_dp) + nf**3*((-43.54458161865569_dp) + z**2*(-21.772290809327846_dp) + &
                        z*(43.54458161865569_dp)) + nf**2*(z*(-2826.417924096937_dp) + z**3*(-267.1787837219936_dp) + &
                        z**2*(817.1486053955191_dp) + (2255.2318244170096_dp)) + nf*((-9786.229586034264_dp) + &
                        z**2*(688.4416770331645_dp) + z**3*(10387.957933241883_dp) + z*(15932.000047854113_dp)) + &
                        z**2*(308347.0218374135_dp) + (1.034366417667588e6_dp))*Hr2(1,1) + &
                        z**2*(164160.89984793545_dp)*Hr3(-1,-1,-1) + (328321.7996958709_dp)*Hr3(-1,-1,-1) + &
                        z*(328321.7996958709_dp)*Hr3(-1,-1,-1) + nf*z*(-45258.44718792867_dp)*Hr3(-1,-1,0) + &
                        nf*(-42210.385002286246_dp)*Hr3(-1,-1,0) + nf*z**2*(-12352.1755829904_dp)*Hr3(-1,-1,0) + &
                        nf*z**3*(-3833.561957018747_dp)*Hr3(-1,-1,0) + nf**2*z**2*(190.94650205761317_dp)*Hr3(-1,-1,0) + &
                        nf**2*(381.89300411522635_dp)*Hr3(-1,-1,0) + nf**2*z*(381.89300411522635_dp)*Hr3(-1,-1,0) + &
                        z**2*(64774.55304540281_dp)*Hr3(-1,-1,0) + z**3*(146782.8148148148_dp)*Hr3(-1,-1,0) + &
                        z*(915751.9702883364_dp)*Hr3(-1,-1,0) + (920084.2172019169_dp)*Hr3(-1,-1,0) + &
                        (-164160.89984793545_dp)*Hr3(-1,0,-1) + z*(-164160.89984793545_dp)*Hr3(-1,0,-1) + &
                        z**2*(-82080.44992396772_dp)*Hr3(-1,0,-1) + (-907898.6147284494_dp)*Hr3(-1,0,0) + &
                        z*(-906645.93983133_dp)*Hr3(-1,0,0) + z**3*(-128448.79012345678_dp)*Hr3(-1,0,0) + &
                        z**2*(-111897.91641772263_dp)*Hr3(-1,0,0) + nf**2*(-300.7736625514404_dp)*Hr3(-1,0,0) + &
                        nf**2*z*(-300.7736625514404_dp)*Hr3(-1,0,0) + nf**2*z**2*(-150.3868312757202_dp)*Hr3(-1,0,0) + &
                        nf*z**3*(3302.6575217192503_dp)*Hr3(-1,0,0) + nf*z**2*(10582.518518518518_dp)*Hr3(-1,0,0) + &
                        nf*(35418.9373571102_dp)*Hr3(-1,0,0) + nf*z*(38431.07818930041_dp)*Hr3(-1,0,0) + &
                        z*(-806581.635213884_dp)*Hr3(-1,0,1) + (-804754.7380945424_dp)*Hr3(-1,0,1) + &
                        z**2*(-113542.14270982267_dp)*Hr3(-1,0,1) + z**3*(-110114.76543209876_dp)*Hr3(-1,0,1) + &
                        nf**2*(-219.65432098765433_dp)*Hr3(-1,0,1) + nf**2*z*(-219.65432098765433_dp)*Hr3(-1,0,1) + &
                        nf**2*z**2*(-109.82716049382717_dp)*Hr3(-1,0,1) + nf*z**3*(2771.753086419753_dp)*Hr3(-1,0,1) + &
                        nf*z**2*(8812.861454046639_dp)*Hr3(-1,0,1) + nf*(28627.48971193416_dp)*Hr3(-1,0,1) + &
                        nf*z*(31603.709190672158_dp)*Hr3(-1,0,1) + z**2*(-298168.7367025839_dp)*Hr3(0,-1,-1) + &
                        (-288117.7157048216_dp)*Hr3(0,-1,-1) + z*(268015.6737092969_dp)*Hr3(0,-1,-1) + &
                        (-530489.9718860553_dp)*Hr3(0,-1,0) + z*(-234285.81743854206_dp)*Hr3(0,-1,0) + &
                        z**3*(-85271.70370370372_dp)*Hr3(0,-1,0) + z**2*(-23134.12992283961_dp)*Hr3(0,-1,0) + &
                        nf*z*(-19226.688614540468_dp)*Hr3(0,-1,0) + nf**2*(-497.7777777777778_dp)*Hr3(0,-1,0) + &
                        nf**2*z**2*(-42.13991769547325_dp)*Hr3(0,-1,0) + nf**2*z*(214.3868312757202_dp)*Hr3(0,-1,0) + &
                        nf*z**3*(1946.7471422039325_dp)*Hr3(0,-1,0) + nf*z**2*(18142.990397805213_dp)*Hr3(0,-1,0) + &
                        nf*(37608.0329218107_dp)*Hr3(0,-1,0) + nf*z*(-7832.8646396354725_dp)*Hr3(0,0,-1) + &
                        nf*z**2*(-3916.4323198177362_dp)*Hr3(0,0,-1) + z*(114460.3339486555_dp)*Hr3(0,0,-1) + &
                        (128670.3981179057_dp)*Hr3(0,0,-1) + z**2*(273318.45375294395_dp)*Hr3(0,0,-1) + &
                        nf*z*(-50885.55949136883_dp)*Hr3(0,0,0) + z**2*(-27126.87687874108_dp)*Hr3(0,0,0) + &
                        nf*z**3*(-2637.1394604481025_dp)*Hr3(0,0,0) + nf**2*z**2*(-1010.304526748971_dp)*Hr3(0,0,0) + &
                        nf**3*z*(-25.28395061728395_dp)*Hr3(0,0,0) + nf**3*z**2*(12.641975308641975_dp)*Hr3(0,0,0) + &
                        (742.9640450091567_dp)*Hr3(0,0,0) + nf**2*z*(1361.119341563786_dp)*Hr3(0,0,0) + &
                        nf*(2318.222222222222_dp)*Hr3(0,0,0) + nf*z**2*(20792.277011462676_dp)*Hr3(0,0,0) + &
                        z**3*(136018.34842249658_dp)*Hr3(0,0,0) + z*(805022.8373155919_dp)*Hr3(0,0,0) + &
                        (-172627.27409288386_dp)*Hr3(0,0,1) + nf*z**3*(-5192.808413351623_dp)*Hr3(0,0,1) + &
                        nf**2*z*(-282.5130315500686_dp)*Hr3(0,0,1) + nf**2*(-207.01234567901238_dp)*Hr3(0,0,1) + &
                        nf**2*z**2*(-60.22496570644719_dp)*Hr3(0,0,1) + nf**3*z**2*(-12.641975308641975_dp)*Hr3(0,0,1) + &
                        nf**3*z*(25.28395061728395_dp)*Hr3(0,0,1) + nf*z*(22554.93212509308_dp)*Hr3(0,0,1) + &
                        nf*z**2*(24684.38853278954_dp)*Hr3(0,0,1) + nf*(27456.526748971195_dp)*Hr3(0,0,1) + &
                        z**2*(124011.06009130736_dp)*Hr3(0,0,1) + z**3*(264698.2057613169_dp)*Hr3(0,0,1) + &
                        z*(695839.2089122888_dp)*Hr3(0,0,1) + (-834356.8062526842_dp)*Hr3(0,1,0) + &
                        z**2*(-93447.95858748356_dp)*Hr3(0,1,0) + nf*z*(-10347.456790123457_dp)*Hr3(0,1,0) + &
                        nf*z**3*(-4453.955189757658_dp)*Hr3(0,1,0) + nf**2*(-268.641975308642_dp)*Hr3(0,1,0) + &
                        nf**2*z**2*(-57.503429355281206_dp)*Hr3(0,1,0) + nf**2*z*(233.52537722908096_dp)*Hr3(0,1,0) + &
                        nf*z**2*(24280.3182441701_dp)*Hr3(0,1,0) + nf*(46499.467764060355_dp)*Hr3(0,1,0) + &
                        z**3*(206967.0452674897_dp)*Hr3(0,1,0) + z*(988253.2409185411_dp)*Hr3(0,1,0) + &
                        (-349342.3544203041_dp)*Hr3(0,1,1) + nf*z*(-50729.26200274349_dp)*Hr3(0,1,1) + &
                        nf*z**3*(-5013.947873799726_dp)*Hr3(0,1,1) + nf**2*z**2*(116.76268861454048_dp)*Hr3(0,1,1) + &
                        nf**2*(362.22770919067216_dp)*Hr3(0,1,1) + nf**2*z*(858.4252400548698_dp)*Hr3(0,1,1) + &
                        nf*z**2*(14846.770919067216_dp)*Hr3(0,1,1) + nf*(23966.37585733882_dp)*Hr3(0,1,1) + &
                        z**2*(158368.76001573587_dp)*Hr3(0,1,1) + z**3*(247290.7325102881_dp)*Hr3(0,1,1) + &
                        z*(1.6781257992193543e6_dp)*Hr3(0,1,1) + (-931868.5397542159_dp)*Hr3(1,0,0) + &
                        z**2*(-121845.58949302016_dp)*Hr3(1,0,0) + nf*z*(-35570.245999085506_dp)*Hr3(1,0,0) + &
                        nf*z**3*(-2735.875628715135_dp)*Hr3(1,0,0) + nf**2*(-176.46090534979427_dp)*Hr3(1,0,0) + &
                        nf**2*z**2*(-88.23045267489714_dp)*Hr3(1,0,0) + nf**2*z*(176.46090534979427_dp)*Hr3(1,0,0) + &
                        nf*z**2*(10473.540009144948_dp)*Hr3(1,0,0) + nf*(32123.888431641513_dp)*Hr3(1,0,0) + &
                        z**3*(111665.1632373114_dp)*Hr3(1,0,0) + z*(941627.9032658758_dp)*Hr3(1,0,0) + &
                        (-984070.4121229998_dp)*Hr3(1,0,1) + z**2*(-68808.73692569742_dp)*Hr3(1,0,1) + &
                        nf*z*(-37706.798353909464_dp)*Hr3(1,0,1) + nf*z**3*(-2499.3653406492913_dp)*Hr3(1,0,1) + &
                        nf**2*z*(-89.54732510288065_dp)*Hr3(1,0,1) + nf**2*z**2*(44.773662551440324_dp)*Hr3(1,0,1) + &
                        nf**2*(89.54732510288065_dp)*Hr3(1,0,1) + nf*z**2*(12963.116598079563_dp)*Hr3(1,0,1) + &
                        nf*(31707.946959304983_dp)*Hr3(1,0,1) + z**3*(107713.31687242798_dp)*Hr3(1,0,1) + &
                        z*(1.0069576549213536e6_dp)*Hr3(1,0,1) + (-1.0505365431624555e6_dp)*Hr3(1,1,0) + &
                        z**2*(-102041.80244542517_dp)*Hr3(1,1,0) + nf*z*(-37706.798353909464_dp)*Hr3(1,1,0) + &
                        nf*z**3*(-2499.3653406492913_dp)*Hr3(1,1,0) + nf**2*z*(-89.54732510288065_dp)*Hr3(1,1,0) + &
                        nf**2*z**2*(44.773662551440324_dp)*Hr3(1,1,0) + nf**2*(89.54732510288065_dp)*Hr3(1,1,0) + &
                        nf*z**2*(12963.116598079563_dp)*Hr3(1,1,0) + nf*(31707.946959304983_dp)*Hr3(1,1,0) + &
                        z**3*(107713.31687242798_dp)*Hr3(1,1,0) + z*(1.073423785960809e6_dp)*Hr3(1,1,0) + &
                        z**3*(-66857.61316872427_dp)*Hr3(1,1,1) + nf*(-21622.167352537723_dp)*Hr3(1,1,1) + &
                        nf**2*z*(-1041.3827160493827_dp)*Hr3(1,1,1) + nf**2*z**2*(520.6913580246913_dp)*Hr3(1,1,1) + &
                        nf**2*(1041.3827160493827_dp)*Hr3(1,1,1) + nf*z**2*(2159.758573388203_dp)*Hr3(1,1,1) + &
                        nf*z**3*(2302.2441700960226_dp)*Hr3(1,1,1) + nf*z*(14072.449931412893_dp)*Hr3(1,1,1) + &
                        z*(18169.254696507065_dp)*Hr3(1,1,1) + (30453.099212957954_dp)*Hr3(1,1,1) + &
                        z**2*(92991.4220344625_dp)*Hr3(1,1,1) + (-538730.4032921811_dp)*Hr4(-1,-1,-1,0) + &
                        z*(-468487.3744855967_dp)*Hr4(-1,-1,-1,0) + z**2*(-59885.56378600822_dp)*Hr4(-1,-1,-1,0) + &
                        z**3*(-45831.374485596716_dp)*Hr4(-1,-1,-1,0) + nf*z**2*(5587.401920438957_dp)*Hr4(-1,-1,-1,0) + &
                        nf*(11174.803840877914_dp)*Hr4(-1,-1,-1,0) + nf*z*(11174.803840877914_dp)*Hr4(-1,-1,-1,0) + &
                        nf*(-9436.707818930041_dp)*Hr4(-1,-1,0,0) + nf*z*(-9436.707818930041_dp)*Hr4(-1,-1,0,0) + &
                        nf*z**2*(-4718.3539094650205_dp)*Hr4(-1,-1,0,0) + z**3*(40650.27160493827_dp)*Hr4(-1,-1,0,0) + &
                        z**2*(49091.81893004116_dp)*Hr4(-1,-1,0,0) + z*(407361.0534979424_dp)*Hr4(-1,-1,0,0) + &
                        (469898.27160493826_dp)*Hr4(-1,-1,0,0) + nf*(-7698.611796982168_dp)*Hr4(-1,-1,0,1) + &
                        nf*z*(-7698.611796982168_dp)*Hr4(-1,-1,0,1) + nf*z**2*(-3849.305898491084_dp)*Hr4(-1,-1,0,1) + &
                        z**3*(35469.168724279836_dp)*Hr4(-1,-1,0,1) + z**2*(38298.07407407407_dp)*Hr4(-1,-1,0,1) + &
                        z*(346234.7325102881_dp)*Hr4(-1,-1,0,1) + (401066.13991769555_dp)*Hr4(-1,-1,0,1) + &
                        nf*(-5587.401920438957_dp)*Hr4(-1,0,-1,0) + nf*z*(-5587.401920438957_dp)*Hr4(-1,0,-1,0) + &
                        nf*z**2*(-2793.7009602194785_dp)*Hr4(-1,0,-1,0) + z**3*(22915.687242798358_dp)*Hr4(-1,0,-1,0) + &
                        z**2*(29942.78189300411_dp)*Hr4(-1,0,-1,0) + z*(234243.68724279836_dp)*Hr4(-1,0,-1,0) + &
                        (269365.20164609054_dp)*Hr4(-1,0,-1,0) + (-200533.06995884777_dp)*Hr4(-1,0,0,0) + &
                        z*(-173117.36625514404_dp)*Hr4(-1,0,0,0) + z**2*(-19149.037037037036_dp)*Hr4(-1,0,0,0) + &
                        z**3*(-17734.584362139918_dp)*Hr4(-1,0,0,0) + nf*z**2*(1924.652949245542_dp)*Hr4(-1,0,0,0) + &
                        nf*(3849.305898491084_dp)*Hr4(-1,0,0,0) + nf*z*(3849.305898491084_dp)*Hr4(-1,0,0,0) + &
                        (-332234.00823045266_dp)*Hr4(-1,0,0,1) + z*(-285108.41152263375_dp)*Hr4(-1,0,0,1) + &
                        z**3*(-30288.065843621403_dp)*Hr4(-1,0,0,1) + z**2*(-27504.329218106996_dp)*Hr4(-1,0,0,1) + &
                        nf*z**2*(2980.2578875171466_dp)*Hr4(-1,0,0,1) + nf*(5960.515775034293_dp)*Hr4(-1,0,0,1) + &
                        nf*z*(5960.515775034293_dp)*Hr4(-1,0,0,1) + (-131700.93827160494_dp)*Hr4(-1,0,1,0) + &
                        z*(-111991.04526748972_dp)*Hr4(-1,0,1,0) + z**3*(-12553.481481481482_dp)*Hr4(-1,0,1,0) + &
                        z**2*(-8355.292181069959_dp)*Hr4(-1,0,1,0) + nf*z**2*(1055.6049382716049_dp)*Hr4(-1,0,1,0) + &
                        nf*(2111.2098765432097_dp)*Hr4(-1,0,1,0) + nf*z*(2111.2098765432097_dp)*Hr4(-1,0,1,0) + &
                        (-263401.8765432099_dp)*Hr4(-1,0,1,1) + z*(-223982.09053497945_dp)*Hr4(-1,0,1,1) + &
                        z**3*(-25106.962962962964_dp)*Hr4(-1,0,1,1) + z**2*(-16710.584362139918_dp)*Hr4(-1,0,1,1) + &
                        nf*z**2*(2111.2098765432097_dp)*Hr4(-1,0,1,1) + nf*(4222.419753086419_dp)*Hr4(-1,0,1,1) + &
                        nf*z*(4222.419753086419_dp)*Hr4(-1,0,1,1) + z*(-135333.39917695473_dp)*Hr4(0,-1,-1,0) + &
                        nf*(-11500.685871056241_dp)*Hr4(0,-1,-1,0) + nf*z**2*(-6910.244170096022_dp)*Hr4(0,-1,-1,0) + &
                        nf*z*(8056.449931412894_dp)*Hr4(0,-1,-1,0) + z**3*(25427.22633744856_dp)*Hr4(0,-1,-1,0) + &
                        z**2*(126193.2510288066_dp)*Hr4(0,-1,-1,0) + (419944.2962962963_dp)*Hr4(0,-1,-1,0) + &
                        (-367208.2962962963_dp)*Hr4(0,-1,0,0) + z**2*(-105536.2633744856_dp)*Hr4(0,-1,0,0) + &
                        z**3*(-22512.197530864192_dp)*Hr4(0,-1,0,0) + nf*z*(-6860.905349794238_dp)*Hr4(0,-1,0,0) + &
                        nf*z**2*(5686.781893004115_dp)*Hr4(0,-1,0,0) + nf*(9907.97256515775_dp)*Hr4(0,-1,0,0) + &
                        z*(115857.11934156377_dp)*Hr4(0,-1,0,0) + (-314472.2962962963_dp)*Hr4(0,-1,0,1) + &
                        z**2*(-84879.2757201646_dp)*Hr4(0,-1,0,1) + z**3*(-19597.168724279836_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(-5665.360768175583_dp)*Hr4(0,-1,0,1) + nf*z**2*(4463.319615912208_dp)*Hr4(0,-1,0,1) + &
                        nf*(8315.25925925926_dp)*Hr4(0,-1,0,1) + z*(96380.83950617284_dp)*Hr4(0,-1,0,1) + &
                        (-178389.33333333334_dp)*Hr4(0,0,-1,0) + z**2*(-126674.69958847738_dp)*Hr4(0,0,-1,0) + &
                        z*(-110851.68724279835_dp)*Hr4(0,0,-1,0) + z**3*(-8655.539094650208_dp)*Hr4(0,0,-1,0) + &
                        nf**2*z*(-341.3333333333333_dp)*Hr4(0,0,-1,0) + nf**2*z**2*(-170.66666666666666_dp)*Hr4(0,0,-1,0) + &
                        nf*(5660.444444444445_dp)*Hr4(0,0,-1,0) + nf*z**2*(7381.157750342935_dp)*Hr4(0,0,-1,0) + &
                        nf*z*(9967.846364883402_dp)*Hr4(0,0,-1,0) + z**2*(-128658.78737997256_dp)*Hr4(0,0,0,0) + &
                        z*(-39420.48834019205_dp)*Hr4(0,0,0,0) + nf**2*z**2*(-488.82304526748965_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(-387.6872427983539_dp)*Hr4(0,0,0,0) + z**3*(12163.687242798353_dp)*Hr4(0,0,0,0) + &
                        nf*z*(12663.045267489711_dp)*Hr4(0,0,0,0) + nf*z**2*(14043.82990397805_dp)*Hr4(0,0,0,0) + &
                        z**2*(-152410.95198902607_dp)*Hr4(0,0,0,1) + (-144384._dp)*Hr4(0,0,0,1) + &
                        nf*z*(-23950.924554183814_dp)*Hr4(0,0,0,1) + nf**2*z**2*(-206.4855967078189_dp)*Hr4(0,0,0,1) + &
                        nf**2*z*(868.0823045267491_dp)*Hr4(0,0,0,1) + nf*(4693.333333333333_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(11119.319615912209_dp)*Hr4(0,0,0,1) + z**3*(37622.51851851852_dp)*Hr4(0,0,0,1) + &
                        z*(335636.7187928669_dp)*Hr4(0,0,0,1) + (-406784._dp)*Hr4(0,0,1,0) + &
                        z**2*(-140862.9465020576_dp)*Hr4(0,0,1,0) + nf**2*z**2*(-37.925925925925924_dp)*Hr4(0,0,1,0) + &
                        nf**2*z*(75.85185185185185_dp)*Hr4(0,0,1,0) + nf*z*(522.1838134430727_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(8652.554183813443_dp)*Hr4(0,0,1,0) + nf*(11415.703703703704_dp)*Hr4(0,0,1,0) + &
                        z*(29590.123456790123_dp)*Hr4(0,0,1,0) + z**3*(49055.07818930041_dp)*Hr4(0,0,1,0) + &
                        (-478080._dp)*Hr4(0,0,1,1) + z*(-151457.1851851852_dp)*Hr4(0,0,1,1) + &
                        z**2*(-32512.877914951994_dp)*Hr4(0,0,1,1) + nf**2*z*(-375.045267489712_dp)*Hr4(0,0,1,1) + &
                        nf**2*z**2*(187.522633744856_dp)*Hr4(0,0,1,1) + nf*z**2*(2672.373113854596_dp)*Hr4(0,0,1,1) + &
                        nf*(12610.370370370372_dp)*Hr4(0,0,1,1) + nf*z*(17522.48010973937_dp)*Hr4(0,0,1,1) + &
                        z**3*(80767.4732510288_dp)*Hr4(0,0,1,1) + (-444875.8518518519_dp)*Hr4(0,1,0,0) + &
                        z**2*(-88151.61591220851_dp)*Hr4(0,1,0,0) + nf*z*(2943.2976680384086_dp)*Hr4(0,1,0,0) + &
                        nf*z**2*(5758.770919067217_dp)*Hr4(0,1,0,0) + nf*(10318.836762688614_dp)*Hr4(0,1,0,0) + &
                        z*(29016.4060356653_dp)*Hr4(0,1,0,0) + z**3*(39401.87654320988_dp)*Hr4(0,1,0,0) + &
                        (-684610.3703703703_dp)*Hr4(0,1,0,1) + z**2*(-65639.76954732511_dp)*Hr4(0,1,0,1) + &
                        nf*z*(-2600.735253772291_dp)*Hr4(0,1,0,1) + nf*z**2*(5939.972565157751_dp)*Hr4(0,1,0,1) + &
                        nf*(13831.725651577504_dp)*Hr4(0,1,0,1) + z**3*(76507.12757201646_dp)*Hr4(0,1,0,1) + &
                        z*(266223.1440329218_dp)*Hr4(0,1,0,1) + (-684610.3703703703_dp)*Hr4(0,1,1,0) + &
                        z**2*(-65639.76954732511_dp)*Hr4(0,1,1,0) + nf*z*(-2600.735253772291_dp)*Hr4(0,1,1,0) + &
                        nf*z**2*(5939.972565157751_dp)*Hr4(0,1,1,0) + nf*(13831.725651577504_dp)*Hr4(0,1,1,0) + &
                        z**3*(76507.12757201646_dp)*Hr4(0,1,1,0) + z*(266223.1440329218_dp)*Hr4(0,1,1,0) + &
                        (-560995.5555555555_dp)*Hr4(0,1,1,1) + nf*z*(-20218.732510288064_dp)*Hr4(0,1,1,1) + &
                        nf*z**2*(-3126.781893004115_dp)*Hr4(0,1,1,1) + nf*(5760.526748971193_dp)*Hr4(0,1,1,1) + &
                        z**3*(102570.66666666667_dp)*Hr4(0,1,1,1) + z**2*(129506.50205761316_dp)*Hr4(0,1,1,1) + &
                        z*(770049.0534979425_dp)*Hr4(0,1,1,1) + (-216047.93415637861_dp)*Hr4(1,0,0,0) + &
                        z**2*(-20183.045267489713_dp)*Hr4(1,0,0,0) + nf*z*(-3927.264746227709_dp)*Hr4(1,0,0,0) + &
                        nf*z**2*(1963.6323731138546_dp)*Hr4(1,0,0,0) + nf*(3927.264746227709_dp)*Hr4(1,0,0,0) + &
                        z**3*(18830.222222222223_dp)*Hr4(1,0,0,0) + z*(189306.46913580247_dp)*Hr4(1,0,0,0) + &
                        (-528707.1604938272_dp)*Hr4(1,0,0,1) + z**2*(-44666.600823045264_dp)*Hr4(1,0,0,1) + &
                        nf*z*(-8994.238683127573_dp)*Hr4(1,0,0,1) + nf*z**2*(4497.119341563786_dp)*Hr4(1,0,0,1) + &
                        nf*(8994.238683127573_dp)*Hr4(1,0,0,1) + z**3*(46328.62551440329_dp)*Hr4(1,0,0,1) + &
                        z*(462827.720164609_dp)*Hr4(1,0,0,1) + (-625318.4526748972_dp)*Hr4(1,0,1,0) + &
                        z**2*(-48967.11111111111_dp)*Hr4(1,0,1,0) + nf*z*(-10133.947873799725_dp)*Hr4(1,0,1,0) + &
                        nf*z**2*(5066.973936899863_dp)*Hr4(1,0,1,0) + nf*(10133.947873799725_dp)*Hr4(1,0,1,0) + &
                        z**3*(54996.80658436214_dp)*Hr4(1,0,1,0) + z*(547042.5020576132_dp)*Hr4(1,0,1,0) + &
                        (-708655.9341563786_dp)*Hr4(1,0,1,1) + z**2*(-29487.40740740741_dp)*Hr4(1,0,1,1) + &
                        nf*z*(-8338.787379972566_dp)*Hr4(1,0,1,1) + nf*z**2*(4169.393689986283_dp)*Hr4(1,0,1,1) + &
                        nf*(8338.787379972566_dp)*Hr4(1,0,1,1) + z**3*(65458.04115226338_dp)*Hr4(1,0,1,1) + &
                        z*(616307.8847736625_dp)*Hr4(1,0,1,1) + (-528707.1604938272_dp)*Hr4(1,1,0,0) + &
                        z**2*(-44666.600823045264_dp)*Hr4(1,1,0,0) + nf*z*(-8994.238683127573_dp)*Hr4(1,1,0,0) + &
                        nf*z**2*(4497.119341563786_dp)*Hr4(1,1,0,0) + nf*(8994.238683127573_dp)*Hr4(1,1,0,0) + &
                        z**3*(46328.62551440329_dp)*Hr4(1,1,0,0) + z*(462827.720164609_dp)*Hr4(1,1,0,0) + &
                        (-708655.9341563786_dp)*Hr4(1,1,0,1) + z**2*(-29487.40740740741_dp)*Hr4(1,1,0,1) + &
                        nf*z*(-8338.787379972566_dp)*Hr4(1,1,0,1) + nf*z**2*(4169.393689986283_dp)*Hr4(1,1,0,1) + &
                        nf*(8338.787379972566_dp)*Hr4(1,1,0,1) + z**3*(65458.04115226338_dp)*Hr4(1,1,0,1) + &
                        z*(616307.8847736625_dp)*Hr4(1,1,0,1) + (-708655.9341563786_dp)*Hr4(1,1,1,0) + &
                        z**2*(-29487.40740740741_dp)*Hr4(1,1,1,0) + nf*z*(-8338.787379972566_dp)*Hr4(1,1,1,0) + &
                        nf*z**2*(4169.393689986283_dp)*Hr4(1,1,1,0) + nf*(8338.787379972566_dp)*Hr4(1,1,1,0) + &
                        z**3*(65458.04115226338_dp)*Hr4(1,1,1,0) + z*(616307.8847736625_dp)*Hr4(1,1,1,0) + &
                        (-333349.9259259259_dp)*Hr4(1,1,1,1) + nf*(-7180.641975308642_dp)*Hr4(1,1,1,1) + &
                        nf*z**2*(-3590.320987654321_dp)*Hr4(1,1,1,1) + nf*z*(7180.641975308642_dp)*Hr4(1,1,1,1) + &
                        z**3*(41844.93827160494_dp)*Hr4(1,1,1,1) + z**2*(77918.81481481482_dp)*Hr4(1,1,1,1) + &
                        z*(277061.5308641975_dp)*Hr4(1,1,1,1) + z**2*(78708.93827160494_dp)*Hr5(-1,-1,-1,-1,0) + &
                        (157417.8765432099_dp)*Hr5(-1,-1,-1,-1,0) + z*(157417.8765432099_dp)*Hr5(-1,-1,-1,-1,0) + &
                        (-139152.32921810698_dp)*Hr5(-1,-1,-1,0,0) + z*(-139152.32921810698_dp)*Hr5(-1,-1,-1,0,0) + &
                        z**2*(-69576.16460905349_dp)*Hr5(-1,-1,-1,0,0) + (-120886.78189300411_dp)*Hr5(-1,-1,-1,0,1) + &
                        z*(-120886.78189300411_dp)*Hr5(-1,-1,-1,0,1) + z**2*(-60443.390946502055_dp)*Hr5(-1,-1,-1,0,1) + &
                        (-78708.93827160494_dp)*Hr5(-1,-1,0,-1,0) + z*(-78708.93827160494_dp)*Hr5(-1,-1,0,-1,0) + &
                        z**2*(-39354.46913580247_dp)*Hr5(-1,-1,0,-1,0) + z**2*(30221.695473251028_dp)*Hr5(-1,-1,0,0,0) + &
                        (60443.390946502055_dp)*Hr5(-1,-1,0,0,0) + z*(60443.390946502055_dp)*Hr5(-1,-1,0,0,0) + &
                        z**2*(51310.61728395062_dp)*Hr5(-1,-1,0,0,1) + (102621.23456790124_dp)*Hr5(-1,-1,0,0,1) + &
                        z*(102621.23456790124_dp)*Hr5(-1,-1,0,0,1) + z**2*(21088.92181069959_dp)*Hr5(-1,-1,0,1,0) + &
                        (42177.84362139918_dp)*Hr5(-1,-1,0,1,0) + z*(42177.84362139918_dp)*Hr5(-1,-1,0,1,0) + &
                        z**2*(42177.84362139918_dp)*Hr5(-1,-1,0,1,1) + (84355.68724279836_dp)*Hr5(-1,-1,0,1,1) + &
                        z*(84355.68724279836_dp)*Hr5(-1,-1,0,1,1) + (-78708.93827160494_dp)*Hr5(-1,0,-1,-1,0) + &
                        z*(-78708.93827160494_dp)*Hr5(-1,0,-1,-1,0) + z**2*(-39354.46913580247_dp)*Hr5(-1,0,-1,-1,0) + &
                        z**2*(34788.082304526746_dp)*Hr5(-1,0,-1,0,0) + (69576.16460905349_dp)*Hr5(-1,0,-1,0,0) + &
                        z*(69576.16460905349_dp)*Hr5(-1,0,-1,0,0) + z**2*(30221.695473251028_dp)*Hr5(-1,0,-1,0,1) + &
                        (60443.390946502055_dp)*Hr5(-1,0,-1,0,1) + z*(60443.390946502055_dp)*Hr5(-1,0,-1,0,1) + &
                        z**2*(13118.156378600823_dp)*Hr5(-1,0,0,-1,0) + (26236.312757201646_dp)*Hr5(-1,0,0,-1,0) + &
                        z*(26236.312757201646_dp)*Hr5(-1,0,0,-1,0) + (-17103.539094650205_dp)*Hr5(-1,0,0,0,0) + &
                        z*(-17103.539094650205_dp)*Hr5(-1,0,0,0,0) + z**2*(-8551.769547325102_dp)*Hr5(-1,0,0,0,0) + &
                        (-42177.84362139918_dp)*Hr5(-1,0,0,0,1) + z*(-42177.84362139918_dp)*Hr5(-1,0,0,0,1) + &
                        z**2*(-21088.92181069959_dp)*Hr5(-1,0,0,0,1) + (-33045.06995884774_dp)*Hr5(-1,0,0,1,0) + &
                        z*(-33045.06995884774_dp)*Hr5(-1,0,0,1,0) + z**2*(-16522.53497942387_dp)*Hr5(-1,0,0,1,0) + &
                        (-66090.13991769547_dp)*Hr5(-1,0,0,1,1) + z*(-66090.13991769547_dp)*Hr5(-1,0,0,1,1) + &
                        z**2*(-33045.06995884774_dp)*Hr5(-1,0,0,1,1) + (-7970.765432098766_dp)*Hr5(-1,0,1,0,0) + &
                        z*(-7970.765432098766_dp)*Hr5(-1,0,1,0,0) + z**2*(-3985.382716049383_dp)*Hr5(-1,0,1,0,0) + &
                        (-23912.296296296296_dp)*Hr5(-1,0,1,0,1) + z*(-23912.296296296296_dp)*Hr5(-1,0,1,0,1) + &
                        z**2*(-11956.148148148148_dp)*Hr5(-1,0,1,0,1) + (-23912.296296296296_dp)*Hr5(-1,0,1,1,0) + &
                        z*(-23912.296296296296_dp)*Hr5(-1,0,1,1,0) + z**2*(-11956.148148148148_dp)*Hr5(-1,0,1,1,0) + &
                        (-47824.59259259259_dp)*Hr5(-1,0,1,1,1) + z*(-47824.59259259259_dp)*Hr5(-1,0,1,1,1) + &
                        z**2*(-23912.296296296296_dp)*Hr5(-1,0,1,1,1) + z**2*(-142475.06172839506_dp)*Hr5(0,-1,-1,-1,0) + &
                        (-137494.12345679014_dp)*Hr5(0,-1,-1,-1,0) + z*(127532.24691358025_dp)*Hr5(0,-1,-1,-1,0) + &
                        z*(-113350.05761316874_dp)*Hr5(0,-1,-1,0,0) + (121950.81481481482_dp)*Hr5(0,-1,-1,0,0) + &
                        z**2*(126251.19341563786_dp)*Hr5(0,-1,-1,0,0) + z*(-99167.8683127572_dp)*Hr5(0,-1,-1,0,1) + &
                        (106407.50617283951_dp)*Hr5(0,-1,-1,0,1) + z**2*(110027.32510288066_dp)*Hr5(0,-1,-1,0,1) + &
                        z*(-63766.12345679013_dp)*Hr5(0,-1,0,-1,0) + (68747.06172839507_dp)*Hr5(0,-1,0,-1,0) + &
                        z**2*(71237.53086419753_dp)*Hr5(0,-1,0,-1,0) + z**2*(-55013.66255144033_dp)*Hr5(0,-1,0,0,0) + &
                        (-53203.753086419754_dp)*Hr5(0,-1,0,0,0) + z*(49583.9341563786_dp)*Hr5(0,-1,0,0,0) + &
                        z**2*(-93803.45679012345_dp)*Hr5(0,-1,0,0,1) + (-90864.1975308642_dp)*Hr5(0,-1,0,0,1) + &
                        z*(84985.67901234567_dp)*Hr5(0,-1,0,0,1) + z**2*(-38789.79423868313_dp)*Hr5(0,-1,0,1,0) + &
                        (-37660.444444444445_dp)*Hr5(0,-1,0,1,0) + z*(35401.74485596709_dp)*Hr5(0,-1,0,1,0) + &
                        z**2*(-77579.58847736625_dp)*Hr5(0,-1,0,1,1) + (-75320.88888888889_dp)*Hr5(0,-1,0,1,1) + &
                        z*(70803.48971193418_dp)*Hr5(0,-1,0,1,1) + nf*z*(-3607.1769547325107_dp)*Hr5(0,0,-1,-1,0) + &
                        nf*z**2*(-1803.5884773662553_dp)*Hr5(0,0,-1,-1,0) + z*(54284.64197530865_dp)*Hr5(0,0,-1,-1,0) + &
                        (61212.444444444445_dp)*Hr5(0,0,-1,-1,0) + z**2*(130262.91358024691_dp)*Hr5(0,0,-1,-1,0) + &
                        z**2*(-115644.57613168725_dp)*Hr5(0,0,-1,0,0) + (-54414.22222222222_dp)*Hr5(0,0,-1,0,0) + &
                        z*(-48362.93004115226_dp)*Hr5(0,0,-1,0,0) + nf*z**2*(1641.3497942386834_dp)*Hr5(0,0,-1,0,0) + &
                        nf*z*(3282.699588477367_dp)*Hr5(0,0,-1,0,0) + z**2*(-101026.23868312757_dp)*Hr5(0,0,-1,0,1) + &
                        (-47616._dp)*Hr5(0,0,-1,0,1) + z*(-42441.21810699588_dp)*Hr5(0,0,-1,0,1) + &
                        nf*z**2*(1479.111111111111_dp)*Hr5(0,0,-1,0,1) + nf*z*(2958.222222222222_dp)*Hr5(0,0,-1,0,1) + &
                        z**2*(-80065.84362139918_dp)*Hr5(0,0,0,-1,0) + (-18432._dp)*Hr5(0,0,0,-1,0) + &
                        nf*z*(-2292.4115226337453_dp)*Hr5(0,0,0,-1,0) + nf*z**2*(2949.7942386831273_dp)*Hr5(0,0,0,-1,0) + &
                        z*(15844.609053497941_dp)*Hr5(0,0,0,-1,0) + z**2*(-85819.34705075447_dp)*Hr5(0,0,0,0,0) + &
                        nf*z*(-3084.6419753086425_dp)*Hr5(0,0,0,0,0) + nf**2*z**2*(-50.5679012345679_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z*(101.1358024691358_dp)*Hr5(0,0,0,0,0) + nf*z**2*(6556.9711934156385_dp)*Hr5(0,0,0,0,0) + &
                        z*(15088.899862825787_dp)*Hr5(0,0,0,0,0) + z**2*(-191846.18930041153_dp)*Hr5(0,0,0,0,1) + &
                        z*(-55108.47736625515_dp)*Hr5(0,0,0,0,1) + (-18432._dp)*Hr5(0,0,0,0,1) + &
                        nf*z*(6978.370370370371_dp)*Hr5(0,0,0,0,1) + nf*z**2*(8002.37037037037_dp)*Hr5(0,0,0,0,1) + &
                        z**2*(-228284.5761316872_dp)*Hr5(0,0,0,1,0) + (-66048._dp)*Hr5(0,0,0,1,0) + &
                        nf*z*(-2739.094650205762_dp)*Hr5(0,0,0,1,0) + nf*z**2*(4896.658436213992_dp)*Hr5(0,0,0,1,0) + &
                        z*(58007.7037037037_dp)*Hr5(0,0,0,1,0) + z**2*(-305090.19478738005_dp)*Hr5(0,0,0,1,1) + &
                        (-97792._dp)*Hr5(0,0,0,1,1) + nf*z*(-9009.514403292183_dp)*Hr5(0,0,0,1,1) + &
                        nf*z**2*(5073.6460905349795_dp)*Hr5(0,0,0,1,1) + z*(144814.52949245545_dp)*Hr5(0,0,0,1,1) + &
                        z**2*(-164815.18792866942_dp)*Hr5(0,0,1,0,0) + (-93212.44444444444_dp)*Hr5(0,0,1,0,0) + &
                        nf*z*(-2903.4403292181073_dp)*Hr5(0,0,1,0,0) + nf*z**2*(1451.7201646090537_dp)*Hr5(0,0,1,0,0) + &
                        z*(7920.548696844993_dp)*Hr5(0,0,1,0,0) + z**2*(-287593.3497942387_dp)*Hr5(0,0,1,0,1) + &
                        (-184149.33333333334_dp)*Hr5(0,0,1,0,1) + z*(-59067.52263374485_dp)*Hr5(0,0,1,0,1) + &
                        nf*z*(-3775.736625514404_dp)*Hr5(0,0,1,0,1) + nf*z**2*(1887.868312757202_dp)*Hr5(0,0,1,0,1) + &
                        z**2*(-287593.3497942387_dp)*Hr5(0,0,1,1,0) + (-184149.33333333334_dp)*Hr5(0,0,1,1,0) + &
                        z*(-59067.52263374485_dp)*Hr5(0,0,1,1,0) + nf*z*(-3775.736625514404_dp)*Hr5(0,0,1,1,0) + &
                        nf*z**2*(1887.868312757202_dp)*Hr5(0,0,1,1,0) + z**2*(-324479.4732510288_dp)*Hr5(0,0,1,1,1) + &
                        (-252416.00000000003_dp)*Hr5(0,0,1,1,1) + z*(-218729.3497942387_dp)*Hr5(0,0,1,1,1) + &
                        nf*z*(-1643.4567901234568_dp)*Hr5(0,0,1,1,1) + nf*z**2*(821.7283950617284_dp)*Hr5(0,0,1,1,1) + &
                        (-64461.43209876544_dp)*Hr5(0,1,0,0,0) + z**2*(-63581.76131687243_dp)*Hr5(0,1,0,0,0) + &
                        z*(-42308.477366255145_dp)*Hr5(0,1,0,0,0) + (-170868.93827160497_dp)*Hr5(0,1,0,0,1) + &
                        z**2*(-162525.23456790124_dp)*Hr5(0,1,0,0,1) + z*(-91907.16049382718_dp)*Hr5(0,1,0,0,1) + &
                        (-212815.01234567902_dp)*Hr5(0,1,0,1,0) + z**2*(-197886.94650205763_dp)*Hr5(0,1,0,1,0) + &
                        z*(-99197.36625514404_dp)*Hr5(0,1,0,1,0) + (-292023.3086419752_dp)*Hr5(0,1,0,1,1) + &
                        z**2*(-255300.47736625513_dp)*Hr5(0,1,0,1,1) + z*(-78521.41563786009_dp)*Hr5(0,1,0,1,1) + &
                        (-170868.93827160497_dp)*Hr5(0,1,1,0,0) + z**2*(-162525.23456790124_dp)*Hr5(0,1,1,0,0) + &
                        z*(-91907.16049382718_dp)*Hr5(0,1,1,0,0) + (-292023.3086419752_dp)*Hr5(0,1,1,0,1) + &
                        z**2*(-255300.47736625513_dp)*Hr5(0,1,1,0,1) + z*(-78521.41563786009_dp)*Hr5(0,1,1,0,1) + &
                        (-292023.3086419752_dp)*Hr5(0,1,1,1,0) + z**2*(-255300.47736625513_dp)*Hr5(0,1,1,1,0) + &
                        z*(-78521.41563786009_dp)*Hr5(0,1,1,1,0) + (-316833.18518518517_dp)*Hr5(0,1,1,1,1) + &
                        z**2*(-229654.1234567901_dp)*Hr5(0,1,1,1,1) + z*(82703.8024691358_dp)*Hr5(0,1,1,1,1) + &
                        (-19926.913580246914_dp)*Hr5(1,0,0,0,0) + z**2*(-9963.456790123457_dp)*Hr5(1,0,0,0,0) + &
                        z*(19926.913580246914_dp)*Hr5(1,0,0,0,0) + (-73341.71742112482_dp)*Hr5(1,0,0,0,1) + &
                        z**2*(-36670.85871056241_dp)*Hr5(1,0,0,0,1) + z*(73341.71742112482_dp)*Hr5(1,0,0,0,1) + &
                        (-120390.5843621399_dp)*Hr5(1,0,0,1,0) + z**2*(-60195.29218106995_dp)*Hr5(1,0,0,1,0) + &
                        z*(120390.5843621399_dp)*Hr5(1,0,0,1,0) + (-187089.9972565158_dp)*Hr5(1,0,0,1,1) + &
                        z**2*(-93544.9986282579_dp)*Hr5(1,0,0,1,1) + z*(187089.9972565158_dp)*Hr5(1,0,0,1,1) + &
                        (-120390.5843621399_dp)*Hr5(1,0,1,0,0) + z**2*(-60195.29218106995_dp)*Hr5(1,0,1,0,0) + &
                        z*(120390.5843621399_dp)*Hr5(1,0,1,0,0) + (-227496.55967078192_dp)*Hr5(1,0,1,0,1) + &
                        z**2*(-113748.27983539096_dp)*Hr5(1,0,1,0,1) + z*(227496.55967078192_dp)*Hr5(1,0,1,0,1) + &
                        (-227496.55967078192_dp)*Hr5(1,0,1,1,0) + z**2*(-113748.27983539096_dp)*Hr5(1,0,1,1,0) + &
                        z*(227496.55967078192_dp)*Hr5(1,0,1,1,0) + (-293919.6049382717_dp)*Hr5(1,0,1,1,1) + &
                        z**2*(-146959.80246913584_dp)*Hr5(1,0,1,1,1) + z*(293919.6049382717_dp)*Hr5(1,0,1,1,1) + &
                        (-73341.71742112482_dp)*Hr5(1,1,0,0,0) + z**2*(-36670.85871056241_dp)*Hr5(1,1,0,0,0) + &
                        z*(73341.71742112482_dp)*Hr5(1,1,0,0,0) + (-187089.9972565158_dp)*Hr5(1,1,0,0,1) + &
                        z**2*(-93544.9986282579_dp)*Hr5(1,1,0,0,1) + z*(187089.9972565158_dp)*Hr5(1,1,0,0,1) + &
                        (-227496.55967078192_dp)*Hr5(1,1,0,1,0) + z**2*(-113748.27983539096_dp)*Hr5(1,1,0,1,0) + &
                        z*(227496.55967078192_dp)*Hr5(1,1,0,1,0) + (-293919.6049382717_dp)*Hr5(1,1,0,1,1) + &
                        z**2*(-146959.80246913584_dp)*Hr5(1,1,0,1,1) + z*(293919.6049382717_dp)*Hr5(1,1,0,1,1) + &
                        (-187089.9972565158_dp)*Hr5(1,1,1,0,0) + z**2*(-93544.9986282579_dp)*Hr5(1,1,1,0,0) + &
                        z*(187089.9972565158_dp)*Hr5(1,1,1,0,0) + (-293919.6049382717_dp)*Hr5(1,1,1,0,1) + &
                        z**2*(-146959.80246913584_dp)*Hr5(1,1,1,0,1) + z*(293919.6049382717_dp)*Hr5(1,1,1,0,1) + &
                        (-293919.6049382717_dp)*Hr5(1,1,1,1,0) + z**2*(-146959.80246913584_dp)*Hr5(1,1,1,1,0) + &
                        z*(293919.6049382717_dp)*Hr5(1,1,1,1,0) + (-265692.1810699589_dp)*Hr5(1,1,1,1,1) + &
                        z**2*(-132846.09053497945_dp)*Hr5(1,1,1,1,1) + z*(265692.1810699589_dp)*Hr5(1,1,1,1,1))/z

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
                        Ibar_select = ((-252912.83048638407_dp) + nf*z*(-52472.41982225307_dp) + nf*z**3*(-5596.990912219252_dp) + &
                        nf*z**2*(-3152.062111814953_dp) + nf**2*(-914.6326481394249_dp) + nf**2*z**2*(-211.31741811183014_dp) + &
                        nf**3*z*(-16.961316872427986_dp) + nf**3*z**3*(-1.1588477366255143_dp) + nf**3*z**2*(8.506995884773664_dp) + &
                        nf**3*(8.981069958847739_dp) + nf**2*z**3*(152.60210333790582_dp) + nf**2*z*(1341.1613673678587_dp) + &
                        nf*(28010.198100491274_dp) + z**2*(31753.466248695815_dp) + z**3*(51592.4561917762_dp) + &
                        z*(566076.2692518136_dp) + (-167100.6787091279_dp)*Hr1(0) + nf*z**3*(-1991.4973936899862_dp)*Hr1(0) + &
                        nf**2*(-434.8283493369913_dp)*Hr1(0) + nf**2*z*(-306.9348581807918_dp)*Hr1(0) + &
                        nf**2*z**2*(-149.44834594389764_dp)*Hr1(0) + nf**3*z**2*(-1.474897119341564_dp)*Hr1(0) + &
                        nf**3*z*(0.10534979423868313_dp)*Hr1(0) + nf**3*(3.476543209876543_dp)*Hr1(0) + &
                        nf**2*z**3*(23.258893461362597_dp)*Hr1(0) + nf*z**2*(1381.2604369529136_dp)*Hr1(0) + &
                        nf*z*(4039.195960041735_dp)*Hr1(0) + nf*(15708.128185365827_dp)*Hr1(0) + z**2*(17171.684866528172_dp)*Hr1(0) &
                        + z**3*(38724.70527877619_dp)*Hr1(0) + z*(76702.3588096652_dp)*Hr1(0) + (-353306.395555647_dp)*Hr1(1) + &
                        z**2*(-48561.40240202713_dp)*Hr1(1) + nf*z*(-25922.28823836733_dp)*Hr1(1) + &
                        nf*z**3*(-3504.1097393689993_dp)*Hr1(1) + nf**2*(-488.70598994055786_dp)*Hr1(1) + &
                        nf**2*z**2*(-230.49657064471882_dp)*Hr1(1) + nf**3*z*(-2.633744855967078_dp)*Hr1(1) + &
                        nf**3*z**2*(1.316872427983539_dp)*Hr1(1) + nf**3*(2.633744855967078_dp)*Hr1(1) + &
                        nf**2*z**3*(46.23685413808872_dp)*Hr1(1) + nf**2*z*(627.2702331961591_dp)*Hr1(1) + &
                        nf*z**2*(3539.536437427834_dp)*Hr1(1) + nf*(24925.679184869383_dp)*Hr1(1) + &
                        z**3*(58709.939587556844_dp)*Hr1(1) + z*(367361.68156790704_dp)*Hr1(1) + (-65038.350062189675_dp)*Hr2(0,0) + &
                        z**2*(-15637.518176694688_dp)*Hr2(0,0) + nf*z*(-5104.337750675571_dp)*Hr2(0,0) + &
                        nf*z**3*(-191.45569272976684_dp)*Hr2(0,0) + nf**2*z**2*(-78.09931412894376_dp)*Hr2(0,0) + &
                        nf**2*(-68.93388203017832_dp)*Hr2(0,0) + nf**3*z*(-3.7925925925925927_dp)*Hr2(0,0) + &
                        nf**3*z**2*(1.8962962962962964_dp)*Hr2(0,0) + nf**2*z*(200.1294924554184_dp)*Hr2(0,0) + &
                        nf*z**2*(3109.1320047141703_dp)*Hr2(0,0) + nf*(4389.925925925926_dp)*Hr2(0,0) + &
                        z**3*(7596.879012345679_dp)*Hr2(0,0) + z*(55732.98823677771_dp)*Hr2(0,0) + (-161196.4854293372_dp)*Hr2(0,1) &
                        + z**2*(-31746.645190491698_dp)*Hr2(0,1) + nf*z**3*(-565.3772290809328_dp)*Hr2(0,1) + &
                        nf**2*(-138.71056241426612_dp)*Hr2(0,1) + nf**2*z*(-59.34705075445816_dp)*Hr2(0,1) + &
                        nf**2*z**2*(-38.27709190672154_dp)*Hr2(0,1) + nf*z*(3971.687242798354_dp)*Hr2(0,1) + &
                        nf*z**2*(4573.761316872428_dp)*Hr2(0,1) + nf*(10034.041152263377_dp)*Hr2(0,1) + &
                        z**3*(21576.69135802469_dp)*Hr2(0,1) + z*(36899.03848254654_dp)*Hr2(0,1) + (-139971.77125812223_dp)*Hr2(1,0) &
                        + z**2*(-17849.721020007615_dp)*Hr2(1,0) + nf*z*(-6503.769547325104_dp)*Hr2(1,0) + &
                        nf*z**3*(-521.8326474622771_dp)*Hr2(1,0) + nf**2*(-68.65294924554183_dp)*Hr2(1,0) + &
                        nf**2*z**2*(-34.32647462277092_dp)*Hr2(1,0) + nf**2*z*(68.65294924554183_dp)*Hr2(1,0) + &
                        nf*z**2*(1682.6995884773664_dp)*Hr2(1,0) + nf*(6309.750342935528_dp)*Hr2(1,0) + &
                        z**3*(18255.012345679013_dp)*Hr2(1,0) + z*(138141.5819577107_dp)*Hr2(1,0) + ((-269217.8322928852_dp) + &
                        z**2*(-30336.586928335608_dp) + nf**2*((-137.30589849108367_dp) + z**2*(-68.65294924554183_dp) + &
                        z*(137.30589849108367_dp)) + nf*(z*(-13007.539094650208_dp) + z**3*(-1043.6652949245542_dp) + &
                        z**2*(3365.399176954733_dp) + (12619.500685871057_dp)) + z**3*(36510.02469135803_dp) + &
                        z*(265557.4536920622_dp))*Hr2(1,1) + (z**2*(-23316.92949245542_dp) + (-12518.4_dp) + &
                        z*(-4260.977777777778_dp) + nf**2*z*(z*(-54.781893004115226_dp) + (-37.08312757201646_dp)) + &
                        z**3*(531.3843621399177_dp) + nf*((443.7333333333334_dp) + z*(1050.5481481481481_dp) + &
                        z**2*(1869.7481481481482_dp)))*Hr3(0,0,0) + (-43434.666666666664_dp)*Hr3(0,0,1) + &
                        z**2*(-37320.5157750343_dp)*Hr3(0,0,1) + nf*z*(-2427.2592592592587_dp)*Hr3(0,0,1) + &
                        nf**2*z**2*(-40.03292181069959_dp)*Hr3(0,0,1) + nf**2*z*(80.06584362139918_dp)*Hr3(0,0,1) + &
                        nf*(1434.8641975308642_dp)*Hr3(0,0,1) + z**3*(2098.567901234568_dp)*Hr3(0,0,1) + &
                        nf*z**2*(2322.9629629629635_dp)*Hr3(0,0,1) + z*(28072.55967078189_dp)*Hr3(0,0,1) + &
                        (-52451.555555555555_dp)*Hr3(0,1,0) + z**2*(-20038.057613168723_dp)*Hr3(0,1,0) + &
                        z*(-17367.9670781893_dp)*Hr3(0,1,0) + nf*z*(1099.851851851852_dp)*Hr3(0,1,0) + &
                        nf*z**2*(1156.7407407407409_dp)*Hr3(0,1,0) + nf*(1565.4979423868313_dp)*Hr3(0,1,0) + &
                        z**3*(3046.716049382716_dp)*Hr3(0,1,0) + (-104903.11111111111_dp)*Hr3(0,1,1) + &
                        z**2*(-40076.115226337446_dp)*Hr3(0,1,1) + z*(-34735.9341563786_dp)*Hr3(0,1,1) + &
                        nf*z*(2199.703703703704_dp)*Hr3(0,1,1) + nf*z**2*(2313.4814814814818_dp)*Hr3(0,1,1) + &
                        nf*(3130.9958847736625_dp)*Hr3(0,1,1) + z**3*(6093.432098765432_dp)*Hr3(0,1,1) + &
                        (-23546.73251028806_dp)*Hr3(1,0,0) + z**2*(-3654.7599451303154_dp)*Hr3(1,0,0) + &
                        nf*z*(-587.851851851852_dp)*Hr3(1,0,0) + nf*z**2*(293.925925925926_dp)*Hr3(1,0,0) + &
                        nf*(587.851851851852_dp)*Hr3(1,0,0) + z**3*(1786.732510288066_dp)*Hr3(1,0,0) + &
                        z*(21000.0768175583_dp)*Hr3(1,0,0) + (-70640.1975308642_dp)*Hr3(1,0,1) + &
                        z**2*(-10964.279835390947_dp)*Hr3(1,0,1) + nf*z*(-1763.5555555555559_dp)*Hr3(1,0,1) + &
                        nf*z**2*(881.7777777777779_dp)*Hr3(1,0,1) + nf*(1763.5555555555559_dp)*Hr3(1,0,1) + &
                        z**3*(5360.197530864198_dp)*Hr3(1,0,1) + z*(63000.2304526749_dp)*Hr3(1,0,1) + &
                        (-70640.1975308642_dp)*Hr3(1,1,0) + z**2*(-10964.279835390947_dp)*Hr3(1,1,0) + &
                        nf*z*(-1763.5555555555559_dp)*Hr3(1,1,0) + nf*z**2*(881.7777777777779_dp)*Hr3(1,1,0) + &
                        nf*(1763.5555555555559_dp)*Hr3(1,1,0) + z**3*(5360.197530864198_dp)*Hr3(1,1,0) + &
                        z*(63000.2304526749_dp)*Hr3(1,1,0) + (-141280.3950617284_dp)*Hr3(1,1,1) + &
                        z**2*(-21928.559670781895_dp)*Hr3(1,1,1) + nf*z*(-3527.1111111111118_dp)*Hr3(1,1,1) + &
                        nf*z**2*(1763.5555555555559_dp)*Hr3(1,1,1) + nf*(3527.1111111111118_dp)*Hr3(1,1,1) + &
                        z**3*(10720.395061728395_dp)*Hr3(1,1,1) + z*(126000.4609053498_dp)*Hr3(1,1,1) + &
                        z**2*(-8294.259533607681_dp)*Hr4(0,0,0,0) + (-921.6_dp)*Hr4(0,0,0,0) + &
                        nf*z*(-154.2320987654321_dp)*Hr4(0,0,0,0) + nf**2*z**2*(-2.528395061728395_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(5.05679012345679_dp)*Hr4(0,0,0,0) + nf*z**2*(475.3382716049383_dp)*Hr4(0,0,0,0) + &
                        z*(754.4449931412894_dp)*Hr4(0,0,0,0) + z**2*(-18142.639231824418_dp)*Hr4(0,0,0,1) + &
                        (-4096._dp)*Hr4(0,0,0,1) + z*(-3423.1659807956103_dp)*Hr4(0,0,0,1) + &
                        nf*z*(429.82716049382714_dp)*Hr4(0,0,0,1) + nf*z**2*(581.5308641975308_dp)*Hr4(0,0,0,1) + &
                        z**2*(-15060.80658436214_dp)*Hr4(0,0,1,0) + (-6940.444444444444_dp)*Hr4(0,0,1,0) + &
                        nf*z*(-366.6172839506173_dp)*Hr4(0,0,1,0) + nf*z**2*(183.30864197530866_dp)*Hr4(0,0,1,0) + &
                        z*(6000.724279835392_dp)*Hr4(0,0,1,0) + z**2*(-30121.61316872428_dp)*Hr4(0,0,1,1) + &
                        (-13880.888888888889_dp)*Hr4(0,0,1,1) + nf*z*(-733.2345679012346_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(366.6172839506173_dp)*Hr4(0,0,1,1) + z*(12001.448559670784_dp)*Hr4(0,0,1,1) + &
                        z**2*(-5650.260631001372_dp)*Hr4(0,1,0,0) + (-5360.197530864198_dp)*Hr4(0,1,0,0) + &
                        z*(-4780.07133058985_dp)*Hr4(0,1,0,0) + z**2*(-16950.781893004114_dp)*Hr4(0,1,0,1) + &
                        (-16080.592592592593_dp)*Hr4(0,1,0,1) + z*(-14340.21399176955_dp)*Hr4(0,1,0,1) + &
                        z**2*(-16950.781893004114_dp)*Hr4(0,1,1,0) + (-16080.592592592593_dp)*Hr4(0,1,1,0) + &
                        z*(-14340.21399176955_dp)*Hr4(0,1,1,0) + z**2*(-33901.56378600823_dp)*Hr4(0,1,1,1) + &
                        (-32161.185185185186_dp)*Hr4(0,1,1,1) + z*(-28680.4279835391_dp)*Hr4(0,1,1,1) + &
                        (-1630.1124828532236_dp)*Hr4(1,0,0,0) + z**2*(-815.0562414266118_dp)*Hr4(1,0,0,0) + &
                        z*(1630.1124828532236_dp)*Hr4(1,0,0,0) + (-6520.449931412894_dp)*Hr4(1,0,0,1) + &
                        z**2*(-3260.224965706447_dp)*Hr4(1,0,0,1) + z*(6520.449931412894_dp)*Hr4(1,0,0,1) + &
                        (-9780.674897119343_dp)*Hr4(1,0,1,0) + z**2*(-4890.337448559671_dp)*Hr4(1,0,1,0) + &
                        z*(9780.674897119343_dp)*Hr4(1,0,1,0) + (-19561.349794238686_dp)*Hr4(1,0,1,1) + &
                        z**2*(-9780.674897119343_dp)*Hr4(1,0,1,1) + z*(19561.349794238686_dp)*Hr4(1,0,1,1) + &
                        (-6520.449931412894_dp)*Hr4(1,1,0,0) + z**2*(-3260.224965706447_dp)*Hr4(1,1,0,0) + &
                        z*(6520.449931412894_dp)*Hr4(1,1,0,0) + (-19561.349794238686_dp)*Hr4(1,1,0,1) + &
                        z**2*(-9780.674897119343_dp)*Hr4(1,1,0,1) + z*(19561.349794238686_dp)*Hr4(1,1,0,1) + &
                        (-19561.349794238686_dp)*Hr4(1,1,1,0) + z**2*(-9780.674897119343_dp)*Hr4(1,1,1,0) + &
                        z*(19561.349794238686_dp)*Hr4(1,1,1,0) + (-39122.69958847737_dp)*Hr4(1,1,1,1) + &
                        z**2*(-19561.349794238686_dp)*Hr4(1,1,1,1) + z*(39122.69958847737_dp)*Hr4(1,1,1,1))/z

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
                        Ibar_select = ((-3.6425099768754384e6_dp) + nf*z*(-1.1693214691542585e6_dp) + &
                        nf*z**3*(-111694.71525304312_dp) + nf*z**2*(-59953.191215402185_dp) + nf**2*(-29933.41811392497_dp) + &
                        nf**3*z*(-1003.797957921116_dp) + nf**3*z**3*(-104.68062795305592_dp) + nf**4*z**2*(-5.068495656149977_dp) + &
                        nf**4*(-4.9592440176802315_dp) + nf**4*z**3*(0.6399024538942234_dp) + nf**4*z*(9.739003200731595_dp) + &
                        nf**2*z**2*(261.3982322402335_dp) + nf**3*z**2*(271.28310116136674_dp) + nf**3*(652.3175502604938_dp) + &
                        nf**2*z**3*(5631.418082517358_dp) + nf**2*z*(52917.63993903681_dp) + z**2*(348492.107544505_dp) + &
                        z**3*(553488.3084905705_dp) + nf*(564476.6865313558_dp) + z*(7.948459623871852e6_dp) + &
                        (-2.391427371320574e6_dp)*Hr1(0) + nf*z**3*(-62398.443249984644_dp)*Hr1(0) + &
                        nf*z*(-19033.502663502848_dp)*Hr1(0) + nf**2*(-15612.142522694416_dp)*Hr1(0) + &
                        nf*z**2*(-10990.190737051858_dp)*Hr1(0) + nf**2*z*(-7328.892315564053_dp)*Hr1(0) + &
                        nf**2*z**2*(-2761.9761636458284_dp)*Hr1(0) + nf**3*z**3*(-13.591424071533812_dp)*Hr1(0) + &
                        nf**4*(-1.91970736168267_dp)*Hr1(0) + nf**4*z*(0.046822130772748057_dp)*Hr1(0) + &
                        nf**4*z**2*(1.0300868770004572_dp)*Hr1(0) + nf**3*z**2*(69.24786778129695_dp)*Hr1(0) + &
                        nf**3*z*(203.4735983804026_dp)*Hr1(0) + nf**3*(299.16024996189606_dp)*Hr1(0) + &
                        nf**2*z**3*(1653.3635726261239_dp)*Hr1(0) + z**2*(297454.2938255803_dp)*Hr1(0) + &
                        nf*(329859.45697906485_dp)*Hr1(0) + z**3*(718097.3971565821_dp)*Hr1(0) + z*(1.9567645969020291e6_dp)*Hr1(0) &
                        + (-5.733856397189344e6_dp)*Hr1(1) + z**2*(-856111.5361314359_dp)*Hr1(1) + &
                        nf*z*(-660576.5919193014_dp)*Hr1(1) + nf*z**3*(-99567.64693141526_dp)*Hr1(1) + &
                        nf**2*(-22801.96522766523_dp)*Hr1(1) + nf**2*z**2*(-5147.557831025117_dp)*Hr1(1) + &
                        nf**3*z*(-423.0067367779302_dp)*Hr1(1) + nf**3*z**3*(-28.473057968805573_dp)*Hr1(1) + &
                        nf**4*(-1.3578417924096935_dp)*Hr1(1) + nf**4*z**2*(-0.6789208962048467_dp)*Hr1(1) + &
                        nf**4*z*(1.3578417924096935_dp)*Hr1(1) + nf**3*z**2*(173.7374180765127_dp)*Hr1(1) + &
                        nf**3*(310.2837575572829_dp)*Hr1(1) + nf**2*z**3*(3047.6154244779764_dp)*Hr1(1) + &
                        nf**2*z*(25378.447593158293_dp)*Hr1(1) + nf*z**2*(84563.60185211537_dp)*Hr1(1) + &
                        nf*(643627.2851825607_dp)*Hr1(1) + z**3*(963156.0076171076_dp)*Hr1(1) + z*(5.890526534415294e6_dp)*Hr1(1) + &
                        (-1.056344224906906e6_dp)*Hr2(0,0) + nf*z*(-128056.12046323087_dp)*Hr2(0,0) + &
                        z**2*(-56633.62576796685_dp)*Hr2(0,0) + nf*z**3*(-10437.68303612254_dp)*Hr2(0,0) + &
                        nf**2*(-4056.265569272977_dp)*Hr2(0,0) + nf**2*z**2*(-2737.8416862398435_dp)*Hr2(0,0) + &
                        nf**3*z*(-178.80591373266267_dp)*Hr2(0,0) + nf**4*z**2*(-1.1237311385459534_dp)*Hr2(0,0) + &
                        nf**4*z*(2.247462277091907_dp)*Hr2(0,0) + nf**3*(44.644901691815264_dp)*Hr2(0,0) + &
                        nf**3*z**2*(61.97689376619417_dp)*Hr2(0,0) + nf**2*z**3*(137.8521566834324_dp)*Hr2(0,0) + &
                        nf**2*z*(6417.160263646571_dp)*Hr2(0,0) + nf*z**2*(37024.2461139215_dp)*Hr2(0,0) + &
                        nf*(116488.81064561542_dp)*Hr2(0,0) + z**3*(191974.69183225045_dp)*Hr2(0,0) + &
                        z*(1.2185530672428133e6_dp)*Hr2(0,0) + (-2.5748317499469365e6_dp)*Hr2(0,1) + &
                        z**2*(-503846.3861492651_dp)*Hr2(0,1) + nf*z**3*(-30166.16442615455_dp)*Hr2(0,1) + &
                        nf**2*(-8751.264990092974_dp)*Hr2(0,1) + nf**2*z**2*(-4059.8416095107455_dp)*Hr2(0,1) + &
                        nf**2*z*(-3804.5907636031093_dp)*Hr2(0,1) + nf**3*z**2*(7.897332723670172_dp)*Hr2(0,1) + &
                        nf**3*z*(28.733180917543056_dp)*Hr2(0,1) + nf**3*(85.4191739064167_dp)*Hr2(0,1) + &
                        nf**2*z**3*(414.30562414266126_dp)*Hr2(0,1) + nf*z*(5835.664404396488_dp)*Hr2(0,1) + &
                        nf*z**2*(78198.38431446634_dp)*Hr2(0,1) + nf*(275407.3621934227_dp)*Hr2(0,1) + &
                        z**3*(507586.2266415312_dp)*Hr2(0,1) + z*(1.7242774644460494e6_dp)*Hr2(0,1) + &
                        (-2.8851616024907255e6_dp)*Hr2(1,0) + z**2*(-405056.1295953139_dp)*Hr2(1,0) + &
                        nf*z*(-217932.25935422644_dp)*Hr2(1,0) + nf*z**3*(-26730.871513488797_dp)*Hr2(1,0) + &
                        nf**2*(-5048.430422191738_dp)*Hr2(1,0) + nf**2*z**2*(-1841.9240969364425_dp)*Hr2(1,0) + &
                        nf**3*z*(-38.90138698369151_dp)*Hr2(1,0) + nf**3*z**2*(19.450693491845755_dp)*Hr2(1,0) + &
                        nf**3*(38.90138698369151_dp)*Hr2(1,0) + nf**2*z**3*(402.9434537418077_dp)*Hr2(1,0) + &
                        nf**2*z*(5656.991312299953_dp)*Hr2(1,0) + nf*z**2*(34945.41332289892_dp)*Hr2(1,0) + &
                        nf*(215113.635363981_dp)*Hr2(1,0) + z**3*(426442.71521215077_dp)*Hr2(1,0) + &
                        z*(2.9173849837141745e6_dp)*Hr2(1,0) + ((-5.28247753582642e6_dp) + z**2*(-705965.4949231313_dp) + &
                        nf**3*(z*(-77.80277396738302_dp) + z**2*(38.90138698369151_dp) + (77.80277396738302_dp)) + &
                        nf**2*((-10096.860844383476_dp) + z**2*(-3683.848193872885_dp) + z**3*(805.8869074836153_dp) + &
                        z*(11313.982624599907_dp)) + nf*(z*(-423098.8234825848_dp) + z**3*(-53461.743026977594_dp) + &
                        z**2*(63507.97903286382_dp) + (417461.57550209406_dp)) + z**3*(822222.4995459174_dp) + &
                        z*(5.392424487300456e6_dp))*Hr2(1,1) + (-268177.55765236466_dp)*Hr3(0,0,0) + &
                        z**2*(-257755.75252890596_dp)*Hr3(0,0,0) + z*(-31824.110266604584_dp)*Hr3(0,0,0) + &
                        nf**2*z**2*(-1854.7418514224746_dp)*Hr3(0,0,0) + nf**2*z*(-1462.137696332005_dp)*Hr3(0,0,0) + &
                        nf*z**3*(-638.7266981659299_dp)*Hr3(0,0,0) + nf**2*(-379.4699588477366_dp)*Hr3(0,0,0) + &
                        nf**3*z*(24.53479652491998_dp)*Hr3(0,0,0) + nf**3*z**2*(37.45770461819845_dp)*Hr3(0,0,0) + &
                        nf*(20544.47407407407_dp)*Hr3(0,0,0) + nf*z*(22324.402070613778_dp)*Hr3(0,0,0) + &
                        z**3*(22734.20466392318_dp)*Hr3(0,0,0) + nf*z**2*(39905.94428214506_dp)*Hr3(0,0,0) + &
                        (-890473.0853798247_dp)*Hr3(0,0,1) + z**2*(-440364.5315108934_dp)*Hr3(0,0,1) + &
                        nf*z*(-61461.39918224429_dp)*Hr3(0,0,1) + nf*z**3*(-2541.255540313976_dp)*Hr3(0,0,1) + &
                        nf**2*z**2*(-1830.792135345222_dp)*Hr3(0,0,1) + nf**2*(-1174.7438500228623_dp)*Hr3(0,0,1) + &
                        nf**3*z*(-60.49419295839049_dp)*Hr3(0,0,1) + nf**3*z**2*(30.247096479195246_dp)*Hr3(0,0,1) + &
                        nf**2*z*(3106.812254229538_dp)*Hr3(0,0,1) + nf*z**2*(65055.906114517966_dp)*Hr3(0,0,1) + &
                        nf*(67909.53086419753_dp)*Hr3(0,0,1) + z**3*(88944.93497942387_dp)*Hr3(0,0,1) + &
                        z*(579184.0094257747_dp)*Hr3(0,0,1) + (-1.1875031579821315e6_dp)*Hr3(0,1,0) + &
                        z**2*(-330876.23233010876_dp)*Hr3(0,1,0) + nf*z**3*(-3755.4470355128797_dp)*Hr3(0,1,0) + &
                        nf**2*(-1208.8303612254228_dp)*Hr3(0,1,0) + nf**2*z*(-718.6026520347508_dp)*Hr3(0,1,0) + &
                        nf**2*z**2*(-630.4599908550525_dp)*Hr3(0,1,0) + z*(24584.76178530643_dp)*Hr3(0,1,0) + &
                        nf*z*(34067.89940557842_dp)*Hr3(0,1,0) + nf*z**2*(36088.918152720624_dp)*Hr3(0,1,0) + &
                        nf*(77001.21810699588_dp)*Hr3(0,1,0) + z**3*(127334.18930041154_dp)*Hr3(0,1,0) + &
                        (-2.2830175233291104e6_dp)*Hr3(0,1,1) + z**2*(-565683.7020200472_dp)*Hr3(0,1,1) + &
                        nf*z**3*(-7510.894071025759_dp)*Hr3(0,1,1) + nf**2*(-2417.6607224508457_dp)*Hr3(0,1,1) + &
                        nf**2*z*(-1437.2053040695016_dp)*Hr3(0,1,1) + nf**2*z**2*(-1260.919981710105_dp)*Hr3(0,1,1) + &
                        nf*z*(68135.79881115684_dp)*Hr3(0,1,1) + nf*z**2*(72177.83630544125_dp)*Hr3(0,1,1) + &
                        z*(132998.3761957303_dp)*Hr3(0,1,1) + nf*(154002.43621399175_dp)*Hr3(0,1,1) + &
                        z**3*(254668.37860082308_dp)*Hr3(0,1,1) + (-618115.0842367429_dp)*Hr3(1,0,0) + &
                        z**2*(-86328.19720143806_dp)*Hr3(1,0,0) + nf*z*(-31110.30147843316_dp)*Hr3(1,0,0) + &
                        nf*z**3*(-2315.614489661129_dp)*Hr3(1,0,0) + nf**2*(-416.05365035817704_dp)*Hr3(1,0,0) + &
                        nf**2*z**2*(-208.02682517908852_dp)*Hr3(1,0,0) + nf**2*z*(416.05365035817704_dp)*Hr3(1,0,0) + &
                        nf*z**2*(8029.858862978203_dp)*Hr3(1,0,0) + nf*(31675.782756693592_dp)*Hr3(1,0,0) + &
                        z**3*(73053.05898491084_dp)*Hr3(1,0,0) + z*(597498.4416517712_dp)*Hr3(1,0,0) + &
                        (-1.7731137482188116e6_dp)*Hr3(1,0,1) + z**2*(-218368.83935860562_dp)*Hr3(1,0,1) + &
                        nf*z*(-93330.90443529948_dp)*Hr3(1,0,1) + nf*z**3*(-6946.843468983387_dp)*Hr3(1,0,1) + &
                        nf**2*(-1248.1609510745313_dp)*Hr3(1,0,1) + nf**2*z**2*(-624.0804755372657_dp)*Hr3(1,0,1) + &
                        nf**2*z*(1248.1609510745313_dp)*Hr3(1,0,1) + nf*z**2*(24089.576588934615_dp)*Hr3(1,0,1) + &
                        nf*(95027.34827008078_dp)*Hr3(1,0,1) + z**3*(219159.17695473254_dp)*Hr3(1,0,1) + &
                        z*(1.7112638204638963e6_dp)*Hr3(1,0,1) + (-1.8272680845464228e6_dp)*Hr3(1,1,0) + &
                        z**2*(-245446.00752241138_dp)*Hr3(1,1,0) + nf*z*(-93330.90443529948_dp)*Hr3(1,1,0) + &
                        nf*z**3*(-6946.843468983387_dp)*Hr3(1,1,0) + nf**2*(-1248.1609510745313_dp)*Hr3(1,1,0) + &
                        nf**2*z**2*(-624.0804755372657_dp)*Hr3(1,1,0) + nf**2*z*(1248.1609510745313_dp)*Hr3(1,1,0) + &
                        nf*z**2*(24089.576588934615_dp)*Hr3(1,1,0) + nf*(95027.34827008078_dp)*Hr3(1,1,0) + &
                        z**3*(219159.17695473254_dp)*Hr3(1,1,0) + z*(1.7654181567915077e6_dp)*Hr3(1,1,0) + &
                        (-3.546227496437623e6_dp)*Hr3(1,1,1) + z**2*(-436737.67871721124_dp)*Hr3(1,1,1) + &
                        nf*z*(-186661.80887059896_dp)*Hr3(1,1,1) + nf*z**3*(-13893.686937966773_dp)*Hr3(1,1,1) + &
                        nf**2*(-2496.3219021490627_dp)*Hr3(1,1,1) + nf**2*z**2*(-1248.1609510745313_dp)*Hr3(1,1,1) + &
                        nf**2*z*(2496.3219021490627_dp)*Hr3(1,1,1) + nf*z**2*(48179.15317786923_dp)*Hr3(1,1,1) + &
                        nf*(190054.69654016156_dp)*Hr3(1,1,1) + z**3*(438318.3539094651_dp)*Hr3(1,1,1) + &
                        z*(3.4225276409277925e6_dp)*Hr3(1,1,1) + z**2*(-160791.27059899404_dp)*Hr4(0,0,0,0) + &
                        (-35174.4_dp)*Hr4(0,0,0,0) + nf*z*(-4457.2639536656_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**2*(-582.9355281207133_dp)*Hr4(0,0,0,0) + nf**3*z*(-6.74238683127572_dp)*Hr4(0,0,0,0) + &
                        nf**3*z**2*(3.37119341563786_dp)*Hr4(0,0,0,0) + nf**2*z*(295.26035665294927_dp)*Hr4(0,0,0,0) + &
                        z**3*(1086.7416552354823_dp)*Hr4(0,0,0,0) + nf*(1388.0888888888887_dp)*Hr4(0,0,0,0) + &
                        nf*z**2*(16762.19795762841_dp)*Hr4(0,0,0,0) + z*(21399.617863130617_dp)*Hr4(0,0,0,0) + &
                        z**2*(-340813.87300716346_dp)*Hr4(0,0,0,1) + (-156962.13333333333_dp)*Hr4(0,0,0,1) + &
                        z*(-66647.18280749886_dp)*Hr4(0,0,0,1) + nf**2*z**2*(-688.2853223593963_dp)*Hr4(0,0,0,1) + &
                        nf**2*z*(-501.18408779149513_dp)*Hr4(0,0,0,1) + z**3*(5409.735345221765_dp)*Hr4(0,0,0,1) + &
                        nf*(5931.614814814814_dp)*Hr4(0,0,0,1) + nf*z*(13953.494406340496_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(27663.107940862672_dp)*Hr4(0,0,0,1) + z**2*(-271123.62627648225_dp)*Hr4(0,0,1,0) + &
                        (-268117.3333333333_dp)*Hr4(0,0,1,0) + nf*z*(-15529.730224051209_dp)*Hr4(0,0,1,0) + &
                        nf**2*z**2*(-235.28120713305896_dp)*Hr4(0,0,1,0) + nf**2*z*(470.5624142661179_dp)*Hr4(0,0,1,0) + &
                        nf*(9574.189300411523_dp)*Hr4(0,0,1,0) + z**3*(10717.58573388203_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(17610.856881572927_dp)*Hr4(0,0,1,0) + z*(181238.40146319158_dp)*Hr4(0,0,1,0) + &
                        z**2*(-542247.2525529645_dp)*Hr4(0,0,1,1) + (-536234.6666666666_dp)*Hr4(0,0,1,1) + &
                        nf*z*(-31059.460448102418_dp)*Hr4(0,0,1,1) + nf**2*z**2*(-470.5624142661179_dp)*Hr4(0,0,1,1) + &
                        nf**2*z*(941.1248285322358_dp)*Hr4(0,0,1,1) + nf*(19148.378600823045_dp)*Hr4(0,0,1,1) + &
                        z**3*(21435.17146776406_dp)*Hr4(0,0,1,1) + nf*z**2*(35221.71376314585_dp)*Hr4(0,0,1,1) + &
                        z*(362476.80292638316_dp)*Hr4(0,0,1,1) + (-210958.74897119342_dp)*Hr4(0,1,0,0) + &
                        z**2*(-95999.37570492302_dp)*Hr4(0,1,0,0) + z*(-88969.93019356806_dp)*Hr4(0,1,0,0) + &
                        nf*z*(5419.973784484073_dp)*Hr4(0,1,0,0) + nf*z**2*(5750.303917085811_dp)*Hr4(0,1,0,0) + &
                        nf*(6946.843468983387_dp)*Hr4(0,1,0,0) + z**3*(10440.39871970736_dp)*Hr4(0,1,0,0) + &
                        (-632876.2469135802_dp)*Hr4(0,1,0,1) + z**2*(-287998.1271147691_dp)*Hr4(0,1,0,1) + &
                        z*(-266909.79058070417_dp)*Hr4(0,1,0,1) + nf*z*(16259.921353452217_dp)*Hr4(0,1,0,1) + &
                        nf*z**2*(17250.91175125743_dp)*Hr4(0,1,0,1) + nf*(20840.53040695016_dp)*Hr4(0,1,0,1) + &
                        z**3*(31321.196159122082_dp)*Hr4(0,1,0,1) + (-632876.2469135802_dp)*Hr4(0,1,1,0) + &
                        z**2*(-287998.1271147691_dp)*Hr4(0,1,1,0) + z*(-266909.79058070417_dp)*Hr4(0,1,1,0) + &
                        nf*z*(16259.921353452217_dp)*Hr4(0,1,1,0) + nf*z**2*(17250.91175125743_dp)*Hr4(0,1,1,0) + &
                        nf*(20840.53040695016_dp)*Hr4(0,1,1,0) + z**3*(31321.196159122082_dp)*Hr4(0,1,1,0) + &
                        (-1.2657524938271604e6_dp)*Hr4(0,1,1,1) + z**2*(-575996.2542295381_dp)*Hr4(0,1,1,1) + &
                        z*(-533819.5811614083_dp)*Hr4(0,1,1,1) + nf*z*(32519.842706904434_dp)*Hr4(0,1,1,1) + &
                        nf*z**2*(34501.82350251486_dp)*Hr4(0,1,1,1) + nf*(41681.06081390032_dp)*Hr4(0,1,1,1) + &
                        z**3*(62642.392318244165_dp)*Hr4(0,1,1,1) + (-68129.12665752172_dp)*Hr4(1,0,0,0) + &
                        z**2*(-12821.147995732357_dp)*Hr4(1,0,0,0) + nf*z*(-1940.1530254534368_dp)*Hr4(1,0,0,0) + &
                        nf*z**2*(970.0765127267184_dp)*Hr4(1,0,0,0) + nf*(1940.1530254534368_dp)*Hr4(1,0,0,0) + &
                        z**3*(4660.206675811614_dp)*Hr4(1,0,0,0) + z*(61213.927145252244_dp)*Hr4(1,0,0,0) + &
                        (-272516.50663008686_dp)*Hr4(1,0,0,1) + z**2*(-51284.59198292943_dp)*Hr4(1,0,0,1) + &
                        nf*z*(-7760.612101813747_dp)*Hr4(1,0,0,1) + nf*z**2*(3880.3060509068737_dp)*Hr4(1,0,0,1) + &
                        nf*(7760.612101813747_dp)*Hr4(1,0,0,1) + z**3*(18640.826703246457_dp)*Hr4(1,0,0,1) + &
                        z*(244855.70858100898_dp)*Hr4(1,0,0,1) + (-408774.7599451303_dp)*Hr4(1,0,1,0) + &
                        z**2*(-76926.88797439415_dp)*Hr4(1,0,1,0) + nf*z*(-11640.918152720622_dp)*Hr4(1,0,1,0) + &
                        nf*z**2*(5820.459076360311_dp)*Hr4(1,0,1,0) + nf*(11640.918152720622_dp)*Hr4(1,0,1,0) + &
                        z**3*(27961.24005486969_dp)*Hr4(1,0,1,0) + z*(367283.5628715135_dp)*Hr4(1,0,1,0) + &
                        (-817549.5198902606_dp)*Hr4(1,0,1,1) + z**2*(-153853.7759487883_dp)*Hr4(1,0,1,1) + &
                        nf*z*(-23281.836305441244_dp)*Hr4(1,0,1,1) + nf*z**2*(11640.918152720622_dp)*Hr4(1,0,1,1) + &
                        nf*(23281.836305441244_dp)*Hr4(1,0,1,1) + z**3*(55922.48010973938_dp)*Hr4(1,0,1,1) + &
                        z*(734567.125743027_dp)*Hr4(1,0,1,1) + (-272516.50663008686_dp)*Hr4(1,1,0,0) + &
                        z**2*(-51284.59198292943_dp)*Hr4(1,1,0,0) + nf*z*(-7760.612101813747_dp)*Hr4(1,1,0,0) + &
                        nf*z**2*(3880.3060509068737_dp)*Hr4(1,1,0,0) + nf*(7760.612101813747_dp)*Hr4(1,1,0,0) + &
                        z**3*(18640.826703246457_dp)*Hr4(1,1,0,0) + z*(244855.70858100898_dp)*Hr4(1,1,0,0) + &
                        (-817549.5198902606_dp)*Hr4(1,1,0,1) + z**2*(-153853.7759487883_dp)*Hr4(1,1,0,1) + &
                        nf*z*(-23281.836305441244_dp)*Hr4(1,1,0,1) + nf*z**2*(11640.918152720622_dp)*Hr4(1,1,0,1) + &
                        nf*(23281.836305441244_dp)*Hr4(1,1,0,1) + z**3*(55922.48010973938_dp)*Hr4(1,1,0,1) + &
                        z*(734567.125743027_dp)*Hr4(1,1,0,1) + (-817549.5198902606_dp)*Hr4(1,1,1,0) + &
                        z**2*(-153853.7759487883_dp)*Hr4(1,1,1,0) + nf*z*(-23281.836305441244_dp)*Hr4(1,1,1,0) + &
                        nf*z**2*(11640.918152720622_dp)*Hr4(1,1,1,0) + nf*(23281.836305441244_dp)*Hr4(1,1,1,0) + &
                        z**3*(55922.48010973938_dp)*Hr4(1,1,1,0) + z*(734567.125743027_dp)*Hr4(1,1,1,0) + &
                        (-1.6350990397805213e6_dp)*Hr4(1,1,1,1) + z**2*(-307707.5518975766_dp)*Hr4(1,1,1,1) + &
                        nf*z*(-46563.67261088249_dp)*Hr4(1,1,1,1) + nf*z**2*(23281.836305441244_dp)*Hr4(1,1,1,1) + &
                        nf*(46563.67261088249_dp)*Hr4(1,1,1,1) + z**3*(111844.96021947876_dp)*Hr4(1,1,1,1) + &
                        z*(1.469134251486054e6_dp)*Hr4(1,1,1,1) + z**2*(-33177.537570492306_dp)*Hr5(0,0,0,0,0) + &
                        (-1843.2_dp)*Hr5(0,0,0,0,0) + z*(-1507.8911141594267_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z**2*(-33.711934156378604_dp)*Hr5(0,0,0,0,0) + nf**2*z*(-23.598353909465022_dp)*Hr5(0,0,0,0,0) + &
                        nf*z*(398.54997713763146_dp)*Hr5(0,0,0,0,0) + nf*z**2*(2592.073159579332_dp)*Hr5(0,0,0,0,0) + &
                        z**2*(-89164.07035512879_dp)*Hr5(0,0,0,0,1) + (-10035.2_dp)*Hr5(0,0,0,0,1) + &
                        nf*z*(-1552.9964334705076_dp)*Hr5(0,0,0,0,1) + nf**2*z**2*(-21.912757201646095_dp)*Hr5(0,0,0,0,1) + &
                        nf**2*z*(43.82551440329219_dp)*Hr5(0,0,0,0,1) + nf*z**2*(4569.090809327847_dp)*Hr5(0,0,0,0,1) + &
                        z*(8365.21066910532_dp)*Hr5(0,0,0,0,1) + z**2*(-96548.48224356043_dp)*Hr5(0,0,0,1,0) + &
                        (-22072.888888888894_dp)*Hr5(0,0,0,1,0) + z*(-18807.82563633592_dp)*Hr5(0,0,0,1,0) + &
                        nf*z*(2112.614540466392_dp)*Hr5(0,0,0,1,0) + nf*z**2*(2736.2853223593966_dp)*Hr5(0,0,0,1,0) + &
                        z**2*(-193096.96448712086_dp)*Hr5(0,0,0,1,1) + (-44145.77777777779_dp)*Hr5(0,0,0,1,1) + &
                        z*(-37615.65127267184_dp)*Hr5(0,0,0,1,1) + nf*z*(4225.229080932784_dp)*Hr5(0,0,0,1,1) + &
                        nf*z**2*(5472.570644718793_dp)*Hr5(0,0,0,1,1) + z**2*(-52762.61057765584_dp)*Hr5(0,0,1,0,0) + &
                        (-24601.28395061728_dp)*Hr5(0,0,1,0,0) + nf*z*(-1119.9853680841336_dp)*Hr5(0,0,1,0,0) + &
                        nf*z**2*(559.9926840420668_dp)*Hr5(0,0,1,0,0) + z*(21641.500990702632_dp)*Hr5(0,0,1,0,0) + &
                        z**2*(-158287.83173296755_dp)*Hr5(0,0,1,0,1) + (-73803.85185185185_dp)*Hr5(0,0,1,0,1) + &
                        nf*z*(-3359.9561042524006_dp)*Hr5(0,0,1,0,1) + nf*z**2*(1679.9780521262003_dp)*Hr5(0,0,1,0,1) + &
                        z*(64924.50297210792_dp)*Hr5(0,0,1,0,1) + z**2*(-158287.83173296755_dp)*Hr5(0,0,1,1,0) + &
                        (-73803.85185185185_dp)*Hr5(0,0,1,1,0) + nf*z*(-3359.9561042524006_dp)*Hr5(0,0,1,1,0) + &
                        nf*z**2*(1679.9780521262003_dp)*Hr5(0,0,1,1,0) + z*(64924.50297210792_dp)*Hr5(0,0,1,1,0) + &
                        z**2*(-316575.6634659351_dp)*Hr5(0,0,1,1,1) + (-147607.7037037037_dp)*Hr5(0,0,1,1,1) + &
                        nf*z*(-6719.912208504801_dp)*Hr5(0,0,1,1,1) + nf*z**2*(3359.9561042524006_dp)*Hr5(0,0,1,1,1) + &
                        z*(129849.00594421584_dp)*Hr5(0,0,1,1,1) + z**2*(-14600.701112635268_dp)*Hr5(0,1,0,0,0) + &
                        (-13980.620027434845_dp)*Hr5(0,1,0,0,0) + z*(-12740.45785703399_dp)*Hr5(0,1,0,0,0) + &
                        z**2*(-58402.80445054107_dp)*Hr5(0,1,0,0,1) + (-55922.48010973938_dp)*Hr5(0,1,0,0,1) + &
                        z*(-50961.83142813596_dp)*Hr5(0,1,0,0,1) + z**2*(-87604.20667581163_dp)*Hr5(0,1,0,1,0) + &
                        (-83883.72016460905_dp)*Hr5(0,1,0,1,0) + z*(-76442.74714220394_dp)*Hr5(0,1,0,1,0) + &
                        z**2*(-175208.41335162325_dp)*Hr5(0,1,0,1,1) + (-167767.4403292181_dp)*Hr5(0,1,0,1,1) + &
                        z*(-152885.49428440788_dp)*Hr5(0,1,0,1,1) + z**2*(-58402.80445054107_dp)*Hr5(0,1,1,0,0) + &
                        (-55922.48010973938_dp)*Hr5(0,1,1,0,0) + z*(-50961.83142813596_dp)*Hr5(0,1,1,0,0) + &
                        z**2*(-175208.41335162325_dp)*Hr5(0,1,1,0,1) + (-167767.4403292181_dp)*Hr5(0,1,1,0,1) + &
                        z*(-152885.49428440788_dp)*Hr5(0,1,1,0,1) + z**2*(-175208.41335162325_dp)*Hr5(0,1,1,1,0) + &
                        (-167767.4403292181_dp)*Hr5(0,1,1,1,0) + z*(-152885.49428440788_dp)*Hr5(0,1,1,1,0) + &
                        z**2*(-350416.8267032465_dp)*Hr5(0,1,1,1,1) + (-335534.8806584362_dp)*Hr5(0,1,1,1,1) + &
                        z*(-305770.98856881575_dp)*Hr5(0,1,1,1,1) + (-3292.1888736473097_dp)*Hr5(1,0,0,0,0) + &
                        z**2*(-1646.0944368236549_dp)*Hr5(1,0,0,0,0) + z*(3292.1888736473097_dp)*Hr5(1,0,0,0,0) + &
                        (-16460.94436823655_dp)*Hr5(1,0,0,0,1) + z**2*(-8230.472184118275_dp)*Hr5(1,0,0,0,1) + &
                        z*(16460.94436823655_dp)*Hr5(1,0,0,0,1) + (-32921.8887364731_dp)*Hr5(1,0,0,1,0) + &
                        z**2*(-16460.94436823655_dp)*Hr5(1,0,0,1,0) + z*(32921.8887364731_dp)*Hr5(1,0,0,1,0) + &
                        (-65843.7774729462_dp)*Hr5(1,0,0,1,1) + z**2*(-32921.8887364731_dp)*Hr5(1,0,0,1,1) + &
                        z*(65843.7774729462_dp)*Hr5(1,0,0,1,1) + (-32921.8887364731_dp)*Hr5(1,0,1,0,0) + &
                        z**2*(-16460.94436823655_dp)*Hr5(1,0,1,0,0) + z*(32921.8887364731_dp)*Hr5(1,0,1,0,0) + &
                        (-98765.6662094193_dp)*Hr5(1,0,1,0,1) + z**2*(-49382.83310470965_dp)*Hr5(1,0,1,0,1) + &
                        z*(98765.6662094193_dp)*Hr5(1,0,1,0,1) + (-98765.6662094193_dp)*Hr5(1,0,1,1,0) + &
                        z**2*(-49382.83310470965_dp)*Hr5(1,0,1,1,0) + z*(98765.6662094193_dp)*Hr5(1,0,1,1,0) + &
                        (-197531.3324188386_dp)*Hr5(1,0,1,1,1) + z**2*(-98765.6662094193_dp)*Hr5(1,0,1,1,1) + &
                        z*(197531.3324188386_dp)*Hr5(1,0,1,1,1) + (-16460.94436823655_dp)*Hr5(1,1,0,0,0) + &
                        z**2*(-8230.472184118275_dp)*Hr5(1,1,0,0,0) + z*(16460.94436823655_dp)*Hr5(1,1,0,0,0) + &
                        (-65843.7774729462_dp)*Hr5(1,1,0,0,1) + z**2*(-32921.8887364731_dp)*Hr5(1,1,0,0,1) + &
                        z*(65843.7774729462_dp)*Hr5(1,1,0,0,1) + (-98765.6662094193_dp)*Hr5(1,1,0,1,0) + &
                        z**2*(-49382.83310470965_dp)*Hr5(1,1,0,1,0) + z*(98765.6662094193_dp)*Hr5(1,1,0,1,0) + &
                        ((-197531.3324188386_dp) + z**2*(-98765.6662094193_dp) + z*(197531.3324188386_dp))*Hr5(1,1,0,1,1) + &
                        (-65843.7774729462_dp)*Hr5(1,1,1,0,0) + z**2*(-32921.8887364731_dp)*Hr5(1,1,1,0,0) + &
                        z*(65843.7774729462_dp)*Hr5(1,1,1,0,0) + (-197531.3324188386_dp)*Hr5(1,1,1,0,1) + &
                        z**2*(-98765.6662094193_dp)*Hr5(1,1,1,0,1) + z*(197531.3324188386_dp)*Hr5(1,1,1,0,1) + &
                        (-197531.3324188386_dp)*Hr5(1,1,1,1,0) + z**2*(-98765.6662094193_dp)*Hr5(1,1,1,1,0) + &
                        z*(197531.3324188386_dp)*Hr5(1,1,1,1,0) + (-395062.6648376772_dp)*Hr5(1,1,1,1,1) + &
                        z**2*(-197531.3324188386_dp)*Hr5(1,1,1,1,1) + z*(395062.6648376772_dp)*Hr5(1,1,1,1,1))/z

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
