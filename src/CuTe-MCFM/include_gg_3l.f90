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
                            Ibar_select = (1._dp)
                        else
                            Ibar_select = (1._dp)
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
                            Ibar_select = (-4.934802200544679_dp)
                        else
                            Ibar_select = (-4.934802200544679_dp)
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
                        Ibar_select = (-12._dp)/z + z*(-12._dp) + z**2*(12._dp) + (24._dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-11._dp) + Nf*(0.6666666666666666_dp)
                        else
                            Ibar_select = (-11._dp) + Nf*(0.6666666666666666_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (-12._dp)
                        else
                            Ibar_select = (-12._dp)
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
                        Ibar_select = z**2*(-1749.3975515463917_dp) + z**4*(-756.6756476683938_dp) + (-437.39410947073037_dp) + &
                        (-99.78215752045685_dp)/z + z**6*(-23.03903743315508_dp) + Nf*(z**2*(-29.111113251155626_dp) + &
                        z*(-18.22222222222222_dp) + z**3*(1.343964992399878e-6_dp) + (5.777777777777779_dp) + &
                        (25.11111111111111_dp)/z) + z**5*(160.2193044712562_dp) + z*(573.9208042850898_dp) + &
                        z**3*(2428.998643147897_dp) + (z**2*(-706.4640151515151_dp) + z**3*(-477.5598548972189_dp) + (-293._dp) + &
                        z*(-269._dp) + Nf*(z**3*(-1.2336221242726255e-6_dp) + z*(22.666666666666668_dp) + &
                        (24.666666666666668_dp)))*Log(z) + (z**2*(-97.54483541430193_dp) + z*(-33._dp) + (3._dp) + &
                        Nf*(z*(3.3333333333333335_dp) + (6._dp)) + z**3*(433.57529794149514_dp))*Log(z)**2 + &
                        (z**3*(-51.513090785447126_dp) + z*(-24._dp) + (-12._dp) + z**2*(-4.948320413436693_dp) + &
                        Nf*((0.8888888888888888_dp) + z*(0.8888888888888888_dp)))*Log(z)**3 + ((-11.3239413203076_dp) + &
                        z**3*(-4.982017982017982_dp) + z**2*(-1.3599053562716357_dp) + z*(Nf*(-2._dp) + &
                        (23.66586465859722_dp)))*Log(z*(-1._dp) + (1._dp)) + ((-88.54241913219718_dp) + z**2*(-76.77727698139373_dp) &
                        + z**3*(3.1174289245982694_dp) + z*(126.20226718899264_dp))*Log(z*(-1._dp) + (1._dp))**2

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-62.104684476395086_dp) + Nf*(7.760195141578497_dp)
                        else
                            Ibar_select = (-62.104684476395086_dp) + Nf*(7.760195141578497_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = Nf*(12.444444444444445_dp) + (33.585006262884406_dp)
                        else
                            Ibar_select = Nf*(12.444444444444445_dp) + (33.585006262884406_dp)
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
                        Ibar_select = z**4*(-542.6627860297069_dp) + (-279.7410074813626_dp) + z*(-132.05787602026794_dp) + &
                        z**6*(-9.846043165467625_dp) + Nf*(z**2*(-37.77777777777778_dp) + (-1.3333333333333357_dp) + &
                        z*(-1.3333333333333357_dp) + (27.11111111111111_dp)/z) + (59.21762640653615_dp)/z + &
                        z**5*(105.14900153609831_dp) + z**3*(368.39046988749175_dp) + z**2*(521.8977288264662_dp) + (z*(-132._dp) + &
                        z**2*(2.2320887991927347_dp) + (12._dp) + Nf*((24._dp) + z*(24._dp)) + z**3*(458.7309417040359_dp))*Log(z) + &
                        (z*(-144._dp) + (-72._dp) + z**2*(0.9825783972125436_dp) + Nf*((5.333333333333333_dp) + &
                        z*(5.333333333333333_dp)) + z**3*(103.00736497545009_dp))*Log(z)**2 + (z**2*(-457.9088724083798_dp) + &
                        (-398.8100218336671_dp) + z**3*(65.54942528735633_dp) + z*(647.1694689546906_dp))*Log(z*(-1._dp) + (1._dp))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-225.82214554123618_dp) + Nf*(10.666666666666668_dp)
                        else
                            Ibar_select = (-225.82214554123618_dp) + Nf*(10.666666666666668_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (-90.34712078039158_dp) + Nf*(13.333333333333334_dp)
                        else
                            Ibar_select = (-90.34712078039158_dp) + Nf*(13.333333333333334_dp)
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
                        Ibar_select = z**2*(-297.33384972889235_dp) + (-198._dp)/z + (-60.00000633565014_dp) + &
                        z**5*(-3.8417329796640143_dp) + z**6*(-1.7164912280701754_dp) + z**4*(0.8431876606683805_dp) + &
                        Nf*(z*(-6.666666666666666_dp) + (-0.44444444444444464_dp)/z + z**2*(0.44444444444444464_dp) + &
                        (10.666666666666666_dp)) + z**3*(261.89724310776944_dp) + z*(304.15165631469984_dp) + (z*(-288._dp) + &
                        (-72._dp) + (-72._dp)/z + z**3*(-44.49895031490553_dp) + z**2*(1.47918188458729_dp) + &
                        Nf*((5.333333333333333_dp) + z*(5.333333333333333_dp)))*Log(z) + (z*(-572.6610814223588_dp) + &
                        z**3*(-96.38817480719794_dp) + (238.1364533075805_dp) + z**2*(286.9128029219764_dp))*Log(z*(-1._dp) + &
                        (1._dp))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-118.43525281307228_dp)
                        else
                            Ibar_select = (-118.43525281307228_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = Nf*(-4._dp) + (66._dp)
                        else
                            Ibar_select = Nf*(-4._dp) + (66._dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (144._dp)
                        else
                            Ibar_select = (144._dp)
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
                        Ibar_select = (z**5*(-1.0589078898860335e6_dp) + z**7*(-6151.961475586699_dp) + &
                        Nf**2*(z*(-80.3831877475492_dp) + (-60.57489132361303_dp) + z**5*(-4.13448275862069_dp) + &
                        z**7*(-0.035789473684210524_dp) + z**6*(0.09665427509293681_dp) + z**2*(3.4523692206836216_dp) + &
                        z**4*(21.63037634408602_dp) + z**3*(148.98202614379085_dp)) + z*(10057.663792648746_dp) + &
                        (23847.83675357499_dp) + Nf*(z**5*(-64513.179313929315_dp) + z*(-5106.295498828374_dp) + &
                        z**7*(-65.33871478315923_dp) + (664.5576023163359_dp) + z**6*(1878.6213307240705_dp) + &
                        z**2*(17389.96411960326_dp) + z**3*(22640.672143402884_dp) + z**4*(25910.682542872746_dp)) + &
                        z**2*(40658.11407642454_dp) + z**6*(65967.73482705458_dp) + z**3*(448303.280245973_dp) + &
                        z**4*(483762.28296938713_dp))/z + (Nf**2*(z*(-155.48293848852651_dp) + (-102.72147365317323_dp) + &
                        z**2*(-20.999877134783144_dp) + z**3*(-9.745923913043478_dp)) + (7128.586925899884_dp)/z + &
                        z*(15184.772999762423_dp) + (46139.85372551423_dp) + Nf*(z**2*(-27711.77269200931_dp) + &
                        (-3956.9032392637387_dp) + z*(162.82994064299692_dp) + (176.95651548462848_dp)/z + &
                        z**3*(136126.95015563903_dp)) + z**2*(618618.8865110949_dp) + z**3*(912867.8787722085_dp))*Log(z) + &
                        (Nf*(-1343.4788758922616_dp) + Nf**2*(-19.224205853928307_dp) + (2077.154328659779_dp)/z + &
                        (2474.086418387495_dp) + z*(Nf**2*(-25.3846996810888_dp) + Nf*(899.1983710183482_dp) + &
                        (7286.177575442274_dp)) + z**2*(Nf*(-13513.01415571284_dp) + Nf**2*(-2.153256704980843_dp) + &
                        (189910.38647334423_dp)) + z**3*(Nf*(-7284.343866171003_dp) + Nf**2*(4.561757719714964_dp) + &
                        (544446.2279675674_dp)))*Log(z)**2 + (z*(-1451.7301965516265_dp) + Nf**2*(z*(-1.6131687242798354_dp) + &
                        (-1.2181069958847737_dp) + z**2*(-0.6197339246119734_dp) + z**3*(-0.3809897879025923_dp)) + &
                        (3114.4958099516366_dp) + Nf*(z**2*(-2113.2833545108006_dp) + (-259.2747952996066_dp) + &
                        z*(372.8898137538914_dp) + z**3*(19161.809097688292_dp)) + z**2*(27151.053881794087_dp) + &
                        z**3*(64215.70796935953_dp))*Log(z)**3 + (z*(-496._dp) + (-62._dp) + Nf**2*(z**2*(-0.0027378507871321013_dp) &
                        + (0.09876543209876543_dp) + z*(0.09876543209876543_dp) + z**3*(0.1307901907356948_dp)) + &
                        Nf*(z**2*(-150.1987534626039_dp) + (-19.358024691358025_dp) + z**3*(47.58798999165972_dp) + &
                        z*(50.23456790123457_dp)) + z**2*(1945.8287910552062_dp) + z**3*(24477.33321799308_dp))*Log(z)**4 + &
                        (z*(-100.8_dp) + (28.8_dp) + z**2*(57.284084084084085_dp) + z**3*(378.8297520661157_dp) + &
                        Nf*(z**2*(-4.214463840399002_dp) + (-3.7925925925925927_dp) + z*(5.674074074074074_dp) + &
                        z**3*(426.68602362204723_dp)))*Log(z)**5 + (z*(-178722.62492810359_dp) + z**3*(-51823.91534627787_dp) + &
                        Nf**2*((-19.74598273265214_dp) + z**2*(-7.258346110868289_dp) + z**3*(-0.4094931617055511_dp) + &
                        z*(33.6360442274482_dp)) + Nf*(z*(-31086.70883343061_dp) + z**3*(-6458.65815324165_dp) + &
                        (12510.219198512514_dp) + z**2*(25001.132899215398_dp)) + (64571.574553930404_dp) + &
                        z**2*(167087.46998763797_dp))*Log(z*(-1._dp) + (1._dp)) + (z**2*(-3027.7046120765613_dp) + &
                        (-1875.029130251597_dp) + Nf**2*((-0.8769988658029368_dp) + z**3*(-0.2931297709923664_dp) + &
                        z*(-0.03295919988031981_dp) + z**2*(2.091976725564512_dp)) + z**3*(953.6201145059462_dp) + &
                        z*(3286.942771107182_dp) + Nf*(z**2*(-8822.627249042573_dp) + (-2216.2369258079734_dp) + &
                        z**3*(3263.8719923002886_dp) + z*(7839.528453618015_dp)))*Log(z*(-1._dp) + (1._dp))**2 + &
                        (z*(-7295.822501797675_dp) + z**3*(-2652.9462169553326_dp) + Nf**2*(z*(-1.513763158838952_dp) + &
                        z**3*(-0.4380911435941531_dp) + (0.6366014397211648_dp) + z**2*(1.3152528627119402_dp)) + &
                        Nf*(z*(-681.2773984447608_dp) + z**3*(-181.90016638935109_dp) + (258.6268876326431_dp) + &
                        z**2*(615.2173438681355_dp)) + (2045.438142421171_dp) + z**2*(7727.330576331837_dp))*Log(z*(-1._dp) + &
                        (1._dp))**3 + Nf*(z**2*(-249.20540500909814_dp) + (-82.5140388720342_dp) + z**3*(83.59259664877888_dp) + &
                        z*(248.12684723235347_dp))*Log(z*(-1._dp) + (1._dp))**4

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-4961.132012905556_dp) + Nf**2*(-19.287855837177126_dp) + Nf*(669.0197491928262_dp)
                        else
                            Ibar_select = (-4961.132012905556_dp) + Nf**2*(-19.287855837177126_dp) + Nf*(669.0197491928262_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = Nf**2*(-20.45980038267353_dp) + (149.57067875345635_dp) + Nf*(285.6679009121497_dp)
                        else
                            Ibar_select = Nf**2*(-20.45980038267353_dp) + (149.57067875345635_dp) + Nf*(285.6679009121497_dp)
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
                        Ibar_select = (z**3*(-1.2554583065334738e6_dp) + z**5*(-853504.1497970934_dp) + (-19466.65114165154_dp) + &
                        z**7*(-6811.911013116697_dp) + Nf**2*(z*(-77.24567376060064_dp) + (-28.691358024691365_dp) + &
                        z**4*(-16.709201558011067_dp) + z**6*(-1.4100801832760597_dp) + z**7*(0.4088717454194793_dp) + &
                        z**5*(6.819786096256685_dp) + z**2*(11.421026980976881_dp) + z**3*(114.59181113469499_dp)) + &
                        z*(45952.07483227205_dp) + z**6*(61335.415549784004_dp) + z**2*(63663.14091758631_dp) + &
                        Nf*(z**3*(-171288.47162680948_dp) + z**5*(-15781.148086720303_dp) + z*(-3726.981267634797_dp) + &
                        z**7*(-13.809068153373534_dp) + z**6*(581.539988291572_dp) + (925.3388560955036_dp) + &
                        z**2*(2807.1680495258897_dp) + z**4*(185784.84368158728_dp)) + z**4*(1.965784058956797e6_dp))/z + &
                        (z**2*(-317560.83570642135_dp) + (-3873.233164028004_dp)/z + Nf**2*(z*(-104.95017339517477_dp) + &
                        (-63.83906228406365_dp) + z**3*(-3.3990107234048486_dp) + z**2*(-2.4098786828422876_dp)) + &
                        Nf*(z**3*(-92485.46657745869_dp) + z**2*(-65679.58010755393_dp) + (-2676.6286425465564_dp) + &
                        (144.06448189823294_dp)/z + z*(223.91888744485428_dp)) + (1128.0389382772778_dp) + &
                        z**3*(2506.2173325119547_dp) + z*(13451.496295308121_dp))*Log(z) + (z**2*(-37112.87219220359_dp) + &
                        z*(-2930.438763701927_dp) + Nf**2*(z*(-16.74074074074074_dp) + (-14.962962962962962_dp) + &
                        z**2*(-4.622988505747126_dp) + z**3*(-1.7004383828543594_dp)) + (11101.97485970982_dp) + &
                        Nf*(z**2*(-9956.312235719159_dp) + (-791.9205325032043_dp) + z*(801.1813869765347_dp) + &
                        z**3*(38165.12544399298_dp)) + z**3*(553273.1420413478_dp))*Log(z)**2 + (z**2*(-1772.9590905441432_dp) + &
                        z*(-1656._dp) + (-288._dp) + Nf**2*((0.5925925925925926_dp) + z*(0.5925925925925926_dp)) + &
                        Nf*(z**3*(-3302.6502615818326_dp) + z**2*(-720.7533922461586_dp) + (-93.33333333333333_dp) + &
                        z*(252.2962962962963_dp)) + z**3*(30589.183357240505_dp))*Log(z)**3 + (z*(-648._dp) + &
                        z**2*(-23.18308227114716_dp) + (180._dp) + Nf*((-23.703703703703702_dp) + z**2*(-21.00261780104712_dp) + &
                        z*(36.96296296296296_dp) + z**3*(948.1076534903281_dp)) + z**3*(17054.803303303303_dp))*Log(z)**4 + &
                        (z*(-268199.10375435144_dp) + z**3*(-83612.82188766198_dp) + Nf**2*((-19.554516026088475_dp) + &
                        z**2*(-19.489470174094443_dp) + z**3*(6.523032629558541_dp) + z*(33.85428690395771_dp)) + &
                        Nf*(z**2*(-2401.1445122102523_dp) + (-538.2676222347458_dp) + z**3*(1007.9971311778443_dp) + &
                        z*(2572.08166993382_dp)) + (89533.79381256434_dp) + z**2*(257284.78349663675_dp))*Log(z*(-1._dp) + (1._dp)) &
                        + (z*(-6822.235668326744_dp) + z**3*(-1629.9576424467502_dp) + Nf*(z**2*(-98.73138900021561_dp) + &
                        (-8.247616523550668_dp) + z**3*(34.79744179388802_dp) + z*(120.18156372987825_dp)) + (1939.1390129399965_dp) &
                        + z**2*(6045.054297833498_dp))*Log(z*(-1._dp) + (1._dp))**2 + (z*(-18021.389615262786_dp) + &
                        z**3*(-5467.282717282717_dp) + Nf*(z**2*(-218.84198432530866_dp) + (-77.03825012450218_dp) + &
                        z**3*(72.38334858188472_dp) + z*(223.49688586792612_dp)) + (6601.053448990035_dp) + &
                        z**2*(17319.61888355547_dp))*Log(z*(-1._dp) + (1._dp))**3

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-3834.620605605593_dp) + Nf**2*(-11.636426390681965_dp) + Nf*(910.99791041049_dp)
                        else
                            Ibar_select = (-3834.620605605593_dp) + Nf**2*(-11.636426390681965_dp) + Nf*(910.99791041049_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (-790.7573854657392_dp) + Nf**2*(-6.518518518518519_dp) + Nf*(367.7508988355743_dp)
                        else
                            Ibar_select = (-790.7573854657392_dp) + Nf**2*(-6.518518518518519_dp) + Nf*(367.7508988355743_dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (-806.0401503092262_dp) + Nf*(-298.6666666666667_dp)
                        else
                            Ibar_select = (-806.0401503092262_dp) + Nf*(-298.6666666666667_dp)
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
                        Ibar_select = (z**5*(-65786.24114205001_dp) + z*(-3074.9274432289267_dp) + z**7*(-1226.5934317276235_dp) + &
                        Nf**2*(z*(-14.067679290179981_dp) + z**4*(-9.551818181818183_dp) + (-7.901234567901238_dp) + &
                        z**6*(-0.6386516853932585_dp) + z**7*(0.18543342269883825_dp) + z**2*(2.353625500849584_dp) + &
                        z**5*(3.0231124807395995_dp) + z**3*(26.597213267426035_dp)) + z**2*(1215.6607819610977_dp) + &
                        Nf*(z**3*(-3278.5966893101418_dp) + z*(-1257.9741077579797_dp) + z**6*(-203.39319233299258_dp) + &
                        z**7*(37.1497005316876_dp) + z**2*(632.7509870824194_dp) + (691.1209035591867_dp) + &
                        z**5*(1619.7007818535747_dp) + z**4*(1670.9808172674916_dp)) + (4795.024945085008_dp) + &
                        z**3*(5731.709702967339_dp) + z**6*(8266.872446492906_dp) + z**4*(46251.95662623072_dp))/z + &
                        (Nf**2*((-18.962962962962962_dp) + z*(-18.962962962962962_dp) + z**2*(-4.134747348721148_dp) + &
                        z**3*(-0.3748558246828143_dp)) + (355.3057584392169_dp)/z + Nf*(z**3*(-3651.3420926003323_dp) + &
                        z**2*(-1555.44970600563_dp) + (-620.0466326947806_dp) + (254.22222222222226_dp)/z + &
                        z*(498.9493723062367_dp)) + z*(1658.115168784338_dp) + z**2*(6286.435327915089_dp) + (10024.197893466688_dp) &
                        + z**3*(68519.53422823384_dp))*Log(z) + (z*(-2448._dp) + (-504._dp) + Nf**2*z**2*(z*(-0.8774002954209749_dp) &
                        + (0.05705394190871369_dp)) + Nf*(z**2*(-204.3363323505762_dp) + (-72.88888888888889_dp) + &
                        z**3*(211.59806204808456_dp) + z*(254.22222222222226_dp)) + z**2*(563.3395498741156_dp) + &
                        z**3*(13954.264456791992_dp))*Log(z)**2 + (z*(-1152._dp) + z**2*(14.483916083916085_dp) + &
                        Nf*(z**3*(-229.83758700696052_dp) + (-37.92592592592592_dp) + z**2*(-9.75536809815951_dp) + &
                        z*(68.74074074074073_dp)) + (288._dp) + z**3*(4431.171881349045_dp))*Log(z)**3 + (z*(-17121.036351811876_dp) &
                        + z**3*(-6506.483101015873_dp) + Nf**2*(z**2*(-7.698287110506561_dp) + (-6.096735434141918_dp) + &
                        z**3*(2.5785536159600997_dp) + z*(11.216468928688379_dp)) + Nf*(z**2*(-969.607784411576_dp) + &
                        (-933.8634129772819_dp) + z**3*(335.2055190504804_dp) + z*(1888.2656783383773_dp)) + (8223.276625398003_dp) &
                        + z**2*(12525.300411821916_dp))*Log(z*(-1._dp) + (1._dp)) + (z*(-10917.885243992832_dp) + &
                        z**3*(-1819.539033457249_dp) + (5845.173105267792_dp) + z**2*(8620.25117218229_dp))*Log(z*(-1._dp) + &
                        (1._dp))**2

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-2928.8352279638225_dp) + Nf**2*(-4.222222222222223_dp) + Nf*(366.8561173623828_dp)
                        else
                            Ibar_select = (-2928.8352279638225_dp) + Nf**2*(-4.222222222222223_dp) + Nf*(366.8561173623828_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = Nf*(-71.73920880217872_dp) + (2423.5626917307827_dp)
                        else
                            Ibar_select = Nf*(-71.73920880217872_dp) + (2423.5626917307827_dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = Nf*(-320._dp) + (2878.942415607831_dp)
                        else
                            Ibar_select = Nf*(-320._dp) + (2878.942415607831_dp)
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
                        Ibar_select = ((-2848.776966243133_dp) + z**3*(-1410.431263657225_dp) + Nf**2*(z*(-2.3703703703703707_dp) + &
                        z**3*(-0.09876543209876498_dp) + (0.09876543209876498_dp) + z**2*(1.4814814814814818_dp)) + &
                        z**6*(19.664018528516635_dp) + z**4*(37.674783792611834_dp) + z**7*(49.363449210621766_dp) + &
                        Nf*(z**4*(-207.14736476711994_dp) + z*(-172.59508244396136_dp) + z**6*(-24.759451214431316_dp) + &
                        z**7*(5.298802632909237_dp) + (88.00000000000004_dp) + z**5*(107.40131324138235_dp) + &
                        z**2*(110.48170115275484_dp) + z**3*(122.6534171464293_dp)) + z**2*(393.9763556136463_dp) + &
                        z**5*(851.3064442484624_dp) + z*(1244.0000687141019_dp))/z + (z*(-1728._dp) + (-1584._dp)/z + (-432._dp) + &
                        z**3*(-249.0346608912705_dp) + z**2*(-68.29966666210295_dp) + Nf**2*((-1.1851851851851851_dp) + &
                        z*(-1.1851851851851851_dp)) + Nf*(z**3*(-16.447871638508445_dp) + z*(-6.518518518518505_dp) + &
                        z**2*(17.393603561750474_dp) + (19.55555555555557_dp) + (28.444444444444443_dp)/z))*Log(z) + (z*(-1152._dp) &
                        + (-144._dp)/z + z**2*(-4.0784447476125525_dp) + Nf*(z**3*(-36.61378413524057_dp) + (-18.962962962962962_dp) &
                        + z**2*(2.121037463976945_dp) + z*(45.03703703703704_dp)) + z**3*(84.60616844602612_dp) + &
                        (144._dp))*Log(z)**2 + (z*(-3488.3481968637952_dp) + z**3*(-911.4251169248644_dp) + &
                        Nf*(z**2*(-72.12186081941573_dp) + (-50.249386470551116_dp) + z**3*(23.97327421146933_dp) + &
                        z*(98.3979730784975_dp)) + z**2*(1167.311773819195_dp) + (2368.461539969464_dp))*Log(z*(-1._dp) + (1._dp)) + &
                        (z**2*(-3241.003370727616_dp) + (-1945.1212990144686_dp) + z**3*(1079.9410358565738_dp) + &
                        z*(4970.183633885511_dp))*Log(z*(-1._dp) + (1._dp))**2

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-692.3847762199265_dp)
                        else
                            Ibar_select = (-692.3847762199265_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = Nf*(-29.333333333333343_dp) + Nf**2*(0.8888888888888888_dp) + (1663.2230337568676_dp)
                        else
                            Ibar_select = Nf*(-29.333333333333343_dp) + Nf**2*(0.8888888888888888_dp) + (1663.2230337568676_dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (0._dp)
                        else
                            Ibar_select = (0._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-864._dp)
                        else
                            Ibar_select = (-864._dp)
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
                        Ibar_select = (z*(-4.525852287590247e6_dp) + z**2*(-1.556739559015692e6_dp) + &
                        z**3*(-1.4757454257594738e6_dp) + nf*z*(-121900.85897687894_dp) + nf*z**4*(-91606.89515101635_dp) + &
                        nf*(-43143.044043592265_dp) + nf**2*z**2*(-5459.766737435285_dp) + nf**2*z**3*(-4632.214358431698_dp) + &
                        nf**2*z*(-1313.1508848269295_dp) + nf**3*z*(-92.37161831927756_dp) + nf**3*z**5*(-30.29659222214686_dp) + &
                        nf**3*(-28.31275720164609_dp) + nf**3*z**2*(-4.453922846026542_dp) + nf**3*z**4*(32.766680047672615_dp) + &
                        nf**3*z**3*(122.66821054142441_dp) + nf**2*(2023.2792405646167_dp) + nf**2*z**4*(3447.2856719670262_dp) + &
                        nf**2*z**5*(5742.237963656736_dp) + nf*z**5*(23020.337174524568_dp) + nf*z**3*(113238.5238815937_dp) + &
                        nf*z**2*(123878.84750828883_dp) + z**4*(487529.52754728484_dp) + (1.0641792490647945e6_dp) + &
                        z**5*(6.239852302666566e6_dp) + ((-689399.3508877959_dp) + z*(-536989.479184669_dp) + &
                        z**3*(-51333.82157901244_dp) + nf**2*(z**2*(-421.2655837370321_dp) + z*(-187.64433058861232_dp) + &
                        z**5*(-180.98336300927645_dp) + (180.98336300927645_dp) + z**4*(187.644330588613_dp) + &
                        z**3*(421.2655837370319_dp)) + z**2*(51333.821579012154_dp) + nf*(z**5*(-60398.705349934724_dp) + &
                        z**2*(-31984.772787023947_dp) + z*(-7916.6558252253135_dp) + z**4*(3210.4385364885748_dp) + &
                        z**3*(41397.2073644974_dp) + (55692.48806119799_dp)) + z**4*(536989.4791846691_dp) + &
                        z**5*(689399.3508877959_dp))*Hr1(-1) + z**5*(-6.529000847165631e6_dp)*Hr1(0) + &
                        z**2*(-243227.22267451117_dp)*Hr1(0) + nf*z*(-162057.24844932597_dp)*Hr1(0) + &
                        z**4*(-118186.96980752655_dp)*Hr1(0) + nf*z**4*(-27572.613978757177_dp)*Hr1(0) + &
                        nf*(-12910.818185117048_dp)*Hr1(0) + nf**2*z**4*(-4844.672531834072_dp)*Hr1(0) + &
                        nf**2*z**3*(-2676.3623859359514_dp)*Hr1(0) + nf**3*z**2*(-106.39079853913235_dp)*Hr1(0) + &
                        nf**3*z*(-73.37845286011999_dp)*Hr1(0) + nf**3*z**5*(22.814814814814813_dp)*Hr1(0) + &
                        nf**3*z**3*(50.56363804530518_dp)*Hr1(0) + nf**3*z**4*(106.39079853913235_dp)*Hr1(0) + &
                        nf**2*(567.4041113735102_dp)*Hr1(0) + nf**2*z**5*(756.6401435086418_dp)*Hr1(0) + &
                        nf**2*z*(2118.2921862466264_dp)*Hr1(0) + nf**2*z**2*(4528.476254419023_dp)*Hr1(0) + &
                        nf*z**2*(29318.81261310885_dp)*Hr1(0) + nf*z**5*(74136.2000831097_dp)*Hr1(0) + &
                        nf*z**3*(101251.23590934378_dp)*Hr1(0) + z**3*(228788.59921246555_dp)*Hr1(0) + (302013.2661103723_dp)*Hr1(0) &
                        + z*(6.179720492284087e6_dp)*Hr1(0) + z**3*(-549882.7359828057_dp)*Hr1(1) + &
                        z**4*(-400313.08806602424_dp)*Hr1(1) + z*(-205566.5672068271_dp)*Hr1(1) + &
                        z**2*(-110750.2151236086_dp)*Hr1(1) + nf*z**4*(-23109.03964717267_dp)*Hr1(1) + &
                        nf*z**5*(-11812.31981875029_dp)*Hr1(1) + nf*z*(-9225.462369810091_dp)*Hr1(1) + &
                        nf*(-9073.416031911664_dp)*Hr1(1) + nf**2*z**2*(-1862.182006496999_dp)*Hr1(1) + &
                        nf**2*z**3*(-1297.540031188357_dp)*Hr1(1) + nf**3*z**3*(-51.18518518518518_dp)*Hr1(1) + &
                        nf**3*z**2*(-42.592592592592595_dp)*Hr1(1) + nf**3*(18.074074074074073_dp)*Hr1(1) + &
                        nf**3*z**5*(22.814814814814813_dp)*Hr1(1) + nf**3*z**4*(24.518518518518515_dp)*Hr1(1) + &
                        nf**3*z*(28.370370370370374_dp)*Hr1(1) + nf**2*z*(333.46030021062506_dp)*Hr1(1) + &
                        nf**2*z**4*(840.670176753835_dp)*Hr1(1) + nf**2*z**5*(964.0797309777321_dp)*Hr1(1) + &
                        nf**2*(1021.5118297431641_dp)*Hr1(1) + nf*z**3*(30038.784279419415_dp)*Hr1(1) + &
                        nf*z**2*(41183.45776994337_dp)*Hr1(1) + (153545.82296790314_dp)*Hr1(1) + z**5*(397931.82296790305_dp)*Hr1(1) &
                        + z**4*(-1.5448108242020824e6_dp)*Hr2(-1,-1) + z**5*(-1.3287849230710384e6_dp)*Hr2(-1,-1) + &
                        z**3*(-856938.8758637585_dp)*Hr2(-1,-1) + nf*(-20230.49577681072_dp)*Hr2(-1,-1) + &
                        nf*z**3*(-8860.711506755779_dp)*Hr2(-1,-1) + nf*z*(-473.7410112522848_dp)*Hr2(-1,-1) + &
                        nf*z**4*(473.74101125228844_dp)*Hr2(-1,-1) + nf*z**2*(8860.71150675578_dp)*Hr2(-1,-1) + &
                        nf*z**5*(20230.49577681072_dp)*Hr2(-1,-1) + z**2*(856938.8758637585_dp)*Hr2(-1,-1) + &
                        (1.3287849230710384e6_dp)*Hr2(-1,-1) + z*(1.5448108242020824e6_dp)*Hr2(-1,-1) + &
                        z*(-2.1198709813133436e6_dp)*Hr2(-1,0) + (-2.019221800410253e6_dp)*Hr2(-1,0) + &
                        z**2*(-457484.51217718376_dp)*Hr2(-1,0) + nf*z*(-53628.58799193395_dp)*Hr2(-1,0) + &
                        nf*z**2*(-48672.96230373865_dp)*Hr2(-1,0) + nf*(-11770.90476960474_dp)*Hr2(-1,0) + &
                        nf**2*z**5*(-647.4403292181071_dp)*Hr2(-1,0) + nf**2*z**2*(-468.9711934156376_dp)*Hr2(-1,0) + &
                        nf**2*z**4*(-29.135802469135797_dp)*Hr2(-1,0) + nf**2*z*(-13.530864197530663_dp)*Hr2(-1,0) + &
                        nf**2*z**3*(554.304526748971_dp)*Hr2(-1,0) + nf**2*(604.7736625514403_dp)*Hr2(-1,0) + &
                        nf*z**5*(22921.52205355536_dp)*Hr2(-1,0) + nf*z**3*(26371.72773583735_dp)*Hr2(-1,0) + &
                        nf*z**4*(64779.20527588459_dp)*Hr2(-1,0) + z**3*(457484.5121771839_dp)*Hr2(-1,0) + &
                        z**5*(2.019221800410253e6_dp)*Hr2(-1,0) + z**4*(2.119870981313344e6_dp)*Hr2(-1,0) + &
                        z*(-1.3737729793130937e6_dp)*Hr2(0,-1) + z**2*(-1.145610845329124e6_dp)*Hr2(0,-1) + &
                        (-863335.4206752796_dp)*Hr2(0,-1) + nf*z**4*(-31811.76422478605_dp)*Hr2(0,-1) + &
                        nf*z*(-29074.59393755061_dp)*Hr2(0,-1) + nf*z**5*(-21166.28049039549_dp)*Hr2(0,-1) + &
                        nf**2*z*(-332.3985284712564_dp)*Hr2(0,-1) + nf**2*z**4*(-332.3985284712564_dp)*Hr2(0,-1) + &
                        nf**2*z**2*(332.3985284712564_dp)*Hr2(0,-1) + nf**2*z**3*(332.3985284712564_dp)*Hr2(0,-1) + &
                        nf*(14838.036364778487_dp)*Hr2(0,-1) + nf*z**2*(27869.77111881022_dp)*Hr2(0,-1) + &
                        nf*z**3*(42187.27723665717_dp)*Hr2(0,-1) + z**3*(60562.89612174796_dp)*Hr2(0,-1) + &
                        z**4*(1.4895468194982861e6_dp)*Hr2(0,-1) + z**5*(1.9601901331909082e6_dp)*Hr2(0,-1) + &
                        z**5*(-354168.7026133328_dp)*Hr2(0,0) + z*(-347133.4485395173_dp)*Hr2(0,0) + &
                        z**4*(-280659.92759531364_dp)*Hr2(0,0) + nf*z**5*(-117740.9105944414_dp)*Hr2(0,0) + &
                        nf*z**2*(-55530.764553861176_dp)*Hr2(0,0) + nf**2*z**4*(-4670.547791564198_dp)*Hr2(0,0) + &
                        nf**2*z**3*(-1927.5316138177227_dp)*Hr2(0,0) + nf*(-717.2182272232396_dp)*Hr2(0,0) + &
                        nf**2*z**5*(-43.65432098765427_dp)*Hr2(0,0) + nf**3*z*(-12.641975308641975_dp)*Hr2(0,0) + &
                        nf**3*z**2*(-12.641975308641975_dp)*Hr2(0,0) + nf**3*z**5*(3.1604938271604937_dp)*Hr2(0,0) + &
                        nf**3*z**3*(9.481481481481481_dp)*Hr2(0,0) + nf**3*z**4*(12.641975308641975_dp)*Hr2(0,0) + &
                        nf**2*z*(2077.8526014720437_dp)*Hr2(0,0) + nf**2*z**2*(4670.547791564198_dp)*Hr2(0,0) + &
                        nf*z*(25905.29167461316_dp)*Hr2(0,0) + (39294.48262686314_dp)*Hr2(0,0) + &
                        nf*z**4*(54116.50075857983_dp)*Hr2(0,0) + nf*z**3*(93091.99307736027_dp)*Hr2(0,0) + &
                        z**2*(227403.19445194132_dp)*Hr2(0,0) + z**3*(455188.6249882117_dp)*Hr2(0,0) + &
                        z**5*(-1.1411721451709713e6_dp)*Hr2(0,1) + (-180197.59854229694_dp)*Hr2(0,1) + &
                        z**2*(-71331.80621106409_dp)*Hr2(0,1) + nf*z**2*(-42273.91194552105_dp)*Hr2(0,1) + &
                        nf*z**5*(-40680.55779356958_dp)*Hr2(0,1) + nf*(-16341.732292926388_dp)*Hr2(0,1) + &
                        nf**2*z**4*(-1105.3305833120646_dp)*Hr2(0,1) + nf**2*z**3*(-983.5528055342867_dp)*Hr2(0,1) + &
                        nf**2*z**5*(-296.09876543209873_dp)*Hr2(0,1) + nf**3*z**4*(-48.592592592592595_dp)*Hr2(0,1) + &
                        nf**3*z**3*(-39.111111111111114_dp)*Hr2(0,1) + nf**3*z**5*(-9.481481481481481_dp)*Hr2(0,1) + &
                        nf**3*z*(48.592592592592595_dp)*Hr2(0,1) + nf**3*z**2*(48.592592592592595_dp)*Hr2(0,1) + &
                        nf**2*(317.03703703703707_dp)*Hr2(0,1) + nf**2*z**2*(894.9602129416946_dp)*Hr2(0,1) + &
                        nf**2*z*(1386.3182376330524_dp)*Hr2(0,1) + nf*z**3*(12139.06085645721_dp)*Hr2(0,1) + &
                        nf*z*(20426.01491460779_dp)*Hr2(0,1) + nf*z**4*(49552.680193438275_dp)*Hr2(0,1) + &
                        z**4*(194248.15340962174_dp)*Hr2(0,1) + z**3*(312031.5063381115_dp)*Hr2(0,1) + &
                        z*(756225.9341177951_dp)*Hr2(0,1) + (-1.4308073374441604e6_dp)*Hr2(1,0) + &
                        z**5*(-1.4217353374441606e6_dp)*Hr2(1,0) + z**3*(-127783.44753633116_dp)*Hr2(1,0) + &
                        z**2*(-109315.44753633105_dp)*Hr2(1,0) + nf*z**2*(-24893.9805132652_dp)*Hr2(1,0) + &
                        nf*z**3*(-3550.128661413333_dp)*Hr2(1,0) + nf*z**5*(-696.0781435291537_dp)*Hr2(1,0) + &
                        nf**2*z**5*(-676.5432098765432_dp)*Hr2(1,0) + nf**2*(-614.9135802469136_dp)*Hr2(1,0) + &
                        nf**2*z**2*(-29.234567901234545_dp)*Hr2(1,0) + nf**3*z**2*(-5.530864197530864_dp)*Hr2(1,0) + &
                        nf**3*z**3*(-5.530864197530864_dp)*Hr2(1,0) + nf**3*z*(2.3703703703703702_dp)*Hr2(1,0) + &
                        nf**3*z**4*(2.3703703703703702_dp)*Hr2(1,0) + nf**3*(3.1604938271604937_dp)*Hr2(1,0) + &
                        nf**3*z**5*(3.1604938271604937_dp)*Hr2(1,0) + nf**2*z**3*(70.02469135802474_dp)*Hr2(1,0) + &
                        nf**2*z*(713.1851851851852_dp)*Hr2(1,0) + nf**2*z**4*(750.8148148148149_dp)*Hr2(1,0) + &
                        nf*z*(1341.875906213092_dp)*Hr2(1,0) + nf*(9885.847782396773_dp)*Hr2(1,0) + &
                        nf*z**4*(12103.801832139026_dp)*Hr2(1,0) + z**4*(1.2877521419334463e6_dp)*Hr2(1,0) + &
                        z*(1.2971481419334468e6_dp)*Hr2(1,0) + ((-246732.74028856884_dp) + z**5*(-228588.74028856884_dp) + &
                        z**3*(-215373.2852250387_dp) + z**2*(-163928.56251947262_dp) + nf**3*((-9.481481481481481_dp) + &
                        z**5*(-9.481481481481481_dp) + z*(-7.11111111111111_dp) + z**4*(-7.11111111111111_dp) + &
                        z**2*(16.59259259259259_dp) + z**3*(16.59259259259259_dp)) + nf**2*(z*(-736.5925925925926_dp) + &
                        z**4*(-561.7777777777778_dp) + z**5*(-207.8024691358025_dp) + (15.012345679012347_dp) + &
                        z**2*(546.7654320987651_dp) + z**3*(944.3950617283948_dp)) + nf*(z**3*(-7494.594543846192_dp) + &
                        z**5*(-5587.570472075438_dp) + z**2*(-5398.495778414097_dp) + (-4127.619854791487_dp) + &
                        z**4*(9526.115633205582_dp) + z*(13082.16501592163_dp)) + z**4*(124014.00545299193_dp) + &
                        z*(157314.72815855814_dp))*Hr2(1,1) + (-613968.3505829668_dp)*Hr3(-1,-1,-1) + &
                        z*(-613968.3505829668_dp)*Hr3(-1,-1,-1) + z**2*(-613968.3505829668_dp)*Hr3(-1,-1,-1) + &
                        z**3*(613968.3505829668_dp)*Hr3(-1,-1,-1) + z**4*(613968.3505829668_dp)*Hr3(-1,-1,-1) + &
                        z**5*(613968.3505829668_dp)*Hr3(-1,-1,-1) + z**4*(-862891.7764312497_dp)*Hr3(-1,-1,0) + &
                        z**3*(-697939.7764312498_dp)*Hr3(-1,-1,0) + z**5*(-691963.7764312497_dp)*Hr3(-1,-1,0) + &
                        nf*z**5*(-28961.97530864198_dp)*Hr3(-1,-1,0) + nf*z**4*(-6935.407407407406_dp)*Hr3(-1,-1,0) + &
                        nf*z**2*(-6000.4938271604915_dp)*Hr3(-1,-1,0) + nf**2*z**2*(-128.1975308641977_dp)*Hr3(-1,-1,0) + &
                        nf**2*z*(-78.81481481481501_dp)*Hr3(-1,-1,0) + nf**2*(-14.617283950617209_dp)*Hr3(-1,-1,0) + &
                        nf**2*z**5*(14.617283950617207_dp)*Hr3(-1,-1,0) + nf**2*z**4*(78.8148148148149_dp)*Hr3(-1,-1,0) + &
                        nf**2*z**3*(128.1975308641977_dp)*Hr3(-1,-1,0) + nf*z*(4541.3333333333285_dp)*Hr3(-1,-1,0) + &
                        nf*z**3*(10788.641975308647_dp)*Hr3(-1,-1,0) + nf*(26567.901234567897_dp)*Hr3(-1,-1,0) + &
                        (691963.7764312497_dp)*Hr3(-1,-1,0) + z**2*(697939.7764312497_dp)*Hr3(-1,-1,0) + &
                        z*(862891.77643125_dp)*Hr3(-1,-1,0) + z**3*(-571331.6595702608_dp)*Hr3(-1,0,-1) + &
                        z**4*(-571331.6595702608_dp)*Hr3(-1,0,-1) + z**5*(-571331.6595702608_dp)*Hr3(-1,0,-1) + &
                        (571331.6595702608_dp)*Hr3(-1,0,-1) + z*(571331.6595702608_dp)*Hr3(-1,0,-1) + &
                        z**2*(571331.6595702608_dp)*Hr3(-1,0,-1) + z**5*(-464802.0825691279_dp)*Hr3(-1,0,0) + &
                        z**2*(-371405.9174308719_dp)*Hr3(-1,0,0) + z**4*(-350178.0825691279_dp)*Hr3(-1,0,0) + &
                        nf*(-41129.28395061728_dp)*Hr3(-1,0,0) + nf*z*(-29042.024691358023_dp)*Hr3(-1,0,0) + &
                        nf*z**3*(-4420.296296296298_dp)*Hr3(-1,0,0) + nf**2*z**5*(-62.419753086419746_dp)*Hr3(-1,0,0) + &
                        nf**2*z**2*(-39.012345679012356_dp)*Hr3(-1,0,0) + nf**2*z*(-8.59259259259261_dp)*Hr3(-1,0,0) + &
                        nf**2*z**4*(8.59259259259261_dp)*Hr3(-1,0,0) + nf**2*z**3*(39.012345679012356_dp)*Hr3(-1,0,0) + &
                        nf**2*(62.419753086419746_dp)*Hr3(-1,0,0) + nf*z**2*(362.2222222222233_dp)*Hr3(-1,0,0) + &
                        nf*z**4*(31071.06172839506_dp)*Hr3(-1,0,0) + nf*z**5*(43158.320987654326_dp)*Hr3(-1,0,0) + &
                        z*(350178.08256912784_dp)*Hr3(-1,0,0) + z**3*(371405.91743087204_dp)*Hr3(-1,0,0) + &
                        (464802.0825691279_dp)*Hr3(-1,0,0) + z**5*(-259465.98539237637_dp)*Hr3(-1,0,1) + &
                        z**4*(-115897.98539237637_dp)*Hr3(-1,0,1) + z**2*(-53374.01460762362_dp)*Hr3(-1,0,1) + &
                        nf*z**3*(-16221.629629629624_dp)*Hr3(-1,0,1) + nf*(-13526.518518518516_dp)*Hr3(-1,0,1) + &
                        nf*z**4*(-5127.111111111108_dp)*Hr3(-1,0,1) + nf**2*z**3*(-192._dp)*Hr3(-1,0,1) + &
                        nf**2*(-117.33333333333334_dp)*Hr3(-1,0,1) + nf**2*z**4*(-74.66666666666667_dp)*Hr3(-1,0,1) + &
                        nf**2*z*(74.66666666666667_dp)*Hr3(-1,0,1) + nf**2*z**5*(117.33333333333334_dp)*Hr3(-1,0,1) + &
                        nf**2*z**2*(192._dp)*Hr3(-1,0,1) + nf*z*(6791.11111111111_dp)*Hr3(-1,0,1) + &
                        nf*z**2*(12893.629629629628_dp)*Hr3(-1,0,1) + nf*z**5*(15190.518518518516_dp)*Hr3(-1,0,1) + &
                        z**3*(53374.014607623634_dp)*Hr3(-1,0,1) + z*(115897.98539237642_dp)*Hr3(-1,0,1) + &
                        (259465.98539237637_dp)*Hr3(-1,0,1) + z**4*(-758933.1000261671_dp)*Hr3(0,-1,-1) + &
                        z**5*(-673659.7180007552_dp)*Hr3(0,-1,-1) + z**3*(-383730.2191143544_dp)*Hr3(0,-1,-1) + &
                        nf*z**2*(-15650.999334705259_dp)*Hr3(0,-1,-1) + nf*z**3*(-15650.999334705259_dp)*Hr3(0,-1,-1) + &
                        nf*z*(15650.999334705259_dp)*Hr3(0,-1,-1) + nf*z**4*(15650.999334705259_dp)*Hr3(0,-1,-1) + &
                        z*(383730.21911435435_dp)*Hr3(0,-1,-1) + (469003.6011397663_dp)*Hr3(0,-1,-1) + &
                        z**2*(758933.100026167_dp)*Hr3(0,-1,-1) + z**2*(-817335.2022795327_dp)*Hr3(0,-1,0) + &
                        (-451394.93934230757_dp)*Hr3(0,-1,0) + z**5*(-253074.2966068684_dp)*Hr3(0,-1,0) + &
                        nf*z*(-50696.8268651662_dp)*Hr3(0,-1,0) + nf*(-16906.666666666664_dp)*Hr3(0,-1,0) + &
                        nf*z**4*(-5102.654025660049_dp)*Hr3(0,-1,0) + nf**2*z**4*(-633.7777777777777_dp)*Hr3(0,-1,0) + &
                        nf**2*z*(-618.9629629629632_dp)*Hr3(0,-1,0) + nf**2*z**5*(-7.901234567901157_dp)*Hr3(0,-1,0) + &
                        nf**2*(184.49382716049382_dp)*Hr3(0,-1,0) + nf**2*z**2*(545.2839506172838_dp)*Hr3(0,-1,0) + &
                        nf**2*z**3*(594.8641975308641_dp)*Hr3(0,-1,0) + nf*z**2*(13081.320692326732_dp)*Hr3(0,-1,0) + &
                        nf*z**3*(22854.851556524256_dp)*Hr3(0,-1,0) + nf*z**5*(34657.97530864197_dp)*Hr3(0,-1,0) + &
                        z**3*(175073.44045590656_dp)*Hr3(0,-1,0) + z*(601590.5595440934_dp)*Hr3(0,-1,0) + &
                        z**4*(933111.2022795326_dp)*Hr3(0,-1,0) + z**2*(-639550.3651905904_dp)*Hr3(0,0,-1) + &
                        z**3*(-486058.2775448489_dp)*Hr3(0,0,-1) + (-196128.7786584477_dp)*Hr3(0,0,-1) + &
                        nf*z*(-30204.401182435045_dp)*Hr3(0,0,-1) + nf*z**4*(-18019.704390966705_dp)*Hr3(0,0,-1) + &
                        nf*z**2*(18019.704390966708_dp)*Hr3(0,0,-1) + nf*z**3*(30204.401182435045_dp)*Hr3(0,0,-1) + &
                        z**5*(417839.5719245192_dp)*Hr3(0,0,-1) + z*(486058.2775448489_dp)*Hr3(0,0,-1) + &
                        z**4*(639550.3651905904_dp)*Hr3(0,0,-1) + z**3*(-1.0410754550635297e6_dp)*Hr3(0,0,0) + &
                        z*(-686564.54493647_dp)*Hr3(0,0,0) + z**4*(-139855.59098653877_dp)*Hr3(0,0,0) + &
                        nf*z**2*(-73722.53213479299_dp)*Hr3(0,0,0) + nf*z**5*(-38609.18518518518_dp)*Hr3(0,0,0) + &
                        nf*z*(-34892.184050826654_dp)*Hr3(0,0,0) + nf**2*z**4*(-2555.5555555555557_dp)*Hr3(0,0,0) + &
                        nf**2*z**3*(-636.641975308642_dp)*Hr3(0,0,0) + nf**2*z**5*(-69.1358024691358_dp)*Hr3(0,0,0) + &
                        nf**2*z*(737.7777777777779_dp)*Hr3(0,0,0) + nf**2*z**2*(2555.5555555555557_dp)*Hr3(0,0,0) + &
                        nf*z**3*(66685.36923601184_dp)*Hr3(0,0,0) + nf*z**4*(73722.53213479299_dp)*Hr3(0,0,0) + &
                        z**2*(139855.59098653882_dp)*Hr3(0,0,0) + z**5*(1.5877805157212226e6_dp)*Hr3(0,0,0) + &
                        z*(-910306.6179745881_dp)*Hr3(0,0,1) + z**3*(-232297.38202541182_dp)*Hr3(0,0,1) + &
                        (-140482.1022533653_dp)*Hr3(0,0,1) + z**4*(-108635.04605006878_dp)*Hr3(0,0,1) + &
                        nf*z**2*(-42316.84725779807_dp)*Hr3(0,0,1) + nf*z**5*(-36149.135802469136_dp)*Hr3(0,0,1) + &
                        nf*(-8021.333333333333_dp)*Hr3(0,0,1) + nf**2*z**3*(-970.0246913580247_dp)*Hr3(0,0,1) + &
                        nf**2*z**4*(-473.3333333333333_dp)*Hr3(0,0,1) + nf**2*z**5*(-253.23456790123453_dp)*Hr3(0,0,1) + &
                        nf**3*z**3*(-4.7407407407407405_dp)*Hr3(0,0,1) + nf**3*z**4*(-4.7407407407407405_dp)*Hr3(0,0,1) + &
                        nf**3*z*(4.7407407407407405_dp)*Hr3(0,0,1) + nf**3*z**2*(4.7407407407407405_dp)*Hr3(0,0,1) + &
                        nf**2*z**2*(489.33333333333326_dp)*Hr3(0,0,1) + nf**2*z*(1239.2592592592594_dp)*Hr3(0,0,1) + &
                        nf*z*(8410.398800143437_dp)*Hr3(0,0,1) + nf*z**3*(15690.737002325695_dp)*Hr3(0,0,1) + &
                        nf*z**4*(44050.1805911314_dp)*Hr3(0,0,1) + z**2*(190175.04605006895_dp)*Hr3(0,0,1) + &
                        z**5*(1.207965235949176e6_dp)*Hr3(0,0,1) + z*(-547594.7056203296_dp)*Hr3(0,1,0) + &
                        (-274270.79326607136_dp)*Hr3(0,1,0) + z**2*(-241937.61574738997_dp)*Hr3(0,1,0) + &
                        nf*z**5*(-25071.209876543213_dp)*Hr3(0,1,0) + nf*z**3*(-15009.175824418277_dp)*Hr3(0,1,0) + &
                        nf*z**2*(-7170.182200273079_dp)*Hr3(0,1,0) + nf*(-5682.962962962965_dp)*Hr3(0,1,0) + &
                        nf**2*(-291.55555555555554_dp)*Hr3(0,1,0) + nf**2*z*(-207.4074074074074_dp)*Hr3(0,1,0) + &
                        nf**2*z**5*(-104.29629629629626_dp)*Hr3(0,1,0) + nf**3*z**3*(-4.7407407407407405_dp)*Hr3(0,1,0) + &
                        nf**3*z**4*(-4.7407407407407405_dp)*Hr3(0,1,0) + nf**3*z*(4.7407407407407405_dp)*Hr3(0,1,0) + &
                        nf**3*z**2*(4.7407407407407405_dp)*Hr3(0,1,0) + nf**2*z**4*(41.48148148148145_dp)*Hr3(0,1,0) + &
                        nf**2*z**2*(186.07407407407405_dp)*Hr3(0,1,0) + nf**2*z**3*(311.70370370370364_dp)*Hr3(0,1,0) + &
                        nf*z**4*(3445.1451632360445_dp)*Hr3(0,1,0) + nf*z*(28560.385700961488_dp)*Hr3(0,1,0) + &
                        z**3*(198178.7056203296_dp)*Hr3(0,1,0) + z**5*(376311.88313901104_dp)*Hr3(0,1,0) + &
                        z**4*(508625.6157473899_dp)*Hr3(0,1,0) + z**5*(-87502.01460762361_dp)*Hr3(0,1,1) + &
                        nf*z**5*(-12326.518518518516_dp)*Hr3(0,1,1) + nf*(-11205.333333333332_dp)*Hr3(0,1,1) + &
                        nf*z**3*(-6465.477930371273_dp)*Hr3(0,1,1) + nf*z**2*(-1724.9418227151455_dp)*Hr3(0,1,1) + &
                        (-1352.7202279532648_dp)*Hr3(0,1,1) + nf**2*z*(-456.29629629629613_dp)*Hr3(0,1,1) + &
                        nf**2*(-174.22222222222223_dp)*Hr3(0,1,1) + nf**2*z**5*(-174.22222222222223_dp)*Hr3(0,1,1) + &
                        nf**2*z**2*(-164.74074074074073_dp)*Hr3(0,1,1) + nf**3*z*(-14.22222222222222_dp)*Hr3(0,1,1) + &
                        nf**3*z**2*(-14.22222222222222_dp)*Hr3(0,1,1) + nf**3*z**3*(14.22222222222222_dp)*Hr3(0,1,1) + &
                        nf**3*z**4*(14.22222222222222_dp)*Hr3(0,1,1) + nf**2*z**4*(306.962962962963_dp)*Hr3(0,1,1) + &
                        nf**2*z**3*(598.5185185185185_dp)*Hr3(0,1,1) + nf*z**4*(2466.2751560484785_dp)*Hr3(0,1,1) + &
                        nf*z*(8327.99644888979_dp)*Hr3(0,1,1) + z**2*(10996.574151717097_dp)*Hr3(0,1,1) + &
                        z*(70547.30898729397_dp)*Hr3(0,1,1) + z**4*(119395.4258482829_dp)*Hr3(0,1,1) + &
                        z**3*(214212.69101270603_dp)*Hr3(0,1,1) + z**4*(-119382.73483557686_dp)*Hr3(1,0,-1) + &
                        z*(-119382.73483557685_dp)*Hr3(1,0,-1) + (119382.73483557686_dp)*Hr3(1,0,-1) + &
                        z**2*(119382.73483557686_dp)*Hr3(1,0,-1) + z**3*(119382.73483557686_dp)*Hr3(1,0,-1) + &
                        z**5*(119382.73483557686_dp)*Hr3(1,0,-1) + z**4*(-201292.18482249323_dp)*Hr3(1,0,0) + &
                        z*(-200860.18482249326_dp)*Hr3(1,0,0) + z**3*(-90055.81517750677_dp)*Hr3(1,0,0) + &
                        z**2*(-89623.81517750675_dp)*Hr3(1,0,0) + nf*z**5*(-30171.456790123455_dp)*Hr3(1,0,0) + &
                        nf*(-28260.938271604937_dp)*Hr3(1,0,0) + nf*z**2*(-12614.938271604946_dp)*Hr3(1,0,0) + &
                        nf*z**3*(-8649.901234567911_dp)*Hr3(1,0,0) + nf**2*z*(-60.888888888888914_dp)*Hr3(1,0,0) + &
                        nf**2*z**4*(-60.888888888888914_dp)*Hr3(1,0,0) + nf**2*(36.74074074074074_dp)*Hr3(1,0,0) + &
                        nf**2*z**5*(36.74074074074074_dp)*Hr3(1,0,0) + nf**2*z**2*(72.1481481481481_dp)*Hr3(1,0,0) + &
                        nf**2*z**3*(72.14814814814815_dp)*Hr3(1,0,0) + nf*z*(28597.358024691355_dp)*Hr3(1,0,0) + &
                        nf*z**4*(30651.876543209873_dp)*Hr3(1,0,0) + (217564.1848224932_dp)*Hr3(1,0,0) + &
                        z**5*(217564.1848224932_dp)*Hr3(1,0,0) + (-265971.3674177884_dp)*Hr3(1,0,1) + &
                        z**5*(-265971.3674177884_dp)*Hr3(1,0,1) + nf*z**2*(-6500.938271604939_dp)*Hr3(1,0,1) + &
                        nf*z**5*(-5148.839506172841_dp)*Hr3(1,0,1) + nf*z*(-4413.629629629635_dp)*Hr3(1,0,1) + &
                        nf*(-2565.1358024691363_dp)*Hr3(1,0,1) + nf*z**4*(-1397.9259259259234_dp)*Hr3(1,0,1) + &
                        nf*z**3*(-901.5308641975298_dp)*Hr3(1,0,1) + nf**2*(-310.12345679012344_dp)*Hr3(1,0,1) + &
                        nf**2*z**5*(-310.12345679012344_dp)*Hr3(1,0,1) + nf**2*z**4*(-176.59259259259258_dp)*Hr3(1,0,1) + &
                        nf**2*z*(-176.59259259259255_dp)*Hr3(1,0,1) + nf**2*z**3*(454.716049382716_dp)*Hr3(1,0,1) + &
                        nf**2*z**2*(454.7160493827161_dp)*Hr3(1,0,1) + z**3*(107996.63258221159_dp)*Hr3(1,0,1) + &
                        z**2*(109292.63258221158_dp)*Hr3(1,0,1) + z**4*(319827.3674177885_dp)*Hr3(1,0,1) + &
                        z*(321123.3674177885_dp)*Hr3(1,0,1) + (-419463.4550635302_dp)*Hr3(1,1,0) + &
                        z**5*(-419463.4550635302_dp)*Hr3(1,1,0) + z**3*(-45495.455063530004_dp)*Hr3(1,1,0) + &
                        z**2*(-44199.45506353019_dp)*Hr3(1,1,0) + nf*z**5*(-6713.283950617284_dp)*Hr3(1,1,0) + &
                        nf*z**2*(-6354.765432098767_dp)*Hr3(1,1,0) + nf*(-4129.580246913578_dp)*Hr3(1,1,0) + &
                        nf*z*(-2995.3580246913584_dp)*Hr3(1,1,0) + nf*z**3*(-755.3580246913584_dp)*Hr3(1,1,0) + &
                        nf**2*z**5*(-278.51851851851853_dp)*Hr3(1,1,0) + nf**2*(-278.5185185185185_dp)*Hr3(1,1,0) + &
                        nf**2*z*(-152.88888888888886_dp)*Hr3(1,1,0) + nf**2*z**4*(-152.88888888888886_dp)*Hr3(1,1,0) + &
                        nf*z**4*(20.34567901234459_dp)*Hr3(1,1,0) + nf**2*z**2*(399.40740740740745_dp)*Hr3(1,1,0) + &
                        nf**2*z**3*(399.40740740740745_dp)*Hr3(1,1,0) + z**4*(473319.45506353024_dp)*Hr3(1,1,0) + &
                        z*(474615.45506352995_dp)*Hr3(1,1,0) + z**4*(-2592._dp)*Hr3(1,1,1) + &
                        nf*z**5*(-2339.5555555555557_dp)*Hr3(1,1,1) + nf*(-1415.111111111111_dp)*Hr3(1,1,1) + &
                        nf*z**2*(-507.6543209876544_dp)*Hr3(1,1,1) + nf*z*(134.32098765432107_dp)*Hr3(1,1,1) + &
                        nf*z**4*(1922.765432098766_dp)*Hr3(1,1,1) + nf*z**3*(2205.234567901235_dp)*Hr3(1,1,1) + &
                        z**2*(2592._dp)*Hr3(1,1,1) + z**4*(-272160._dp)*Hr4(-1,-1,-1,0) + &
                        z**3*(-196128.00000000003_dp)*Hr4(-1,-1,-1,0) + z**5*(-161568._dp)*Hr4(-1,-1,-1,0) + &
                        nf*(-7715.555555555556_dp)*Hr4(-1,-1,-1,0) + nf*z*(-1685.3333333333337_dp)*Hr4(-1,-1,-1,0) + &
                        nf*z**3*(-846.2222222222199_dp)*Hr4(-1,-1,-1,0) + nf*z**2*(846.2222222222199_dp)*Hr4(-1,-1,-1,0) + &
                        nf*z**4*(1685.3333333333335_dp)*Hr4(-1,-1,-1,0) + nf*z**5*(7715.555555555556_dp)*Hr4(-1,-1,-1,0) + &
                        (161568._dp)*Hr4(-1,-1,-1,0) + z**2*(196128.00000000003_dp)*Hr4(-1,-1,-1,0) + z*(272160._dp)*Hr4(-1,-1,-1,0) &
                        + z*(-808704.0000000001_dp)*Hr4(-1,-1,0,0) + (-798336._dp)*Hr4(-1,-1,0,0) + &
                        z**2*(-124415.99999999996_dp)*Hr4(-1,-1,0,0) + nf*z**5*(-14882.765432098764_dp)*Hr4(-1,-1,0,0) + &
                        nf*z**2*(-8476.839506172844_dp)*Hr4(-1,-1,0,0) + nf*z*(-506.074074074073_dp)*Hr4(-1,-1,0,0) + &
                        nf*z**4*(506.0740740740748_dp)*Hr4(-1,-1,0,0) + nf*z**3*(8476.839506172842_dp)*Hr4(-1,-1,0,0) + &
                        nf*(14882.765432098764_dp)*Hr4(-1,-1,0,0) + z**3*(124415.99999999997_dp)*Hr4(-1,-1,0,0) + &
                        z**5*(798336._dp)*Hr4(-1,-1,0,0) + z**4*(808704.0000000001_dp)*Hr4(-1,-1,0,0) + &
                        z*(-522720.00000000006_dp)*Hr4(-1,-1,0,1) + (-446688.00000000006_dp)*Hr4(-1,-1,0,1) + &
                        z**2*(-142560._dp)*Hr4(-1,-1,0,1) + nf*z**5*(-8440.888888888889_dp)*Hr4(-1,-1,0,1) + &
                        nf*z**2*(-4963.555555555551_dp)*Hr4(-1,-1,0,1) + nf*z*(-554.6666666666679_dp)*Hr4(-1,-1,0,1) + &
                        nf*z**4*(554.6666666666661_dp)*Hr4(-1,-1,0,1) + nf*z**3*(4963.555555555557_dp)*Hr4(-1,-1,0,1) + &
                        nf*(8440.888888888889_dp)*Hr4(-1,-1,0,1) + z**3*(142560._dp)*Hr4(-1,-1,0,1) + &
                        z**5*(446688.00000000006_dp)*Hr4(-1,-1,0,1) + z**4*(522720._dp)*Hr4(-1,-1,0,1) + &
                        z*(-371520._dp)*Hr4(-1,0,-1,0) + (-361152._dp)*Hr4(-1,0,-1,0) + z**2*(-86399.99999999994_dp)*Hr4(-1,0,-1,0) &
                        + nf*z**5*(-8095.604938271604_dp)*Hr4(-1,0,-1,0) + nf*z**2*(-2295.3086419753117_dp)*Hr4(-1,0,-1,0) + &
                        nf*z**4*(-1192.2962962962956_dp)*Hr4(-1,0,-1,0) + nf*z*(1192.2962962962956_dp)*Hr4(-1,0,-1,0) + &
                        nf*z**3*(2295.30864197531_dp)*Hr4(-1,0,-1,0) + nf*(8095.604938271604_dp)*Hr4(-1,0,-1,0) + &
                        z**3*(86399.99999999997_dp)*Hr4(-1,0,-1,0) + z**5*(361152._dp)*Hr4(-1,0,-1,0) + &
                        z**4*(371520._dp)*Hr4(-1,0,-1,0) + z**5*(-639144._dp)*Hr4(-1,0,0,0) + z**4*(-606312._dp)*Hr4(-1,0,0,0) + &
                        z**3*(-26568.000000000015_dp)*Hr4(-1,0,0,0) + nf*(-9646.024691358025_dp)*Hr4(-1,0,0,0) + &
                        nf*z**3*(-7620.54320987654_dp)*Hr4(-1,0,0,0) + nf*z**4*(-1574.518518518518_dp)*Hr4(-1,0,0,0) + &
                        nf*z*(1574.5185185185146_dp)*Hr4(-1,0,0,0) + nf*z**2*(7620.543209876542_dp)*Hr4(-1,0,0,0) + &
                        nf*z**5*(9646.024691358025_dp)*Hr4(-1,0,0,0) + z**2*(26568.00000000011_dp)*Hr4(-1,0,0,0) + &
                        z*(606312._dp)*Hr4(-1,0,0,0) + (639144._dp)*Hr4(-1,0,0,0) + z**4*(-658800._dp)*Hr4(-1,0,0,1) + &
                        z**5*(-603504._dp)*Hr4(-1,0,0,1) + z**3*(-107568.00000000003_dp)*Hr4(-1,0,0,1) + &
                        nf*(-8771.555555555557_dp)*Hr4(-1,0,0,1) + nf*z**3*(-7758.222222222222_dp)*Hr4(-1,0,0,1) + &
                        nf*z**4*(-2154.6666666666633_dp)*Hr4(-1,0,0,1) + nf*z*(2154.666666666668_dp)*Hr4(-1,0,0,1) + &
                        nf*z**2*(7758.22222222222_dp)*Hr4(-1,0,0,1) + nf*z**5*(8771.555555555557_dp)*Hr4(-1,0,0,1) + &
                        z**2*(107568.00000000006_dp)*Hr4(-1,0,0,1) + (603504._dp)*Hr4(-1,0,0,1) + z*(658800._dp)*Hr4(-1,0,0,1) + &
                        z**4*(-293327.99999999994_dp)*Hr4(-1,0,1,0) + z**5*(-251856.00000000003_dp)*Hr4(-1,0,1,0) + &
                        z**3*(-65232._dp)*Hr4(-1,0,1,0) + nf*(-3047.111111111111_dp)*Hr4(-1,0,1,0) + &
                        nf*z**3*(-1852.4444444444443_dp)*Hr4(-1,0,1,0) + nf*z**4*(-245.33333333333394_dp)*Hr4(-1,0,1,0) + &
                        nf*z*(245.333333333333_dp)*Hr4(-1,0,1,0) + nf*z**2*(1852.444444444445_dp)*Hr4(-1,0,1,0) + &
                        nf*z**5*(3047.1111111111113_dp)*Hr4(-1,0,1,0) + z**2*(65231.99999999997_dp)*Hr4(-1,0,1,0) + &
                        (251856.00000000003_dp)*Hr4(-1,0,1,0) + z*(293328._dp)*Hr4(-1,0,1,0) + &
                        z**4*(-259200.00000000006_dp)*Hr4(-1,0,1,1) + z**3*(-145152.00000000003_dp)*Hr4(-1,0,1,1) + &
                        z**5*(-114048._dp)*Hr4(-1,0,1,1) + nf*z**4*(-768._dp)*Hr4(-1,0,1,1) + nf*z**3*(-512._dp)*Hr4(-1,0,1,1) + &
                        nf*z**5*(-256._dp)*Hr4(-1,0,1,1) + nf*(256._dp)*Hr4(-1,0,1,1) + nf*z**2*(512._dp)*Hr4(-1,0,1,1) + &
                        nf*z*(768._dp)*Hr4(-1,0,1,1) + (114048._dp)*Hr4(-1,0,1,1) + z**2*(145152.00000000003_dp)*Hr4(-1,0,1,1) + &
                        z*(259200.00000000006_dp)*Hr4(-1,0,1,1) + z*(-396576._dp)*Hr4(0,-1,-1,0) + z**2*(-270432._dp)*Hr4(0,-1,-1,0) &
                        + (-161568._dp)*Hr4(0,-1,-1,0) + z**3*(-97632.00000000003_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**4*(-13616._dp)*Hr4(0,-1,-1,0) + nf*z*(-10757.333333333332_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**5*(-8853.333333333334_dp)*Hr4(0,-1,-1,0) + nf**2*z*(-62.814814814814824_dp)*Hr4(0,-1,-1,0) + &
                        nf**2*z**4*(-62.814814814814824_dp)*Hr4(0,-1,-1,0) + nf**2*z**2*(62.814814814814824_dp)*Hr4(0,-1,-1,0) + &
                        nf**2*z**3*(62.814814814814824_dp)*Hr4(0,-1,-1,0) + nf*(6705.777777777779_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**2*(12094.222222222224_dp)*Hr4(0,-1,-1,0) + nf*z**3*(15578.666666666664_dp)*Hr4(0,-1,-1,0) + &
                        z**4*(346464._dp)*Hr4(0,-1,-1,0) + z**5*(560736.0000000001_dp)*Hr4(0,-1,-1,0) + &
                        z**5*(-826848.0000000001_dp)*Hr4(0,-1,0,0) + z**4*(-756864._dp)*Hr4(0,-1,0,0) + &
                        nf*z**3*(-18208.098765432103_dp)*Hr4(0,-1,0,0) + nf*(-11434.666666666666_dp)*Hr4(0,-1,0,0) + &
                        nf*z**2*(-7413.037037037039_dp)*Hr4(0,-1,0,0) + nf**2*z*(-53.925925925925924_dp)*Hr4(0,-1,0,0) + &
                        nf**2*z**4*(-53.925925925925924_dp)*Hr4(0,-1,0,0) + nf**2*z**2*(53.925925925925924_dp)*Hr4(0,-1,0,0) + &
                        nf**2*z**3*(53.925925925925924_dp)*Hr4(0,-1,0,0) + nf*z**5*(10203.654320987655_dp)*Hr4(0,-1,0,0) + &
                        nf*z*(12612.444444444449_dp)*Hr4(0,-1,0,0) + nf*z**4*(13087.703703703704_dp)*Hr4(0,-1,0,0) + &
                        z**3*(299807.99999999994_dp)*Hr4(0,-1,0,0) + z**2*(405216._dp)*Hr4(0,-1,0,0) + &
                        (446688.00000000006_dp)*Hr4(0,-1,0,0) + z*(451008._dp)*Hr4(0,-1,0,0) + z**5*(-598752._dp)*Hr4(0,-1,0,1) + &
                        z**4*(-431136._dp)*Hr4(0,-1,0,1) + nf*z**3*(-12206.222222222219_dp)*Hr4(0,-1,0,1) + &
                        nf*(-5667.555555555556_dp)*Hr4(0,-1,0,1) + nf*z**2*(-5244.444444444445_dp)*Hr4(0,-1,0,1) + &
                        nf**2*z**2*(-170.66666666666666_dp)*Hr4(0,-1,0,1) + nf**2*z**3*(-170.66666666666666_dp)*Hr4(0,-1,0,1) + &
                        nf**2*z*(170.66666666666666_dp)*Hr4(0,-1,0,1) + nf**2*z**4*(170.66666666666666_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(6645.3333333333285_dp)*Hr4(0,-1,0,1) + nf*z**4*(6880._dp)*Hr4(0,-1,0,1) + &
                        nf*z**5*(8440.888888888889_dp)*Hr4(0,-1,0,1) + z**3*(132191.99999999994_dp)*Hr4(0,-1,0,1) + &
                        (237599.99999999997_dp)*Hr4(0,-1,0,1) + z**2*(260064.00000000006_dp)*Hr4(0,-1,0,1) + &
                        z*(419040._dp)*Hr4(0,-1,0,1) + z**5*(-389664._dp)*Hr4(0,0,-1,0) + z**4*(-312768._dp)*Hr4(0,0,-1,0) + &
                        nf*z**2*(-26262.814814814814_dp)*Hr4(0,0,-1,0) + nf*z**3*(-3087.703703703702_dp)*Hr4(0,0,-1,0) + &
                        nf*(-2880._dp)*Hr4(0,0,-1,0) + nf**2*z**4*(-232.88888888888889_dp)*Hr4(0,0,-1,0) + &
                        nf**2*z*(-216.2962962962963_dp)*Hr4(0,0,-1,0) + nf**2*z**3*(216.2962962962963_dp)*Hr4(0,0,-1,0) + &
                        nf**2*z**2*(232.88888888888889_dp)*Hr4(0,0,-1,0) + nf*z**5*(1879.703703703704_dp)*Hr4(0,0,-1,0) + &
                        nf*z*(2936.0000000000005_dp)*Hr4(0,0,-1,0) + nf*z**4*(26262.814814814814_dp)*Hr4(0,0,-1,0) + &
                        (47520.00000000001_dp)*Hr4(0,0,-1,0) + z*(50112._dp)*Hr4(0,0,-1,0) + z**3*(311040._dp)*Hr4(0,0,-1,0) + &
                        z**2*(312768._dp)*Hr4(0,0,-1,0) + z**4*(-394848._dp)*Hr4(0,0,0,0) + z*(-294624._dp)*Hr4(0,0,0,0) + &
                        z**3*(-123552.00000000012_dp)*Hr4(0,0,0,0) + nf*z**2*(-54250.370370370365_dp)*Hr4(0,0,0,0) + &
                        nf*z**5*(-2268.8395061728397_dp)*Hr4(0,0,0,0) + nf**2*z**4*(-727.111111111111_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**3*(-469.33333333333337_dp)*Hr4(0,0,0,0) + nf*z**3*(-247.3086419753081_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(469.33333333333337_dp)*Hr4(0,0,0,0) + nf*z*(500.1481481481478_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**2*(727.111111111111_dp)*Hr4(0,0,0,0) + nf*z**4*(54250.370370370365_dp)*Hr4(0,0,0,0) + &
                        z**2*(394847.99999999994_dp)*Hr4(0,0,0,0) + z**5*(451440.0000000001_dp)*Hr4(0,0,0,0) + &
                        z**3*(-430271.99999999994_dp)*Hr4(0,0,0,1) + z*(-301536._dp)*Hr4(0,0,0,1) + &
                        z**4*(-257903.99999999994_dp)*Hr4(0,0,0,1) + nf*z**2*(-37680.740740740745_dp)*Hr4(0,0,0,1) + &
                        nf*z**5*(-5646.617283950617_dp)*Hr4(0,0,0,1) + nf*z**3*(-1655.753086419756_dp)*Hr4(0,0,0,1) + &
                        nf**2*z**2*(-496.2962962962963_dp)*Hr4(0,0,0,1) + nf**2*z**3*(-357.03703703703707_dp)*Hr4(0,0,0,1) + &
                        nf**2*z*(357.037037037037_dp)*Hr4(0,0,0,1) + nf**2*z**4*(496.2962962962963_dp)*Hr4(0,0,0,1) + &
                        nf*z*(3558.3703703703704_dp)*Hr4(0,0,0,1) + nf*z**4*(37104.740740740745_dp)*Hr4(0,0,0,1) + &
                        z**2*(267408._dp)*Hr4(0,0,0,1) + z**5*(793584._dp)*Hr4(0,0,0,1) + z**3*(-270000._dp)*Hr4(0,0,1,0) + &
                        z*(-262224._dp)*Hr4(0,0,1,0) + z**4*(-79056.00000000003_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(-8104.888888888889_dp)*Hr4(0,0,1,0) + nf*z*(-5316.148148148147_dp)*Hr4(0,0,1,0) + &
                        nf*z**5*(-2668.246913580247_dp)*Hr4(0,0,1,0) + nf**2*z**2*(-696.8888888888889_dp)*Hr4(0,0,1,0) + &
                        nf*(-447.99999999999994_dp)*Hr4(0,0,1,0) + nf**2*z*(-14.22222222222222_dp)*Hr4(0,0,1,0) + &
                        nf**2*z**3*(14.222222222222248_dp)*Hr4(0,0,1,0) + nf**2*z**4*(696.8888888888889_dp)*Hr4(0,0,1,0) + &
                        nf*z**3*(5680.395061728394_dp)*Hr4(0,0,1,0) + nf*z**4*(8552.888888888889_dp)*Hr4(0,0,1,0) + &
                        z**2*(31535.999999999975_dp)*Hr4(0,0,1,0) + (47520.00000000001_dp)*Hr4(0,0,1,0) + &
                        z**5*(570240._dp)*Hr4(0,0,1,0) + z*(-371087.99999999994_dp)*Hr4(0,0,1,1) + &
                        z**3*(-85104.00000000003_dp)*Hr4(0,0,1,1) + z**2*(-60912.00000000003_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(-8272.592592592591_dp)*Hr4(0,0,1,1) + nf*z*(-3496.296296296297_dp)*Hr4(0,0,1,1) + &
                        nf*(-1056._dp)*Hr4(0,0,1,1) + nf*z**5*(-1005.4320987654321_dp)*Hr4(0,0,1,1) + &
                        nf**2*z*(-649.4814814814815_dp)*Hr4(0,0,1,1) + nf**2*z**2*(-649.4814814814815_dp)*Hr4(0,0,1,1) + &
                        nf**2*z**3*(649.4814814814815_dp)*Hr4(0,0,1,1) + nf**2*z**4*(649.4814814814815_dp)*Hr4(0,0,1,1) + &
                        nf*z**3*(3061.728395061728_dp)*Hr4(0,0,1,1) + nf*z**4*(7888.592592592591_dp)*Hr4(0,0,1,1) + &
                        z**4*(22895.99999999997_dp)*Hr4(0,0,1,1) + (61776.00000000001_dp)*Hr4(0,0,1,1) + &
                        z**5*(479951.99999999994_dp)*Hr4(0,0,1,1) + z**4*(-323568._dp)*Hr4(0,1,0,0) + &
                        z*(-201744.00000000003_dp)*Hr4(0,1,0,0) + z**3*(-168912.00000000003_dp)*Hr4(0,1,0,0) + &
                        nf*z*(-6634.222222222222_dp)*Hr4(0,1,0,0) + nf*(-5070.222222222223_dp)*Hr4(0,1,0,0) + &
                        nf*z**5*(-4223.20987654321_dp)*Hr4(0,1,0,0) + nf*z**2*(-3239.259259259261_dp)*Hr4(0,1,0,0) + &
                        nf**2*z*(-19.85185185185185_dp)*Hr4(0,1,0,0) + nf**2*z**2*(-19.85185185185185_dp)*Hr4(0,1,0,0) + &
                        nf**2*z**3*(19.85185185185185_dp)*Hr4(0,1,0,0) + nf**2*z**4*(19.85185185185185_dp)*Hr4(0,1,0,0) + &
                        nf*z**4*(7157.481481481482_dp)*Hr4(0,1,0,0) + nf*z**3*(7977.432098765432_dp)*Hr4(0,1,0,0) + &
                        z**2*(114479.99999999997_dp)*Hr4(0,1,0,0) + (228096._dp)*Hr4(0,1,0,0) + z**5*(418176._dp)*Hr4(0,1,0,0) + &
                        z*(-227232._dp)*Hr4(0,1,0,1) + z**2*(-95904._dp)*Hr4(0,1,0,1) + z**3*(-10368._dp)*Hr4(0,1,0,1) + &
                        (-9504._dp)*Hr4(0,1,0,1) + nf*z**4*(-3249.7777777777783_dp)*Hr4(0,1,0,1) + &
                        nf*z**5*(-808.2962962962963_dp)*Hr4(0,1,0,1) + nf*z*(-791.7037037037035_dp)*Hr4(0,1,0,1) + &
                        nf**2*z*(-417.18518518518505_dp)*Hr4(0,1,0,1) + nf**2*z**2*(-417.18518518518505_dp)*Hr4(0,1,0,1) + &
                        nf**2*z**3*(417.18518518518505_dp)*Hr4(0,1,0,1) + nf**2*z**4*(417.18518518518505_dp)*Hr4(0,1,0,1) + &
                        nf*(1564.444444444444_dp)*Hr4(0,1,0,1) + nf*z**3*(2175.999999999999_dp)*Hr4(0,1,0,1) + &
                        nf*z**2*(3989.333333333333_dp)*Hr4(0,1,0,1) + z**4*(67392.00000000003_dp)*Hr4(0,1,0,1) + &
                        z**5*(228096._dp)*Hr4(0,1,0,1) + z*(-227232._dp)*Hr4(0,1,1,0) + z**2*(-95904._dp)*Hr4(0,1,1,0) + &
                        z**3*(-10368._dp)*Hr4(0,1,1,0) + (-9504._dp)*Hr4(0,1,1,0) + nf*z**4*(-2775.703703703704_dp)*Hr4(0,1,1,0) + &
                        nf*z*(-1573.9259259259252_dp)*Hr4(0,1,1,0) + nf*z**5*(-934.7160493827159_dp)*Hr4(0,1,1,0) + &
                        nf**2*z*(-369.77777777777777_dp)*Hr4(0,1,1,0) + nf**2*z**2*(-369.77777777777777_dp)*Hr4(0,1,1,0) + &
                        nf**2*z**3*(369.77777777777777_dp)*Hr4(0,1,1,0) + nf**2*z**4*(369.77777777777777_dp)*Hr4(0,1,1,0) + &
                        nf*(1280._dp)*Hr4(0,1,1,0) + nf*z**3*(3084.641975308642_dp)*Hr4(0,1,1,0) + &
                        nf*z**2*(3799.7037037037026_dp)*Hr4(0,1,1,0) + z**4*(67392.00000000003_dp)*Hr4(0,1,1,0) + &
                        z**5*(228096._dp)*Hr4(0,1,1,0) + (-99792._dp)*Hr4(0,1,1,1) + z**5*(-99792._dp)*Hr4(0,1,1,1) + &
                        z**2*(-32400.000000000007_dp)*Hr4(0,1,1,1) + z**3*(-32400.000000000007_dp)*Hr4(0,1,1,1) + &
                        nf*z*(-2961.777777777777_dp)*Hr4(0,1,1,1) + nf*z**4*(-1544.2962962962963_dp)*Hr4(0,1,1,1) + &
                        nf*z**2*(1686.518518518519_dp)*Hr4(0,1,1,1) + nf*(2449.777777777778_dp)*Hr4(0,1,1,1) + &
                        nf*z**3*(2531.9506172839506_dp)*Hr4(0,1,1,1) + nf*z**5*(3021.8271604938273_dp)*Hr4(0,1,1,1) + &
                        z*(89423.99999999999_dp)*Hr4(0,1,1,1) + z**4*(89424._dp)*Hr4(0,1,1,1) + (-142560._dp)*Hr4(1,0,-1,0) + &
                        z**5*(-142560._dp)*Hr4(1,0,-1,0) + nf*z**2*(-3275.061728395063_dp)*Hr4(1,0,-1,0) + &
                        nf*z**3*(-3275.061728395063_dp)*Hr4(1,0,-1,0) + nf*(966.3209876543211_dp)*Hr4(1,0,-1,0) + &
                        nf*z**5*(966.3209876543211_dp)*Hr4(1,0,-1,0) + nf*z*(1732.7407407407402_dp)*Hr4(1,0,-1,0) + &
                        nf*z**4*(1732.7407407407402_dp)*Hr4(1,0,-1,0) + z**3*(37151.999999999985_dp)*Hr4(1,0,-1,0) + &
                        z**2*(37152.000000000015_dp)*Hr4(1,0,-1,0) + z**4*(114911.99999999999_dp)*Hr4(1,0,-1,0) + &
                        z*(114912._dp)*Hr4(1,0,-1,0) + z**4*(-319464._dp)*Hr4(1,0,0,0) + z*(-319463.99999999994_dp)*Hr4(1,0,0,0) + &
                        z**3*(-22679.999999999985_dp)*Hr4(1,0,0,0) + z**2*(-22679.999999999913_dp)*Hr4(1,0,0,0) + &
                        nf*(-6795.25925925926_dp)*Hr4(1,0,0,0) + nf*z**5*(-6795.25925925926_dp)*Hr4(1,0,0,0) + &
                        nf*z**4*(355.55555555555566_dp)*Hr4(1,0,0,0) + nf*z*(355.55555555555657_dp)*Hr4(1,0,0,0) + &
                        nf*z**3*(3415.7037037037016_dp)*Hr4(1,0,0,0) + nf*z**2*(3415.7037037037035_dp)*Hr4(1,0,0,0) + &
                        (392040._dp)*Hr4(1,0,0,0) + z**5*(392040._dp)*Hr4(1,0,0,0) + z*(-185760._dp)*Hr4(1,0,0,1) + &
                        z**4*(-185760._dp)*Hr4(1,0,0,1) + z**2*(-13824._dp)*Hr4(1,0,0,1) + z**3*(-13824._dp)*Hr4(1,0,0,1) + &
                        nf*(-4853.728395061728_dp)*Hr4(1,0,0,1) + nf*z**5*(-4853.728395061728_dp)*Hr4(1,0,0,1) + &
                        nf*z*(551.7037037037035_dp)*Hr4(1,0,0,1) + nf*z**4*(551.7037037037039_dp)*Hr4(1,0,0,1) + &
                        nf*z**3*(1998.024691358024_dp)*Hr4(1,0,0,1) + nf*z**2*(1998.024691358025_dp)*Hr4(1,0,0,1) + &
                        (237599.99999999997_dp)*Hr4(1,0,0,1) + z**5*(237599.99999999997_dp)*Hr4(1,0,0,1) + &
                        z**2*(-51408.000000000015_dp)*Hr4(1,0,1,0) + z**3*(-51407.999999999985_dp)*Hr4(1,0,1,0) + &
                        z*(-43632._dp)*Hr4(1,0,1,0) + z**4*(-43631.999999999985_dp)*Hr4(1,0,1,0) + &
                        nf*z*(-2262.518518518519_dp)*Hr4(1,0,1,0) + nf*z**4*(-2262.5185185185187_dp)*Hr4(1,0,1,0) + &
                        nf*(343.30864197530855_dp)*Hr4(1,0,1,0) + nf*z**5*(343.30864197530855_dp)*Hr4(1,0,1,0) + &
                        nf*z**2*(3359.20987654321_dp)*Hr4(1,0,1,0) + nf*z**3*(3359.209876543211_dp)*Hr4(1,0,1,0) + &
                        (71280._dp)*Hr4(1,0,1,0) + z**5*(71280._dp)*Hr4(1,0,1,0) + (-175824.00000000003_dp)*Hr4(1,0,1,1) + &
                        z**5*(-175824._dp)*Hr4(1,0,1,1) + z**2*(-28944.000000000004_dp)*Hr4(1,0,1,1) + &
                        z**3*(-28944.000000000004_dp)*Hr4(1,0,1,1) + nf*z*(-1999.4074074074074_dp)*Hr4(1,0,1,1) + &
                        nf*z**4*(-1999.4074074074074_dp)*Hr4(1,0,1,1) + nf*z**2*(1209.283950617284_dp)*Hr4(1,0,1,1) + &
                        nf*z**3*(1209.283950617284_dp)*Hr4(1,0,1,1) + nf*(3382.1234567901233_dp)*Hr4(1,0,1,1) + &
                        nf*z**5*(3382.1234567901233_dp)*Hr4(1,0,1,1) + z*(162000._dp)*Hr4(1,0,1,1) + z**4*(162000._dp)*Hr4(1,0,1,1) &
                        + z*(-215568._dp)*Hr4(1,1,0,0) + z**4*(-215568._dp)*Hr4(1,1,0,0) + &
                        nf*z**5*(-6068.543209876544_dp)*Hr4(1,1,0,0) + nf*(-6068.543209876542_dp)*Hr4(1,1,0,0) + &
                        z**3*(-3024.0000000000146_dp)*Hr4(1,1,0,0) + z**2*(-3024._dp)*Hr4(1,1,0,0) + &
                        nf*z**4*(1152.5925925925922_dp)*Hr4(1,1,0,0) + nf*z*(1152.592592592593_dp)*Hr4(1,1,0,0) + &
                        nf*z**3*(1747.950617283951_dp)*Hr4(1,1,0,0) + nf*z**2*(1747.9506172839513_dp)*Hr4(1,1,0,0) + &
                        (270864._dp)*Hr4(1,1,0,0) + z**5*(270864._dp)*Hr4(1,1,0,0) + (-175824.00000000003_dp)*Hr4(1,1,0,1) + &
                        z**5*(-175824._dp)*Hr4(1,1,0,1) + z**2*(-28944.000000000004_dp)*Hr4(1,1,0,1) + &
                        z**3*(-28944.000000000004_dp)*Hr4(1,1,0,1) + nf*z**4*(-2307.5555555555557_dp)*Hr4(1,1,0,1) + &
                        nf*z*(-2307.555555555555_dp)*Hr4(1,1,0,1) + nf*z**2*(1928.2962962962963_dp)*Hr4(1,1,0,1) + &
                        nf*z**3*(1928.2962962962963_dp)*Hr4(1,1,0,1) + nf*(2971.259259259259_dp)*Hr4(1,1,0,1) + &
                        nf*z**5*(2971.259259259259_dp)*Hr4(1,1,0,1) + z*(162000._dp)*Hr4(1,1,0,1) + z**4*(162000._dp)*Hr4(1,1,0,1) + &
                        (-175824.00000000003_dp)*Hr4(1,1,1,0) + z**5*(-175824._dp)*Hr4(1,1,1,0) + &
                        z**2*(-28944.000000000004_dp)*Hr4(1,1,1,0) + z**3*(-28944.000000000004_dp)*Hr4(1,1,1,0) + &
                        nf*z*(-2615.703703703704_dp)*Hr4(1,1,1,0) + nf*z**4*(-2615.703703703704_dp)*Hr4(1,1,1,0) + &
                        nf*(2560.395061728395_dp)*Hr4(1,1,1,0) + nf*z**5*(2560.395061728395_dp)*Hr4(1,1,1,0) + &
                        nf*z**2*(2647.3086419753085_dp)*Hr4(1,1,1,0) + nf*z**3*(2647.3086419753085_dp)*Hr4(1,1,1,0) + &
                        z*(162000._dp)*Hr4(1,1,1,0) + z**4*(162000._dp)*Hr4(1,1,1,0) + nf*(-821.7283950617284_dp)*Hr4(1,1,1,1) + &
                        nf*z**5*(-821.7283950617284_dp)*Hr4(1,1,1,1) + nf*z*(-616.2962962962963_dp)*Hr4(1,1,1,1) + &
                        nf*z**4*(-616.2962962962963_dp)*Hr4(1,1,1,1) + nf*z**2*(1438.0246913580245_dp)*Hr4(1,1,1,1) + &
                        nf*z**3*(1438.0246913580245_dp)*Hr4(1,1,1,1) + z**2*(-124416.00000000001_dp)*Hr5(-1,-1,-1,-1,0) + &
                        (-124416._dp)*Hr5(-1,-1,-1,-1,0) + z*(-124415.99999999997_dp)*Hr5(-1,-1,-1,-1,0) + &
                        z**3*(124416._dp)*Hr5(-1,-1,-1,-1,0) + z**4*(124416._dp)*Hr5(-1,-1,-1,-1,0) + &
                        z**5*(124416._dp)*Hr5(-1,-1,-1,-1,0) + z**3*(-466560._dp)*Hr5(-1,-1,-1,0,0) + &
                        z**4*(-466560._dp)*Hr5(-1,-1,-1,0,0) + z**5*(-466560._dp)*Hr5(-1,-1,-1,0,0) + (466560._dp)*Hr5(-1,-1,-1,0,0) &
                        + z*(466560._dp)*Hr5(-1,-1,-1,0,0) + z**2*(466560._dp)*Hr5(-1,-1,-1,0,0) + &
                        z**3*(-311040._dp)*Hr5(-1,-1,-1,0,1) + z**4*(-311040._dp)*Hr5(-1,-1,-1,0,1) + &
                        z**5*(-311040._dp)*Hr5(-1,-1,-1,0,1) + (311040._dp)*Hr5(-1,-1,-1,0,1) + z*(311040._dp)*Hr5(-1,-1,-1,0,1) + &
                        z**2*(311040._dp)*Hr5(-1,-1,-1,0,1) + z**3*(-196992._dp)*Hr5(-1,-1,0,-1,0) + &
                        z**4*(-196992._dp)*Hr5(-1,-1,0,-1,0) + z**5*(-196992._dp)*Hr5(-1,-1,0,-1,0) + (196992._dp)*Hr5(-1,-1,0,-1,0) &
                        + z*(196992._dp)*Hr5(-1,-1,0,-1,0) + z**2*(196992._dp)*Hr5(-1,-1,0,-1,0) + &
                        z*(-368064.0000000001_dp)*Hr5(-1,-1,0,0,0) + (-368064.00000000006_dp)*Hr5(-1,-1,0,0,0) + &
                        z**2*(-368064._dp)*Hr5(-1,-1,0,0,0) + z**3*(368064.00000000006_dp)*Hr5(-1,-1,0,0,0) + &
                        z**4*(368064.00000000006_dp)*Hr5(-1,-1,0,0,0) + z**5*(368064.00000000006_dp)*Hr5(-1,-1,0,0,0) + &
                        (-414720._dp)*Hr5(-1,-1,0,0,1) + z*(-414720._dp)*Hr5(-1,-1,0,0,1) + z**2*(-414720._dp)*Hr5(-1,-1,0,0,1) + &
                        z**3*(414720._dp)*Hr5(-1,-1,0,0,1) + z**4*(414720._dp)*Hr5(-1,-1,0,0,1) + z**5*(414720._dp)*Hr5(-1,-1,0,0,1) &
                        + (-186624._dp)*Hr5(-1,-1,0,1,0) + z*(-186624._dp)*Hr5(-1,-1,0,1,0) + z**2*(-186624._dp)*Hr5(-1,-1,0,1,0) + &
                        z**3*(186624._dp)*Hr5(-1,-1,0,1,0) + z**4*(186624._dp)*Hr5(-1,-1,0,1,0) + z**5*(186624._dp)*Hr5(-1,-1,0,1,0) &
                        + (-186624._dp)*Hr5(-1,-1,0,1,1) + z*(-186624._dp)*Hr5(-1,-1,0,1,1) + z**2*(-186624._dp)*Hr5(-1,-1,0,1,1) + &
                        z**3*(186624._dp)*Hr5(-1,-1,0,1,1) + z**4*(186624._dp)*Hr5(-1,-1,0,1,1) + z**5*(186624._dp)*Hr5(-1,-1,0,1,1) &
                        + z**3*(-196992._dp)*Hr5(-1,0,-1,-1,0) + z**4*(-196992._dp)*Hr5(-1,0,-1,-1,0) + &
                        z**5*(-196992._dp)*Hr5(-1,0,-1,-1,0) + (196992._dp)*Hr5(-1,0,-1,-1,0) + z*(196992._dp)*Hr5(-1,0,-1,-1,0) + &
                        z**2*(196992._dp)*Hr5(-1,0,-1,-1,0) + z*(-347328.00000000006_dp)*Hr5(-1,0,-1,0,0) + &
                        (-347328._dp)*Hr5(-1,0,-1,0,0) + z**2*(-347327.99999999994_dp)*Hr5(-1,0,-1,0,0) + &
                        z**3*(347328._dp)*Hr5(-1,0,-1,0,0) + z**4*(347328._dp)*Hr5(-1,0,-1,0,0) + z**5*(347328._dp)*Hr5(-1,0,-1,0,0) &
                        + z**2*(-248832.00000000003_dp)*Hr5(-1,0,-1,0,1) + (-248832._dp)*Hr5(-1,0,-1,0,1) + &
                        z*(-248831.99999999994_dp)*Hr5(-1,0,-1,0,1) + z**3*(248832._dp)*Hr5(-1,0,-1,0,1) + &
                        z**4*(248832._dp)*Hr5(-1,0,-1,0,1) + z**5*(248832._dp)*Hr5(-1,0,-1,0,1) + &
                        z**2*(-124416.00000000001_dp)*Hr5(-1,0,0,-1,0) + (-124416._dp)*Hr5(-1,0,0,-1,0) + &
                        z*(-124415.99999999997_dp)*Hr5(-1,0,0,-1,0) + z**3*(124416._dp)*Hr5(-1,0,0,-1,0) + &
                        z**4*(124416._dp)*Hr5(-1,0,0,-1,0) + z**5*(124416._dp)*Hr5(-1,0,0,-1,0) + z**3*(-134784._dp)*Hr5(-1,0,0,0,0) &
                        + z**4*(-134784._dp)*Hr5(-1,0,0,0,0) + z**5*(-134784._dp)*Hr5(-1,0,0,0,0) + (134784._dp)*Hr5(-1,0,0,0,0) + &
                        z*(134784._dp)*Hr5(-1,0,0,0,0) + z**2*(134784._dp)*Hr5(-1,0,0,0,0) + &
                        z**3*(-243647.99999999997_dp)*Hr5(-1,0,0,0,1) + z**4*(-243647.99999999997_dp)*Hr5(-1,0,0,0,1) + &
                        z**5*(-243647.99999999997_dp)*Hr5(-1,0,0,0,1) + (243647.99999999997_dp)*Hr5(-1,0,0,0,1) + &
                        z*(243647.99999999997_dp)*Hr5(-1,0,0,0,1) + z**2*(243647.99999999997_dp)*Hr5(-1,0,0,0,1) + &
                        z**3*(-160704._dp)*Hr5(-1,0,0,1,0) + z**4*(-160704._dp)*Hr5(-1,0,0,1,0) + z**5*(-160704._dp)*Hr5(-1,0,0,1,0) &
                        + (160704._dp)*Hr5(-1,0,0,1,0) + z*(160704._dp)*Hr5(-1,0,0,1,0) + z**2*(160704._dp)*Hr5(-1,0,0,1,0) + &
                        z**3*(-176256._dp)*Hr5(-1,0,0,1,1) + z**4*(-176256._dp)*Hr5(-1,0,0,1,1) + z**5*(-176256._dp)*Hr5(-1,0,0,1,1) &
                        + z*(176255.99999999997_dp)*Hr5(-1,0,0,1,1) + (176256._dp)*Hr5(-1,0,0,1,1) + &
                        z**2*(176256.00000000003_dp)*Hr5(-1,0,0,1,1) + z**3*(-41472._dp)*Hr5(-1,0,1,0,0) + &
                        z**4*(-41472._dp)*Hr5(-1,0,1,0,0) + z**5*(-41472._dp)*Hr5(-1,0,1,0,0) + (41472._dp)*Hr5(-1,0,1,0,0) + &
                        z*(41472._dp)*Hr5(-1,0,1,0,0) + z**2*(41472._dp)*Hr5(-1,0,1,0,0) + z**3*(-51840._dp)*Hr5(-1,0,1,0,1) + &
                        z**4*(-51840._dp)*Hr5(-1,0,1,0,1) + z**5*(-51840._dp)*Hr5(-1,0,1,0,1) + (51840._dp)*Hr5(-1,0,1,0,1) + &
                        z*(51840._dp)*Hr5(-1,0,1,0,1) + z**2*(51840._dp)*Hr5(-1,0,1,0,1) + z**3*(-51840._dp)*Hr5(-1,0,1,1,0) + &
                        z**4*(-51840._dp)*Hr5(-1,0,1,1,0) + z**5*(-51840._dp)*Hr5(-1,0,1,1,0) + (51840._dp)*Hr5(-1,0,1,1,0) + &
                        z*(51840._dp)*Hr5(-1,0,1,1,0) + z**2*(51840._dp)*Hr5(-1,0,1,1,0) + z**5*(-238464._dp)*Hr5(0,-1,-1,-1,0) + &
                        z**3*(-217728._dp)*Hr5(0,-1,-1,-1,0) + z**4*(-176256._dp)*Hr5(0,-1,-1,-1,0) + &
                        nf*z**2*(-4608._dp)*Hr5(0,-1,-1,-1,0) + nf*z**3*(-4608._dp)*Hr5(0,-1,-1,-1,0) + &
                        nf*z*(4608._dp)*Hr5(0,-1,-1,-1,0) + nf*z**4*(4608._dp)*Hr5(0,-1,-1,-1,0) + (155520._dp)*Hr5(0,-1,-1,-1,0) + &
                        z**2*(176255.99999999997_dp)*Hr5(0,-1,-1,-1,0) + z*(217728._dp)*Hr5(0,-1,-1,-1,0) + &
                        z**2*(-648000._dp)*Hr5(0,-1,-1,0,0) + (-305856._dp)*Hr5(0,-1,-1,0,0) + &
                        z*(-46656.00000000006_dp)*Hr5(0,-1,-1,0,0) + nf*z*(-12660.148148148146_dp)*Hr5(0,-1,-1,0,0) + &
                        nf*z**4*(-12660.148148148146_dp)*Hr5(0,-1,-1,0,0) + nf*z**2*(12660.148148148146_dp)*Hr5(0,-1,-1,0,0) + &
                        nf*z**3*(12660.148148148146_dp)*Hr5(0,-1,-1,0,0) + z**3*(46656.000000000015_dp)*Hr5(0,-1,-1,0,0) + &
                        z**5*(388800._dp)*Hr5(0,-1,-1,0,0) + z**4*(648000._dp)*Hr5(0,-1,-1,0,0) + &
                        z**2*(-373248._dp)*Hr5(0,-1,-1,0,1) + (-207360._dp)*Hr5(0,-1,-1,0,1) + &
                        z*(-124415.99999999997_dp)*Hr5(0,-1,-1,0,1) + nf*z*(-7210.666666666667_dp)*Hr5(0,-1,-1,0,1) + &
                        nf*z**4*(-7210.666666666667_dp)*Hr5(0,-1,-1,0,1) + nf*z**2*(7210.666666666667_dp)*Hr5(0,-1,-1,0,1) + &
                        nf*z**3*(7210.666666666667_dp)*Hr5(0,-1,-1,0,1) + z**3*(124415.99999999997_dp)*Hr5(0,-1,-1,0,1) + &
                        z**5*(290304._dp)*Hr5(0,-1,-1,0,1) + z**4*(373248._dp)*Hr5(0,-1,-1,0,1) + &
                        z**2*(-290304._dp)*Hr5(0,-1,0,-1,0) + (-145152._dp)*Hr5(0,-1,0,-1,0) + z*(-41472._dp)*Hr5(0,-1,0,-1,0) + &
                        nf*z*(-5636.740740740741_dp)*Hr5(0,-1,0,-1,0) + nf*z**4*(-5636.740740740741_dp)*Hr5(0,-1,0,-1,0) + &
                        nf*z**2*(5636.740740740741_dp)*Hr5(0,-1,0,-1,0) + nf*z**3*(5636.740740740741_dp)*Hr5(0,-1,0,-1,0) + &
                        z**3*(41472.00000000001_dp)*Hr5(0,-1,0,-1,0) + z**5*(186624._dp)*Hr5(0,-1,0,-1,0) + &
                        z**4*(290304._dp)*Hr5(0,-1,0,-1,0) + z**4*(-505439.99999999994_dp)*Hr5(0,-1,0,0,0) + &
                        z**5*(-251424._dp)*Hr5(0,-1,0,0,0) + z*(-44064._dp)*Hr5(0,-1,0,0,0) + &
                        nf*z**2*(-9367.703703703704_dp)*Hr5(0,-1,0,0,0) + nf*z**3*(-9367.703703703704_dp)*Hr5(0,-1,0,0,0) + &
                        nf*z*(9367.703703703704_dp)*Hr5(0,-1,0,0,0) + nf*z**4*(9367.703703703704_dp)*Hr5(0,-1,0,0,0) + &
                        z**3*(44064.00000000001_dp)*Hr5(0,-1,0,0,0) + (209952._dp)*Hr5(0,-1,0,0,0) + &
                        z**2*(505439.99999999994_dp)*Hr5(0,-1,0,0,0) + z**4*(-487295.99999999994_dp)*Hr5(0,-1,0,0,1) + &
                        z**5*(-321408._dp)*Hr5(0,-1,0,0,1) + z**3*(-72576.00000000001_dp)*Hr5(0,-1,0,0,1) + &
                        nf*z**2*(-8896._dp)*Hr5(0,-1,0,0,1) + nf*z**3*(-8896._dp)*Hr5(0,-1,0,0,1) + nf*z*(8896._dp)*Hr5(0,-1,0,0,1) &
                        + nf*z**4*(8896._dp)*Hr5(0,-1,0,0,1) + z*(72576._dp)*Hr5(0,-1,0,0,1) + (238464._dp)*Hr5(0,-1,0,0,1) + &
                        z**2*(487295.99999999994_dp)*Hr5(0,-1,0,0,1) + z**4*(-191808._dp)*Hr5(0,-1,0,1,0) + &
                        z**5*(-139968._dp)*Hr5(0,-1,0,1,0) + z**3*(-46655.99999999998_dp)*Hr5(0,-1,0,1,0) + &
                        nf*z**2*(-2602.6666666666665_dp)*Hr5(0,-1,0,1,0) + nf*z**3*(-2602.6666666666665_dp)*Hr5(0,-1,0,1,0) + &
                        nf*z*(2602.6666666666665_dp)*Hr5(0,-1,0,1,0) + nf*z**4*(2602.6666666666665_dp)*Hr5(0,-1,0,1,0) + &
                        z*(46656._dp)*Hr5(0,-1,0,1,0) + (98496._dp)*Hr5(0,-1,0,1,0) + z**2*(191808._dp)*Hr5(0,-1,0,1,0) + &
                        z**3*(-186624.00000000003_dp)*Hr5(0,-1,0,1,1) + z**5*(-186624._dp)*Hr5(0,-1,0,1,1) + &
                        z**4*(-103680.00000000001_dp)*Hr5(0,-1,0,1,1) + (103680._dp)*Hr5(0,-1,0,1,1) + &
                        z**2*(103680._dp)*Hr5(0,-1,0,1,1) + z*(186624._dp)*Hr5(0,-1,0,1,1) + &
                        z**2*(-238464.00000000006_dp)*Hr5(0,0,-1,-1,0) + z**3*(-217728._dp)*Hr5(0,0,-1,-1,0) + &
                        (-72576._dp)*Hr5(0,0,-1,-1,0) + nf*z*(-12105.481481481482_dp)*Hr5(0,0,-1,-1,0) + &
                        nf*z**4*(-8725.333333333334_dp)*Hr5(0,0,-1,-1,0) + nf*z**2*(8725.333333333336_dp)*Hr5(0,0,-1,-1,0) + &
                        nf*z**3*(12105.481481481485_dp)*Hr5(0,0,-1,-1,0) + z**5*(176256._dp)*Hr5(0,0,-1,-1,0) + &
                        z*(217728._dp)*Hr5(0,0,-1,-1,0) + z**4*(238464._dp)*Hr5(0,0,-1,-1,0) + z**4*(-544320._dp)*Hr5(0,0,-1,0,0) + &
                        z**5*(-212544._dp)*Hr5(0,0,-1,0,0) + z*(-160704._dp)*Hr5(0,0,-1,0,0) + &
                        nf*z**2*(-20250.074074074077_dp)*Hr5(0,0,-1,0,0) + nf*z**3*(-6366.814814814816_dp)*Hr5(0,0,-1,0,0) + &
                        nf*z*(6366.814814814815_dp)*Hr5(0,0,-1,0,0) + nf*z**4*(20250.074074074073_dp)*Hr5(0,0,-1,0,0) + &
                        (119232._dp)*Hr5(0,0,-1,0,0) + z**3*(160703.99999999997_dp)*Hr5(0,0,-1,0,0) + &
                        z**2*(544320._dp)*Hr5(0,0,-1,0,0) + z**4*(-269568._dp)*Hr5(0,0,-1,0,1) + z*(-186624._dp)*Hr5(0,0,-1,0,1) + &
                        z**5*(-165888._dp)*Hr5(0,0,-1,0,1) + nf*z**3*(-12309.333333333334_dp)*Hr5(0,0,-1,0,1) + &
                        nf*z**2*(-6592._dp)*Hr5(0,0,-1,0,1) + nf*z**4*(6592._dp)*Hr5(0,0,-1,0,1) + &
                        nf*z*(12309.333333333334_dp)*Hr5(0,0,-1,0,1) + (82944._dp)*Hr5(0,0,-1,0,1) + &
                        z**3*(186624._dp)*Hr5(0,0,-1,0,1) + z**2*(269568._dp)*Hr5(0,0,-1,0,1) + z**4*(-196992._dp)*Hr5(0,0,0,-1,0) + &
                        z**3*(-93312.00000000001_dp)*Hr5(0,0,0,-1,0) + z**5*(-72576._dp)*Hr5(0,0,0,-1,0) + &
                        nf*z*(-13575.111111111111_dp)*Hr5(0,0,0,-1,0) + nf*z**2*(-12373.333333333334_dp)*Hr5(0,0,0,-1,0) + &
                        (10368._dp)*Hr5(0,0,0,-1,0) + nf*z**4*(12373.333333333334_dp)*Hr5(0,0,0,-1,0) + &
                        nf*z**3*(13575.111111111111_dp)*Hr5(0,0,0,-1,0) + z*(93312._dp)*Hr5(0,0,0,-1,0) + &
                        z**2*(196992._dp)*Hr5(0,0,0,-1,0) + z**4*(-373248._dp)*Hr5(0,0,0,0,0) + &
                        z**3*(-57023.999999999956_dp)*Hr5(0,0,0,0,0) + nf*z**2*(-30539.85185185185_dp)*Hr5(0,0,0,0,0) + &
                        nf*z*(-11003.25925925926_dp)*Hr5(0,0,0,0,0) + nf**2*z**3*(-312.8888888888889_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z**4*(-312.8888888888889_dp)*Hr5(0,0,0,0,0) + nf**2*z*(312.8888888888889_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z**2*(312.8888888888889_dp)*Hr5(0,0,0,0,0) + nf*z**3*(11003.25925925926_dp)*Hr5(0,0,0,0,0) + &
                        nf*z**4*(30539.85185185185_dp)*Hr5(0,0,0,0,0) + z*(57024._dp)*Hr5(0,0,0,0,0) + &
                        z**5*(57024._dp)*Hr5(0,0,0,0,0) + z**2*(373248._dp)*Hr5(0,0,0,0,0) + z**4*(-539136._dp)*Hr5(0,0,0,0,1) + &
                        z*(-124416.00000000004_dp)*Hr5(0,0,0,0,1) + nf*z**2*(-24877.037037037036_dp)*Hr5(0,0,0,0,1) + &
                        nf*z**3*(-13362.96296296296_dp)*Hr5(0,0,0,0,1) + nf*z*(13362.962962962962_dp)*Hr5(0,0,0,0,1) + &
                        nf*z**4*(24877.037037037036_dp)*Hr5(0,0,0,0,1) + z**3*(124416.00000000007_dp)*Hr5(0,0,0,0,1) + &
                        z**5*(155520._dp)*Hr5(0,0,0,0,1) + z**2*(539136._dp)*Hr5(0,0,0,0,1) + &
                        z**4*(-476928.0000000002_dp)*Hr5(0,0,0,1,0) + z**3*(-25920.000000000127_dp)*Hr5(0,0,0,1,0) + &
                        nf*z**2*(-15493.925925925927_dp)*Hr5(0,0,0,1,0) + nf*z**3*(-2852.740740740742_dp)*Hr5(0,0,0,1,0) + &
                        nf*z*(2852.740740740741_dp)*Hr5(0,0,0,1,0) + nf*z**4*(15493.925925925927_dp)*Hr5(0,0,0,1,0) + &
                        z*(25920._dp)*Hr5(0,0,0,1,0) + (41472._dp)*Hr5(0,0,0,1,0) + z**5*(181440._dp)*Hr5(0,0,0,1,0) + &
                        z**2*(476928._dp)*Hr5(0,0,0,1,0) + z**4*(-440640.00000000006_dp)*Hr5(0,0,0,1,1) + &
                        z**3*(-181439.99999999997_dp)*Hr5(0,0,0,1,1) + nf*z*(-9204.148148148148_dp)*Hr5(0,0,0,1,1) + &
                        nf*z**2*(-7540.148148148148_dp)*Hr5(0,0,0,1,1) + nf*z**4*(7540.148148148148_dp)*Hr5(0,0,0,1,1) + &
                        nf*z**3*(9204.148148148148_dp)*Hr5(0,0,0,1,1) + (57024._dp)*Hr5(0,0,0,1,1) + &
                        z*(181440.00000000003_dp)*Hr5(0,0,0,1,1) + z**5*(233280._dp)*Hr5(0,0,0,1,1) + &
                        z**2*(440640.00000000006_dp)*Hr5(0,0,0,1,1) + z**4*(-440640.00000000006_dp)*Hr5(0,0,1,0,0) + &
                        z**3*(-46656.000000000015_dp)*Hr5(0,0,1,0,0) + nf*z**2*(-13757.62962962963_dp)*Hr5(0,0,1,0,0) + &
                        nf*z*(-835.5555555555557_dp)*Hr5(0,0,1,0,0) + nf*z**3*(835.5555555555563_dp)*Hr5(0,0,1,0,0) + &
                        nf*z**4*(13757.62962962963_dp)*Hr5(0,0,1,0,0) + z*(46655.999999999985_dp)*Hr5(0,0,1,0,0) + &
                        (98496._dp)*Hr5(0,0,1,0,0) + z**5*(150336._dp)*Hr5(0,0,1,0,0) + z**2*(440640.00000000006_dp)*Hr5(0,0,1,0,0) &
                        + z**4*(-264383.99999999994_dp)*Hr5(0,0,1,0,1) + z**3*(-57023.999999999956_dp)*Hr5(0,0,1,0,1) + &
                        nf*z*(-7329.185185185185_dp)*Hr5(0,0,1,0,1) + nf*z**2*(-3915.8518518518517_dp)*Hr5(0,0,1,0,1) + &
                        nf*z**4*(3915.8518518518517_dp)*Hr5(0,0,1,0,1) + nf*z**3*(7329.185185185185_dp)*Hr5(0,0,1,0,1) + &
                        z*(57024.00000000003_dp)*Hr5(0,0,1,0,1) + (108864._dp)*Hr5(0,0,1,0,1) + z**5*(171072._dp)*Hr5(0,0,1,0,1) + &
                        z**2*(264384._dp)*Hr5(0,0,1,0,1) + z**4*(-264383.99999999994_dp)*Hr5(0,0,1,1,0) + &
                        z**3*(-57023.999999999956_dp)*Hr5(0,0,1,1,0) + nf*z*(-6997.333333333333_dp)*Hr5(0,0,1,1,0) + &
                        nf*z**2*(-4864.000000000001_dp)*Hr5(0,0,1,1,0) + nf*z**4*(4864._dp)*Hr5(0,0,1,1,0) + &
                        nf*z**3*(6997.333333333333_dp)*Hr5(0,0,1,1,0) + z*(57024.00000000003_dp)*Hr5(0,0,1,1,0) + &
                        (108864._dp)*Hr5(0,0,1,1,0) + z**5*(171072._dp)*Hr5(0,0,1,1,0) + z**2*(264384._dp)*Hr5(0,0,1,1,0) + &
                        z*(-155520._dp)*Hr5(0,0,1,1,1) + z**4*(-31104.000000000004_dp)*Hr5(0,0,1,1,1) + &
                        nf*z**3*(-905.4814814814815_dp)*Hr5(0,0,1,1,1) + nf*z**4*(-265.48148148148147_dp)*Hr5(0,0,1,1,1) + &
                        nf*z**2*(265.48148148148147_dp)*Hr5(0,0,1,1,1) + nf*z*(905.4814814814815_dp)*Hr5(0,0,1,1,1) + &
                        z**2*(31104.00000000003_dp)*Hr5(0,0,1,1,1) + (93312._dp)*Hr5(0,0,1,1,1) + z**5*(93312._dp)*Hr5(0,0,1,1,1) + &
                        z**3*(155520.00000000003_dp)*Hr5(0,0,1,1,1) + z**2*(-124415.99999999997_dp)*Hr5(0,1,0,-1,0) + &
                        (-41472._dp)*Hr5(0,1,0,-1,0) + z*(-41472._dp)*Hr5(0,1,0,-1,0) + z**5*(-41472._dp)*Hr5(0,1,0,-1,0) + &
                        nf*z**3*(-2313.4814814814813_dp)*Hr5(0,1,0,-1,0) + nf*z**4*(-2313.4814814814813_dp)*Hr5(0,1,0,-1,0) + &
                        nf*z*(2313.4814814814813_dp)*Hr5(0,1,0,-1,0) + nf*z**2*(2313.4814814814813_dp)*Hr5(0,1,0,-1,0) + &
                        z**3*(41472.00000000001_dp)*Hr5(0,1,0,-1,0) + z**4*(124416._dp)*Hr5(0,1,0,-1,0) + &
                        z**4*(-329184._dp)*Hr5(0,1,0,0,0) + z**3*(-64799.99999999997_dp)*Hr5(0,1,0,0,0) + &
                        nf*z*(-5550.222222222223_dp)*Hr5(0,1,0,0,0) + nf*z**2*(-5550.222222222223_dp)*Hr5(0,1,0,0,0) + &
                        nf*z**3*(5550.222222222223_dp)*Hr5(0,1,0,0,0) + nf*z**4*(5550.222222222223_dp)*Hr5(0,1,0,0,0) + &
                        z*(64799.999999999985_dp)*Hr5(0,1,0,0,0) + (127008._dp)*Hr5(0,1,0,0,0) + z**5*(137376._dp)*Hr5(0,1,0,0,0) + &
                        z**2*(329184.00000000006_dp)*Hr5(0,1,0,0,0) + z**4*(-274752._dp)*Hr5(0,1,0,0,1) + &
                        z*(-25920._dp)*Hr5(0,1,0,0,1) + nf*z*(-3717.925925925926_dp)*Hr5(0,1,0,0,1) + &
                        nf*z**2*(-3717.925925925926_dp)*Hr5(0,1,0,0,1) + nf*z**3*(3717.925925925926_dp)*Hr5(0,1,0,0,1) + &
                        nf*z**4*(3717.925925925926_dp)*Hr5(0,1,0,0,1) + z**3*(25919.999999999975_dp)*Hr5(0,1,0,0,1) + &
                        (139968._dp)*Hr5(0,1,0,0,1) + z**5*(160704._dp)*Hr5(0,1,0,0,1) + z**2*(274752._dp)*Hr5(0,1,0,0,1) + &
                        z**4*(-202176.00000000003_dp)*Hr5(0,1,0,1,0) + z*(-77760._dp)*Hr5(0,1,0,1,0) + &
                        nf*z*(-1645.037037037037_dp)*Hr5(0,1,0,1,0) + nf*z**2*(-1645.037037037037_dp)*Hr5(0,1,0,1,0) + &
                        nf*z**3*(1645.037037037037_dp)*Hr5(0,1,0,1,0) + nf*z**4*(1645.037037037037_dp)*Hr5(0,1,0,1,0) + &
                        z**3*(77760._dp)*Hr5(0,1,0,1,0) + (129600.00000000001_dp)*Hr5(0,1,0,1,0) + z**5*(150336._dp)*Hr5(0,1,0,1,0) &
                        + z**2*(202176.00000000003_dp)*Hr5(0,1,0,1,0) + z*(-124416._dp)*Hr5(0,1,0,1,1) + &
                        nf*z**3*(-1185.1851851851852_dp)*Hr5(0,1,0,1,1) + nf*z**4*(-1185.1851851851852_dp)*Hr5(0,1,0,1,1) + &
                        z**2*(1.4551915228366852e-11_dp)*Hr5(0,1,0,1,1) + nf*z*(1185.1851851851852_dp)*Hr5(0,1,0,1,1) + &
                        nf*z**2*(1185.1851851851852_dp)*Hr5(0,1,0,1,1) + (51840._dp)*Hr5(0,1,0,1,1) + &
                        z**5*(72576._dp)*Hr5(0,1,0,1,1) + z**3*(124416.00000000001_dp)*Hr5(0,1,0,1,1) + &
                        z**4*(-290303.99999999994_dp)*Hr5(0,1,1,0,0) + z*(-20735.999999999913_dp)*Hr5(0,1,1,0,0) + &
                        nf*z*(-4244.148148148148_dp)*Hr5(0,1,1,0,0) + nf*z**2*(-4244.148148148148_dp)*Hr5(0,1,1,0,0) + &
                        nf*z**3*(4244.148148148148_dp)*Hr5(0,1,1,0,0) + nf*z**4*(4244.148148148148_dp)*Hr5(0,1,1,0,0) + &
                        z**3*(20736.00000000004_dp)*Hr5(0,1,1,0,0) + (145152._dp)*Hr5(0,1,1,0,0) + z**5*(165888._dp)*Hr5(0,1,1,0,0) &
                        + z**2*(290304.0000000001_dp)*Hr5(0,1,1,0,0) + z*(-124416._dp)*Hr5(0,1,1,0,1) + &
                        nf*z**3*(-568.8888888888889_dp)*Hr5(0,1,1,0,1) + nf*z**4*(-568.8888888888889_dp)*Hr5(0,1,1,0,1) + &
                        z**2*(1.4551915228366852e-11_dp)*Hr5(0,1,1,0,1) + nf*z*(568.8888888888889_dp)*Hr5(0,1,1,0,1) + &
                        nf*z**2*(568.8888888888889_dp)*Hr5(0,1,1,0,1) + (51840._dp)*Hr5(0,1,1,0,1) + z**5*(72576._dp)*Hr5(0,1,1,0,1) &
                        + z**3*(124416.00000000001_dp)*Hr5(0,1,1,0,1) + z*(-124416._dp)*Hr5(0,1,1,1,0) + &
                        nf*z*(-47.407407407407405_dp)*Hr5(0,1,1,1,0) + nf*z**2*(-47.407407407407405_dp)*Hr5(0,1,1,1,0) + &
                        z**2*(1.4551915228366852e-11_dp)*Hr5(0,1,1,1,0) + nf*z**3*(47.407407407407405_dp)*Hr5(0,1,1,1,0) + &
                        nf*z**4*(47.407407407407405_dp)*Hr5(0,1,1,1,0) + (51840._dp)*Hr5(0,1,1,1,0) + &
                        z**5*(72576._dp)*Hr5(0,1,1,1,0) + z**3*(124416.00000000001_dp)*Hr5(0,1,1,1,0) + (-62208._dp)*Hr5(0,1,1,1,1) &
                        + z**2*(-62208._dp)*Hr5(0,1,1,1,1) + z**5*(-62208._dp)*Hr5(0,1,1,1,1) + &
                        z**3*(-62207.999999999985_dp)*Hr5(0,1,1,1,1) + nf*z*(-1232.5925925925926_dp)*Hr5(0,1,1,1,1) + &
                        nf*z**2*(-1232.5925925925926_dp)*Hr5(0,1,1,1,1) + nf*z**3*(1232.5925925925926_dp)*Hr5(0,1,1,1,1) + &
                        nf*z**4*(1232.5925925925926_dp)*Hr5(0,1,1,1,1) + z*(62207.99999999997_dp)*Hr5(0,1,1,1,1) + &
                        z**4*(62208._dp)*Hr5(0,1,1,1,1) + z**4*(-62208._dp)*Hr5(1,0,-1,-1,0) + &
                        z*(-62207.99999999997_dp)*Hr5(1,0,-1,-1,0) + z**3*(62207.999999999985_dp)*Hr5(1,0,-1,-1,0) + &
                        (62208._dp)*Hr5(1,0,-1,-1,0) + z**2*(62208._dp)*Hr5(1,0,-1,-1,0) + z**5*(62208._dp)*Hr5(1,0,-1,-1,0) + &
                        z**3*(-51840.00000000001_dp)*Hr5(1,0,-1,0,0) + (-51840._dp)*Hr5(1,0,-1,0,0) + &
                        z**2*(-51840._dp)*Hr5(1,0,-1,0,0) + z**5*(-51840._dp)*Hr5(1,0,-1,0,0) + z**4*(51840._dp)*Hr5(1,0,-1,0,0) + &
                        z*(51840.000000000015_dp)*Hr5(1,0,-1,0,0) + (-41472._dp)*Hr5(1,0,-1,0,1) + z**2*(-41472._dp)*Hr5(1,0,-1,0,1) &
                        + z**5*(-41472._dp)*Hr5(1,0,-1,0,1) + z**3*(-41471.99999999999_dp)*Hr5(1,0,-1,0,1) + &
                        z*(41471.999999999985_dp)*Hr5(1,0,-1,0,1) + z**4*(41472._dp)*Hr5(1,0,-1,0,1) + &
                        z**3*(-72576.00000000001_dp)*Hr5(1,0,0,-1,0) + (-72576._dp)*Hr5(1,0,0,-1,0) + &
                        z**2*(-72576._dp)*Hr5(1,0,0,-1,0) + z**5*(-72576._dp)*Hr5(1,0,0,-1,0) + z**4*(72576._dp)*Hr5(1,0,0,-1,0) + &
                        z*(72576.00000000003_dp)*Hr5(1,0,0,-1,0) + z*(-93312.00000000003_dp)*Hr5(1,0,0,0,0) + &
                        z**4*(-93312._dp)*Hr5(1,0,0,0,0) + (93312._dp)*Hr5(1,0,0,0,0) + z**2*(93312._dp)*Hr5(1,0,0,0,0) + &
                        z**5*(93312._dp)*Hr5(1,0,0,0,0) + z**3*(93312.00000000003_dp)*Hr5(1,0,0,0,0) + &
                        z*(-155520._dp)*Hr5(1,0,0,0,1) + z**4*(-155520._dp)*Hr5(1,0,0,0,1) + (155520._dp)*Hr5(1,0,0,0,1) + &
                        z**2*(155520._dp)*Hr5(1,0,0,0,1) + z**3*(155520._dp)*Hr5(1,0,0,0,1) + z**5*(155520._dp)*Hr5(1,0,0,0,1) + &
                        z*(-145152.00000000006_dp)*Hr5(1,0,0,1,0) + z**4*(-145152._dp)*Hr5(1,0,0,1,0) + (145152._dp)*Hr5(1,0,0,1,0) &
                        + z**2*(145152._dp)*Hr5(1,0,0,1,0) + z**5*(145152._dp)*Hr5(1,0,0,1,0) + &
                        z**3*(145152.00000000003_dp)*Hr5(1,0,0,1,0) + z*(-103680.00000000003_dp)*Hr5(1,0,0,1,1) + &
                        z**4*(-103680._dp)*Hr5(1,0,0,1,1) + (103680._dp)*Hr5(1,0,0,1,1) + z**2*(103680._dp)*Hr5(1,0,0,1,1) + &
                        z**5*(103680._dp)*Hr5(1,0,0,1,1) + z**3*(103680.00000000001_dp)*Hr5(1,0,0,1,1) + &
                        z*(-155520._dp)*Hr5(1,0,1,0,0) + z**4*(-155520._dp)*Hr5(1,0,1,0,0) + (155520._dp)*Hr5(1,0,1,0,0) + &
                        z**2*(155520._dp)*Hr5(1,0,1,0,0) + z**3*(155520._dp)*Hr5(1,0,1,0,0) + z**5*(155520._dp)*Hr5(1,0,1,0,0) + &
                        z**4*(-62208._dp)*Hr5(1,0,1,0,1) + z*(-62207.99999999997_dp)*Hr5(1,0,1,0,1) + &
                        z**3*(62207.999999999985_dp)*Hr5(1,0,1,0,1) + (62208._dp)*Hr5(1,0,1,0,1) + z**2*(62208._dp)*Hr5(1,0,1,0,1) + &
                        z**5*(62208._dp)*Hr5(1,0,1,0,1) + z**4*(-62208._dp)*Hr5(1,0,1,1,0) + &
                        z*(-62207.99999999997_dp)*Hr5(1,0,1,1,0) + z**3*(62207.999999999985_dp)*Hr5(1,0,1,1,0) + &
                        (62208._dp)*Hr5(1,0,1,1,0) + z**2*(62208._dp)*Hr5(1,0,1,1,0) + z**5*(62208._dp)*Hr5(1,0,1,1,0) + &
                        (-62208._dp)*Hr5(1,0,1,1,1) + z**2*(-62208._dp)*Hr5(1,0,1,1,1) + z**5*(-62208._dp)*Hr5(1,0,1,1,1) + &
                        z**3*(-62207.999999999985_dp)*Hr5(1,0,1,1,1) + z*(62207.99999999997_dp)*Hr5(1,0,1,1,1) + &
                        z**4*(62208._dp)*Hr5(1,0,1,1,1) + (-82944._dp)*Hr5(1,1,0,-1,0) + z**2*(-82944._dp)*Hr5(1,1,0,-1,0) + &
                        z**5*(-82944._dp)*Hr5(1,1,0,-1,0) + z**3*(-82943.99999999999_dp)*Hr5(1,1,0,-1,0) + &
                        z*(82943.99999999997_dp)*Hr5(1,1,0,-1,0) + z**4*(82944._dp)*Hr5(1,1,0,-1,0) + &
                        z*(-186624.00000000006_dp)*Hr5(1,1,0,0,0) + z**4*(-186624._dp)*Hr5(1,1,0,0,0) + (186624._dp)*Hr5(1,1,0,0,0) &
                        + z**2*(186624._dp)*Hr5(1,1,0,0,0) + z**5*(186624._dp)*Hr5(1,1,0,0,0) + &
                        z**3*(186624.00000000006_dp)*Hr5(1,1,0,0,0) + z*(-114048._dp)*Hr5(1,1,0,0,1) + &
                        z**4*(-114048._dp)*Hr5(1,1,0,0,1) + (114048._dp)*Hr5(1,1,0,0,1) + z**2*(114048._dp)*Hr5(1,1,0,0,1) + &
                        z**3*(114048._dp)*Hr5(1,1,0,0,1) + z**5*(114048._dp)*Hr5(1,1,0,0,1) + z**4*(-62208._dp)*Hr5(1,1,0,1,0) + &
                        z*(-62207.99999999997_dp)*Hr5(1,1,0,1,0) + z**3*(62207.999999999985_dp)*Hr5(1,1,0,1,0) + &
                        (62208._dp)*Hr5(1,1,0,1,0) + z**2*(62208._dp)*Hr5(1,1,0,1,0) + z**5*(62208._dp)*Hr5(1,1,0,1,0) + &
                        (-62208._dp)*Hr5(1,1,0,1,1) + z**2*(-62208._dp)*Hr5(1,1,0,1,1) + z**5*(-62208._dp)*Hr5(1,1,0,1,1) + &
                        z**3*(-62207.999999999985_dp)*Hr5(1,1,0,1,1) + z*(62207.99999999997_dp)*Hr5(1,1,0,1,1) + &
                        z**4*(62208._dp)*Hr5(1,1,0,1,1) + z**4*(-124416._dp)*Hr5(1,1,1,0,0) + &
                        z*(-124415.99999999994_dp)*Hr5(1,1,1,0,0) + z**3*(124415.99999999997_dp)*Hr5(1,1,1,0,0) + &
                        (124416._dp)*Hr5(1,1,1,0,0) + z**2*(124416._dp)*Hr5(1,1,1,0,0) + z**5*(124416._dp)*Hr5(1,1,1,0,0) + &
                        (-62208._dp)*Hr5(1,1,1,0,1) + z**2*(-62208._dp)*Hr5(1,1,1,0,1) + z**5*(-62208._dp)*Hr5(1,1,1,0,1) + &
                        z**3*(-62207.999999999985_dp)*Hr5(1,1,1,0,1) + z*(62207.99999999997_dp)*Hr5(1,1,1,0,1) + &
                        z**4*(62208._dp)*Hr5(1,1,1,0,1) + (-62208._dp)*Hr5(1,1,1,1,0) + z**2*(-62208._dp)*Hr5(1,1,1,1,0) + &
                        z**5*(-62208._dp)*Hr5(1,1,1,1,0) + z**3*(-62207.999999999985_dp)*Hr5(1,1,1,1,0) + &
                        z*(62207.99999999997_dp)*Hr5(1,1,1,1,0) + z**4*(62208._dp)*Hr5(1,1,1,1,0))/(z*(z + (1._dp))*((-1._dp) + &
                        z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-62428.9691251482_dp) + nf**2*(-921.0798595122442_dp) + nf**3*(7.609469445639825_dp) + &
                            nf*(24230.497556856026_dp)

                        else
                            Ibar_select = (-62428.9691251482_dp) + nf**2*(-921.0798595122442_dp) + nf**3*(7.609469445639825_dp) + &
                            nf*(24230.497556856026_dp)

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-10270.141656414533_dp) + nf**2*(-143.43961619616834_dp) + &
                            nf**3*(4.938271604938271_dp) + (56748.451098868834_dp)

                        else
                            Ibar_select = nf*(-10270.141656414533_dp) + nf**2*(-143.43961619616834_dp) + &
                            nf**3*(4.938271604938271_dp) + (56748.451098868834_dp)

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-15305.86954732275_dp) + nf**2*(433.77777777777777_dp) + (45857.644590167445_dp)
                        else
                            Ibar_select = nf*(-15305.86954732275_dp) + nf**2*(433.77777777777777_dp) + (45857.644590167445_dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = nf*(2688._dp) + (7254.361352783031_dp)
                        else
                            Ibar_select = nf*(2688._dp) + (7254.361352783031_dp)
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
                        Ibar_select = (z**5*(-1.3052591732880324e6_dp) + (-169070.83120486594_dp) + nf*z*(-25319.910347967227_dp) + &
                        nf*(-8515.925588693572_dp) + nf**2*z**2*(-674.2428664828883_dp) + nf**2*z**3*(-403.9658353157456_dp) + &
                        nf**2*z*(-399.8622629214053_dp) + nf**3*z*(-12.538946687280479_dp) + nf**3*(-5.267489711934156_dp) + &
                        nf**3*z**5*(-5.267489711934156_dp) + nf**3*z**2*(2.2100245061351576_dp) + nf**3*z**4*(3.0574652057989966_dp) &
                        + nf**2*z**4*(13.823912164138164_dp) + nf**3*z**3*(17.806436399214636_dp) + nf**2*(660.4189543187501_dp) + &
                        nf**2*z**5*(803.8280982371507_dp) + nf*z**2*(919.8013513631277_dp) + nf*z**4*(5692.066102725642_dp) + &
                        nf*z**5*(5962.457074307651_dp) + z**4*(6700.774589389937_dp) + nf*z**3*(19876.74185582451_dp) + &
                        z**3*(24793.579681110714_dp) + z**2*(195189.7067473448_dp) + z*(1.2985484593081032e6_dp) + &
                        ((-337375.7357449194_dp) + z*(-297581.4907997271_dp) + z**2*(-78713.14360116946_dp) + &
                        nf*(z**5*(-4227.494090351932_dp) + z**2*(-2186.963534340701_dp) + z**4*(-145.566511002075_dp) + &
                        z*(145.56651100207452_dp) + z**3*(2186.9635343407003_dp) + (4227.494090351932_dp)) + &
                        z**3*(78713.14360116946_dp) + z**4*(297581.4907997271_dp) + z**5*(337375.7357449194_dp))*Hr1(-1) + &
                        z*(-306704.50824196974_dp)*Hr1(0) + nf*z**5*(-28138.84557244405_dp)*Hr1(0) + (-27052.331341925725_dp)*Hr1(0) &
                        + nf*z**2*(-13503.458160557462_dp)*Hr1(0) + nf*(-7838.454837709223_dp)*Hr1(0) + &
                        nf**2*z**4*(-809.147053841506_dp)*Hr1(0) + nf**2*z**3*(-472.651170595414_dp)*Hr1(0) + &
                        nf**3*z*(-12.641975308641975_dp)*Hr1(0) + nf**3*z**2*(-12.641975308641975_dp)*Hr1(0) + &
                        nf**3*z**5*(3.1604938271604937_dp)*Hr1(0) + nf**3*z**3*(9.481481481481481_dp)*Hr1(0) + &
                        nf**3*z**4*(12.641975308641975_dp)*Hr1(0) + nf**2*z**5*(83.2702331961591_dp)*Hr1(0) + &
                        nf**2*(190.26611796982166_dp)*Hr1(0) + nf**2*z*(442.71427073258826_dp)*Hr1(0) + &
                        nf**2*z**2*(672.2142692050177_dp)*Hr1(0) + z**4*(6334.2370314133295_dp)*Hr1(0) + &
                        nf*z*(11382.211790976311_dp)*Hr1(0) + nf*z**3*(16792.468332103042_dp)*Hr1(0) + &
                        nf*z**4*(19482.78350389283_dp)*Hr1(0) + z**2*(59855.58291159835_dp)*Hr1(0) + &
                        z**3*(86912.48886884123_dp)*Hr1(0) + z**5*(171579.43435774941_dp)*Hr1(0) + (-181252.89067028777_dp)*Hr1(1) + &
                        z**5*(-181252.89067028777_dp)*Hr1(1) + nf*z**5*(-17943.823154364458_dp)*Hr1(1) + &
                        nf*(-14036.135911566105_dp)*Hr1(1) + nf*z**2*(-4416.197977372031_dp)*Hr1(1) + &
                        nf**2*z**2*(-307.3580246913581_dp)*Hr1(1) + nf**2*z**3*(-254.41975308641972_dp)*Hr1(1) + &
                        nf**3*z**2*(-5.530864197530864_dp)*Hr1(1) + nf**3*z**3*(-5.530864197530864_dp)*Hr1(1) + &
                        nf**3*z*(2.3703703703703702_dp)*Hr1(1) + nf**3*z**4*(2.3703703703703702_dp)*Hr1(1) + &
                        nf**3*(3.1604938271604937_dp)*Hr1(1) + nf**3*z**5*(3.1604938271604937_dp)*Hr1(1) + &
                        nf**2*(68.80658436213992_dp)*Hr1(1) + nf**2*z**5*(95.67078189300412_dp)*Hr1(1) + &
                        nf**2*z*(158.74897119341563_dp)*Hr1(1) + nf**2*z**4*(238.5514403292181_dp)*Hr1(1) + &
                        nf*z**3*(628.3995419815378_dp)*Hr1(1) + nf*z*(17315.42361238293_dp)*Hr1(1) + &
                        nf*z**4*(18452.33388893814_dp)*Hr1(1) + z**3*(27986.405401587206_dp)*Hr1(1) + &
                        z**2*(33890.97441651346_dp)*Hr1(1) + z**4*(172287.76819769168_dp)*Hr1(1) + z*(178192.33721261792_dp)*Hr1(1) &
                        + z**3*(-119382.73483557686_dp)*Hr2(-1,-1) + z**4*(-119382.73483557686_dp)*Hr2(-1,-1) + &
                        z**5*(-119382.73483557686_dp)*Hr2(-1,-1) + z*(119382.73483557685_dp)*Hr2(-1,-1) + &
                        (119382.73483557686_dp)*Hr2(-1,-1) + z**2*(119382.73483557688_dp)*Hr2(-1,-1) + &
                        z**4*(-128021.95617712916_dp)*Hr2(-1,0) + z**5*(-120533.95617712915_dp)*Hr2(-1,0) + &
                        z**2*(-69258.04382287088_dp)*Hr2(-1,0) + nf*z*(-8108.115226337449_dp)*Hr2(-1,0) + &
                        nf*(-6642.567901234568_dp)*Hr2(-1,0) + nf*z**2*(-1465.5473251028807_dp)*Hr2(-1,0) + &
                        nf**2*z**2*(-82.96296296296296_dp)*Hr2(-1,0) + nf**2*z**5*(-47.407407407407405_dp)*Hr2(-1,0) + &
                        nf**2*z*(-35.55555555555556_dp)*Hr2(-1,0) + nf**2*z**4*(35.55555555555556_dp)*Hr2(-1,0) + &
                        nf**2*(47.407407407407405_dp)*Hr2(-1,0) + nf**2*z**3*(82.96296296296296_dp)*Hr2(-1,0) + &
                        nf*z**3*(1465.5473251028807_dp)*Hr2(-1,0) + nf*z**5*(6642.567901234568_dp)*Hr2(-1,0) + &
                        nf*z**4*(8108.115226337449_dp)*Hr2(-1,0) + z**3*(69258.04382287087_dp)*Hr2(-1,0) + &
                        (120533.95617712913_dp)*Hr2(-1,0) + z*(128021.95617712915_dp)*Hr2(-1,0) + &
                        z**2*(-179074.1022533653_dp)*Hr2(0,-1) + (-59691.36741778843_dp)*Hr2(0,-1) + &
                        z**3*(-59691.36741778843_dp)*Hr2(0,-1) + nf*z*(-3498.7950680141635_dp)*Hr2(0,-1) + &
                        nf*z**4*(-3498.7950680141635_dp)*Hr2(0,-1) + nf*z**2*(3498.7950680141635_dp)*Hr2(0,-1) + &
                        nf*z**3*(3498.7950680141635_dp)*Hr2(0,-1) + z**5*(59691.36741778843_dp)*Hr2(0,-1) + &
                        z*(59691.36741778844_dp)*Hr2(0,-1) + z**4*(179074.1022533653_dp)*Hr2(0,-1) + &
                        z**3*(-225457.25224900412_dp)*Hr2(0,0) + z*(-136126.7477509959_dp)*Hr2(0,0) + &
                        z**4*(-21662.02752304262_dp)*Hr2(0,0) + nf*z**2*(-12950.667908696938_dp)*Hr2(0,0) + &
                        nf*z**5*(-8055.440329218107_dp)*Hr2(0,0) + nf*z*(-7783.67391399428_dp)*Hr2(0,0) + &
                        (-1421.2230337568676_dp)*Hr2(0,0) + nf*(-1383.111111111111_dp)*Hr2(0,0) + &
                        nf**2*z**4*(-591.0123456790124_dp)*Hr2(0,0) + nf**2*z**3*(-77.03703703703704_dp)*Hr2(0,0) + &
                        nf**2*z**5*(-12.641975308641975_dp)*Hr2(0,0) + nf**2*z*(89.67901234567901_dp)*Hr2(0,0) + &
                        nf**2*z**2*(591.0123456790124_dp)*Hr2(0,0) + nf*z**4*(13373.779019808047_dp)*Hr2(0,0) + &
                        nf*z**3*(14879.114243212387_dp)*Hr2(0,0) + z**2*(40958.027523042634_dp)*Hr2(0,0) + &
                        z**5*(328294.7477509959_dp)*Hr2(0,0) + z*(-252065.654493647_dp)*Hr2(0,1) + &
                        z**3*(-70926.34550635303_dp)*Hr2(0,1) + z**4*(-41242.919658070125_dp)*Hr2(0,1) + &
                        nf*z**5*(-11230.683127572016_dp)*Hr2(0,1) + nf*z**2*(-9748.721545556473_dp)*Hr2(0,1) + &
                        nf*(-6642.172839506173_dp)*Hr2(0,1) + nf**2*z**4*(-223.20987654320987_dp)*Hr2(0,1) + &
                        nf**2*z**3*(-193.4485596707819_dp)*Hr2(0,1) + nf**2*z**5*(-30.551440329218106_dp)*Hr2(0,1) + &
                        nf**2*(-9.481481481481481_dp)*Hr2(0,1) + nf**3*z**3*(-4.7407407407407405_dp)*Hr2(0,1) + &
                        nf**3*z**4*(-4.7407407407407405_dp)*Hr2(0,1) + nf**3*z*(4.7407407407407405_dp)*Hr2(0,1) + &
                        nf**3*z**2*(4.7407407407407405_dp)*Hr2(0,1) + nf**2*z*(224._dp)*Hr2(0,1) + &
                        nf**2*z**2*(232.69135802469134_dp)*Hr2(0,1) + nf*z**3*(2328.803850083221_dp)*Hr2(0,1) + &
                        nf*z*(6021.879277488795_dp)*Hr2(0,1) + z**2*(8266.919658070125_dp)*Hr2(0,1) + &
                        (9854.287075858549_dp)*Hr2(0,1) + nf*z**4*(13510.894385062646_dp)*Hr2(0,1) + &
                        z**5*(299870.28707585856_dp)*Hr2(0,1) + z*(-163708.96348094096_dp)*Hr2(1,0) + &
                        z**4*(-163708.96348094096_dp)*Hr2(1,0) + z**2*(-14275.03651905909_dp)*Hr2(1,0) + &
                        z**3*(-14275.036519059071_dp)*Hr2(1,0) + nf*z**5*(-9770.798353909464_dp)*Hr2(1,0) + &
                        nf*(-9350.452674897118_dp)*Hr2(1,0) + nf*z**2*(-2572.773662551441_dp)*Hr2(1,0) + &
                        nf*z**3*(-1732.0823045267502_dp)*Hr2(1,0) + nf**2*(-5.267489711934156_dp)*Hr2(1,0) + &
                        nf**2*z**5*(-5.267489711934156_dp)*Hr2(1,0) + nf**2*z*(-3.9506172839506175_dp)*Hr2(1,0) + &
                        nf**2*z**4*(-3.9506172839506175_dp)*Hr2(1,0) + nf**2*z**2*(9.218106995884774_dp)*Hr2(1,0) + &
                        nf**2*z**3*(9.218106995884774_dp)*Hr2(1,0) + nf*z*(8622.880658436216_dp)*Hr2(1,0) + &
                        nf*z**4*(9043.22633744856_dp)*Hr2(1,0) + (171916.96348094096_dp)*Hr2(1,0) + &
                        z**5*(171916.96348094096_dp)*Hr2(1,0) + (z*(-42185.26516442312_dp) + z**2*(-31392._dp) + &
                        nf**2*((-87.440329218107_dp) + z**5*(-87.440329218107_dp) + z*(-65.58024691358025_dp) + &
                        z**4*(-65.58024691358025_dp) + z**2*(153.02057613168725_dp) + z**3*(153.02057613168725_dp)) + &
                        nf*(z**5*(-7929.942386831276_dp) + (-7089.251028806584_dp) + z**3*(-1302.2551440329216_dp) + &
                        z**2*(2776.3621399176945_dp) + z**4*(4312.88888888889_dp) + z*(9232.197530864198_dp)) + &
                        z**4*(13899.367417788439_dp) + (17492.63258221156_dp) + z**5*(17492.63258221156_dp) + &
                        z**3*(24692.632582211565_dp))*Hr2(1,1) + (-133056._dp)*Hr3(-1,-1,0) + z*(-112320._dp)*Hr3(-1,-1,0) + &
                        nf*z**5*(-2319.275720164609_dp)*Hr3(-1,-1,0) + z**3*(-1728.0000000000073_dp)*Hr3(-1,-1,0) + &
                        nf*z**2*(-890.7325102880659_dp)*Hr3(-1,-1,0) + nf*z**4*(-276.5432098765432_dp)*Hr3(-1,-1,0) + &
                        nf*z*(276.5432098765432_dp)*Hr3(-1,-1,0) + nf*z**3*(890.7325102880659_dp)*Hr3(-1,-1,0) + &
                        z**2*(1728.0000000000146_dp)*Hr3(-1,-1,0) + nf*(2319.275720164609_dp)*Hr3(-1,-1,0) + &
                        z**4*(112320._dp)*Hr3(-1,-1,0) + z**5*(133056._dp)*Hr3(-1,-1,0) + z**5*(-109296.00000000001_dp)*Hr3(-1,0,0) &
                        + z**4*(-92016._dp)*Hr3(-1,0,0) + z**2*(-3024.0000000000036_dp)*Hr3(-1,0,0) + &
                        nf*(-1864.8230452674898_dp)*Hr3(-1,0,0) + nf*z**3*(-887.4403292181069_dp)*Hr3(-1,0,0) + &
                        nf*z*(-113.3827160493829_dp)*Hr3(-1,0,0) + nf*z**4*(113.38271604938245_dp)*Hr3(-1,0,0) + &
                        nf*z**2*(887.4403292181075_dp)*Hr3(-1,0,0) + nf*z**5*(1864.8230452674898_dp)*Hr3(-1,0,0) + &
                        z**3*(3024.0000000000036_dp)*Hr3(-1,0,0) + z*(92016._dp)*Hr3(-1,0,0) + (109296.00000000001_dp)*Hr3(-1,0,0) + &
                        z**5*(-85536._dp)*Hr3(-1,0,1) + z**4*(-71712._dp)*Hr3(-1,0,1) + z**2*(-4320.000000000015_dp)*Hr3(-1,0,1) + &
                        nf*(-1410.3703703703702_dp)*Hr3(-1,0,1) + nf*z**3*(-884.1481481481483_dp)*Hr3(-1,0,1) + &
                        nf*z**4*(-49.77777777777783_dp)*Hr3(-1,0,1) + nf*z*(49.77777777777828_dp)*Hr3(-1,0,1) + &
                        nf*z**2*(884.1481481481483_dp)*Hr3(-1,0,1) + nf*z**5*(1410.3703703703702_dp)*Hr3(-1,0,1) + &
                        z**3*(4320.000000000007_dp)*Hr3(-1,0,1) + z*(71712._dp)*Hr3(-1,0,1) + (85536._dp)*Hr3(-1,0,1) + &
                        z**4*(-106272._dp)*Hr3(0,-1,0) + z**5*(-66528._dp)*Hr3(0,-1,0) + &
                        nf*z**3*(-1883.1275720164608_dp)*Hr3(0,-1,0) + nf*(-1600._dp)*Hr3(0,-1,0) + &
                        nf*z**2*(-992.3950617283951_dp)*Hr3(0,-1,0) + nf**2*z*(-71.11111111111111_dp)*Hr3(0,-1,0) + &
                        nf**2*z**4*(-71.11111111111111_dp)*Hr3(0,-1,0) + nf**2*z**2*(71.11111111111111_dp)*Hr3(0,-1,0) + &
                        nf**2*z**3*(71.11111111111111_dp)*Hr3(0,-1,0) + nf*z**5*(719.2757201646091_dp)*Hr3(0,-1,0) + &
                        nf*z*(1739.8518518518522_dp)*Hr3(0,-1,0) + nf*z**4*(2016.3950617283951_dp)*Hr3(0,-1,0) + &
                        z*(6047.999999999996_dp)*Hr3(0,-1,0) + z**2*(49247.999999999985_dp)*Hr3(0,-1,0) + &
                        z**3*(50976._dp)*Hr3(0,-1,0) + (66528._dp)*Hr3(0,-1,0) + z**4*(-117504._dp)*Hr3(0,0,0) + &
                        z*(-62208.00000000001_dp)*Hr3(0,0,0) + z**3*(-13824._dp)*Hr3(0,0,0) + &
                        nf*z**2*(-11922.172839506173_dp)*Hr3(0,0,0) + nf*z**3*(-1816.230452674897_dp)*Hr3(0,0,0) + &
                        nf*z**5*(-643.4238683127571_dp)*Hr3(0,0,0) + nf**2*z**4*(-203.85185185185185_dp)*Hr3(0,0,0) + &
                        nf**2*z**3*(-9.481481481481476_dp)*Hr3(0,0,0) + nf**2*z*(9.481481481481481_dp)*Hr3(0,0,0) + &
                        nf**2*z**2*(203.85185185185185_dp)*Hr3(0,0,0) + nf*z*(1883.6543209876545_dp)*Hr3(0,0,0) + &
                        nf*z**4*(11922.172839506173_dp)*Hr3(0,0,0) + z**5*(85536._dp)*Hr3(0,0,0) + z**2*(117504._dp)*Hr3(0,0,0) + &
                        z**4*(-158112.00000000003_dp)*Hr3(0,0,1) + z**3*(-110592._dp)*Hr3(0,0,1) + z*(-41472._dp)*Hr3(0,0,1) + &
                        nf*z**2*(-7909.925925925927_dp)*Hr3(0,0,1) + nf*z**5*(-1446.9794238683128_dp)*Hr3(0,0,1) + &
                        nf*z*(-1360.5925925925926_dp)*Hr3(0,0,1) + nf*(-1258.6666666666667_dp)*Hr3(0,0,1) + &
                        nf**2*z**3*(-148.54320987654324_dp)*Hr3(0,0,1) + nf**2*z**2*(-22.123456790123488_dp)*Hr3(0,0,1) + &
                        nf**2*z**4*(22.123456790123456_dp)*Hr3(0,0,1) + nf**2*z*(148.54320987654322_dp)*Hr3(0,0,1) + &
                        nf*z**3*(1655.5720164609052_dp)*Hr3(0,0,1) + nf*z**4*(8592.592592592591_dp)*Hr3(0,0,1) + &
                        (47520._dp)*Hr3(0,0,1) + z**2*(120096._dp)*Hr3(0,0,1) + z**5*(171072._dp)*Hr3(0,0,1) + &
                        z**4*(-156383.99999999997_dp)*Hr3(0,1,0) + z*(-68256._dp)*Hr3(0,1,0) + z**3*(-64800._dp)*Hr3(0,1,0) + &
                        nf*(-2403.5555555555557_dp)*Hr3(0,1,0) + nf*z*(-2153.8765432098767_dp)*Hr3(0,1,0) + &
                        nf*z**2*(-1847.3086419753083_dp)*Hr3(0,1,0) + nf*z**5*(-1712.4609053497943_dp)*Hr3(0,1,0) + &
                        nf**2*z*(-7.901234567901235_dp)*Hr3(0,1,0) + nf**2*z**2*(-7.901234567901235_dp)*Hr3(0,1,0) + &
                        nf**2*z**3*(7.901234567901235_dp)*Hr3(0,1,0) + nf**2*z**4*(7.901234567901235_dp)*Hr3(0,1,0) + &
                        nf*z**3*(2714.337448559671_dp)*Hr3(0,1,0) + nf*z**4*(3098.864197530864_dp)*Hr3(0,1,0) + &
                        z**2*(61344.00000000001_dp)*Hr3(0,1,0) + (114048.00000000001_dp)*Hr3(0,1,0) + z**5*(152064._dp)*Hr3(0,1,0) + &
                        z*(-124416.00000000003_dp)*Hr3(0,1,1) + z**4*(-100223.99999999999_dp)*Hr3(0,1,1) + &
                        z**3*(-27648._dp)*Hr3(0,1,1) + nf*(-1986.3703703703704_dp)*Hr3(0,1,1) + &
                        nf*z**5*(-1986.3703703703704_dp)*Hr3(0,1,1) + nf*z**2*(-1470.4197530864194_dp)*Hr3(0,1,1) + &
                        nf*z*(-1252.3456790123455_dp)*Hr3(0,1,1) + nf**2*z*(-131.1604938271605_dp)*Hr3(0,1,1) + &
                        nf**2*z**2*(-131.1604938271605_dp)*Hr3(0,1,1) + nf**2*z**3*(131.1604938271605_dp)*Hr3(0,1,1) + &
                        nf**2*z**4*(131.1604938271605_dp)*Hr3(0,1,1) + nf*z**3*(2086.716049382716_dp)*Hr3(0,1,1) + &
                        nf*z**4*(2304.7901234567903_dp)*Hr3(0,1,1) + z**2*(24191.999999999978_dp)*Hr3(0,1,1) + &
                        (95040._dp)*Hr3(0,1,1) + z**5*(171072._dp)*Hr3(0,1,1) + z*(-92016._dp)*Hr3(1,0,0) + &
                        z**4*(-92016._dp)*Hr3(1,0,0) + z**2*(-3024.0000000000073_dp)*Hr3(1,0,0) + &
                        z**3*(-3024.0000000000036_dp)*Hr3(1,0,0) + nf*(-1675.1934156378602_dp)*Hr3(1,0,0) + &
                        nf*z**5*(-1675.1934156378602_dp)*Hr3(1,0,0) + nf*z*(255.60493827160462_dp)*Hr3(1,0,0) + &
                        nf*z**4*(255.60493827160485_dp)*Hr3(1,0,0) + nf*z**2*(555.5884773662551_dp)*Hr3(1,0,0) + &
                        nf*z**3*(555.5884773662551_dp)*Hr3(1,0,0) + (109296.00000000001_dp)*Hr3(1,0,0) + &
                        z**5*(109296.00000000001_dp)*Hr3(1,0,0) + z**4*(-112320._dp)*Hr3(1,0,1) + &
                        z*(-112319.99999999999_dp)*Hr3(1,0,1) + nf*(-2129.6460905349795_dp)*Hr3(1,0,1) + &
                        nf*z**5*(-2129.6460905349795_dp)*Hr3(1,0,1) + z**3*(-1728._dp)*Hr3(1,0,1) + &
                        z**2*(-1727.9999999999964_dp)*Hr3(1,0,1) + nf*z*(418.7654320987649_dp)*Hr3(1,0,1) + &
                        nf*z**4*(418.7654320987656_dp)*Hr3(1,0,1) + nf*z**3*(558.8806584362137_dp)*Hr3(1,0,1) + &
                        nf*z**2*(558.8806584362139_dp)*Hr3(1,0,1) + (133056._dp)*Hr3(1,0,1) + z**5*(133056._dp)*Hr3(1,0,1) + &
                        z**4*(-112320._dp)*Hr3(1,1,0) + z*(-112319.99999999999_dp)*Hr3(1,1,0) + &
                        nf*(-2129.6460905349795_dp)*Hr3(1,1,0) + nf*z**5*(-2129.6460905349795_dp)*Hr3(1,1,0) + &
                        z**3*(-1728._dp)*Hr3(1,1,0) + z**2*(-1727.9999999999964_dp)*Hr3(1,1,0) + &
                        nf*z*(418.7654320987649_dp)*Hr3(1,1,0) + nf*z**4*(418.7654320987656_dp)*Hr3(1,1,0) + &
                        nf*z**3*(558.8806584362137_dp)*Hr3(1,1,0) + nf*z**2*(558.8806584362139_dp)*Hr3(1,1,0) + &
                        (133056._dp)*Hr3(1,1,0) + z**5*(133056._dp)*Hr3(1,1,0) + z**3*(-62208.00000000001_dp)*Hr4(-1,-1,-1,0) + &
                        z**4*(-62208.00000000001_dp)*Hr4(-1,-1,-1,0) + z**5*(-62208.00000000001_dp)*Hr4(-1,-1,-1,0) + &
                        z**2*(62207.99999999999_dp)*Hr4(-1,-1,-1,0) + (62208.00000000001_dp)*Hr4(-1,-1,-1,0) + &
                        z*(62208.00000000001_dp)*Hr4(-1,-1,-1,0) + (-51840._dp)*Hr4(-1,-1,0,0) + z*(-51840._dp)*Hr4(-1,-1,0,0) + &
                        z**2*(-51840._dp)*Hr4(-1,-1,0,0) + z**3*(51840._dp)*Hr4(-1,-1,0,0) + z**4*(51840._dp)*Hr4(-1,-1,0,0) + &
                        z**5*(51840._dp)*Hr4(-1,-1,0,0) + (-41472._dp)*Hr4(-1,-1,0,1) + z*(-41472._dp)*Hr4(-1,-1,0,1) + &
                        z**2*(-41472._dp)*Hr4(-1,-1,0,1) + z**3*(41472._dp)*Hr4(-1,-1,0,1) + z**4*(41472._dp)*Hr4(-1,-1,0,1) + &
                        z**5*(41472._dp)*Hr4(-1,-1,0,1) + (-31104.000000000004_dp)*Hr4(-1,0,-1,0) + &
                        z*(-31104.000000000004_dp)*Hr4(-1,0,-1,0) + z**2*(-31103.999999999996_dp)*Hr4(-1,0,-1,0) + &
                        z**3*(31104.000000000004_dp)*Hr4(-1,0,-1,0) + z**4*(31104.000000000004_dp)*Hr4(-1,0,-1,0) + &
                        z**5*(31104.000000000004_dp)*Hr4(-1,0,-1,0) + z**3*(-20736._dp)*Hr4(-1,0,0,0) + &
                        z**4*(-20736._dp)*Hr4(-1,0,0,0) + z**5*(-20736._dp)*Hr4(-1,0,0,0) + (20736._dp)*Hr4(-1,0,0,0) + &
                        z*(20736._dp)*Hr4(-1,0,0,0) + z**2*(20736._dp)*Hr4(-1,0,0,0) + z**3*(-31104.000000000004_dp)*Hr4(-1,0,0,1) + &
                        z**4*(-31104.000000000004_dp)*Hr4(-1,0,0,1) + z**5*(-31104.000000000004_dp)*Hr4(-1,0,0,1) + &
                        z**2*(31103.999999999996_dp)*Hr4(-1,0,0,1) + (31104.000000000004_dp)*Hr4(-1,0,0,1) + &
                        z*(31104.000000000004_dp)*Hr4(-1,0,0,1) + z**3*(-10368._dp)*Hr4(-1,0,1,0) + z**4*(-10368._dp)*Hr4(-1,0,1,0) &
                        + z**5*(-10368._dp)*Hr4(-1,0,1,0) + (10368._dp)*Hr4(-1,0,1,0) + z*(10368._dp)*Hr4(-1,0,1,0) + &
                        z**2*(10368._dp)*Hr4(-1,0,1,0) + z**3*(-20736._dp)*Hr4(-1,0,1,1) + z**4*(-20736._dp)*Hr4(-1,0,1,1) + &
                        z**5*(-20736._dp)*Hr4(-1,0,1,1) + (20736._dp)*Hr4(-1,0,1,1) + z*(20736._dp)*Hr4(-1,0,1,1) + &
                        z**2*(20736._dp)*Hr4(-1,0,1,1) + z**2*(-93312._dp)*Hr4(0,-1,-1,0) + (-31104.000000000004_dp)*Hr4(0,-1,-1,0) &
                        + z**3*(-31103.999999999993_dp)*Hr4(0,-1,-1,0) + nf*z*(-1750.9135802469134_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**4*(-1750.9135802469134_dp)*Hr4(0,-1,-1,0) + nf*z**2*(1750.9135802469134_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**3*(1750.9135802469134_dp)*Hr4(0,-1,-1,0) + z*(31103.999999999996_dp)*Hr4(0,-1,-1,0) + &
                        z**5*(31104.000000000004_dp)*Hr4(0,-1,-1,0) + z**4*(93312.00000000001_dp)*Hr4(0,-1,-1,0) + &
                        z**4*(-77760._dp)*Hr4(0,-1,0,0) + z*(-25920._dp)*Hr4(0,-1,0,0) + z**5*(-25920._dp)*Hr4(0,-1,0,0) + &
                        nf*z**2*(-1501.2345679012346_dp)*Hr4(0,-1,0,0) + nf*z**3*(-1501.2345679012346_dp)*Hr4(0,-1,0,0) + &
                        nf*z*(1501.2345679012346_dp)*Hr4(0,-1,0,0) + nf*z**4*(1501.2345679012346_dp)*Hr4(0,-1,0,0) + &
                        (25920._dp)*Hr4(0,-1,0,0) + z**3*(25920._dp)*Hr4(0,-1,0,0) + z**2*(77760._dp)*Hr4(0,-1,0,0) + &
                        z**4*(-62208.00000000001_dp)*Hr4(0,-1,0,1) + z*(-20736._dp)*Hr4(0,-1,0,1) + z**5*(-20736._dp)*Hr4(0,-1,0,1) &
                        + nf*z**2*(-1251.5555555555559_dp)*Hr4(0,-1,0,1) + nf*z**3*(-1251.5555555555559_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(1251.5555555555559_dp)*Hr4(0,-1,0,1) + nf*z**4*(1251.5555555555559_dp)*Hr4(0,-1,0,1) + &
                        z**3*(20735.999999999996_dp)*Hr4(0,-1,0,1) + (20736._dp)*Hr4(0,-1,0,1) + &
                        z**2*(62208.00000000001_dp)*Hr4(0,-1,0,1) + z**4*(-72576._dp)*Hr4(0,0,-1,0) + &
                        z**3*(-10368.000000000004_dp)*Hr4(0,0,-1,0) + z**5*(-10368._dp)*Hr4(0,0,-1,0) + &
                        nf*z**2*(-3179.4567901234573_dp)*Hr4(0,0,-1,0) + nf*z*(-1428.5432098765432_dp)*Hr4(0,0,-1,0) + &
                        nf*z**3*(1428.543209876543_dp)*Hr4(0,0,-1,0) + nf*z**4*(3179.456790123457_dp)*Hr4(0,0,-1,0) + &
                        (10368._dp)*Hr4(0,0,-1,0) + z*(10368._dp)*Hr4(0,0,-1,0) + z**2*(72576._dp)*Hr4(0,0,-1,0) + &
                        z**4*(-82944._dp)*Hr4(0,0,0,0) + z**3*(-10367.999999999998_dp)*Hr4(0,0,0,0) + &
                        nf*z**2*(-6757.135802469136_dp)*Hr4(0,0,0,0) + nf*z*(-2000.5925925925928_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**3*(-56.888888888888886_dp)*Hr4(0,0,0,0) + nf**2*z**4*(-56.888888888888886_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(56.888888888888886_dp)*Hr4(0,0,0,0) + nf**2*z**2*(56.888888888888886_dp)*Hr4(0,0,0,0) + &
                        nf*z**3*(2000.5925925925928_dp)*Hr4(0,0,0,0) + nf*z**4*(6757.135802469136_dp)*Hr4(0,0,0,0) + &
                        z*(10368._dp)*Hr4(0,0,0,0) + z**5*(10368._dp)*Hr4(0,0,0,0) + z**2*(82944._dp)*Hr4(0,0,0,0) + &
                        z**4*(-155520._dp)*Hr4(0,0,0,1) + z*(-31104.000000000004_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(-7212.246913580248_dp)*Hr4(0,0,0,1) + nf*z**3*(-3653.5308641975316_dp)*Hr4(0,0,0,1) + &
                        nf*z*(3653.5308641975307_dp)*Hr4(0,0,0,1) + nf*z**4*(7212.246913580248_dp)*Hr4(0,0,0,1) + &
                        (10368._dp)*Hr4(0,0,0,1) + z**3*(31104.000000000004_dp)*Hr4(0,0,0,1) + &
                        z**5*(31104.000000000004_dp)*Hr4(0,0,0,1) + z**2*(155520.00000000003_dp)*Hr4(0,0,0,1) + &
                        z**4*(-155520._dp)*Hr4(0,0,1,0) + z**3*(-20735.999999999996_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(-4282.469135802469_dp)*Hr4(0,0,1,0) + nf*z*(-527.8024691358024_dp)*Hr4(0,0,1,0) + &
                        nf*z**3*(527.8024691358027_dp)*Hr4(0,0,1,0) + nf*z**4*(4282.469135802469_dp)*Hr4(0,0,1,0) + &
                        z*(20736._dp)*Hr4(0,0,1,0) + (31104.000000000004_dp)*Hr4(0,0,1,0) + z**5*(41472._dp)*Hr4(0,0,1,0) + &
                        z**2*(155520._dp)*Hr4(0,0,1,0) + z**4*(-165888._dp)*Hr4(0,0,1,1) + z**3*(-62208.00000000001_dp)*Hr4(0,0,1,1) &
                        + nf*z*(-3343.8024691358023_dp)*Hr4(0,0,1,1) + nf*z**2*(-3343.8024691358023_dp)*Hr4(0,0,1,1) + &
                        nf*z**3*(3343.8024691358023_dp)*Hr4(0,0,1,1) + nf*z**4*(3343.8024691358023_dp)*Hr4(0,0,1,1) + &
                        (41472._dp)*Hr4(0,0,1,1) + z*(62208.00000000001_dp)*Hr4(0,0,1,1) + z**5*(62208.00000000001_dp)*Hr4(0,0,1,1) &
                        + z**2*(165888._dp)*Hr4(0,0,1,1) + z**4*(-88128._dp)*Hr4(0,1,0,0) + &
                        z**3*(-15551.99999999999_dp)*Hr4(0,1,0,0) + nf*z*(-1216.79012345679_dp)*Hr4(0,1,0,0) + &
                        nf*z**2*(-1216.79012345679_dp)*Hr4(0,1,0,0) + nf*z**3*(1216.79012345679_dp)*Hr4(0,1,0,0) + &
                        nf*z**4*(1216.79012345679_dp)*Hr4(0,1,0,0) + z*(15551.999999999993_dp)*Hr4(0,1,0,0) + &
                        (36288._dp)*Hr4(0,1,0,0) + z**5*(36288._dp)*Hr4(0,1,0,0) + z**2*(88128._dp)*Hr4(0,1,0,0) + &
                        z**4*(-124415.99999999997_dp)*Hr4(0,1,0,1) + nf*z*(-1466.469135802469_dp)*Hr4(0,1,0,1) + &
                        nf*z**2*(-1466.469135802469_dp)*Hr4(0,1,0,1) + z*(-7.275957614183426e-12_dp)*Hr4(0,1,0,1) + &
                        nf*z**3*(1466.469135802469_dp)*Hr4(0,1,0,1) + nf*z**4*(1466.469135802469_dp)*Hr4(0,1,0,1) + &
                        (62208.00000000001_dp)*Hr4(0,1,0,1) + z**5*(62208.00000000001_dp)*Hr4(0,1,0,1) + &
                        z**2*(124415.99999999999_dp)*Hr4(0,1,0,1) + z**4*(-124415.99999999997_dp)*Hr4(0,1,1,0) + &
                        nf*z*(-1466.469135802469_dp)*Hr4(0,1,1,0) + nf*z**2*(-1466.469135802469_dp)*Hr4(0,1,1,0) + &
                        z*(-7.275957614183426e-12_dp)*Hr4(0,1,1,0) + nf*z**3*(1466.469135802469_dp)*Hr4(0,1,1,0) + &
                        nf*z**4*(1466.469135802469_dp)*Hr4(0,1,1,0) + (62208.00000000001_dp)*Hr4(0,1,1,0) + &
                        z**5*(62208.00000000001_dp)*Hr4(0,1,1,0) + z**2*(124415.99999999999_dp)*Hr4(0,1,1,0) + &
                        z**4*(-62208.00000000001_dp)*Hr4(0,1,1,1) + z*(-62207.99999999997_dp)*Hr4(0,1,1,1) + &
                        z**3*(62207.999999999985_dp)*Hr4(0,1,1,1) + (62208.00000000001_dp)*Hr4(0,1,1,1) + &
                        z**2*(62208.00000000001_dp)*Hr4(0,1,1,1) + z**5*(62208.00000000001_dp)*Hr4(0,1,1,1) + &
                        z*(-20736.000000000007_dp)*Hr4(1,0,0,0) + z**4*(-20736._dp)*Hr4(1,0,0,0) + (20736._dp)*Hr4(1,0,0,0) + &
                        z**2*(20736._dp)*Hr4(1,0,0,0) + z**5*(20736._dp)*Hr4(1,0,0,0) + z**3*(20736.000000000007_dp)*Hr4(1,0,0,0) + &
                        z*(-51840._dp)*Hr4(1,0,0,1) + z**4*(-51840._dp)*Hr4(1,0,0,1) + (51840._dp)*Hr4(1,0,0,1) + &
                        z**2*(51840._dp)*Hr4(1,0,0,1) + z**3*(51840._dp)*Hr4(1,0,0,1) + z**5*(51840._dp)*Hr4(1,0,0,1) + &
                        z**4*(-62208.00000000001_dp)*Hr4(1,0,1,0) + z*(-62207.99999999997_dp)*Hr4(1,0,1,0) + &
                        z**3*(62207.999999999985_dp)*Hr4(1,0,1,0) + (62208.00000000001_dp)*Hr4(1,0,1,0) + &
                        z**2*(62208.00000000001_dp)*Hr4(1,0,1,0) + z**5*(62208.00000000001_dp)*Hr4(1,0,1,0) + &
                        z**4*(-62208.00000000001_dp)*Hr4(1,0,1,1) + z*(-62207.99999999997_dp)*Hr4(1,0,1,1) + &
                        z**3*(62207.999999999985_dp)*Hr4(1,0,1,1) + (62208.00000000001_dp)*Hr4(1,0,1,1) + &
                        z**2*(62208.00000000001_dp)*Hr4(1,0,1,1) + z**5*(62208.00000000001_dp)*Hr4(1,0,1,1) + &
                        z*(-51840._dp)*Hr4(1,1,0,0) + z**4*(-51840._dp)*Hr4(1,1,0,0) + (51840._dp)*Hr4(1,1,0,0) + &
                        z**2*(51840._dp)*Hr4(1,1,0,0) + z**3*(51840._dp)*Hr4(1,1,0,0) + z**5*(51840._dp)*Hr4(1,1,0,0) + &
                        z**4*(-62208.00000000001_dp)*Hr4(1,1,0,1) + z*(-62207.99999999997_dp)*Hr4(1,1,0,1) + &
                        z**3*(62207.999999999985_dp)*Hr4(1,1,0,1) + (62208.00000000001_dp)*Hr4(1,1,0,1) + &
                        z**2*(62208.00000000001_dp)*Hr4(1,1,0,1) + z**5*(62208.00000000001_dp)*Hr4(1,1,0,1) + &
                        z**4*(-62208.00000000001_dp)*Hr4(1,1,1,0) + z*(-62207.99999999997_dp)*Hr4(1,1,1,0) + &
                        z**3*(62207.999999999985_dp)*Hr4(1,1,1,0) + (62208.00000000001_dp)*Hr4(1,1,1,0) + &
                        z**2*(62208.00000000001_dp)*Hr4(1,1,1,0) + z**5*(62208.00000000001_dp)*Hr4(1,1,1,0))/(z*(z + &
                        (1._dp))*((-1._dp) + z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-21729.993224184578_dp) + nf**2*(-164.7668539356091_dp) + nf**3*(1.876543209876543_dp) + &
                            nf*(5596.3775213278295_dp)

                        else
                            Ibar_select = (-21729.993224184578_dp) + nf**2*(-164.7668539356091_dp) + nf**3*(1.876543209876543_dp) + &
                            nf*(5596.3775213278295_dp)

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-7048.873757119791_dp) + nf**2*(68.32853724541276_dp) + (65664.46725135812_dp)
                        else
                            Ibar_select = nf*(-7048.873757119791_dp) + nf**2*(68.32853724541276_dp) + (65664.46725135812_dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (-5904.569014926325_dp) + nf*(-2770.7769662431324_dp) + nf**2*(106.66666666666667_dp)
                        else
                            Ibar_select = (-5904.569014926325_dp) + nf*(-2770.7769662431324_dp) + nf**2*(106.66666666666667_dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-28042.31629110578_dp) + nf*(2880._dp)
                        else
                            Ibar_select = (-28042.31629110578_dp) + nf*(2880._dp)
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
                        Ibar_select = (z*(-102586.79373986671_dp) + z**3*(-29741.026243435335_dp) + nf*z**2*(-5059.707602809469_dp) &
                        + nf*z**4*(-3472.664905217326_dp) + nf*(-3079.2887385807658_dp) + nf**2*z*(-108.4415237142116_dp) + &
                        nf**2*z**3*(-32.12912100458128_dp) + nf**3*z*(-0.8230452674897119_dp) + nf**3*z**3*(-0.5267489711934156_dp) &
                        + nf**3*(0.03292181069958847_dp) + nf**3*z**4*(0.03292181069958847_dp) + nf**3*z**2*(1.2839506172839505_dp) &
                        + nf**2*(44.12071330589849_dp) + nf**2*z**4*(44.12071330589849_dp) + nf**2*z**2*(52.32921810699589_dp) + &
                        nf*z**3*(4563.055400716628_dp) + nf*z*(7048.605845890932_dp) + (40865.05791436658_dp) + &
                        z**4*(40865.05791436658_dp) + z**2*(54752.01281188843_dp) + z**3*(-47856.86291978024_dp)*Hr1(0) + &
                        z*(-36966.215729945056_dp)*Hr1(0) + nf*(-1378.6666666666667_dp)*Hr1(0) + &
                        nf*z**4*(-883.2592592592592_dp)*Hr1(0) + nf*z**2*(-174.88490305533003_dp)*Hr1(0) + &
                        nf**2*z**2*(-21.827160493827158_dp)*Hr1(0) + nf**2*z**3*(-15.539094650205762_dp)*Hr1(0) + &
                        nf**2*z**4*(-0.5596707818930042_dp)*Hr1(0) + nf**3*z*(-0.3950617283950617_dp)*Hr1(0) + &
                        nf**3*z**3*(0.3950617283950617_dp)*Hr1(0) + nf**2*(13.662551440329215_dp)*Hr1(0) + &
                        nf**2*z*(21.596707818930042_dp)*Hr1(0) + nf*z*(44.60644254140357_dp)*Hr1(0) + &
                        nf*z**3*(2480.2043864398524_dp)*Hr1(0) + z**2*(19665.970784752764_dp)*Hr1(0) + (26157.10786497253_dp)*Hr1(0) &
                        + z**4*(32589.107864972535_dp)*Hr1(0) + z**3*(-107562.64718983517_dp)*Hr1(1) + &
                        z*(-89055.97078475278_dp)*Hr1(1) + nf*z**2*(-4701.168724279835_dp)*Hr1(1) + &
                        nf*(-2261.925925925926_dp)*Hr1(1) + nf*z**4*(-2261.925925925926_dp)*Hr1(1) + &
                        nf**2*z**2*(-43.654320987654316_dp)*Hr1(1) + nf**2*z**3*(6.057613168724281_dp)*Hr1(1) + &
                        nf**2*z*(11.390946502057613_dp)*Hr1(1) + nf**2*(13.102880658436213_dp)*Hr1(1) + &
                        nf**2*z**4*(13.102880658436213_dp)*Hr1(1) + nf*z*(4524.510288065843_dp)*Hr1(1) + &
                        nf*z**3*(4700.510288065843_dp)*Hr1(1) + (53061.32359491759_dp)*Hr1(1) + z**4*(53061.32359491759_dp)*Hr1(1) + &
                        z**2*(90495.97078475278_dp)*Hr1(1) + z**3*(-34560._dp)*Hr2(0,0) + z*(-14688._dp)*Hr2(0,0) + &
                        nf*z**2*(-1955.9506172839504_dp)*Hr2(0,0) + nf*(-224._dp)*Hr2(0,0) + &
                        nf*z**4*(-104.42798353909464_dp)*Hr2(0,0) + nf**2*z**3*(-41.08641975308642_dp)*Hr2(0,0) + &
                        nf**2*z*(-13.432098765432098_dp)*Hr2(0,0) + nf**2*z**2*(54.51851851851852_dp)*Hr2(0,0) + &
                        nf*z*(703.2098765432099_dp)*Hr2(0,0) + nf*z**3*(1485.1687242798355_dp)*Hr2(0,0) + &
                        (7919.999999999998_dp)*Hr2(0,0) + z**4*(7920._dp)*Hr2(0,0) + z**2*(34992._dp)*Hr2(0,0) + &
                        z**3*(-63936._dp)*Hr2(0,1) + z*(-24191.999999999996_dp)*Hr2(0,1) + nf*z**2*(-602.4691358024692_dp)*Hr2(0,1) &
                        + nf*(-581.9259259259259_dp)*Hr2(0,1) + nf*z**4*(-342.7818930041152_dp)*Hr2(0,1) + &
                        nf*z*(-162.37037037037038_dp)*Hr2(0,1) + nf**2*z**3*(-27.65432098765432_dp)*Hr2(0,1) + &
                        nf**2*z*(27.65432098765432_dp)*Hr2(0,1) + nf*z**3*(1401.5473251028807_dp)*Hr2(0,1) + (23760._dp)*Hr2(0,1) + &
                        z**4*(23760._dp)*Hr2(0,1) + z**2*(45360._dp)*Hr2(0,1) + z*(-44064._dp)*Hr2(1,0) + z**3*(-44064._dp)*Hr2(1,0) &
                        + nf*z**2*(-602.4691358024692_dp)*Hr2(1,0) + nf*(-462.3539094650206_dp)*Hr2(1,0) + &
                        nf*z**4*(-462.3539094650206_dp)*Hr2(1,0) + nf*z**3*(619.5884773662551_dp)*Hr2(1,0) + &
                        nf*z*(619.5884773662553_dp)*Hr2(1,0) + (23760._dp)*Hr2(1,0) + z**4*(23760._dp)*Hr2(1,0) + &
                        z**2*(45360._dp)*Hr2(1,0) + (z*(-97632._dp) + z**3*(-88128._dp) + nf*(z**2*(-1204.9382716049383_dp) + &
                        (-924.7078189300412_dp) + z**4*(-924.7078189300412_dp) + z**3*(1239.1769547325102_dp) + &
                        z*(1815.1769547325102_dp)) + (47520._dp) + z**4*(47520._dp) + z**2*(90720._dp))*Hr2(1,1) + &
                        (z**3*(-13824._dp) + nf**2*(z**3*(-4.7407407407407405_dp) + z*(4.7407407407407405_dp)) + &
                        nf*z*(z*(-661.3333333333334_dp) + (-166.71604938271605_dp) + z**2*(828.0493827160494_dp)) + (864._dp) + &
                        z**4*(864._dp) + z**2*(12960._dp))*Hr3(0,0,0) + z**3*(-27648._dp)*Hr3(0,0,1) + z*(-6912._dp)*Hr3(0,0,1) + &
                        nf*z**2*(-1322.6666666666665_dp)*Hr3(0,0,1) + nf*z*(399.8024691358025_dp)*Hr3(0,0,1) + &
                        nf*z**3*(922.8641975308642_dp)*Hr3(0,0,1) + (3456._dp)*Hr3(0,0,1) + z**4*(3456._dp)*Hr3(0,0,1) + &
                        z**2*(31104._dp)*Hr3(0,0,1) + z**3*(-20736._dp)*Hr3(0,1,0) + nf*z*(-261.5308641975309_dp)*Hr3(0,1,0) + &
                        nf*z**3*(261.5308641975309_dp)*Hr3(0,1,0) + (5184._dp)*Hr3(0,1,0) + z**4*(5184._dp)*Hr3(0,1,0) + &
                        z**2*(15552._dp)*Hr3(0,1,0) + z**3*(-41472._dp)*Hr3(0,1,1) + nf*z*(-523.0617283950618_dp)*Hr3(0,1,1) + &
                        nf*z**3*(523.0617283950618_dp)*Hr3(0,1,1) + (10368._dp)*Hr3(0,1,1) + z**4*(10368._dp)*Hr3(0,1,1) + &
                        z**2*(31104._dp)*Hr3(0,1,1) + z*(-6912._dp)*Hr3(1,0,0) + z**3*(-6912._dp)*Hr3(1,0,0) + (3456._dp)*Hr3(1,0,0) &
                        + z**4*(3456._dp)*Hr3(1,0,0) + z**2*(10368._dp)*Hr3(1,0,0) + z*(-20736._dp)*Hr3(1,0,1) + &
                        z**3*(-20736._dp)*Hr3(1,0,1) + (10368._dp)*Hr3(1,0,1) + z**4*(10368._dp)*Hr3(1,0,1) + &
                        z**2*(31104._dp)*Hr3(1,0,1) + z*(-20736._dp)*Hr3(1,1,0) + z**3*(-20736._dp)*Hr3(1,1,0) + &
                        (10368._dp)*Hr3(1,1,0) + z**4*(10368._dp)*Hr3(1,1,0) + z**2*(31104._dp)*Hr3(1,1,0) + &
                        z*(-62208._dp)*Hr3(1,1,1) + z**3*(-41472._dp)*Hr3(1,1,1) + (20736._dp)*Hr3(1,1,1) + &
                        z**4*(20736._dp)*Hr3(1,1,1) + z**2*(62208._dp)*Hr3(1,1,1))/(z*((-1._dp) + z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-1211.203225788147_dp) + nf**2*(4.386490844928604_dp) + nf*(86.04072752399819_dp)

                        else
                            Ibar_select = (z*(-1211.2032257881465_dp) + nf**2*((-4.386490844928604_dp) + z*(4.386490844928604_dp)) + &
                            nf*((-86.04072752399816_dp) + z*(86.04072752399816_dp)) + (1211.2032257881465_dp))/((-1._dp) + &
                            z*(1._dp))

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-715.7410112522891_dp) + nf**3*(-0.2962962962962963_dp) + &
                            nf**2*(14.666666666666666_dp) + (17456.344000301888_dp)

                        else
                            Ibar_select = ((-17456.344000301888_dp) + nf**3*(z*(-0.2962962962962963_dp) + (0.2962962962962963_dp)) + &
                            nf**2*((-14.666666666666668_dp) + z*(14.666666666666668_dp)) + nf*(z*(-715.7410112522892_dp) + &
                            (715.7410112522892_dp)) + z*(17456.344000301888_dp))/((-1._dp) + z*(1._dp))

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (-18506.67640508241_dp) + nf**2*(-5.333333333333333_dp) + nf*(176._dp)
                        else
                            Ibar_select = (-18506.67640508241_dp) + nf**2*(-5.333333333333333_dp) + nf*(176._dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-4752._dp) + nf*(288._dp)
                        else
                            Ibar_select = (z*(-4752._dp) + nf*((-288._dp) + z*(288._dp)) + (4752._dp))/((-1._dp) + z*(1._dp))

                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (3456._dp)
                        else
                            Ibar_select = (3456._dp)
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
                        Ibar_select = (z**5*(-3.595299775629447e7_dp) + (-4.238413990861842e6_dp) + nf*z*(-1.098909341318743e6_dp) + &
                        nf*z**3*(-287452.4479656426_dp) + nf*z**2*(-128789.92166623939_dp) + nf**2*z**4*(-15156.185003047027_dp) + &
                        nf**2*z*(-10564.926706413065_dp) + nf**2*z**3*(-5422.272663333591_dp) + nf**2*z**2*(-3666.813808814596_dp) + &
                        nf**3*z**5*(-644.4604807918503_dp) + nf**3*(-558.1252825066533_dp) + nf**3*z*(-180.67148884355728_dp) + &
                        nf**3*z**4*(-172.75802262412248_dp) + nf**4*z**3*(-11.870957599476425_dp) + &
                        nf**4*z**4*(-2.0383101371993324_dp) + nf**4*z**2*(-1.4733496707567704_dp) + nf**4*(3.511659807956104_dp) + &
                        nf**4*z**5*(3.511659807956104_dp) + nf**4*z*(8.359297791520321_dp) + nf**3*z**2*(730.8833051307759_dp) + &
                        nf**3*z**3*(825.1319696354077_dp) + nf**2*z**5*(15756.404444340013_dp) + nf**2*(19669.24687168598_dp) + &
                        nf*(34376.329789389216_dp) + nf*z**4*(47088.68459618879_dp) + z**4*(97217.26131592109_dp) + &
                        nf*z**5*(1.3577470095856714e6_dp) + z**3*(3.3063145545618325e6_dp) + z**2*(3.8472467740044985e6_dp) + &
                        z*(3.3439012623164866e7_dp) + ((-1.1415247292583121e7_dp) + z*(-1.0908859078165937e7_dp) + &
                        z**2*(-278572.7525335288_dp) + nf**2*((-2847.2116378143983_dp) + z**3*(-2666.5532000528947_dp) + &
                        z**4*(-661.547804464788_dp) + z*(661.5478044647889_dp) + z**2*(2666.5532000528938_dp) + &
                        z**5*(2847.2116378143983_dp)) + z**3*(278572.7525335285_dp) + nf*(z**5*(-394526.3015504059_dp) + &
                        z**4*(-357210.2594829997_dp) + z**3*(-48637.08512853526_dp) + z**2*(48637.085128535226_dp) + &
                        z*(357210.2594829997_dp) + (394526.3015504059_dp)) + z**4*(1.090885907816594e7_dp) + &
                        z**5*(1.1415247292583121e7_dp))*Hr1(-1) + z*(-3.1074898525837217e6_dp)*Hr1(0) + &
                        z**5*(-3.083689676024273e6_dp)*Hr1(0) + (-941493.0765290941_dp)*Hr1(0) + &
                        nf*z**5*(-734518.1102242261_dp)*Hr1(0) + nf*(-138460.77605292882_dp)*Hr1(0) + &
                        nf*z**3*(-54875.57345538099_dp)*Hr1(0) + nf*z**2*(-49551.32707640482_dp)*Hr1(0) + &
                        nf**2*z**4*(-18852.770884740392_dp)*Hr1(0) + nf**2*z**3*(-13194.353130819776_dp)*Hr1(0) + &
                        nf**2*z*(-6678.220197785403_dp)*Hr1(0) + nf**3*z**2*(-788.7884969655465_dp)*Hr1(0) + &
                        nf**3*z*(-747.2365377277679_dp)*Hr1(0) + nf**3*(-139.20804755372654_dp)*Hr1(0) + &
                        nf**4*z**4*(-8.427983539094651_dp)*Hr1(0) + nf**4*z**3*(-6.320987654320988_dp)*Hr1(0) + &
                        nf**4*z**5*(-2.1069958847736627_dp)*Hr1(0) + nf**4*z*(8.427983539094651_dp)*Hr1(0) + &
                        nf**4*z**2*(8.427983539094651_dp)*Hr1(0) + nf**3*z**5*(53.34796524919983_dp)*Hr1(0) + &
                        nf**3*z**3*(670.1848687748645_dp)*Hr1(0) + nf**3*z**4*(904.2928408155694_dp)*Hr1(0) + &
                        nf**2*z**2*(9374.707029353356_dp)*Hr1(0) + nf**2*(10485.309657504396_dp)*Hr1(0) + &
                        nf**2*z**5*(20037.61288849624_dp)*Hr1(0) + nf*z**4*(182365.94429697137_dp)*Hr1(0) + &
                        z**4*(314216.2271188573_dp)*Hr1(0) + z**2*(663725.0175254318_dp)*Hr1(0) + nf*z*(848929.1087565884_dp)*Hr1(0) &
                        + z**3*(5.774175666714225e6_dp)*Hr1(0) + (-8.951406301137175e6_dp)*Hr1(1) + &
                        z**5*(-8.951406301137175e6_dp)*Hr1(1) + z**3*(-446898.5305006087_dp)*Hr1(1) + &
                        nf*z**5*(-206086.18101545802_dp)*Hr1(1) + nf*(-96349.00038653224_dp)*Hr1(1) + &
                        nf*z**2*(-87850.71737396951_dp)*Hr1(1) + nf*z*(-15722.365327015541_dp)*Hr1(1) + &
                        nf**2*z**3*(-9969.733545853494_dp)*Hr1(1) + nf**2*z**4*(-6069.061901887954_dp)*Hr1(1) + &
                        nf**2*z**2*(-5268.296240223345_dp)*Hr1(1) + nf**2*z*(-4515.037490633664_dp)*Hr1(1) + &
                        nf**3*z**4*(-147.57750342935532_dp)*Hr1(1) + nf**3*z*(-111.23182441700962_dp)*Hr1(1) + &
                        nf**4*(-2.1069958847736627_dp)*Hr1(1) + nf**4*z**5*(-2.1069958847736627_dp)*Hr1(1) + &
                        nf**4*z*(-1.5802469135802466_dp)*Hr1(1) + nf**4*z**4*(-1.5802469135802466_dp)*Hr1(1) + &
                        nf**4*z**2*(3.6872427983539096_dp)*Hr1(1) + nf**4*z**3*(3.6872427983539096_dp)*Hr1(1) + &
                        nf**3*z**5*(36.17009602194788_dp)*Hr1(1) + nf**3*(47.231824417009605_dp)*Hr1(1) + &
                        nf**3*z**3*(75.0617283950617_dp)*Hr1(1) + nf**3*z**2*(100.34567901234568_dp)*Hr1(1) + &
                        nf**2*(11337.358142111298_dp)*Hr1(1) + nf**2*z**5*(14484.77103648716_dp)*Hr1(1) + &
                        nf*z**4*(167582.48313122356_dp)*Hr1(1) + nf*z**3*(205191.31171319535_dp)*Hr1(1) + &
                        z**2*(503229.73381381895_dp)*Hr1(1) + z**4*(8.890683848013204e6_dp)*Hr1(1) + &
                        z*(9.84081211232763e6_dp)*Hr1(1) + z**5*(-6.635567564776269e6_dp)*Hr2(-1,-1) + &
                        z**4*(-5.953380508572971e6_dp)*Hr2(-1,-1) + z**3*(-2.201351699454841e6_dp)*Hr2(-1,-1) + &
                        nf*(-111601.56941314253_dp)*Hr2(-1,-1) + nf*z*(-55578.6802483161_dp)*Hr2(-1,-1) + &
                        nf*z**2*(-23565.600725558175_dp)*Hr2(-1,-1) + nf*z**3*(23565.60072555815_dp)*Hr2(-1,-1) + &
                        nf*z**4*(55578.68024831612_dp)*Hr2(-1,-1) + nf*z**5*(111601.56941314253_dp)*Hr2(-1,-1) + &
                        z**2*(2.2013516994548417e6_dp)*Hr2(-1,-1) + z*(5.953380508572971e6_dp)*Hr2(-1,-1) + &
                        (6.635567564776269e6_dp)*Hr2(-1,-1) + z**4*(-1.5039286221698965e6_dp)*Hr2(-1,0) + &
                        z**5*(-1.0482856828275891e6_dp)*Hr2(-1,0) + z**2*(-949395.2114474123_dp)*Hr2(-1,0) + &
                        nf*z*(-260026.42100138124_dp)*Hr2(-1,0) + nf*(-227943.25161159842_dp)*Hr2(-1,0) + &
                        nf*z**3*(-19080.85982546444_dp)*Hr2(-1,0) + nf**2*z**5*(-5370.732510288066_dp)*Hr2(-1,0) + &
                        nf**2*z**4*(-4357.0041152263375_dp)*Hr2(-1,0) + nf**2*z**2*(-1013.7283950617284_dp)*Hr2(-1,0) + &
                        nf**3*z**3*(-49.77777777777777_dp)*Hr2(-1,0) + nf**3*(-28.444444444444443_dp)*Hr2(-1,0) + &
                        nf**3*z**4*(-21.333333333333336_dp)*Hr2(-1,0) + nf**3*z*(21.333333333333336_dp)*Hr2(-1,0) + &
                        nf**3*z**5*(28.444444444444443_dp)*Hr2(-1,0) + nf**3*z**2*(49.77777777777777_dp)*Hr2(-1,0) + &
                        nf**2*z**3*(1013.7283950617284_dp)*Hr2(-1,0) + nf**2*z*(4357.0041152263375_dp)*Hr2(-1,0) + &
                        nf**2*(5370.732510288066_dp)*Hr2(-1,0) + nf*z**2*(19080.859825464453_dp)*Hr2(-1,0) + &
                        nf*z**5*(227943.25161159842_dp)*Hr2(-1,0) + nf*z**4*(260026.42100138124_dp)*Hr2(-1,0) + &
                        z**3*(949395.2114474123_dp)*Hr2(-1,0) + (1.048285682827589e6_dp)*Hr2(-1,0) + &
                        z*(1.5039286221698965e6_dp)*Hr2(-1,0) + z**2*(-5.007495189690503e6_dp)*Hr2(0,-1) + &
                        (-3.3177837823881344e6_dp)*Hr2(0,-1) + z**3*(-2.8061434902356616e6_dp)*Hr2(0,-1) + &
                        nf*z**4*(-182100.51960476278_dp)*Hr2(0,-1) + nf*z*(-126521.83935644673_dp)*Hr2(0,-1) + &
                        nf*z**5*(-43213.9026639712_dp)*Hr2(0,-1) + nf**2*z**2*(-3007.50809338216_dp)*Hr2(0,-1) + &
                        nf**2*z**3*(-3007.50809338216_dp)*Hr2(0,-1) + nf**2*z*(3007.50809338216_dp)*Hr2(0,-1) + &
                        nf**2*z**4*(3007.50809338216_dp)*Hr2(0,-1) + nf*(68387.66674917135_dp)*Hr2(0,-1) + &
                        nf*z**3*(129941.4970752257_dp)*Hr2(0,-1) + nf*z**2*(153507.0978007838_dp)*Hr2(0,-1) + &
                        z*(930129.0856765972_dp)*Hr2(0,-1) + z**5*(3.3177837823881344e6_dp)*Hr2(0,-1) + &
                        z**4*(6.883509594249569e6_dp)*Hr2(0,-1) + z*(-4.961245011573579e6_dp)*Hr2(0,0) + &
                        z**3*(-3.8186736771247042e6_dp)*Hr2(0,0) + nf*z**5*(-492890.2703616504_dp)*Hr2(0,0) + &
                        nf*z**2*(-295929.39825709944_dp)*Hr2(0,0) + z**4*(-248234.1570332076_dp)*Hr2(0,0) + &
                        (-101338.0520533657_dp)*Hr2(0,0) + nf*(-50265.74101125229_dp)*Hr2(0,0) + &
                        nf**2*z**4*(-23024.725992287647_dp)*Hr2(0,0) + nf**2*z**3*(-4769.173980108522_dp)*Hr2(0,0) + &
                        nf**3*z**2*(-431.5716218931046_dp)*Hr2(0,0) + nf**3*z*(-141.06956427993586_dp)*Hr2(0,0) + &
                        nf**3*z**5*(12.641975308641973_dp)*Hr2(0,0) + nf**3*z**3*(128.42758897129386_dp)*Hr2(0,0) + &
                        nf**3*z**4*(431.5716218931046_dp)*Hr2(0,0) + nf**2*z*(889.6092796055491_dp)*Hr2(0,0) + &
                        nf**2*(1568.9218106995886_dp)*Hr2(0,0) + nf**2*z**5*(4519.564700502972_dp)*Hr2(0,0) + &
                        nf**2*z**2*(22095.804181588057_dp)*Hr2(0,0) + nf*z*(142158.5080448354_dp)*Hr2(0,0) + &
                        nf*z**4*(321792.36230210867_dp)*Hr2(0,0) + nf*z**3*(360438.3381607368_dp)*Hr2(0,0) + &
                        z**2*(654497.436805254_dp)*Hr2(0,0) + z**5*(8.247855223666153e6_dp)*Hr2(0,0) + &
                        z*(-5.4084992185497675e6_dp)*Hr2(0,1) + (-1.7410781846450255e6_dp)*Hr2(0,1) + &
                        z**3*(-958159.8948395886_dp)*Hr2(0,1) + nf*z**5*(-562902.891138973_dp)*Hr2(0,1) + &
                        nf*z**2*(-312285.934382299_dp)*Hr2(0,1) + nf*(-187700.6868618072_dp)*Hr2(0,1) + &
                        nf**2*z**4*(-11778.741868457164_dp)*Hr2(0,1) + nf**2*z**3*(-9018.635787067135_dp)*Hr2(0,1) + &
                        nf**3*z**2*(-110.26611796982166_dp)*Hr2(0,1) + nf**3*z*(-76.72976680384087_dp)*Hr2(0,1) + &
                        nf**4*z*(-3.160493827160493_dp)*Hr2(0,1) + nf**4*z**2*(-3.160493827160493_dp)*Hr2(0,1) + &
                        nf**4*z**3*(3.160493827160493_dp)*Hr2(0,1) + nf**4*z**4*(3.160493827160493_dp)*Hr2(0,1) + &
                        nf**3*(21.421124828532232_dp)*Hr2(0,1) + nf**3*z**5*(34.06310013717421_dp)*Hr2(0,1) + &
                        nf**3*z**3*(42.666666666666664_dp)*Hr2(0,1) + nf**3*z**4*(88.84499314128944_dp)*Hr2(0,1) + &
                        nf**2*(4726.108824874257_dp)*Hr2(0,1) + nf**2*z**5*(5397.128486511203_dp)*Hr2(0,1) + &
                        nf**2*z*(5541.507300555932_dp)*Hr2(0,1) + nf**2*z**2*(8972.633043582908_dp)*Hr2(0,1) + &
                        nf*z**3*(215953.20329117915_dp)*Hr2(0,1) + nf*z*(324905.3861643117_dp)*Hr2(0,1) + &
                        nf*z**4*(477942.31956062425_dp)*Hr2(0,1) + z**2*(527791.5953776948_dp)*Hr2(0,1) + &
                        z**4*(822727.5294093239_dp)*Hr2(0,1) + z**5*(5.976100053531348e6_dp)*Hr2(0,1) + &
                        z**4*(-2.862249016149364e6_dp)*Hr2(1,0) + z*(-2.8622490161493635e6_dp)*Hr2(1,0) + &
                        nf*z**5*(-412491.52845684346_dp)*Hr2(1,0) + nf*(-391250.37619346904_dp)*Hr2(1,0) + &
                        nf*z**2*(-40770.88256040506_dp)*Hr2(1,0) + z**3*(-21023.859063617885_dp)*Hr2(1,0) + &
                        z**2*(-21023.859063617652_dp)*Hr2(1,0) + nf**2*z**4*(-2797.3004115226336_dp)*Hr2(1,0) + &
                        nf**2*z*(-2551.308641975308_dp)*Hr2(1,0) + nf**2*z**3*(-1477.9990855052583_dp)*Hr2(1,0) + &
                        nf**2*z**2*(-986.0155464106083_dp)*Hr2(1,0) + nf**3*z**2*(-30.112482853223593_dp)*Hr2(1,0) + &
                        nf**3*z**3*(-30.112482853223593_dp)*Hr2(1,0) + nf**3*z*(12.905349794238683_dp)*Hr2(1,0) + &
                        nf**3*z**4*(12.905349794238683_dp)*Hr2(1,0) + nf**3*(17.207133058984912_dp)*Hr2(1,0) + &
                        nf**3*z**5*(17.207133058984912_dp)*Hr2(1,0) + nf*z**3*(1711.4219663438816_dp)*Hr2(1,0) + &
                        nf**2*(5703.315957933242_dp)*Hr2(1,0) + nf**2*z**5*(5949.307727480567_dp)*Hr2(1,0) + &
                        nf*z*(377366.0205369626_dp)*Hr2(1,0) + nf*z**4*(398607.172800337_dp)*Hr2(1,0) + &
                        (2.7550928116426333e6_dp)*Hr2(1,0) + z**5*(2.7550928116426333e6_dp)*Hr2(1,0) + ((-1.3984208099810563e6_dp) + &
                        z**5*(-1.398420809981056e6_dp) + nf**3*(z**2*(-146.87517146776406_dp) + z**3*(-146.87517146776406_dp) + &
                        z*(62.946502057613166_dp) + z**4*(62.946502057613166_dp) + (83.9286694101509_dp) + &
                        z**5*(83.9286694101509_dp)) + nf**2*(z**2*(-4298.827617741198_dp) + z**3*(-1442.7946959304982_dp) + &
                        z*(-1124.5212620027437_dp) + (2075.332418838592_dp) + z**4*(2223.4951989026063_dp) + &
                        z**5*(2567.315957933242_dp)) + z**2*(147886.39861849174_dp) + nf*(z**5*(-329427.4885231026_dp) + &
                        (-286945.1839963536_dp) + z**3*(-27792.920725554948_dp) + z**2*(-450.22079175886785_dp) + &
                        z**4*(287395.40478811244_dp) + z*(357220.40924865747_dp)) + z**3*(592080.3294135854_dp) + &
                        z*(1.1054507038944783e6_dp) + z**4*(1.5496446346895723e6_dp))*Hr2(1,1) + &
                        (-2.0465611686098895e6_dp)*Hr3(-1,-1,-1) + z*(-2.0465611686098895e6_dp)*Hr3(-1,-1,-1) + &
                        z**2*(-2.0465611686098895e6_dp)*Hr3(-1,-1,-1) + z**3*(2.0465611686098895e6_dp)*Hr3(-1,-1,-1) + &
                        z**4*(2.0465611686098895e6_dp)*Hr3(-1,-1,-1) + z**5*(2.0465611686098895e6_dp)*Hr3(-1,-1,-1) + &
                        (-2.5602391819730774e6_dp)*Hr3(-1,-1,0) + z*(-2.553327181973078e6_dp)*Hr3(-1,-1,0) + &
                        z**3*(-1.300112818026922e6_dp)*Hr3(-1,-1,0) + nf*z**4*(-159624.51577503426_dp)*Hr3(-1,-1,0) + &
                        nf*z**5*(-159194.07407407407_dp)*Hr3(-1,-1,0) + nf*z**3*(-17326.44170096021_dp)*Hr3(-1,-1,0) + &
                        nf**2*(-1508.6090534979423_dp)*Hr3(-1,-1,0) + nf**2*z**3*(-1232.0658436213994_dp)*Hr3(-1,-1,0) + &
                        nf**2*z**4*(-235.4567901234568_dp)*Hr3(-1,-1,0) + nf**2*z*(235.4567901234568_dp)*Hr3(-1,-1,0) + &
                        nf**2*z**2*(1232.065843621399_dp)*Hr3(-1,-1,0) + nf**2*z**5*(1508.6090534979423_dp)*Hr3(-1,-1,0) + &
                        nf*z**2*(17326.441700960208_dp)*Hr3(-1,-1,0) + nf*(159194.07407407407_dp)*Hr3(-1,-1,0) + &
                        nf*z*(159624.5157750343_dp)*Hr3(-1,-1,0) + z**2*(1.3001128180269229e6_dp)*Hr3(-1,-1,0) + &
                        z**4*(2.5533271819730774e6_dp)*Hr3(-1,-1,0) + z**5*(2.5602391819730774e6_dp)*Hr3(-1,-1,0) + &
                        z**3*(-1.0232805843049447e6_dp)*Hr3(-1,0,-1) + z**4*(-1.0232805843049447e6_dp)*Hr3(-1,0,-1) + &
                        z**5*(-1.0232805843049447e6_dp)*Hr3(-1,0,-1) + (1.0232805843049447e6_dp)*Hr3(-1,0,-1) + &
                        z*(1.0232805843049447e6_dp)*Hr3(-1,0,-1) + z**2*(1.0232805843049447e6_dp)*Hr3(-1,0,-1) + &
                        z**4*(-2.881738355037363e6_dp)*Hr3(-1,0,0) + z**5*(-2.878282355037363e6_dp)*Hr3(-1,0,0) + &
                        z**2*(-437749.64496263734_dp)*Hr3(-1,0,0) + nf*z*(-136476.31275720167_dp)*Hr3(-1,0,0) + &
                        nf*(-134481.64609053498_dp)*Hr3(-1,0,0) + nf*z**2*(-14666.66666666668_dp)*Hr3(-1,0,0) + &
                        nf**2*z**5*(-1242.6008230452676_dp)*Hr3(-1,0,0) + nf**2*z**2*(-1118.551440329218_dp)*Hr3(-1,0,0) + &
                        nf**2*z*(-259.9506172839506_dp)*Hr3(-1,0,0) + nf**2*z**4*(259.9506172839506_dp)*Hr3(-1,0,0) + &
                        nf**2*z**3*(1118.5514403292182_dp)*Hr3(-1,0,0) + nf**2*(1242.6008230452676_dp)*Hr3(-1,0,0) + &
                        nf*z**3*(14666.666666666679_dp)*Hr3(-1,0,0) + nf*z**5*(134481.64609053498_dp)*Hr3(-1,0,0) + &
                        nf*z**4*(136476.31275720167_dp)*Hr3(-1,0,0) + z**3*(437749.64496263704_dp)*Hr3(-1,0,0) + &
                        (2.878282355037363e6_dp)*Hr3(-1,0,0) + z*(2.8817383550373632e6_dp)*Hr3(-1,0,0) + &
                        z**4*(-2.664399883139011e6_dp)*Hr3(-1,0,1) + z**5*(-2.650575883139011e6_dp)*Hr3(-1,0,1) + &
                        z**2*(-121136.11686098891_dp)*Hr3(-1,0,1) + nf*z*(-113328.10973936899_dp)*Hr3(-1,0,1) + &
                        nf*(-109769.21810699589_dp)*Hr3(-1,0,1) + nf*z**2*(-12006.891632373112_dp)*Hr3(-1,0,1) + &
                        nf**2*z**2*(-1005.0370370370371_dp)*Hr3(-1,0,1) + nf**2*z**5*(-976.5925925925926_dp)*Hr3(-1,0,1) + &
                        nf**2*z*(-284.44444444444434_dp)*Hr3(-1,0,1) + nf**2*z**4*(284.44444444444434_dp)*Hr3(-1,0,1) + &
                        nf**2*(976.5925925925926_dp)*Hr3(-1,0,1) + nf**2*z**3*(1005.0370370370368_dp)*Hr3(-1,0,1) + &
                        nf*z**3*(12006.8916323731_dp)*Hr3(-1,0,1) + nf*z**5*(109769.21810699589_dp)*Hr3(-1,0,1) + &
                        nf*z**4*(113328.10973936899_dp)*Hr3(-1,0,1) + z**3*(121136.11686098887_dp)*Hr3(-1,0,1) + &
                        (2.650575883139011e6_dp)*Hr3(-1,0,1) + z*(2.664399883139011e6_dp)*Hr3(-1,0,1) + &
                        z**4*(-3.069841752914834e6_dp)*Hr3(0,-1,-1) + z*(-1.0232805843049447e6_dp)*Hr3(0,-1,-1) + &
                        z**5*(-1.0232805843049447e6_dp)*Hr3(0,-1,-1) + nf*z**2*(-48019.619284136934_dp)*Hr3(0,-1,-1) + &
                        nf*z**3*(-48019.619284136934_dp)*Hr3(0,-1,-1) + nf*z*(48019.619284136934_dp)*Hr3(0,-1,-1) + &
                        nf*z**4*(48019.619284136934_dp)*Hr3(0,-1,-1) + (1.0232805843049447e6_dp)*Hr3(0,-1,-1) + &
                        z**3*(1.0232805843049448e6_dp)*Hr3(0,-1,-1) + z**2*(3.069841752914834e6_dp)*Hr3(0,-1,-1) + &
                        z**5*(-1.4344875909865387e6_dp)*Hr3(0,-1,0) + z**4*(-1.092838772959616e6_dp)*Hr3(0,-1,0) + &
                        z**2*(-679513.227040384_dp)*Hr3(0,-1,0) + nf*(-94298.07407407407_dp)*Hr3(0,-1,0) + &
                        nf*z**2*(-89688.98159692141_dp)*Hr3(0,-1,0) + nf*z**3*(-72362.53989596119_dp)*Hr3(0,-1,0) + &
                        nf**2*z*(-2849.1851851851848_dp)*Hr3(0,-1,0) + nf**2*z**4*(-2613.7283950617284_dp)*Hr3(0,-1,0) + &
                        nf**2*z**5*(-298.13991769547323_dp)*Hr3(0,-1,0) + nf**3*z**2*(-42.666666666666664_dp)*Hr3(0,-1,0) + &
                        nf**3*z**3*(-42.666666666666664_dp)*Hr3(0,-1,0) + nf**3*z*(42.666666666666664_dp)*Hr3(0,-1,0) + &
                        nf**3*z**4*(42.666666666666664_dp)*Hr3(0,-1,0) + nf**2*(1210.469135802469_dp)*Hr3(0,-1,0) + &
                        nf**2*z**2*(1659.2592592592594_dp)*Hr3(0,-1,0) + nf**2*z**3*(2891.325102880659_dp)*Hr3(0,-1,0) + &
                        nf*z*(15914.539895961207_dp)*Hr3(0,-1,0) + nf*z**5*(64896._dp)*Hr3(0,-1,0) + &
                        nf*z**4*(175539.0556709955_dp)*Hr3(0,-1,0) + z**3*(620599.5909865387_dp)*Hr3(0,-1,0) + &
                        (1.1257515909865387e6_dp)*Hr3(0,-1,0) + z*(1.4604884090134613e6_dp)*Hr3(0,-1,0) + &
                        z**2*(-2.3876546967115374e6_dp)*Hr3(0,0,-1) + (-341093.52810164815_dp)*Hr3(0,0,-1) + &
                        z*(-341093.5281016481_dp)*Hr3(0,0,-1) + nf*z**4*(-88345.00870102132_dp)*Hr3(0,0,-1) + &
                        nf*z**3*(-40325.38941688438_dp)*Hr3(0,0,-1) + nf*z*(40325.38941688438_dp)*Hr3(0,0,-1) + &
                        nf*z**2*(88345.00870102132_dp)*Hr3(0,0,-1) + z**3*(341093.52810164815_dp)*Hr3(0,0,-1) + &
                        z**5*(341093.52810164815_dp)*Hr3(0,0,-1) + z**4*(2.3876546967115374e6_dp)*Hr3(0,0,-1) + &
                        z*(-2.7778424915825897e6_dp)*Hr3(0,0,0) + z**4*(-1.3665720285930296e6_dp)*Hr3(0,0,0) + &
                        z**3*(-393613.50841741054_dp)*Hr3(0,0,0) + nf*z**2*(-292433.0914834052_dp)*Hr3(0,0,0) + &
                        nf*z**3*(-88879.77450423062_dp)*Hr3(0,0,0) + nf*z**5*(-85204.45541838135_dp)*Hr3(0,0,0) + &
                        nf**2*z**4*(-14696.083047396518_dp)*Hr3(0,0,0) + nf*(-5248._dp)*Hr3(0,0,0) + &
                        (-4263.669101270602_dp)*Hr3(0,0,0) + nf**2*z*(-3526.066228983893_dp)*Hr3(0,0,0) + &
                        nf**3*z**2*(-144.32921810699588_dp)*Hr3(0,0,0) + nf**3*z*(-71.63786008230453_dp)*Hr3(0,0,0) + &
                        nf**3*z**3*(71.63786008230453_dp)*Hr3(0,0,0) + nf**3*z**4*(144.32921810699588_dp)*Hr3(0,0,0) + &
                        nf**2*z**5*(250.38134430727015_dp)*Hr3(0,0,0) + nf**2*z**3*(3531.6848846766225_dp)*Hr3(0,0,0) + &
                        nf**2*z**2*(14696.083047396518_dp)*Hr3(0,0,0) + nf*z*(161796.22992261196_dp)*Hr3(0,0,0) + &
                        nf*z**4*(293841.0914834052_dp)*Hr3(0,0,0) + z**2*(1.4437560285930294e6_dp)*Hr3(0,0,0) + &
                        z**5*(3.041197508417411e6_dp)*Hr3(0,0,0) + z*(-3.5785403797982135e6_dp)*Hr3(0,0,1) + &
                        z**3*(-2.532819620201787e6_dp)*Hr3(0,0,1) + z**4*(-1.5616981010089302e6_dp)*Hr3(0,0,1) + &
                        nf*z**2*(-331881.80348196946_dp)*Hr3(0,0,1) + nf*z**5*(-183603.62139917695_dp)*Hr3(0,0,1) + &
                        nf*(-80564.14814814815_dp)*Hr3(0,0,1) + nf**2*z**4*(-9486.046639231823_dp)*Hr3(0,0,1) + &
                        nf**2*z**3*(-3331.775034293552_dp)*Hr3(0,0,1) + nf**3*z*(-85.86008230452674_dp)*Hr3(0,0,1) + &
                        nf**3*z**4*(-45.300411522633745_dp)*Hr3(0,0,1) + nf**3*z**2*(45.30041152263375_dp)*Hr3(0,0,1) + &
                        nf**3*z**3*(85.86008230452674_dp)*Hr3(0,0,1) + nf**2*z**5*(520.4279835390946_dp)*Hr3(0,0,1) + &
                        nf**2*(733.2345679012345_dp)*Hr3(0,0,1) + nf*z*(770.8479830362559_dp)*Hr3(0,0,1) + &
                        nf**2*z*(3323.347050754458_dp)*Hr3(0,0,1) + nf**2*z**2*(9008.812071330589_dp)*Hr3(0,0,1) + &
                        nf*z**3*(150576.77341614064_dp)*Hr3(0,0,1) + (370139.56177129154_dp)*Hr3(0,0,1) + &
                        nf*z**4*(388637.95163011755_dp)*Hr3(0,0,1) + z**2*(802530.1010089298_dp)*Hr3(0,0,1) + &
                        z**5*(5.587371444910303e6_dp)*Hr3(0,0,1) + z**4*(-3.527266257238393e6_dp)*Hr3(0,1,0) + &
                        z*(-3.330683484278778e6_dp)*Hr3(0,1,0) + z**3*(-1.5359405157212226e6_dp)*Hr3(0,1,0) + &
                        nf*z**2*(-186177.71108006244_dp)*Hr3(0,1,0) + nf*z**5*(-182011.78600823044_dp)*Hr3(0,1,0) + &
                        nf*(-169382.45267489713_dp)*Hr3(0,1,0) + nf**2*z**4*(-3395.511659807956_dp)*Hr3(0,1,0) + &
                        nf**2*z**3*(-3384.537722908093_dp)*Hr3(0,1,0) + nf**3*z**3*(-25.810699588477366_dp)*Hr3(0,1,0) + &
                        nf**3*z**4*(-25.810699588477366_dp)*Hr3(0,1,0) + nf**3*z*(25.810699588477366_dp)*Hr3(0,1,0) + &
                        nf**3*z**2*(25.810699588477366_dp)*Hr3(0,1,0) + nf**2*z**5*(659.8408779149521_dp)*Hr3(0,1,0) + &
                        nf**2*(1330.9190672153634_dp)*Hr3(0,1,0) + nf**2*z**2*(2576.5925925925926_dp)*Hr3(0,1,0) + &
                        nf**2*z*(3236.696844993141_dp)*Hr3(0,1,0) + nf*z*(28774.1901545054_dp)*Hr3(0,1,0) + &
                        nf*z**3*(113301.59585372506_dp)*Hr3(0,1,0) + nf*z**4*(315624.1637549596_dp)*Hr3(0,1,0) + &
                        z**2*(918562.2572383931_dp)*Hr3(0,1,0) + (2.4011028707585856e6_dp)*Hr3(0,1,0) + &
                        z**5*(4.659022870758586e6_dp)*Hr3(0,1,0) + z*(-4.28613979549327e6_dp)*Hr3(0,1,1) + &
                        z**4*(-2.1401067437443696e6_dp)*Hr3(0,1,1) + z**3*(-1.2849322045067304e6_dp)*Hr3(0,1,1) + &
                        nf*z**5*(-236955.9176954733_dp)*Hr3(0,1,1) + nf*z**2*(-172792.038851038_dp)*Hr3(0,1,1) + &
                        nf*(-163860.01646090535_dp)*Hr3(0,1,1) + nf**2*z**4*(-2391.0891632373114_dp)*Hr3(0,1,1) + &
                        nf**2*z**3*(-1917.5418381344307_dp)*Hr3(0,1,1) + nf**3*z**3*(-125.89300411522633_dp)*Hr3(0,1,1) + &
                        nf**3*z**4*(-125.89300411522633_dp)*Hr3(0,1,1) + nf**3*z*(125.89300411522633_dp)*Hr3(0,1,1) + &
                        nf**3*z**2*(125.89300411522633_dp)*Hr3(0,1,1) + nf**2*(562.5679012345679_dp)*Hr3(0,1,1) + &
                        nf**2*z**5*(644.7407407407408_dp)*Hr3(0,1,1) + nf**2*z*(1784.80109739369_dp)*Hr3(0,1,1) + &
                        nf**2*z**2*(2340.5212620027432_dp)*Hr3(0,1,1) + nf*z**3*(56816.92773992689_dp)*Hr3(0,1,1) + &
                        nf*z*(117162.98995554636_dp)*Hr3(0,1,1) + nf*z**4*(273676.0553119433_dp)*Hr3(0,1,1) + &
                        z**2*(467402.74374436896_dp)*Hr3(0,1,1) + (1.4592032696188195e6_dp)*Hr3(0,1,1) + &
                        z**5*(5.3575712696188195e6_dp)*Hr3(0,1,1) + z*(-3.1016703842526106e6_dp)*Hr3(1,0,0) + &
                        z**4*(-3.10167038425261e6_dp)*Hr3(1,0,0) + nf*z**5*(-145447.1550068587_dp)*Hr3(1,0,0) + &
                        nf*(-144118.6941015089_dp)*Hr3(1,0,0) + z**2*(-140633.61574739066_dp)*Hr3(1,0,0) + &
                        z**3*(-140633.6157473899_dp)*Hr3(1,0,0) + nf*z**2*(-22428.181069958864_dp)*Hr3(1,0,0) + &
                        nf*z**3*(-19771.25925925924_dp)*Hr3(1,0,0) + nf**2*z**2*(-539.0397805212617_dp)*Hr3(1,0,0) + &
                        nf**2*z**3*(-539.0397805212617_dp)*Hr3(1,0,0) + nf**2*z**4*(11.58847736625512_dp)*Hr3(1,0,0) + &
                        nf**2*z*(11.588477366255574_dp)*Hr3(1,0,0) + nf**2*(911.4513031550067_dp)*Hr3(1,0,0) + &
                        nf**2*z**5*(911.4513031550067_dp)*Hr3(1,0,0) + nf*z*(137186.41426611796_dp)*Hr3(1,0,0) + &
                        nf*z**4*(138514.87517146778_dp)*Hr3(1,0,0) + (3.1609983842526106e6_dp)*Hr3(1,0,0) + &
                        z**5*(3.1609983842526106e6_dp)*Hr3(1,0,0) + z*(-3.21312326961882e6_dp)*Hr3(1,0,1) + &
                        z**4*(-3.21312326961882e6_dp)*Hr3(1,0,1) + z**2*(-408764.730381181_dp)*Hr3(1,0,1) + &
                        z**3*(-408764.7303811806_dp)*Hr3(1,0,1) + nf*z**5*(-208506.20576131687_dp)*Hr3(1,0,1) + &
                        nf*(-204520.8230452675_dp)*Hr3(1,0,1) + nf*z**2*(-45486.74897119346_dp)*Hr3(1,0,1) + &
                        nf*z**3*(-37515.98353909464_dp)*Hr3(1,0,1) + nf**2*z**4*(-327.37448559670815_dp)*Hr3(1,0,1) + &
                        nf**2*z*(-327.374485596707_dp)*Hr3(1,0,1) + nf**2*z**3*(81.20713305898471_dp)*Hr3(1,0,1) + &
                        nf**2*z**2*(81.20713305898511_dp)*Hr3(1,0,1) + nf**2*(758.1673525377229_dp)*Hr3(1,0,1) + &
                        nf**2*z**5*(758.1673525377229_dp)*Hr3(1,0,1) + nf*z*(183046.18930041147_dp)*Hr3(1,0,1) + &
                        nf*z**4*(187031.5720164609_dp)*Hr3(1,0,1) + (3.4083872696188195e6_dp)*Hr3(1,0,1) + &
                        z**5*(3.4083872696188195e6_dp)*Hr3(1,0,1) + z**4*(-3.417779386479808e6_dp)*Hr3(1,1,0) + &
                        z*(-3.417779386479807e6_dp)*Hr3(1,1,0) + nf*z**5*(-208506.20576131687_dp)*Hr3(1,1,0) + &
                        nf*(-204520.8230452675_dp)*Hr3(1,1,0) + z**3*(-204108.61352019233_dp)*Hr3(1,1,0) + &
                        z**2*(-204108.6135201922_dp)*Hr3(1,1,0) + nf*z**2*(-45486.74897119346_dp)*Hr3(1,1,0) + &
                        nf*z**3*(-37515.98353909464_dp)*Hr3(1,1,0) + nf**2*z**4*(-327.37448559670815_dp)*Hr3(1,1,0) + &
                        nf**2*z*(-327.374485596707_dp)*Hr3(1,1,0) + nf**2*z**3*(81.20713305898471_dp)*Hr3(1,1,0) + &
                        nf**2*z**2*(81.20713305898511_dp)*Hr3(1,1,0) + nf**2*(758.1673525377229_dp)*Hr3(1,1,0) + &
                        nf**2*z**5*(758.1673525377229_dp)*Hr3(1,1,0) + nf*z*(183046.18930041147_dp)*Hr3(1,1,0) + &
                        nf*z**4*(187031.5720164609_dp)*Hr3(1,1,0) + (3.613043386479809e6_dp)*Hr3(1,1,0) + &
                        z**5*(3.613043386479809e6_dp)*Hr3(1,1,0) + z*(-1.0231189482510997e6_dp)*Hr3(1,1,1) + &
                        z**2*(-376703.9999999999_dp)*Hr3(1,1,1) + nf*z**5*(-115039.8683127572_dp)*Hr3(1,1,1) + &
                        nf*(-107069.10288065845_dp)*Hr3(1,1,1) + z**4*(-91655.4741255498_dp)*Hr3(1,1,1) + &
                        nf*z**3*(-45254.84773662552_dp)*Hr3(1,1,1) + nf**2*(-1257.8765432098764_dp)*Hr3(1,1,1) + &
                        nf**2*z**5*(-1257.8765432098764_dp)*Hr3(1,1,1) + nf**2*z*(-943.4074074074074_dp)*Hr3(1,1,1) + &
                        nf**2*z**4*(-943.4074074074074_dp)*Hr3(1,1,1) + nf**2*z**2*(2201.283950617284_dp)*Hr3(1,1,1) + &
                        nf**2*z**3*(2201.283950617284_dp)*Hr3(1,1,1) + nf*z**2*(30963.621399176965_dp)*Hr3(1,1,1) + &
                        nf*z**4*(76105.48148148147_dp)*Hr3(1,1,1) + nf*z*(160294.71604938273_dp)*Hr3(1,1,1) + &
                        (468359.4741255498_dp)*Hr3(1,1,1) + z**5*(468359.4741255498_dp)*Hr3(1,1,1) + &
                        z**3*(554759.4741255499_dp)*Hr3(1,1,1) + z**5*(-2.509056e6_dp)*Hr4(-1,-1,-1,0) + &
                        z**4*(-2.17728e6_dp)*Hr4(-1,-1,-1,0) + z**3*(-352511.9999999998_dp)*Hr4(-1,-1,-1,0) + &
                        nf*(-56749.12482853225_dp)*Hr4(-1,-1,-1,0) + nf*z*(-30014.156378600816_dp)*Hr4(-1,-1,-1,0) + &
                        nf*z**2*(-14737.031550068583_dp)*Hr4(-1,-1,-1,0) + nf*z**3*(14737.031550068576_dp)*Hr4(-1,-1,-1,0) + &
                        nf*z**4*(30014.156378600823_dp)*Hr4(-1,-1,-1,0) + nf*z**5*(56749.12482853225_dp)*Hr4(-1,-1,-1,0) + &
                        z**2*(352511.9999999999_dp)*Hr4(-1,-1,-1,0) + z*(2.1772799999999995e6_dp)*Hr4(-1,-1,-1,0) + &
                        (2.509056e6_dp)*Hr4(-1,-1,-1,0) + (-2.166912e6_dp)*Hr4(-1,-1,0,0) + &
                        z*(-1.8766079999999998e6_dp)*Hr4(-1,-1,0,0) + z**2*(-279936.0000000001_dp)*Hr4(-1,-1,0,0) + &
                        nf*z**5*(-48110.09053497943_dp)*Hr4(-1,-1,0,0) + nf*z**4*(-24397.432098765432_dp)*Hr4(-1,-1,0,0) + &
                        nf*z**3*(-10847.341563786014_dp)*Hr4(-1,-1,0,0) + nf*z**2*(10847.341563786009_dp)*Hr4(-1,-1,0,0) + &
                        nf*z*(24397.432098765436_dp)*Hr4(-1,-1,0,0) + nf*(48110.09053497943_dp)*Hr4(-1,-1,0,0) + &
                        z**3*(279936.0000000001_dp)*Hr4(-1,-1,0,0) + z**4*(1.8766079999999998e6_dp)*Hr4(-1,-1,0,0) + &
                        z**5*(2.166912e6_dp)*Hr4(-1,-1,0,0) + (-1.8247680000000002e6_dp)*Hr4(-1,-1,0,1) + &
                        z*(-1.575936e6_dp)*Hr4(-1,-1,0,1) + z**2*(-207359.99999999997_dp)*Hr4(-1,-1,0,1) + &
                        nf*z**5*(-39471.05624142661_dp)*Hr4(-1,-1,0,1) + nf*z**4*(-18780.707818930037_dp)*Hr4(-1,-1,0,1) + &
                        nf*z**3*(-6957.651577503427_dp)*Hr4(-1,-1,0,1) + nf*z**2*(6957.651577503431_dp)*Hr4(-1,-1,0,1) + &
                        nf*z*(18780.707818930045_dp)*Hr4(-1,-1,0,1) + nf*(39471.05624142661_dp)*Hr4(-1,-1,0,1) + &
                        z**3*(207359.99999999983_dp)*Hr4(-1,-1,0,1) + z**4*(1.575936e6_dp)*Hr4(-1,-1,0,1) + &
                        z**5*(1.8247680000000002e6_dp)*Hr4(-1,-1,0,1) + (-1.254528e6_dp)*Hr4(-1,0,-1,0) + &
                        z*(-1.0886399999999998e6_dp)*Hr4(-1,0,-1,0) + z**2*(-176255.99999999994_dp)*Hr4(-1,0,-1,0) + &
                        nf*z**5*(-28374.562414266125_dp)*Hr4(-1,0,-1,0) + nf*z**4*(-15007.078189300411_dp)*Hr4(-1,0,-1,0) + &
                        nf*z**3*(-7368.515775034288_dp)*Hr4(-1,0,-1,0) + nf*z**2*(7368.515775034291_dp)*Hr4(-1,0,-1,0) + &
                        nf*z*(15007.078189300408_dp)*Hr4(-1,0,-1,0) + nf*(28374.562414266125_dp)*Hr4(-1,0,-1,0) + &
                        z**3*(176255.9999999999_dp)*Hr4(-1,0,-1,0) + z**4*(1.08864e6_dp)*Hr4(-1,0,-1,0) + &
                        z**5*(1.254528e6_dp)*Hr4(-1,0,-1,0) + z**5*(-912384.0000000001_dp)*Hr4(-1,0,0,0) + &
                        z**4*(-787968._dp)*Hr4(-1,0,0,0) + z**3*(-103679.99999999991_dp)*Hr4(-1,0,0,0) + &
                        nf*(-19735.528120713305_dp)*Hr4(-1,0,0,0) + nf*z*(-9390.353909465022_dp)*Hr4(-1,0,0,0) + &
                        nf*z**2*(-3478.8257887517157_dp)*Hr4(-1,0,0,0) + nf*z**3*(3478.8257887517134_dp)*Hr4(-1,0,0,0) + &
                        nf*z**4*(9390.353909465019_dp)*Hr4(-1,0,0,0) + nf*z**5*(19735.528120713305_dp)*Hr4(-1,0,0,0) + &
                        z**2*(103679.99999999999_dp)*Hr4(-1,0,0,0) + z*(787968._dp)*Hr4(-1,0,0,0) + &
                        (912384.0000000001_dp)*Hr4(-1,0,0,0) + z**5*(-1.482624e6_dp)*Hr4(-1,0,0,1) + &
                        z**4*(-1.275264e6_dp)*Hr4(-1,0,0,1) + z**3*(-134783.9999999999_dp)*Hr4(-1,0,0,1) + &
                        nf*(-30832.0219478738_dp)*Hr4(-1,0,0,1) + nf*z*(-13163.983539094654_dp)*Hr4(-1,0,0,1) + &
                        nf*z**2*(-3067.96159122085_dp)*Hr4(-1,0,0,1) + nf*z**3*(3067.961591220857_dp)*Hr4(-1,0,0,1) + &
                        nf*z**4*(13163.98353909465_dp)*Hr4(-1,0,0,1) + nf*z**5*(30832.0219478738_dp)*Hr4(-1,0,0,1) + &
                        z**2*(134783.99999999994_dp)*Hr4(-1,0,0,1) + z*(1.275264e6_dp)*Hr4(-1,0,0,1) + (1.482624e6_dp)*Hr4(-1,0,0,1) &
                        + z**5*(-570240._dp)*Hr4(-1,0,1,0) + z**4*(-487296._dp)*Hr4(-1,0,1,0) + &
                        z**3*(-31103.999999999964_dp)*Hr4(-1,0,1,0) + nf*(-11096.493827160493_dp)*Hr4(-1,0,1,0) + &
                        nf*z*(-3773.6296296296296_dp)*Hr4(-1,0,1,0) + nf*z**3*(-410.8641975308634_dp)*Hr4(-1,0,1,0) + &
                        nf*z**2*(410.86419753086284_dp)*Hr4(-1,0,1,0) + nf*z**4*(3773.6296296296296_dp)*Hr4(-1,0,1,0) + &
                        nf*z**5*(11096.493827160493_dp)*Hr4(-1,0,1,0) + z**2*(31104._dp)*Hr4(-1,0,1,0) + &
                        z*(487296._dp)*Hr4(-1,0,1,0) + (570240._dp)*Hr4(-1,0,1,0) + z**5*(-1.14048e6_dp)*Hr4(-1,0,1,1) + &
                        z**4*(-974592._dp)*Hr4(-1,0,1,1) + z**3*(-62207.99999999993_dp)*Hr4(-1,0,1,1) + &
                        nf*(-22192.987654320987_dp)*Hr4(-1,0,1,1) + nf*z*(-7547.259259259259_dp)*Hr4(-1,0,1,1) + &
                        nf*z**3*(-821.7283950617268_dp)*Hr4(-1,0,1,1) + nf*z**2*(821.7283950617257_dp)*Hr4(-1,0,1,1) + &
                        nf*z**4*(7547.259259259259_dp)*Hr4(-1,0,1,1) + nf*z**5*(22192.987654320987_dp)*Hr4(-1,0,1,1) + &
                        z**2*(62208._dp)*Hr4(-1,0,1,1) + z*(974592._dp)*Hr4(-1,0,1,1) + (1.14048e6_dp)*Hr4(-1,0,1,1) + &
                        z**2*(-1.358208e6_dp)*Hr4(0,-1,-1,0) + (-1.254528e6_dp)*Hr4(0,-1,-1,0) + &
                        z**3*(-1.005696e6_dp)*Hr4(0,-1,-1,0) + nf*z**4*(-75334.05761316873_dp)*Hr4(0,-1,-1,0) + &
                        nf*z*(-45319.90123456791_dp)*Hr4(0,-1,-1,0) + nf*z**5*(-22410.359396433472_dp)*Hr4(0,-1,-1,0) + &
                        nf**2*z**2*(-1494.9135802469136_dp)*Hr4(0,-1,-1,0) + nf**2*z**3*(-1494.9135802469136_dp)*Hr4(0,-1,-1,0) + &
                        nf**2*z*(1494.9135802469136_dp)*Hr4(0,-1,-1,0) + nf**2*z**4*(1494.9135802469136_dp)*Hr4(0,-1,-1,0) + &
                        nf*(34338.765432098764_dp)*Hr4(0,-1,-1,0) + nf*z**3*(46994.26063100137_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**2*(61731.29218106996_dp)*Hr4(0,-1,-1,0) + z*(93311.99999999991_dp)*Hr4(0,-1,-1,0) + &
                        z**5*(1.254528e6_dp)*Hr4(0,-1,-1,0) + z**4*(2.270592e6_dp)*Hr4(0,-1,-1,0) + &
                        z**4*(-1.9440000000000002e6_dp)*Hr4(0,-1,0,0) + z**5*(-1.083456e6_dp)*Hr4(0,-1,0,0) + &
                        z*(-67391.99999999994_dp)*Hr4(0,-1,0,0) + nf*z**2*(-50753.58024691356_dp)*Hr4(0,-1,0,0) + &
                        nf*z**3*(-39906.23868312756_dp)*Hr4(0,-1,0,0) + nf*(-29372.04938271605_dp)*Hr4(0,-1,0,0) + &
                        nf**2*z*(-1287.9012345679012_dp)*Hr4(0,-1,0,0) + nf**2*z**4*(-1287.9012345679012_dp)*Hr4(0,-1,0,0) + &
                        nf**2*z**2*(1287.9012345679012_dp)*Hr4(0,-1,0,0) + nf**2*z**3*(1287.9012345679012_dp)*Hr4(0,-1,0,0) + &
                        nf*z**5*(18738.041152263373_dp)*Hr4(0,-1,0,0) + nf*z*(38448.1975308642_dp)*Hr4(0,-1,0,0) + &
                        nf*z**4*(62845.62962962963_dp)*Hr4(0,-1,0,0) + z**3*(865728._dp)*Hr4(0,-1,0,0) + &
                        (1.083456e6_dp)*Hr4(0,-1,0,0) + z**2*(1.145664e6_dp)*Hr4(0,-1,0,0) + &
                        z**4*(-1.6174079999999998e6_dp)*Hr4(0,-1,0,1) + z**5*(-912384.0000000001_dp)*Hr4(0,-1,0,1) + &
                        z*(-41471.99999999991_dp)*Hr4(0,-1,0,1) + nf*z**2*(-39775.86831275721_dp)*Hr4(0,-1,0,1) + &
                        nf*z**3*(-32818.21673525377_dp)*Hr4(0,-1,0,1) + nf*(-24405.333333333332_dp)*Hr4(0,-1,0,1) + &
                        nf**2*z*(-1080.888888888889_dp)*Hr4(0,-1,0,1) + nf**2*z**4*(-1080.888888888889_dp)*Hr4(0,-1,0,1) + &
                        nf**2*z**2*(1080.888888888889_dp)*Hr4(0,-1,0,1) + nf**2*z**3*(1080.888888888889_dp)*Hr4(0,-1,0,1) + &
                        nf*z**5*(15065.722908093281_dp)*Hr4(0,-1,0,1) + nf*z*(31576.49382716049_dp)*Hr4(0,-1,0,1) + &
                        nf*z**4*(50357.201646090536_dp)*Hr4(0,-1,0,1) + z**3*(725760._dp)*Hr4(0,-1,0,1) + &
                        (912384.0000000001_dp)*Hr4(0,-1,0,1) + z**2*(933119.9999999998_dp)*Hr4(0,-1,0,1) + &
                        z**4*(-1.544832e6_dp)*Hr4(0,0,-1,0) + z**5*(-418175.99999999994_dp)*Hr4(0,0,-1,0) + &
                        nf*z**2*(-84546.37037037038_dp)*Hr4(0,0,-1,0) + z**3*(-58751.99999999992_dp)*Hr4(0,0,-1,0) + &
                        nf*z*(-30363.390946502055_dp)*Hr4(0,0,-1,0) + nf*(-13055.999999999998_dp)*Hr4(0,0,-1,0) + &
                        nf**2*z**4*(-2439.901234567901_dp)*Hr4(0,0,-1,0) + nf**2*z**3*(-944.9876543209875_dp)*Hr4(0,0,-1,0) + &
                        nf**2*z*(944.9876543209875_dp)*Hr4(0,0,-1,0) + nf**2*z**2*(2439.901234567901_dp)*Hr4(0,0,-1,0) + &
                        nf*z**5*(7091.796982167352_dp)*Hr4(0,0,-1,0) + nf*z**3*(30183.5939643347_dp)*Hr4(0,0,-1,0) + &
                        nf*z**4*(90690.37037037038_dp)*Hr4(0,0,-1,0) + z*(362880.00000000006_dp)*Hr4(0,0,-1,0) + &
                        (418175.99999999994_dp)*Hr4(0,0,-1,0) + z**2*(1.240704e6_dp)*Hr4(0,0,-1,0) + &
                        z**4*(-1.7141760000000002e6_dp)*Hr4(0,0,0,0) + z**3*(-352512.0000000001_dp)*Hr4(0,0,0,0) + &
                        nf*z**2*(-208249.1522633744_dp)*Hr4(0,0,0,0) + nf*z*(-35523.95061728395_dp)*Hr4(0,0,0,0) + &
                        z*(-27648._dp)*Hr4(0,0,0,0) + nf*z**5*(-7046.847736625516_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**4*(-6584.3621399176955_dp)*Hr4(0,0,0,0) + nf**2*z**3*(-1729.843621399177_dp)*Hr4(0,0,0,0) + &
                        nf**3*z*(-50.56790123456789_dp)*Hr4(0,0,0,0) + nf**3*z**2*(-50.56790123456789_dp)*Hr4(0,0,0,0) + &
                        nf**3*z**3*(50.56790123456789_dp)*Hr4(0,0,0,0) + nf**3*z**4*(50.56790123456789_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(1729.843621399177_dp)*Hr4(0,0,0,0) + nf**2*z**2*(6584.3621399176955_dp)*Hr4(0,0,0,0) + &
                        nf*z**3*(35658.79835390948_dp)*Hr4(0,0,0,0) + nf*z**4*(208249.15226337445_dp)*Hr4(0,0,0,0) + &
                        z**5*(494208._dp)*Hr4(0,0,0,0) + z**2*(1.7141760000000002e6_dp)*Hr4(0,0,0,0) + &
                        z**4*(-3.224448e6_dp)*Hr4(0,0,0,1) + z*(-1.1612160000000002e6_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(-257705.08641975312_dp)*Hr4(0,0,0,1) + nf*z**3*(-66720.92181069961_dp)*Hr4(0,0,0,1) + &
                        nf*z**5*(-21511.37448559671_dp)*Hr4(0,0,0,1) + nf*(-11520.000000000002_dp)*Hr4(0,0,0,1) + &
                        nf**2*z**4*(-4592.197530864198_dp)*Hr4(0,0,0,1) + nf**2*z*(-1961.0864197530864_dp)*Hr4(0,0,0,1) + &
                        nf**3*z**3*(-12.641975308641973_dp)*Hr4(0,0,0,1) + nf**3*z**4*(-12.641975308641973_dp)*Hr4(0,0,0,1) + &
                        nf**3*z*(12.641975308641973_dp)*Hr4(0,0,0,1) + nf**3*z**2*(12.641975308641973_dp)*Hr4(0,0,0,1) + &
                        nf**2*z**3*(1961.086419753086_dp)*Hr4(0,0,0,1) + nf**2*z**2*(4592.197530864198_dp)*Hr4(0,0,0,1) + &
                        z**3*(20736.000000000073_dp)*Hr4(0,0,0,1) + nf*z*(67496.29629629629_dp)*Hr4(0,0,0,1) + &
                        nf*z**4*(262313.0864197531_dp)*Hr4(0,0,0,1) + (342144._dp)*Hr4(0,0,0,1) + z**5*(1.482624e6_dp)*Hr4(0,0,0,1) &
                        + z**2*(2.9963519999999995e6_dp)*Hr4(0,0,0,1) + z**4*(-3.431808e6_dp)*Hr4(0,0,1,0) + &
                        z**3*(-998784.0000000007_dp)*Hr4(0,0,1,0) + z*(-445824.0000000001_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(-166227.4897119342_dp)*Hr4(0,0,1,0) + nf*(-33080.88888888889_dp)*Hr4(0,0,1,0) + &
                        nf*z**5*(-29763.42386831276_dp)*Hr4(0,0,1,0) + nf*z*(-9139.621399176958_dp)*Hr4(0,0,1,0) + &
                        nf**2*z**4*(-2262.38683127572_dp)*Hr4(0,0,1,0) + nf**2*z**3*(-1040.855967078189_dp)*Hr4(0,0,1,0) + &
                        nf**2*z*(1040.8559670781892_dp)*Hr4(0,0,1,0) + nf**2*z**2*(2262.3868312757204_dp)*Hr4(0,0,1,0) + &
                        nf*z**3*(11255.0452674897_dp)*Hr4(0,0,1,0) + nf*z**4*(178572.37860082305_dp)*Hr4(0,0,1,0) + &
                        (1.102464e6_dp)*Hr4(0,0,1,0) + z**5*(1.9008e6_dp)*Hr4(0,0,1,0) + z**2*(2.6714879999999995e6_dp)*Hr4(0,0,1,0) &
                        + z**4*(-3.773952000000001e6_dp)*Hr4(0,0,1,1) + z**3*(-2.1150720000000005e6_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(-174329.1522633745_dp)*Hr4(0,0,1,1) + z*(-165888.00000000023_dp)*Hr4(0,0,1,1) + &
                        nf*z*(-70986.79835390946_dp)*Hr4(0,0,1,1) + nf*z**5*(-45230.88065843622_dp)*Hr4(0,0,1,1) + &
                        nf*(-41756.444444444445_dp)*Hr4(0,0,1,1) + nf**2*z**3*(-3683.028806584362_dp)*Hr4(0,0,1,1) + &
                        nf**2*z**4*(-298.13991769547323_dp)*Hr4(0,0,1,1) + nf**2*z**2*(298.1399176954734_dp)*Hr4(0,0,1,1) + &
                        nf**2*z*(3683.028806584362_dp)*Hr4(0,0,1,1) + nf*z**3*(74745.67901234569_dp)*Hr4(0,0,1,1) + &
                        nf*z**4*(188437.59670781894_dp)*Hr4(0,0,1,1) + (1.368576e6_dp)*Hr4(0,0,1,1) + &
                        z**2*(2.861568e6_dp)*Hr4(0,0,1,1) + z**5*(2.965248e6_dp)*Hr4(0,0,1,1) + &
                        z**4*(-2.2826879999999995e6_dp)*Hr4(0,1,0,0) + z**3*(-858816.0000000001_dp)*Hr4(0,1,0,0) + &
                        z*(-319680.00000000023_dp)*Hr4(0,1,0,0) + nf*z**2*(-56925.23456790124_dp)*Hr4(0,1,0,0) + &
                        nf*(-36669.62962962963_dp)*Hr4(0,1,0,0) + nf*z*(-33124.082304526746_dp)*Hr4(0,1,0,0) + &
                        nf*z**5*(-27854.661179698214_dp)*Hr4(0,1,0,0) + nf**2*z**3*(-791.1769547325102_dp)*Hr4(0,1,0,0) + &
                        nf**2*z**4*(-791.1769547325102_dp)*Hr4(0,1,0,0) + nf**2*z*(791.1769547325102_dp)*Hr4(0,1,0,0) + &
                        nf**2*z**2*(791.1769547325102_dp)*Hr4(0,1,0,0) + nf*z**3*(36786.743484224964_dp)*Hr4(0,1,0,0) + &
                        nf*z**4*(69402.86419753087_dp)*Hr4(0,1,0,0) + z**2*(1.2562559999999995e6_dp)*Hr4(0,1,0,0) + &
                        (1.4255999999999998e6_dp)*Hr4(0,1,0,0) + z**5*(1.577664e6_dp)*Hr4(0,1,0,0) + &
                        z**4*(-3.2866559999999995e6_dp)*Hr4(0,1,0,1) + z*(-1.0679039999999998e6_dp)*Hr4(0,1,0,1) + &
                        z**3*(-984959.9999999993_dp)*Hr4(0,1,0,1) + nf*z**2*(-80291.55555555558_dp)*Hr4(0,1,0,1) + &
                        nf*(-59644.83950617284_dp)*Hr4(0,1,0,1) + nf*z**5*(-49535.47325102881_dp)*Hr4(0,1,0,1) + &
                        nf*z*(-35197.629629629635_dp)*Hr4(0,1,0,1) + nf**2*z**3*(-369.25102880658443_dp)*Hr4(0,1,0,1) + &
                        nf**2*z**4*(-369.25102880658443_dp)*Hr4(0,1,0,1) + nf**2*z*(369.25102880658443_dp)*Hr4(0,1,0,1) + &
                        nf**2*z**2*(369.25102880658443_dp)*Hr4(0,1,0,1) + nf*z**3*(43261.102880658415_dp)*Hr4(0,1,0,1) + &
                        nf*z**4*(98464.39506172837_dp)*Hr4(0,1,0,1) + z**2*(1.689984e6_dp)*Hr4(0,1,0,1) + &
                        (2.28096e6_dp)*Hr4(0,1,0,1) + z**5*(2.737152e6_dp)*Hr4(0,1,0,1) + &
                        z**4*(-3.2866559999999995e6_dp)*Hr4(0,1,1,0) + z*(-1.0679039999999998e6_dp)*Hr4(0,1,1,0) + &
                        z**3*(-984959.9999999993_dp)*Hr4(0,1,1,0) + nf*z**2*(-80291.55555555558_dp)*Hr4(0,1,1,0) + &
                        nf*(-59644.83950617284_dp)*Hr4(0,1,1,0) + nf*z**5*(-49535.47325102881_dp)*Hr4(0,1,1,0) + &
                        nf*z*(-35197.629629629635_dp)*Hr4(0,1,1,0) + nf**2*z**3*(-369.25102880658443_dp)*Hr4(0,1,1,0) + &
                        nf**2*z**4*(-369.25102880658443_dp)*Hr4(0,1,1,0) + nf**2*z*(369.25102880658443_dp)*Hr4(0,1,1,0) + &
                        nf**2*z**2*(369.25102880658443_dp)*Hr4(0,1,1,0) + nf*z**3*(43261.102880658415_dp)*Hr4(0,1,1,0) + &
                        nf*z**4*(98464.39506172837_dp)*Hr4(0,1,1,0) + z**2*(1.689984e6_dp)*Hr4(0,1,1,0) + &
                        (2.28096e6_dp)*Hr4(0,1,1,0) + z**5*(2.737152e6_dp)*Hr4(0,1,1,0) + z*(-2.3224320000000014e6_dp)*Hr4(0,1,1,1) &
                        + z**4*(-2.0321279999999998e6_dp)*Hr4(0,1,1,1) + nf*(-54025.48148148148_dp)*Hr4(0,1,1,1) + &
                        nf*z**5*(-54025.48148148148_dp)*Hr4(0,1,1,1) + nf*z**2*(-37165.82716049383_dp)*Hr4(0,1,1,1) + &
                        nf**2*z*(-1886.8148148148148_dp)*Hr4(0,1,1,1) + nf**2*z**2*(-1886.8148148148148_dp)*Hr4(0,1,1,1) + &
                        nf*z**3*(-1840.9876543209996_dp)*Hr4(0,1,1,1) + nf**2*z**3*(1886.8148148148148_dp)*Hr4(0,1,1,1) + &
                        nf**2*z**4*(1886.8148148148148_dp)*Hr4(0,1,1,1) + nf*z*(14394.469135802481_dp)*Hr4(0,1,1,1) + &
                        z**3*(41472.00000000132_dp)*Hr4(0,1,1,1) + nf*z**4*(49719.30864197531_dp)*Hr4(0,1,1,1) + &
                        z**2*(663551.9999999999_dp)*Hr4(0,1,1,1) + (2.052864e6_dp)*Hr4(0,1,1,1) + z**5*(2.965248e6_dp)*Hr4(0,1,1,1) &
                        + z**4*(-787968.0000000002_dp)*Hr4(1,0,0,0) + z*(-787967.9999999999_dp)*Hr4(1,0,0,0) + &
                        nf*(-18938.38134430727_dp)*Hr4(1,0,0,0) + nf*z**5*(-18938.38134430727_dp)*Hr4(1,0,0,0) + &
                        nf*z**3*(-4873.832647462276_dp)*Hr4(1,0,0,0) + nf*z**2*(-4873.832647462272_dp)*Hr4(1,0,0,0) + &
                        nf*z**4*(9988.213991769546_dp)*Hr4(1,0,0,0) + nf*z*(9988.21399176955_dp)*Hr4(1,0,0,0) + &
                        z**3*(103679.99999999978_dp)*Hr4(1,0,0,0) + z**2*(103680.00000000015_dp)*Hr4(1,0,0,0) + &
                        (912384.0000000001_dp)*Hr4(1,0,0,0) + z**5*(912384.0000000001_dp)*Hr4(1,0,0,0) + &
                        z*(-1.8766080000000005e6_dp)*Hr4(1,0,0,1) + z**4*(-1.8766079999999998e6_dp)*Hr4(1,0,0,1) + &
                        nf*(-46515.79698216737_dp)*Hr4(1,0,0,1) + nf*z**5*(-46515.79698216737_dp)*Hr4(1,0,0,1) + &
                        nf*z**3*(-13637.355281207132_dp)*Hr4(1,0,0,1) + nf*z**2*(-13637.355281207125_dp)*Hr4(1,0,0,1) + &
                        nf*z*(25593.152263374475_dp)*Hr4(1,0,0,1) + nf*z**4*(25593.152263374475_dp)*Hr4(1,0,0,1) + &
                        z**2*(279935.99999999953_dp)*Hr4(1,0,0,1) + z**3*(279935.9999999998_dp)*Hr4(1,0,0,1) + &
                        (2.166912e6_dp)*Hr4(1,0,0,1) + z**5*(2.166912e6_dp)*Hr4(1,0,0,1) + z**4*(-2.17728e6_dp)*Hr4(1,0,1,0) + &
                        z*(-2.1772799999999995e6_dp)*Hr4(1,0,1,0) + nf*(-55154.831275720164_dp)*Hr4(1,0,1,0) + &
                        nf*z**5*(-55154.831275720164_dp)*Hr4(1,0,1,0) + nf*z**2*(-17527.04526748971_dp)*Hr4(1,0,1,0) + &
                        nf*z**3*(-17527.04526748971_dp)*Hr4(1,0,1,0) + nf*z*(31209.87654320987_dp)*Hr4(1,0,1,0) + &
                        nf*z**4*(31209.87654320987_dp)*Hr4(1,0,1,0) + z**3*(352511.99999999895_dp)*Hr4(1,0,1,0) + &
                        z**2*(352511.9999999994_dp)*Hr4(1,0,1,0) + (2.509056e6_dp)*Hr4(1,0,1,0) + z**5*(2.509056e6_dp)*Hr4(1,0,1,0) &
                        + z**4*(-2.17728e6_dp)*Hr4(1,0,1,1) + z*(-2.1772799999999995e6_dp)*Hr4(1,0,1,1) + &
                        nf*(-55154.831275720164_dp)*Hr4(1,0,1,1) + nf*z**5*(-55154.831275720164_dp)*Hr4(1,0,1,1) + &
                        nf*z**2*(-17527.04526748971_dp)*Hr4(1,0,1,1) + nf*z**3*(-17527.04526748971_dp)*Hr4(1,0,1,1) + &
                        nf*z*(31209.87654320987_dp)*Hr4(1,0,1,1) + nf*z**4*(31209.87654320987_dp)*Hr4(1,0,1,1) + &
                        z**3*(352511.99999999895_dp)*Hr4(1,0,1,1) + z**2*(352511.9999999994_dp)*Hr4(1,0,1,1) + &
                        (2.509056e6_dp)*Hr4(1,0,1,1) + z**5*(2.509056e6_dp)*Hr4(1,0,1,1) + z*(-1.8766080000000005e6_dp)*Hr4(1,1,0,0) &
                        + z**4*(-1.8766079999999998e6_dp)*Hr4(1,1,0,0) + nf*(-46515.79698216737_dp)*Hr4(1,1,0,0) + &
                        nf*z**5*(-46515.79698216737_dp)*Hr4(1,1,0,0) + nf*z**3*(-13637.355281207132_dp)*Hr4(1,1,0,0) + &
                        nf*z**2*(-13637.355281207125_dp)*Hr4(1,1,0,0) + nf*z*(25593.152263374475_dp)*Hr4(1,1,0,0) + &
                        nf*z**4*(25593.152263374475_dp)*Hr4(1,1,0,0) + z**2*(279935.99999999953_dp)*Hr4(1,1,0,0) + &
                        z**3*(279935.9999999998_dp)*Hr4(1,1,0,0) + (2.166912e6_dp)*Hr4(1,1,0,0) + z**5*(2.166912e6_dp)*Hr4(1,1,0,0) &
                        + z**4*(-2.17728e6_dp)*Hr4(1,1,0,1) + z*(-2.1772799999999995e6_dp)*Hr4(1,1,0,1) + &
                        nf*(-55154.831275720164_dp)*Hr4(1,1,0,1) + nf*z**5*(-55154.831275720164_dp)*Hr4(1,1,0,1) + &
                        nf*z**2*(-17527.04526748971_dp)*Hr4(1,1,0,1) + nf*z**3*(-17527.04526748971_dp)*Hr4(1,1,0,1) + &
                        nf*z*(31209.87654320987_dp)*Hr4(1,1,0,1) + nf*z**4*(31209.87654320987_dp)*Hr4(1,1,0,1) + &
                        z**3*(352511.99999999895_dp)*Hr4(1,1,0,1) + z**2*(352511.9999999994_dp)*Hr4(1,1,0,1) + &
                        (2.509056e6_dp)*Hr4(1,1,0,1) + z**5*(2.509056e6_dp)*Hr4(1,1,0,1) + z**4*(-2.17728e6_dp)*Hr4(1,1,1,0) + &
                        z*(-2.1772799999999995e6_dp)*Hr4(1,1,1,0) + nf*(-55154.831275720164_dp)*Hr4(1,1,1,0) + &
                        nf*z**5*(-55154.831275720164_dp)*Hr4(1,1,1,0) + nf*z**2*(-17527.04526748971_dp)*Hr4(1,1,1,0) + &
                        nf*z**3*(-17527.04526748971_dp)*Hr4(1,1,1,0) + nf*z*(31209.87654320987_dp)*Hr4(1,1,1,0) + &
                        nf*z**4*(31209.87654320987_dp)*Hr4(1,1,1,0) + z**3*(352511.99999999895_dp)*Hr4(1,1,1,0) + &
                        z**2*(352511.9999999994_dp)*Hr4(1,1,1,0) + (2.509056e6_dp)*Hr4(1,1,1,0) + z**5*(2.509056e6_dp)*Hr4(1,1,1,0) &
                        + (-995328.0000000002_dp)*Hr5(-1,-1,-1,-1,0) + z*(-995328.0000000002_dp)*Hr5(-1,-1,-1,-1,0) + &
                        z**2*(-995328.0000000002_dp)*Hr5(-1,-1,-1,-1,0) + z**3*(995328.0000000002_dp)*Hr5(-1,-1,-1,-1,0) + &
                        z**4*(995328.0000000002_dp)*Hr5(-1,-1,-1,-1,0) + z**5*(995328.0000000002_dp)*Hr5(-1,-1,-1,-1,0) + &
                        z**3*(-870912.0000000001_dp)*Hr5(-1,-1,-1,0,0) + z**4*(-870912.0000000001_dp)*Hr5(-1,-1,-1,0,0) + &
                        z**5*(-870912.0000000001_dp)*Hr5(-1,-1,-1,0,0) + z**2*(870912._dp)*Hr5(-1,-1,-1,0,0) + &
                        (870912.0000000001_dp)*Hr5(-1,-1,-1,0,0) + z*(870912.0000000002_dp)*Hr5(-1,-1,-1,0,0) + &
                        z**3*(-746496._dp)*Hr5(-1,-1,-1,0,1) + z**4*(-746496._dp)*Hr5(-1,-1,-1,0,1) + &
                        z**5*(-746496._dp)*Hr5(-1,-1,-1,0,1) + (746496._dp)*Hr5(-1,-1,-1,0,1) + z*(746496._dp)*Hr5(-1,-1,-1,0,1) + &
                        z**2*(746496._dp)*Hr5(-1,-1,-1,0,1) + z**3*(-497664.0000000001_dp)*Hr5(-1,-1,0,-1,0) + &
                        z**4*(-497664.0000000001_dp)*Hr5(-1,-1,0,-1,0) + z**5*(-497664.0000000001_dp)*Hr5(-1,-1,0,-1,0) + &
                        (497664.0000000001_dp)*Hr5(-1,-1,0,-1,0) + z*(497664.0000000001_dp)*Hr5(-1,-1,0,-1,0) + &
                        z**2*(497664.0000000001_dp)*Hr5(-1,-1,0,-1,0) + (-373248._dp)*Hr5(-1,-1,0,0,0) + &
                        z*(-373248._dp)*Hr5(-1,-1,0,0,0) + z**2*(-373248._dp)*Hr5(-1,-1,0,0,0) + z**3*(373248._dp)*Hr5(-1,-1,0,0,0) &
                        + z**4*(373248._dp)*Hr5(-1,-1,0,0,0) + z**5*(373248._dp)*Hr5(-1,-1,0,0,0) + (-622080._dp)*Hr5(-1,-1,0,0,1) + &
                        z*(-622080._dp)*Hr5(-1,-1,0,0,1) + z**2*(-622080._dp)*Hr5(-1,-1,0,0,1) + z**3*(622080._dp)*Hr5(-1,-1,0,0,1) &
                        + z**4*(622080._dp)*Hr5(-1,-1,0,0,1) + z**5*(622080._dp)*Hr5(-1,-1,0,0,1) + &
                        (-248832.00000000006_dp)*Hr5(-1,-1,0,1,0) + z*(-248832.00000000006_dp)*Hr5(-1,-1,0,1,0) + &
                        z**2*(-248832.00000000006_dp)*Hr5(-1,-1,0,1,0) + z**3*(248832.00000000006_dp)*Hr5(-1,-1,0,1,0) + &
                        z**4*(248832.00000000006_dp)*Hr5(-1,-1,0,1,0) + z**5*(248832.00000000006_dp)*Hr5(-1,-1,0,1,0) + &
                        (-497664.0000000001_dp)*Hr5(-1,-1,0,1,1) + z*(-497664.0000000001_dp)*Hr5(-1,-1,0,1,1) + &
                        z**2*(-497664.0000000001_dp)*Hr5(-1,-1,0,1,1) + z**3*(497664.0000000001_dp)*Hr5(-1,-1,0,1,1) + &
                        z**4*(497664.0000000001_dp)*Hr5(-1,-1,0,1,1) + z**5*(497664.0000000001_dp)*Hr5(-1,-1,0,1,1) + &
                        z**3*(-497664.0000000001_dp)*Hr5(-1,0,-1,-1,0) + z**4*(-497664.0000000001_dp)*Hr5(-1,0,-1,-1,0) + &
                        z**5*(-497664.0000000001_dp)*Hr5(-1,0,-1,-1,0) + (497664.0000000001_dp)*Hr5(-1,0,-1,-1,0) + &
                        z*(497664.0000000001_dp)*Hr5(-1,0,-1,-1,0) + z**2*(497664.0000000001_dp)*Hr5(-1,0,-1,-1,0) + &
                        z*(-435456.0000000001_dp)*Hr5(-1,0,-1,0,0) + (-435456.00000000006_dp)*Hr5(-1,0,-1,0,0) + &
                        z**2*(-435456._dp)*Hr5(-1,0,-1,0,0) + z**3*(435456.00000000006_dp)*Hr5(-1,0,-1,0,0) + &
                        z**4*(435456.00000000006_dp)*Hr5(-1,0,-1,0,0) + z**5*(435456.00000000006_dp)*Hr5(-1,0,-1,0,0) + &
                        (-373248._dp)*Hr5(-1,0,-1,0,1) + z*(-373248._dp)*Hr5(-1,0,-1,0,1) + z**2*(-373248._dp)*Hr5(-1,0,-1,0,1) + &
                        z**3*(373248._dp)*Hr5(-1,0,-1,0,1) + z**4*(373248._dp)*Hr5(-1,0,-1,0,1) + z**5*(373248._dp)*Hr5(-1,0,-1,0,1) &
                        + (-165888._dp)*Hr5(-1,0,0,-1,0) + z*(-165888._dp)*Hr5(-1,0,0,-1,0) + z**2*(-165888._dp)*Hr5(-1,0,0,-1,0) + &
                        z**3*(165888._dp)*Hr5(-1,0,0,-1,0) + z**4*(165888._dp)*Hr5(-1,0,0,-1,0) + z**5*(165888._dp)*Hr5(-1,0,0,-1,0) &
                        + z**3*(-103680._dp)*Hr5(-1,0,0,0,0) + z**4*(-103680._dp)*Hr5(-1,0,0,0,0) + &
                        z**5*(-103680._dp)*Hr5(-1,0,0,0,0) + (103680._dp)*Hr5(-1,0,0,0,0) + z*(103680._dp)*Hr5(-1,0,0,0,0) + &
                        z**2*(103680._dp)*Hr5(-1,0,0,0,0) + z**3*(-248832.00000000006_dp)*Hr5(-1,0,0,0,1) + &
                        z**4*(-248832.00000000006_dp)*Hr5(-1,0,0,0,1) + z**5*(-248832.00000000006_dp)*Hr5(-1,0,0,0,1) + &
                        (248832.00000000006_dp)*Hr5(-1,0,0,0,1) + z*(248832.00000000006_dp)*Hr5(-1,0,0,0,1) + &
                        z**2*(248832.00000000006_dp)*Hr5(-1,0,0,0,1) + z**3*(-186624._dp)*Hr5(-1,0,0,1,0) + &
                        z**4*(-186624._dp)*Hr5(-1,0,0,1,0) + z**5*(-186624._dp)*Hr5(-1,0,0,1,0) + (186624._dp)*Hr5(-1,0,0,1,0) + &
                        z*(186624._dp)*Hr5(-1,0,0,1,0) + z**2*(186624._dp)*Hr5(-1,0,0,1,0) + z**3*(-373248._dp)*Hr5(-1,0,0,1,1) + &
                        z**4*(-373248._dp)*Hr5(-1,0,0,1,1) + z**5*(-373248._dp)*Hr5(-1,0,0,1,1) + (373248._dp)*Hr5(-1,0,0,1,1) + &
                        z*(373248._dp)*Hr5(-1,0,0,1,1) + z**2*(373248._dp)*Hr5(-1,0,0,1,1) + z**3*(-41472._dp)*Hr5(-1,0,1,0,0) + &
                        z**4*(-41472._dp)*Hr5(-1,0,1,0,0) + z**5*(-41472._dp)*Hr5(-1,0,1,0,0) + (41472._dp)*Hr5(-1,0,1,0,0) + &
                        z*(41472._dp)*Hr5(-1,0,1,0,0) + z**2*(41472._dp)*Hr5(-1,0,1,0,0) + &
                        z**3*(-124416.00000000003_dp)*Hr5(-1,0,1,0,1) + z**4*(-124416.00000000003_dp)*Hr5(-1,0,1,0,1) + &
                        z**5*(-124416.00000000003_dp)*Hr5(-1,0,1,0,1) + (124416.00000000003_dp)*Hr5(-1,0,1,0,1) + &
                        z*(124416.00000000003_dp)*Hr5(-1,0,1,0,1) + z**2*(124416.00000000003_dp)*Hr5(-1,0,1,0,1) + &
                        z**3*(-124416.00000000003_dp)*Hr5(-1,0,1,1,0) + z**4*(-124416.00000000003_dp)*Hr5(-1,0,1,1,0) + &
                        z**5*(-124416.00000000003_dp)*Hr5(-1,0,1,1,0) + (124416.00000000003_dp)*Hr5(-1,0,1,1,0) + &
                        z*(124416.00000000003_dp)*Hr5(-1,0,1,1,0) + z**2*(124416.00000000003_dp)*Hr5(-1,0,1,1,0) + &
                        z**3*(-248832.00000000006_dp)*Hr5(-1,0,1,1,1) + z**4*(-248832.00000000006_dp)*Hr5(-1,0,1,1,1) + &
                        z**5*(-248832.00000000006_dp)*Hr5(-1,0,1,1,1) + (248832.00000000006_dp)*Hr5(-1,0,1,1,1) + &
                        z*(248832.00000000006_dp)*Hr5(-1,0,1,1,1) + z**2*(248832.00000000006_dp)*Hr5(-1,0,1,1,1) + &
                        z**4*(-1.492992e6_dp)*Hr5(0,-1,-1,-1,0) + z*(-497664.0000000001_dp)*Hr5(0,-1,-1,-1,0) + &
                        z**5*(-497664.0000000001_dp)*Hr5(0,-1,-1,-1,0) + nf*z**2*(-22915.687242798354_dp)*Hr5(0,-1,-1,-1,0) + &
                        nf*z**3*(-22915.687242798354_dp)*Hr5(0,-1,-1,-1,0) + nf*z*(22915.687242798354_dp)*Hr5(0,-1,-1,-1,0) + &
                        nf*z**4*(22915.687242798354_dp)*Hr5(0,-1,-1,-1,0) + (497664.0000000001_dp)*Hr5(0,-1,-1,-1,0) + &
                        z**3*(497664.0000000001_dp)*Hr5(0,-1,-1,-1,0) + z**2*(1.492992e6_dp)*Hr5(0,-1,-1,-1,0) + &
                        z**2*(-1.306368e6_dp)*Hr5(0,-1,-1,0,0) + z**3*(-435456.0000000001_dp)*Hr5(0,-1,-1,0,0) + &
                        (-435456.00000000006_dp)*Hr5(0,-1,-1,0,0) + nf*z*(-20325.135802469136_dp)*Hr5(0,-1,-1,0,0) + &
                        nf*z**4*(-20325.135802469136_dp)*Hr5(0,-1,-1,0,0) + nf*z**2*(20325.135802469136_dp)*Hr5(0,-1,-1,0,0) + &
                        nf*z**3*(20325.135802469136_dp)*Hr5(0,-1,-1,0,0) + z*(435456._dp)*Hr5(0,-1,-1,0,0) + &
                        z**5*(435456.00000000006_dp)*Hr5(0,-1,-1,0,0) + z**4*(1.306368e6_dp)*Hr5(0,-1,-1,0,0) + &
                        z**2*(-1.119744e6_dp)*Hr5(0,-1,-1,0,1) + (-373248._dp)*Hr5(0,-1,-1,0,1) + &
                        z**3*(-373248._dp)*Hr5(0,-1,-1,0,1) + nf*z*(-17734.584362139918_dp)*Hr5(0,-1,-1,0,1) + &
                        nf*z**4*(-17734.584362139918_dp)*Hr5(0,-1,-1,0,1) + nf*z**2*(17734.584362139918_dp)*Hr5(0,-1,-1,0,1) + &
                        nf*z**3*(17734.584362139918_dp)*Hr5(0,-1,-1,0,1) + z*(373248._dp)*Hr5(0,-1,-1,0,1) + &
                        z**5*(373248._dp)*Hr5(0,-1,-1,0,1) + z**4*(1.119744e6_dp)*Hr5(0,-1,-1,0,1) + &
                        z**2*(-746496._dp)*Hr5(0,-1,0,-1,0) + (-248832.00000000006_dp)*Hr5(0,-1,0,-1,0) + &
                        z**3*(-248832.00000000006_dp)*Hr5(0,-1,0,-1,0) + nf*z*(-11457.843621399177_dp)*Hr5(0,-1,0,-1,0) + &
                        nf*z**4*(-11457.843621399177_dp)*Hr5(0,-1,0,-1,0) + nf*z**2*(11457.843621399177_dp)*Hr5(0,-1,0,-1,0) + &
                        nf*z**3*(11457.843621399177_dp)*Hr5(0,-1,0,-1,0) + z*(248832.00000000006_dp)*Hr5(0,-1,0,-1,0) + &
                        z**5*(248832.00000000006_dp)*Hr5(0,-1,0,-1,0) + z**4*(746496._dp)*Hr5(0,-1,0,-1,0) + &
                        z**4*(-559872._dp)*Hr5(0,-1,0,0,0) + z*(-186624._dp)*Hr5(0,-1,0,0,0) + z**5*(-186624._dp)*Hr5(0,-1,0,0,0) + &
                        nf*z**2*(-8867.292181069959_dp)*Hr5(0,-1,0,0,0) + nf*z**3*(-8867.292181069959_dp)*Hr5(0,-1,0,0,0) + &
                        nf*z*(8867.292181069959_dp)*Hr5(0,-1,0,0,0) + nf*z**4*(8867.292181069959_dp)*Hr5(0,-1,0,0,0) + &
                        (186624._dp)*Hr5(0,-1,0,0,0) + z**3*(186624._dp)*Hr5(0,-1,0,0,0) + z**2*(559872._dp)*Hr5(0,-1,0,0,0) + &
                        z**4*(-933120._dp)*Hr5(0,-1,0,0,1) + z*(-311040._dp)*Hr5(0,-1,0,0,1) + z**5*(-311040._dp)*Hr5(0,-1,0,0,1) + &
                        nf*z**2*(-15144.032921810702_dp)*Hr5(0,-1,0,0,1) + nf*z**3*(-15144.032921810702_dp)*Hr5(0,-1,0,0,1) + &
                        nf*z*(15144.032921810702_dp)*Hr5(0,-1,0,0,1) + nf*z**4*(15144.032921810702_dp)*Hr5(0,-1,0,0,1) + &
                        (311040._dp)*Hr5(0,-1,0,0,1) + z**3*(311040._dp)*Hr5(0,-1,0,0,1) + z**2*(933120._dp)*Hr5(0,-1,0,0,1) + &
                        z**4*(-373248._dp)*Hr5(0,-1,0,1,0) + z*(-124416.00000000003_dp)*Hr5(0,-1,0,1,0) + &
                        z**5*(-124416.00000000003_dp)*Hr5(0,-1,0,1,0) + nf*z**2*(-6276.740740740741_dp)*Hr5(0,-1,0,1,0) + &
                        nf*z**3*(-6276.740740740741_dp)*Hr5(0,-1,0,1,0) + nf*z*(6276.740740740741_dp)*Hr5(0,-1,0,1,0) + &
                        nf*z**4*(6276.740740740741_dp)*Hr5(0,-1,0,1,0) + (124416.00000000003_dp)*Hr5(0,-1,0,1,0) + &
                        z**3*(124416.00000000003_dp)*Hr5(0,-1,0,1,0) + z**2*(373248._dp)*Hr5(0,-1,0,1,0) + &
                        z**4*(-746496._dp)*Hr5(0,-1,0,1,1) + z*(-248832.00000000006_dp)*Hr5(0,-1,0,1,1) + &
                        z**5*(-248832.00000000006_dp)*Hr5(0,-1,0,1,1) + nf*z**2*(-12553.481481481482_dp)*Hr5(0,-1,0,1,1) + &
                        nf*z**3*(-12553.481481481482_dp)*Hr5(0,-1,0,1,1) + nf*z*(12553.481481481482_dp)*Hr5(0,-1,0,1,1) + &
                        nf*z**4*(12553.481481481482_dp)*Hr5(0,-1,0,1,1) + (248832.00000000006_dp)*Hr5(0,-1,0,1,1) + &
                        z**3*(248832.00000000006_dp)*Hr5(0,-1,0,1,1) + z**2*(746496._dp)*Hr5(0,-1,0,1,1) + &
                        z**2*(-1.161216e6_dp)*Hr5(0,0,-1,-1,0) + (-165888._dp)*Hr5(0,0,-1,-1,0) + z*(-165888._dp)*Hr5(0,0,-1,-1,0) + &
                        nf*z**4*(-42064.0658436214_dp)*Hr5(0,0,-1,-1,0) + nf*z**3*(-19148.37860082305_dp)*Hr5(0,0,-1,-1,0) + &
                        nf*z*(19148.378600823045_dp)*Hr5(0,0,-1,-1,0) + nf*z**2*(42064.0658436214_dp)*Hr5(0,0,-1,-1,0) + &
                        z**3*(165888._dp)*Hr5(0,0,-1,-1,0) + z**5*(165888._dp)*Hr5(0,0,-1,-1,0) + &
                        z**4*(1.161216e6_dp)*Hr5(0,0,-1,-1,0) + z**4*(-1.0160639999999999e6_dp)*Hr5(0,0,-1,0,0) + &
                        z**3*(-145152.0000000001_dp)*Hr5(0,0,-1,0,0) + z**5*(-145152._dp)*Hr5(0,0,-1,0,0) + &
                        nf*z**2*(-37369.67901234568_dp)*Hr5(0,0,-1,0,0) + nf*z*(-17044.543209876545_dp)*Hr5(0,0,-1,0,0) + &
                        nf*z**3*(17044.543209876545_dp)*Hr5(0,0,-1,0,0) + nf*z**4*(37369.67901234568_dp)*Hr5(0,0,-1,0,0) + &
                        (145152._dp)*Hr5(0,0,-1,0,0) + z*(145152._dp)*Hr5(0,0,-1,0,0) + &
                        z**2*(1.0160639999999999e6_dp)*Hr5(0,0,-1,0,0) + z**4*(-870912.0000000001_dp)*Hr5(0,0,-1,0,1) + &
                        z**3*(-124416.00000000009_dp)*Hr5(0,0,-1,0,1) + z**5*(-124416.00000000003_dp)*Hr5(0,0,-1,0,1) + &
                        nf*z**2*(-32675.292181069955_dp)*Hr5(0,0,-1,0,1) + nf*z*(-14940.707818930041_dp)*Hr5(0,0,-1,0,1) + &
                        nf*z**3*(14940.707818930041_dp)*Hr5(0,0,-1,0,1) + nf*z**4*(32675.292181069955_dp)*Hr5(0,0,-1,0,1) + &
                        (124416.00000000003_dp)*Hr5(0,0,-1,0,1) + z*(124416.00000000003_dp)*Hr5(0,0,-1,0,1) + &
                        z**2*(870912._dp)*Hr5(0,0,-1,0,1) + z**4*(-622080._dp)*Hr5(0,0,0,-1,0) + z*(-41472._dp)*Hr5(0,0,0,-1,0) + &
                        z**5*(-41472._dp)*Hr5(0,0,0,-1,0) + nf*z**2*(-38903.5720164609_dp)*Hr5(0,0,0,-1,0) + &
                        nf*z**3*(-8297.34979423868_dp)*Hr5(0,0,0,-1,0) + nf**2*z*(-227.55555555555554_dp)*Hr5(0,0,0,-1,0) + &
                        nf**2*z**4*(-227.55555555555554_dp)*Hr5(0,0,0,-1,0) + nf**2*z**2*(227.55555555555554_dp)*Hr5(0,0,0,-1,0) + &
                        nf**2*z**3*(227.55555555555554_dp)*Hr5(0,0,0,-1,0) + nf*z*(8297.349794238684_dp)*Hr5(0,0,0,-1,0) + &
                        nf*z**4*(38903.57201646091_dp)*Hr5(0,0,0,-1,0) + (41472._dp)*Hr5(0,0,0,-1,0) + &
                        z**3*(41472._dp)*Hr5(0,0,0,-1,0) + z**2*(622080._dp)*Hr5(0,0,0,-1,0) + z**4*(-663552._dp)*Hr5(0,0,0,0,0) + &
                        nf*z**2*(-68232.9547325103_dp)*Hr5(0,0,0,0,0) + z*(-41472._dp)*Hr5(0,0,0,0,0) + &
                        nf*z**3*(-10509.695473251037_dp)*Hr5(0,0,0,0,0) + nf**2*z**4*(-1238.9135802469136_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z*(-581.5308641975308_dp)*Hr5(0,0,0,0,0) + nf**2*z**3*(581.530864197531_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z**2*(1238.9135802469136_dp)*Hr5(0,0,0,0,0) + nf*z*(10509.69547325103_dp)*Hr5(0,0,0,0,0) + &
                        z**3*(41472._dp)*Hr5(0,0,0,0,0) + z**5*(41472._dp)*Hr5(0,0,0,0,0) + &
                        nf*z**4*(68232.95473251029_dp)*Hr5(0,0,0,0,0) + z**2*(663552._dp)*Hr5(0,0,0,0,0) + &
                        z**4*(-1.6174079999999998e6_dp)*Hr5(0,0,0,0,1) + z**3*(-165888._dp)*Hr5(0,0,0,0,1) + &
                        nf*z**2*(-108485.00411522631_dp)*Hr5(0,0,0,0,1) + nf*z*(-29117.62962962963_dp)*Hr5(0,0,0,0,1) + &
                        nf**2*z**3*(-720.5925925925927_dp)*Hr5(0,0,0,0,1) + nf**2*z**4*(-720.5925925925927_dp)*Hr5(0,0,0,0,1) + &
                        nf**2*z*(720.5925925925927_dp)*Hr5(0,0,0,0,1) + nf**2*z**2*(720.5925925925927_dp)*Hr5(0,0,0,0,1) + &
                        nf*z**3*(29117.629629629646_dp)*Hr5(0,0,0,0,1) + (41472._dp)*Hr5(0,0,0,0,1) + &
                        nf*z**4*(108485.00411522634_dp)*Hr5(0,0,0,0,1) + z*(165888._dp)*Hr5(0,0,0,0,1) + &
                        z**5*(165888._dp)*Hr5(0,0,0,0,1) + z**2*(1.6174079999999998e6_dp)*Hr5(0,0,0,0,1) + &
                        z**4*(-1.9906560000000005e6_dp)*Hr5(0,0,0,1,0) + z*(-207360.0000000001_dp)*Hr5(0,0,0,1,0) + &
                        nf*z**2*(-86704.98765432097_dp)*Hr5(0,0,0,1,0) + nf*z**3*(-19060.938271604937_dp)*Hr5(0,0,0,1,0) + &
                        nf**2*z**3*(-164.3456790123457_dp)*Hr5(0,0,0,1,0) + nf**2*z**4*(-164.3456790123457_dp)*Hr5(0,0,0,1,0) + &
                        nf**2*z*(164.3456790123457_dp)*Hr5(0,0,0,1,0) + nf**2*z**2*(164.3456790123457_dp)*Hr5(0,0,0,1,0) + &
                        nf*z*(19060.938271604937_dp)*Hr5(0,0,0,1,0) + nf*z**4*(86704.98765432097_dp)*Hr5(0,0,0,1,0) + &
                        (165888._dp)*Hr5(0,0,0,1,0) + z**3*(207359.99999999983_dp)*Hr5(0,0,0,1,0) + z**5*(290304._dp)*Hr5(0,0,0,1,0) &
                        + z**2*(1.9906559999999993e6_dp)*Hr5(0,0,0,1,0) + z**4*(-2.737152e6_dp)*Hr5(0,0,0,1,1) + &
                        z*(-497664._dp)*Hr5(0,0,0,1,1) + nf*z**2*(-105800.6913580247_dp)*Hr5(0,0,0,1,1) + &
                        nf*z**3*(-52198.716049382725_dp)*Hr5(0,0,0,1,1) + nf*z*(52198.71604938272_dp)*Hr5(0,0,0,1,1) + &
                        nf*z**4*(105800.69135802468_dp)*Hr5(0,0,0,1,1) + (248832.00000000006_dp)*Hr5(0,0,0,1,1) + &
                        z**3*(497664.0000000001_dp)*Hr5(0,0,0,1,1) + z**5*(497664.0000000001_dp)*Hr5(0,0,0,1,1) + &
                        z**2*(2.737152e6_dp)*Hr5(0,0,0,1,1) + z**4*(-1.472256e6_dp)*Hr5(0,0,1,0,0) + &
                        nf*z**2*(-39939.16049382716_dp)*Hr5(0,0,1,0,0) + z**3*(-20736.000000000073_dp)*Hr5(0,0,1,0,0) + &
                        nf*z**3*(-6795.061728395059_dp)*Hr5(0,0,1,0,0) + nf*z*(6795.061728395061_dp)*Hr5(0,0,1,0,0) + &
                        z*(20736._dp)*Hr5(0,0,1,0,0) + nf*z**4*(39939.16049382716_dp)*Hr5(0,0,1,0,0) + (269568._dp)*Hr5(0,0,1,0,0) + &
                        z**5*(311040._dp)*Hr5(0,0,1,0,0) + z**2*(1.472256e6_dp)*Hr5(0,0,1,0,0) + &
                        z**4*(-2.5297919999999995e6_dp)*Hr5(0,0,1,0,1) + z**3*(-331776._dp)*Hr5(0,0,1,0,1) + &
                        nf*z**2*(-59843.95061728395_dp)*Hr5(0,0,1,0,1) + nf*z*(-6311.506172839506_dp)*Hr5(0,0,1,0,1) + &
                        nf*z**3*(6311.506172839506_dp)*Hr5(0,0,1,0,1) + nf*z**4*(59843.95061728395_dp)*Hr5(0,0,1,0,1) + &
                        z*(331776.00000000006_dp)*Hr5(0,0,1,0,1) + (539136._dp)*Hr5(0,0,1,0,1) + z**5*(663552._dp)*Hr5(0,0,1,0,1) + &
                        z**2*(2.529792e6_dp)*Hr5(0,0,1,0,1) + z**4*(-2.5297919999999995e6_dp)*Hr5(0,0,1,1,0) + &
                        z**3*(-331776._dp)*Hr5(0,0,1,1,0) + nf*z**2*(-59843.95061728395_dp)*Hr5(0,0,1,1,0) + &
                        nf*z*(-6311.506172839506_dp)*Hr5(0,0,1,1,0) + nf*z**3*(6311.506172839506_dp)*Hr5(0,0,1,1,0) + &
                        nf*z**4*(59843.95061728395_dp)*Hr5(0,0,1,1,0) + z*(331776.00000000006_dp)*Hr5(0,0,1,1,0) + &
                        (539136._dp)*Hr5(0,0,1,1,0) + z**5*(663552._dp)*Hr5(0,0,1,1,0) + z**2*(2.529792e6_dp)*Hr5(0,0,1,1,0) + &
                        z**4*(-2.737152e6_dp)*Hr5(0,0,1,1,1) + z**3*(-995328.0000000002_dp)*Hr5(0,0,1,1,1) + &
                        nf*z*(-45631.20987654322_dp)*Hr5(0,0,1,1,1) + nf*z**2*(-45631.20987654322_dp)*Hr5(0,0,1,1,1) + &
                        nf*z**3*(45631.20987654322_dp)*Hr5(0,0,1,1,1) + nf*z**4*(45631.20987654322_dp)*Hr5(0,0,1,1,1) + &
                        (746496._dp)*Hr5(0,0,1,1,1) + z*(995328._dp)*Hr5(0,0,1,1,1) + z**5*(995328.0000000002_dp)*Hr5(0,0,1,1,1) + &
                        z**2*(2.7371519999999995e6_dp)*Hr5(0,0,1,1,1) + z**4*(-601344._dp)*Hr5(0,1,0,0,0) + &
                        z**3*(-145152.0000000001_dp)*Hr5(0,1,0,0,0) + nf*z*(-7671.572016460906_dp)*Hr5(0,1,0,0,0) + &
                        nf*z**2*(-7671.572016460906_dp)*Hr5(0,1,0,0,0) + nf*z**3*(7671.572016460906_dp)*Hr5(0,1,0,0,0) + &
                        nf*z**4*(7671.572016460906_dp)*Hr5(0,1,0,0,0) + z*(145152.00000000003_dp)*Hr5(0,1,0,0,0) + &
                        (228096.00000000003_dp)*Hr5(0,1,0,0,0) + z**5*(228096.00000000003_dp)*Hr5(0,1,0,0,0) + &
                        z**2*(601344._dp)*Hr5(0,1,0,0,0) + z**4*(-1.472256e6_dp)*Hr5(0,1,0,0,1) + &
                        z**3*(-269567.9999999998_dp)*Hr5(0,1,0,0,1) + nf*z*(-17933.695473251028_dp)*Hr5(0,1,0,0,1) + &
                        nf*z**2*(-17933.695473251028_dp)*Hr5(0,1,0,0,1) + nf*z**3*(17933.695473251028_dp)*Hr5(0,1,0,0,1) + &
                        nf*z**4*(17933.695473251028_dp)*Hr5(0,1,0,0,1) + z*(269568._dp)*Hr5(0,1,0,0,1) + (601344._dp)*Hr5(0,1,0,0,1) &
                        + z**5*(601344._dp)*Hr5(0,1,0,0,1) + z**2*(1.472256e6_dp)*Hr5(0,1,0,0,1) + &
                        z**4*(-1.741824e6_dp)*Hr5(0,1,0,1,0) + z**3*(-248832.00000000006_dp)*Hr5(0,1,0,1,0) + &
                        nf*z*(-20524.246913580246_dp)*Hr5(0,1,0,1,0) + nf*z**2*(-20524.246913580246_dp)*Hr5(0,1,0,1,0) + &
                        nf*z**3*(20524.246913580246_dp)*Hr5(0,1,0,1,0) + nf*z**4*(20524.246913580246_dp)*Hr5(0,1,0,1,0) + &
                        z*(248832.00000000017_dp)*Hr5(0,1,0,1,0) + (746496._dp)*Hr5(0,1,0,1,0) + z**5*(746496._dp)*Hr5(0,1,0,1,0) + &
                        z**2*(1.741824e6_dp)*Hr5(0,1,0,1,0) + z**4*(-1.9906560000000005e6_dp)*Hr5(0,1,0,1,1) + &
                        nf*z*(-20524.246913580246_dp)*Hr5(0,1,0,1,1) + nf*z**2*(-20524.246913580246_dp)*Hr5(0,1,0,1,1) + &
                        z*(-2.3283064365386963e-10_dp)*Hr5(0,1,0,1,1) + nf*z**3*(20524.246913580246_dp)*Hr5(0,1,0,1,1) + &
                        nf*z**4*(20524.246913580246_dp)*Hr5(0,1,0,1,1) + (995328.0000000002_dp)*Hr5(0,1,0,1,1) + &
                        z**5*(995328.0000000002_dp)*Hr5(0,1,0,1,1) + z**2*(1.9906560000000005e6_dp)*Hr5(0,1,0,1,1) + &
                        z**4*(-1.472256e6_dp)*Hr5(0,1,1,0,0) + z**3*(-269567.9999999998_dp)*Hr5(0,1,1,0,0) + &
                        nf*z*(-17933.695473251028_dp)*Hr5(0,1,1,0,0) + nf*z**2*(-17933.695473251028_dp)*Hr5(0,1,1,0,0) + &
                        nf*z**3*(17933.695473251028_dp)*Hr5(0,1,1,0,0) + nf*z**4*(17933.695473251028_dp)*Hr5(0,1,1,0,0) + &
                        z*(269568._dp)*Hr5(0,1,1,0,0) + (601344._dp)*Hr5(0,1,1,0,0) + z**5*(601344._dp)*Hr5(0,1,1,0,0) + &
                        z**2*(1.472256e6_dp)*Hr5(0,1,1,0,0) + z**4*(-1.9906560000000005e6_dp)*Hr5(0,1,1,0,1) + &
                        nf*z*(-20524.246913580246_dp)*Hr5(0,1,1,0,1) + nf*z**2*(-20524.246913580246_dp)*Hr5(0,1,1,0,1) + &
                        z*(-2.3283064365386963e-10_dp)*Hr5(0,1,1,0,1) + nf*z**3*(20524.246913580246_dp)*Hr5(0,1,1,0,1) + &
                        nf*z**4*(20524.246913580246_dp)*Hr5(0,1,1,0,1) + (995328.0000000002_dp)*Hr5(0,1,1,0,1) + &
                        z**5*(995328.0000000002_dp)*Hr5(0,1,1,0,1) + z**2*(1.9906560000000005e6_dp)*Hr5(0,1,1,0,1) + &
                        z**4*(-1.9906560000000005e6_dp)*Hr5(0,1,1,1,0) + nf*z*(-20524.246913580246_dp)*Hr5(0,1,1,1,0) + &
                        nf*z**2*(-20524.246913580246_dp)*Hr5(0,1,1,1,0) + z*(-2.3283064365386963e-10_dp)*Hr5(0,1,1,1,0) + &
                        nf*z**3*(20524.246913580246_dp)*Hr5(0,1,1,1,0) + nf*z**4*(20524.246913580246_dp)*Hr5(0,1,1,1,0) + &
                        (995328.0000000002_dp)*Hr5(0,1,1,1,0) + z**5*(995328.0000000002_dp)*Hr5(0,1,1,1,0) + &
                        z**2*(1.9906560000000005e6_dp)*Hr5(0,1,1,1,0) + z*(-995328.0000000005_dp)*Hr5(0,1,1,1,1) + &
                        z**4*(-995328.0000000002_dp)*Hr5(0,1,1,1,1) + z**2*(995328._dp)*Hr5(0,1,1,1,1) + &
                        (995328.0000000002_dp)*Hr5(0,1,1,1,1) + z**3*(995328.0000000002_dp)*Hr5(0,1,1,1,1) + &
                        z**5*(995328.0000000002_dp)*Hr5(0,1,1,1,1) + z**4*(-103680._dp)*Hr5(1,0,0,0,0) + &
                        z*(-103679.99999999999_dp)*Hr5(1,0,0,0,0) + z**3*(103679.99999999999_dp)*Hr5(1,0,0,0,0) + &
                        (103680._dp)*Hr5(1,0,0,0,0) + z**2*(103680._dp)*Hr5(1,0,0,0,0) + z**5*(103680._dp)*Hr5(1,0,0,0,0) + &
                        z**4*(-373248._dp)*Hr5(1,0,0,0,1) + z*(-373247.9999999999_dp)*Hr5(1,0,0,0,1) + &
                        z**3*(373247.9999999999_dp)*Hr5(1,0,0,0,1) + (373248._dp)*Hr5(1,0,0,0,1) + z**2*(373248._dp)*Hr5(1,0,0,0,1) &
                        + z**5*(373248._dp)*Hr5(1,0,0,0,1) + z*(-601344._dp)*Hr5(1,0,0,1,0) + z**4*(-601344._dp)*Hr5(1,0,0,1,0) + &
                        (601344._dp)*Hr5(1,0,0,1,0) + z**2*(601344._dp)*Hr5(1,0,0,1,0) + z**5*(601344._dp)*Hr5(1,0,0,1,0) + &
                        z**3*(601344.0000000001_dp)*Hr5(1,0,0,1,0) + z**4*(-870912.0000000001_dp)*Hr5(1,0,0,1,1) + &
                        z*(-870912._dp)*Hr5(1,0,0,1,1) + z**2*(870912._dp)*Hr5(1,0,0,1,1) + z**3*(870912._dp)*Hr5(1,0,0,1,1) + &
                        (870912.0000000001_dp)*Hr5(1,0,0,1,1) + z**5*(870912.0000000001_dp)*Hr5(1,0,0,1,1) + &
                        z*(-601344._dp)*Hr5(1,0,1,0,0) + z**4*(-601344._dp)*Hr5(1,0,1,0,0) + (601344._dp)*Hr5(1,0,1,0,0) + &
                        z**2*(601344._dp)*Hr5(1,0,1,0,0) + z**5*(601344._dp)*Hr5(1,0,1,0,0) + &
                        z**3*(601344.0000000001_dp)*Hr5(1,0,1,0,0) + z*(-995328.0000000005_dp)*Hr5(1,0,1,0,1) + &
                        z**4*(-995328.0000000002_dp)*Hr5(1,0,1,0,1) + z**2*(995328._dp)*Hr5(1,0,1,0,1) + &
                        (995328.0000000002_dp)*Hr5(1,0,1,0,1) + z**3*(995328.0000000002_dp)*Hr5(1,0,1,0,1) + &
                        z**5*(995328.0000000002_dp)*Hr5(1,0,1,0,1) + z*(-995328.0000000005_dp)*Hr5(1,0,1,1,0) + &
                        z**4*(-995328.0000000002_dp)*Hr5(1,0,1,1,0) + z**2*(995328._dp)*Hr5(1,0,1,1,0) + &
                        (995328.0000000002_dp)*Hr5(1,0,1,1,0) + z**3*(995328.0000000002_dp)*Hr5(1,0,1,1,0) + &
                        z**5*(995328.0000000002_dp)*Hr5(1,0,1,1,0) + z*(-995328.0000000005_dp)*Hr5(1,0,1,1,1) + &
                        z**4*(-995328.0000000002_dp)*Hr5(1,0,1,1,1) + z**2*(995328._dp)*Hr5(1,0,1,1,1) + &
                        (995328.0000000002_dp)*Hr5(1,0,1,1,1) + z**3*(995328.0000000002_dp)*Hr5(1,0,1,1,1) + &
                        z**5*(995328.0000000002_dp)*Hr5(1,0,1,1,1) + z**4*(-373248._dp)*Hr5(1,1,0,0,0) + &
                        z*(-373247.9999999999_dp)*Hr5(1,1,0,0,0) + z**3*(373247.9999999999_dp)*Hr5(1,1,0,0,0) + &
                        (373248._dp)*Hr5(1,1,0,0,0) + z**2*(373248._dp)*Hr5(1,1,0,0,0) + z**5*(373248._dp)*Hr5(1,1,0,0,0) + &
                        z**4*(-870912.0000000001_dp)*Hr5(1,1,0,0,1) + z*(-870912._dp)*Hr5(1,1,0,0,1) + &
                        z**2*(870912._dp)*Hr5(1,1,0,0,1) + z**3*(870912._dp)*Hr5(1,1,0,0,1) + (870912.0000000001_dp)*Hr5(1,1,0,0,1) &
                        + z**5*(870912.0000000001_dp)*Hr5(1,1,0,0,1) + z*(-995328.0000000005_dp)*Hr5(1,1,0,1,0) + &
                        z**4*(-995328.0000000002_dp)*Hr5(1,1,0,1,0) + z**2*(995328._dp)*Hr5(1,1,0,1,0) + &
                        (995328.0000000002_dp)*Hr5(1,1,0,1,0) + z**3*(995328.0000000002_dp)*Hr5(1,1,0,1,0) + &
                        z**5*(995328.0000000002_dp)*Hr5(1,1,0,1,0) + z*(-995328.0000000005_dp)*Hr5(1,1,0,1,1) + &
                        z**4*(-995328.0000000002_dp)*Hr5(1,1,0,1,1) + z**2*(995328._dp)*Hr5(1,1,0,1,1) + &
                        (995328.0000000002_dp)*Hr5(1,1,0,1,1) + z**3*(995328.0000000002_dp)*Hr5(1,1,0,1,1) + &
                        z**5*(995328.0000000002_dp)*Hr5(1,1,0,1,1) + z**4*(-870912.0000000001_dp)*Hr5(1,1,1,0,0) + &
                        z*(-870912._dp)*Hr5(1,1,1,0,0) + z**2*(870912._dp)*Hr5(1,1,1,0,0) + z**3*(870912._dp)*Hr5(1,1,1,0,0) + &
                        (870912.0000000001_dp)*Hr5(1,1,1,0,0) + z**5*(870912.0000000001_dp)*Hr5(1,1,1,0,0) + &
                        z*(-995328.0000000005_dp)*Hr5(1,1,1,0,1) + z**4*(-995328.0000000002_dp)*Hr5(1,1,1,0,1) + &
                        z**2*(995328._dp)*Hr5(1,1,1,0,1) + (995328.0000000002_dp)*Hr5(1,1,1,0,1) + &
                        z**3*(995328.0000000002_dp)*Hr5(1,1,1,0,1) + z**5*(995328.0000000002_dp)*Hr5(1,1,1,0,1) + &
                        z*(-995328.0000000005_dp)*Hr5(1,1,1,1,0) + z**4*(-995328.0000000002_dp)*Hr5(1,1,1,1,0) + &
                        z**2*(995328._dp)*Hr5(1,1,1,1,0) + (995328.0000000002_dp)*Hr5(1,1,1,1,0) + &
                        z**3*(995328.0000000002_dp)*Hr5(1,1,1,1,0) + z**5*(995328.0000000002_dp)*Hr5(1,1,1,1,0))/(z*(z + &
                        (1._dp))*((-1._dp) + z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = nf**2*(-3989.475668182756_dp) + nf**4*(-0.9382716049382716_dp) + &
                            nf**3*(92.9910297326987_dp) + (2445.2333096632065_dp) + nf*(50986.413601230975_dp)

                        else
                            Ibar_select = ((-2445.233309663221_dp) + nf**4*(z*(-0.9382716049382716_dp) + (0.9382716049382716_dp)) + &
                            nf**3*((-92.9910297326987_dp) + z*(92.9910297326987_dp)) + z*(2445.233309663221_dp) + &
                            nf**2*(z*(-3989.475668182756_dp) + (3989.475668182756_dp)) + nf*((-50986.413601230975_dp) + &
                            z*(50986.413601230975_dp)))/((-1._dp) + z*(1._dp))

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-163404.47958903562_dp) + nf**3*(-45.86797232641008_dp) + nf**2*(5327.088519468091_dp) &
                            + (916929.4811926486_dp)

                        else
                            Ibar_select = ((-916929.4811926486_dp) + nf**3*(z*(-45.86797232641008_dp) + (45.86797232641008_dp)) + &
                            nf**2*((-5327.088519468091_dp) + z*(5327.088519468091_dp)) + nf*(z*(-163404.47958903562_dp) + &
                            (163404.47958903562_dp)) + z*(916929.4811926486_dp))/((-1._dp) + z*(1._dp))

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (-950128.2643144309_dp) + nf**3*(-47.407407407407405_dp) + nf**2*(1593.3884831215662_dp) + &
                            nf*(73567.66782931326_dp)

                        else
                            Ibar_select = (-950128.2643144309_dp) + nf**3*(-47.407407407407405_dp) + nf**2*(1593.3884831215662_dp) + &
                            nf*(73567.66782931326_dp)

                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-222096.9653975472_dp) + nf**2*(-1920._dp) + nf*(56153.65449364699_dp)
                        else
                            Ibar_select = (z*(-222096.9653975472_dp) + nf**2*(z*(-1920._dp) + (1920._dp)) + &
                            nf*((-56153.65449364699_dp) + z*(56153.65449364699_dp)) + (222096.9653975472_dp))/((-1._dp) + z*(1._dp))

                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-15360._dp) + (155243.9123542583_dp)
                        else
                            Ibar_select = nf*(-15360._dp) + (155243.9123542583_dp)
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
                        Ibar_select = (z*(-1.63451709595519e6_dp) + z**3*(-189061.60522967763_dp) + nf*z**2*(-127319.58383841898_dp) &
                        + nf*z**4*(-84575.10375143704_dp) + nf*(-71844.4459345784_dp) + nf**2*z*(-6065.161495692442_dp) + &
                        nf**2*z**3*(-3585.8188976549654_dp) + nf**3*(-23.57933241883859_dp) + nf**3*z**4*(-23.57933241883859_dp) + &
                        nf**3*z**3*(-0.6054088775707491_dp) + nf**4*z**2*(-0.5135802469135803_dp) + nf**4*(-0.01316872427983539_dp) &
                        + nf**4*z**4*(-0.01316872427983539_dp) + nf**4*z**3*(0.21069958847736625_dp) + &
                        nf**4*z*(0.3292181069958848_dp) + nf**3*z**2*(9.40246913580247_dp) + nf**3*z*(38.36160457944547_dp) + &
                        nf**2*(2453.9443410726553_dp) + nf**2*z**4*(2825.947646966847_dp) + nf**2*z**2*(4371.088405307904_dp) + &
                        nf*z**3*(82553.7737411744_dp) + nf*z*(198415.82067838032_dp) + z**4*(461186.1963162266_dp) + &
                        (588148.78733184_dp) + z**2*(847994.930985109_dp) + z*(-719688.6876282968_dp)*Hr1(0) + &
                        z**3*(-667118.874304156_dp)*Hr1(0) + nf*z**4*(-37806.74927726053_dp)*Hr1(0) + &
                        nf*(-37240.01081960046_dp)*Hr1(0) + nf*z**2*(-10278.255634812002_dp)*Hr1(0) + &
                        nf**2*z**3*(-1244.311411240329_dp)*Hr1(0) + nf**2*z**2*(-526.0431005017153_dp)*Hr1(0) + &
                        nf**3*z*(-22.132235939643344_dp)*Hr1(0) + nf**3*(-7.141838134430727_dp)*Hr1(0) + &
                        nf**4*z**3*(-0.1580246913580247_dp)*Hr1(0) + nf**4*z*(0.1580246913580247_dp)*Hr1(0) + &
                        nf**3*z**4*(1.2159122085048013_dp)*Hr1(0) + nf**3*z**2*(12.444444444444443_dp)*Hr1(0) + &
                        nf**3*z**3*(17.391495198902607_dp)*Hr1(0) + nf**2*z*(240.63526162773235_dp)*Hr1(0) + &
                        nf**2*z**4*(445.1826245999086_dp)*Hr1(0) + nf**2*(996.5366255144032_dp)*Hr1(0) + &
                        nf*z*(30622.613275172458_dp)*Hr1(0) + nf*z**3*(59944.33054651886_dp)*Hr1(0) + &
                        z**2*(269149.5553309225_dp)*Hr1(0) + (391385.6335721348_dp)*Hr1(0) + z**4*(643289.6335721348_dp)*Hr1(0) + &
                        z**3*(-1.7139066551092213e6_dp)*Hr1(1) + z*(-1.4106298068776452e6_dp)*Hr1(1) + &
                        nf*z**2*(-135972.61566766168_dp)*Hr1(1) + nf*(-70674.56596283437_dp)*Hr1(1) + &
                        nf*z**4*(-70674.56596283437_dp)*Hr1(1) + nf**2*z**3*(-2636.100594421582_dp)*Hr1(1) + &
                        nf**2*z*(-2460.100594421582_dp)*Hr1(1) + nf**3*z*(-8.296296296296296_dp)*Hr1(1) + &
                        nf**3*(-5.925925925925926_dp)*Hr1(1) + nf**3*z**4*(-5.925925925925926_dp)*Hr1(1) + &
                        nf**3*z**3*(-4.7407407407407405_dp)*Hr1(1) + nf**3*z**2*(24.88888888888889_dp)*Hr1(1) + &
                        nf**2*(1441.7192501143118_dp)*Hr1(1) + nf**2*z**4*(1441.7192501143118_dp)*Hr1(1) + &
                        nf**2*z**2*(2212.7626886145404_dp)*Hr1(1) + nf*z*(131523.98166163775_dp)*Hr1(1) + &
                        nf*z**3*(145797.7659316927_dp)*Hr1(1) + (822147.9747444456_dp)*Hr1(1) + z**4*(822147.9747444456_dp)*Hr1(1) + &
                        z**2*(1.5300922163858097e6_dp)*Hr1(1) + z**3*(-466709.1775186815_dp)*Hr2(0,0) + z*(-342864._dp)*Hr2(0,0) + &
                        nf*z**2*(-49431.37820894347_dp)*Hr2(0,0) + nf*(-9973.333333333334_dp)*Hr2(0,0) + &
                        nf*z**4*(-6577.909465020576_dp)*Hr2(0,0) + nf**2*z**3*(-1245.8450023864796_dp)*Hr2(0,0) + &
                        nf**2*z*(-720.9478645545355_dp)*Hr2(0,0) + nf**3*z**2*(-22.439506172839504_dp)*Hr2(0,0) + &
                        nf**3*z*(4.108641975308641_dp)*Hr2(0,0) + nf**3*z**3*(18.330864197530865_dp)*Hr2(0,0) + &
                        nf**2*z**4*(36.49492455418381_dp)*Hr2(0,0) + nf**2*(149.57037037037037_dp)*Hr2(0,0) + &
                        nf**2*z**2*(1812.7275720164607_dp)*Hr2(0,0) + nf*z*(25584.16086075511_dp)*Hr2(0,0) + &
                        nf*z**3*(39342.46014654226_dp)*Hr2(0,0) + (150345.3235949176_dp)*Hr2(0,0) + &
                        z**4*(188937.3235949176_dp)*Hr2(0,0) + z**2*(461947.85392376385_dp)*Hr2(0,0) + &
                        z**3*(-1.0972157662780222e6_dp)*Hr2(0,1) + z*(-713088._dp)*Hr2(0,1) + &
                        nf*z**2*(-47496.82304526749_dp)*Hr2(0,1) + nf*(-27312.592592592595_dp)*Hr2(0,1) + &
                        nf*z**4*(-20521.744855967077_dp)*Hr2(0,1) + nf**2*z**3*(-899.1604938271605_dp)*Hr2(0,1) + &
                        nf**3*z*(-14.222222222222221_dp)*Hr2(0,1) + nf**3*z**3*(14.222222222222221_dp)*Hr2(0,1) + &
                        nf**2*z**2*(67.81893004115238_dp)*Hr2(0,1) + nf**2*z**4*(129.7119341563786_dp)*Hr2(0,1) + &
                        nf**2*(355.8628257887517_dp)*Hr2(0,1) + nf**2*z*(441.7668038408779_dp)*Hr2(0,1) + &
                        nf*z*(19367.861481941745_dp)*Hr2(0,1) + nf*z**3*(72795.29901188542_dp)*Hr2(0,1) + &
                        (419167.94156950555_dp)*Hr2(0,1) + z**4*(496351.94156950555_dp)*Hr2(0,1) + &
                        z**2*(818591.8247085167_dp)*Hr2(0,1) + z*(-973370.5887593407_dp)*Hr2(1,0) + &
                        z**3*(-973370.5887593407_dp)*Hr2(1,0) + nf*z**2*(-47496.82304526749_dp)*Hr2(1,0) + &
                        nf*(-23917.168724279836_dp)*Hr2(1,0) + nf*z**4*(-23917.168724279836_dp)*Hr2(1,0) + &
                        nf**2*z*(-228.69684499314133_dp)*Hr2(1,0) + nf**2*z**3*(-228.69684499314127_dp)*Hr2(1,0) + &
                        nf**2*z**2*(67.81893004115238_dp)*Hr2(1,0) + nf**2*(242.78737997256516_dp)*Hr2(1,0) + &
                        nf**2*z**4*(242.78737997256516_dp)*Hr2(1,0) + nf*z**3*(46081.580246913574_dp)*Hr2(1,0) + &
                        nf*z*(46081.58024691358_dp)*Hr2(1,0) + z**4*(491869.29437967035_dp)*Hr2(1,0) + &
                        (491869.2943796704_dp)*Hr2(1,0) + z**2*(920919.883139011_dp)*Hr2(1,0) + (z**3*(-1.810303766278022e6_dp) + &
                        z*(-1.657919649417033e6_dp) + nf**2*(z*(-649.3936899862827_dp) + z**3*(-457.39368998628254_dp) + &
                        z**2*(135.63786008230477_dp) + (485.5747599451303_dp) + z**4*(485.5747599451303_dp)) + &
                        nf*(z**2*(-94993.64609053498_dp) + (-47834.33744855967_dp) + z**4*(-47834.33744855967_dp) + &
                        z**3*(92163.16049382715_dp) + z*(98499.16049382716_dp)) + (915519.8831390111_dp) + &
                        z**4*(915519.8831390111_dp) + z**2*(1.6371836494170334e6_dp))*Hr2(1,1) + (z**3*(-259200.00000000006_dp) + &
                        z*(-24191.999999999993_dp) + nf**3*(z*(-4.424691358024691_dp) + z**3*(4.424691358024691_dp)) + &
                        nf**2*z*(z**2*(-685.0897119341564_dp) + (188.89218106995884_dp) + z*(496.1975308641976_dp)) + &
                        nf*(z**2*(-18014.551440329218_dp) + z*(-2267.865020576132_dp) + (-985.5999999999999_dp) + &
                        z**4*(-584.9898491083676_dp) + z**3*(21277.00631001372_dp)) + (28512._dp) + z**4*(28512._dp) + &
                        z**2*(235872.00000000003_dp))*Hr3(0,0,0) + z**3*(-566784._dp)*Hr3(0,0,1) + &
                        z*(-214272.00000000003_dp)*Hr3(0,0,1) + nf*z**2*(-39038.41975308642_dp)*Hr3(0,0,1) + &
                        nf*(-3669.333333333334_dp)*Hr3(0,0,1) + nf*z**4*(-2371.4238683127573_dp)*Hr3(0,0,1) + &
                        nf**2*z**3*(-716.3786008230453_dp)*Hr3(0,0,1) + nf**2*z*(-276.0164609053498_dp)*Hr3(0,0,1) + &
                        nf**2*z**2*(992.3950617283951_dp)*Hr3(0,0,1) + nf*z*(12161.58024691358_dp)*Hr3(0,0,1) + &
                        nf*z**3*(30613.596707818924_dp)*Hr3(0,0,1) + (114048._dp)*Hr3(0,0,1) + z**4*(114048._dp)*Hr3(0,0,1) + &
                        z**2*(590976._dp)*Hr3(0,0,1) + z**3*(-497664._dp)*Hr3(0,1,0) + z*(-145151.99999999997_dp)*Hr3(0,1,0) + &
                        nf*z**2*(-9027.950617283952_dp)*Hr3(0,1,0) + nf*(-4998.320987654321_dp)*Hr3(0,1,0) + &
                        nf*z**4*(-3700.4115226337453_dp)*Hr3(0,1,0) + nf*z*(-2090.666666666667_dp)*Hr3(0,1,0) + &
                        nf**2*z**3*(-220.18106995884773_dp)*Hr3(0,1,0) + nf**2*z*(220.18106995884773_dp)*Hr3(0,1,0) + &
                        nf*z**3*(16361.349794238682_dp)*Hr3(0,1,0) + (171072._dp)*Hr3(0,1,0) + z**4*(171072._dp)*Hr3(0,1,0) + &
                        z**2*(357696._dp)*Hr3(0,1,0) + z**3*(-995328._dp)*Hr3(0,1,1) + z*(-290303.99999999994_dp)*Hr3(0,1,1) + &
                        nf*z**2*(-18055.901234567904_dp)*Hr3(0,1,1) + nf*(-9996.641975308641_dp)*Hr3(0,1,1) + &
                        nf*z**4*(-7400.823045267491_dp)*Hr3(0,1,1) + nf*z*(-4181.333333333334_dp)*Hr3(0,1,1) + &
                        nf**2*z**3*(-440.36213991769546_dp)*Hr3(0,1,1) + nf**2*z*(440.36213991769546_dp)*Hr3(0,1,1) + &
                        nf*z**3*(32722.699588477364_dp)*Hr3(0,1,1) + (342144._dp)*Hr3(0,1,1) + z**4*(342144._dp)*Hr3(0,1,1) + &
                        z**2*(715392._dp)*Hr3(0,1,1) + z*(-214272.00000000003_dp)*Hr3(1,0,0) + z**3*(-214272._dp)*Hr3(1,0,0) + &
                        nf*z**2*(-6018.633744855968_dp)*Hr3(1,0,0) + nf*(-2899.577503429355_dp)*Hr3(1,0,0) + &
                        nf*z**4*(-2899.577503429355_dp)*Hr3(1,0,0) + nf*z*(4756.894375857339_dp)*Hr3(1,0,0) + &
                        nf*z**3*(4756.894375857339_dp)*Hr3(1,0,0) + (114048._dp)*Hr3(1,0,0) + z**4*(114048._dp)*Hr3(1,0,0) + &
                        z**2*(238464.00000000003_dp)*Hr3(1,0,0) + z*(-642816._dp)*Hr3(1,0,1) + z**3*(-642816._dp)*Hr3(1,0,1) + &
                        nf*z**2*(-18055.9012345679_dp)*Hr3(1,0,1) + nf*(-8698.732510288066_dp)*Hr3(1,0,1) + &
                        nf*z**4*(-8698.732510288066_dp)*Hr3(1,0,1) + nf*z*(14270.683127572014_dp)*Hr3(1,0,1) + &
                        nf*z**3*(14270.683127572018_dp)*Hr3(1,0,1) + (342144._dp)*Hr3(1,0,1) + z**4*(342144._dp)*Hr3(1,0,1) + &
                        z**2*(715392._dp)*Hr3(1,0,1) + z*(-642816._dp)*Hr3(1,1,0) + z**3*(-642816._dp)*Hr3(1,1,0) + &
                        nf*z**2*(-18055.9012345679_dp)*Hr3(1,1,0) + nf*(-8698.732510288066_dp)*Hr3(1,1,0) + &
                        nf*z**4*(-8698.732510288066_dp)*Hr3(1,1,0) + nf*z*(14270.683127572014_dp)*Hr3(1,1,0) + &
                        nf*z**3*(14270.683127572018_dp)*Hr3(1,1,0) + (342144._dp)*Hr3(1,1,0) + z**4*(342144._dp)*Hr3(1,1,0) + &
                        z**2*(715392._dp)*Hr3(1,1,0) + z*(-1.513728e6_dp)*Hr3(1,1,1) + z**3*(-1.285632e6_dp)*Hr3(1,1,1) + &
                        nf*z**2*(-36111.8024691358_dp)*Hr3(1,1,1) + nf*(-17397.46502057613_dp)*Hr3(1,1,1) + &
                        nf*z**4*(-17397.46502057613_dp)*Hr3(1,1,1) + nf*z**3*(28541.366255144036_dp)*Hr3(1,1,1) + &
                        nf*z*(42365.36625514403_dp)*Hr3(1,1,1) + (684288._dp)*Hr3(1,1,1) + z**4*(684288._dp)*Hr3(1,1,1) + &
                        z**2*(1.430784e6_dp)*Hr3(1,1,1) + z**3*(-66355.2_dp)*Hr4(0,0,0,0) + &
                        nf*z**2*(-5882.31111111111_dp)*Hr4(0,0,0,0) + z*(-4147.200000000001_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**3*(-73.32345679012346_dp)*Hr4(0,0,0,0) + nf**2*z*(-29.076543209876544_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**2*(102.39999999999999_dp)*Hr4(0,0,0,0) + nf*z*(525.4847736625514_dp)*Hr4(0,0,0,0) + &
                        (2073.6_dp)*Hr4(0,0,0,0) + z**4*(2073.6_dp)*Hr4(0,0,0,0) + nf*z**3*(5356.82633744856_dp)*Hr4(0,0,0,0) + &
                        z**2*(68428.8_dp)*Hr4(0,0,0,0) + z**3*(-165888._dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(-6940.444444444445_dp)*Hr4(0,0,0,1) + nf*z*(-1807.8024691358023_dp)*Hr4(0,0,0,1) + &
                        nf**2*z**3*(-44.24691358024691_dp)*Hr4(0,0,0,1) + nf**2*z*(44.24691358024691_dp)*Hr4(0,0,0,1) + &
                        nf*z**3*(8748.246913580248_dp)*Hr4(0,0,0,1) + (10368._dp)*Hr4(0,0,0,1) + z**4*(10368._dp)*Hr4(0,0,0,1) + &
                        z**2*(155520._dp)*Hr4(0,0,0,1) + z**3*(-165888._dp)*Hr4(0,0,1,0) + z*(-41471.99999999999_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(-6940.444444444443_dp)*Hr4(0,0,1,0) + nf*z*(2130.1728395061727_dp)*Hr4(0,0,1,0) + &
                        nf*z**3*(4810.271604938272_dp)*Hr4(0,0,1,0) + (20736._dp)*Hr4(0,0,1,0) + z**4*(20736._dp)*Hr4(0,0,1,0) + &
                        z**2*(186624._dp)*Hr4(0,0,1,0) + z**3*(-331776._dp)*Hr4(0,0,1,1) + z*(-82943.99999999999_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(-13880.888888888887_dp)*Hr4(0,0,1,1) + nf*z*(4260.3456790123455_dp)*Hr4(0,0,1,1) + &
                        nf*z**3*(9620.543209876543_dp)*Hr4(0,0,1,1) + (41472._dp)*Hr4(0,0,1,1) + z**4*(41472._dp)*Hr4(0,0,1,1) + &
                        z**2*(373248._dp)*Hr4(0,0,1,1) + z**3*(-82944._dp)*Hr4(0,1,0,0) + nf*z*(-893.366255144033_dp)*Hr4(0,1,0,0) + &
                        nf*z**3*(893.366255144033_dp)*Hr4(0,1,0,0) + (20736._dp)*Hr4(0,1,0,0) + z**4*(20736._dp)*Hr4(0,1,0,0) + &
                        z**2*(62208._dp)*Hr4(0,1,0,0) + z**3*(-248832._dp)*Hr4(0,1,0,1) + nf*z*(-2680.098765432099_dp)*Hr4(0,1,0,1) &
                        + z*(7.275957614183426e-12_dp)*Hr4(0,1,0,1) + nf*z**3*(2680.098765432099_dp)*Hr4(0,1,0,1) + &
                        (62208._dp)*Hr4(0,1,0,1) + z**4*(62208._dp)*Hr4(0,1,0,1) + z**2*(186624._dp)*Hr4(0,1,0,1) + &
                        z**3*(-248832._dp)*Hr4(0,1,1,0) + nf*z*(-2680.098765432099_dp)*Hr4(0,1,1,0) + &
                        z*(7.275957614183426e-12_dp)*Hr4(0,1,1,0) + nf*z**3*(2680.098765432099_dp)*Hr4(0,1,1,0) + &
                        (62208._dp)*Hr4(0,1,1,0) + z**4*(62208._dp)*Hr4(0,1,1,0) + z**2*(186624._dp)*Hr4(0,1,1,0) + &
                        z**3*(-497664._dp)*Hr4(0,1,1,1) + nf*z*(-5360.197530864198_dp)*Hr4(0,1,1,1) + &
                        z*(1.4551915228366852e-11_dp)*Hr4(0,1,1,1) + nf*z**3*(5360.197530864198_dp)*Hr4(0,1,1,1) + &
                        (124416._dp)*Hr4(0,1,1,1) + z**4*(124416._dp)*Hr4(0,1,1,1) + z**2*(373248._dp)*Hr4(0,1,1,1) + &
                        z**3*(-20736._dp)*Hr4(1,0,0,0) + z*(-20735.999999999996_dp)*Hr4(1,0,0,0) + (10368._dp)*Hr4(1,0,0,0) + &
                        z**4*(10368._dp)*Hr4(1,0,0,0) + z**2*(31104._dp)*Hr4(1,0,0,0) + z**3*(-82944._dp)*Hr4(1,0,0,1) + &
                        z*(-82943.99999999999_dp)*Hr4(1,0,0,1) + (41472._dp)*Hr4(1,0,0,1) + z**4*(41472._dp)*Hr4(1,0,0,1) + &
                        z**2*(124416._dp)*Hr4(1,0,0,1) + z**3*(-124416._dp)*Hr4(1,0,1,0) + z*(-124415.99999999997_dp)*Hr4(1,0,1,0) + &
                        (62208._dp)*Hr4(1,0,1,0) + z**4*(62208._dp)*Hr4(1,0,1,0) + z**2*(186623.99999999997_dp)*Hr4(1,0,1,0) + &
                        z**3*(-248832._dp)*Hr4(1,0,1,1) + z*(-248831.99999999994_dp)*Hr4(1,0,1,1) + (124416._dp)*Hr4(1,0,1,1) + &
                        z**4*(124416._dp)*Hr4(1,0,1,1) + z**2*(373247.99999999994_dp)*Hr4(1,0,1,1) + z**3*(-82944._dp)*Hr4(1,1,0,0) &
                        + z*(-82943.99999999999_dp)*Hr4(1,1,0,0) + (41472._dp)*Hr4(1,1,0,0) + z**4*(41472._dp)*Hr4(1,1,0,0) + &
                        z**2*(124416._dp)*Hr4(1,1,0,0) + z**3*(-248832._dp)*Hr4(1,1,0,1) + z*(-248831.99999999994_dp)*Hr4(1,1,0,1) + &
                        (124416._dp)*Hr4(1,1,0,1) + z**4*(124416._dp)*Hr4(1,1,0,1) + z**2*(373247.99999999994_dp)*Hr4(1,1,0,1) + &
                        z**3*(-248832._dp)*Hr4(1,1,1,0) + z*(-248831.99999999994_dp)*Hr4(1,1,1,0) + (124416._dp)*Hr4(1,1,1,0) + &
                        z**4*(124416._dp)*Hr4(1,1,1,0) + z**2*(373247.99999999994_dp)*Hr4(1,1,1,0) + &
                        z*(-746495.9999999999_dp)*Hr4(1,1,1,1) + z**3*(-497664._dp)*Hr4(1,1,1,1) + (248832._dp)*Hr4(1,1,1,1) + &
                        z**4*(248832._dp)*Hr4(1,1,1,1) + z**2*(746495.9999999999_dp)*Hr4(1,1,1,1))/(z*((-1._dp) + z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-784.8273595169849_dp) + nf**3*(-2.9243272299524024_dp) + nf**2*(67.82255608042989_dp) &
                            + (38020.00857894558_dp)

                        else
                            Ibar_select = ((-38020.008578945584_dp) + nf**3*(z*(-2.9243272299524024_dp) + (2.9243272299524024_dp)) + &
                            nf**2*((-67.8225560804299_dp) + z*(67.8225560804299_dp)) + nf*(z*(-784.8273595169848_dp) + &
                            (784.8273595169848_dp)) + z*(38020.008578945584_dp))/((-1._dp) + z*(1._dp))

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-12879.829333534592_dp) + nf**3*(-7.822222222222222_dp) + &
                            nf**4*(0.11851851851851852_dp) + nf**2*(351.5136704174297_dp) + (126339.09630149989_dp)

                        else
                            Ibar_select = ((-126339.0963014999_dp) + nf**4*((-0.11851851851851852_dp) + z*(0.11851851851851852_dp)) &
                            + nf**3*(z*(-7.822222222222225_dp) + (7.822222222222225_dp)) + nf**2*((-351.5136704174298_dp) + &
                            z*(351.5136704174298_dp)) + nf*(z*(-12879.82933353459_dp) + (12879.82933353459_dp)) + &
                            z*(126339.0963014999_dp))/((-1._dp) + z*(1._dp))

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (-303276.8482315759_dp) + nf**2*(-176._dp) + nf**3*(3.5555555555555554_dp) + &
                            nf*(14273.78427005494_dp)

                        else
                            Ibar_select = (-303276.8482315759_dp) + nf**2*(-176._dp) + nf**3*(3.5555555555555554_dp) + &
                            nf*(14273.78427005494_dp)

                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = nf**2*(-96._dp) + nf*(3168._dp) + (76192.05843049446_dp)
                        else
                            Ibar_select = ((-76192.05843049448_dp) + nf**2*(z*(-96._dp) + (96._dp)) + nf*((-3168._dp) + &
                            z*(3168._dp)) + z*(76192.05843049448_dp))/((-1._dp) + z*(1._dp))

                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-2304._dp) + (38016._dp)
                        else
                            Ibar_select = nf*(-2304._dp) + (38016._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (-10368._dp)
                        else
                            Ibar_select = (-10368._dp)
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
                        Ibar_select = (z*(-2.1881388304818448e7_dp) + nf*z**2*(-2.476525489938832e6_dp) + &
                        nf*z**4*(-1.5115659312363262e6_dp) + nf*(-1.422823727192426e6_dp) + nf**2*z*(-218719.94593548917_dp) + &
                        nf**2*z**3*(-106500.99467388017_dp) + nf**3*z**2*(-2797.349598781436_dp) + nf**3*z**4*(-2053.62548280808_dp) &
                        + nf**3*(-1778.3841134159727_dp) + nf**4*z**2*(-20.979240969364422_dp) + nf**4*z*(-12.818441466198026_dp) + &
                        nf**5*z*(-0.14631915866483763_dp) + nf**5*z**3*(-0.0936442615454961_dp) + nf**5*(0.005852766346593506_dp) + &
                        nf**5*z**4*(0.005852766346593506_dp) + nf**5*z**2*(0.22825788751714676_dp) + &
                        nf**4*z**3*(7.559080596411915_dp) + nf**4*(13.119300919575268_dp) + nf**4*z**4*(13.119300919575268_dp) + &
                        nf**3*z**3*(2240.0399155751347_dp) + nf**3*z*(4389.319279430354_dp) + nf**2*(79655.52870939834_dp) + &
                        nf**2*z**4*(97494.05296380713_dp) + nf**2*z**2*(149609.99177220813_dp) + nf*z**3*(1.2400203863026542e6_dp) + &
                        z**3*(1.727327365206176e6_dp) + z**4*(2.1880952548908265e6_dp) + nf*z*(4.092066060257676e6_dp) + &
                        (8.415452519661546e6_dp) + z**2*(1.0287716917612791e7_dp) + z*(-1.211266511228937e7_dp)*Hr1(0) + &
                        z**3*(-1.0074060168700002e7_dp)*Hr1(0) + nf*z**4*(-1.129980175238819e6_dp)*Hr1(0) + &
                        nf*(-789495.651853655_dp)*Hr1(0) + nf*z**2*(-237966.6961954545_dp)*Hr1(0) + &
                        nf**2*z**3*(-52259.046851133186_dp)*Hr1(0) + nf**2*z*(-16052.58298198683_dp)*Hr1(0) + &
                        nf**2*z**2*(-4250.437182527716_dp)*Hr1(0) + nf**3*(-681.8477671086724_dp)*Hr1(0) + &
                        nf**3*z*(-357.3242349723403_dp)*Hr1(0) + nf**3*z**4*(-219.3387593354671_dp)*Hr1(0) + &
                        nf**4*z**3*(-13.793019356805365_dp)*Hr1(0) + nf**4*z**2*(-6.838957475994513_dp)*Hr1(0) + &
                        nf**4*z**4*(-0.9647309861301631_dp)*Hr1(0) + nf**5*z*(-0.07023319615912206_dp)*Hr1(0) + &
                        nf**5*z**3*(0.07023319615912206_dp)*Hr1(0) + nf**4*(3.9203779911598846_dp)*Hr1(0) + &
                        nf**4*z*(16.649169333942996_dp)*Hr1(0) + nf**3*z**2*(592.3403251079803_dp)*Hr1(0) + &
                        nf**3*z**3*(733.9630289010923_dp)*Hr1(0) + nf**2*z**4*(31623.69101619516_dp)*Hr1(0) + &
                        nf**2*(37154.99372722018_dp)*Hr1(0) + nf*z*(992406.4799372845_dp)*Hr1(0) + &
                        nf*z**3*(1.2654375176396053e6_dp)*Hr1(0) + z**2*(4.41077545766107e6_dp)*Hr1(0) + &
                        (5.590625655543536e6_dp)*Hr1(0) + z**4*(1.1471717309449855e7_dp)*Hr1(0) + &
                        z**3*(-2.3971806965120688e7_dp)*Hr1(1) + z*(-2.065674195714444e7_dp)*Hr1(1) + &
                        nf*z**2*(-3.2111751141209053e6_dp)*Hr1(1) + nf*(-1.6604595806288861e6_dp)*Hr1(1) + &
                        nf*z**4*(-1.6604595806288861e6_dp)*Hr1(1) + nf**2*z**3*(-135228.68482571546_dp)*Hr1(1) + &
                        nf**2*z*(-125556.40467568496_dp)*Hr1(1) + nf**3*(-901.1865264441395_dp)*Hr1(1) + &
                        nf**3*z**4*(-901.1865264441395_dp)*Hr1(1) + nf**3*z**2*(-799.8761164456638_dp)*Hr1(1) + &
                        nf**4*z**2*(-13.677914951989026_dp)*Hr1(1) + nf**4*z**3*(2.856149977137632_dp)*Hr1(1) + &
                        nf**4*(2.9556470050297214_dp)*Hr1(1) + nf**4*z**4*(2.9556470050297214_dp)*Hr1(1) + &
                        nf**4*z*(4.910470964791953_dp)*Hr1(1) + nf**3*z*(1233.3319920743788_dp)*Hr1(1) + &
                        nf**3*z**3*(1368.9171772595641_dp)*Hr1(1) + nf**2*(66126.64005738891_dp)*Hr1(1) + &
                        nf**2*z**4*(66126.64005738891_dp)*Hr1(1) + nf**2*z**2*(128531.8093866226_dp)*Hr1(1) + &
                        nf*z*(3.0685900260487497e6_dp)*Hr1(1) + nf*z**3*(3.413652545442093e6_dp)*Hr1(1) + &
                        (1.1584186575898275e7_dp)*Hr1(1) + z**4*(1.1584186575898275e7_dp)*Hr1(1) + &
                        z**2*(2.2619374703231364e7_dp)*Hr1(1) + z*(-6.615458215992291e6_dp)*Hr2(0,0) + &
                        z**3*(-6.416788040181599e6_dp)*Hr2(0,0) + nf*z**2*(-894552.2757606346_dp)*Hr2(0,0) + &
                        nf*(-269152.771440413_dp)*Hr2(0,0) + nf*z**4*(-261301.05541731758_dp)*Hr2(0,0) + &
                        nf**2*z**3*(-35674.2846918305_dp)*Hr2(0,0) + nf**2*z*(-28418.31509448524_dp)*Hr2(0,0) + &
                        nf**3*z**2*(-1381.9756744398721_dp)*Hr2(0,0) + nf**3*(-96.07803688462124_dp)*Hr2(0,0) + &
                        nf**3*z**4*(-12.414692882182592_dp)*Hr2(0,0) + nf**4*z**3*(-8.427983539094651_dp)*Hr2(0,0) + &
                        nf**4*z*(-0.9130315500685864_dp)*Hr2(0,0) + nf**4*z**2*(9.341015089163236_dp)*Hr2(0,0) + &
                        nf**3*z*(546.9313665220056_dp)*Hr2(0,0) + nf**3*z**3*(932.8703710180035_dp)*Hr2(0,0) + &
                        nf**2*z**4*(4325.372351775644_dp)*Hr2(0,0) + nf**2*(9153.303703703703_dp)*Hr2(0,0) + &
                        nf**2*z**2*(51141.923730836395_dp)*Hr2(0,0) + nf*z*(678567.2186434069_dp)*Hr2(0,0) + &
                        nf*z**3*(754781.5603800406_dp)*Hr2(0,0) + (2.4518415360716265e6_dp)*Hr2(0,0) + &
                        z**4*(4.175521536071626e6_dp)*Hr2(0,0) + z**2*(6.141486001014076e6_dp)*Hr2(0,0) + &
                        z**3*(-1.7437280211893618e7_dp)*Hr2(0,1) + z*(-1.4166558363946155e7_dp)*Hr2(0,1) + &
                        nf*z**2*(-1.3782108353958824e6_dp)*Hr2(0,1) + nf*(-726365.8519008246_dp)*Hr2(0,1) + &
                        nf*z**4*(-719506.4967932728_dp)*Hr2(0,1) + nf**2*z**3*(-58062.14412762053_dp)*Hr2(0,1) + &
                        nf**2*z*(-9958.4039595613_dp)*Hr2(0,1) + nf**3*z*(-478.6860539551898_dp)*Hr2(0,1) + &
                        nf**3*(-215.12135345221765_dp)*Hr2(0,1) + nf**3*z**4*(-47.79466544734035_dp)*Hr2(0,1) + &
                        nf**4*z**3*(-7.514951989026063_dp)*Hr2(0,1) + nf**4*z*(7.514951989026063_dp)*Hr2(0,1) + &
                        nf**3*z**2*(53.18701417466849_dp)*Hr2(0,1) + nf**3*z**3*(656.4150586800793_dp)*Hr2(0,1) + &
                        nf**2*z**4*(13849.287364730986_dp)*Hr2(0,1) + nf**2*(23505.150068587107_dp)*Hr2(0,1) + &
                        nf**2*z**2*(32250.110653863743_dp)*Hr2(0,1) + nf*z*(972164.5585268016_dp)*Hr2(0,1) + &
                        nf*z**3*(1.9281106839936727e6_dp)*Hr2(0,1) + (6.014383527098954e6_dp)*Hr2(0,1) + &
                        z**4*(9.461743527098954e6_dp)*Hr2(0,1) + z**2*(1.4283936334211698e7_dp)*Hr2(0,1) + &
                        z**3*(-1.8591012266388934e7_dp)*Hr2(1,0) + z*(-1.8591012266388927e7_dp)*Hr2(1,0) + &
                        nf*z**2*(-1.472873161106781e6_dp)*Hr2(1,0) + nf*(-762156.0156369442_dp)*Hr2(1,0) + &
                        nf*z**4*(-762156.0156369442_dp)*Hr2(1,0) + nf**2*z*(-34010.27404359092_dp)*Hr2(1,0) + &
                        nf**2*z**3*(-34010.27404359092_dp)*Hr2(1,0) + nf**3*(-131.45800944977898_dp)*Hr2(1,0) + &
                        nf**3*z**4*(-131.45800944977898_dp)*Hr2(1,0) + nf**3*z**2*(53.18701417466854_dp)*Hr2(1,0) + &
                        nf**3*z*(88.86450236244474_dp)*Hr2(1,0) + nf**3*z**3*(88.86450236244474_dp)*Hr2(1,0) + &
                        nf**2*(18677.218716659045_dp)*Hr2(1,0) + nf**2*z**4*(18677.218716659045_dp)*Hr2(1,0) + &
                        nf**2*z**2*(32250.110653863743_dp)*Hr2(1,0) + nf*z*(1.5196339490004994e6_dp)*Hr2(1,0) + &
                        nf*z**3*(1.5196339490004994e6_dp)*Hr2(1,0) + (9.200828721953806e6_dp)*Hr2(1,0) + &
                        z**4*(9.200828721953806e6_dp)*Hr2(1,0) + z**2*(1.7648951334471308e7_dp)*Hr2(1,0) + &
                        (z**3*(-3.1603838575839777e7_dp) + z*(-2.731806775432542e7_dp) + nf**3*((-262.91601889955797_dp) + &
                        z**4*(-262.91601889955797_dp) + z**2*(106.37402834933708_dp) + z**3*(177.72900472488948_dp) + &
                        z*(241.72900472488948_dp)) + nf**2*(z*(-71188.54808718184_dp) + z**3*(-68020.54808718184_dp) + &
                        (37354.43743331809_dp) + z**4*(37354.43743331809_dp) + z**2*(64500.22130772749_dp)) + &
                        nf*(z**2*(-2.7564216707917647e6_dp) + (-1.4458723486940975e6_dp) + z**4*(-1.4458723486940975e6_dp) + &
                        z*(2.7478911256594853e6_dp) + z**3*(2.900275242520475e6_dp)) + (1.547612705419791e7_dp) + &
                        z**4*(1.547612705419791e7_dp) + z**2*(2.85678726684234e7_dp))*Hr2(1,1) + &
                        z**3*(-3.802431252089671e6_dp)*Hr3(0,0,0) + z*(-779027.1532556044_dp)*Hr3(0,0,0) + &
                        nf*z**2*(-389473.99067238206_dp)*Hr3(0,0,0) + nf*(-46464._dp)*Hr3(0,0,0) + &
                        nf*z**4*(-32681.437585733882_dp)*Hr3(0,0,0) + nf**2*z**3*(-23124.72441104388_dp)*Hr3(0,0,0) + &
                        nf*z*(-12144.194616135988_dp)*Hr3(0,0,0) + nf**3*z**2*(-301.6164609053498_dp)*Hr3(0,0,0) + &
                        nf**3*z*(-185.91897576588934_dp)*Hr3(0,0,0) + nf**4*z**3*(-3.230727023319616_dp)*Hr3(0,0,0) + &
                        nf**4*z*(3.230727023319616_dp)*Hr3(0,0,0) + nf**2*z**4*(326.9921048620637_dp)*Hr3(0,0,0) + &
                        nf**3*z**3*(487.5354366712391_dp)*Hr3(0,0,0) + nf**2*(835.7925925925925_dp)*Hr3(0,0,0) + &
                        nf**2*z*(4756.902885401481_dp)*Hr3(0,0,0) + nf**2*z**2*(17525.036828187745_dp)*Hr3(0,0,0) + &
                        nf*z**3*(470203.6228742519_dp)*Hr3(0,0,0) + (615219.1766278022_dp)*Hr3(0,0,0) + &
                        z**4*(754150.3766278022_dp)*Hr3(0,0,0) + z**2*(3.258277628717473e6_dp)*Hr3(0,0,0) + &
                        z**3*(-8.996606130224176e6_dp)*Hr3(0,0,1) + z*(-4.543935532556044e6_dp)*Hr3(0,0,1) + &
                        nf*z**2*(-1.0210500154819434e6_dp)*Hr3(0,0,1) + nf*(-176260.74074074073_dp)*Hr3(0,0,1) + &
                        nf*z**4*(-131714.63374485596_dp)*Hr3(0,0,1) + nf**2*z**3*(-29661.936899862827_dp)*Hr3(0,0,1) + &
                        nf**2*z*(-13099.895747599452_dp)*Hr3(0,0,1) + nf**3*z**2*(-603.2329218106996_dp)*Hr3(0,0,1) + &
                        nf**3*z*(152.4294467306813_dp)*Hr3(0,0,1) + nf**3*z**3*(450.8034750800183_dp)*Hr3(0,0,1) + &
                        nf**2*z**4*(1333.9390946502058_dp)*Hr3(0,0,1) + nf**2*(2960.5399176954734_dp)*Hr3(0,0,1) + &
                        nf**2*z**2*(39747.3536351166_dp)*Hr3(0,0,1) + nf*z*(367853.72853350476_dp)*Hr3(0,0,1) + &
                        nf*z**3*(918931.6614340354_dp)*Hr3(0,0,1) + (2.261599766278022e6_dp)*Hr3(0,0,1) + &
                        z**4*(2.724703766278022e6_dp)*Hr3(0,0,1) + z**2*(8.493405896502199e6_dp)*Hr3(0,0,1) + &
                        z**3*(-1.0348159065112088e7_dp)*Hr3(0,1,0) + z*(-5.076864e6_dp)*Hr3(0,1,0) + &
                        nf*z**2*(-462714.2642889804_dp)*Hr3(0,1,0) + nf*(-246799.80246913582_dp)*Hr3(0,1,0) + &
                        nf*z**4*(-202253.69547325102_dp)*Hr3(0,1,0) + nf**2*z**3*(-12379.361682670324_dp)*Hr3(0,1,0) + &
                        nf**3*z*(-149.1870141746685_dp)*Hr3(0,1,0) + nf**3*z**3*(149.1870141746685_dp)*Hr3(0,1,0) + &
                        nf**2*z**4*(2140.064014631916_dp)*Hr3(0,1,0) + nf**2*(3766.6648376771827_dp)*Hr3(0,1,0) + &
                        nf**2*z*(4182.679469593049_dp)*Hr3(0,1,0) + nf**2*z**2*(4209.953360768176_dp)*Hr3(0,1,0) + &
                        nf*z*(140999.18194582223_dp)*Hr3(0,1,0) + nf*z**3*(707408.580285545_dp)*Hr3(0,1,0) + &
                        (3.712831766278022e6_dp)*Hr3(0,1,0) + z**4*(4.175935766278022e6_dp)*Hr3(0,1,0) + &
                        z**2*(7.649663298834067e6_dp)*Hr3(0,1,0) + z**3*(-1.9059069195336267e7_dp)*Hr3(0,1,1) + &
                        z*(-1.0153728e7_dp)*Hr3(0,1,1) + nf*z**2*(-925428.5285779608_dp)*Hr3(0,1,1) + &
                        nf*(-493599.60493827163_dp)*Hr3(0,1,1) + nf*z**4*(-404507.39094650204_dp)*Hr3(0,1,1) + &
                        nf**2*z**3*(-24758.723365340647_dp)*Hr3(0,1,1) + nf**3*z*(-298.374028349337_dp)*Hr3(0,1,1) + &
                        nf**3*z**3*(298.374028349337_dp)*Hr3(0,1,1) + nf**2*z**4*(4280.128029263832_dp)*Hr3(0,1,1) + &
                        nf**2*(7533.3296753543655_dp)*Hr3(0,1,1) + nf**2*z*(8365.358939186099_dp)*Hr3(0,1,1) + &
                        nf**2*z**2*(8419.906721536352_dp)*Hr3(0,1,1) + nf*z*(297329.82933083654_dp)*Hr3(0,1,1) + &
                        nf*z**3*(1.399485695131898e6_dp)*Hr3(0,1,1) + (7.016351298834067e6_dp)*Hr3(0,1,1) + &
                        z**4*(7.942559298834068e6_dp)*Hr3(0,1,1) + z**2*(1.40713898965022e7_dp)*Hr3(0,1,1) + &
                        z**3*(-5.278111766278023e6_dp)*Hr3(1,0,0) + z*(-5.278111766278021e6_dp)*Hr3(1,0,0) + &
                        nf*z**2*(-308476.17619265354_dp)*Hr3(1,0,0) + nf*(-149684.49931412895_dp)*Hr3(1,0,0) + &
                        nf*z**4*(-149684.49931412895_dp)*Hr3(1,0,0) + nf**2*z*(-2732.227404359092_dp)*Hr3(1,0,0) + &
                        nf**2*z**3*(-2732.227404359092_dp)*Hr3(1,0,0) + nf**2*(1968.9096174363667_dp)*Hr3(1,0,0) + &
                        nf**2*z**4*(1968.9096174363667_dp)*Hr3(1,0,0) + nf**2*z**2*(2806.6355738454504_dp)*Hr3(1,0,0) + &
                        nf*z*(282802.58741045574_dp)*Hr3(1,0,0) + nf*z**3*(282802.58741045574_dp)*Hr3(1,0,0) + &
                        (2.697807883139011e6_dp)*Hr3(1,0,0) + z**4*(2.697807883139011e6_dp)*Hr3(1,0,0) + &
                        z**2*(5.304431649417033e6_dp)*Hr3(1,0,0) + z*(-1.4606398597668132e7_dp)*Hr3(1,0,1) + &
                        z**3*(-1.4606398597668132e7_dp)*Hr3(1,0,1) + nf*z**2*(-925428.5285779606_dp)*Hr3(1,0,1) + &
                        nf*(-449053.49794238684_dp)*Hr3(1,0,1) + nf*z**4*(-449053.49794238684_dp)*Hr3(1,0,1) + &
                        nf**2*z*(-8196.682213077274_dp)*Hr3(1,0,1) + nf**2*z**3*(-8196.682213077274_dp)*Hr3(1,0,1) + &
                        nf**2*(5906.728852309099_dp)*Hr3(1,0,1) + nf**2*z**4*(5906.728852309099_dp)*Hr3(1,0,1) + &
                        nf**2*z**2*(8419.90672153635_dp)*Hr3(1,0,1) + nf*z*(848407.762231367_dp)*Hr3(1,0,1) + &
                        nf*z**3*(848407.7622313672_dp)*Hr3(1,0,1) + (7.479455298834067e6_dp)*Hr3(1,0,1) + &
                        z**4*(7.479455298834067e6_dp)*Hr3(1,0,1) + z**2*(1.40713898965022e7_dp)*Hr3(1,0,1) + &
                        z*(-1.5425023065112086e7_dp)*Hr3(1,1,0) + z**3*(-1.5425023065112086e7_dp)*Hr3(1,1,0) + &
                        nf*z**2*(-925428.5285779606_dp)*Hr3(1,1,0) + nf*(-449053.49794238684_dp)*Hr3(1,1,0) + &
                        nf*z**4*(-449053.49794238684_dp)*Hr3(1,1,0) + nf**2*z*(-8196.682213077274_dp)*Hr3(1,1,0) + &
                        nf**2*z**3*(-8196.682213077274_dp)*Hr3(1,1,0) + nf**2*(5906.728852309099_dp)*Hr3(1,1,0) + &
                        nf**2*z**4*(5906.728852309099_dp)*Hr3(1,1,0) + nf**2*z**2*(8419.90672153635_dp)*Hr3(1,1,0) + &
                        nf*z*(848407.762231367_dp)*Hr3(1,1,0) + nf*z**3*(848407.7622313672_dp)*Hr3(1,1,0) + &
                        (7.888767532556043e6_dp)*Hr3(1,1,0) + z**4*(7.888767532556043e6_dp)*Hr3(1,1,0) + &
                        z**2*(1.529932659766813e7_dp)*Hr3(1,1,0) + z**3*(-2.9212797195336264e7_dp)*Hr3(1,1,1) + &
                        z*(-2.8847803793004397e7_dp)*Hr3(1,1,1) + nf*z**2*(-1.8508570571559211e6_dp)*Hr3(1,1,1) + &
                        nf*(-898106.9958847737_dp)*Hr3(1,1,1) + nf*z**4*(-898106.9958847737_dp)*Hr3(1,1,1) + &
                        nf**2*z*(-24073.36442615455_dp)*Hr3(1,1,1) + nf**2*z**3*(-16393.36442615455_dp)*Hr3(1,1,1) + &
                        nf**2*(11813.457704618198_dp)*Hr3(1,1,1) + nf**2*z**4*(11813.457704618198_dp)*Hr3(1,1,1) + &
                        nf**2*z**2*(16839.8134430727_dp)*Hr3(1,1,1) + nf*z**3*(1.6968155244627343e6_dp)*Hr3(1,1,1) + &
                        nf*z*(1.950255524462734e6_dp)*Hr3(1,1,1) + (1.4958910597668134e7_dp)*Hr3(1,1,1) + &
                        z**4*(1.4958910597668134e7_dp)*Hr3(1,1,1) + z**2*(2.81427797930044e7_dp)*Hr3(1,1,1) + &
                        z**3*(-1.5510528e6_dp)*Hr4(0,0,0,0) + nf*z**2*(-190489.27224508458_dp)*Hr4(0,0,0,0) + &
                        z*(-151372.80000000002_dp)*Hr4(0,0,0,0) + nf**2*z**3*(-7137.144215820759_dp)*Hr4(0,0,0,0) + &
                        nf*(-3097.6_dp)*Hr4(0,0,0,0) + nf*z**4*(-2081.5909769852155_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(-1152.4096936442616_dp)*Hr4(0,0,0,0) + nf**3*z**2*(-120.09876543209877_dp)*Hr4(0,0,0,0) + &
                        nf**3*z*(32.4477366255144_dp)*Hr4(0,0,0,0) + nf**3*z**3*(87.65102880658436_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**2*(8289.55390946502_dp)*Hr4(0,0,0,0) + nf*z*(16617.384910836765_dp)*Hr4(0,0,0,0) + &
                        (79833.6_dp)*Hr4(0,0,0,0) + z**4*(79833.6_dp)*Hr4(0,0,0,0) + nf*z**3*(176977.47831123305_dp)*Hr4(0,0,0,0) + &
                        z**2*(1.5769727999999998e6_dp)*Hr4(0,0,0,0) + z**3*(-4.0227840000000005e6_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(-290897.7631458619_dp)*Hr4(0,0,0,1) + z*(-290304._dp)*Hr4(0,0,0,1) + &
                        nf*z*(-34379.336808413354_dp)*Hr4(0,0,0,1) + nf*(-14828.088888888888_dp)*Hr4(0,0,0,1) + &
                        nf*z**4*(-10435.923304374333_dp)*Hr4(0,0,0,1) + nf**2*z**3*(-10235.013443072703_dp)*Hr4(0,0,0,1) + &
                        nf**3*z*(-55.20329218106996_dp)*Hr4(0,0,0,1) + nf**3*z**3*(55.20329218106996_dp)*Hr4(0,0,0,1) + &
                        nf**2*z*(2767.8200274348424_dp)*Hr4(0,0,0,1) + nf**2*z**2*(7467.193415637862_dp)*Hr4(0,0,0,1) + &
                        nf*z**3*(340173.11214753846_dp)*Hr4(0,0,0,1) + (399167.99999999994_dp)*Hr4(0,0,0,1) + &
                        z**4*(399167.99999999994_dp)*Hr4(0,0,0,1) + z**2*(3.6858239999999995e6_dp)*Hr4(0,0,0,1) + &
                        z**3*(-4.313088e6_dp)*Hr4(0,0,1,0) + z*(-1.513728e6_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(-319671.6598079561_dp)*Hr4(0,0,1,0) + nf*(-28093.62962962963_dp)*Hr4(0,0,1,0) + &
                        nf*z**4*(-21013.18701417467_dp)*Hr4(0,0,1,0) + nf**2*z**3*(-5283.6433470507545_dp)*Hr4(0,0,1,0) + &
                        nf**2*z*(-2183.5500685871057_dp)*Hr4(0,0,1,0) + nf**2*z**2*(7467.19341563786_dp)*Hr4(0,0,1,0) + &
                        nf*z*(93272.8449931413_dp)*Hr4(0,0,1,0) + nf*z**3*(254769.63145861911_dp)*Hr4(0,0,1,0) + &
                        (798335.9999999999_dp)*Hr4(0,0,1,0) + z**4*(798335.9999999999_dp)*Hr4(0,0,1,0) + &
                        z**2*(4.572288e6_dp)*Hr4(0,0,1,0) + z**3*(-8.626176e6_dp)*Hr4(0,0,1,1) + z*(-3.027456e6_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(-639343.3196159122_dp)*Hr4(0,0,1,1) + nf*(-56187.25925925926_dp)*Hr4(0,0,1,1) + &
                        nf*z**4*(-42026.37402834934_dp)*Hr4(0,0,1,1) + nf**2*z**3*(-10567.286694101509_dp)*Hr4(0,0,1,1) + &
                        nf**2*z*(-4367.100137174211_dp)*Hr4(0,0,1,1) + nf**2*z**2*(14934.38683127572_dp)*Hr4(0,0,1,1) + &
                        nf*z*(186545.6899862826_dp)*Hr4(0,0,1,1) + nf*z**3*(509539.26291723823_dp)*Hr4(0,0,1,1) + &
                        (1.5966719999999998e6_dp)*Hr4(0,0,1,1) + z**4*(1.5966719999999998e6_dp)*Hr4(0,0,1,1) + &
                        z**2*(9.144576e6_dp)*Hr4(0,0,1,1) + z**3*(-2.446848e6_dp)*Hr4(0,1,0,0) + z*(-580608._dp)*Hr4(0,1,0,0) + &
                        nf*z**2*(-57547.79332418839_dp)*Hr4(0,1,0,0) + nf*(-26202.951989026064_dp)*Hr4(0,1,0,0) + &
                        nf*z**4*(-21482.65691205609_dp)*Hr4(0,1,0,0) + nf*z*(-11583.561042524005_dp)*Hr4(0,1,0,0) + &
                        nf**2*z**3*(-1033.3644261545496_dp)*Hr4(0,1,0,0) + nf**2*z*(1033.3644261545496_dp)*Hr4(0,1,0,0) + &
                        nf*z**3*(96080.96326779455_dp)*Hr4(0,1,0,0) + (798335.9999999999_dp)*Hr4(0,1,0,0) + &
                        z**4*(798335.9999999999_dp)*Hr4(0,1,0,0) + z**2*(1.772928e6_dp)*Hr4(0,1,0,0) + &
                        z**3*(-7.340544e6_dp)*Hr4(0,1,0,1) + z*(-1.7418239999999995e6_dp)*Hr4(0,1,0,1) + &
                        nf*z**2*(-172643.37997256516_dp)*Hr4(0,1,0,1) + nf*(-78608.85596707818_dp)*Hr4(0,1,0,1) + &
                        nf*z**4*(-64447.97073616827_dp)*Hr4(0,1,0,1) + nf*z*(-34750.683127572025_dp)*Hr4(0,1,0,1) + &
                        nf**2*z**3*(-3100.093278463649_dp)*Hr4(0,1,0,1) + nf**2*z*(3100.093278463649_dp)*Hr4(0,1,0,1) + &
                        nf*z**3*(288242.8898033836_dp)*Hr4(0,1,0,1) + (2.395008e6_dp)*Hr4(0,1,0,1) + &
                        z**4*(2.395008e6_dp)*Hr4(0,1,0,1) + z**2*(5.318784e6_dp)*Hr4(0,1,0,1) + z**3*(-7.340544e6_dp)*Hr4(0,1,1,0) + &
                        z*(-1.7418239999999995e6_dp)*Hr4(0,1,1,0) + nf*z**2*(-172643.37997256516_dp)*Hr4(0,1,1,0) + &
                        nf*(-78608.85596707818_dp)*Hr4(0,1,1,0) + nf*z**4*(-64447.97073616827_dp)*Hr4(0,1,1,0) + &
                        nf*z*(-34750.683127572025_dp)*Hr4(0,1,1,0) + nf**2*z**3*(-3100.093278463649_dp)*Hr4(0,1,1,0) + &
                        nf**2*z*(3100.093278463649_dp)*Hr4(0,1,1,0) + nf*z**3*(288242.8898033836_dp)*Hr4(0,1,1,0) + &
                        (2.395008e6_dp)*Hr4(0,1,1,0) + z**4*(2.395008e6_dp)*Hr4(0,1,1,0) + z**2*(5.318784e6_dp)*Hr4(0,1,1,0) + &
                        z**3*(-1.4681088e7_dp)*Hr4(0,1,1,1) + z*(-3.483647999999999e6_dp)*Hr4(0,1,1,1) + &
                        nf*z**2*(-345286.7599451303_dp)*Hr4(0,1,1,1) + nf*(-157217.71193415637_dp)*Hr4(0,1,1,1) + &
                        nf*z**4*(-128895.94147233653_dp)*Hr4(0,1,1,1) + nf*z*(-69501.36625514405_dp)*Hr4(0,1,1,1) + &
                        nf**2*z**3*(-6200.186556927298_dp)*Hr4(0,1,1,1) + nf**2*z*(6200.186556927298_dp)*Hr4(0,1,1,1) + &
                        nf*z**3*(576485.7796067672_dp)*Hr4(0,1,1,1) + (4.790016e6_dp)*Hr4(0,1,1,1) + &
                        z**4*(4.790016e6_dp)*Hr4(0,1,1,1) + z**2*(1.0637568e7_dp)*Hr4(0,1,1,1) + z*(-756864._dp)*Hr4(1,0,0,0) + &
                        z**3*(-756864._dp)*Hr4(1,0,0,0) + nf*z**2*(-28773.896662094194_dp)*Hr4(1,0,0,0) + &
                        nf*(-11921.402225270538_dp)*Hr4(1,0,0,0) + nf*z**4*(-11921.402225270538_dp)*Hr4(1,0,0,0) + &
                        nf*z*(21124.350556317633_dp)*Hr4(1,0,0,0) + nf*z**3*(21124.350556317637_dp)*Hr4(1,0,0,0) + &
                        (399167.99999999994_dp)*Hr4(1,0,0,0) + z**4*(399167.99999999994_dp)*Hr4(1,0,0,0) + &
                        z**2*(886464.0000000001_dp)*Hr4(1,0,0,0) + z*(-3.027456e6_dp)*Hr4(1,0,0,1) + &
                        z**3*(-3.027456e6_dp)*Hr4(1,0,0,1) + nf*z**2*(-115095.58664837677_dp)*Hr4(1,0,0,1) + &
                        nf*(-47685.60890108215_dp)*Hr4(1,0,0,1) + nf*z**4*(-47685.60890108215_dp)*Hr4(1,0,0,1) + &
                        nf*z*(84497.40222527053_dp)*Hr4(1,0,0,1) + nf*z**3*(84497.40222527055_dp)*Hr4(1,0,0,1) + &
                        (1.5966719999999998e6_dp)*Hr4(1,0,0,1) + z**4*(1.5966719999999998e6_dp)*Hr4(1,0,0,1) + &
                        z**2*(3.5458560000000005e6_dp)*Hr4(1,0,0,1) + z*(-4.541184e6_dp)*Hr4(1,0,1,0) + &
                        z**3*(-4.541184e6_dp)*Hr4(1,0,1,0) + nf*z**2*(-172643.37997256516_dp)*Hr4(1,0,1,0) + &
                        nf*(-71528.41335162323_dp)*Hr4(1,0,1,0) + nf*z**4*(-71528.41335162323_dp)*Hr4(1,0,1,0) + &
                        nf*z**3*(126746.10333790579_dp)*Hr4(1,0,1,0) + nf*z*(126746.1033379058_dp)*Hr4(1,0,1,0) + &
                        (2.395008e6_dp)*Hr4(1,0,1,0) + z**4*(2.395008e6_dp)*Hr4(1,0,1,0) + z**2*(5.318784e6_dp)*Hr4(1,0,1,0) + &
                        z*(-9.082368e6_dp)*Hr4(1,0,1,1) + z**3*(-9.082368e6_dp)*Hr4(1,0,1,1) + &
                        nf*z**2*(-345286.7599451303_dp)*Hr4(1,0,1,1) + nf*(-143056.82670324645_dp)*Hr4(1,0,1,1) + &
                        nf*z**4*(-143056.82670324645_dp)*Hr4(1,0,1,1) + nf*z**3*(253492.20667581158_dp)*Hr4(1,0,1,1) + &
                        nf*z*(253492.2066758116_dp)*Hr4(1,0,1,1) + (4.790016e6_dp)*Hr4(1,0,1,1) + z**4*(4.790016e6_dp)*Hr4(1,0,1,1) &
                        + z**2*(1.0637568e7_dp)*Hr4(1,0,1,1) + z*(-3.027456e6_dp)*Hr4(1,1,0,0) + z**3*(-3.027456e6_dp)*Hr4(1,1,0,0) &
                        + nf*z**2*(-115095.58664837677_dp)*Hr4(1,1,0,0) + nf*(-47685.60890108215_dp)*Hr4(1,1,0,0) + &
                        nf*z**4*(-47685.60890108215_dp)*Hr4(1,1,0,0) + nf*z*(84497.40222527053_dp)*Hr4(1,1,0,0) + &
                        nf*z**3*(84497.40222527055_dp)*Hr4(1,1,0,0) + (1.5966719999999998e6_dp)*Hr4(1,1,0,0) + &
                        z**4*(1.5966719999999998e6_dp)*Hr4(1,1,0,0) + z**2*(3.5458560000000005e6_dp)*Hr4(1,1,0,0) + &
                        z*(-9.082368e6_dp)*Hr4(1,1,0,1) + z**3*(-9.082368e6_dp)*Hr4(1,1,0,1) + &
                        nf*z**2*(-345286.7599451303_dp)*Hr4(1,1,0,1) + nf*(-143056.82670324645_dp)*Hr4(1,1,0,1) + &
                        nf*z**4*(-143056.82670324645_dp)*Hr4(1,1,0,1) + nf*z**3*(253492.20667581158_dp)*Hr4(1,1,0,1) + &
                        nf*z*(253492.2066758116_dp)*Hr4(1,1,0,1) + (4.790016e6_dp)*Hr4(1,1,0,1) + z**4*(4.790016e6_dp)*Hr4(1,1,0,1) &
                        + z**2*(1.0637568e7_dp)*Hr4(1,1,0,1) + z*(-9.082368e6_dp)*Hr4(1,1,1,0) + z**3*(-9.082368e6_dp)*Hr4(1,1,1,0) &
                        + nf*z**2*(-345286.7599451303_dp)*Hr4(1,1,1,0) + nf*(-143056.82670324645_dp)*Hr4(1,1,1,0) + &
                        nf*z**4*(-143056.82670324645_dp)*Hr4(1,1,1,0) + nf*z**3*(253492.20667581158_dp)*Hr4(1,1,1,0) + &
                        nf*z*(253492.2066758116_dp)*Hr4(1,1,1,0) + (4.790016e6_dp)*Hr4(1,1,1,0) + z**4*(4.790016e6_dp)*Hr4(1,1,1,0) &
                        + z**2*(1.0637568e7_dp)*Hr4(1,1,1,0) + z*(-2.2270464e7_dp)*Hr4(1,1,1,1) + &
                        z**3*(-1.8164736e7_dp)*Hr4(1,1,1,1) + nf*z**2*(-690573.5198902606_dp)*Hr4(1,1,1,1) + &
                        nf*(-286113.6534064929_dp)*Hr4(1,1,1,1) + nf*z**4*(-286113.6534064929_dp)*Hr4(1,1,1,1) + &
                        nf*z**3*(506984.41335162317_dp)*Hr4(1,1,1,1) + nf*z*(755816.4133516232_dp)*Hr4(1,1,1,1) + &
                        (9.580032e6_dp)*Hr4(1,1,1,1) + z**4*(9.580032e6_dp)*Hr4(1,1,1,1) + z**2*(2.1275136e7_dp)*Hr4(1,1,1,1) + &
                        z**3*(-265420.8_dp)*Hr5(0,0,0,0,0) + nf*z**2*(-25654.360493827164_dp)*Hr5(0,0,0,0,0) + &
                        nf*z*(-1302.4512117055326_dp)*Hr5(0,0,0,0,0) + nf**2*z**3*(-610.1860082304527_dp)*Hr5(0,0,0,0,0) + &
                        nf**3*z*(-1.68559670781893_dp)*Hr5(0,0,0,0,0) + nf**3*z**3*(1.68559670781893_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z*(109.56378600823045_dp)*Hr5(0,0,0,0,0) + nf**2*z**2*(500.62222222222226_dp)*Hr5(0,0,0,0,0) + &
                        (4147.2_dp)*Hr5(0,0,0,0,0) + z**4*(4147.2_dp)*Hr5(0,0,0,0,0) + &
                        nf*z**3*(26956.811705532695_dp)*Hr5(0,0,0,0,0) + z**2*(261273.60000000003_dp)*Hr5(0,0,0,0,0) + &
                        z**3*(-796262.4_dp)*Hr5(0,0,0,0,1) + nf*z**2*(-63609.362962962965_dp)*Hr5(0,0,0,0,1) + &
                        z*(-49766.40000000001_dp)*Hr5(0,0,0,0,1) + nf**2*z**3*(-711.3218106995885_dp)*Hr5(0,0,0,0,1) + &
                        nf**2*z*(-289.92263374485594_dp)*Hr5(0,0,0,0,1) + nf**2*z**2*(1001.2444444444443_dp)*Hr5(0,0,0,0,1) + &
                        nf*z*(5807.6298125285775_dp)*Hr5(0,0,0,0,1) + (24883.2_dp)*Hr5(0,0,0,0,1) + z**4*(24883.2_dp)*Hr5(0,0,0,0,1) &
                        + nf*z**3*(57801.73315043438_dp)*Hr5(0,0,0,0,1) + z**2*(821145.6_dp)*Hr5(0,0,0,0,1) + &
                        z**3*(-995327.9999999999_dp)*Hr5(0,0,0,1,0) + nf*z**2*(-36901.92592592593_dp)*Hr5(0,0,0,1,0) + &
                        nf*z*(-9876.192043895748_dp)*Hr5(0,0,0,1,0) + nf**2*z**3*(-210.69958847736626_dp)*Hr5(0,0,0,1,0) + &
                        z*(7.275957614183426e-12_dp)*Hr5(0,0,0,1,0) + nf**2*z*(210.69958847736626_dp)*Hr5(0,0,0,1,0) + &
                        nf*z**3*(46778.11796982167_dp)*Hr5(0,0,0,1,0) + (62207.99999999999_dp)*Hr5(0,0,0,1,0) + &
                        z**4*(62207.99999999999_dp)*Hr5(0,0,0,1,0) + z**2*(933119.9999999999_dp)*Hr5(0,0,0,1,0) + &
                        z**3*(-1.9906559999999998e6_dp)*Hr5(0,0,0,1,1) + nf*z**2*(-73803.85185185185_dp)*Hr5(0,0,0,1,1) + &
                        nf*z*(-19752.384087791495_dp)*Hr5(0,0,0,1,1) + nf**2*z**3*(-421.3991769547325_dp)*Hr5(0,0,0,1,1) + &
                        z*(1.4551915228366852e-11_dp)*Hr5(0,0,0,1,1) + nf**2*z*(421.3991769547325_dp)*Hr5(0,0,0,1,1) + &
                        nf*z**3*(93556.23593964335_dp)*Hr5(0,0,0,1,1) + (124415.99999999999_dp)*Hr5(0,0,0,1,1) + &
                        z**4*(124415.99999999999_dp)*Hr5(0,0,0,1,1) + z**2*(1.8662399999999998e6_dp)*Hr5(0,0,0,1,1) + &
                        z**3*(-663552._dp)*Hr5(0,0,1,0,0) + z*(-165888.00000000003_dp)*Hr5(0,0,1,0,0) + &
                        nf*z**2*(-24601.283950617286_dp)*Hr5(0,0,1,0,0) + nf*z*(7640.435299497028_dp)*Hr5(0,0,1,0,0) + &
                        nf*z**3*(16960.848651120257_dp)*Hr5(0,0,1,0,0) + (82944._dp)*Hr5(0,0,1,0,0) + &
                        z**4*(82944._dp)*Hr5(0,0,1,0,0) + z**2*(746496._dp)*Hr5(0,0,1,0,0) + &
                        z**3*(-1.9906559999999998e6_dp)*Hr5(0,0,1,0,1) + z*(-497663.9999999998_dp)*Hr5(0,0,1,0,1) + &
                        nf*z**2*(-73803.85185185185_dp)*Hr5(0,0,1,0,1) + nf*z*(22921.305898491082_dp)*Hr5(0,0,1,0,1) + &
                        nf*z**3*(50882.545953360765_dp)*Hr5(0,0,1,0,1) + (248831.99999999997_dp)*Hr5(0,0,1,0,1) + &
                        z**4*(248831.99999999997_dp)*Hr5(0,0,1,0,1) + z**2*(2.2394879999999995e6_dp)*Hr5(0,0,1,0,1) + &
                        z**3*(-1.9906559999999998e6_dp)*Hr5(0,0,1,1,0) + z*(-497663.9999999998_dp)*Hr5(0,0,1,1,0) + &
                        nf*z**2*(-73803.85185185185_dp)*Hr5(0,0,1,1,0) + nf*z*(22921.305898491082_dp)*Hr5(0,0,1,1,0) + &
                        nf*z**3*(50882.545953360765_dp)*Hr5(0,0,1,1,0) + (248831.99999999997_dp)*Hr5(0,0,1,1,0) + &
                        z**4*(248831.99999999997_dp)*Hr5(0,0,1,1,0) + z**2*(2.2394879999999995e6_dp)*Hr5(0,0,1,1,0) + &
                        z**3*(-3.9813119999999995e6_dp)*Hr5(0,0,1,1,1) + z*(-995327.9999999997_dp)*Hr5(0,0,1,1,1) + &
                        nf*z**2*(-147607.7037037037_dp)*Hr5(0,0,1,1,1) + nf*z*(45842.611796982164_dp)*Hr5(0,0,1,1,1) + &
                        nf*z**3*(101765.09190672153_dp)*Hr5(0,0,1,1,1) + (497663.99999999994_dp)*Hr5(0,0,1,1,1) + &
                        z**4*(497663.99999999994_dp)*Hr5(0,0,1,1,1) + z**2*(4.478975999999999e6_dp)*Hr5(0,0,1,1,1) + &
                        z**3*(-248831.99999999997_dp)*Hr5(0,1,0,0,0) + nf*z*(-2330.103337905807_dp)*Hr5(0,1,0,0,0) + &
                        z*(7.275957614183426e-12_dp)*Hr5(0,1,0,0,0) + nf*z**3*(2330.103337905807_dp)*Hr5(0,1,0,0,0) + &
                        (62207.99999999999_dp)*Hr5(0,1,0,0,0) + z**4*(62207.99999999999_dp)*Hr5(0,1,0,0,0) + &
                        z**2*(186624._dp)*Hr5(0,1,0,0,0) + z**3*(-995327.9999999999_dp)*Hr5(0,1,0,0,1) + &
                        nf*z*(-9320.413351623229_dp)*Hr5(0,1,0,0,1) + z*(2.9103830456733704e-11_dp)*Hr5(0,1,0,0,1) + &
                        nf*z**3*(9320.413351623229_dp)*Hr5(0,1,0,0,1) + (248831.99999999997_dp)*Hr5(0,1,0,0,1) + &
                        z**4*(248831.99999999997_dp)*Hr5(0,1,0,0,1) + z**2*(746496._dp)*Hr5(0,1,0,0,1) + &
                        z**3*(-1.492992e6_dp)*Hr5(0,1,0,1,0) + nf*z*(-13980.620027434843_dp)*Hr5(0,1,0,1,0) + &
                        nf*z**3*(13980.620027434843_dp)*Hr5(0,1,0,1,0) + (373248._dp)*Hr5(0,1,0,1,0) + &
                        z**4*(373248._dp)*Hr5(0,1,0,1,0) + z**2*(1.119744e6_dp)*Hr5(0,1,0,1,0) + &
                        z**3*(-2.985984e6_dp)*Hr5(0,1,0,1,1) + nf*z*(-27961.240054869686_dp)*Hr5(0,1,0,1,1) + &
                        nf*z**3*(27961.240054869686_dp)*Hr5(0,1,0,1,1) + (746496._dp)*Hr5(0,1,0,1,1) + &
                        z**4*(746496._dp)*Hr5(0,1,0,1,1) + z**2*(2.239488e6_dp)*Hr5(0,1,0,1,1) + &
                        z**3*(-995327.9999999999_dp)*Hr5(0,1,1,0,0) + nf*z*(-9320.413351623229_dp)*Hr5(0,1,1,0,0) + &
                        z*(2.9103830456733704e-11_dp)*Hr5(0,1,1,0,0) + nf*z**3*(9320.413351623229_dp)*Hr5(0,1,1,0,0) + &
                        (248831.99999999997_dp)*Hr5(0,1,1,0,0) + z**4*(248831.99999999997_dp)*Hr5(0,1,1,0,0) + &
                        z**2*(746496._dp)*Hr5(0,1,1,0,0) + z**3*(-2.985984e6_dp)*Hr5(0,1,1,0,1) + &
                        nf*z*(-27961.240054869686_dp)*Hr5(0,1,1,0,1) + nf*z**3*(27961.240054869686_dp)*Hr5(0,1,1,0,1) + &
                        (746496._dp)*Hr5(0,1,1,0,1) + z**4*(746496._dp)*Hr5(0,1,1,0,1) + z**2*(2.239488e6_dp)*Hr5(0,1,1,0,1) + &
                        z**3*(-2.985984e6_dp)*Hr5(0,1,1,1,0) + nf*z*(-27961.240054869686_dp)*Hr5(0,1,1,1,0) + &
                        nf*z**3*(27961.240054869686_dp)*Hr5(0,1,1,1,0) + (746496._dp)*Hr5(0,1,1,1,0) + &
                        z**4*(746496._dp)*Hr5(0,1,1,1,0) + z**2*(2.239488e6_dp)*Hr5(0,1,1,1,0) + &
                        z**3*(-5.971968e6_dp)*Hr5(0,1,1,1,1) + nf*z*(-55922.48010973937_dp)*Hr5(0,1,1,1,1) + &
                        nf*z**3*(55922.48010973937_dp)*Hr5(0,1,1,1,1) + (1.492992e6_dp)*Hr5(0,1,1,1,1) + &
                        z**4*(1.492992e6_dp)*Hr5(0,1,1,1,1) + z**2*(4.478976e6_dp)*Hr5(0,1,1,1,1) + &
                        z*(-49766.40000000001_dp)*Hr5(1,0,0,0,0) + z**3*(-49766.4_dp)*Hr5(1,0,0,0,0) + (24883.2_dp)*Hr5(1,0,0,0,0) + &
                        z**4*(24883.2_dp)*Hr5(1,0,0,0,0) + z**2*(74649.6_dp)*Hr5(1,0,0,0,0) + &
                        z**3*(-248831.99999999997_dp)*Hr5(1,0,0,0,1) + z*(-248831.9999999999_dp)*Hr5(1,0,0,0,1) + &
                        (124415.99999999999_dp)*Hr5(1,0,0,0,1) + z**4*(124415.99999999999_dp)*Hr5(1,0,0,0,1) + &
                        z**2*(373247.99999999994_dp)*Hr5(1,0,0,0,1) + z**3*(-497663.99999999994_dp)*Hr5(1,0,0,1,0) + &
                        z*(-497663.9999999998_dp)*Hr5(1,0,0,1,0) + (248831.99999999997_dp)*Hr5(1,0,0,1,0) + &
                        z**4*(248831.99999999997_dp)*Hr5(1,0,0,1,0) + z**2*(746495.9999999999_dp)*Hr5(1,0,0,1,0) + &
                        z**3*(-995327.9999999999_dp)*Hr5(1,0,0,1,1) + z*(-995327.9999999997_dp)*Hr5(1,0,0,1,1) + &
                        (497663.99999999994_dp)*Hr5(1,0,0,1,1) + z**4*(497663.99999999994_dp)*Hr5(1,0,0,1,1) + &
                        z**2*(1.4929919999999998e6_dp)*Hr5(1,0,0,1,1) + z**3*(-497663.99999999994_dp)*Hr5(1,0,1,0,0) + &
                        z*(-497663.9999999998_dp)*Hr5(1,0,1,0,0) + (248831.99999999997_dp)*Hr5(1,0,1,0,0) + &
                        z**4*(248831.99999999997_dp)*Hr5(1,0,1,0,0) + z**2*(746495.9999999999_dp)*Hr5(1,0,1,0,0) + &
                        z*(-1.492992e6_dp)*Hr5(1,0,1,0,1) + z**3*(-1.492992e6_dp)*Hr5(1,0,1,0,1) + (746496._dp)*Hr5(1,0,1,0,1) + &
                        z**4*(746496._dp)*Hr5(1,0,1,0,1) + z**2*(2.239488e6_dp)*Hr5(1,0,1,0,1) + z*(-1.492992e6_dp)*Hr5(1,0,1,1,0) + &
                        z**3*(-1.492992e6_dp)*Hr5(1,0,1,1,0) + (746496._dp)*Hr5(1,0,1,1,0) + z**4*(746496._dp)*Hr5(1,0,1,1,0) + &
                        z**2*(2.239488e6_dp)*Hr5(1,0,1,1,0) + z*(-2.985984e6_dp)*Hr5(1,0,1,1,1) + &
                        z**3*(-2.985984e6_dp)*Hr5(1,0,1,1,1) + (1.492992e6_dp)*Hr5(1,0,1,1,1) + z**4*(1.492992e6_dp)*Hr5(1,0,1,1,1) &
                        + z**2*(4.478976e6_dp)*Hr5(1,0,1,1,1) + z**3*(-248831.99999999997_dp)*Hr5(1,1,0,0,0) + &
                        z*(-248831.9999999999_dp)*Hr5(1,1,0,0,0) + (124415.99999999999_dp)*Hr5(1,1,0,0,0) + &
                        z**4*(124415.99999999999_dp)*Hr5(1,1,0,0,0) + z**2*(373247.99999999994_dp)*Hr5(1,1,0,0,0) + &
                        z**3*(-995327.9999999999_dp)*Hr5(1,1,0,0,1) + z*(-995327.9999999997_dp)*Hr5(1,1,0,0,1) + &
                        (497663.99999999994_dp)*Hr5(1,1,0,0,1) + z**4*(497663.99999999994_dp)*Hr5(1,1,0,0,1) + &
                        z**2*(1.4929919999999998e6_dp)*Hr5(1,1,0,0,1) + z*(-1.492992e6_dp)*Hr5(1,1,0,1,0) + &
                        z**3*(-1.492992e6_dp)*Hr5(1,1,0,1,0) + (746496._dp)*Hr5(1,1,0,1,0) + z**4*(746496._dp)*Hr5(1,1,0,1,0) + &
                        z**2*(2.239488e6_dp)*Hr5(1,1,0,1,0) + (z*(-1._dp) + (0.9999999999999998_dp) + &
                        z**2*(1._dp))**2*(1.4929919999999998e6_dp)*Hr5(1,1,0,1,1) + z**3*(-995327.9999999999_dp)*Hr5(1,1,1,0,0) + &
                        z*(-995327.9999999997_dp)*Hr5(1,1,1,0,0) + (497663.99999999994_dp)*Hr5(1,1,1,0,0) + &
                        z**4*(497663.99999999994_dp)*Hr5(1,1,1,0,0) + z**2*(1.4929919999999998e6_dp)*Hr5(1,1,1,0,0) + &
                        z*(-2.985984e6_dp)*Hr5(1,1,1,0,1) + z**3*(-2.985984e6_dp)*Hr5(1,1,1,0,1) + (1.492992e6_dp)*Hr5(1,1,1,0,1) + &
                        z**4*(1.492992e6_dp)*Hr5(1,1,1,0,1) + z**2*(4.478976e6_dp)*Hr5(1,1,1,0,1) + &
                        z*(-2.985984e6_dp)*Hr5(1,1,1,1,0) + z**3*(-2.985984e6_dp)*Hr5(1,1,1,1,0) + (1.492992e6_dp)*Hr5(1,1,1,1,0) + &
                        z**4*(1.492992e6_dp)*Hr5(1,1,1,1,0) + z**2*(4.478976e6_dp)*Hr5(1,1,1,1,0) + &
                        z*(-8.957952e6_dp)*Hr5(1,1,1,1,1) + z**3*(-5.971968e6_dp)*Hr5(1,1,1,1,1) + (2.985984e6_dp)*Hr5(1,1,1,1,1) + &
                        z**4*(2.985984e6_dp)*Hr5(1,1,1,1,1) + z**2*(8.957952e6_dp)*Hr5(1,1,1,1,1))/(z*((-1._dp) + z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-56957.74555921394_dp) + nf**3*(-85.87046443478027_dp) + nf**4*(1.6896112884169436_dp) &
                            + nf**2*(2010.1231020738924_dp) + (774088.9973590237_dp)

                        else
                            Ibar_select = ((-774088.9973590238_dp) + nf**4*((-1.6896112884169436_dp) + z*(1.6896112884169436_dp)) + &
                            nf**3*(z*(-85.87046443478027_dp) + (85.87046443478027_dp)) + nf**2*((-2010.1231020738924_dp) + &
                            z*(2010.1231020738924_dp)) + nf*(z*(-56957.74555921393_dp) + (56957.74555921393_dp)) + &
                            z*(774088.9973590238_dp))/((-1._dp) + z*(1._dp))

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-147230.80635272549_dp) + nf**3*(-196.04529754655067_dp) + &
                            nf**5*(-0.05267489711934156_dp) + nf**4*(4.345679012345679_dp) + nf**2*(8049.063456198373_dp) + &
                            (496145.90612771263_dp)

                        else
                            Ibar_select = ((-496145.9061277126_dp) + nf**5*(z*(-0.05267489711934158_dp) + (0.05267489711934158_dp)) &
                            + nf**4*((-4.34567901234568_dp) + z*(4.34567901234568_dp)) + nf**3*(z*(-196.04529754655061_dp) + &
                            (196.04529754655061_dp)) + nf**2*((-8049.063456198376_dp) + z*(8049.063456198376_dp)) + &
                            nf*(z*(-147230.80635272549_dp) + (147230.80635272549_dp)) + z*(496145.9061277126_dp))/((-1._dp) + &
                            z*(1._dp))

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (-3.315065007976247e6_dp) + nf**2*(-9672.280150030523_dp) + nf**4*(-2.054320987654321_dp) &
                            + nf**3*(135.5851851851852_dp) + nf*(345062.5193933433_dp)

                        else
                            Ibar_select = (-3.315065007976247e6_dp) + nf**2*(-9672.280150030523_dp) + nf**4*(-2.054320987654321_dp) &
                            + nf**3*(135.5851851851852_dp) + nf*(345062.5193933433_dp)

                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-76192.05843049446_dp) + nf**2*(-1584._dp) + nf**3*(32._dp) + (2.1428854107571747e6_dp)
                        else
                            Ibar_select = ((-2.1428854107571747e6_dp) + nf**3*((-32._dp) + z*(32._dp)) + nf**2*(z*(-1584._dp) + &
                            (1584._dp)) + nf*(z*(-76192.05843049445_dp) + (76192.05843049445_dp)) + &
                            z*(2.1428854107571747e6_dp))/((-1._dp) + z*(1._dp))

                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (-60832.23372197786_dp) + nf*(-42240._dp) + nf**2*(1280._dp)
                        else
                            Ibar_select = (-60832.23372197786_dp) + nf*(-42240._dp) + nf**2*(1280._dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (-171072._dp) + nf*(10368._dp)
                        else
                            Ibar_select = (-171072._dp) + nf*(10368._dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (24883.2_dp)
                        else
                            Ibar_select = (24883.2_dp)
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
