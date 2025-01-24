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
                        Ibar_select = z**3*(-482.71622869222256_dp) + z**4*(-60.49891979174238_dp) + (-6.575689489654321_dp) + &
                        z**5*(1.668274502516088_dp) + (2.6451727005890224_dp)/z + z*(7.993737974147562_dp) + &
                        z**2*(537.483652803029_dp) + (z*(-6.666666666666666_dp) + (1.3333333333333333_dp) + &
                        z**2*(213.40596468377385_dp) + z**3*(392.525592453522_dp))*Log(z) + (z**3*(-68.67155311692338_dp) + &
                        z*(-1.5555555555555554_dp) + (-0.6666666666666666_dp) + z**2*(27.80137308718255_dp))*Log(z)**2 + &
                        (z*(0.2962962962962963_dp) + (0.5925925925925926_dp) + z**2*(1.7401934484961763_dp) + &
                        z**3*(23.858585392707216_dp))*Log(z)**3

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
                        Ibar_select = z**2*(-40.60694027360694_dp) + z*(-13.054680678027468_dp) + (-7.703703703703703_dp)/z + &
                        z**5*(-0.627328431372549_dp) + z**4*(1.8567738905846443_dp) + (11.297892423978965_dp) + &
                        z**3*(48.83798581122764_dp) + (z**3*(-36.29506406523005_dp) + z**2*(-17.27835736192615_dp) + &
                        z*(-9.777777777777779_dp) + (-0.8888888888888888_dp))*Log(z) + (z**2*(0.16551702514579544_dp) + &
                        z*(1.7777777777777777_dp) + (3.5555555555555554_dp) + z**3*(9.853247118358167_dp))*Log(z)**2

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
                        Ibar_select = (z**3*(-3.4341099262359136e8_dp) + z**6*(-39560.364185309394_dp) + &
                        z**2*(-3295.3625557732435_dp) + (-472.185801152152_dp) + z*(504.2200309897908_dp) + &
                        Nf*(z**4*(-4943.92914118205_dp) + z*(-57.61279662683297_dp) + z**6*(-16.560320644171973_dp) + &
                        (5.88282387997262_dp) + z**2*(50.50110896509391_dp) + z**5*(687.1501547509737_dp) + &
                        z**3*(4274.56817078873_dp)) + z**5*(3.2795734938256973e6_dp) + z**4*(3.4017424282250994e8_dp))/z + &
                        (z**3*(-2.0680195860704923e8_dp) + z**2*(-1.398157629035882e8_dp) + (-410.79394755522935_dp) + &
                        z*(-318.5923450744127_dp) + (-78.98467140066984_dp)/z + Nf*((-35.457848078043234_dp) + &
                        z*(-0.515555083643898_dp) + z**2*(1433.424565505143_dp) + z**3*(2183.0567985987045_dp)))*Log(z) + &
                        (z**2*(-2.4322447525880173e7_dp) + z*(-133.9743881253827_dp) + (106.97016215684792_dp) + &
                        Nf*(z**3*(-1186.0167155711883_dp) + (-13.966979547350224_dp) + z*(-3.9053004732635905_dp) + &
                        z**2*(180.08534114311618_dp)) + z**3*(5.4647073401290886e7_dp))*Log(z)**2 + (z**3*(-9.882034596600426e6_dp) &
                        + z**2*(-2.264770837914974e6_dp) + (-36.21099191941935_dp) + z*(-17.07029459145344_dp) + &
                        Nf*((-2.5569272976680386_dp) + z*(-0.7352537722908093_dp) + z**2*(9.383835657232206_dp) + &
                        z**3*(69.21788306131498_dp)))*Log(z)**3 + (z**2*(-113156.93113263487_dp) + Nf*(z**3*(-37.21414424673223_dp) &
                        + (-0.3786008230452675_dp) + z*(-0.11522633744855965_dp) + z**2*(-0.056455745460574316_dp)) + &
                        z*(0.8641975308641978_dp) + (7.506172839506172_dp) + z**3*(903434.533209043_dp))*Log(z)**4 + &
                        (z**3*(-86524.9654140194_dp) + z**2*(-2437.510647811278_dp) + (-0.45925925925925926_dp) + &
                        z*(1.288888888888889_dp))*Log(z)**5 + ((-2795.883307118071_dp) + z**2*(-1070.5460224804667_dp) + &
                        z**3*(-878.5694657716556_dp) + Nf*(z**2*(-7.379410689810273_dp) + (-5.030944165540118_dp) + &
                        z**3*(3.083698282711209_dp) + z*(9.326656572639182_dp)) + z*(4744.998795370194_dp))*Log(z*(-1._dp) + &
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
                        Ibar_select = (z**3*(-4.137570544800112e6_dp) + z**5*(-143224.0584354666_dp) + z*(-455.71095448386416_dp) + &
                        z**2*(-216.886626563896_dp) + (483.9673467888446_dp) + Nf*(z**3*(-1464.9706389304206_dp) + &
                        z*(-53.202465875626004_dp) + z**6*(-5.3557525871830025_dp) + (-1.1565265637483262_dp) + &
                        z**2*(56.37673348167625_dp) + z**5*(187.3190234295046_dp) + z**4*(1280.9896270460613_dp)) + &
                        z**6*(3284.794493857353_dp) + z**4*(4.277698438734031e6_dp))/z + (z**3*(-2.3497334571323595e6_dp) + &
                        z**2*(-1.651476910119308e6_dp) + z*(-248.1425173541941_dp) + Nf*(z**3*(-1081.7001034429313_dp) + &
                        z**2*(-555.5684605091913_dp) + (-38.14068756693837_dp) + z*(-8.562572247989484_dp)) + &
                        (98.11244277743864_dp)/z + (497.14886135146355_dp))*Log(z) + (z**2*(-267989.4145365464_dp) + &
                        (-117.43267732966508_dp) + z*(-33.550538103564136_dp) + Nf*(z**2*(-80.05862177563473_dp) + &
                        (-11.45679012345679_dp) + z*(-5.135802469135802_dp) + z**3*(176.27122444255033_dp)) + &
                        z**3*(750337.2465155071_dp))*Log(z)**2 + (z**3*(-90958.87649645458_dp) + z**2*(-21011.775038243_dp) + &
                        z*(-13.234567901234568_dp) + Nf*(z**3*(-67.89203291934805_dp) + z**2*(-4.974774569535159_dp) + &
                        (-1.7777777777777777_dp) + z*(-0.5925925925925926_dp)) + (36.345679012345684_dp))*Log(z)**3 + &
                        (z**2*(-676.3207898614029_dp) + (-2.9382716049382718_dp) + z*(8.271604938271604_dp) + &
                        z**3*(14692.230963005815_dp))*Log(z)**4 + (z**2*(-1447.6040371123736_dp) + (-544.4791642015681_dp) + &
                        Nf*(z*(-4.559282020884582_dp) + z**3*(-0.17257441490570324_dp) + (0.7118723215079577_dp) + &
                        z**2*(4.019984114282327_dp)) + z**3*(437.7395614033088_dp) + z*(1554.3436399106329_dp))*Log(z*(-1._dp) + &
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
                        Ibar_select = (z**4*(-9731.273258243622_dp) + z**5*(-887.5293603025209_dp) + (-220.95532436198914_dp) + &
                        z**2*(-74.60853927232975_dp) + z**6*(11.140806991105725_dp) + Nf*(z**4*(-33.032735542313304_dp) + &
                        z*(-16.343545776367677_dp) + z**5*(-2.0156321469592604_dp) + z**6*(0.5980532083920125_dp) + &
                        (1.1851851851851851_dp) + z**2*(8.222169068877047_dp) + z**3*(41.38650653945485_dp)) + &
                        z*(251.62454557814533_dp) + z**3*(10651.601134270764_dp))/z + ((-137.5464974390064_dp) + &
                        z*(-116.44444444444444_dp) + (-46.22222222222222_dp)/z + Nf*((-8.88888888888889_dp) + &
                        z*(-6.518518518518518_dp) + z**2*(8.687208138840761_dp) + z**3*(24.043713577571747_dp)) + &
                        z**2*(3991.944090625632_dp) + z**3*(7496.693974079974_dp))*Log(z) + (z**3*(-1490.003176703887_dp) + &
                        z*(-25.185185185185187_dp) + Nf*(z**3*(-6.0637677494352396_dp) + (-2.3703703703703702_dp) + &
                        z*(-1.1851851851851851_dp) + z**2*(-0.15664547357115216_dp)) + (49.185185185185176_dp) + &
                        z**2*(576.4586405141953_dp))*Log(z)**2 + ((-5.135802469135802_dp) + z*(14.617283950617283_dp) + &
                        z**2*(32.0389454382639_dp) + z**3*(441.5034897636785_dp))*Log(z)**3 + (z*(-177.77963268425538_dp) + &
                        z**3*(-30.466535499847346_dp) + Nf*(z**2*(-3.511929408639212_dp) + (-2.9291277214703895_dp) + &
                        z**3*(1.1802897324733002_dp) + z*(5.260767397636301_dp)) + (87.4343263699818_dp) + &
                        z**2*(120.81184181412095_dp))*Log(z*(-1._dp) + (1._dp))

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
                        Ibar_select = (-50.62385710459352_dp) + z**4*(-4.407436322826194_dp) + z**3*(-2.6864439571778687_dp) + &
                        z**2*(-0.5265500760136046_dp) + Nf*((-0.7901234567901235_dp)/z + (-0.5925925925925926_dp) + &
                        z*(0.5925925925925926_dp) + z**2*(0.7901234567901235_dp)) + z**5*(1.019060773480663_dp) + &
                        z*(18.11411506631996_dp) + (39.11111111111111_dp)/z + (z**2*(-8.93517064115379_dp) + &
                        Nf*((-1.1851851851851851_dp) + z*(-1.1851851851851851_dp)) + z**3*(-0.866976222063479_dp) + &
                        (7.111111111111111_dp)/z + z*(9.481481481481486_dp) + (19.55555555555555_dp))*Log(z) + &
                        ((-2.9629629629629624_dp) + z**2*(-0.2623712088231026_dp) + z**3*(2.862027514358221_dp) + &
                        z*(13.037037037037035_dp))*Log(z)**2 + (z**2*(-19.900933497371888_dp) + (-16.598390495119958_dp) + &
                        z**3*(6.688308538163001_dp) + z*(29.81101545432885_dp))*Log(z*(-1._dp) + (1._dp))

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
                        Ibar_select = (z**3*(-140008.74876212107_dp) + z**4*(-120491.36783272112_dp) + z**2*(-68932.65964084199_dp) &
                        + z*(-32946.70923411154_dp) + nf*(-386.143015602376_dp) + nf*z*(-376.2017650436285_dp) + &
                        nf**2*z**3*(-17.696731935969694_dp) + nf**2*z**4*(-10.534979423868313_dp) + nf**2*(0.3664031069582032_dp) + &
                        nf**2*z**2*(56.95781199614854_dp) + nf**2*z*(57.187597837805896_dp) + nf*z**4*(1197.1006285558358_dp) + &
                        nf*z**2*(1217.9304463257856_dp) + nf*z**3*(2882.63236867176_dp) + (22375.471060762604_dp) + &
                        ((-2837.889696814507_dp) + z**4*(-484.7810524461406_dp) + nf**2*(z*(-3.8991029732698697_dp) + &
                        z**3*(-3.8991029732698697_dp)) + nf*((-23.39461783961922_dp) + z**4*(-23.39461783961922_dp) + &
                        z**3*(290.3823156154507_dp) + z*(290.3823156154508_dp) + z**2*(697.9394322153069_dp)) + &
                        z*(6930.288296909876_dp) + z**3*(11636.505585646612_dp) + z**2*(16690.19972495895_dp))*Hr1(-1) + &
                        nf*z**2*(-2392.7235752187225_dp)*Hr1(0) + nf*z*(-2348.9598555561392_dp)*Hr1(0) + &
                        nf*z**4*(-1159.7330744041965_dp)*Hr1(0) + nf*z**3*(-1096.2998583183735_dp)*Hr1(0) + &
                        nf*(-19.342312845824782_dp)*Hr1(0) + nf**2*z**4*(2.3703703703703702_dp)*Hr1(0) + &
                        nf**2*z**3*(32.21476164101245_dp)*Hr1(0) + nf**2*z*(47.622169048419856_dp)*Hr1(0) + &
                        nf**2*z**2*(73.56745734579206_dp)*Hr1(0) + (5458.699640174869_dp)*Hr1(0) + &
                        z**2*(51355.66190018232_dp)*Hr1(0) + z**3*(68859.21165189499_dp)*Hr1(0) + z*(74319.4053365387_dp)*Hr1(0) + &
                        z**4*(89022.34445378877_dp)*Hr1(0) + z**4*(-12247.672701091975_dp)*Hr1(1) + &
                        z**3*(-4643.8800003542_dp)*Hr1(1) + nf*z**3*(-684.4764792513485_dp)*Hr1(1) + &
                        nf*(-169.28655538377376_dp)*Hr1(1) + nf*z**2*(-64.98765432098755_dp)*Hr1(1) + &
                        nf**2*z*(-27.25925925925926_dp)*Hr1(1) + nf**2*z**4*(-2.3703703703703702_dp)*Hr1(1) + &
                        nf**2*(2.3703703703703702_dp)*Hr1(1) + nf**2*z**3*(27.25925925925926_dp)*Hr1(1) + &
                        nf*z**4*(234.2742097047614_dp)*Hr1(1) + nf*z*(684.4764792513487_dp)*Hr1(1) + (2588.7836896544595_dp)*Hr1(1) &
                        + z*(4643.880000354197_dp)*Hr1(1) + z**2*(9658.889011437512_dp)*Hr1(1) + &
                        z**2*(-753.8265748321749_dp)*Hr2(-1,-1) + nf*z*(-140.36770703771532_dp)*Hr2(-1,-1) + &
                        nf*z**3*(-140.36770703771532_dp)*Hr2(-1,-1) + (4356.597722133534_dp)*Hr2(-1,-1) + &
                        z**4*(4356.597722133534_dp)*Hr2(-1,-1) + z**3*(7061.960740459517_dp)*Hr2(-1,-1) + &
                        z*(7061.960740459519_dp)*Hr2(-1,-1) + z**3*(-44077.727696879505_dp)*Hr2(-1,0) + &
                        z**2*(-32810.541238505284_dp)*Hr2(-1,0) + z*(-32769.08572157087_dp)*Hr2(-1,0) + &
                        z**4*(-22232.42778950975_dp)*Hr2(-1,0) + (-16578.106801855432_dp)*Hr2(-1,0) + &
                        nf**2*z*(-10.271604938271604_dp)*Hr2(-1,0) + nf**2*z**3*(-10.271604938271604_dp)*Hr2(-1,0) + &
                        nf**2*z**2*(-4.7407407407407405_dp)*Hr2(-1,0) + nf*(79.53909465020575_dp)*Hr2(-1,0) + &
                        nf*z**4*(100.8724279835391_dp)*Hr2(-1,0) + nf*z*(820.3772984398004_dp)*Hr2(-1,0) + &
                        nf*z**2*(835.9506172839507_dp)*Hr2(-1,0) + nf*z**3*(863.043965106467_dp)*Hr2(-1,0) + &
                        z**4*(-7403.096845248393_dp)*Hr2(0,-1) + z**3*(-5678.314585179259_dp)*Hr2(0,-1) + &
                        (-2105.51560556573_dp)*Hr2(0,-1) + nf*z**3*(-246.2933378115468_dp)*Hr2(0,-1) + &
                        nf*z*(100.72682680947165_dp)*Hr2(0,-1) + z*(4064.9839243484353_dp)*Hr2(0,-1) + &
                        z**2*(5530.227717087766_dp)*Hr2(0,-1) + z*(-9977.170765967558_dp)*Hr2(0,0) + &
                        z**2*(-6654.2445136044_dp)*Hr2(0,0) + z**4*(-4288.541551294237_dp)*Hr2(0,0) + &
                        nf*z**2*(-1131.8415452756865_dp)*Hr2(0,0) + nf*z*(-662.6758114718295_dp)*Hr2(0,0) + &
                        nf**2*z**4*(-4.7407407407407405_dp)*Hr2(0,0) + nf**2*z**3*(18.17283950617284_dp)*Hr2(0,0) + &
                        nf**2*z*(22.91358024691358_dp)*Hr2(0,0) + nf**2*z**2*(37.925925925925924_dp)*Hr2(0,0) + &
                        nf*z**3*(100.82852337762745_dp)*Hr2(0,0) + z**3*(447.0099776172025_dp)*Hr2(0,0) + &
                        (588.6746566646319_dp)*Hr2(0,0) + nf*z**4*(785.2510288065843_dp)*Hr2(0,0) + (-2153.499958396644_dp)*Hr2(0,1) &
                        + nf*z**4*(-218.8641975308642_dp)*Hr2(0,1) + nf**2*z**2*(-35.55555555555556_dp)*Hr2(0,1) + &
                        nf**2*z*(-17.77777777777778_dp)*Hr2(0,1) + nf*(-16.59259259259259_dp)*Hr2(0,1) + &
                        nf**2*z**3*(-13.03703703703704_dp)*Hr2(0,1) + nf**2*z**4*(4.7407407407407405_dp)*Hr2(0,1) + &
                        nf*z*(96.7734267392672_dp)*Hr2(0,1) + nf*z**3*(348.62527859111907_dp)*Hr2(0,1) + &
                        nf*z**2*(634.0662597746045_dp)*Hr2(0,1) + z**2*(9677.716748526997_dp)*Hr2(0,1) + &
                        z*(10372.619878208488_dp)*Hr2(0,1) + z**3*(13312.73752595613_dp)*Hr2(0,1) + &
                        z**4*(17159.872264625257_dp)*Hr2(0,1) + (-13913.798714631414_dp)*Hr2(1,0) + &
                        z**2*(-4097.777777777779_dp)*Hr2(1,0) + nf*z**4*(-312.2304526748971_dp)*Hr2(1,0) + &
                        nf*z*(-309.59670781893004_dp)*Hr2(1,0) + z**3*(-306.8332977604946_dp)*Hr2(1,0) + &
                        nf*z**2*(-11.851851851851885_dp)*Hr2(1,0) + nf**2*(-4.7407407407407405_dp)*Hr2(1,0) + &
                        nf**2*z**3*(-4.7407407407407405_dp)*Hr2(1,0) + nf**2*z*(4.7407407407407405_dp)*Hr2(1,0) + &
                        nf**2*z**4*(4.7407407407407405_dp)*Hr2(1,0) + z*(306.8332977604946_dp)*Hr2(1,0) + &
                        nf*z**3*(309.59670781893004_dp)*Hr2(1,0) + nf*(324.082304526749_dp)*Hr2(1,0) + &
                        z**4*(18011.57649240919_dp)*Hr2(1,0) + (z*(-2792.87513526076_dp) + z**4*(-1820.4637074001841_dp) + &
                        nf**2*(z**3*(-8.296296296296296_dp) + z**4*(-4.7407407407407405_dp) + (4.7407407407407405_dp) + &
                        z*(8.296296296296296_dp)) + (43.574818511295554_dp) + nf*((-202.2716049382716_dp) + &
                        z**2*(-73.48148148148145_dp) + z**3*(-71.50617283950618_dp) + z*(71.50617283950618_dp) + &
                        z**4*(275.75308641975306_dp)) + z**2*(1776.8888888888887_dp) + z**3*(2792.8751352607596_dp))*Hr2(1,1) + &
                        z*(-1996.3407223141733_dp)*Hr3(-1,-1,-1) + z**3*(-1996.3407223141733_dp)*Hr3(-1,-1,-1) + &
                        nf*z*(-312.0987654320988_dp)*Hr3(-1,-1,0) + nf*z**3*(-312.0987654320988_dp)*Hr3(-1,-1,0) + &
                        nf*(-47.407407407407405_dp)*Hr3(-1,-1,0) + nf*z**4*(-47.407407407407405_dp)*Hr3(-1,-1,0) + &
                        nf*z**2*(-11.061728395061728_dp)*Hr3(-1,-1,0) + nf**2*z*(4.7407407407407405_dp)*Hr3(-1,-1,0) + &
                        nf**2*z**3*(4.7407407407407405_dp)*Hr3(-1,-1,0) + (2189.958847736626_dp)*Hr3(-1,-1,0) + &
                        z**4*(3386.9958847736625_dp)*Hr3(-1,-1,0) + z*(12854.744951553886_dp)*Hr3(-1,-1,0) + &
                        z**2*(14008.098765432098_dp)*Hr3(-1,-1,0) + z**3*(15248.819025627961_dp)*Hr3(-1,-1,0) + &
                        z*(-291.1330220041503_dp)*Hr3(-1,0,-1) + z**3*(-291.1330220041503_dp)*Hr3(-1,0,-1) + &
                        z**2*(-1515.8518518518517_dp)*Hr3(-1,0,0) + nf*(-94.81481481481481_dp)*Hr3(-1,0,0) + &
                        nf*z**4*(-94.81481481481481_dp)*Hr3(-1,0,0) + nf**2*z*(-4.7407407407407405_dp)*Hr3(-1,0,0) + &
                        nf**2*z**3*(-4.7407407407407405_dp)*Hr3(-1,0,0) + nf*z**2*(87.7037037037037_dp)*Hr3(-1,0,0) + &
                        nf*z**3*(164.74074074074073_dp)*Hr3(-1,0,0) + nf*z*(164.7407407407408_dp)*Hr3(-1,0,0) + &
                        z**3*(695.5829164362839_dp)*Hr3(-1,0,0) + z*(2724.6199534733205_dp)*Hr3(-1,0,0) + &
                        z**4*(4059.3909465020574_dp)*Hr3(-1,0,0) + (5073.909465020576_dp)*Hr3(-1,0,0) + &
                        z**2*(-4168.691358024692_dp)*Hr3(-1,0,1) + z**3*(-2744.166976460303_dp)*Hr3(-1,0,1) + &
                        z*(-1080.1669764603025_dp)*Hr3(-1,0,1) + nf*z**2*(-429.82716049382714_dp)*Hr3(-1,0,1) + &
                        nf*z*(-270.22222222222223_dp)*Hr3(-1,0,1) + nf*z**3*(-270.22222222222223_dp)*Hr3(-1,0,1) + &
                        nf*(-9.481481481481481_dp)*Hr3(-1,0,1) + nf*z**4*(-9.481481481481481_dp)*Hr3(-1,0,1) + &
                        nf**2*z*(4.7407407407407405_dp)*Hr3(-1,0,1) + nf**2*z**3*(4.7407407407407405_dp)*Hr3(-1,0,1) + &
                        z**4*(279.7037037037037_dp)*Hr3(-1,0,1) + (1111.7037037037037_dp)*Hr3(-1,0,1) + &
                        z*(-6638.8726624875_dp)*Hr3(0,-1,-1) + z**3*(6056.606618479198_dp)*Hr3(0,-1,-1) + &
                        z**3*(-13142.429868728268_dp)*Hr3(0,-1,0) + z**2*(-5018.074074074075_dp)*Hr3(0,-1,0) + &
                        (-2090.6666666666665_dp)*Hr3(0,-1,0) + nf*z**4*(-121.67901234567901_dp)*Hr3(0,-1,0) + &
                        nf**2*z*(-7.111111111111111_dp)*Hr3(0,-1,0) + nf**2*z**3*(-7.111111111111111_dp)*Hr3(0,-1,0) + &
                        nf*(-1.5802469135802468_dp)*Hr3(0,-1,0) + nf*z**2*(191.20987654320993_dp)*Hr3(0,-1,0) + &
                        nf*z**3*(218.8641975308642_dp)*Hr3(0,-1,0) + nf*z*(350.0246913580247_dp)*Hr3(0,-1,0) + &
                        z**4*(8073.744855967078_dp)*Hr3(0,-1,0) + z*(10003.164222422674_dp)*Hr3(0,-1,0) + &
                        z**3*(-4638.632837200055_dp)*Hr3(0,0,-1) + z**2*(5408.05582392531_dp)*Hr3(0,0,-1) + &
                        z*(11169.630317427089_dp)*Hr3(0,0,-1) + z**4*(-10975.604938271606_dp)*Hr3(0,0,0) + &
                        z**3*(-4385.404143951035_dp)*Hr3(0,0,0) + nf*z**2*(-1367.3086419753085_dp)*Hr3(0,0,0) + &
                        nf*z*(-974.8148148148148_dp)*Hr3(0,0,0) + nf*z**3*(-565.1358024691358_dp)*Hr3(0,0,0) + &
                        nf**2*z*(10.666666666666666_dp)*Hr3(0,0,0) + nf**2*z**3*(10.666666666666666_dp)*Hr3(0,0,0) + &
                        nf**2*z**2*(14.222222222222221_dp)*Hr3(0,0,0) + nf*z**4*(83.75308641975309_dp)*Hr3(0,0,0) + &
                        z*(8106.8807347987695_dp)*Hr3(0,0,0) + z**2*(12783.565473594437_dp)*Hr3(0,0,0) + &
                        z**3*(-8745.50566121914_dp)*Hr3(0,0,1) + z**4*(-8290.502057613168_dp)*Hr3(0,0,1) + &
                        z*(-1605.4984078211285_dp)*Hr3(0,0,1) + z**2*(-1225.1510889561564_dp)*Hr3(0,0,1) + &
                        (-1109.3333333333333_dp)*Hr3(0,0,1) + nf*z**2*(-101.53086419753086_dp)*Hr3(0,0,1) + &
                        nf*z*(-83.35802469135803_dp)*Hr3(0,0,1) + nf**2*z**2*(-14.222222222222221_dp)*Hr3(0,0,1) + &
                        nf**2*z*(-9.481481481481481_dp)*Hr3(0,0,1) + nf**2*z**3*(-9.481481481481481_dp)*Hr3(0,0,1) + &
                        nf*z**4*(61.629629629629626_dp)*Hr3(0,0,1) + nf*z**3*(89.28395061728396_dp)*Hr3(0,0,1) + &
                        z**2*(-13881.504268011273_dp)*Hr3(0,1,0) + z**3*(-9982.47188555651_dp)*Hr3(0,1,0) + &
                        z*(-8769.434848519471_dp)*Hr3(0,1,0) + (-3320.8888888888896_dp)*Hr3(0,1,0) + &
                        z**4*(-2144.7901234567903_dp)*Hr3(0,1,0) + nf*z**3*(-153.679012345679_dp)*Hr3(0,1,0) + &
                        nf*z**2*(-118.91358024691358_dp)*Hr3(0,1,0) + nf*z*(-41.08641975308642_dp)*Hr3(0,1,0) + &
                        nf*z**4*(-26.864197530864196_dp)*Hr3(0,1,0) + nf**2*z*(2.3703703703703702_dp)*Hr3(0,1,0) + &
                        nf**2*z**3*(2.3703703703703702_dp)*Hr3(0,1,0) + nf*(83.75308641975309_dp)*Hr3(0,1,0) + &
                        z**2*(-2649.8516761533565_dp)*Hr3(0,1,1) + z*(-2530.7153419461965_dp)*Hr3(0,1,1) + &
                        (-2042.0740740740741_dp)*Hr3(0,1,1) + nf*z**3*(-218.4691358024692_dp)*Hr3(0,1,1) + &
                        nf*z**4*(-124.8395061728395_dp)*Hr3(0,1,1) + nf*z**2*(-81.77777777777779_dp)*Hr3(0,1,1) + &
                        nf**2*z*(7.111111111111111_dp)*Hr3(0,1,1) + nf**2*z**3*(7.111111111111111_dp)*Hr3(0,1,1) + &
                        nf**2*z**2*(14.222222222222221_dp)*Hr3(0,1,1) + nf*(56.888888888888886_dp)*Hr3(0,1,1) + &
                        nf*z*(68.74074074074073_dp)*Hr3(0,1,1) + z**3*(372.7249872719101_dp)*Hr3(0,1,1) + &
                        z**4*(3033.0205761316874_dp)*Hr3(0,1,1) + z**2*(-801.1851851851851_dp)*Hr3(1,0,0) + &
                        z**4*(-449.3168724279836_dp)*Hr3(1,0,0) + z*(-262.3868312757201_dp)*Hr3(1,0,0) + &
                        nf*z*(-175.01234567901236_dp)*Hr3(1,0,0) + nf*(-44.24691358024691_dp)*Hr3(1,0,0) + &
                        nf*z**4*(44.24691358024691_dp)*Hr3(1,0,0) + nf*z**3*(175.01234567901236_dp)*Hr3(1,0,0) + &
                        z**3*(262.3868312757201_dp)*Hr3(1,0,0) + (1250.5020576131687_dp)*Hr3(1,0,0) + &
                        (-3921.514403292181_dp)*Hr3(1,0,1) + z*(-1314.6995884773667_dp)*Hr3(1,0,1) + &
                        z**2*(-829.6296296296294_dp)*Hr3(1,0,1) + nf*z**4*(-93.23456790123457_dp)*Hr3(1,0,1) + &
                        nf*z**3*(-76.24691358024691_dp)*Hr3(1,0,1) + nf*z*(76.24691358024691_dp)*Hr3(1,0,1) + &
                        nf*(93.23456790123457_dp)*Hr3(1,0,1) + z**3*(1314.6995884773662_dp)*Hr3(1,0,1) + &
                        z**4*(4751.14403292181_dp)*Hr3(1,0,1) + (-4703.736625514403_dp)*Hr3(1,1,0) + &
                        z*(-1387.786008230453_dp)*Hr3(1,1,0) + z**2*(-829.6296296296298_dp)*Hr3(1,1,0) + &
                        nf*z**4*(-109.03703703703704_dp)*Hr3(1,1,0) + nf*z**3*(-103.90123456790123_dp)*Hr3(1,1,0) + &
                        nf*z*(103.90123456790123_dp)*Hr3(1,1,0) + nf*(109.03703703703704_dp)*Hr3(1,1,0) + &
                        z**3*(1387.7860082304524_dp)*Hr3(1,1,0) + z**4*(5533.366255144033_dp)*Hr3(1,1,0) + &
                        z**3*(-1031.3086419753085_dp)*Hr3(1,1,1) + z**4*(-247.7037037037037_dp)*Hr3(1,1,1) + &
                        (-214.51851851851848_dp)*Hr3(1,1,1) + nf*z*(-107.85185185185185_dp)*Hr3(1,1,1) + &
                        nf*(-61.629629629629626_dp)*Hr3(1,1,1) + nf*z**4*(61.629629629629626_dp)*Hr3(1,1,1) + &
                        nf*z**3*(107.85185185185185_dp)*Hr3(1,1,1) + z**2*(462.22222222222234_dp)*Hr3(1,1,1) + &
                        z*(1031.3086419753088_dp)*Hr3(1,1,1) + z*(-1205.7283950617284_dp)*Hr4(-1,-1,-1,0) + &
                        z**3*(-1205.7283950617284_dp)*Hr4(-1,-1,-1,0) + z**2*(-786.9629629629629_dp)*Hr4(-1,-1,-1,0) + &
                        nf*z*(94.81481481481481_dp)*Hr4(-1,-1,-1,0) + nf*z**3*(94.81481481481481_dp)*Hr4(-1,-1,-1,0) + &
                        (524.6419753086421_dp)*Hr4(-1,-1,-1,0) + z**4*(524.6419753086421_dp)*Hr4(-1,-1,-1,0) + &
                        (-3136.7901234567903_dp)*Hr4(-1,-1,0,0) + z**4*(-3136.7901234567903_dp)*Hr4(-1,-1,0,0) + &
                        nf*z*(-107.45679012345678_dp)*Hr4(-1,-1,0,0) + nf*z**3*(-107.45679012345678_dp)*Hr4(-1,-1,0,0) + &
                        z*(337.7777777777776_dp)*Hr4(-1,-1,0,0) + z**3*(337.7777777777776_dp)*Hr4(-1,-1,0,0) + &
                        z**2*(4161.580246913581_dp)*Hr4(-1,-1,0,0) + z*(-3861.3333333333335_dp)*Hr4(-1,-1,0,1) + &
                        z**3*(-3861.3333333333335_dp)*Hr4(-1,-1,0,1) + (-2386.1728395061727_dp)*Hr4(-1,-1,0,1) + &
                        z**4*(-2386.1728395061727_dp)*Hr4(-1,-1,0,1) + z**2*(64.79012345679013_dp)*Hr4(-1,-1,0,1) + &
                        nf*z*(132.74074074074073_dp)*Hr4(-1,-1,0,1) + nf*z**3*(132.74074074074073_dp)*Hr4(-1,-1,0,1) + &
                        (-1238.9135802469136_dp)*Hr4(-1,0,-1,0) + z**4*(-1238.9135802469136_dp)*Hr4(-1,0,-1,0) + &
                        nf*z*(-126.41975308641975_dp)*Hr4(-1,0,-1,0) + nf*z**3*(-126.41975308641975_dp)*Hr4(-1,0,-1,0) + &
                        z*(1649.7777777777778_dp)*Hr4(-1,0,-1,0) + z**3*(1649.7777777777778_dp)*Hr4(-1,0,-1,0) + &
                        z**2*(2515.7530864197533_dp)*Hr4(-1,0,-1,0) + z**2*(-3628.246913580247_dp)*Hr4(-1,0,0,0) + &
                        z*(-619.8518518518522_dp)*Hr4(-1,0,0,0) + z**3*(-619.851851851852_dp)*Hr4(-1,0,0,0) + &
                        nf*z**3*(101.92592592592594_dp)*Hr4(-1,0,0,0) + nf*z*(101.92592592592595_dp)*Hr4(-1,0,0,0) + &
                        (2477.8271604938273_dp)*Hr4(-1,0,0,0) + z**4*(2477.8271604938273_dp)*Hr4(-1,0,0,0) + &
                        z**2*(-914.172839506173_dp)*Hr4(-1,0,0,1) + nf*z*(-28.444444444444443_dp)*Hr4(-1,0,0,1) + &
                        nf*z**3*(-28.444444444444443_dp)*Hr4(-1,0,0,1) + z**3*(2564.3456790123455_dp)*Hr4(-1,0,0,1) + &
                        z*(2564.345679012346_dp)*Hr4(-1,0,0,1) + (2741.7283950617284_dp)*Hr4(-1,0,0,1) + &
                        z**4*(2741.7283950617284_dp)*Hr4(-1,0,0,1) + nf*z*(22.123456790123456_dp)*Hr4(-1,0,1,0) + &
                        nf*z**3*(22.123456790123456_dp)*Hr4(-1,0,1,0) + z**2*(164.34567901234567_dp)*Hr4(-1,0,1,0) + &
                        z*(842.2716049382718_dp)*Hr4(-1,0,1,0) + z**3*(842.2716049382718_dp)*Hr4(-1,0,1,0) + &
                        (1087.20987654321_dp)*Hr4(-1,0,1,0) + z**4*(1087.20987654321_dp)*Hr4(-1,0,1,0) + &
                        nf*z*(-50.5679012345679_dp)*Hr4(-1,0,1,1) + nf*z**3*(-50.5679012345679_dp)*Hr4(-1,0,1,1) + &
                        (707.9506172839506_dp)*Hr4(-1,0,1,1) + z**4*(707.9506172839506_dp)*Hr4(-1,0,1,1) + &
                        z*(2958.2222222222226_dp)*Hr4(-1,0,1,1) + z**3*(2958.2222222222226_dp)*Hr4(-1,0,1,1) + &
                        z**2*(3438.617283950617_dp)*Hr4(-1,0,1,1) + z**4*(-2408.296296296296_dp)*Hr4(0,-1,-1,0) + &
                        (-512.0000000000001_dp)*Hr4(0,-1,-1,0) + nf*z**3*(-192.00000000000003_dp)*Hr4(0,-1,-1,0) + &
                        nf*z*(-60.839506172839506_dp)*Hr4(0,-1,-1,0) + z**3*(1585.382716049383_dp)*Hr4(0,-1,-1,0) + &
                        z*(1603.5555555555554_dp)*Hr4(0,-1,-1,0) + z**2*(2847.604938271605_dp)*Hr4(0,-1,-1,0) + &
                        z**2*(-4260.740740740742_dp)*Hr4(0,-1,0,0) + z*(-2826.074074074074_dp)*Hr4(0,-1,0,0) + &
                        nf*z**3*(47.01234567901235_dp)*Hr4(0,-1,0,0) + nf*z*(186.86419753086423_dp)*Hr4(0,-1,0,0) + &
                        z**3*(918.3209876543211_dp)*Hr4(0,-1,0,0) + (1792._dp)*Hr4(0,-1,0,0) + &
                        z**4*(3509.7283950617284_dp)*Hr4(0,-1,0,0) + z**2*(-1938.1728395061725_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(-91.65432098765432_dp)*Hr4(0,-1,0,1) + nf*z**3*(53.72839506172839_dp)*Hr4(0,-1,0,1) + &
                        z*(521.4814814814815_dp)*Hr4(0,-1,0,1) + (1024.0000000000002_dp)*Hr4(0,-1,0,1) + &
                        z**3*(2183.1111111111104_dp)*Hr4(0,-1,0,1) + z**4*(3296.3950617283945_dp)*Hr4(0,-1,0,1) + &
                        z**2*(-5527.703703703703_dp)*Hr4(0,0,-1,0) + z*(-4896.197530864199_dp)*Hr4(0,0,-1,0) + &
                        z**3*(-1409.7777777777778_dp)*Hr4(0,0,-1,0) + nf*z**3*(-41.48148148148148_dp)*Hr4(0,0,-1,0) + &
                        nf*z**2*(71.11111111111113_dp)*Hr4(0,0,-1,0) + nf*z*(240.5925925925926_dp)*Hr4(0,0,-1,0) + &
                        z**4*(821.7283950617284_dp)*Hr4(0,0,-1,0) + z*(-1482.962962962963_dp)*Hr4(0,0,0,0) + &
                        z**4*(-1112.4938271604938_dp)*Hr4(0,0,0,0) + nf*z**2*(-303.4074074074074_dp)*Hr4(0,0,0,0) + &
                        nf*z**3*(-275.55555555555554_dp)*Hr4(0,0,0,0) + nf*z*(-95.8024691358025_dp)*Hr4(0,0,0,0) + &
                        z**2*(382.81481481481455_dp)*Hr4(0,0,0,0) + z**3*(1646.9135802469134_dp)*Hr4(0,0,0,0) + &
                        z**3*(-4129.777777777778_dp)*Hr4(0,0,0,1) + z**4*(-2607.4074074074074_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(-265.48148148148147_dp)*Hr4(0,0,0,1) + nf*z*(-246.5185185185185_dp)*Hr4(0,0,0,1) + &
                        nf*z**3*(-63.20987654320987_dp)*Hr4(0,0,0,1) + z**2*(3499.061728395061_dp)*Hr4(0,0,0,1) + &
                        z*(5694.617283950617_dp)*Hr4(0,0,0,1) + z**3*(-6183.9012345679_dp)*Hr4(0,0,1,0) + &
                        z**2*(-3667.753086419753_dp)*Hr4(0,0,1,0) + z**4*(-1670.3209876543208_dp)*Hr4(0,0,1,0) + &
                        nf*z*(-27.65432098765432_dp)*Hr4(0,0,1,0) + nf*z**3*(98.76543209876546_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(109.03703703703705_dp)*Hr4(0,0,1,0) + (256.00000000000006_dp)*Hr4(0,0,1,0) + &
                        z*(1651.753086419753_dp)*Hr4(0,0,1,0) + z**2*(-6582.123456790124_dp)*Hr4(0,0,1,1) + &
                        z**3*(-5530.864197530865_dp)*Hr4(0,0,1,1) + z*(-2901.7283950617284_dp)*Hr4(0,0,1,1) + &
                        z**4*(-1703.506172839506_dp)*Hr4(0,0,1,1) + nf*z**3*(207.80246913580248_dp)*Hr4(0,0,1,1) + &
                        nf*z*(271.8024691358024_dp)*Hr4(0,0,1,1) + (383.99999999999994_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(429.03703703703707_dp)*Hr4(0,0,1,1) + z**3*(-1821.8271604938275_dp)*Hr4(0,1,0,0) + &
                        z**4*(-1286.320987654321_dp)*Hr4(0,1,0,0) + nf*z**2*(-241.77777777777783_dp)*Hr4(0,1,0,0) + &
                        nf*z*(-136.69135802469137_dp)*Hr4(0,1,0,0) + nf*z**3*(-136.69135802469137_dp)*Hr4(0,1,0,0) + &
                        z**2*(237.03703703703698_dp)*Hr4(0,1,0,0) + (625.7777777777779_dp)*Hr4(0,1,0,0) + &
                        z*(1843.9506172839506_dp)*Hr4(0,1,0,0) + z**2*(-4895.6049382716055_dp)*Hr4(0,1,0,1) + &
                        z**3*(-3594.666666666667_dp)*Hr4(0,1,0,1) + z*(-2037.3333333333337_dp)*Hr4(0,1,0,1) + &
                        z**4*(-982.9135802469136_dp)*Hr4(0,1,0,1) + (-199.11111111111111_dp)*Hr4(0,1,0,1) + &
                        nf*z*(76.64197530864199_dp)*Hr4(0,1,0,1) + nf*z**3*(76.64197530864199_dp)*Hr4(0,1,0,1) + &
                        nf*z**2*(184.8888888888889_dp)*Hr4(0,1,0,1) + z**2*(-5523.753086419753_dp)*Hr4(0,1,1,0) + &
                        z**3*(-3768.4938271604933_dp)*Hr4(0,1,1,0) + z*(-2570.666666666667_dp)*Hr4(0,1,1,0) + &
                        z**4*(-919.7037037037036_dp)*Hr4(0,1,1,0) + (-341.3333333333333_dp)*Hr4(0,1,1,0) + &
                        nf*z*(100.34567901234568_dp)*Hr4(0,1,1,0) + nf*z**3*(100.34567901234568_dp)*Hr4(0,1,1,0) + &
                        nf*z**2*(232.2962962962963_dp)*Hr4(0,1,1,0) + (-625.7777777777779_dp)*Hr4(0,1,1,1) + &
                        nf*z**2*(-184.8888888888889_dp)*Hr4(0,1,1,1) + nf*z*(-92.44444444444446_dp)*Hr4(0,1,1,1) + &
                        nf*z**3*(-92.44444444444446_dp)*Hr4(0,1,1,1) + z*(-24.888888888888903_dp)*Hr4(0,1,1,1) + &
                        z**4*(339.7530864197531_dp)*Hr4(0,1,1,1) + z**3*(1895.901234567901_dp)*Hr4(0,1,1,1) + &
                        z**2*(2157.037037037037_dp)*Hr4(0,1,1,1) + z*(-1349.5308641975307_dp)*Hr4(1,0,-1,0) + &
                        (-771.1604938271605_dp)*Hr4(1,0,-1,0) + z**4*(771.1604938271605_dp)*Hr4(1,0,-1,0) + &
                        z**3*(1349.5308641975307_dp)*Hr4(1,0,-1,0) + z**3*(-2585.283950617284_dp)*Hr4(1,0,0,0) + &
                        z**4*(-1365.3333333333333_dp)*Hr4(1,0,0,0) + (1365.3333333333333_dp)*Hr4(1,0,0,0) + &
                        z*(2585.283950617284_dp)*Hr4(1,0,0,0) + z**3*(-1629.2345679012344_dp)*Hr4(1,0,0,1) + &
                        z**4*(-739.5555555555557_dp)*Hr4(1,0,0,1) + (739.5555555555557_dp)*Hr4(1,0,0,1) + &
                        z*(1629.2345679012344_dp)*Hr4(1,0,0,1) + z**3*(-576.7901234567901_dp)*Hr4(1,0,1,0) + &
                        z**4*(-170.66666666666666_dp)*Hr4(1,0,1,0) + (170.66666666666666_dp)*Hr4(1,0,1,0) + &
                        z*(576.7901234567901_dp)*Hr4(1,0,1,0) + z*(-769.9753086419754_dp)*Hr4(1,0,1,1) + &
                        (-598.9135802469136_dp)*Hr4(1,0,1,1) + z**4*(598.9135802469136_dp)*Hr4(1,0,1,1) + &
                        z**3*(769.9753086419754_dp)*Hr4(1,0,1,1) + z**3*(-1037.432098765432_dp)*Hr4(1,1,0,0) + &
                        z**4*(-401.38271604938274_dp)*Hr4(1,1,0,0) + (401.38271604938274_dp)*Hr4(1,1,0,0) + &
                        z*(1037.432098765432_dp)*Hr4(1,1,0,0) + z*(-1074.172839506173_dp)*Hr4(1,1,0,1) + &
                        (-772.740740740741_dp)*Hr4(1,1,0,1) + z**4*(772.740740740741_dp)*Hr4(1,1,0,1) + &
                        z**3*(1074.172839506173_dp)*Hr4(1,1,0,1) + z*(-1433.6790123456788_dp)*Hr4(1,1,1,0) + &
                        (-978.1728395061727_dp)*Hr4(1,1,1,0) + z**4*(978.1728395061727_dp)*Hr4(1,1,1,0) + &
                        z**3*(1433.6790123456788_dp)*Hr4(1,1,1,0) + z*(-719.0123456790124_dp)*Hr4(1,1,1,1) + &
                        (-410.8641975308642_dp)*Hr4(1,1,1,1) + z**4*(410.8641975308642_dp)*Hr4(1,1,1,1) + &
                        z**3*(719.0123456790124_dp)*Hr4(1,1,1,1) + z*(606.8148148148148_dp)*Hr5(-1,-1,-1,-1,0) + &
                        z**3*(606.8148148148148_dp)*Hr5(-1,-1,-1,-1,0) + z*(-758.5185185185185_dp)*Hr5(-1,-1,-1,0,0) + &
                        z**3*(-758.5185185185185_dp)*Hr5(-1,-1,-1,0,0) + z*(1517.037037037037_dp)*Hr5(-1,-1,-1,0,1) + &
                        z**3*(1517.037037037037_dp)*Hr5(-1,-1,-1,0,1) + z*(-960.7901234567901_dp)*Hr5(-1,-1,0,-1,0) + &
                        z**3*(-960.7901234567901_dp)*Hr5(-1,-1,0,-1,0) + z*(910.2222222222222_dp)*Hr5(-1,-1,0,0,0) + &
                        z**3*(910.2222222222222_dp)*Hr5(-1,-1,0,0,0) + z*(-505.679012345679_dp)*Hr5(-1,-1,0,0,1) + &
                        z**3*(-505.679012345679_dp)*Hr5(-1,-1,0,0,1) + z*(-101.1358024691358_dp)*Hr5(-1,-1,0,1,0) + &
                        z**3*(-101.1358024691358_dp)*Hr5(-1,-1,0,1,0) + z*(-1112.4938271604938_dp)*Hr5(-1,-1,0,1,1) + &
                        z**3*(-1112.4938271604938_dp)*Hr5(-1,-1,0,1,1) + z*(-960.7901234567901_dp)*Hr5(-1,0,-1,-1,0) + &
                        z**3*(-960.7901234567901_dp)*Hr5(-1,0,-1,-1,0) + z*(935.5061728395061_dp)*Hr5(-1,0,-1,0,0) + &
                        z**3*(935.5061728395061_dp)*Hr5(-1,0,-1,0,0) + z*(-303.4074074074074_dp)*Hr5(-1,0,-1,0,1) + &
                        z**3*(-303.4074074074074_dp)*Hr5(-1,0,-1,0,1) + z*(606.8148148148148_dp)*Hr5(-1,0,0,-1,0) + &
                        z**3*(606.8148148148148_dp)*Hr5(-1,0,0,-1,0) + z*(-467.75308641975306_dp)*Hr5(-1,0,0,0,0) + &
                        z**3*(-467.75308641975306_dp)*Hr5(-1,0,0,0,0) + z*(-366.61728395061726_dp)*Hr5(-1,0,0,0,1) + &
                        z**3*(-366.61728395061726_dp)*Hr5(-1,0,0,0,1) + z*(-278.12345679012344_dp)*Hr5(-1,0,0,1,0) + &
                        z**3*(-278.12345679012344_dp)*Hr5(-1,0,0,1,0) + z*(151.7037037037037_dp)*Hr5(-1,0,0,1,1) + &
                        z**3*(151.7037037037037_dp)*Hr5(-1,0,0,1,1) + z*(-265.48148148148147_dp)*Hr5(-1,0,1,0,0) + &
                        z**3*(-265.48148148148147_dp)*Hr5(-1,0,1,0,0) + z*(-252.8395061728395_dp)*Hr5(-1,0,1,0,1) + &
                        z**3*(-252.8395061728395_dp)*Hr5(-1,0,1,0,1) + z*(-252.8395061728395_dp)*Hr5(-1,0,1,1,0) + &
                        z**3*(-252.8395061728395_dp)*Hr5(-1,0,1,1,0) + z*(-1848.8888888888887_dp)*Hr5(0,-1,-1,-1,0) + &
                        z**3*(-72.69135802469135_dp)*Hr5(0,-1,-1,-1,0) + z**3*(-3434.666666666666_dp)*Hr5(0,-1,-1,0,0) + &
                        z*(5305.6790123456785_dp)*Hr5(0,-1,-1,0,0) + z**3*(-3718.320987654321_dp)*Hr5(0,-1,-1,0,1) + &
                        z*(3111.5061728395062_dp)*Hr5(0,-1,-1,0,1) + z**3*(-935.5061728395061_dp)*Hr5(0,-1,0,-1,0) + &
                        z*(2553.6790123456794_dp)*Hr5(0,-1,0,-1,0) + z*(-3900.049382716049_dp)*Hr5(0,-1,0,0,0) + &
                        z**3*(2534.716049382716_dp)*Hr5(0,-1,0,0,0) + z*(-3398.320987654321_dp)*Hr5(0,-1,0,0,1) + &
                        z**3*(3196.0493827160494_dp)*Hr5(0,-1,0,0,1) + z*(-1061.9259259259259_dp)*Hr5(0,-1,0,1,0) + &
                        z**3*(910.2222222222222_dp)*Hr5(0,-1,0,1,0) + z*(303.4074074074074_dp)*Hr5(0,-1,0,1,1) + &
                        z**3*(303.4074074074074_dp)*Hr5(0,-1,0,1,1) + z**3*(-409.2839506172839_dp)*Hr5(0,0,-1,-1,0) + &
                        z**2*(2991.4074074074074_dp)*Hr5(0,0,-1,-1,0) + z*(4614.320987654321_dp)*Hr5(0,0,-1,-1,0) + &
                        z*(-3556.345679012346_dp)*Hr5(0,0,-1,0,0) + z**2*(2216.296296296296_dp)*Hr5(0,0,-1,0,0) + &
                        z**3*(4723.358024691358_dp)*Hr5(0,0,-1,0,0) + z*(-4483.16049382716_dp)*Hr5(0,0,-1,0,1) + &
                        z**2*(-1792.0000000000007_dp)*Hr5(0,0,-1,0,1) + z**3*(2615.3086419753085_dp)*Hr5(0,0,-1,0,1) + &
                        z**3*(1895.5061728395062_dp)*Hr5(0,0,0,-1,0) + z*(2535.5061728395067_dp)*Hr5(0,0,0,-1,0) + &
                        z**2*(4835.555555555557_dp)*Hr5(0,0,0,-1,0) + nf*z**2*(-312.88888888888897_dp)*Hr5(0,0,0,0,0) + &
                        nf*z*(-156.44444444444449_dp)*Hr5(0,0,0,0,0) + nf*z**3*(-156.44444444444449_dp)*Hr5(0,0,0,0,0) + &
                        z*(1547.8518518518522_dp)*Hr5(0,0,0,0,0) + z**3*(5100.246913580247_dp)*Hr5(0,0,0,0,0) + &
                        z**2*(6433.185185185184_dp)*Hr5(0,0,0,0,0) + z*(-1212.4444444444446_dp)*Hr5(0,0,0,0,1) + &
                        z**2*(3342.222222222222_dp)*Hr5(0,0,0,0,1) + z**3*(4933.925925925927_dp)*Hr5(0,0,0,0,1) + &
                        z*(429.82716049382714_dp)*Hr5(0,0,0,1,0) + z**2*(3607.703703703704_dp)*Hr5(0,0,0,1,0) + &
                        z**3*(3607.703703703704_dp)*Hr5(0,0,0,1,0) + z**3*(2433.58024691358_dp)*Hr5(0,0,0,1,1) + &
                        z*(3009.58024691358_dp)*Hr5(0,0,0,1,1) + z**2*(5089.185185185185_dp)*Hr5(0,0,0,1,1) + &
                        z*(515.5555555555554_dp)*Hr5(0,0,1,0,0) + z**3*(3491.160493827161_dp)*Hr5(0,0,1,0,0) + &
                        z**2*(3690.666666666667_dp)*Hr5(0,0,1,0,0) + z**3*(520.6913580246913_dp)*Hr5(0,0,1,0,1) + &
                        z*(1600.79012345679_dp)*Hr5(0,0,1,0,1) + z**2*(1818.074074074074_dp)*Hr5(0,0,1,0,1) + &
                        z**3*(46.617283950617285_dp)*Hr5(0,0,1,1,0) + z**2*(1509.925925925926_dp)*Hr5(0,0,1,1,0) + &
                        z*(1766.7160493827157_dp)*Hr5(0,0,1,1,0) + z**2*(-1818.074074074074_dp)*Hr5(0,0,1,1,1) + &
                        z**3*(-1069.037037037037_dp)*Hr5(0,0,1,1,1) + z*(-749.0370370370371_dp)*Hr5(0,0,1,1,1) + &
                        z**2*(-2313.4814814814818_dp)*Hr5(0,1,0,-1,0) + z*(-1156.7407407407409_dp)*Hr5(0,1,0,-1,0) + &
                        z**3*(-1156.7407407407409_dp)*Hr5(0,1,0,-1,0) + z*(2187.0617283950614_dp)*Hr5(0,1,0,0,0) + &
                        z**3*(2187.061728395062_dp)*Hr5(0,1,0,0,0) + z**2*(4323.555555555556_dp)*Hr5(0,1,0,0,0) + &
                        z*(1349.530864197531_dp)*Hr5(0,1,0,0,1) + z**3*(1349.530864197531_dp)*Hr5(0,1,0,0,1) + &
                        z**2*(2597.925925925926_dp)*Hr5(0,1,0,0,1) + z*(458.27160493827165_dp)*Hr5(0,1,0,1,0) + &
                        z**3*(458.2716049382717_dp)*Hr5(0,1,0,1,0) + z**2*(815.4074074074076_dp)*Hr5(0,1,0,1,0) + &
                        z**2*(-1493.3333333333333_dp)*Hr5(0,1,0,1,1) + z*(-696.0987654320987_dp)*Hr5(0,1,0,1,1) + &
                        z**3*(-696.0987654320987_dp)*Hr5(0,1,0,1,1) + z*(842.2716049382716_dp)*Hr5(0,1,1,0,0) + &
                        z**3*(842.2716049382716_dp)*Hr5(0,1,1,0,0) + z**2*(1583.4074074074074_dp)*Hr5(0,1,1,0,0) + &
                        z**2*(-2014.814814814815_dp)*Hr5(0,1,1,0,1) + z*(-956.8395061728396_dp)*Hr5(0,1,1,0,1) + &
                        z**3*(-956.8395061728396_dp)*Hr5(0,1,1,0,1) + z**2*(-2631.1111111111118_dp)*Hr5(0,1,1,1,0) + &
                        z*(-1264.987654320988_dp)*Hr5(0,1,1,1,0) + z**3*(-1264.987654320988_dp)*Hr5(0,1,1,1,0) + &
                        z**2*(-1232.5925925925926_dp)*Hr5(0,1,1,1,1) + z*(-616.2962962962963_dp)*Hr5(0,1,1,1,1) + &
                        z**3*(-616.2962962962963_dp)*Hr5(0,1,1,1,1))/(z*(z + (1._dp)))

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
                        Ibar_select = ((-5155.828741530614_dp) + nf*z*(-352.8381326682064_dp) + nf*z**4*(-319.11655903545085_dp) + &
                        nf**2*z**3*(-6.487017770906506_dp) + nf**2*(-0.7901234567901234_dp) + nf**2*z**4*(0.7901234567901234_dp) + &
                        nf**2*z**2*(7.798205946539738_dp) + nf**2*z*(11.685821735266329_dp) + nf*(133.03829986033438_dp) + &
                        nf*z**2*(158.6519681080079_dp) + nf*z**3*(467.57845104237845_dp) + z**2*(4189.821114673114_dp) + &
                        z**3*(6827.728970392585_dp) + z*(7751.589223472723_dp) + z**4*(16787.908742689935_dp) + &
                        ((-858.6691214467646_dp) + z**4*(-858.6691214467646_dp) + z*(nf*(-62.385647572317914_dp) + &
                        (1281.7847995975671_dp)) + z**3*(nf*(-62.385647572317914_dp) + (1281.7847995975674_dp)) + &
                        z**2*(1870.2697261784474_dp))*Hr1(-1) + z**3*(-7148.515643143709_dp)*Hr1(0) + &
                        z**2*(-6642.755006077893_dp)*Hr1(0) + z*(-5812.01531866995_dp)*Hr1(0) + z**4*(-3702.673411300039_dp)*Hr1(0) &
                        + (-1561.747223373883_dp)*Hr1(0) + nf*z**2*(-183.86937778076123_dp)*Hr1(0) + &
                        nf*z*(-57.53911314058605_dp)*Hr1(0) + nf**2*z**4*(-1.5802469135802468_dp)*Hr1(0) + &
                        nf**2*z**3*(4.345679012345679_dp)*Hr1(0) + nf**2*z*(5.925925925925926_dp)*Hr1(0) + &
                        nf**2*z**2*(11.851851851851851_dp)*Hr1(0) + nf*(25.371742112482853_dp)*Hr1(0) + &
                        nf*z**3*(164.79612367288345_dp)*Hr1(0) + nf*z**4*(203.36899862825783_dp)*Hr1(0) + &
                        (-2554.700222468202_dp)*Hr1(1) + z**2*(-1953.843621399176_dp)*Hr1(1) + z**3*(-1244.5174090345547_dp)*Hr1(1) &
                        + nf*z**4*(-92.3127572016461_dp)*Hr1(1) + nf*z*(-65.44855967078189_dp)*Hr1(1) + &
                        nf**2*z*(-2.765432098765432_dp)*Hr1(1) + nf**2*(-1.5802469135802468_dp)*Hr1(1) + &
                        nf**2*z**4*(1.5802469135802468_dp)*Hr1(1) + nf**2*z**3*(2.765432098765432_dp)*Hr1(1) + &
                        nf*z**2*(13.432098765432103_dp)*Hr1(1) + nf*z**3*(65.44855967078189_dp)*Hr1(1) + &
                        nf*(78.88065843621398_dp)*Hr1(1) + z*(1244.5174090345558_dp)*Hr1(1) + z**4*(4508.543843867379_dp)*Hr1(1) + &
                        z*(-582.2660440083005_dp)*Hr2(-1,-1) + z**3*(-582.2660440083005_dp)*Hr2(-1,-1) + &
                        nf*(-33.18518518518518_dp)*Hr2(-1,0) + nf*z**4*(-33.18518518518518_dp)*Hr2(-1,0) + &
                        nf**2*z*(-1.5802469135802468_dp)*Hr2(-1,0) + nf**2*z**3*(-1.5802469135802468_dp)*Hr2(-1,0) + &
                        nf*z*(48.592592592592595_dp)*Hr2(-1,0) + nf*z**3*(48.592592592592595_dp)*Hr2(-1,0) + &
                        nf*z**2*(87.70370370370371_dp)*Hr2(-1,0) + (1479.9012345679012_dp)*Hr2(-1,0) + &
                        z**4*(1479.9012345679012_dp)*Hr2(-1,0) + z**2*(1734.4526748971193_dp)*Hr2(-1,0) + &
                        z**3*(2497.441457450368_dp)*Hr2(-1,0) + z*(2497.4414574503685_dp)*Hr2(-1,0) + &
                        z**3*(-851.3041491639216_dp)*Hr2(0,-1) + z*(1433.570193172222_dp)*Hr2(0,-1) + &
                        z**4*(-2090.1399176954733_dp)*Hr2(0,0) + z**3*(-2063.388942872342_dp)*Hr2(0,0) + &
                        nf*z**2*(-397.4320987654321_dp)*Hr2(0,0) + nf*z*(-232.29629629629628_dp)*Hr2(0,0) + &
                        (-184.88888888888886_dp)*Hr2(0,0) + nf*z**3*(-177.77777777777777_dp)*Hr2(0,0) + &
                        nf**2*z*(3.1604938271604937_dp)*Hr2(0,0) + nf**2*z**3*(3.1604938271604937_dp)*Hr2(0,0) + &
                        nf**2*z**2*(4.7407407407407405_dp)*Hr2(0,0) + nf*z**4*(25.28395061728395_dp)*Hr2(0,0) + &
                        z*(2180.8011436595243_dp)*Hr2(0,0) + z**2*(2417.9835976608156_dp)*Hr2(0,0) + &
                        z**3*(-2722.486547694072_dp)*Hr2(0,1) + z**4*(-2136.230452674897_dp)*Hr2(0,1) + &
                        z**2*(-1884.3475809848528_dp)*Hr2(0,1) + z*(-1199.3260538669122_dp)*Hr2(0,1) + &
                        nf*z**2*(-47.80246913580247_dp)*Hr2(0,1) + nf*z*(-33.1851851851852_dp)*Hr2(0,1) + &
                        nf*(-4.7407407407407405_dp)*Hr2(0,1) + nf**2*z**2*(-4.7407407407407405_dp)*Hr2(0,1) + &
                        nf**2*z*(-2.3703703703703702_dp)*Hr2(0,1) + nf**2*z**3*(-2.3703703703703702_dp)*Hr2(0,1) + &
                        nf*z**3*(8.559670781893006_dp)*Hr2(0,1) + nf*z**4*(27.91769547325103_dp)*Hr2(0,1) + &
                        (98.76543209876543_dp)*Hr2(0,1) + z**4*(-700.312757201646_dp)*Hr2(1,0) + &
                        z**2*(-210.17283950617292_dp)*Hr2(1,0) + z*(-17.185185185185325_dp)*Hr2(1,0) + &
                        nf*z*(-15.670781893004113_dp)*Hr2(1,0) + nf*(-8.954732510288066_dp)*Hr2(1,0) + &
                        nf*z**4*(8.954732510288066_dp)*Hr2(1,0) + nf*z**3*(15.670781893004113_dp)*Hr2(1,0) + &
                        z**3*(17.185185185185176_dp)*Hr2(1,0) + (910.4855967078189_dp)*Hr2(1,0) + (z*(-683.1275720164609_dp) + &
                        (-638.6831275720165_dp) + z**2*(-420.34567901234567_dp) + nf*(z**3*(-67.29218106995884_dp) + &
                        z**4*(-38.452674897119344_dp) + (38.452674897119344_dp) + z*(67.29218106995884_dp)) + &
                        z**3*(683.1275720164609_dp) + z**4*(1059.0288065843622_dp))*Hr2(1,1) + (-399.275720164609_dp)*Hr3(-1,-1,0) + &
                        z**4*(-399.275720164609_dp)*Hr3(-1,-1,0) + nf*z*(-37.925925925925924_dp)*Hr3(-1,-1,0) + &
                        nf*z**3*(-37.925925925925924_dp)*Hr3(-1,-1,0) + z*(525.9588477366254_dp)*Hr3(-1,-1,0) + &
                        z**3*(525.9588477366254_dp)*Hr3(-1,-1,0) + z**2*(902.3209876543209_dp)*Hr3(-1,-1,0) + &
                        z**2*(-794.0740740740741_dp)*Hr3(-1,0,0) + z*(-391.76954732510296_dp)*Hr3(-1,0,0) + &
                        z**3*(-391.7695473251029_dp)*Hr3(-1,0,0) + nf*z*(28.44444444444444_dp)*Hr3(-1,0,0) + &
                        nf*z**3*(28.44444444444444_dp)*Hr3(-1,0,0) + (360.82304526748965_dp)*Hr3(-1,0,0) + &
                        z**4*(360.82304526748965_dp)*Hr3(-1,0,0) + z**2*(-685.8271604938273_dp)*Hr3(-1,0,1) + &
                        z*(-257.5802469135802_dp)*Hr3(-1,0,1) + z**3*(-257.5802469135802_dp)*Hr3(-1,0,1) + &
                        nf*z*(18.962962962962962_dp)*Hr3(-1,0,1) + nf*z**3*(18.962962962962962_dp)*Hr3(-1,0,1) + &
                        (322.3703703703703_dp)*Hr3(-1,0,1) + z**4*(322.3703703703703_dp)*Hr3(-1,0,1) + &
                        z*(-712.2962962962964_dp)*Hr3(0,-1,0) + z**2*(-451.16049382716045_dp)*Hr3(0,-1,0) + &
                        nf*z**3*(-21.333333333333332_dp)*Hr3(0,-1,0) + nf*z*(59.25925925925925_dp)*Hr3(0,-1,0) + &
                        z**4*(143.27572016460906_dp)*Hr3(0,-1,0) + z**3*(186.33744855967075_dp)*Hr3(0,-1,0) + &
                        (255.99999999999997_dp)*Hr3(0,-1,0) + z*(-381.6296296296296_dp)*Hr3(0,0,0) + &
                        z**4*(-168.559670781893_dp)*Hr3(0,0,0) + nf*z**3*(-71.11111111111111_dp)*Hr3(0,0,0) + &
                        nf*z**2*(-61.629629629629626_dp)*Hr3(0,0,0) + nf*z*(-9.481481481481481_dp)*Hr3(0,0,0) + &
                        z**2*(371.358024691358_dp)*Hr3(0,0,0) + z**3*(821.4650205761318_dp)*Hr3(0,0,0) + &
                        z**3*(-426.5349794238683_dp)*Hr3(0,0,1) + z**4*(-383.4732510288066_dp)*Hr3(0,0,1) + &
                        nf*z**2*(-97.97530864197529_dp)*Hr3(0,0,1) + nf*z*(-79.80246913580247_dp)*Hr3(0,0,1) + &
                        nf*z**3*(-37.1358024691358_dp)*Hr3(0,0,1) + (170.66666666666663_dp)*Hr3(0,0,1) + z**2*(1152._dp)*Hr3(0,0,1) &
                        + z*(1602.7654320987654_dp)*Hr3(0,0,1) + z**3*(-563.2263374485597_dp)*Hr3(0,1,0) + &
                        z**4*(-405.59670781893_dp)*Hr3(0,1,0) + nf*z**2*(-26.864197530864196_dp)*Hr3(0,1,0) + &
                        nf*z*(-13.432098765432098_dp)*Hr3(0,1,0) + nf*z**3*(-13.432098765432098_dp)*Hr3(0,1,0) + &
                        z**2*(266.2716049382716_dp)*Hr3(0,1,0) + (369.7777777777777_dp)*Hr3(0,1,0) + &
                        z*(793.679012345679_dp)*Hr3(0,1,0) + z**2*(-1729.5802469135804_dp)*Hr3(0,1,1) + &
                        z**3*(-1707.061728395062_dp)*Hr3(0,1,1) + z**4*(-524.641975308642_dp)*Hr3(0,1,1) + &
                        z*(-224.79012345679016_dp)*Hr3(0,1,1) + nf*z*(57.679012345679006_dp)*Hr3(0,1,1) + &
                        nf*z**3*(57.679012345679006_dp)*Hr3(0,1,1) + nf*z**2*(115.35802469135801_dp)*Hr3(0,1,1) + &
                        (322.3703703703703_dp)*Hr3(0,1,1) + z**3*(-548.4773662551439_dp)*Hr3(1,0,0) + &
                        z**4*(-313.4156378600823_dp)*Hr3(1,0,0) + (313.4156378600823_dp)*Hr3(1,0,0) + &
                        z*(548.4773662551439_dp)*Hr3(1,0,0) + z**3*(-615.7695473251028_dp)*Hr3(1,0,1) + &
                        z**4*(-351.86831275720164_dp)*Hr3(1,0,1) + (351.86831275720164_dp)*Hr3(1,0,1) + &
                        z*(615.7695473251028_dp)*Hr3(1,0,1) + z**3*(-615.7695473251028_dp)*Hr3(1,1,0) + &
                        z**4*(-351.86831275720164_dp)*Hr3(1,1,0) + (351.86831275720164_dp)*Hr3(1,1,0) + &
                        z*(615.7695473251028_dp)*Hr3(1,1,0) + z*(-303.4074074074074_dp)*Hr4(-1,-1,-1,0) + &
                        z**3*(-303.4074074074074_dp)*Hr4(-1,-1,-1,0) + z*(252.83950617283952_dp)*Hr4(-1,-1,0,0) + &
                        z**3*(252.83950617283952_dp)*Hr4(-1,-1,0,0) + z*(202.2716049382716_dp)*Hr4(-1,-1,0,1) + &
                        z**3*(202.2716049382716_dp)*Hr4(-1,-1,0,1) + z*(151.7037037037037_dp)*Hr4(-1,0,-1,0) + &
                        z**3*(151.7037037037037_dp)*Hr4(-1,0,-1,0) + z*(-101.1358024691358_dp)*Hr4(-1,0,0,0) + &
                        z**3*(-101.1358024691358_dp)*Hr4(-1,0,0,0) + z*(-151.7037037037037_dp)*Hr4(-1,0,0,1) + &
                        z**3*(-151.7037037037037_dp)*Hr4(-1,0,0,1) + z*(-50.5679012345679_dp)*Hr4(-1,0,1,0) + &
                        z**3*(-50.5679012345679_dp)*Hr4(-1,0,1,0) + z*(-101.1358024691358_dp)*Hr4(-1,0,1,1) + &
                        z**3*(-101.1358024691358_dp)*Hr4(-1,0,1,1) + z**3*(-371.358024691358_dp)*Hr4(0,-1,-1,0) + &
                        z*(674.7654320987652_dp)*Hr4(0,-1,-1,0) + z*(-604.4444444444443_dp)*Hr4(0,-1,0,0) + &
                        z**3*(351.6049382716049_dp)*Hr4(0,-1,0,0) + z*(-534.1234567901234_dp)*Hr4(0,-1,0,1) + &
                        z**3*(331.85185185185185_dp)*Hr4(0,-1,0,1) + z*(263.9012345679012_dp)*Hr4(0,0,-1,0) + &
                        z**3*(786.9629629629629_dp)*Hr4(0,0,-1,0) + z**2*(1151.9999999999995_dp)*Hr4(0,0,-1,0) + &
                        nf*z**2*(-56.88888888888888_dp)*Hr4(0,0,0,0) + nf*z*(-28.44444444444444_dp)*Hr4(0,0,0,0) + &
                        nf*z**3*(-28.44444444444444_dp)*Hr4(0,0,0,0) + z*(276.5432098765432_dp)*Hr4(0,0,0,0) + &
                        z**3*(1103.0123456790122_dp)*Hr4(0,0,0,0) + z**2*(1328.9876543209875_dp)*Hr4(0,0,0,0) + &
                        z*(-480.39506172839504_dp)*Hr4(0,0,0,1) + z**2*(894.4197530864199_dp)*Hr4(0,0,0,1) + &
                        z**3*(1475.9506172839504_dp)*Hr4(0,0,0,1) + z*(244.93827160493828_dp)*Hr4(0,0,1,0) + &
                        z**3*(1183.6049382716049_dp)*Hr4(0,0,1,0) + z**2*(1377.9753086419753_dp)*Hr4(0,0,1,0) + &
                        z*(900.7407407407408_dp)*Hr4(0,0,1,1) + z**3*(900.7407407407408_dp)*Hr4(0,0,1,1) + &
                        z**2*(1700.3456790123457_dp)*Hr4(0,0,1,1) + z*(470.12345679012344_dp)*Hr4(0,1,0,0) + &
                        z**3*(470.12345679012344_dp)*Hr4(0,1,0,0) + z**2*(940.2469135802469_dp)*Hr4(0,1,0,0) + &
                        z*(527.8024691358024_dp)*Hr4(0,1,0,1) + z**3*(527.8024691358024_dp)*Hr4(0,1,0,1) + &
                        z**2*(1055.6049382716049_dp)*Hr4(0,1,0,1) + z*(527.8024691358024_dp)*Hr4(0,1,1,0) + &
                        z**3*(527.8024691358024_dp)*Hr4(0,1,1,0) + z**2*(1055.6049382716049_dp)*Hr4(0,1,1,0))/(z*(z + (1._dp)))

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
                        Ibar_select = ((-107237.23285764939_dp) + nf*z**4*(-19010.08237291654_dp) + nf*z*(-9120.939057232137_dp) + &
                        nf**2*z**3*(-576.3472052566617_dp) + nf**2*(-90.31270205667104_dp) + nf**3*z*(-7.790547823510888_dp) + &
                        nf**3*z**2*(-5.198803964359825_dp) + nf**3*z**4*(-0.5267489711934157_dp) + nf**3*(0.5267489711934157_dp) + &
                        nf**3*z**3*(4.3246785139376716_dp) + nf**2*z**2*(16.90377296554736_dp) + nf**2*z**4*(218.97174416874213_dp) &
                        + nf**2*z*(590.1853264539891_dp) + nf*z**3*(863.5568654147571_dp) + nf*z**2*(2417.6998964923_dp) + &
                        nf*(5930.765024741579_dp) + z**2*(187759.47053503798_dp) + z*(197263.05344895637_dp) + &
                        z**3*(301314.66213063063_dp) + z**4*(428328.00975198246_dp) + ((-43386.493766560634_dp) + &
                        z**4*(-43386.493766560634_dp) + z**3*(-40654.873813341226_dp) + z*(-40654.873813341204_dp) + &
                        z**2*(-25800.047406724985_dp) + nf**2*(z*(50.83274987374052_dp) + z**3*(50.832749873740525_dp)) + &
                        nf*(z**2*(-2156.203944218238_dp) + z**3*(-2070.2910693758154_dp) + z*(-2070.291069375815_dp) + &
                        (919.8994792492247_dp) + z**4*(919.8994792492247_dp)))*Hr1(-1) + z**2*(-171573.67841174387_dp)*Hr1(0) + &
                        z*(-145190.82554733325_dp)*Hr1(0) + z**3*(-107756.0243689796_dp)*Hr1(0) + (-37135.867149334656_dp)*Hr1(0) + &
                        z**4*(-31316.683156337316_dp)*Hr1(0) + nf**2*z**4*(-165.17969821673526_dp)*Hr1(0) + &
                        nf**2*(-16.095107453132144_dp)*Hr1(0) + nf**3*z**2*(-7.901234567901233_dp)*Hr1(0) + &
                        nf**3*z*(-3.9506172839506166_dp)*Hr1(0) + nf**3*z**3*(-2.8971193415637853_dp)*Hr1(0) + &
                        nf**3*z**4*(1.0534979423868314_dp)*Hr1(0) + nf**2*z**3*(28.011235690797406_dp)*Hr1(0) + &
                        nf**2*z*(246.73450820906498_dp)*Hr1(0) + nf**2*z**2*(506.85329944347035_dp)*Hr1(0) + &
                        nf*(1574.4864368705594_dp)*Hr1(0) + nf*z*(6339.649829899339_dp)*Hr1(0) + &
                        nf*z**4*(6780.612470334416_dp)*Hr1(0) + nf*z**2*(7199.837304874719_dp)*Hr1(0) + &
                        nf*z**3*(10924.756642838594_dp)*Hr1(0) + (-84150.20390796897_dp)*Hr1(1) + &
                        z**2*(-54868.590314462905_dp)*Hr1(1) + z**3*(-34724.92929054642_dp)*Hr1(1) + &
                        nf*z**4*(-4521.370452228104_dp)*Hr1(1) + nf*z*(-3975.0996509876495_dp)*Hr1(1) + &
                        nf**2*(-79.2318244170096_dp)*Hr1(1) + nf**2*z**2*(-5.530864197530864_dp)*Hr1(1) + &
                        nf**2*z*(-3.5994513031550133_dp)*Hr1(1) + nf**3*z**3*(-1.8436213991769548_dp)*Hr1(1) + &
                        nf**3*z**4*(-1.0534979423868314_dp)*Hr1(1) + nf**3*(1.0534979423868314_dp)*Hr1(1) + &
                        nf**3*z*(1.8436213991769548_dp)*Hr1(1) + nf**2*z**3*(3.5994513031550133_dp)*Hr1(1) + &
                        nf**2*z**4*(84.76268861454047_dp)*Hr1(1) + nf*z**2*(1573.7064471879287_dp)*Hr1(1) + &
                        nf*(2947.664005040175_dp)*Hr1(1) + nf*z**3*(3975.0996509876495_dp)*Hr1(1) + z*(34724.92929054642_dp)*Hr1(1) &
                        + z**4*(139018.79422243187_dp)*Hr1(1) + z**2*(-17199.376448757095_dp)*Hr2(-1,-1) + &
                        z**3*(-13975.755567934784_dp)*Hr2(-1,-1) + z*(-13975.75556793478_dp)*Hr2(-1,-1) + &
                        nf*z*(776.3547253444007_dp)*Hr2(-1,-1) + nf*z**3*(776.3547253444008_dp)*Hr2(-1,-1) + &
                        (8508.709155002249_dp)*Hr2(-1,-1) + z**4*(8508.709155002249_dp)*Hr2(-1,-1) + &
                        nf*(-1898.9300411522634_dp)*Hr2(-1,0) + nf*z**4*(-1898.930041152263_dp)*Hr2(-1,0) + &
                        nf*z*(-1603.1949199475418_dp)*Hr2(-1,0) + nf*z**3*(-1603.1949199475418_dp)*Hr2(-1,0) + &
                        nf**2*z**2*(-66.45816186556928_dp)*Hr2(-1,0) + nf**2*z*(-50.085048010973935_dp)*Hr2(-1,0) + &
                        nf**2*z**3*(-50.085048010973935_dp)*Hr2(-1,0) + nf**3*z*(1.0534979423868314_dp)*Hr2(-1,0) + &
                        nf**3*z**3*(1.0534979423868314_dp)*Hr2(-1,0) + nf**2*(23.703703703703702_dp)*Hr2(-1,0) + &
                        nf**2*z**4*(23.703703703703702_dp)*Hr2(-1,0) + nf*z**2*(557.7393689986283_dp)*Hr2(-1,0) + &
                        (32947.38054981388_dp)*Hr2(-1,0) + z**4*(32947.38054981388_dp)*Hr2(-1,0) + &
                        z**2*(69011.2506829668_dp)*Hr2(-1,0) + z*(73241.03239092037_dp)*Hr2(-1,0) + &
                        z**3*(73241.03239092037_dp)*Hr2(-1,0) + z**3*(-16944.503177854596_dp)*Hr2(0,-1) + &
                        (-5874.648479726603_dp)*Hr2(0,-1) + z**4*(-2634.060675275645_dp)*Hr2(0,-1) + &
                        nf*z*(-1573.9379002099372_dp)*Hr2(0,-1) + nf*z**3*(797.5831748655365_dp)*Hr2(0,-1) + &
                        z**2*(8599.688224378548_dp)*Hr2(0,-1) + z*(30920.258745789375_dp)*Hr2(0,-1) + &
                        z**3*(-72871.2915123544_dp)*Hr2(0,0) + z**4*(-61870.87779773989_dp)*Hr2(0,0) + &
                        nf*z**2*(-7317.712461562962_dp)*Hr2(0,0) + (-7227.463892343869_dp)*Hr2(0,0) + &
                        nf*z*(-4331.5362379658645_dp)*Hr2(0,0) + nf*z**3*(-732.1944251135506_dp)*Hr2(0,0) + &
                        nf**2*z**4*(-16.153635116598075_dp)*Hr2(0,0) + nf**3*z**2*(-3.160493827160493_dp)*Hr2(0,0) + &
                        nf**3*z*(-2.1069958847736623_dp)*Hr2(0,0) + nf**3*z**3*(-2.1069958847736623_dp)*Hr2(0,0) + &
                        nf*(168.8230452674897_dp)*Hr2(0,0) + nf**2*z**3*(172.5046038134933_dp)*Hr2(0,0) + &
                        nf**2*z*(274.2988424966208_dp)*Hr2(0,0) + nf**2*z**2*(422.39741064481933_dp)*Hr2(0,0) + &
                        nf*z**4*(2268.3859167809787_dp)*Hr2(0,0) + z**2*(7655.96265504701_dp)*Hr2(0,0) + &
                        z*(8508.132125594882_dp)*Hr2(0,0) + z**3*(-87822.34084218701_dp)*Hr2(0,1) + &
                        z**2*(-72139.16751907964_dp)*Hr2(0,1) + z**4*(-48264.867526517206_dp)*Hr2(0,1) + &
                        z*(-43429.184550526035_dp)*Hr2(0,1) + (-10847.490347116198_dp)*Hr2(0,1) + &
                        nf*z*(-1374.3064615755395_dp)*Hr2(0,1) + nf*z**2*(-1287.1479025749472_dp)*Hr2(0,1) + &
                        nf*(-29.936899862825793_dp)*Hr2(0,1) + nf**2*z**3*(-29.761316872427987_dp)*Hr2(0,1) + &
                        nf**2*z**2*(-22.56241426611797_dp)*Hr2(0,1) + nf**2*z**4*(-19.13854595336077_dp)*Hr2(0,1) + &
                        nf**2*z*(-13.168724279835388_dp)*Hr2(0,1) + nf**2*(-1.2290809327846366_dp)*Hr2(0,1) + &
                        nf**3*z*(1.5802469135802466_dp)*Hr2(0,1) + nf**3*z**3*(1.5802469135802466_dp)*Hr2(0,1) + &
                        nf**3*z**2*(3.160493827160493_dp)*Hr2(0,1) + nf*z**4*(2378.534979423868_dp)*Hr2(0,1) + &
                        nf*z**3*(2435.756638561634_dp)*Hr2(0,1) + z**2*(-10620.576131687243_dp)*Hr2(1,0) + &
                        z**4*(-7778.877662037416_dp)*Hr2(1,0) + z**3*(-5658.660051226658_dp)*Hr2(1,0) + &
                        nf*z*(-1300.5285779606766_dp)*Hr2(1,0) + nf*(-883.3872885230909_dp)*Hr2(1,0) + &
                        nf**2*z**3*(-6.4526748971193415_dp)*Hr2(1,0) + nf**2*z**4*(-3.6872427983539096_dp)*Hr2(1,0) + &
                        nf**2*(3.6872427983539087_dp)*Hr2(1,0) + nf**2*z*(6.45267489711934_dp)*Hr2(1,0) + &
                        nf*z**2*(122.99588477366255_dp)*Hr2(1,0) + nf*z**4*(760.3914037494285_dp)*Hr2(1,0) + &
                        nf*z**3*(1300.5285779606766_dp)*Hr2(1,0) + z*(5658.660051226656_dp)*Hr2(1,0) + &
                        (18399.453793724657_dp)*Hr2(1,0) + ((-27184.746729252103_dp) + z**2*(-21241.15226337448_dp) + &
                        z**3*(-5066.649328061223_dp) + nf**2*(z*(-64.21947873799724_dp) + (-36.69684499314129_dp) + &
                        z**4*(36.69684499314129_dp) + z**3*(64.21947873799724_dp)) + nf*(z**4*(-1823.7805212620028_dp) + &
                        z**3*(-334.26611796982166_dp) + z**2*(245.9917695473251_dp) + z*(334.26611796982166_dp) + &
                        (1577.7887517146776_dp)) + z*(5066.649328061218_dp) + z**4*(48425.89899262658_dp))*Hr2(1,1) + &
                        z*(4436.3127162537185_dp)*Hr3(-1,-1,-1) + z**3*(4436.3127162537185_dp)*Hr3(-1,-1,-1) + &
                        z**3*(-23327.706967221227_dp)*Hr3(-1,-1,0) + z*(-23327.706967221224_dp)*Hr3(-1,-1,0) + &
                        z**2*(-17487.802469135804_dp)*Hr3(-1,-1,0) + (-16857.54732510288_dp)*Hr3(-1,-1,0) + &
                        z**4*(-16857.54732510288_dp)*Hr3(-1,-1,0) + nf*z**2*(-1039.2757201646089_dp)*Hr3(-1,-1,0) + &
                        nf*z*(-913.997256515775_dp)*Hr3(-1,-1,0) + nf*z**3*(-913.997256515775_dp)*Hr3(-1,-1,0) + &
                        nf**2*z*(30.902606310013716_dp)*Hr3(-1,-1,0) + nf**2*z**3*(30.902606310013716_dp)*Hr3(-1,-1,0) + &
                        nf*(423.1550068587105_dp)*Hr3(-1,-1,0) + nf*z**4*(423.1550068587105_dp)*Hr3(-1,-1,0) + &
                        z*(-2218.1563581268592_dp)*Hr3(-1,0,-1) + z**3*(-2218.1563581268592_dp)*Hr3(-1,0,-1) + &
                        nf*(-385.40466392318245_dp)*Hr3(-1,0,0) + nf*z**4*(-385.40466392318245_dp)*Hr3(-1,0,0) + &
                        nf**2*z*(-23.17695473251029_dp)*Hr3(-1,0,0) + nf**2*z**3*(-23.17695473251029_dp)*Hr3(-1,0,0) + &
                        nf*z**3*(685.344307270233_dp)*Hr3(-1,0,0) + nf*z*(685.3443072702331_dp)*Hr3(-1,0,0) + &
                        nf*z**2*(915.2263374485597_dp)*Hr3(-1,0,0) + (15342.353909465019_dp)*Hr3(-1,0,0) + &
                        z**4*(15342.353909465019_dp)*Hr3(-1,0,0) + z**2*(16338.831275720166_dp)*Hr3(-1,0,0) + &
                        z*(20660.415901289052_dp)*Hr3(-1,0,0) + z**3*(20660.415901289052_dp)*Hr3(-1,0,0) + &
                        nf*(-347.6543209876543_dp)*Hr3(-1,0,1) + nf*z**4*(-347.6543209876543_dp)*Hr3(-1,0,1) + &
                        nf**2*z*(-15.451303155006858_dp)*Hr3(-1,0,1) + nf**2*z**3*(-15.451303155006858_dp)*Hr3(-1,0,1) + &
                        nf*z**3*(456.6913580246913_dp)*Hr3(-1,0,1) + nf*z*(456.6913580246914_dp)*Hr3(-1,0,1) + &
                        nf*z**2*(791.1769547325102_dp)*Hr3(-1,0,1) + (13827.160493827161_dp)*Hr3(-1,0,1) + &
                        z**4*(13827.160493827161_dp)*Hr3(-1,0,1) + z**2*(15189.860082304524_dp)*Hr3(-1,0,1) + &
                        z*(19176.14155969122_dp)*Hr3(-1,0,1) + z**3*(19176.14155969122_dp)*Hr3(-1,0,1) + &
                        z*(-13872.141911566803_dp)*Hr3(0,-1,-1) + z**3*(9435.829195313085_dp)*Hr3(0,-1,-1) + &
                        nf*z**3*(-793.986282578875_dp)*Hr3(0,-1,0) + nf*(-306.5679012345679_dp)*Hr3(0,-1,0) + &
                        nf*z**4*(-116.58710562414265_dp)*Hr3(0,-1,0) + nf**2*z*(-43.28120713305898_dp)*Hr3(0,-1,0) + &
                        nf**2*z**3*(12.378600823045266_dp)*Hr3(0,-1,0) + nf*z**2*(519.6378600823048_dp)*Hr3(0,-1,0) + &
                        nf*z*(1707.9835390946505_dp)*Hr3(0,-1,0) + z*(3646.927247674335_dp)*Hr3(0,-1,0) + &
                        z**4*(7020.510288065843_dp)*Hr3(0,-1,0) + z**2*(8743.901234567897_dp)*Hr3(0,-1,0) + &
                        (9837.037037037036_dp)*Hr3(0,-1,0) + z**3*(19680.77971954689_dp)*Hr3(0,-1,0) + &
                        z**2*(-26435.91815876972_dp)*Hr3(0,0,-1) + z**3*(-18305.56640339588_dp)*Hr3(0,0,-1) + &
                        z*(-6651.580849955933_dp)*Hr3(0,0,-1) + z*(-18797.597045726143_dp)*Hr3(0,0,0) + &
                        z**4*(-9823.458619112938_dp)*Hr3(0,0,0) + z**2*(-7455.518731875785_dp)*Hr3(0,0,0) + &
                        nf*z**3*(-2696.880991342825_dp)*Hr3(0,0,0) + nf*z**2*(-1436.4870701470225_dp)*Hr3(0,0,0) + &
                        (-554.6666666666666_dp)*Hr3(0,0,0) + nf**2*z*(34.41426611796981_dp)*Hr3(0,0,0) + &
                        nf**2*z**3*(65.84362139917694_dp)*Hr3(0,0,0) + nf**2*z**2*(84.8065843621399_dp)*Hr3(0,0,0) + &
                        nf*z**4*(137.65706447187927_dp)*Hr3(0,0,0) + nf*z*(989.2937840215908_dp)*Hr3(0,0,0) + &
                        z**3*(3062.0249128171695_dp)*Hr3(0,0,0) + z**4*(-20674.19478737997_dp)*Hr3(0,0,1) + &
                        z**3*(-15600.601952787802_dp)*Hr3(0,0,1) + nf*z**2*(-4543.604938271605_dp)*Hr3(0,0,1) + &
                        nf*z*(-3148.246913580246_dp)*Hr3(0,0,1) + nf*z**3*(-1663.8244170096018_dp)*Hr3(0,0,1) + &
                        nf*(-153.28395061728395_dp)*Hr3(0,0,1) + nf**2*z**3*(30.3758573388203_dp)*Hr3(0,0,1) + &
                        nf**2*z*(53.289437585733886_dp)*Hr3(0,0,1) + nf**2*z**2*(68.21399176954732_dp)*Hr3(0,0,1) + &
                        nf*z**4*(293.5747599451303_dp)*Hr3(0,0,1) + (3496.2962962962956_dp)*Hr3(0,0,1) + &
                        z**2*(20331.412566941402_dp)*Hr3(0,0,1) + z*(21000.06500393373_dp)*Hr3(0,0,1) + &
                        z**4*(-18438.93552812071_dp)*Hr3(0,1,0) + z**3*(-3281.3350844473075_dp)*Hr3(0,1,0) + &
                        nf*z**2*(-2081.8436213991768_dp)*Hr3(0,1,0) + nf*z*(-1567.5171467764062_dp)*Hr3(0,1,0) + &
                        nf*z**3*(-514.5020576131686_dp)*Hr3(0,1,0) + nf*(-320.9657064471879_dp)*Hr3(0,1,0) + &
                        nf**2*z*(5.530864197530864_dp)*Hr3(0,1,0) + nf**2*z**3*(5.530864197530864_dp)*Hr3(0,1,0) + &
                        nf**2*z**2*(11.061728395061728_dp)*Hr3(0,1,0) + nf*z**4*(320.7901234567901_dp)*Hr3(0,1,0) + &
                        z*(10971.197151492337_dp)*Hr3(0,1,0) + (11371.456790123457_dp)*Hr3(0,1,0) + &
                        z**2*(14757.34080504229_dp)*Hr3(0,1,0) + z**3*(-30331.205141889128_dp)*Hr3(0,1,1) + &
                        z**2*(-26337.120846192534_dp)*Hr3(0,1,1) + z**4*(-22751.517146776405_dp)*Hr3(0,1,1) + &
                        z*(-13014.815567129192_dp)*Hr3(0,1,1) + nf*(-126.41975308641973_dp)*Hr3(0,1,1) + &
                        nf**2*z**2*(-110.09053497942385_dp)*Hr3(0,1,1) + nf**2*z*(-55.045267489711925_dp)*Hr3(0,1,1) + &
                        nf**2*z**3*(-55.045267489711925_dp)*Hr3(0,1,1) + nf*z**4*(346.6008230452674_dp)*Hr3(0,1,1) + &
                        nf*z*(424.82304526748976_dp)*Hr3(0,1,1) + nf*z**3*(1112.4938271604938_dp)*Hr3(0,1,1) + &
                        nf*z**2*(1317.1358024691358_dp)*Hr3(0,1,1) + (5742.617283950617_dp)*Hr3(0,1,1) + &
                        z**4*(-12259.672610882486_dp)*Hr3(1,0,0) + z**3*(-3449.752171925011_dp)*Hr3(1,0,0) + &
                        z**2*(-664.2304526748966_dp)*Hr3(1,0,0) + nf*z*(-499.6213991769547_dp)*Hr3(1,0,0) + &
                        nf*(-285.4979423868312_dp)*Hr3(1,0,0) + nf*z**4*(285.4979423868312_dp)*Hr3(1,0,0) + &
                        nf*z**3*(499.6213991769547_dp)*Hr3(1,0,0) + z*(3449.752171925012_dp)*Hr3(1,0,0) + &
                        (12923.903063557384_dp)*Hr3(1,0,0) + z**4*(-10197.947873799723_dp)*Hr3(1,0,1) + &
                        z**2*(-1992.6913580246921_dp)*Hr3(1,0,1) + z**3*(-415.1001371742113_dp)*Hr3(1,0,1) + &
                        nf*z*(-278.69410150891633_dp)*Hr3(1,0,1) + nf*(-159.25377229080934_dp)*Hr3(1,0,1) + &
                        nf*z**4*(159.25377229080934_dp)*Hr3(1,0,1) + nf*z**3*(278.69410150891633_dp)*Hr3(1,0,1) + &
                        z*(415.10013717421_dp)*Hr3(1,0,1) + (12190.639231824416_dp)*Hr3(1,0,1) + &
                        z**4*(-10197.947873799723_dp)*Hr3(1,1,0) + z**2*(-1992.6913580246921_dp)*Hr3(1,1,0) + &
                        z**3*(-415.1001371742113_dp)*Hr3(1,1,0) + nf*z*(-278.69410150891633_dp)*Hr3(1,1,0) + &
                        nf*(-159.25377229080934_dp)*Hr3(1,1,0) + nf*z**4*(159.25377229080934_dp)*Hr3(1,1,0) + &
                        nf*z**3*(278.69410150891633_dp)*Hr3(1,1,0) + z*(415.10013717421_dp)*Hr3(1,1,0) + &
                        (12190.639231824416_dp)*Hr3(1,1,0) + z*(-10358.47462277092_dp)*Hr3(1,1,1) + &
                        (-6745.371742112483_dp)*Hr3(1,1,1) + z**2*(-3985.382716049383_dp)*Hr3(1,1,1) + &
                        nf*z**3*(-860.9711934156377_dp)*Hr3(1,1,1) + nf*z**4*(-491.9835390946502_dp)*Hr3(1,1,1) + &
                        nf*(491.9835390946502_dp)*Hr3(1,1,1) + nf*z*(860.9711934156377_dp)*Hr3(1,1,1) + &
                        z**3*(10358.474622770918_dp)*Hr3(1,1,1) + z**4*(10730.754458161866_dp)*Hr3(1,1,1) + &
                        z**2*(-8010.798353909467_dp)*Hr4(-1,-1,-1,0) + z**3*(-5564.576131687243_dp)*Hr4(-1,-1,-1,0) + &
                        z*(-5564.576131687241_dp)*Hr4(-1,-1,-1,0) + nf*z**3*(404.54320987654313_dp)*Hr4(-1,-1,-1,0) + &
                        nf*z*(404.54320987654324_dp)*Hr4(-1,-1,-1,0) + z**4*(3902.1563786008232_dp)*Hr4(-1,-1,-1,0) + &
                        (3902.1563786008237_dp)*Hr4(-1,-1,-1,0) + (-3561.876543209876_dp)*Hr4(-1,-1,0,0) + &
                        z**4*(-3561.876543209876_dp)*Hr4(-1,-1,0,0) + nf*z*(-337.11934156378595_dp)*Hr4(-1,-1,0,0) + &
                        nf*z**3*(-337.11934156378595_dp)*Hr4(-1,-1,0,0) + z*(4604.576131687243_dp)*Hr4(-1,-1,0,0) + &
                        z**3*(4604.576131687243_dp)*Hr4(-1,-1,0,0) + z**2*(7230.683127572017_dp)*Hr4(-1,-1,0,0) + &
                        (-3221.5967078189306_dp)*Hr4(-1,-1,0,1) + z**4*(-3221.5967078189306_dp)*Hr4(-1,-1,0,1) + &
                        nf*z*(-269.69547325102883_dp)*Hr4(-1,-1,0,1) + nf*z**3*(-269.69547325102883_dp)*Hr4(-1,-1,0,1) + &
                        z*(3644.576131687241_dp)*Hr4(-1,-1,0,1) + z**3*(3644.576131687242_dp)*Hr4(-1,-1,0,1) + &
                        z**2*(6450.567901234567_dp)*Hr4(-1,-1,0,1) + (-1951.0781893004119_dp)*Hr4(-1,0,-1,0) + &
                        z**4*(-1951.0781893004116_dp)*Hr4(-1,0,-1,0) + nf*z*(-202.27160493827162_dp)*Hr4(-1,0,-1,0) + &
                        nf*z**3*(-202.27160493827157_dp)*Hr4(-1,0,-1,0) + z*(2782.2880658436206_dp)*Hr4(-1,0,-1,0) + &
                        z**3*(2782.2880658436216_dp)*Hr4(-1,0,-1,0) + z**2*(4005.3991769547333_dp)*Hr4(-1,0,-1,0) + &
                        z**2*(-3225.2839506172836_dp)*Hr4(-1,0,0,0) + z**3*(-1822.288065843621_dp)*Hr4(-1,0,0,0) + &
                        z*(-1822.2880658436204_dp)*Hr4(-1,0,0,0) + nf*z*(134.84773662551441_dp)*Hr4(-1,0,0,0) + &
                        nf*z**3*(134.84773662551441_dp)*Hr4(-1,0,0,0) + (1610.7983539094653_dp)*Hr4(-1,0,0,0) + &
                        z**4*(1610.7983539094653_dp)*Hr4(-1,0,0,0) + z**2*(-5670.452674897119_dp)*Hr4(-1,0,0,1) + &
                        z*(-2684.5761316872427_dp)*Hr4(-1,0,0,1) + z**3*(-2684.5761316872427_dp)*Hr4(-1,0,0,1) + &
                        nf*z**3*(202.27160493827157_dp)*Hr4(-1,0,0,1) + nf*z*(202.27160493827162_dp)*Hr4(-1,0,0,1) + &
                        (2881.3168724279835_dp)*Hr4(-1,0,0,1) + z**4*(2881.3168724279835_dp)*Hr4(-1,0,0,1) + &
                        z**2*(-2445.1687242798353_dp)*Hr4(-1,0,1,0) + z**3*(-862.2880658436213_dp)*Hr4(-1,0,1,0) + &
                        z*(-862.2880658436211_dp)*Hr4(-1,0,1,0) + nf*z*(67.42386831275721_dp)*Hr4(-1,0,1,0) + &
                        nf*z**3*(67.42386831275721_dp)*Hr4(-1,0,1,0) + (1270.5185185185182_dp)*Hr4(-1,0,1,0) + &
                        z**4*(1270.5185185185182_dp)*Hr4(-1,0,1,0) + z**2*(-4890.3374485596705_dp)*Hr4(-1,0,1,1) + &
                        z**3*(-1724.5761316872427_dp)*Hr4(-1,0,1,1) + z*(-1724.5761316872422_dp)*Hr4(-1,0,1,1) + &
                        nf*z*(134.84773662551441_dp)*Hr4(-1,0,1,1) + nf*z**3*(134.84773662551441_dp)*Hr4(-1,0,1,1) + &
                        (2541.0370370370365_dp)*Hr4(-1,0,1,1) + z**4*(2541.0370370370365_dp)*Hr4(-1,0,1,1) + &
                        z**3*(-3771.7860082304533_dp)*Hr4(0,-1,-1,0) + (-2705.382716049383_dp)*Hr4(0,-1,-1,0) + &
                        z**4*(-1196.77366255144_dp)*Hr4(0,-1,-1,0) + nf*z*(-735.8683127572017_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**3*(331.3251028806584_dp)*Hr4(0,-1,-1,0) + z**2*(4005.3991769547324_dp)*Hr4(0,-1,-1,0) + &
                        z*(9336.362139917694_dp)*Hr4(0,-1,-1,0) + z*(-8384.1316872428_dp)*Hr4(0,-1,0,0) + &
                        z**2*(-3615.3415637860076_dp)*Hr4(0,-1,0,0) + nf*z**3*(-325.2674897119341_dp)*Hr4(0,-1,0,0) + &
                        nf*z*(662.3868312757202_dp)*Hr4(0,-1,0,0) + z**4*(1099.8518518518517_dp)*Hr4(0,-1,0,0) + &
                        (2462.0246913580245_dp)*Hr4(0,-1,0,0) + z**3*(3779.5555555555557_dp)*Hr4(0,-1,0,0) + &
                        z*(-7431.901234567901_dp)*Hr4(0,-1,0,1) + z**2*(-3225.2839506172845_dp)*Hr4(0,-1,0,1) + &
                        nf*z**3*(-319.2098765432099_dp)*Hr4(0,-1,0,1) + nf*z*(588.9053497942386_dp)*Hr4(0,-1,0,1) + &
                        z**4*(1002.9300411522632_dp)*Hr4(0,-1,0,1) + (2218.6666666666665_dp)*Hr4(0,-1,0,1) + &
                        z**3*(3787.325102880658_dp)*Hr4(0,-1,0,1) + nf*z**2*(-924.4444444444442_dp)*Hr4(0,0,-1,0) + &
                        nf*z**3*(-661.59670781893_dp)*Hr4(0,0,-1,0) + nf*z*(-128._dp)*Hr4(0,0,-1,0) + &
                        z**4*(269.69547325102883_dp)*Hr4(0,0,-1,0) + (1024._dp)*Hr4(0,0,-1,0) + &
                        z*(5581.51989026063_dp)*Hr4(0,0,-1,0) + z**3*(12135.593964334706_dp)*Hr4(0,0,-1,0) + &
                        z**2*(18243.862825788754_dp)*Hr4(0,0,-1,0) + nf*z**2*(-2118.584362139917_dp)*Hr4(0,0,0,0) + &
                        nf*z**3*(-1466.469135802469_dp)*Hr4(0,0,0,0) + nf*z*(-719.5390946502058_dp)*Hr4(0,0,0,0) + &
                        z**4*(-382.06858710562415_dp)*Hr4(0,0,0,0) + nf**2*z*(31.604938271604933_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**3*(31.604938271604933_dp)*Hr4(0,0,0,0) + nf**2*z**2*(63.209876543209866_dp)*Hr4(0,0,0,0) + &
                        z*(5157.750342935527_dp)*Hr4(0,0,0,0) + z**3*(17313.53635116598_dp)*Hr4(0,0,0,0) + &
                        z**2*(21943.13305898491_dp)*Hr4(0,0,0,0) + z*(-7919.802469135803_dp)*Hr4(0,0,0,1) + &
                        z**4*(-1297.9094650205761_dp)*Hr4(0,0,0,1) + nf*z**3*(-1098.2716049382716_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(-719.8024691358024_dp)*Hr4(0,0,0,1) + nf**2*z**2*(-12.641975308641973_dp)*Hr4(0,0,0,1) + &
                        nf**2*z*(-6.320987654320986_dp)*Hr4(0,0,0,1) + nf**2*z**3*(-6.320987654320986_dp)*Hr4(0,0,0,1) + &
                        nf*z*(243.6213991769547_dp)*Hr4(0,0,0,1) + (768._dp)*Hr4(0,0,0,1) + z**2*(7365.66255144033_dp)*Hr4(0,0,0,1) &
                        + z**3*(16576._dp)*Hr4(0,0,0,1) + z**4*(-1942.6502057613166_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(-1300.2798353909463_dp)*Hr4(0,0,1,0) + nf*z**3*(-896.7901234567901_dp)*Hr4(0,0,1,0) + &
                        nf*z*(-470.91358024691357_dp)*Hr4(0,0,1,0) + (2360.8888888888887_dp)*Hr4(0,0,1,0) + &
                        z*(8493.563786008228_dp)*Hr4(0,0,1,0) + z**3*(10313.61316872428_dp)*Hr4(0,0,1,0) + &
                        z**2*(17478.716049382714_dp)*Hr4(0,0,1,0) + z**3*(-3243.9835390946505_dp)*Hr4(0,0,1,1) + &
                        z**4*(-3177.349794238683_dp)*Hr4(0,0,1,1) + nf*z**2*(-1384.8230452674895_dp)*Hr4(0,0,1,1) + &
                        nf*z*(-1165.1687242798355_dp)*Hr4(0,0,1,1) + nf*z**3*(-354.50205761316874_dp)*Hr4(0,0,1,1) + &
                        (2958.222222222222_dp)*Hr4(0,0,1,1) + z**2*(17737.48148148148_dp)*Hr4(0,0,1,1) + &
                        z*(22582.781893004118_dp)*Hr4(0,0,1,1) + z**4*(-2083.467764060357_dp)*Hr4(0,1,0,0) + &
                        nf*z**2*(-856.4938271604938_dp)*Hr4(0,1,0,0) + nf*z*(-428.2469135802469_dp)*Hr4(0,1,0,0) + &
                        nf*z**3*(-428.2469135802469_dp)*Hr4(0,1,0,0) + z**3*(2134.9135802469136_dp)*Hr4(0,1,0,0) + &
                        (2812.8395061728397_dp)*Hr4(0,1,0,0) + z*(10703.451303155005_dp)*Hr4(0,1,0,0) + &
                        z**2*(12108.99314128944_dp)*Hr4(0,1,0,0) + z**4*(-3810.5020576131683_dp)*Hr4(0,1,0,1) + &
                        z**3*(-3498.271604938271_dp)*Hr4(0,1,0,1) + nf*z**2*(-477.76131687242787_dp)*Hr4(0,1,0,1) + &
                        nf*z*(-238.88065843621393_dp)*Hr4(0,1,0,1) + nf*z**3*(-238.88065843621393_dp)*Hr4(0,1,0,1) + &
                        (4326.716049382715_dp)*Hr4(0,1,0,1) + z**2*(6727.374485596706_dp)*Hr4(0,1,0,1) + &
                        z*(10741.860082304524_dp)*Hr4(0,1,0,1) + z**4*(-3810.5020576131683_dp)*Hr4(0,1,1,0) + &
                        z**3*(-3498.271604938271_dp)*Hr4(0,1,1,0) + nf*z**2*(-477.76131687242787_dp)*Hr4(0,1,1,0) + &
                        nf*z*(-238.88065843621393_dp)*Hr4(0,1,1,0) + nf*z**3*(-238.88065843621393_dp)*Hr4(0,1,1,0) + &
                        (4326.716049382715_dp)*Hr4(0,1,1,0) + z**2*(6727.374485596706_dp)*Hr4(0,1,1,0) + &
                        z*(10741.860082304524_dp)*Hr4(0,1,1,0) + z**2*(-20985.152263374486_dp)*Hr4(0,1,1,1) + &
                        z**3*(-18646.123456790123_dp)*Hr4(0,1,1,1) + z**4*(-4890.3374485596705_dp)*Hr4(0,1,1,1) + &
                        z*(-3417.810699588478_dp)*Hr4(0,1,1,1) + nf*z*(737.9753086419752_dp)*Hr4(0,1,1,1) + &
                        nf*z**3*(737.9753086419752_dp)*Hr4(0,1,1,1) + nf*z**2*(1475.9506172839504_dp)*Hr4(0,1,1,1) + &
                        (3811.5555555555547_dp)*Hr4(0,1,1,1) + z**3*(-2717.497942386831_dp)*Hr4(1,0,0,0) + &
                        z**4*(-1552.8559670781892_dp)*Hr4(1,0,0,0) + (1552.8559670781892_dp)*Hr4(1,0,0,0) + &
                        z*(2717.497942386831_dp)*Hr4(1,0,0,0) + z**3*(-6030.485596707817_dp)*Hr4(1,0,0,1) + &
                        z**4*(-3445.991769547325_dp)*Hr4(1,0,0,1) + (3445.991769547325_dp)*Hr4(1,0,0,1) + &
                        z*(6030.485596707817_dp)*Hr4(1,0,0,1) + z**3*(-6625.975308641975_dp)*Hr4(1,0,1,0) + &
                        z**4*(-3786.271604938272_dp)*Hr4(1,0,1,0) + (3786.271604938272_dp)*Hr4(1,0,1,0) + &
                        z*(6625.975308641975_dp)*Hr4(1,0,1,0) + z**3*(-6625.975308641975_dp)*Hr4(1,0,1,1) + &
                        z**4*(-3786.271604938272_dp)*Hr4(1,0,1,1) + (3786.271604938272_dp)*Hr4(1,0,1,1) + &
                        z*(6625.975308641975_dp)*Hr4(1,0,1,1) + z**3*(-6030.485596707817_dp)*Hr4(1,1,0,0) + &
                        z**4*(-3445.991769547325_dp)*Hr4(1,1,0,0) + (3445.991769547325_dp)*Hr4(1,1,0,0) + &
                        z*(6030.485596707817_dp)*Hr4(1,1,0,0) + z**3*(-6625.975308641975_dp)*Hr4(1,1,0,1) + &
                        z**4*(-3786.271604938272_dp)*Hr4(1,1,0,1) + (3786.271604938272_dp)*Hr4(1,1,0,1) + &
                        z*(6625.975308641975_dp)*Hr4(1,1,0,1) + z**3*(-6625.975308641975_dp)*Hr4(1,1,1,0) + &
                        z**4*(-3786.271604938272_dp)*Hr4(1,1,1,0) + (3786.271604938272_dp)*Hr4(1,1,1,0) + &
                        z*(6625.975308641975_dp)*Hr4(1,1,1,0) + z*(2157.5637860082306_dp)*Hr5(-1,-1,-1,-1,0) + &
                        z**3*(2157.5637860082306_dp)*Hr5(-1,-1,-1,-1,0) + z*(-1887.8683127572021_dp)*Hr5(-1,-1,-1,0,0) + &
                        z**3*(-1887.8683127572017_dp)*Hr5(-1,-1,-1,0,0) + z*(-1618.172839506173_dp)*Hr5(-1,-1,-1,0,1) + &
                        z**3*(-1618.1728395061725_dp)*Hr5(-1,-1,-1,0,1) + z*(-1078.7818930041153_dp)*Hr5(-1,-1,0,-1,0) + &
                        z**3*(-1078.7818930041153_dp)*Hr5(-1,-1,0,-1,0) + z**3*(809.0864197530863_dp)*Hr5(-1,-1,0,0,0) + &
                        z*(809.0864197530865_dp)*Hr5(-1,-1,0,0,0) + z*(1348.4773662551438_dp)*Hr5(-1,-1,0,0,1) + &
                        z**3*(1348.4773662551438_dp)*Hr5(-1,-1,0,0,1) + z*(539.3909465020577_dp)*Hr5(-1,-1,0,1,0) + &
                        z**3*(539.3909465020577_dp)*Hr5(-1,-1,0,1,0) + z*(1078.7818930041153_dp)*Hr5(-1,-1,0,1,1) + &
                        z**3*(1078.7818930041153_dp)*Hr5(-1,-1,0,1,1) + z*(-1078.7818930041153_dp)*Hr5(-1,0,-1,-1,0) + &
                        z**3*(-1078.7818930041153_dp)*Hr5(-1,0,-1,-1,0) + z**3*(943.9341563786008_dp)*Hr5(-1,0,-1,0,0) + &
                        z*(943.9341563786011_dp)*Hr5(-1,0,-1,0,0) + z**3*(809.0864197530863_dp)*Hr5(-1,0,-1,0,1) + &
                        z*(809.0864197530865_dp)*Hr5(-1,0,-1,0,1) + z*(359.59396433470505_dp)*Hr5(-1,0,0,-1,0) + &
                        z**3*(359.59396433470505_dp)*Hr5(-1,0,0,-1,0) + z*(-224.74622770919066_dp)*Hr5(-1,0,0,0,0) + &
                        z**3*(-224.74622770919063_dp)*Hr5(-1,0,0,0,0) + z*(-539.3909465020577_dp)*Hr5(-1,0,0,0,1) + &
                        z**3*(-539.3909465020577_dp)*Hr5(-1,0,0,0,1) + z*(-404.54320987654324_dp)*Hr5(-1,0,0,1,0) + &
                        z**3*(-404.54320987654313_dp)*Hr5(-1,0,0,1,0) + z*(-809.0864197530865_dp)*Hr5(-1,0,0,1,1) + &
                        z**3*(-809.0864197530863_dp)*Hr5(-1,0,0,1,1) + z*(-89.89849108367626_dp)*Hr5(-1,0,1,0,0) + &
                        z**3*(-89.89849108367626_dp)*Hr5(-1,0,1,0,0) + z*(-269.69547325102883_dp)*Hr5(-1,0,1,0,1) + &
                        z**3*(-269.69547325102883_dp)*Hr5(-1,0,1,0,1) + z*(-269.69547325102883_dp)*Hr5(-1,0,1,1,0) + &
                        z**3*(-269.69547325102883_dp)*Hr5(-1,0,1,1,0) + z*(-539.3909465020577_dp)*Hr5(-1,0,1,1,1) + &
                        z**3*(-539.3909465020577_dp)*Hr5(-1,0,1,1,1) + z*(-6392.625514403292_dp)*Hr5(0,-1,-1,-1,0) + &
                        z**3*(4235.061728395061_dp)*Hr5(0,-1,-1,-1,0) + z**3*(-3926.9135802469136_dp)*Hr5(0,-1,-1,0,0) + &
                        z*(5814.781893004116_dp)*Hr5(0,-1,-1,0,0) + z**3*(-3618.765432098765_dp)*Hr5(0,-1,-1,0,1) + &
                        z*(5236.938271604939_dp)*Hr5(0,-1,-1,0,1) + z**3*(-2117.5308641975307_dp)*Hr5(0,-1,0,-1,0) + &
                        z*(3196.312757201646_dp)*Hr5(0,-1,0,-1,0) + z*(-2618.4691358024693_dp)*Hr5(0,-1,0,0,0) + &
                        z**3*(1809.3827160493825_dp)*Hr5(0,-1,0,0,0) + z*(-4659.09465020576_dp)*Hr5(0,-1,0,0,1) + &
                        z**3*(3310.6172839506166_dp)*Hr5(0,-1,0,0,1) + z*(-2040.6255144032923_dp)*Hr5(0,-1,0,1,0) + &
                        z**3*(1501.2345679012346_dp)*Hr5(0,-1,0,1,0) + z*(-4081.2510288065846_dp)*Hr5(0,-1,0,1,1) + &
                        z**3*(3002.4691358024693_dp)*Hr5(0,-1,0,1,1) + z**2*(-12174.222222222223_dp)*Hr5(0,0,-1,-1,0) + &
                        z**3*(-8384.438957475995_dp)*Hr5(0,0,-1,-1,0) + z*(-3070.595336076817_dp)*Hr5(0,0,-1,-1,0) + &
                        z*(2789.4869684499313_dp)*Hr5(0,0,-1,0,0) + z**3*(7660.3347050754455_dp)*Hr5(0,0,-1,0,0) + &
                        z**2*(11079.111111111111_dp)*Hr5(0,0,-1,0,0) + z*(2508.3786008230454_dp)*Hr5(0,0,-1,0,1) + &
                        z**3*(6936.2304526748985_dp)*Hr5(0,0,-1,0,1) + z**2*(9983.999999999998_dp)*Hr5(0,0,-1,0,1) + &
                        z*(-1155.3360768175583_dp)*Hr5(0,0,0,-1,0) + nf*z**3*(-113.77777777777777_dp)*Hr5(0,0,0,-1,0) + &
                        nf*z*(113.77777777777777_dp)*Hr5(0,0,0,-1,0) + z**2*(6087.111111111111_dp)*Hr5(0,0,0,-1,0) + &
                        z**3*(7062.650205761317_dp)*Hr5(0,0,0,-1,0) + z*(-1007.846364883402_dp)*Hr5(0,0,0,0,0) + &
                        nf*z**3*(-417.18518518518516_dp)*Hr5(0,0,0,0,0) + nf*z**2*(-265.4814814814815_dp)*Hr5(0,0,0,0,0) + &
                        nf*z*(151.7037037037037_dp)*Hr5(0,0,0,0,0) + z**2*(7296.526748971194_dp)*Hr5(0,0,0,0,0) + &
                        z**3*(8394.271604938273_dp)*Hr5(0,0,0,0,0) + nf*z**2*(-594.1728395061729_dp)*Hr5(0,0,0,0,1) + &
                        nf*z*(-297.08641975308643_dp)*Hr5(0,0,0,0,1) + nf*z**3*(-297.08641975308643_dp)*Hr5(0,0,0,0,1) + &
                        z*(3401.218106995884_dp)*Hr5(0,0,0,0,1) + z**3*(15494.847736625512_dp)*Hr5(0,0,0,0,1) + &
                        z**2*(18626.37037037037_dp)*Hr5(0,0,0,0,1) + z*(-2209.7119341563784_dp)*Hr5(0,0,0,1,0) + &
                        nf*z**2*(-164.34567901234567_dp)*Hr5(0,0,0,1,0) + nf*z*(-82.17283950617283_dp)*Hr5(0,0,0,1,0) + &
                        nf*z**3*(-82.17283950617283_dp)*Hr5(0,0,0,1,0) + z**2*(12839.506172839507_dp)*Hr5(0,0,0,1,0) + &
                        z**3*(15318.913580246914_dp)*Hr5(0,0,0,1,0) + z*(-6344.164609053497_dp)*Hr5(0,0,0,1,1) + &
                        z**2*(10935.308641975309_dp)*Hr5(0,0,0,1,1) + z**3*(17818.864197530864_dp)*Hr5(0,0,0,1,1) + &
                        z*(-758.869684499314_dp)*Hr5(0,0,1,0,0) + z**2*(8191.473251028808_dp)*Hr5(0,0,1,0,0) + &
                        z**3*(9040.241426611798_dp)*Hr5(0,0,1,0,0) + z*(2061.1687242798353_dp)*Hr5(0,0,1,0,1) + &
                        z**3*(12955.390946502055_dp)*Hr5(0,0,1,0,1) + z**2*(14746.864197530862_dp)*Hr5(0,0,1,0,1) + &
                        z*(2061.1687242798353_dp)*Hr5(0,0,1,1,0) + z**3*(12955.390946502055_dp)*Hr5(0,0,1,1,0) + &
                        z**2*(14746.864197530862_dp)*Hr5(0,0,1,1,0) + z*(9337.152263374484_dp)*Hr5(0,0,1,1,1) + &
                        z**3*(9337.152263374486_dp)*Hr5(0,0,1,1,1) + z**2*(18134.91358024691_dp)*Hr5(0,0,1,1,1) + &
                        z*(2329.283950617284_dp)*Hr5(0,1,0,0,0) + z**3*(2329.283950617284_dp)*Hr5(0,1,0,0,0) + &
                        z**2*(4658.567901234568_dp)*Hr5(0,1,0,0,0) + z*(5168.987654320987_dp)*Hr5(0,1,0,0,1) + &
                        z**3*(5168.987654320987_dp)*Hr5(0,1,0,0,1) + z**2*(10337.975308641973_dp)*Hr5(0,1,0,0,1) + &
                        z*(5679.407407407408_dp)*Hr5(0,1,0,1,0) + z**3*(5679.407407407408_dp)*Hr5(0,1,0,1,0) + &
                        z**2*(11358.814814814816_dp)*Hr5(0,1,0,1,0) + z*(5679.407407407408_dp)*Hr5(0,1,0,1,1) + &
                        z**3*(5679.407407407408_dp)*Hr5(0,1,0,1,1) + z**2*(11358.814814814816_dp)*Hr5(0,1,0,1,1) + &
                        z*(5168.987654320987_dp)*Hr5(0,1,1,0,0) + z**3*(5168.987654320987_dp)*Hr5(0,1,1,0,0) + &
                        z**2*(10337.975308641973_dp)*Hr5(0,1,1,0,0) + z*(5679.407407407408_dp)*Hr5(0,1,1,0,1) + &
                        z**3*(5679.407407407408_dp)*Hr5(0,1,1,0,1) + z**2*(11358.814814814816_dp)*Hr5(0,1,1,0,1) + &
                        z*(5679.407407407408_dp)*Hr5(0,1,1,1,0) + z**3*(5679.407407407408_dp)*Hr5(0,1,1,1,0) + &
                        z**2*(11358.814814814816_dp)*Hr5(0,1,1,1,0))/(z*(z + (1._dp)))

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
