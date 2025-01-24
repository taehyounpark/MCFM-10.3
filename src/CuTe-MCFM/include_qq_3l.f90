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
                        Ibar_select = z*(-2.6666666666666665_dp) + (2.6666666666666665_dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-2.1932454224643014_dp)
                        else
                            Ibar_select = (-2.1932454224643014_dp)
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
                        Ibar_select = (2.6666666666666665_dp) + z*(2.6666666666666665_dp)
                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-4._dp)
                        else
                            Ibar_select = (-4._dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (-5.333333333333333_dp)
                        else
                            Ibar_select = (-5.333333333333333_dp)
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
                        Ibar_select = (-19.478221870606756_dp) + z**3*((-52.22838795753341_dp) + Nf*(-4.392485119047619_dp)) + &
                        Nf*(-3.654320928943654_dp) + z**4*((-2.2025221668310504_dp) + Nf*(-0.24639893160354862_dp)) + &
                        z**5*(Nf*(0.017839871037076842_dp) + (0.06960260542677314_dp)) + (2.6451727005890224_dp)/z + &
                        z*(Nf*(-1.8770055143244497_dp) + (8.390702182275046_dp)) + z**2*(Nf*(1.658543476725295_dp) + &
                        (63.36166396971877_dp)) + (z*(-16.88888888888889_dp) + (-5.333333333333334_dp) + &
                        z**2*(-1.857655468080786_dp) + Nf*((1.4814814814814814_dp) + z*(1.4814814814814814_dp) + &
                        z**3*(1.9203081586259156_dp) + z**2*(3.503592902185071_dp)) + z**3*(61.37866450247033_dp))*Log(z) + &
                        (z**3*(-18.786707737231765_dp) + z**2*(-12.803702472061493_dp) + z*(-10.666666666666666_dp) + &
                        (-2.6666666666666665_dp) + Nf*(z**3*(-0.053204192649595326_dp) + (0.4444444444444444_dp) + &
                        z*(0.4444444444444444_dp) + z**2*(0.9395466949666464_dp)))*Log(z)**2 + (z**2*(-2.496058052317927_dp) + &
                        (-0.2962962962962963_dp) + z*(-0.2962962962962963_dp) + z**3*(3.4789361919552912_dp))*Log(z)**3 + &
                        (z*(-18.739650524376295_dp) + z**3*(-3.1972624798711755_dp) + (21.104527355585894_dp) + &
                        z**2*(23.0546078708838_dp))*Log(z*(-1._dp) + (1._dp)) + ((-7.259164154264592_dp) + &
                        z**2*(-6.356310832246821_dp) + z**3*(1.3263511167688924_dp) + z*(5.178012758631411_dp))*Log(z*(-1._dp) + &
                        (1._dp))**2

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-30.608535416484305_dp) + Nf*(3.4489756184793317_dp)
                        else
                            Ibar_select = (-30.608535416484305_dp) + Nf*(3.4489756184793317_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = Nf*(5.530864197530864_dp) + (14.926669450170854_dp)
                        else
                            Ibar_select = Nf*(5.530864197530864_dp) + (14.926669450170854_dp)
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
                        Ibar_select = (-7.703703703703703_dp)/z + Nf*(-1.1851842266201185_dp) + z**2*((-15.681439494386987_dp) + &
                        Nf*(-0.46283902782669945_dp)) + (-0.1676001386040955_dp) + z**5*((-0.12440984013915349_dp) + &
                        Nf*(0.07365718799368089_dp)) + z**4*(Nf*(-0.6510300044782803_dp) + (3.674698435447271_dp)) + &
                        z*(Nf*(-4.742561980915629_dp) + (22.721853107547204_dp)) + z**3*(Nf*(-2.5135230066736916_dp) + &
                        (89.38984667302472_dp)) + (z**2*(-75.21947494222778_dp) + z**3*(-43.591358772461135_dp) + &
                        z*(-39.111111111111114_dp) + (-19.555555555555557_dp) + Nf*((1.7777777777777775_dp) + &
                        z*(1.7777777777777775_dp) + z**3*(1.9345473739773027_dp) + z**2*(3.436207859936674_dp)))*Log(z) + &
                        (z**2*(-16.857399697199092_dp) + (-1.7777777777777781_dp) + z*(-1.7777777777777781_dp) + &
                        z**3*(3.654691276642496_dp))*Log(z)**2 + ((-41.11719493959703_dp) + z**2*(-19.226359593565867_dp) + &
                        z**3*(3.2768045492066404_dp) + z*(28.622305539511814_dp))*Log(z*(-1._dp) + (1._dp))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-84.35261963606607_dp) + Nf*(7.75526251932545_dp)
                        else
                            Ibar_select = (-84.35261963606607_dp) + Nf*(7.75526251932545_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (-54.77591205215827_dp) + Nf*(5.9259259259259265_dp)
                        else
                            Ibar_select = (-54.77591205215827_dp) + Nf*(5.9259259259259265_dp)
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
                        Ibar_select = z**2*(-2.2092456459304928_dp) + (-1.7777794164880518_dp) + Nf*((-0.8888888888888888_dp) + &
                        z*(-0.8888888888888888_dp)) + z**5*(-0.5918134963442406_dp) + (1.7777777777777777_dp)/z + &
                        z**4*(3.431416450624558_dp) + z*(10.536329359013553_dp) + z**3*(11.055536391859262_dp) + &
                        (z**2*(-13.889402138049359_dp) + z**3*(-8.607642691980042_dp) + (-0.8888888888888888_dp) + &
                        z*(-0.8888888888888888_dp))*Log(z) + (z*(-16.399098675437834_dp) + (-13.467620810972969_dp) + &
                        z**3*(-0.667673630717109_dp) + z**2*(2.0899486726834704_dp))*Log(z*(-1._dp) + (1._dp))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-37.39461783961921_dp) + Nf*(1.333333333333333_dp)
                        else
                            Ibar_select = (-37.39461783961921_dp) + Nf*(1.333333333333333_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (-8._dp) + Nf*(1.7777777777777777_dp)
                        else
                            Ibar_select = (-8._dp) + Nf*(1.7777777777777777_dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (28.444444444444443_dp)
                        else
                            Ibar_select = (28.444444444444443_dp)
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
                        Ibar_select = z**3*(-3.3427217363336796e8_dp) + z**4*(-3.2431051381686223e6_dp) + (-1080.5982614131724_dp) + &
                        (-472.185801152152_dp)/z + Nf**2*(z**3*(-13.450264046402909_dp) + z**4*(-0.40541711809317443_dp) + &
                        z**5*(-0.00043641441913240813_dp) + (1.6495029510255073_dp) + z*(7.444117645315636_dp) + &
                        z**2*(20.703478227151113_dp)) + Nf*(z**3*(-2275.457084325846_dp) + z*(-219.92491205988716_dp) + &
                        z**5*(-1.27764543675966_dp) + (5.88282387997262_dp)/z + z**4*(82.20918803777823_dp) + (84.97942723038443_dp) &
                        + z**2*(1984.0321035087022_dp)) + z*(6263.120863556335_dp) + z**5*(39341.60322033864_dp) + &
                        z**2*(3.3747238190769726e8_dp) + ((-1429.24047488185_dp) + (-78.98467140066984_dp)/z + &
                        Nf**2*(-5.662551440329218_dp) + Nf*(153.53624971840847_dp) + z*(Nf**2*(-1.1851851851851851_dp) + &
                        Nf*(28.742781690295573_dp) + (164.6569218640053_dp)) + z**2*(Nf**2*(0.6223743940909441_dp) + &
                        Nf*(979.951826373415_dp) + (1.3738623388470685e8_dp)) + z**3*(Nf**2*(13.50872115760657_dp) + &
                        Nf*(1168.676727125948_dp) + (2.032582522664056e8_dp)))*Log(z) + (z**3*(-5.370075316682352e7_dp) + &
                        (-196.22980905904754_dp) + Nf**2*(z**3*(-4.080640889830509_dp) + (-2.4362139917695473_dp) + &
                        z**2*(-2.3073629620695892_dp) + z*(-0.8559670781893004_dp)) + z*(26.241691401506486_dp) + &
                        Nf*(z**3*(-402.55244147369024_dp) + z*(6.6863345960316956_dp) + (47.65684715008643_dp) + &
                        z**2*(210.75077176993884_dp)) + z**2*(2.3897672173141126e7_dp))*Log(z)**2 + ((-115.53991840076135_dp) + &
                        z*(-47.96790193985601_dp) + Nf**2*(z**2*(-0.6104142911440517_dp) + (-0.3292181069958848_dp) + &
                        z*(-0.3292181069958848_dp) + z**3*(0.6660089450144698_dp)) + Nf*((10.326474622770919_dp) + &
                        z*(11.050754458161865_dp) + z**2*(38.63021707930867_dp) + z**3*(56.80190606582247_dp)) + &
                        z**2*(2.2248225933966124e6_dp) + z**3*(9.716442709638055e6_dp))*Log(z)**3 + (z**3*(-887702.7570951866_dp) + &
                        z*(-10.691358024691358_dp) + (-1.530864197530864_dp) + Nf*(z**3*(-10.936522351273572_dp) + &
                        (0.6090534979423868_dp) + z*(0.6090534979423868_dp) + z**2*(3.378615094127529_dp)) + &
                        z**2*(111106.8459163747_dp))*Log(z)**4 + ((-1.0419753086419754_dp) + z*(0.9728395061728395_dp) + &
                        z**2*(2391.7673040141467_dp) + z**3*(85158.56322382427_dp))*Log(z)**5 + (z*(-8229.839494390406_dp) + &
                        Nf*(z**2*(-30.59983085491485_dp) + z*(-0.11155794282490561_dp) + z**3*(5.808624467512806_dp) + &
                        (14.456967794684449_dp)) + z**3*(102.93734727964954_dp) + z**2*(4165.354484787279_dp) + &
                        (4550.651748188859_dp))*Log(z*(-1._dp) + (1._dp)) + (z**2*(-13798.89473059282_dp) + (-4342.16176940399_dp) + &
                        Nf*(z*(-9.224322621215682_dp) + z**3*(-1.7788024634202817_dp) + (10.043747733218687_dp) + &
                        z**2*(11.329747721787648_dp)) + z**3*(4719.296429028051_dp) + z*(13432.266111454075_dp))*Log(z*(-1._dp) + &
                        (1._dp))**2 + ((-48.032152469655344_dp) + z**3*(-3.6747150447818466_dp) + z**2*(1.4407344100318031_dp) + &
                        Nf*(z*(-0.1234260607052864_dp) + z**3*(0.03754220270917316_dp) + z**2*(0.5815314135316819_dp) + &
                        (1.611348329238094_dp)) + z*(15.500701005639957_dp))*Log(z*(-1._dp) + (1._dp))**3 + &
                        (z**2*(-7.504455857694885_dp) + (-2.2065175424656394_dp) + z**3*(2.525512367491166_dp) + &
                        z*(7.185461032669358_dp))*Log(z*(-1._dp) + (1._dp))**4

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-2288.313108178154_dp) + Nf**2*(-8.57238037207872_dp) + Nf*(306.7976732366431_dp)
                        else
                            Ibar_select = (-2288.313108178154_dp) + Nf**2*(-8.57238037207872_dp) + Nf*(306.7976732366431_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = Nf**2*(-9.09324461452157_dp) + (107.39816902904101_dp) + Nf*(142.12668974589386_dp)
                        else
                            Ibar_select = Nf**2*(-9.09324461452157_dp) + (107.39816902904101_dp) + Nf*(142.12668974589386_dp)
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
                        Ibar_select = z**3*(-4.22972933300531e6_dp) + z**5*(-3089.5349177331004_dp) + (-1757.6987405291425_dp) + &
                        Nf**2*(z**2*(-4.422782608695652_dp) + (-0.6584363525049426_dp) + z**5*(-0.04757297210127399_dp) + &
                        z**4*(0.657063817609463_dp) + z*(7.244031151067413_dp) + z**3*(11.713293621097375_dp)) + &
                        Nf*(z**2*(-436.09391015522965_dp) + z*(-237.36245015658233_dp) + (-1.1565265637483262_dp)/z + &
                        z**5*(-0.4200837414838903_dp) + z**4*(25.05741279537468_dp) + (58.80384996366821_dp) + &
                        z**3*(144.25386156284176_dp)) + (483.9673467888446_dp)/z + z*(2185.746709484631_dp) + &
                        z**4*(135352.07357320288_dp) + z**2*(4.098196640506589e6_dp) + ((-801.1344497497917_dp) + &
                        Nf**2*(-6.320987654320988_dp) + z*((-17.561133205221694_dp) + Nf**2*(-1.580246913580247_dp) + &
                        Nf*(9.34091182062968_dp)) + (98.11244277743864_dp)/z + Nf*(147.6125167589013_dp) + &
                        z**2*(Nf**2*(-9.342914380280464_dp) + Nf*(116.84643790984825_dp) + (1.6390928823203878e6_dp)) + &
                        z**3*(Nf*(-284.68715341647203_dp) + Nf**2*(-5.120821746779411_dp) + (2.3298399152559782e6_dp)))*Log(z) + &
                        (z**3*(-737901.8415141026_dp) + (-497.1782563419991_dp) + z*(-89.92539186663392_dp) + &
                        Nf**2*(z**2*(-2.5054578258908435_dp) + (-1.1851851851851851_dp) + z*(-1.1851851851851851_dp) + &
                        z**3*(0.14187787498872553_dp)) + Nf*(z*(33.77777777777778_dp) + (39.7037037037037_dp) + &
                        z**3*(93.75245473222748_dp) + z**2*(95.34625293470691_dp)) + z**2*(265879.8115445273_dp))*Log(z)**2 + &
                        (z*(-66.37037037037037_dp) + (-11.654320987654323_dp) + Nf*(z**3*(-21.23221015364744_dp) + &
                        (2.9629629629629632_dp) + z*(2.9629629629629632_dp) + z**2*(14.542821676253126_dp)) + &
                        z**2*(20689.513723999866_dp) + z**3*(90225.0765892128_dp))*Log(z)**3 + (z**3*(-14263.124546787465_dp) + &
                        (-6.444444444444445_dp) + z*(6.296296296296296_dp) + z**2*(660.4845220177932_dp))*Log(z)**4 + &
                        ((-416.0510248494995_dp) + z**3*(-388.0581448811905_dp) + z*(-293.1852857091508_dp) + &
                        Nf*(z*(-12.011371315711287_dp) + z**3*(-3.8729492332335544_dp) + z**2*(42.8603769974179_dp) + &
                        (81.40254437457219_dp)) + z**2*(802.4939709767361_dp))*Log(z*(-1._dp) + (1._dp)) + ((-283.3406638471292_dp) &
                        + z**2*(-187.35895209869057_dp) + Nf*(z*(-0.05433895247430164_dp) + z**3*(0.2376982680130455_dp) + &
                        z**2*(2.4854368890323975_dp) + (6.812685276910341_dp)) + z**3*(19.63877375541033_dp) + &
                        z*(204.54232367189087_dp))*Log(z*(-1._dp) + (1._dp))**2 + (z*(-61.54102392563138_dp) + &
                        z**3*(-17.881838821881402_dp) + (51.01477773706017_dp) + z**2*(66.33401093637853_dp))*Log(z*(-1._dp) + &
                        (1._dp))**3

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-1788.0038383692386_dp) + Nf**2*(-5.729300528155038_dp) + Nf*(339.3314757224809_dp)
                        else
                            Ibar_select = (-1788.0038383692386_dp) + Nf**2*(-5.729300528155038_dp) + Nf*(339.3314757224809_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (-597.180611891856_dp) + Nf**2*(-6.584362139917696_dp) + Nf*(231.44921352134406_dp)
                        else
                            Ibar_select = (-597.180611891856_dp) + Nf**2*(-6.584362139917696_dp) + Nf*(231.44921352134406_dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (-159.21780746848913_dp) + Nf*(-58.99588477366255_dp)
                        else
                            Ibar_select = (-159.21780746848913_dp) + Nf*(-58.99588477366255_dp)
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
                        Ibar_select = z**2*(-858.5631726043488_dp) + (-316.407615223112_dp) + (-220.95532436198914_dp)/z + &
                        z**5*(-16.93670096591297_dp) + Nf**2*(z**5*(-0.04910479199578725_dp) + z**2*(0.30855935188446637_dp) + &
                        z**4*(0.43402001177163035_dp) + (0.7901228150561073_dp) + z**3*(1.6756820044491276_dp) + &
                        z*(3.1617079872770857_dp)) + Nf*(z**3*(-162.59090205685936_dp) + z*(-55.49584817664684_dp) + &
                        z**4*(-0.31576611913438773_dp) + z**5*(-0.07353941357773612_dp) + (1.1851851851851851_dp)/z + &
                        (11.315211727782213_dp) + z**2*(50.27870465679118_dp)) + z**4*(124.86727654204465_dp) + &
                        z*(346.28439705409136_dp) + z**3*(1499.9296187480459_dp) + ((-476.40998289081455_dp) + &
                        (-46.22222222222222_dp)/z + Nf**2*(-1.1851851851851851_dp) + z*((-247.28235076067622_dp) + &
                        Nf**2*(-1.1851851851851851_dp) + Nf*(31.01234567901234_dp)) + Nf*(35.75308641975309_dp) + &
                        z**3*((-603.2013327681985_dp) + Nf**2*(-1.2896982493182019_dp) + Nf*(75.03343321959325_dp)) + &
                        z**2*((-1056.7021973040094_dp) + Nf**2*(-2.2908052449521357_dp) + Nf*(126.1268812243271_dp)))*Log(z) + &
                        (z**2*(-418.9234294120876_dp) + z*(-93.92592592592591_dp) + z**3*(-35.0728346050264_dp) + (-32._dp) + &
                        Nf*(z**3*(-10.969063474679714_dp) + (3.5555555555555554_dp) + z*(3.5555555555555554_dp) + &
                        z**2*(21.353051658168518_dp)))*Log(z)**2 + (z**2*(-32.233904473757036_dp) + (-9.876543209876543_dp) + &
                        z**3*(-1.5127359104039044_dp) + z*(11.45679012345679_dp))*Log(z)**3 + ((-412.82558000246803_dp) + &
                        z*(-400.5102010164616_dp) + z**2*(-239.38018717024332_dp) + z**3*(7.8317031716114105_dp) + &
                        Nf*(z**3*(-3.2004963392531387_dp) + z**2*(22.1593548666296_dp) + z*(22.36919903758664_dp) + &
                        (59.8077449041727_dp)))*Log(z*(-1._dp) + (1._dp)) + (z*(-140.08680428408184_dp) + &
                        z**3*(-25.516642384848144_dp) + z**2*(127.28132868227607_dp) + (190.02582169035756_dp))*Log(z*(-1._dp) + &
                        (1._dp))**2

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-1357.031985332822_dp) + Nf**2*(-5.170175012883632_dp) + Nf*(190.76761759953948_dp)
                        else
                            Ibar_select = (-1357.031985332822_dp) + Nf**2*(-5.170175012883632_dp) + Nf*(190.76761759953948_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (-223.096709685803_dp) + Nf**2*(-3.950617283950617_dp) + Nf*(74.31423683089888_dp)
                        else
                            Ibar_select = (-223.096709685803_dp) + Nf**2*(-3.950617283950617_dp) + Nf*(74.31423683089888_dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = Nf*(-63.20987654320987_dp) + (646.6620427953394_dp)
                        else
                            Ibar_select = Nf*(-63.20987654320987_dp) + (646.6620427953394_dp)
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
                        Ibar_select = (-201.3875074710237_dp) + z**2*(-38.23462815087049_dp) + z**4*(-22.37170466331597_dp) + &
                        z*(-12.51025086819331_dp) + Nf**2*((0.3950617283950617_dp) + z*(0.3950617283950617_dp)) + &
                        Nf*(z*(-10.579775121412132_dp) + z**3*(-7.370357594572841_dp) + z**4*(-2.2876109733663923_dp) + &
                        (-1.777776679831051_dp) + (-0.7901234567901234_dp)/z + z**5*(0.3945423271814421_dp) + &
                        z**2*(1.0777687022252667_dp)) + z**5*(2.6164131874537233_dp) + (39.111111111111114_dp)/z + &
                        z**3*(230.6719409454128_dp) + (z**2*(-179.05897772016493_dp) + z**3*(-105.1648184234266_dp) + &
                        z*(-24.8888888888889_dp) + (-14.814814814814824_dp) + (7.111111111111111_dp)/z + Nf*((1.1851851851851851_dp) &
                        + z*(1.1851851851851851_dp) + z**3*(5.738428469070934_dp) + z**2*(9.25960142536624_dp)))*Log(z) + &
                        (z**2*(-13.947536596824612_dp) + (-4.54320987654321_dp) + z*(11.456790123456791_dp) + &
                        z**3*(30.073287001340642_dp))*Log(z)**2 + ((-114.52294112178362_dp) + z*(-101.8766124927582_dp) + &
                        z**2*(-74.30758127876247_dp) + Nf*(z**2*(-1.3932991005450936_dp) + z**3*(0.445115753811406_dp) + &
                        (8.9784138885592_dp) + z*(10.932732421137453_dp)) + z**3*(15.744171930341277_dp))*Log(z*(-1._dp) + (1._dp)) &
                        + (z**3*(-2.6110727599061163_dp) + z**2*(7.681222264632603_dp) + z*(30.396699676379065_dp) + &
                        (40.385002670746296_dp))*Log(z*(-1._dp) + (1._dp))**2

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-308.5478196346396_dp) + Nf**2*(-0.5925925925925923_dp) + Nf*(29.81863411530169_dp)
                        else
                            Ibar_select = (-308.5478196346396_dp) + Nf**2*(-0.5925925925925923_dp) + Nf*(29.81863411530169_dp)
                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = Nf**2*(-0.7901234567901234_dp) + Nf*(11.851851851851851_dp) + (101.66018403352467_dp)
                        else
                            Ibar_select = Nf**2*(-0.7901234567901234_dp) + Nf*(11.851851851851851_dp) + (101.66018403352467_dp)
                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = Nf*(-18.962962962962962_dp) + (199.11111111111114_dp)
                        else
                            Ibar_select = Nf*(-18.962962962962962_dp) + (199.11111111111114_dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-75.85185185185185_dp)
                        else
                            Ibar_select = (-75.85185185185185_dp)
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
                        Ibar_select = (z**2*(-147035.5065162744_dp) + z**4*(-122055.30143098741_dp) + (-22375.471060762604_dp) + &
                        nf*z**3*(-10164.535188631842_dp) + nf*z*(-5952.682145979667_dp) + nf**2*z**2*(-651.9106286860859_dp) + &
                        nf**2*z**4*(-10.534979423868313_dp) + nf**3*z**3*(-7.242798353909465_dp) + nf**3*z*(-0.6584362139917695_dp) &
                        + nf**2*(-0.3664031069582032_dp) + nf**3*z**2*(7.901234567901234_dp) + nf**2*z*(94.80031125221302_dp) + &
                        nf*(386.14301560237607_dp) + nf**2*z**3*(552.8153262753732_dp) + nf*z**4*(1339.5666318932692_dp) + &
                        nf*z**2*(15392.049477971586_dp) + z*(113519.71512152723_dp) + z**3*(164453.32488050626_dp) + &
                        (z*(-13041.405403074034_dp) + z**2*(-2353.108644368366_dp) + (-155.912433206915_dp) + &
                        nf*(z*(-269.038105155621_dp) + z**4*(-140.36770703771532_dp) + (140.36770703771532_dp) + &
                        z**3*(269.038105155621_dp)) + z**4*(2509.0210775752817_dp) + z**3*(13041.405403074034_dp))*Hr1(-1) + &
                        z**3*(-70980.43336036279_dp)*Hr1(0) + z*(-33546.43826627107_dp)*Hr1(0) + (-5458.699640174869_dp)*Hr1(0) + &
                        nf*z**3*(-5413.4364780884525_dp)*Hr1(0) + nf*z*(-3298.6142760600687_dp)*Hr1(0) + &
                        nf*z**4*(-1003.7689554734017_dp)*Hr1(0) + nf**2*z**2*(-329.08641975308643_dp)*Hr1(0) + &
                        nf**3*z*(-6.320987654320987_dp)*Hr1(0) + nf**3*z**3*(-6.320987654320987_dp)*Hr1(0) + &
                        nf**2*z**4*(2.3703703703703702_dp)*Hr1(0) + nf**3*z**2*(4.7407407407407405_dp)*Hr1(0) + &
                        nf*(19.34231284582478_dp)*Hr1(0) + nf**2*z*(283.1663352876262_dp)*Hr1(0) + &
                        nf**2*z**3*(389.1727501094128_dp)*Hr1(0) + nf*z**2*(5625.885751835184_dp)*Hr1(0) + &
                        z**2*(31976.021666070552_dp)*Hr1(0) + z**4*(85385.19033317783_dp)*Hr1(0) + &
                        z**2*(-54318.00541585909_dp)*Hr1(1) + z**4*(-12891.024691681501_dp)*Hr1(1) + &
                        nf*z**3*(-5268.0193418161625_dp)*Hr1(1) + (-3232.1356802439873_dp)*Hr1(1) + &
                        nf*z*(-1677.8386219948734_dp)*Hr1(1) + nf**2*z**2*(-125.49794238683127_dp)*Hr1(1) + &
                        nf**2*(-2.3703703703703702_dp)*Hr1(1) + nf**2*z**4*(-2.3703703703703702_dp)*Hr1(1) + &
                        nf**2*z*(5.574759945130316_dp)*Hr1(1) + nf**2*z**3*(124.6639231824417_dp)*Hr1(1) + &
                        nf*(208.27758511647247_dp)*Hr1(1) + nf*z**4*(273.2652394374601_dp)*Hr1(1) + &
                        nf*z**2*(6010.323475288482_dp)*Hr1(1) + z*(25336.93221056996_dp)*Hr1(1) + z**3*(54590.19894864441_dp)*Hr1(1) &
                        + z**3*(-4941.463168124015_dp)*Hr2(-1,-1) + (-3108.8847706871766_dp)*Hr2(-1,-1) + &
                        z**4*(3108.8847706871766_dp)*Hr2(-1,-1) + z*(4941.463168124015_dp)*Hr2(-1,-1) + &
                        z**4*(-21192.666996637785_dp)*Hr2(-1,0) + z*(-5745.04852213961_dp)*Hr2(-1,0) + &
                        nf*z*(-162.23868312757202_dp)*Hr2(-1,0) + nf*(-79.53909465020575_dp)*Hr2(-1,0) + &
                        nf*z**2*(-21.333333333333304_dp)*Hr2(-1,0) + nf*z**4*(100.87242798353908_dp)*Hr2(-1,0) + &
                        nf*z**3*(162.238683127572_dp)*Hr2(-1,0) + z**2*(5654.320987654321_dp)*Hr2(-1,0) + &
                        z**3*(5745.048522139608_dp)*Hr2(-1,0) + (15538.346008983466_dp)*Hr2(-1,0) + &
                        z*(-6908.6145196428515_dp)*Hr2(0,-1) + z**4*(-6155.383893802034_dp)*Hr2(0,-1) + &
                        nf*z*(-210.55156055657298_dp)*Hr2(0,-1) + nf*z**3*(-210.55156055657298_dp)*Hr2(0,-1) + &
                        nf*z**2*(514.6815924716228_dp)*Hr2(0,-1) + (2105.51560556573_dp)*Hr2(0,-1) + &
                        z**2*(3262.9598869764436_dp)*Hr2(0,-1) + z**3*(7048.874423672193_dp)*Hr2(0,-1) + &
                        z**2*(-34875.024794928715_dp)*Hr2(0,0) + z**4*(-5536.254502740596_dp)*Hr2(0,0) + &
                        nf*z**3*(-5269.793668337065_dp)*Hr2(0,0) + nf*z*(-2387.535221386207_dp)*Hr2(0,0) + &
                        (-588.6746566646319_dp)*Hr2(0,0) + nf**2*z**2*(-63.60493827160494_dp)*Hr2(0,0) + &
                        nf**2*z**4*(-4.7407407407407405_dp)*Hr2(0,0) + nf**3*z*(-2.3703703703703702_dp)*Hr2(0,0) + &
                        nf**3*z**3*(-2.3703703703703702_dp)*Hr2(0,0) + nf**2*z*(140.37860082304528_dp)*Hr2(0,0) + &
                        nf**2*z**3*(238.22222222222223_dp)*Hr2(0,0) + nf*z**4*(785.2510288065844_dp)*Hr2(0,0) + &
                        nf*z**2*(2574.643931204778_dp)*Hr2(0,0) + z**3*(24725.82777622362_dp)*Hr2(0,0) + &
                        z*(28053.730371909456_dp)*Hr2(0,0) + z**2*(-8610.67553670462_dp)*Hr2(0,1) + &
                        nf*z*(-3354.326489244371_dp)*Hr2(0,1) + nf*z**3*(-1870.712076703011_dp)*Hr2(0,1) + &
                        nf*z**4*(-218.8641975308642_dp)*Hr2(0,1) + nf**2*z**2*(-44.24691358024691_dp)*Hr2(0,1) + &
                        nf**2*z**4*(4.7407407407407405_dp)*Hr2(0,1) + nf*(16.59259259259259_dp)*Hr2(0,1) + &
                        nf**2*z**3*(57.15226337448559_dp)*Hr2(0,1) + nf**2*z*(114.0411522633745_dp)*Hr2(0,1) + &
                        nf*z**2*(1664.1975308641975_dp)*Hr2(0,1) + z**3*(1780.7397019990058_dp)*Hr2(0,1) + &
                        (2153.499958396644_dp)*Hr2(0,1) + z*(12674.158044883017_dp)*Hr2(0,1) + z**4*(16743.967947476474_dp)*Hr2(0,1) &
                        + z**3*(-20748.417665450404_dp)*Hr2(1,0) + z*(-11995.825072857811_dp)*Hr2(1,0) + &
                        nf*z*(-1310.6957891748884_dp)*Hr2(1,0) + nf*z**3*(-1294.8933200390866_dp)*Hr2(1,0) + &
                        nf*(-324.08230452674906_dp)*Hr2(1,0) + nf*z**4*(-312.2304526748971_dp)*Hr2(1,0) + &
                        nf**2*z**2*(-9.481481481481488_dp)*Hr2(1,0) + nf**2*(4.7407407407407405_dp)*Hr2(1,0) + &
                        nf**2*z**4*(4.7407407407407405_dp)*Hr2(1,0) + nf**2*z*(65.84362139917695_dp)*Hr2(1,0) + &
                        nf**2*z**3*(65.84362139917695_dp)*Hr2(1,0) + nf*z**2*(515.9506172839508_dp)*Hr2(1,0) + &
                        z**2*(10354.174257474235_dp)*Hr2(1,0) + (13705.846556057018_dp)*Hr2(1,0) + &
                        z**4*(17803.624333834796_dp)*Hr2(1,0) + (z**2*(-12234.215621028887_dp) + z**4*(-2236.3680245489704_dp) + &
                        (-459.4791356600816_dp) + nf**2*((-4.7407407407407405_dp) + z**4*(-4.7407407407407405_dp) + &
                        z*(1.1851851851851858_dp) + z**3*(1.1851851851851858_dp) + z**2*(7.111111111111111_dp)) + &
                        nf*(z**3*(-763.1275720164608_dp) + z*(-223.20987654320987_dp) + (202.2716049382716_dp) + &
                        z**4*(275.75308641975306_dp) + z**2*(508.3127572016461_dp)) + z*(11015.464412137808_dp) + &
                        z**3*(12181.4256560936_dp))*Hr2(1,1) + (-2972.181069958848_dp)*Hr3(-1,-1,0) + &
                        z**2*(-1197.0370370370367_dp)*Hr3(-1,-1,0) + z*(-975.6707818930036_dp)*Hr3(-1,-1,0) + &
                        nf*z*(-194.37037037037038_dp)*Hr3(-1,-1,0) + nf*z**4*(-94.81481481481482_dp)*Hr3(-1,-1,0) + &
                        nf*(94.81481481481482_dp)*Hr3(-1,-1,0) + nf*z**3*(194.37037037037038_dp)*Hr3(-1,-1,0) + &
                        z**3*(975.6707818930031_dp)*Hr3(-1,-1,0) + z**4*(4169.218106995885_dp)*Hr3(-1,-1,0) + &
                        (-4291.6872427983535_dp)*Hr3(-1,0,0) + z**3*(-4159.736625514403_dp)*Hr3(-1,0,0) + &
                        nf*z*(-68.74074074074075_dp)*Hr3(-1,0,0) + nf*z**4*(-47.40740740740741_dp)*Hr3(-1,0,0) + &
                        nf*(47.40740740740741_dp)*Hr3(-1,0,0) + nf*z**3*(68.74074074074075_dp)*Hr3(-1,0,0) + &
                        z**2*(1014.5185185185184_dp)*Hr3(-1,0,0) + z**4*(3277.1687242798353_dp)*Hr3(-1,0,0) + &
                        z*(4159.736625514403_dp)*Hr3(-1,0,0) + z**3*(-5776.197530864198_dp)*Hr3(-1,0,1) + &
                        z**4*(-502.51851851851853_dp)*Hr3(-1,0,1) + (-329.48148148148147_dp)*Hr3(-1,0,1) + &
                        nf*z**3*(-66.37037037037037_dp)*Hr3(-1,0,1) + nf*(-37.925925925925924_dp)*Hr3(-1,0,1) + &
                        nf*z**4*(37.925925925925924_dp)*Hr3(-1,0,1) + nf*z*(66.37037037037037_dp)*Hr3(-1,0,1) + &
                        z**2*(832.0000000000006_dp)*Hr3(-1,0,1) + z*(5776.197530864197_dp)*Hr3(-1,0,1) + &
                        z**2*(-8827.569131482986_dp)*Hr3(0,-1,-1) + z*(4289.013270596856_dp)*Hr3(0,-1,-1) + &
                        z**3*(4289.013270596857_dp)*Hr3(0,-1,-1) + z**3*(-21330.43377742551_dp)*Hr3(0,-1,0) + &
                        z*(-14860.24447701399_dp)*Hr3(0,-1,0) + nf*z*(-156.8395061728395_dp)*Hr3(0,-1,0) + &
                        nf*z**4*(-74.27160493827161_dp)*Hr3(0,-1,0) + nf*z**3*(-73.8765432098765_dp)*Hr3(0,-1,0) + &
                        nf*(1.5802469135802468_dp)*Hr3(0,-1,0) + nf*z**2*(367.4074074074074_dp)*Hr3(0,-1,0) + &
                        (2090.6666666666665_dp)*Hr3(0,-1,0) + z**4*(7291.522633744856_dp)*Hr3(0,-1,0) + &
                        z**2*(25384.6324783048_dp)*Hr3(0,-1,0) + z*(-10680.942744777265_dp)*Hr3(0,0,-1) + &
                        z**3*(-5771.972101430498_dp)*Hr3(0,0,-1) + z**2*(15163.611463046527_dp)*Hr3(0,0,-1) + &
                        z**4*(-10193.382716049382_dp)*Hr3(0,0,0) + z**2*(-5117.923724693246_dp)*Hr3(0,0,0) + &
                        nf*z**3*(-3540.4773662551443_dp)*Hr3(0,0,0) + z*(-1151.1556165865686_dp)*Hr3(0,0,0) + &
                        nf*z*(-92.77366255144034_dp)*Hr3(0,0,0) + nf**2*z*(24.098765432098773_dp)*Hr3(0,0,0) + &
                        nf*z**4*(36.34567901234568_dp)*Hr3(0,0,0) + nf**2*z**3*(93.62962962962962_dp)*Hr3(0,0,0) + &
                        nf*z**2*(166.32098765432107_dp)*Hr3(0,0,0) + z**3*(33514.139016285466_dp)*Hr3(0,0,0) + &
                        z**2*(-14479.788939633738_dp)*Hr3(0,0,1) + z**4*(-7508.279835390947_dp)*Hr3(0,0,1) + &
                        nf*z**3*(-2284.8395061728397_dp)*Hr3(0,0,1) + nf*z*(-1515.3909465020577_dp)*Hr3(0,0,1) + &
                        nf*z**4*(14.222222222222221_dp)*Hr3(0,0,1) + nf**2*z**3*(36.74074074074075_dp)*Hr3(0,0,1) + &
                        nf**2*z*(50.96296296296296_dp)*Hr3(0,0,1) + nf*z**2*(766.4197530864199_dp)*Hr3(0,0,1) + &
                        (1109.3333333333333_dp)*Hr3(0,0,1) + z*(12799.82491926513_dp)*Hr3(0,0,1) + &
                        z**3*(26825.58749673797_dp)*Hr3(0,0,1) + z**2*(-5868.444444444445_dp)*Hr3(0,1,0) + &
                        z**4*(-2144.7901234567903_dp)*Hr3(0,1,0) + nf*z**3*(-1541.6625514403293_dp)*Hr3(0,1,0) + &
                        nf*z*(-1188.3456790123457_dp)*Hr3(0,1,0) + nf*(-83.75308641975309_dp)*Hr3(0,1,0) + &
                        nf*z**4*(-26.864197530864196_dp)*Hr3(0,1,0) + nf**2*z**3*(27.65432098765432_dp)*Hr3(0,1,0) + &
                        nf**2*z*(27.654320987654323_dp)*Hr3(0,1,0) + nf*z**2*(474.46913580246905_dp)*Hr3(0,1,0) + &
                        (3320.8888888888887_dp)*Hr3(0,1,0) + z**3*(3797.9542419568484_dp)*Hr3(0,1,0) + &
                        z*(13830.373064413969_dp)*Hr3(0,1,0) + z**3*(-2116.622464320164_dp)*Hr3(0,1,1) + &
                        nf*z*(-370.0411522633744_dp)*Hr3(0,1,1) + nf*z**3*(-206.88065843621396_dp)*Hr3(0,1,1) + &
                        nf*z**4*(-124.8395061728395_dp)*Hr3(0,1,1) + nf*(-56.888888888888886_dp)*Hr3(0,1,1) + &
                        nf**2*z*(-16.59259259259259_dp)*Hr3(0,1,1) + nf**2*z**3*(-2.3703703703703702_dp)*Hr3(0,1,1) + &
                        z**2*(28.839506172839783_dp)*Hr3(0,1,1) + nf*z**2*(193.9753086419753_dp)*Hr3(0,1,1) + &
                        (2042.0740740740741_dp)*Hr3(0,1,1) + z**4*(3033.0205761316874_dp)*Hr3(0,1,1) + &
                        z*(6589.753815282146_dp)*Hr3(0,1,1) + z*(-1414.0746783058728_dp)*Hr3(1,0,-1) + &
                        z**3*(-1414.0746783058728_dp)*Hr3(1,0,-1) + z**2*(-5067.654320987653_dp)*Hr3(1,0,0) + &
                        nf*z*(-1342.8806584362142_dp)*Hr3(1,0,0) + nf*z**3*(-1342.8806584362142_dp)*Hr3(1,0,0) + &
                        (-1250.5020576131687_dp)*Hr3(1,0,0) + z**4*(-449.3168724279835_dp)*Hr3(1,0,0) + &
                        nf*z**2*(-147.35802469135797_dp)*Hr3(1,0,0) + nf**2*z**3*(29.62962962962963_dp)*Hr3(1,0,0) + &
                        nf**2*z*(29.629629629629633_dp)*Hr3(1,0,0) + nf*(44.24691358024691_dp)*Hr3(1,0,0) + &
                        nf*z**4*(44.24691358024691_dp)*Hr3(1,0,0) + z**3*(9119.717118689863_dp)*Hr3(1,0,0) + &
                        z*(10785.297365603446_dp)*Hr3(1,0,0) + z**3*(-560.3936809006236_dp)*Hr3(1,0,1) + &
                        nf*z**3*(-538.7325102880659_dp)*Hr3(1,0,1) + nf*z*(-538.7325102880656_dp)*Hr3(1,0,1) + &
                        nf*(-93.23456790123457_dp)*Hr3(1,0,1) + nf*z**4*(-93.23456790123457_dp)*Hr3(1,0,1) + &
                        nf**2*z*(2.3703703703703702_dp)*Hr3(1,0,1) + nf**2*z**3*(2.3703703703703702_dp)*Hr3(1,0,1) + &
                        nf*z**2*(35.55555555555554_dp)*Hr3(1,0,1) + z*(1288.495207988266_dp)*Hr3(1,0,1) + &
                        z**2*(3568.5925925925912_dp)*Hr3(1,0,1) + (3921.514403292181_dp)*Hr3(1,0,1) + &
                        z**4*(4751.14403292181_dp)*Hr3(1,0,1) + z**3*(-4214.454154765671_dp)*Hr3(1,1,0) + &
                        z*(-2365.5652658767813_dp)*Hr3(1,1,0) + nf*z**3*(-534.7818930041152_dp)*Hr3(1,1,0) + &
                        nf*z*(-534.7818930041151_dp)*Hr3(1,1,0) + nf*(-109.03703703703704_dp)*Hr3(1,1,0) + &
                        nf*z**4*(-109.03703703703704_dp)*Hr3(1,1,0) + nf**2*z*(2.3703703703703702_dp)*Hr3(1,1,0) + &
                        nf**2*z**3*(2.3703703703703702_dp)*Hr3(1,1,0) + nf*z**2*(59.25925925925924_dp)*Hr3(1,1,0) + &
                        (4703.736625514403_dp)*Hr3(1,1,0) + z**2*(4986.864197530864_dp)*Hr3(1,1,0) + &
                        z**4*(5533.366255144033_dp)*Hr3(1,1,0) + z*(-1460.345679012346_dp)*Hr3(1,1,1) + &
                        z**3*(-915.1604938271605_dp)*Hr3(1,1,1) + z**4*(-247.7037037037037_dp)*Hr3(1,1,1) + &
                        nf*z**2*(-92.44444444444444_dp)*Hr3(1,1,1) + nf*z*(-15.4074074074074_dp)*Hr3(1,1,1) + &
                        nf*z**3*(-15.4074074074074_dp)*Hr3(1,1,1) + nf*(61.629629629629626_dp)*Hr3(1,1,1) + &
                        nf*z**4*(61.629629629629626_dp)*Hr3(1,1,1) + (214.5185185185185_dp)*Hr3(1,1,1) + &
                        z**2*(2408.6913580246915_dp)*Hr3(1,1,1) + z*(-270.22222222222223_dp)*Hr4(-1,-1,-1,0) + &
                        (-18.962962962962962_dp)*Hr4(-1,-1,-1,0) + z**4*(18.962962962962962_dp)*Hr4(-1,-1,-1,0) + &
                        z**3*(270.22222222222223_dp)*Hr4(-1,-1,-1,0) + z*(-4376.888888888889_dp)*Hr4(-1,-1,0,0) + &
                        z**4*(-2631.1111111111113_dp)*Hr4(-1,-1,0,0) + (2631.1111111111113_dp)*Hr4(-1,-1,0,0) + &
                        z**3*(4376.888888888889_dp)*Hr4(-1,-1,0,0) + z*(-3139.1604938271603_dp)*Hr4(-1,-1,0,1) + &
                        z**4*(-1880.4938271604938_dp)*Hr4(-1,-1,0,1) + (1880.4938271604938_dp)*Hr4(-1,-1,0,1) + &
                        z**3*(3139.1604938271603_dp)*Hr4(-1,-1,0,1) + z*(-1573.9259259259256_dp)*Hr4(-1,0,-1,0) + &
                        z**4*(-986.074074074074_dp)*Hr4(-1,0,-1,0) + (986.074074074074_dp)*Hr4(-1,0,-1,0) + &
                        z**3*(1573.9259259259256_dp)*Hr4(-1,0,-1,0) + z**3*(-3817.876543209876_dp)*Hr4(-1,0,0,0) + &
                        (-2224.9876543209875_dp)*Hr4(-1,0,0,0) + z**4*(2224.9876543209875_dp)*Hr4(-1,0,0,0) + &
                        z*(3817.876543209876_dp)*Hr4(-1,0,0,0) + z**3*(-3837.2345679012346_dp)*Hr4(-1,0,0,1) + &
                        (-2236.0493827160494_dp)*Hr4(-1,0,0,1) + z**4*(2236.0493827160494_dp)*Hr4(-1,0,0,1) + &
                        z*(3837.2345679012346_dp)*Hr4(-1,0,0,1) + z**3*(-1460.148148148148_dp)*Hr4(-1,0,1,0) + &
                        (-834.3703703703702_dp)*Hr4(-1,0,1,0) + z**4*(834.3703703703702_dp)*Hr4(-1,0,1,0) + &
                        z*(1460.148148148148_dp)*Hr4(-1,0,1,0) + z**3*(-353.9753086419753_dp)*Hr4(-1,0,1,1) + &
                        (-202.2716049382716_dp)*Hr4(-1,0,1,1) + z**4*(202.2716049382716_dp)*Hr4(-1,0,1,1) + &
                        z*(353.9753086419753_dp)*Hr4(-1,0,1,1) + z**2*(-2782.814814814815_dp)*Hr4(0,-1,-1,0) + &
                        z**4*(-1902.6172839506173_dp)*Hr4(0,-1,-1,0) + nf*z*(-170.66666666666666_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**3*(-170.66666666666666_dp)*Hr4(0,-1,-1,0) + nf*z**2*(284.44444444444446_dp)*Hr4(0,-1,-1,0) + &
                        (512._dp)*Hr4(0,-1,-1,0) + z*(651.8518518518518_dp)*Hr4(0,-1,-1,0) + &
                        z**3*(4574.024691358024_dp)*Hr4(0,-1,-1,0) + z**2*(-3211.456790123457_dp)*Hr4(0,-1,0,0) + &
                        z**3*(-3164.8395061728393_dp)*Hr4(0,-1,0,0) + (-1792._dp)*Hr4(0,-1,0,0) + nf*z*(-64._dp)*Hr4(0,-1,0,0) + &
                        nf*z**3*(-64._dp)*Hr4(0,-1,0,0) + nf*z**2*(113.77777777777777_dp)*Hr4(0,-1,0,0) + &
                        z**4*(3004.0493827160494_dp)*Hr4(0,-1,0,0) + z*(5000.691358024692_dp)*Hr4(0,-1,0,0) + &
                        z**3*(-3493.925925925926_dp)*Hr4(0,-1,0,1) + (-1024._dp)*Hr4(0,-1,0,1) + &
                        z**2*(-577.58024691358_dp)*Hr4(0,-1,0,1) + nf*z**2*(-170.66666666666666_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(42.666666666666664_dp)*Hr4(0,-1,0,1) + nf*z**3*(42.666666666666664_dp)*Hr4(0,-1,0,1) + &
                        z**4*(2790.716049382716_dp)*Hr4(0,-1,0,1) + z*(3030.1234567901233_dp)*Hr4(0,-1,0,1) + &
                        z**2*(-3958.5185185185173_dp)*Hr4(0,0,-1,0) + z**3*(-615.9012345679014_dp)*Hr4(0,0,-1,0) + &
                        nf*z*(-199.11111111111111_dp)*Hr4(0,0,-1,0) + nf*z**3*(-146.96296296296296_dp)*Hr4(0,0,-1,0) + &
                        nf*z**2*(284.44444444444446_dp)*Hr4(0,0,-1,0) + z**4*(568.8888888888889_dp)*Hr4(0,0,-1,0) + &
                        z*(4434.567901234567_dp)*Hr4(0,0,-1,0) + z**2*(-6246.3209876543215_dp)*Hr4(0,0,0,0) + &
                        nf*z**3*(-1027.7530864197531_dp)*Hr4(0,0,0,0) + z**4*(-859.6543209876543_dp)*Hr4(0,0,0,0) + &
                        nf*z*(-70.51851851851859_dp)*Hr4(0,0,0,0) + nf*z**2*(170.66666666666666_dp)*Hr4(0,0,0,0) + &
                        z*(3650.4691358024693_dp)*Hr4(0,0,0,0) + z**3*(14437.432098765432_dp)*Hr4(0,0,0,0) + &
                        z**4*(-2101.7283950617284_dp)*Hr4(0,0,0,1) + nf*z**3*(-1190.9135802469136_dp)*Hr4(0,0,0,1) + &
                        nf*z*(-283.4567901234568_dp)*Hr4(0,0,0,1) + nf*z**2*(-170.66666666666666_dp)*Hr4(0,0,0,1) + &
                        z*(1347.25925925926_dp)*Hr4(0,0,0,1) + z**2*(2074.8641975308633_dp)*Hr4(0,0,0,1) + &
                        z**3*(17137.679012345678_dp)*Hr4(0,0,0,1) + z**4*(-1417.4814814814815_dp)*Hr4(0,0,1,0) + &
                        nf*z**3*(-691.358024691358_dp)*Hr4(0,0,1,0) + nf*z*(-430.61728395061726_dp)*Hr4(0,0,1,0) + &
                        (-256._dp)*Hr4(0,0,1,0) + nf*z**2*(-142.22222222222217_dp)*Hr4(0,0,1,0) + &
                        z**2*(1492.5432098765432_dp)*Hr4(0,0,1,0) + z*(4990.024691358026_dp)*Hr4(0,0,1,0) + &
                        z**3*(8692.543209876543_dp)*Hr4(0,0,1,0) + z**2*(-4047.4074074074074_dp)*Hr4(0,0,1,1) + &
                        z**4*(-1197.8271604938273_dp)*Hr4(0,0,1,1) + nf*z*(-707.9506172839506_dp)*Hr4(0,0,1,1) + &
                        (-384._dp)*Hr4(0,0,1,1) + nf*z**3*(-215.7037037037037_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(64.00000000000001_dp)*Hr4(0,0,1,1) + z**3*(4149.728395061728_dp)*Hr4(0,0,1,1) + &
                        z*(8988.839506172839_dp)*Hr4(0,0,1,1) + z**4*(-1286.320987654321_dp)*Hr4(0,1,0,0) + &
                        nf*z**3*(-633.8765432098766_dp)*Hr4(0,1,0,0) + (-625.7777777777778_dp)*Hr4(0,1,0,0) + &
                        nf*z*(-266.8641975308642_dp)*Hr4(0,1,0,0) + z**2*(833.975308641975_dp)*Hr4(0,1,0,0) + &
                        z*(4749.135802469135_dp)*Hr4(0,1,0,0) + z**3*(7512.3950617283945_dp)*Hr4(0,1,0,0) + &
                        z**4*(-982.9135802469136_dp)*Hr4(0,1,0,1) + nf*z*(-406.9135802469136_dp)*Hr4(0,1,0,1) + &
                        nf*z**3*(-237.8271604938272_dp)*Hr4(0,1,0,1) + z**2*(-86.12345679012351_dp)*Hr4(0,1,0,1) + &
                        (199.11111111111111_dp)*Hr4(0,1,0,1) + z**3*(2270.41975308642_dp)*Hr4(0,1,0,1) + &
                        z*(5748.543209876543_dp)*Hr4(0,1,0,1) + z**4*(-919.7037037037037_dp)*Hr4(0,1,1,0) + &
                        nf*z*(-430.61728395061726_dp)*Hr4(0,1,1,0) + z**2*(-240.19753086419763_dp)*Hr4(0,1,1,0) + &
                        nf*z**3*(-214.12345679012353_dp)*Hr4(0,1,1,0) + (341.3333333333333_dp)*Hr4(0,1,1,0) + &
                        z**3*(1970.1728395061725_dp)*Hr4(0,1,1,0) + z*(5997.432098765432_dp)*Hr4(0,1,1,0) + &
                        z*(-4810.666666666666_dp)*Hr4(0,1,1,1) + z**3*(-2367.6049382716046_dp)*Hr4(0,1,1,1) + &
                        nf*z**3*(97.18518518518518_dp)*Hr4(0,1,1,1) + nf*z*(282.0740740740741_dp)*Hr4(0,1,1,1) + &
                        z**4*(339.75308641975306_dp)*Hr4(0,1,1,1) + (625.7777777777778_dp)*Hr4(0,1,1,1) + &
                        z**2*(865.1851851851852_dp)*Hr4(0,1,1,1) + z**2*(-2066.9629629629626_dp)*Hr4(1,0,-1,0) + &
                        nf*z*(-90.07407407407408_dp)*Hr4(1,0,-1,0) + nf*z**3*(-90.07407407407408_dp)*Hr4(1,0,-1,0) + &
                        (771.1604938271605_dp)*Hr4(1,0,-1,0) + z**4*(771.1604938271605_dp)*Hr4(1,0,-1,0) + &
                        z*(1217.5802469135801_dp)*Hr4(1,0,-1,0) + z**3*(1217.5802469135801_dp)*Hr4(1,0,-1,0) + &
                        (-1365.3333333333333_dp)*Hr4(1,0,0,0) + z**4*(-1365.3333333333333_dp)*Hr4(1,0,0,0) + &
                        nf*z*(-617.0864197530864_dp)*Hr4(1,0,0,0) + nf*z**3*(-617.0864197530864_dp)*Hr4(1,0,0,0) + &
                        z**2*(878.6172839506168_dp)*Hr4(1,0,0,0) + z*(8699.654320987655_dp)*Hr4(1,0,0,0) + &
                        z**3*(8699.654320987655_dp)*Hr4(1,0,0,0) + (-739.5555555555555_dp)*Hr4(1,0,0,1) + &
                        z**4*(-739.5555555555555_dp)*Hr4(1,0,0,1) + nf*z*(-682.6666666666666_dp)*Hr4(1,0,0,1) + &
                        nf*z**3*(-682.6666666666666_dp)*Hr4(1,0,0,1) + z**2*(-211.75308641975306_dp)*Hr4(1,0,0,1) + &
                        z*(8999.506172839507_dp)*Hr4(1,0,0,1) + z**3*(8999.506172839507_dp)*Hr4(1,0,0,1) + &
                        nf*z*(-401.3827160493828_dp)*Hr4(1,0,1,0) + nf*z**3*(-401.3827160493827_dp)*Hr4(1,0,1,0) + &
                        (-170.66666666666666_dp)*Hr4(1,0,1,0) + z**4*(-170.66666666666666_dp)*Hr4(1,0,1,0) + &
                        z**2*(-47.40740740740778_dp)*Hr4(1,0,1,0) + z*(4882.9629629629635_dp)*Hr4(1,0,1,0) + &
                        z**3*(4882.9629629629635_dp)*Hr4(1,0,1,0) + z*(-1783.703703703704_dp)*Hr4(1,0,1,1) + &
                        z**3*(-1783.703703703704_dp)*Hr4(1,0,1,1) + nf*z*(31.604938271604937_dp)*Hr4(1,0,1,1) + &
                        nf*z**3*(31.604938271604937_dp)*Hr4(1,0,1,1) + (598.9135802469136_dp)*Hr4(1,0,1,1) + &
                        z**4*(598.9135802469136_dp)*Hr4(1,0,1,1) + z**2*(1478.320987654321_dp)*Hr4(1,0,1,1) + &
                        (-401.38271604938274_dp)*Hr4(1,1,0,0) + z**4*(-401.38271604938274_dp)*Hr4(1,1,0,0) + &
                        nf*z*(-145.3827160493827_dp)*Hr4(1,1,0,0) + nf*z**3*(-145.3827160493827_dp)*Hr4(1,1,0,0) + &
                        z*(884.1481481481478_dp)*Hr4(1,1,0,0) + z**3*(884.1481481481482_dp)*Hr4(1,1,0,0) + &
                        z**2*(4135.506172839507_dp)*Hr4(1,1,0,0) + z*(-586.6666666666669_dp)*Hr4(1,1,0,1) + &
                        z**3*(-586.6666666666667_dp)*Hr4(1,1,0,1) + nf*z*(-47.40740740740741_dp)*Hr4(1,1,0,1) + &
                        nf*z**3*(-47.40740740740741_dp)*Hr4(1,1,0,1) + (772.7407407407409_dp)*Hr4(1,1,0,1) + &
                        z**4*(772.7407407407409_dp)*Hr4(1,1,0,1) + z**2*(964.7407407407409_dp)*Hr4(1,1,0,1) + &
                        z**3*(-638.0246913580247_dp)*Hr4(1,1,1,0) + z*(-638.0246913580246_dp)*Hr4(1,1,1,0) + &
                        nf*z*(-47.40740740740741_dp)*Hr4(1,1,1,0) + nf*z**3*(-47.40740740740741_dp)*Hr4(1,1,1,0) + &
                        z**2*(656.5925925925925_dp)*Hr4(1,1,1,0) + (978.1728395061729_dp)*Hr4(1,1,1,0) + &
                        z**4*(978.1728395061729_dp)*Hr4(1,1,1,0) + z**2*(-616.2962962962963_dp)*Hr4(1,1,1,1) + &
                        z*(-102.71604938271608_dp)*Hr4(1,1,1,1) + z**3*(-102.71604938271608_dp)*Hr4(1,1,1,1) + &
                        (410.8641975308642_dp)*Hr4(1,1,1,1) + z**4*(410.8641975308642_dp)*Hr4(1,1,1,1) + &
                        z**3*(-274.962962962963_dp)*Hr5(0,-1,-1,-1,0) + z*(-274.96296296296293_dp)*Hr5(0,-1,-1,-1,0) + &
                        z**2*(-56.88888888888886_dp)*Hr5(0,-1,-1,-1,0) + z*(-3757.037037037037_dp)*Hr5(0,-1,-1,0,0) + &
                        z**3*(-3757.037037037037_dp)*Hr5(0,-1,-1,0,0) + z**2*(7741.62962962963_dp)*Hr5(0,-1,-1,0,0) + &
                        z*(-2744.8888888888887_dp)*Hr5(0,-1,-1,0,1) + z**3*(-2744.8888888888887_dp)*Hr5(0,-1,-1,0,1) + &
                        z**2*(5338.074074074074_dp)*Hr5(0,-1,-1,0,1) + z*(-1327.4074074074074_dp)*Hr5(0,-1,0,-1,0) + &
                        z**3*(-1327.4074074074074_dp)*Hr5(0,-1,0,-1,0) + z**2*(2958.222222222222_dp)*Hr5(0,-1,0,-1,0) + &
                        z**2*(-6523.259259259259_dp)*Hr5(0,-1,0,0,0) + z*(3299.5555555555557_dp)*Hr5(0,-1,0,0,0) + &
                        z**3*(3299.5555555555557_dp)*Hr5(0,-1,0,0,0) + z**2*(-6253.0370370370365_dp)*Hr5(0,-1,0,0,1) + &
                        z*(3392._dp)*Hr5(0,-1,0,0,1) + z**3*(3392._dp)*Hr5(0,-1,0,0,1) + &
                        z**2*(-2199.703703703704_dp)*Hr5(0,-1,0,1,0) + z*(1327.4074074074074_dp)*Hr5(0,-1,0,1,0) + &
                        z**3*(1327.4074074074074_dp)*Hr5(0,-1,0,1,0) + z*(455.1111111111111_dp)*Hr5(0,-1,0,1,1) + &
                        z**3*(455.1111111111111_dp)*Hr5(0,-1,0,1,1) + z*(-3931.6543209876545_dp)*Hr5(0,0,-1,-1,0) + &
                        z**3*(-889.6790123456789_dp)*Hr5(0,0,-1,-1,0) + z**2*(4821.333333333334_dp)*Hr5(0,0,-1,-1,0) + &
                        z**2*(-8336.592592592593_dp)*Hr5(0,0,-1,0,0) + z*(3186.5679012345677_dp)*Hr5(0,0,-1,0,0) + &
                        z**3*(5541.925925925926_dp)*Hr5(0,0,-1,0,0) + z**2*(-6807.703703703704_dp)*Hr5(0,0,-1,0,1) + &
                        z**3*(3064.098765432099_dp)*Hr5(0,0,-1,0,1) + z*(4527.407407407408_dp)*Hr5(0,0,-1,0,1) + &
                        z*(-2620.8395061728397_dp)*Hr5(0,0,0,-1,0) + z**2*(507.2592592592592_dp)*Hr5(0,0,0,-1,0) + &
                        z**3*(2416.987654320988_dp)*Hr5(0,0,0,-1,0) + z**2*(-3460.740740740741_dp)*Hr5(0,0,0,0,0) + &
                        z*(-1285.5308641975307_dp)*Hr5(0,0,0,0,0) + nf*z**3*(-156.44444444444446_dp)*Hr5(0,0,0,0,0) + &
                        nf*z*(156.44444444444446_dp)*Hr5(0,0,0,0,0) + z**3*(7186.172839506173_dp)*Hr5(0,0,0,0,0) + &
                        z**2*(-6089.481481481481_dp)*Hr5(0,0,0,0,1) + z*(2535.1111111111113_dp)*Hr5(0,0,0,0,1) + &
                        z**3*(10381.037037037036_dp)*Hr5(0,0,0,0,1) + z**2*(-3038.8148148148157_dp)*Hr5(0,0,0,1,0) + &
                        z*(1810.9629629629626_dp)*Hr5(0,0,0,1,0) + z**3*(9116.444444444445_dp)*Hr5(0,0,0,1,0) + &
                        z**2*(576.0000000000017_dp)*Hr5(0,0,0,1,1) + z*(612.3456790123456_dp)*Hr5(0,0,0,1,1) + &
                        z**3*(9026.37037037037_dp)*Hr5(0,0,0,1,1) + z**2*(-2925.0370370370374_dp)*Hr5(0,0,1,0,0) + &
                        z*(1306.469135802469_dp)*Hr5(0,0,1,0,0) + z**3*(7547.6543209876545_dp)*Hr5(0,0,1,0,0) + &
                        z*(984.4938271604951_dp)*Hr5(0,0,1,0,1) + z**2*(1130.6666666666665_dp)*Hr5(0,0,1,0,1) + &
                        z**3*(5116.049382716049_dp)*Hr5(0,0,1,0,1) + z*(818.5679012345681_dp)*Hr5(0,0,1,1,0) + &
                        z**2*(1770.6666666666672_dp)*Hr5(0,0,1,1,0) + z**3*(4641.975308641975_dp)*Hr5(0,0,1,1,0) + &
                        z**2*(320.00000000000006_dp)*Hr5(0,0,1,1,1) + z**3*(448._dp)*Hr5(0,0,1,1,1) + &
                        z*(2872.8888888888887_dp)*Hr5(0,0,1,1,1) + z**3*(-474.0740740740741_dp)*Hr5(0,1,0,-1,0) + &
                        z*(1384.2962962962963_dp)*Hr5(0,1,0,-1,0) + z*(-268.64197530864203_dp)*Hr5(0,1,0,0,0) + &
                        z**3*(5805.827160493827_dp)*Hr5(0,1,0,0,0) + z*(2066.9629629629635_dp)*Hr5(0,1,0,0,1) + &
                        z**3*(7085.827160493827_dp)*Hr5(0,1,0,0,1) + z*(2537.8765432098767_dp)*Hr5(0,1,0,1,0) + &
                        z**3*(4693.333333333333_dp)*Hr5(0,1,0,1,0) + z**3*(1225.4814814814815_dp)*Hr5(0,1,0,1,1) + &
                        z*(2718.814814814815_dp)*Hr5(0,1,0,1,1) + z*(1360.5925925925922_dp)*Hr5(0,1,1,0,0) + &
                        z**3*(2937.6790123456794_dp)*Hr5(0,1,1,0,0) + z**3*(1154.3703703703704_dp)*Hr5(0,1,1,0,1) + &
                        z*(3042.765432098766_dp)*Hr5(0,1,1,0,1) + z**3*(846.2222222222222_dp)*Hr5(0,1,1,1,0) + &
                        z*(3350.913580246914_dp)*Hr5(0,1,1,1,0) + z**3*(-1829.9259259259259_dp)*Hr5(0,1,1,1,1) + &
                        z*(-597.3333333333333_dp)*Hr5(0,1,1,1,1) + z*(-303.4074074074074_dp)*Hr5(1,0,-1,-1,0) + &
                        z**3*(-303.4074074074074_dp)*Hr5(1,0,-1,-1,0) + z*(505.679012345679_dp)*Hr5(1,0,-1,0,0) + &
                        z**3*(505.679012345679_dp)*Hr5(1,0,-1,0,0) + z*(707.9506172839506_dp)*Hr5(1,0,-1,0,1) + &
                        z**3*(707.9506172839506_dp)*Hr5(1,0,-1,0,1) + z*(606.8148148148148_dp)*Hr5(1,0,0,-1,0) + &
                        z**3*(606.8148148148148_dp)*Hr5(1,0,0,-1,0) + z*(2768.5925925925926_dp)*Hr5(1,0,0,0,0) + &
                        z**3*(2768.5925925925926_dp)*Hr5(1,0,0,0,0) + z*(5941.728395061728_dp)*Hr5(1,0,0,0,1) + &
                        z**3*(5941.728395061728_dp)*Hr5(1,0,0,0,1) + z**3*(5233.7777777777765_dp)*Hr5(1,0,0,1,0) + &
                        z*(5233.777777777778_dp)*Hr5(1,0,0,1,0) + z*(5056.79012345679_dp)*Hr5(1,0,0,1,1) + &
                        z**3*(5056.79012345679_dp)*Hr5(1,0,0,1,1) + z*(2275.5555555555557_dp)*Hr5(1,0,1,0,0) + &
                        z**3*(2275.5555555555557_dp)*Hr5(1,0,1,0,0) + z*(2477.8271604938273_dp)*Hr5(1,0,1,0,1) + &
                        z**3*(2477.8271604938273_dp)*Hr5(1,0,1,0,1) + z*(2477.8271604938273_dp)*Hr5(1,0,1,1,0) + &
                        z**3*(2477.8271604938273_dp)*Hr5(1,0,1,1,0) + z*(-455.1111111111111_dp)*Hr5(1,0,1,1,1) + &
                        z**3*(-455.1111111111111_dp)*Hr5(1,0,1,1,1) + z*(910.2222222222222_dp)*Hr5(1,1,0,-1,0) + &
                        z**3*(910.2222222222222_dp)*Hr5(1,1,0,-1,0) + z*(2882.3703703703704_dp)*Hr5(1,1,0,0,0) + &
                        z**3*(2882.3703703703704_dp)*Hr5(1,1,0,0,0) + z*(4247.7037037037035_dp)*Hr5(1,1,0,0,1) + &
                        z**3*(4247.7037037037035_dp)*Hr5(1,1,0,0,1) + z*(2730.6666666666665_dp)*Hr5(1,1,0,1,0) + &
                        z**3*(2730.6666666666665_dp)*Hr5(1,1,0,1,0) + z*(50.5679012345679_dp)*Hr5(1,1,0,1,1) + &
                        z**3*(50.5679012345679_dp)*Hr5(1,1,0,1,1) + z*(-606.8148148148148_dp)*Hr5(1,1,1,0,0) + &
                        z**3*(-606.8148148148148_dp)*Hr5(1,1,1,0,0) + z*(303.4074074074074_dp)*Hr5(1,1,1,0,1) + &
                        z**3*(303.4074074074074_dp)*Hr5(1,1,1,0,1) + z*(303.4074074074074_dp)*Hr5(1,1,1,1,0) + &
                        z**3*(303.4074074074074_dp)*Hr5(1,1,1,1,0))/(z*((-1._dp) + z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-42143.16101826339_dp) + nf**2*(-560.3116650738498_dp) + nf**3*(5.72930052815504_dp) + &
                            nf*(10845.048595100448_dp)

                        else
                            Ibar_select = (-42143.16101826339_dp) + nf**2*(-560.3116650738498_dp) + nf**3*(5.72930052815504_dp) + &
                            nf*(10845.048595100448_dp)

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = (-3487.836384444294_dp) + nf**2*(-359.6958029239425_dp) + nf**3*(6.584362139917696_dp) + &
                            nf*(3118.727692513168_dp)

                        else
                            Ibar_select = (-3487.836384444294_dp) + nf**2*(-359.6958029239425_dp) + nf**3*(6.584362139917696_dp) + &
                            nf*(3118.727692513168_dp)

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-3626.131337105237_dp) + nf**2*(125.0150891632373_dp) + (10783.560477337836_dp)
                        else
                            Ibar_select = nf*(-3626.131337105237_dp) + nf**2*(125.0150891632373_dp) + (10783.560477337836_dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = nf*(235.9835390946502_dp) + (636.8712298739565_dp)
                        else
                            Ibar_select = nf*(235.9835390946502_dp) + (636.8712298739565_dp)
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
                        Ibar_select = (z**3*(-19164.335695824557_dp) + z*(-6949.2465180887875_dp) + nf*z*(-697.3221031459686_dp) + &
                        nf*z**3*(-453.8871120376542_dp) + nf*z**4*(-319.11655903545085_dp) + nf*(-133.0382998603344_dp) + &
                        nf**2*z**2*(-107.98353909465021_dp) + nf**3*z**3*(-2.1069958847736627_dp) + nf**3*z*(0.5267489711934157_dp) &
                        + nf**2*(0.7901234567901234_dp) + nf**2*z**4*(0.7901234567901234_dp) + nf**3*z**2*(1.580246913580247_dp) + &
                        nf**2*z*(5.707513243163054_dp) + nf**2*z**3*(100.69577893790691_dp) + nf*z**2*(1648.9531951473869_dp) + &
                        z**2*(3463.153680377601_dp) + (5155.828741530615_dp) + z**4*(16787.908742689935_dp) + &
                        (z*(-1502.670962531838_dp) + z**4*(-858.6691214467646_dp) + (858.6691214467646_dp) + &
                        z**3*(1502.670962531838_dp))*Hr1(-1) + z**2*(-8870.98211338786_dp)*Hr1(0) + &
                        z**4*(-3702.67341130004_dp)*Hr1(0) + nf*z**3*(-1411.021893910288_dp)*Hr1(0) + &
                        nf*z*(-721.2058692040127_dp)*Hr1(0) + nf*(-25.371742112482853_dp)*Hr1(0) + &
                        nf**2*z**2*(-23.17695473251029_dp)*Hr1(0) + nf**2*z**4*(-1.5802469135802468_dp)*Hr1(0) + &
                        nf**3*z*(-0.7901234567901234_dp)*Hr1(0) + nf**3*z**3*(-0.7901234567901234_dp)*Hr1(0) + &
                        nf**2*z*(40.82304526748971_dp)*Hr1(0) + nf**2*z**3*(70.05761316872427_dp)*Hr1(0) + &
                        nf*z**4*(203.36899862825788_dp)*Hr1(0) + nf*z**2*(897.2831229534097_dp)*Hr1(0) + &
                        (1561.7472233738827_dp)*Hr1(0) + z**3*(5621.596885273833_dp)*Hr1(0) + z*(7030.491200080555_dp)*Hr1(0) + &
                        z*(-5630.423334614203_dp)*Hr1(1) + nf*z**3*(-1003.6784759651185_dp)*Hr1(1) + &
                        z**3*(-859.2866407604407_dp)*Hr1(1) + z**2*(-816.6760699901585_dp)*Hr1(1) + &
                        nf*z**4*(-92.3127572016461_dp)*Hr1(1) + nf*(-78.88065843621399_dp)*Hr1(1) + &
                        nf**2*z*(-24.62551440329218_dp)*Hr1(1) + nf**2*z**2*(-17.11934156378601_dp)*Hr1(1) + &
                        nf**2*(1.5802469135802468_dp)*Hr1(1) + nf**2*z**4*(1.5802469135802468_dp)*Hr1(1) + &
                        nf**2*z**3*(38.584362139917694_dp)*Hr1(1) + nf*z*(566.7402043601805_dp)*Hr1(1) + &
                        nf*z**2*(608.1316872427983_dp)*Hr1(1) + (2554.700222468202_dp)*Hr1(1) + z**4*(4508.543843867379_dp)*Hr1(1) + &
                        (-1479.9012345679012_dp)*Hr2(-1,0) + z**3*(-299.78600823045264_dp)*Hr2(-1,0) + &
                        nf*z*(-58.07407407407407_dp)*Hr2(-1,0) + nf*z**4*(-33.18518518518518_dp)*Hr2(-1,0) + &
                        nf*(33.18518518518518_dp)*Hr2(-1,0) + nf*z**3*(58.07407407407407_dp)*Hr2(-1,0) + &
                        z*(299.78600823045264_dp)*Hr2(-1,0) + z**4*(1479.9012345679012_dp)*Hr2(-1,0) + &
                        z*(-1288.003682170147_dp)*Hr2(0,-1) + z**3*(-1288.003682170147_dp)*Hr2(0,-1) + &
                        z**2*(2576.007364340294_dp)*Hr2(0,-1) + z**4*(-2090.1399176954733_dp)*Hr2(0,0) + &
                        nf*z**3*(-998.9794238683129_dp)*Hr2(0,0) + z**2*(-861.7714861986041_dp)*Hr2(0,0) + &
                        z*(-454.2721279868255_dp)*Hr2(0,0) + nf*z*(-61.10288065843622_dp)*Hr2(0,0) + &
                        nf**2*z*(6.320987654320988_dp)*Hr2(0,0) + nf*z**4*(25.28395061728395_dp)*Hr2(0,0) + &
                        nf**2*z**3*(26.864197530864196_dp)*Hr2(0,0) + nf*z**2*(173.03703703703707_dp)*Hr2(0,0) + &
                        (184.88888888888889_dp)*Hr2(0,0) + z**3*(7118.8710016823115_dp)*Hr2(0,0) + &
                        z**2*(-4099.884773662552_dp)*Hr2(0,1) + z**4*(-2136.230452674897_dp)*Hr2(0,1) + &
                        nf*z**3*(-707.8189300411523_dp)*Hr2(0,1) + nf*z*(-510.4197530864197_dp)*Hr2(0,1) + &
                        (-98.76543209876543_dp)*Hr2(0,1) + nf*(4.7407407407407405_dp)*Hr2(0,1) + &
                        nf**2*z**3*(13.432098765432098_dp)*Hr2(0,1) + nf**2*z*(18.17283950617284_dp)*Hr2(0,1) + &
                        nf*z**4*(27.91769547325103_dp)*Hr2(0,1) + nf*z**2*(193.1851851851852_dp)*Hr2(0,1) + &
                        z*(5001.74875383174_dp)*Hr2(0,1) + z**3*(6287.509105790305_dp)*Hr2(0,1) + &
                        z**2*(-5033.218106995885_dp)*Hr2(1,0) + (-910.4855967078189_dp)*Hr2(1,0) + &
                        z**4*(-700.312757201646_dp)*Hr2(1,0) + nf*z*(-599.5720164609054_dp)*Hr2(1,0) + &
                        nf*z**3*(-599.5720164609054_dp)*Hr2(1,0) + nf*(8.954732510288066_dp)*Hr2(1,0) + &
                        nf*z**4*(8.954732510288066_dp)*Hr2(1,0) + nf**2*z*(15.802469135802466_dp)*Hr2(1,0) + &
                        nf**2*z**3*(15.802469135802466_dp)*Hr2(1,0) + nf*z**2*(188.8395061728395_dp)*Hr2(1,0) + &
                        z**3*(5505.84312810966_dp)*Hr2(1,0) + z*(5926.1888071220055_dp)*Hr2(1,0) + (z*(-2910.0924090967114_dp) + &
                        z**2*(-376.49382716049377_dp) + nf*(z**3*(-243.22633744855966_dp) + (-38.452674897119344_dp) + &
                        z**4*(-38.452674897119344_dp) + z**2*(57.679012345679006_dp) + z*(262.4526748971193_dp)) + &
                        (638.6831275720165_dp) + z**4*(1059.0288065843622_dp) + z**3*(1588.8743021008263_dp))*Hr2(1,1) + &
                        z*(-698.7325102880659_dp)*Hr3(-1,-1,0) + z**4*(-399.275720164609_dp)*Hr3(-1,-1,0) + &
                        (399.275720164609_dp)*Hr3(-1,-1,0) + z**3*(698.7325102880659_dp)*Hr3(-1,-1,0) + &
                        z**3*(-631.440329218107_dp)*Hr3(-1,0,0) + (-360.8230452674897_dp)*Hr3(-1,0,0) + &
                        z**4*(360.8230452674897_dp)*Hr3(-1,0,0) + z*(631.440329218107_dp)*Hr3(-1,0,0) + &
                        z**3*(-564.1481481481483_dp)*Hr3(-1,0,1) + (-322.3703703703704_dp)*Hr3(-1,0,1) + &
                        z**4*(322.3703703703704_dp)*Hr3(-1,0,1) + z*(564.1481481481483_dp)*Hr3(-1,0,1) + &
                        z**2*(-1324.2469135802469_dp)*Hr3(0,-1,0) + (-255.99999999999997_dp)*Hr3(0,-1,0) + &
                        nf*z*(-49.77777777777778_dp)*Hr3(0,-1,0) + nf*z**3*(-49.77777777777778_dp)*Hr3(0,-1,0) + &
                        nf*z**2*(99.55555555555556_dp)*Hr3(0,-1,0) + z**4*(143.27572016460906_dp)*Hr3(0,-1,0) + &
                        z**3*(369.11934156378595_dp)*Hr3(0,-1,0) + z*(1067.851851851852_dp)*Hr3(0,-1,0) + &
                        z**2*(-1830.452674897119_dp)*Hr3(0,0,0) + nf*z**3*(-248.09876543209876_dp)*Hr3(0,0,0) + &
                        z**4*(-168.559670781893_dp)*Hr3(0,0,0) + nf*z*(-25.28395061728395_dp)*Hr3(0,0,0) + &
                        nf*z**2*(52.14814814814816_dp)*Hr3(0,0,0) + z*(896.2633744855967_dp)*Hr3(0,0,0) + &
                        z**3*(3767.045267489711_dp)*Hr3(0,0,0) + z**4*(-383.4732510288066_dp)*Hr3(0,0,1) + &
                        nf*z**3*(-340.5432098765432_dp)*Hr3(0,0,1) + (-170.66666666666666_dp)*Hr3(0,0,1) + &
                        nf*z**2*(-42.66666666666667_dp)*Hr3(0,0,1) + nf*z*(-40.296296296296276_dp)*Hr3(0,0,1) + &
                        z**2*(362.6666666666667_dp)*Hr3(0,0,1) + z*(492.64197530864203_dp)*Hr3(0,0,1) + &
                        z**3*(4790.38683127572_dp)*Hr3(0,0,1) + z**2*(-429.82716049382714_dp)*Hr3(0,1,0) + &
                        z**4*(-405.59670781893004_dp)*Hr3(0,1,0) + (-369.77777777777777_dp)*Hr3(0,1,0) + &
                        nf*z**3*(-266.2716049382716_dp)*Hr3(0,1,0) + nf*z*(-138.2716049382716_dp)*Hr3(0,1,0) + &
                        z*(2120.2962962962965_dp)*Hr3(0,1,0) + z**3*(3939.423868312757_dp)*Hr3(0,1,0) + &
                        z**4*(-524.641975308642_dp)*Hr3(0,1,1) + (-322.3703703703704_dp)*Hr3(0,1,1) + &
                        nf*z*(-259.9506172839506_dp)*Hr3(0,1,1) + z**2*(-173.82716049382717_dp)*Hr3(0,1,1) + &
                        nf*z**3*(-144.59259259259258_dp)*Hr3(0,1,1) + z**3*(2174.0246913580245_dp)*Hr3(0,1,1) + &
                        z*(3701.3333333333335_dp)*Hr3(0,1,1) + z**2*(-659.2263374485595_dp)*Hr3(1,0,0) + &
                        (-313.4156378600823_dp)*Hr3(1,0,0) + z**4*(-313.4156378600823_dp)*Hr3(1,0,0) + &
                        nf*z*(-211.75308641975306_dp)*Hr3(1,0,0) + nf*z**3*(-211.75308641975306_dp)*Hr3(1,0,0) + &
                        z**3*(3188.80658436214_dp)*Hr3(1,0,0) + z*(3188.8065843621407_dp)*Hr3(1,0,0) + &
                        (-351.86831275720164_dp)*Hr3(1,0,1) + z**4*(-351.86831275720164_dp)*Hr3(1,0,1) + &
                        z**2*(-281.283950617284_dp)*Hr3(1,0,1) + nf*z*(-202.2716049382716_dp)*Hr3(1,0,1) + &
                        nf*z**3*(-202.2716049382716_dp)*Hr3(1,0,1) + z*(2919.769547325103_dp)*Hr3(1,0,1) + &
                        z**3*(2919.769547325103_dp)*Hr3(1,0,1) + (-351.86831275720164_dp)*Hr3(1,1,0) + &
                        z**4*(-351.86831275720164_dp)*Hr3(1,1,0) + z**2*(-281.283950617284_dp)*Hr3(1,1,0) + &
                        nf*z*(-202.2716049382716_dp)*Hr3(1,1,0) + nf*z**3*(-202.2716049382716_dp)*Hr3(1,1,0) + &
                        z*(2919.769547325103_dp)*Hr3(1,1,0) + z**3*(2919.769547325103_dp)*Hr3(1,1,0) + &
                        z*(-404.5432098765432_dp)*Hr3(1,1,1) + z**3*(-404.5432098765432_dp)*Hr3(1,1,1) + &
                        z**2*(809.0864197530864_dp)*Hr3(1,1,1) + z*(-598.9135802469136_dp)*Hr4(0,-1,-1,0) + &
                        z**3*(-598.9135802469136_dp)*Hr4(0,-1,-1,0) + z**2*(1197.8271604938273_dp)*Hr4(0,-1,-1,0) + &
                        z**2*(-1082.469135802469_dp)*Hr4(0,-1,0,0) + z*(541.2345679012345_dp)*Hr4(0,-1,0,0) + &
                        z**3*(541.2345679012345_dp)*Hr4(0,-1,0,0) + z**2*(-967.1111111111111_dp)*Hr4(0,-1,0,1) + &
                        z*(483.55555555555554_dp)*Hr4(0,-1,0,1) + z**3*(483.55555555555554_dp)*Hr4(0,-1,0,1) + &
                        z**2*(-598.9135802469136_dp)*Hr4(0,0,-1,0) + z*(-276.5432098765432_dp)*Hr4(0,0,-1,0) + &
                        z**3*(875.4567901234568_dp)*Hr4(0,0,-1,0) + z**2*(-788.5432098765432_dp)*Hr4(0,0,0,0) + &
                        z*(-238.61728395061732_dp)*Hr4(0,0,0,0) + nf*z**3*(-28.44444444444444_dp)*Hr4(0,0,0,0) + &
                        nf*z*(28.44444444444444_dp)*Hr4(0,0,0,0) + z**3*(1482.2716049382716_dp)*Hr4(0,0,0,0) + &
                        z**2*(-1905.777777777778_dp)*Hr4(0,0,0,1) + z*(682.6666666666666_dp)*Hr4(0,0,0,1) + &
                        z**3*(2537.8765432098767_dp)*Hr4(0,0,0,1) + z**2*(-938.6666666666665_dp)*Hr4(0,0,1,0) + &
                        z*(159.60493827160485_dp)*Hr4(0,0,1,0) + z**3*(2447.8024691358023_dp)*Hr4(0,0,1,0) + &
                        z*(-192.7901234567903_dp)*Hr4(0,0,1,1) + z**3*(2721.185185185185_dp)*Hr4(0,0,1,1) + &
                        z*(48.19753086419746_dp)*Hr4(0,1,0,0) + z**3*(1620.5432098765432_dp)*Hr4(0,1,0,0) + &
                        z*(382.41975308641963_dp)*Hr4(0,1,0,1) + z**3*(2044.8395061728395_dp)*Hr4(0,1,0,1) + &
                        z*(382.41975308641963_dp)*Hr4(0,1,1,0) + z**3*(2044.8395061728395_dp)*Hr4(0,1,1,0) + &
                        z*(1213.6296296296296_dp)*Hr4(0,1,1,1) + z**3*(1213.6296296296296_dp)*Hr4(0,1,1,1) + &
                        z*(657.3827160493827_dp)*Hr4(1,0,0,0) + z**3*(657.3827160493827_dp)*Hr4(1,0,0,0) + &
                        z*(1264.1975308641975_dp)*Hr4(1,0,0,1) + z**3*(1264.1975308641975_dp)*Hr4(1,0,0,1) + &
                        z*(1213.6296296296296_dp)*Hr4(1,0,1,0) + z**3*(1213.6296296296296_dp)*Hr4(1,0,1,0) + &
                        z*(1213.6296296296296_dp)*Hr4(1,0,1,1) + z**3*(1213.6296296296296_dp)*Hr4(1,0,1,1) + &
                        z*(1264.1975308641975_dp)*Hr4(1,1,0,0) + z**3*(1264.1975308641975_dp)*Hr4(1,1,0,0) + &
                        z*(1213.6296296296296_dp)*Hr4(1,1,0,1) + z**3*(1213.6296296296296_dp)*Hr4(1,1,0,1) + &
                        z*(1213.6296296296296_dp)*Hr4(1,1,1,0) + z**3*(1213.6296296296296_dp)*Hr4(1,1,1,0))/(z*((-1._dp) + &
                        z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-17884.251125708088_dp) + nf**2*(-202.94844315164937_dp) + nf**3*(3.4467833419224223_dp) &
                            + nf*(3556.8370661328086_dp)

                        else
                            Ibar_select = (-17884.251125708088_dp) + nf**2*(-202.94844315164937_dp) + nf**3*(3.4467833419224223_dp) &
                            + nf*(3556.8370661328086_dp)

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf**2*(-87.1900597482041_dp) + nf**3*(2.633744855967078_dp) + nf*(278.0545709630437_dp) + &
                            (3652.2749041724173_dp)

                        else
                            Ibar_select = nf**2*(-87.1900597482041_dp) + nf**3*(2.633744855967078_dp) + nf*(278.0545709630437_dp) + &
                            (3652.2749041724173_dp)

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-1597.2828778561632_dp) + nf**2*(63.20987654320987_dp) + (8678.823936652116_dp)
                        else
                            Ibar_select = nf*(-1597.2828778561632_dp) + nf**2*(63.20987654320987_dp) + (8678.823936652116_dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-2669.8290346111144_dp) + nf*(252.8395061728395_dp)
                        else
                            Ibar_select = (-2669.8290346111144_dp) + nf*(252.8395061728395_dp)
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
                        Ibar_select = (z**2*(-1967.5219328495803_dp) + z**4*(-690.3278556100735_dp) + (-661.7344338060943_dp) + &
                        nf*z*(-199.84755013431635_dp) + z**3*(-130.47648483896614_dp) + nf*z**3*(-41.82048827446275_dp) + &
                        nf**2*z**2*(-5.201646090534979_dp) + nf**2*z*(-2.4362139917695473_dp) + nf**2*(-0.3950617283950617_dp) + &
                        nf**2*z**4*(-0.3950617283950617_dp) + nf**3*z**3*(-0.19753086419753085_dp) + &
                        nf**3*z*(0.19753086419753085_dp) + nf**2*z**3*(8.42798353909465_dp) + nf*(34.28257887517147_dp) + &
                        nf*z**4*(34.28257887517147_dp) + nf*z**2*(173.10288065843625_dp) + z*(3612.155359790861_dp) + &
                        (z**2*(-362.96544595458573_dp) + (-195.55555555555551_dp) + z**4*(-195.55555555555551_dp) + &
                        nf**2*(z*(0.8559670781893004_dp) + z**3*(4.938271604938271_dp)) + nf*(z**3*(-157.89300411522632_dp) + &
                        z*(-16.855967078189302_dp) + z**4*(4.2139917695473255_dp) + (5.794238683127572_dp) + &
                        z**2*(30.419753086419753_dp)) + z*(156.65555349619424_dp) + z**3*(1128.1979603494092_dp))*Hr1(0) + &
                        z**2*(-1638.3209876543208_dp)*Hr1(1) + (-391.11111111111103_dp)*Hr1(1) + &
                        z**4*(-391.11111111111103_dp)*Hr1(1) + nf*z**3*(-174.74897119341563_dp)*Hr1(1) + &
                        nf**2*z*(-5.794238683127572_dp)*Hr1(1) + nf**2*z**3*(5.794238683127572_dp)*Hr1(1) + &
                        nf*(10.008230452674896_dp)*Hr1(1) + nf*z**4*(10.008230452674896_dp)*Hr1(1) + &
                        nf*z**2*(60.839506172839506_dp)*Hr1(1) + nf*z*(93.89300411522633_dp)*Hr1(1) + &
                        z*(790.4024660647077_dp)*Hr1(1) + z**3*(1630.1407438118354_dp)*Hr1(1) + z**2*(-402.238683127572_dp)*Hr2(0,0) &
                        + nf*z**3*(-37.925925925925924_dp)*Hr2(0,0) + (-21.333333333333332_dp)*Hr2(0,0) + &
                        z**4*(-12.641975308641975_dp)*Hr2(0,0) + nf*z*(-3.950617283950617_dp)*Hr2(0,0) + &
                        nf*z**2*(16.59259259259259_dp)*Hr2(0,0) + z*(132.67489711934158_dp)*Hr2(0,0) + &
                        z**3*(619.5884773662551_dp)*Hr2(0,0) + z**2*(-122.86419753086422_dp)*Hr2(0,1) + &
                        nf*z**3*(-71.90123456790124_dp)*Hr2(0,1) + (-61.62962962962963_dp)*Hr2(0,1) + &
                        z**4*(-44.24691358024691_dp)*Hr2(0,1) + nf*z*(-3.9506172839506135_dp)*Hr2(0,1) + &
                        z*(101.53086419753089_dp)*Hr2(0,1) + z**3*(1075.358024691358_dp)*Hr2(0,1) + &
                        z**2*(-122.86419753086419_dp)*Hr2(1,0) + (-52.93827160493827_dp)*Hr2(1,0) + &
                        z**4*(-52.93827160493827_dp)*Hr2(1,0) + nf*z*(-37.925925925925924_dp)*Hr2(1,0) + &
                        nf*z**3*(-37.925925925925924_dp)*Hr2(1,0) + z*(588.4444444444445_dp)*Hr2(1,0) + &
                        z**3*(588.4444444444445_dp)*Hr2(1,0) + (z**2*(-245.72839506172838_dp) + (-105.87654320987654_dp) + &
                        z**4*(-105.87654320987654_dp) + z*((-719.4074074074073_dp) + nf*(75.85185185185185_dp)) + &
                        z**3*(nf*(-75.85185185185185_dp) + (1176.888888888889_dp)))*Hr2(1,1) + &
                        z**2*(-138.66666666666669_dp)*Hr3(0,0,0) + z*(-20.41152263374486_dp)*Hr3(0,0,0) + &
                        nf*z**3*(-2.3703703703703702_dp)*Hr3(0,0,0) + nf*z*(2.3703703703703702_dp)*Hr3(0,0,0) + &
                        z**3*(192.79012345679013_dp)*Hr3(0,0,0) + z**2*(-277.3333333333333_dp)*Hr3(0,0,1) + &
                        z*(76.11522633744858_dp)*Hr3(0,0,1) + z**3*(336.06584362139915_dp)*Hr3(0,0,1) + &
                        z*(-28.83950617283949_dp)*Hr3(0,1,0) + z**3*(231.1111111111111_dp)*Hr3(0,1,0) + &
                        z*(-57.67901234567898_dp)*Hr3(0,1,1) + z**3*(462.2222222222222_dp)*Hr3(0,1,1) + &
                        z*(67.42386831275721_dp)*Hr3(1,0,0) + z**3*(67.42386831275721_dp)*Hr3(1,0,0) + &
                        z*(202.2716049382716_dp)*Hr3(1,0,1) + z**3*(202.2716049382716_dp)*Hr3(1,0,1) + &
                        z*(202.2716049382716_dp)*Hr3(1,1,0) + z**3*(202.2716049382716_dp)*Hr3(1,1,0) + &
                        z*(-404.5432098765432_dp)*Hr3(1,1,1) + z**3*(404.5432098765432_dp)*Hr3(1,1,1))/(z*((-1._dp) + z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-2552.053064093312_dp) + nf**2*(-20.938548008733754_dp) + nf**3*(0.2962962962962963_dp) + &
                            nf*(426.84577435370863_dp)

                        else
                            Ibar_select = (-2552.053064093312_dp) + nf**2*(-20.938548008733754_dp) + nf**3*(0.2962962962962963_dp) + &
                            nf*(426.84577435370863_dp)

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-46.25277662611731_dp) + nf**2*(-10.864197530864198_dp) + &
                            nf**3*(0.3950617283950617_dp) + (1828.274939124685_dp)

                        else
                            Ibar_select = nf*(-46.25277662611731_dp) + nf**2*(-10.864197530864198_dp) + &
                            nf**3*(0.3950617283950617_dp) + (1828.274939124685_dp)

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-268.641975308642_dp) + nf**2*(11.588477366255145_dp) + (839.7382777471274_dp)
                        else
                            Ibar_select = nf*(-268.641975308642_dp) + nf**2*(11.588477366255145_dp) + (839.7382777471274_dp)
                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-948.1481481481482_dp) + nf*(75.85185185185185_dp)
                        else
                            Ibar_select = (-948.1481481481482_dp) + nf*(75.85185185185185_dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = (134.84773662551441_dp)
                        else
                            Ibar_select = (134.84773662551441_dp)
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
                        Ibar_select = (z**3*(-522028.688047478_dp) + z*(-239415.98359990568_dp) + nf*z**4*(-19010.08237291654_dp) + &
                        nf*z*(-10528.674771623722_dp) + nf*(-5930.765024741579_dp) + nf**2*z**2*(-2172.6732486038095_dp) + &
                        nf**3*z**3*(-92.16539836104523_dp) + nf**3*z*(-1.712973838314615_dp) + nf**4*z**2*(-1.0534979423868314_dp) + &
                        nf**3*(-0.5267489711934157_dp) + nf**3*z**4*(-0.5267489711934157_dp) + nf**4*z*(-0.3511659807956104_dp) + &
                        nf**4*z**3*(1.4046639231824416_dp) + nf**2*(90.31270205667106_dp) + nf**3*z**2*(94.93187014174667_dp) + &
                        nf**2*z**4*(218.97174416874213_dp) + nf**2*z*(750.6434785078349_dp) + nf**2*z**3*(1075.5986326299856_dp) + &
                        nf*z**2*(12627.533035586508_dp) + nf*z**3*(23286.703334953945_dp) + (107237.23285764939_dp) + &
                        z**2*(225437.03321293375_dp) + z**4*(428328.0097519824_dp) + (z**4*(-43386.49376656064_dp) + &
                        z*(-22435.003634855475_dp) + nf*(z**3*(-1609.8240886861433_dp) + (-919.8994792492248_dp) + &
                        z**4*(919.8994792492248_dp) + z*(1609.8240886861433_dp)) + z**3*(22435.003634855475_dp) + &
                        (43386.49376656064_dp))*Hr1(-1) + z**2*(-145696.23985653347_dp)*Hr1(0) + &
                        z**4*(-31316.683156337313_dp)*Hr1(0) + nf*z**3*(-20714.025030014418_dp)*Hr1(0) + &
                        nf*z*(-16194.366420785766_dp)*Hr1(0) + nf*(-1574.4864368705591_dp)*Hr1(0) + &
                        nf**2*z**2*(-1063.8157334229613_dp)*Hr1(0) + nf**2*z**4*(-165.17969821673526_dp)*Hr1(0) + &
                        nf**3*z**3*(-59.434842249657066_dp)*Hr1(0) + nf**3*z*(-37.604023776863286_dp)*Hr1(0) + &
                        nf**4*z*(0.5267489711934157_dp)*Hr1(0) + nf**4*z**3*(0.5267489711934157_dp)*Hr1(0) + &
                        nf**3*z**4*(1.0534979423868314_dp)*Hr1(0) + nf**2*(16.095107453132144_dp)*Hr1(0) + &
                        nf**3*z**2*(18.611796982167355_dp)*Hr1(0) + nf**2*z*(870.6964571923204_dp)*Hr1(0) + &
                        nf**2*z**3*(2014.9521675333667_dp)*Hr1(0) + nf*z**4*(6780.612470334416_dp)*Hr1(0) + &
                        z**3*(15923.611081655856_dp)*Hr1(0) + nf*z**2*(21431.318865538724_dp)*Hr1(0) + (37135.86714933466_dp)*Hr1(0) &
                        + z*(121330.59858639313_dp)*Hr1(0) + z**3*(-192246.82341870404_dp)*Hr1(1) + z*(-152118.8251463271_dp)*Hr1(1) &
                        + nf*z**3*(-12092.013842650538_dp)*Hr1(1) + nf*z**4*(-4521.370452228104_dp)*Hr1(1) + &
                        nf*(-2947.664005040175_dp)*Hr1(1) + nf**2*z**2*(-998.1893004115226_dp)*Hr1(1) + &
                        nf**2*z*(-647.9480744959533_dp)*Hr1(1) + nf**3*z**3*(-32.160951074531326_dp)*Hr1(1) + &
                        nf**3*(-1.0534979423868314_dp)*Hr1(1) + nf**3*z**4*(-1.0534979423868314_dp)*Hr1(1) + &
                        nf**3*z**2*(14.924554183813445_dp)*Hr1(1) + nf**3*z*(19.34339277549154_dp)*Hr1(1) + &
                        nf**2*(79.2318244170096_dp)*Hr1(1) + nf**2*z**4*(84.76268861454047_dp)*Hr1(1) + &
                        nf**2*z**3*(1482.142861875926_dp)*Hr1(1) + nf*z*(8576.547738323861_dp)*Hr1(1) + &
                        nf*z**2*(11308.68986696725_dp)*Hr1(1) + (84150.20390796897_dp)*Hr1(1) + z**2*(116090.35062028855_dp)*Hr1(1) &
                        + z**4*(139018.79422243187_dp)*Hr1(1) + z**3*(-14890.241021253938_dp)*Hr2(-1,-1) + &
                        (-8508.709155002249_dp)*Hr2(-1,-1) + z**4*(8508.709155002249_dp)*Hr2(-1,-1) + &
                        z*(14890.241021253938_dp)*Hr2(-1,-1) + (-32947.38054981388_dp)*Hr2(-1,0) + &
                        z*(-4285.195148936809_dp)*Hr2(-1,0) + nf*z**4*(-1898.9300411522636_dp)*Hr2(-1,0) + &
                        nf*z*(-1667.9067215363511_dp)*Hr2(-1,0) + nf**2*z**3*(-41.48148148148148_dp)*Hr2(-1,0) + &
                        nf**2*(-23.703703703703702_dp)*Hr2(-1,0) + nf**2*z**4*(23.703703703703702_dp)*Hr2(-1,0) + &
                        nf**2*z*(41.48148148148148_dp)*Hr2(-1,0) + nf*z**3*(1667.9067215363511_dp)*Hr2(-1,0) + &
                        nf*(1898.9300411522636_dp)*Hr2(-1,0) + z**3*(4285.195148936809_dp)*Hr2(-1,0) + &
                        z**4*(32947.38054981388_dp)*Hr2(-1,0) + z*(-38643.822491071245_dp)*Hr2(0,-1) + &
                        z**3*(-23753.58146981731_dp)*Hr2(0,-1) + nf*z**2*(-2759.6984377476747_dp)*Hr2(0,-1) + &
                        z**4*(-2634.0606752756453_dp)*Hr2(0,-1) + nf*z*(1379.8492188738373_dp)*Hr2(0,-1) + &
                        nf*z**3*(1379.8492188738373_dp)*Hr2(0,-1) + (5874.648479726604_dp)*Hr2(0,-1) + &
                        z**2*(59156.8161564376_dp)*Hr2(0,-1) + z**4*(-61870.877797739886_dp)*Hr2(0,0) + &
                        z**2*(-24389.900020953926_dp)*Hr2(0,0) + nf*z**3*(-19081.327229861497_dp)*Hr2(0,0) + &
                        nf*z*(-555.320494255282_dp)*Hr2(0,0) + nf*(-168.8230452674897_dp)*Hr2(0,0) + &
                        nf**2*z**2*(-114.52400548696842_dp)*Hr2(0,0) + nf**3*z**3*(-20.016460905349795_dp)*Hr2(0,0) + &
                        nf**2*z**4*(-16.15363511659808_dp)*Hr2(0,0) + nf**3*z*(-4.916323731138546_dp)*Hr2(0,0) + &
                        nf**2*z*(64.22516299034754_dp)*Hr2(0,0) + nf**2*z**3*(1114.8585132785138_dp)*Hr2(0,0) + &
                        nf*z**4*(2268.3859167809787_dp)*Hr2(0,0) + z*(2436.6967042356937_dp)*Hr2(0,0) + &
                        nf*z**2*(5197.669593678116_dp)*Hr2(0,0) + (7227.46389234387_dp)*Hr2(0,0) + &
                        z**3*(109945.76117249667_dp)*Hr2(0,0) + z**2*(-62490.2632643357_dp)*Hr2(0,1) + &
                        z**4*(-48264.867526517206_dp)*Hr2(0,1) + nf*z**3*(-17482.818607251873_dp)*Hr2(0,1) + &
                        nf*z*(-9056.504226624633_dp)*Hr2(0,1) + nf**2*z**2*(-228.87242798353913_dp)*Hr2(0,1) + &
                        nf**2*z**4*(-19.13854595336077_dp)*Hr2(0,1) + nf**3*z*(-13.519890260631001_dp)*Hr2(0,1) + &
                        nf**3*z**3*(-10.359396433470508_dp)*Hr2(0,1) + nf**2*(1.2290809327846366_dp)*Hr2(0,1) + &
                        nf*(29.93689986282579_dp)*Hr2(0,1) + nf**2*z*(656.5048010973937_dp)*Hr2(0,1) + &
                        nf**2*z**3*(822.8696844993142_dp)*Hr2(0,1) + nf*z**4*(2378.534979423868_dp)*Hr2(0,1) + &
                        nf*z**2*(7107.862825788753_dp)*Hr2(0,1) + (10847.4903471162_dp)*Hr2(0,1) + z*(69139.52761326992_dp)*Hr2(0,1) &
                        + z**3*(85146.20900662764_dp)*Hr2(0,1) + z**2*(-83291.91409611185_dp)*Hr2(1,0) + &
                        (-18399.453793724657_dp)*Hr2(1,0) + nf*z*(-13121.559268899307_dp)*Hr2(1,0) + &
                        nf*z**3*(-12875.567499351982_dp)*Hr2(1,0) + z**4*(-7778.8776620374165_dp)*Hr2(1,0) + &
                        nf**2*z**2*(-222.02469135802468_dp)*Hr2(1,0) + nf**3*z*(-11.939643347050755_dp)*Hr2(1,0) + &
                        nf**3*z**3*(-11.939643347050755_dp)*Hr2(1,0) + nf**2*(-3.687242798353909_dp)*Hr2(1,0) + &
                        nf**2*z**4*(-3.687242798353909_dp)*Hr2(1,0) + nf**2*z*(730.9958847736625_dp)*Hr2(1,0) + &
                        nf**2*z**3*(730.9958847736626_dp)*Hr2(1,0) + nf*z**4*(760.3914037494285_dp)*Hr2(1,0) + &
                        nf*(883.387288523091_dp)*Hr2(1,0) + nf*z**2*(7552.175582990397_dp)*Hr2(1,0) + &
                        z**3*(69644.24453790838_dp)*Hr2(1,0) + z*(90885.39680128286_dp)*Hr2(1,0) + (z*(-82244.38784493068_dp) + &
                        z**3*(-17699.700094335294_dp) + nf**2*(z*(-309.7722908093279_dp) + z**2*(-128.08779149519887_dp) + &
                        (36.69684499314129_dp) + z**4*(36.69684499314129_dp) + z**3*(364.4663923182441_dp)) + &
                        nf*(z**3*(-8383.265740741077_dp) + z**4*(-1823.7805212620028_dp) + (-1577.7887517146776_dp) + &
                        z**2*(2621.5418381344316_dp) + z*(9163.293175583327_dp)) + z**2*(23036.684995898104_dp) + &
                        (27184.746729252107_dp) + z**4*(48425.898992626586_dp))*Hr2(1,1) + z**4*(-16857.54732510288_dp)*Hr3(-1,-1,0) &
                        + z*(-4843.588477366257_dp)*Hr3(-1,-1,0) + nf*z**3*(-740.5212620027436_dp)*Hr3(-1,-1,0) + &
                        nf*(-423.1550068587106_dp)*Hr3(-1,-1,0) + nf*z**4*(423.1550068587106_dp)*Hr3(-1,-1,0) + &
                        nf*z*(740.5212620027436_dp)*Hr3(-1,-1,0) + z**3*(4843.588477366257_dp)*Hr3(-1,-1,0) + &
                        (16857.54732510288_dp)*Hr3(-1,-1,0) + (-15342.35390946502_dp)*Hr3(-1,0,0) + &
                        z**3*(-4425.415637860083_dp)*Hr3(-1,0,0) + nf*z*(-674.4581618655693_dp)*Hr3(-1,0,0) + &
                        nf*z**4*(-385.40466392318245_dp)*Hr3(-1,0,0) + nf*(385.40466392318245_dp)*Hr3(-1,0,0) + &
                        nf*z**3*(674.4581618655693_dp)*Hr3(-1,0,0) + z*(4425.415637860083_dp)*Hr3(-1,0,0) + &
                        z**4*(15342.35390946502_dp)*Hr3(-1,0,0) + (-13827.160493827161_dp)*Hr3(-1,0,1) + &
                        z**3*(-4007.2427983539096_dp)*Hr3(-1,0,1) + nf*z*(-608.395061728395_dp)*Hr3(-1,0,1) + &
                        nf*z**4*(-347.65432098765433_dp)*Hr3(-1,0,1) + nf*(347.65432098765433_dp)*Hr3(-1,0,1) + &
                        nf*z**3*(608.395061728395_dp)*Hr3(-1,0,1) + z*(4007.2427983539096_dp)*Hr3(-1,0,1) + &
                        z**4*(13827.160493827161_dp)*Hr3(-1,0,1) + z**2*(-25526.12746500675_dp)*Hr3(0,-1,-1) + &
                        z*(12763.063732503375_dp)*Hr3(0,-1,-1) + z**3*(12763.063732503375_dp)*Hr3(0,-1,-1) + &
                        z**2*(-21667.228069194724_dp)*Hr3(0,-1,0) + (-9837.037037037036_dp)*Hr3(0,-1,0) + &
                        nf*z*(-2066.1728395061727_dp)*Hr3(0,-1,0) + nf*z**3*(-1325.6515775034293_dp)*Hr3(0,-1,0) + &
                        nf*z**4*(-116.58710562414267_dp)*Hr3(0,-1,0) + nf**2*z**2*(-71.11111111111111_dp)*Hr3(0,-1,0) + &
                        nf**2*z*(35.55555555555556_dp)*Hr3(0,-1,0) + nf**2*z**3*(35.55555555555556_dp)*Hr3(0,-1,0) + &
                        nf*(306.5679012345679_dp)*Hr3(0,-1,0) + nf*z**2*(3201.843621399177_dp)*Hr3(0,-1,0) + &
                        z**4*(7020.510288065843_dp)*Hr3(0,-1,0) + z**3*(9820.08317039983_dp)*Hr3(0,-1,0) + &
                        z*(14663.671647766088_dp)*Hr3(0,-1,0) + z**3*(-19599.490945636546_dp)*Hr3(0,0,-1) + &
                        z*(6836.427213133172_dp)*Hr3(0,0,-1) + z**2*(12763.063732503373_dp)*Hr3(0,0,-1) + &
                        z**2*(-44365.42145130579_dp)*Hr3(0,0,0) + z**4*(-9823.45861911294_dp)*Hr3(0,0,0) + &
                        nf*z**3*(-8162.428316445705_dp)*Hr3(0,0,0) + nf*z*(-1962.9599934408868_dp)*Hr3(0,0,0) + &
                        nf**2*z**2*(-23.70370370370372_dp)*Hr3(0,0,0) + nf**2*z*(-2.80932784636488_dp)*Hr3(0,0,0) + &
                        nf*z**4*(137.6570644718793_dp)*Hr3(0,0,0) + nf**2*z**3*(232.99862825788756_dp)*Hr3(0,0,0) + &
                        (554.6666666666666_dp)*Hr3(0,0,0) + nf*z**2*(4540.912772620932_dp)*Hr3(0,0,0) + &
                        z*(24952.718140227953_dp)*Hr3(0,0,0) + z**3*(60059.13646560747_dp)*Hr3(0,0,0) + &
                        z**4*(-20674.194787379973_dp)*Hr3(0,0,1) + z**2*(-16118.213701714683_dp)*Hr3(0,0,1) + &
                        nf*z**3*(-11439.348879743942_dp)*Hr3(0,0,1) + (-3496.296296296296_dp)*Hr3(0,0,1) + &
                        nf*z*(-608.6438042981245_dp)*Hr3(0,0,1) + nf**2*z**2*(22.91358024691358_dp)*Hr3(0,0,1) + &
                        nf**2*z*(57.67901234567899_dp)*Hr3(0,0,1) + nf*(153.28395061728395_dp)*Hr3(0,0,1) + &
                        nf*z**4*(293.5747599451303_dp)*Hr3(0,0,1) + nf**2*z**3*(316.9272976680384_dp)*Hr3(0,0,1) + &
                        nf*z**2*(816.5925925925917_dp)*Hr3(0,0,1) + z*(9258.318852845896_dp)*Hr3(0,0,1) + &
                        z**3*(87655.1531977299_dp)*Hr3(0,0,1) + z**2*(-44381.563786008235_dp)*Hr3(0,1,0) + &
                        z**4*(-18438.935528120714_dp)*Hr3(0,1,0) + (-11371.456790123457_dp)*Hr3(0,1,0) + &
                        nf*z**3*(-9086.814814814814_dp)*Hr3(0,1,0) + nf*z*(-3639.5720164609065_dp)*Hr3(0,1,0) + &
                        nf**2*z*(137.74485596707817_dp)*Hr3(0,1,0) + nf**2*z**3*(244.32373113854598_dp)*Hr3(0,1,0) + &
                        nf*z**4*(320.7901234567901_dp)*Hr3(0,1,0) + nf*(320.9657064471879_dp)*Hr3(0,1,0) + &
                        nf*z**2*(1409.1851851851852_dp)*Hr3(0,1,0) + z*(44203.49810610229_dp)*Hr3(0,1,0) + &
                        z**3*(89503.2126474023_dp)*Hr3(0,1,0) + z**2*(-39832.88888888889_dp)*Hr3(0,1,1) + &
                        z**4*(-22751.517146776405_dp)*Hr3(0,1,1) + nf*z*(-6891.89574759945_dp)*Hr3(0,1,1) + &
                        nf*z**3*(-6693.925925925925_dp)*Hr3(0,1,1) + (-5742.617283950617_dp)*Hr3(0,1,1) + &
                        nf*(126.41975308641975_dp)*Hr3(0,1,1) + nf**2*z**3*(135.98902606310014_dp)*Hr3(0,1,1) + &
                        nf**2*z*(246.07956104252403_dp)*Hr3(0,1,1) + nf*z**4*(346.6008230452675_dp)*Hr3(0,1,1) + &
                        nf*z**2*(1538.3703703703698_dp)*Hr3(0,1,1) + z**3*(67396.78900466385_dp)*Hr3(0,1,1) + &
                        z*(68163.18933409784_dp)*Hr3(0,1,1) + z**2*(-48209.03155006859_dp)*Hr3(1,0,0) + &
                        (-12923.903063557385_dp)*Hr3(1,0,0) + z**4*(-12259.672610882488_dp)*Hr3(1,0,0) + &
                        nf*z**3*(-6418.816643804299_dp)*Hr3(1,0,0) + nf*z*(-6418.816643804297_dp)*Hr3(1,0,0) + &
                        nf**2*z**3*(198.7599451303155_dp)*Hr3(1,0,0) + nf**2*z*(198.75994513031551_dp)*Hr3(1,0,0) + &
                        nf*(285.49794238683126_dp)*Hr3(1,0,0) + nf*z**4*(285.49794238683126_dp)*Hr3(1,0,0) + &
                        nf*z**2*(1482.0960219478743_dp)*Hr3(1,0,0) + z**3*(66488.67460502771_dp)*Hr3(1,0,0) + &
                        z*(67817.13551037751_dp)*Hr3(1,0,0) + z**2*(-50213.53086419753_dp)*Hr3(1,0,1) + &
                        (-12190.639231824416_dp)*Hr3(1,0,1) + z**4*(-10197.947873799725_dp)*Hr3(1,0,1) + &
                        nf*z**3*(-6703.539094650207_dp)*Hr3(1,0,1) + nf*z*(-6703.539094650205_dp)*Hr3(1,0,1) + &
                        nf*(159.25377229080934_dp)*Hr3(1,0,1) + nf*z**4*(159.25377229080934_dp)*Hr3(1,0,1) + &
                        nf**2*z**3*(191.03429355281207_dp)*Hr3(1,0,1) + nf**2*z*(191.03429355281213_dp)*Hr3(1,0,1) + &
                        nf*z**2*(1514.1399176954735_dp)*Hr3(1,0,1) + z**3*(67924.84513645903_dp)*Hr3(1,0,1) + &
                        z*(71910.22785250841_dp)*Hr3(1,0,1) + z**2*(-50213.53086419753_dp)*Hr3(1,1,0) + &
                        (-12190.639231824416_dp)*Hr3(1,1,0) + z**4*(-10197.947873799725_dp)*Hr3(1,1,0) + &
                        nf*z**3*(-6703.539094650207_dp)*Hr3(1,1,0) + nf*z*(-6703.539094650205_dp)*Hr3(1,1,0) + &
                        nf*(159.25377229080934_dp)*Hr3(1,1,0) + nf*z**4*(159.25377229080934_dp)*Hr3(1,1,0) + &
                        nf**2*z**3*(191.03429355281207_dp)*Hr3(1,1,0) + nf**2*z*(191.03429355281213_dp)*Hr3(1,1,0) + &
                        nf*z**2*(1514.1399176954735_dp)*Hr3(1,1,0) + z**3*(67481.21386483367_dp)*Hr3(1,1,0) + &
                        z*(71466.59658088305_dp)*Hr3(1,1,0) + z*(-22750.593296891562_dp)*Hr3(1,1,1) + &
                        z**2*(-2566.5843621399194_dp)*Hr3(1,1,1) + nf*z**3*(-1540.1262002743488_dp)*Hr3(1,1,1) + &
                        nf*(-491.98353909465015_dp)*Hr3(1,1,1) + nf*z**4*(-491.98353909465015_dp)*Hr3(1,1,1) + &
                        nf*z**2*(468.27983539094663_dp)*Hr3(1,1,1) + nf*z*(2055.8134430727023_dp)*Hr3(1,1,1) + &
                        (6745.371742112483_dp)*Hr3(1,1,1) + z**3*(7841.051458757133_dp)*Hr3(1,1,1) + &
                        z**4*(10730.754458161866_dp)*Hr3(1,1,1) + z**3*(-6828.773662551441_dp)*Hr4(-1,-1,-1,0) + &
                        (-3902.1563786008237_dp)*Hr4(-1,-1,-1,0) + z**4*(3902.1563786008237_dp)*Hr4(-1,-1,-1,0) + &
                        z*(6828.773662551441_dp)*Hr4(-1,-1,-1,0) + z*(-6233.283950617284_dp)*Hr4(-1,-1,0,0) + &
                        z**4*(-3561.8765432098767_dp)*Hr4(-1,-1,0,0) + (3561.8765432098767_dp)*Hr4(-1,-1,0,0) + &
                        z**3*(6233.283950617284_dp)*Hr4(-1,-1,0,0) + z*(-5637.794238683127_dp)*Hr4(-1,-1,0,1) + &
                        z**4*(-3221.59670781893_dp)*Hr4(-1,-1,0,1) + (3221.59670781893_dp)*Hr4(-1,-1,0,1) + &
                        z**3*(5637.794238683127_dp)*Hr4(-1,-1,0,1) + z*(-3414.3868312757204_dp)*Hr4(-1,0,-1,0) + &
                        z**4*(-1951.0781893004119_dp)*Hr4(-1,0,-1,0) + (1951.0781893004119_dp)*Hr4(-1,0,-1,0) + &
                        z**3*(3414.3868312757204_dp)*Hr4(-1,0,-1,0) + z**3*(-2818.8971193415637_dp)*Hr4(-1,0,0,0) + &
                        (-1610.798353909465_dp)*Hr4(-1,0,0,0) + z**4*(1610.798353909465_dp)*Hr4(-1,0,0,0) + &
                        z*(2818.8971193415637_dp)*Hr4(-1,0,0,0) + z**3*(-5042.3045267489715_dp)*Hr4(-1,0,0,1) + &
                        (-2881.3168724279835_dp)*Hr4(-1,0,0,1) + z**4*(2881.3168724279835_dp)*Hr4(-1,0,0,1) + &
                        z*(5042.3045267489715_dp)*Hr4(-1,0,0,1) + z**3*(-2223.4074074074074_dp)*Hr4(-1,0,1,0) + &
                        (-1270.5185185185185_dp)*Hr4(-1,0,1,0) + z**4*(1270.5185185185185_dp)*Hr4(-1,0,1,0) + &
                        z*(2223.4074074074074_dp)*Hr4(-1,0,1,0) + z**3*(-4446.814814814815_dp)*Hr4(-1,0,1,1) + &
                        (-2541.037037037037_dp)*Hr4(-1,0,1,1) + z**4*(2541.037037037037_dp)*Hr4(-1,0,1,1) + &
                        z*(4446.814814814815_dp)*Hr4(-1,0,1,1) + z*(-13112.098765432098_dp)*Hr4(0,-1,-1,0) + &
                        z**3*(-6283.325102880658_dp)*Hr4(0,-1,-1,0) + nf*z**2*(-1269.4650205761318_dp)*Hr4(0,-1,-1,0) + &
                        z**4*(-1196.7736625514403_dp)*Hr4(0,-1,-1,0) + nf*z*(634.7325102880659_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**3*(634.7325102880659_dp)*Hr4(0,-1,-1,0) + (2705.382716049383_dp)*Hr4(0,-1,-1,0) + &
                        z**2*(17886.814814814814_dp)*Hr4(0,-1,-1,0) + z**2*(-16273.382716049382_dp)*Hr4(0,-1,0,0) + &
                        (-2462.0246913580245_dp)*Hr4(0,-1,0,0) + nf*z*(-578.1069958847737_dp)*Hr4(0,-1,0,0) + &
                        nf*z**3*(-578.1069958847737_dp)*Hr4(0,-1,0,0) + z**4*(1099.851851851852_dp)*Hr4(0,-1,0,0) + &
                        nf*z**2*(1156.2139917695474_dp)*Hr4(0,-1,0,0) + z**3*(5701.135802469135_dp)*Hr4(0,-1,0,0) + &
                        z*(11934.41975308642_dp)*Hr4(0,-1,0,0) + z**2*(-14659.950617283952_dp)*Hr4(0,-1,0,1) + &
                        (-2218.6666666666665_dp)*Hr4(0,-1,0,1) + nf*z*(-521.4814814814815_dp)*Hr4(0,-1,0,1) + &
                        nf*z**3*(-521.4814814814815_dp)*Hr4(0,-1,0,1) + z**4*(1002.9300411522634_dp)*Hr4(0,-1,0,1) + &
                        nf*z**2*(1042.962962962963_dp)*Hr4(0,-1,0,1) + z**3*(5118.9465020576135_dp)*Hr4(0,-1,0,1) + &
                        z*(10756.740740740743_dp)*Hr4(0,-1,0,1) + z**2*(-8943.407407407409_dp)*Hr4(0,0,-1,0) + &
                        z*(-3806.0246913580245_dp)*Hr4(0,0,-1,0) + (-1024._dp)*Hr4(0,0,-1,0) + &
                        nf*z**3*(-779.5884773662551_dp)*Hr4(0,0,-1,0) + nf*z*(144.8559670781893_dp)*Hr4(0,0,-1,0) + &
                        z**4*(269.69547325102883_dp)*Hr4(0,0,-1,0) + nf*z**2*(634.7325102880658_dp)*Hr4(0,0,-1,0) + &
                        z**3*(13503.736625514404_dp)*Hr4(0,0,-1,0) + z**2*(-13237.377229080934_dp)*Hr4(0,0,0,0) + &
                        z*(-4298.096021947875_dp)*Hr4(0,0,0,0) + nf*z**3*(-2056.427983539095_dp)*Hr4(0,0,0,0) + &
                        z**4*(-382.06858710562415_dp)*Hr4(0,0,0,0) + nf**2*z*(-31.604938271604937_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**3*(31.604938271604937_dp)*Hr4(0,0,0,0) + nf*z*(663.3525377229081_dp)*Hr4(0,0,0,0) + &
                        nf*z**2*(696.3621399176958_dp)*Hr4(0,0,0,0) + z**3*(27053.475994513035_dp)*Hr4(0,0,0,0) + &
                        z**2*(-29809.294924554186_dp)*Hr4(0,0,0,1) + nf*z**3*(-2750.1563786008232_dp)*Hr4(0,0,0,1) + &
                        z**4*(-1297.9094650205761_dp)*Hr4(0,0,0,1) + (-768._dp)*Hr4(0,0,0,1) + &
                        nf*z*(-547.028806584362_dp)*Hr4(0,0,0,1) + nf**2*z**3*(-6.320987654320987_dp)*Hr4(0,0,0,1) + &
                        nf**2*z*(6.320987654320987_dp)*Hr4(0,0,0,1) + nf*z**2*(1274.469135802469_dp)*Hr4(0,0,0,1) + &
                        z*(14079.956104252402_dp)*Hr4(0,0,0,1) + z**3*(44292.82853223594_dp)*Hr4(0,0,0,1) + &
                        z**2*(-11607.396433470502_dp)*Hr4(0,0,1,0) + nf*z**3*(-2863.3196159122085_dp)*Hr4(0,0,1,0) + &
                        (-2360.8888888888887_dp)*Hr4(0,0,1,0) + z**4*(-1942.650205761317_dp)*Hr4(0,0,1,0) + &
                        nf*z*(-147.1385459533608_dp)*Hr4(0,0,1,0) + nf*z**2*(425.8765432098766_dp)*Hr4(0,0,1,0) + &
                        z*(5594.513031550065_dp)*Hr4(0,0,1,0) + z**3*(44129.49245541838_dp)*Hr4(0,0,1,0) + &
                        nf*z**3*(-3186.304526748971_dp)*Hr4(0,0,1,1) + z**4*(-3177.349794238683_dp)*Hr4(0,0,1,1) + &
                        (-2958.222222222222_dp)*Hr4(0,0,1,1) + nf*z**2*(-810.6666666666669_dp)*Hr4(0,0,1,1) + &
                        z*(-349.7613168724274_dp)*Hr4(0,0,1,1) + nf*z*(86.3868312757204_dp)*Hr4(0,0,1,1) + &
                        z**2*(11600.329218106996_dp)*Hr4(0,0,1,1) + z**3*(46059.720164609054_dp)*Hr4(0,0,1,1) + &
                        (-2812.8395061728397_dp)*Hr4(0,1,0,0) + nf*z**3*(-2209.360768175583_dp)*Hr4(0,1,0,0) + &
                        z**4*(-2083.467764060357_dp)*Hr4(0,1,0,0) + z**2*(-1856.526748971195_dp)*Hr4(0,1,0,0) + &
                        nf*z*(-375.2208504801097_dp)*Hr4(0,1,0,0) + z*(7083.807956104254_dp)*Hr4(0,1,0,0) + &
                        z**3*(33482.09602194787_dp)*Hr4(0,1,0,0) + (-4326.716049382716_dp)*Hr4(0,1,0,1) + &
                        z**4*(-3810.5020576131687_dp)*Hr4(0,1,0,1) + nf*z**3*(-2598.716049382716_dp)*Hr4(0,1,0,1) + &
                        z**2*(-1448.559670781893_dp)*Hr4(0,1,0,1) + nf*z*(-1177.0205761316874_dp)*Hr4(0,1,0,1) + &
                        z*(19116.378600823045_dp)*Hr4(0,1,0,1) + z**3*(39823.670781893_dp)*Hr4(0,1,0,1) + &
                        (-4326.716049382716_dp)*Hr4(0,1,1,0) + z**4*(-3810.5020576131687_dp)*Hr4(0,1,1,0) + &
                        nf*z**3*(-2598.716049382716_dp)*Hr4(0,1,1,0) + z**2*(-1448.559670781893_dp)*Hr4(0,1,1,0) + &
                        nf*z*(-1177.0205761316874_dp)*Hr4(0,1,1,0) + z*(19116.378600823045_dp)*Hr4(0,1,1,0) + &
                        z**3*(39823.670781893_dp)*Hr4(0,1,1,0) + z**4*(-4890.3374485596705_dp)*Hr4(0,1,1,1) + &
                        (-3811.555555555555_dp)*Hr4(0,1,1,1) + nf*z*(-2625.843621399177_dp)*Hr4(0,1,1,1) + &
                        nf*z**3*(-1149.8930041152264_dp)*Hr4(0,1,1,1) + z**2*(53.72839506172817_dp)*Hr4(0,1,1,1) + &
                        z**3*(19048.0329218107_dp)*Hr4(0,1,1,1) + z*(38954.40329218107_dp)*Hr4(0,1,1,1) + &
                        z**2*(-2075.742112482852_dp)*Hr4(1,0,0,0) + (-1552.8559670781892_dp)*Hr4(1,0,0,0) + &
                        z**4*(-1552.8559670781892_dp)*Hr4(1,0,0,0) + nf*z*(-1011.358024691358_dp)*Hr4(1,0,0,0) + &
                        nf*z**3*(-1011.358024691358_dp)*Hr4(1,0,0,0) + z*(15839.517146776407_dp)*Hr4(1,0,0,0) + &
                        z**3*(15839.517146776407_dp)*Hr4(1,0,0,0) + (-3445.991769547325_dp)*Hr4(1,0,0,1) + &
                        z**4*(-3445.991769547325_dp)*Hr4(1,0,0,1) + z**2*(-3011.7750342935524_dp)*Hr4(1,0,0,1) + &
                        nf*z*(-1955.2921810699588_dp)*Hr4(1,0,0,1) + nf*z**3*(-1955.2921810699588_dp)*Hr4(1,0,0,1) + &
                        z*(30539.237311385463_dp)*Hr4(1,0,0,1) + z**3*(30539.237311385463_dp)*Hr4(1,0,0,1) + &
                        (-3786.2716049382716_dp)*Hr4(1,0,1,0) + z**4*(-3786.2716049382716_dp)*Hr4(1,0,1,0) + &
                        nf*z*(-1887.868312757202_dp)*Hr4(1,0,1,0) + nf*z**3*(-1887.8683127572015_dp)*Hr4(1,0,1,0) + &
                        z**2*(-1872.0658436214003_dp)*Hr4(1,0,1,0) + z*(29399.440329218105_dp)*Hr4(1,0,1,0) + &
                        z**3*(29399.44032921811_dp)*Hr4(1,0,1,0) + (-3786.2716049382716_dp)*Hr4(1,0,1,1) + &
                        z**4*(-3786.2716049382716_dp)*Hr4(1,0,1,1) + nf*z*(-1887.868312757202_dp)*Hr4(1,0,1,1) + &
                        nf*z**3*(-1887.8683127572015_dp)*Hr4(1,0,1,1) + z**2*(-793.2839506172816_dp)*Hr4(1,0,1,1) + &
                        z**3*(28860.049382716046_dp)*Hr4(1,0,1,1) + z*(28860.04938271605_dp)*Hr4(1,0,1,1) + &
                        (-3445.991769547325_dp)*Hr4(1,1,0,0) + z**4*(-3445.991769547325_dp)*Hr4(1,1,0,0) + &
                        z**2*(-3011.7750342935524_dp)*Hr4(1,1,0,0) + nf*z*(-1955.2921810699588_dp)*Hr4(1,1,0,0) + &
                        nf*z**3*(-1955.2921810699588_dp)*Hr4(1,1,0,0) + z*(30539.237311385463_dp)*Hr4(1,1,0,0) + &
                        z**3*(30539.237311385463_dp)*Hr4(1,1,0,0) + (-3786.2716049382716_dp)*Hr4(1,1,0,1) + &
                        z**4*(-3786.2716049382716_dp)*Hr4(1,1,0,1) + nf*z*(-1887.868312757202_dp)*Hr4(1,1,0,1) + &
                        nf*z**3*(-1887.8683127572015_dp)*Hr4(1,1,0,1) + z**2*(-793.2839506172816_dp)*Hr4(1,1,0,1) + &
                        z**3*(28860.049382716046_dp)*Hr4(1,1,0,1) + z*(28860.04938271605_dp)*Hr4(1,1,0,1) + &
                        (-3786.2716049382716_dp)*Hr4(1,1,1,0) + z**4*(-3786.2716049382716_dp)*Hr4(1,1,1,0) + &
                        nf*z*(-1887.868312757202_dp)*Hr4(1,1,1,0) + nf*z**3*(-1887.8683127572015_dp)*Hr4(1,1,1,0) + &
                        z**2*(-793.2839506172816_dp)*Hr4(1,1,1,0) + z**3*(28860.049382716046_dp)*Hr4(1,1,1,0) + &
                        z*(28860.04938271605_dp)*Hr4(1,1,1,0) + z*(-2157.5637860082306_dp)*Hr4(1,1,1,1) + &
                        z**3*(-2157.5637860082306_dp)*Hr4(1,1,1,1) + z**2*(4315.127572016461_dp)*Hr4(1,1,1,1) + &
                        z**2*(-11706.469135802468_dp)*Hr5(0,-1,-1,-1,0) + z*(5853.234567901234_dp)*Hr5(0,-1,-1,-1,0) + &
                        z**3*(5853.234567901234_dp)*Hr5(0,-1,-1,-1,0) + z*(-5342.814814814815_dp)*Hr5(0,-1,-1,0,0) + &
                        z**3*(-5342.814814814815_dp)*Hr5(0,-1,-1,0,0) + z**2*(10685.62962962963_dp)*Hr5(0,-1,-1,0,0) + &
                        z*(-4832.395061728395_dp)*Hr5(0,-1,-1,0,1) + z**3*(-4832.395061728395_dp)*Hr5(0,-1,-1,0,1) + &
                        z**2*(9664.79012345679_dp)*Hr5(0,-1,-1,0,1) + z*(-2926.617283950617_dp)*Hr5(0,-1,0,-1,0) + &
                        z**3*(-2926.617283950617_dp)*Hr5(0,-1,0,-1,0) + z**2*(5853.234567901234_dp)*Hr5(0,-1,0,-1,0) + &
                        z**2*(-4832.395061728395_dp)*Hr5(0,-1,0,0,0) + z*(2416.1975308641977_dp)*Hr5(0,-1,0,0,0) + &
                        z**3*(2416.1975308641977_dp)*Hr5(0,-1,0,0,0) + z**2*(-8643.95061728395_dp)*Hr5(0,-1,0,0,1) + &
                        z*(4321.975308641975_dp)*Hr5(0,-1,0,0,1) + z**3*(4321.975308641975_dp)*Hr5(0,-1,0,0,1) + &
                        z**2*(-3811.555555555555_dp)*Hr5(0,-1,0,1,0) + z*(1905.7777777777776_dp)*Hr5(0,-1,0,1,0) + &
                        z**3*(1905.7777777777776_dp)*Hr5(0,-1,0,1,0) + z**2*(-7623.11111111111_dp)*Hr5(0,-1,0,1,1) + &
                        z*(3811.555555555555_dp)*Hr5(0,-1,0,1,1) + z**3*(3811.555555555555_dp)*Hr5(0,-1,0,1,1) + &
                        z**3*(-9013.728395061727_dp)*Hr5(0,0,-1,-1,0) + z*(3160.4938271604938_dp)*Hr5(0,0,-1,-1,0) + &
                        z**2*(5853.234567901235_dp)*Hr5(0,0,-1,-1,0) + z**2*(-5342.814814814815_dp)*Hr5(0,0,-1,0,0) + &
                        z*(-2868.1481481481483_dp)*Hr5(0,0,-1,0,0) + z**3*(8210.962962962964_dp)*Hr5(0,0,-1,0,0) + &
                        z**2*(-4832.395061728395_dp)*Hr5(0,0,-1,0,1) + z*(-2575.8024691358023_dp)*Hr5(0,0,-1,0,1) + &
                        z**3*(7408.197530864198_dp)*Hr5(0,0,-1,0,1) + z**2*(-8375.308641975309_dp)*Hr5(0,0,0,-1,0) + &
                        nf*z*(-113.77777777777777_dp)*Hr5(0,0,0,-1,0) + nf*z**3*(-113.77777777777777_dp)*Hr5(0,0,0,-1,0) + &
                        nf*z**2*(227.55555555555554_dp)*Hr5(0,0,0,-1,0) + z*(1144.0987654320988_dp)*Hr5(0,0,0,-1,0) + &
                        z**3*(7231.20987654321_dp)*Hr5(0,0,0,-1,0) + z**2*(-9323.456790123457_dp)*Hr5(0,0,0,0,0) + &
                        nf*z**3*(-417.18518518518516_dp)*Hr5(0,0,0,0,0) + nf*z*(-151.7037037037037_dp)*Hr5(0,0,0,0,0) + &
                        nf*z**2*(568.8888888888889_dp)*Hr5(0,0,0,0,0) + z*(1041.5582990397804_dp)*Hr5(0,0,0,0,0) + &
                        z**3*(9090.984910836763_dp)*Hr5(0,0,0,0,0) + z**2*(-11891.358024691363_dp)*Hr5(0,0,0,0,1) + &
                        z*(-3153.997256515775_dp)*Hr5(0,0,0,0,1) + nf*z**3*(-297.08641975308643_dp)*Hr5(0,0,0,0,1) + &
                        nf*z*(297.08641975308643_dp)*Hr5(0,0,0,0,1) + z**3*(18191.802469135804_dp)*Hr5(0,0,0,0,1) + &
                        z**2*(-17393.777777777777_dp)*Hr5(0,0,0,1,0) + nf*z**3*(-82.17283950617283_dp)*Hr5(0,0,0,1,0) + &
                        nf*z*(82.17283950617283_dp)*Hr5(0,0,0,1,0) + z*(2928.8998628257887_dp)*Hr5(0,0,0,1,0) + &
                        z**3*(19768.88888888889_dp)*Hr5(0,0,0,1,0) + z**2*(-23893.333333333336_dp)*Hr5(0,0,0,1,1) + &
                        z*(7692.641975308642_dp)*Hr5(0,0,0,1,1) + z**3*(25370.337448559672_dp)*Hr5(0,0,0,1,1) + &
                        z**2*(-9799.111111111111_dp)*Hr5(0,0,1,0,0) + z*(1938.7873799725648_dp)*Hr5(0,0,1,0,0) + &
                        z**3*(13793.62414266118_dp)*Hr5(0,0,1,0,0) + z**2*(-10894.222222222223_dp)*Hr5(0,0,1,0,1) + &
                        z*(725.6844993141294_dp)*Hr5(0,0,1,0,1) + z**3*(21945.240054869686_dp)*Hr5(0,0,1,0,1) + &
                        z**2*(-10894.222222222223_dp)*Hr5(0,0,1,1,0) + z*(725.6844993141294_dp)*Hr5(0,0,1,1,0) + &
                        z**3*(21945.240054869686_dp)*Hr5(0,0,1,1,0) + z*(-4482.633744855968_dp)*Hr5(0,0,1,1,1) + &
                        z**3*(22282.53497942387_dp)*Hr5(0,0,1,1,1) + z*(-823.4842249657067_dp)*Hr5(0,1,0,0,0) + &
                        z**3*(6127.495198902607_dp)*Hr5(0,1,0,0,0) + z*(-1505.6241426611796_dp)*Hr5(0,1,0,0,1) + &
                        z**3*(13282.32647462277_dp)*Hr5(0,1,0,0,1) + z*(-1364.2798353909475_dp)*Hr5(0,1,0,1,0) + &
                        z**3*(14309.66255144033_dp)*Hr5(0,1,0,1,0) + z*(793.2839506172859_dp)*Hr5(0,1,0,1,1) + &
                        z**3*(16467.22633744856_dp)*Hr5(0,1,0,1,1) + z*(-1505.6241426611796_dp)*Hr5(0,1,1,0,0) + &
                        z**3*(13282.32647462277_dp)*Hr5(0,1,1,0,0) + z*(793.2839506172859_dp)*Hr5(0,1,1,0,1) + &
                        z**3*(16467.22633744856_dp)*Hr5(0,1,1,0,1) + z*(793.2839506172859_dp)*Hr5(0,1,1,1,0) + &
                        z**3*(16467.22633744856_dp)*Hr5(0,1,1,1,0) + z*(8630.255144032923_dp)*Hr5(0,1,1,1,1) + &
                        z**3*(8630.255144032923_dp)*Hr5(0,1,1,1,1) + z*(1573.2235939643347_dp)*Hr5(1,0,0,0,0) + &
                        z**3*(1573.2235939643347_dp)*Hr5(1,0,0,0,0) + z*(4584.82304526749_dp)*Hr5(1,0,0,0,1) + &
                        z**3*(4584.82304526749_dp)*Hr5(1,0,0,0,1) + z*(5888.3511659807955_dp)*Hr5(1,0,0,1,0) + &
                        z**3*(5888.3511659807955_dp)*Hr5(1,0,0,1,0) + z*(8899.95061728395_dp)*Hr5(1,0,0,1,1) + &
                        z**3*(8899.95061728395_dp)*Hr5(1,0,0,1,1) + z*(5888.3511659807955_dp)*Hr5(1,0,1,0,0) + &
                        z**3*(5888.3511659807955_dp)*Hr5(1,0,1,0,0) + z*(8630.255144032923_dp)*Hr5(1,0,1,0,1) + &
                        z**3*(8630.255144032923_dp)*Hr5(1,0,1,0,1) + z*(8630.255144032923_dp)*Hr5(1,0,1,1,0) + &
                        z**3*(8630.255144032923_dp)*Hr5(1,0,1,1,0) + z*(8630.255144032923_dp)*Hr5(1,0,1,1,1) + &
                        z**3*(8630.255144032923_dp)*Hr5(1,0,1,1,1) + z*(4584.82304526749_dp)*Hr5(1,1,0,0,0) + &
                        z**3*(4584.82304526749_dp)*Hr5(1,1,0,0,0) + z*(8899.95061728395_dp)*Hr5(1,1,0,0,1) + &
                        z**3*(8899.95061728395_dp)*Hr5(1,1,0,0,1) + z*(8630.255144032923_dp)*Hr5(1,1,0,1,0) + &
                        z**3*(8630.255144032923_dp)*Hr5(1,1,0,1,0) + z*(8630.255144032923_dp)*Hr5(1,1,0,1,1) + &
                        z**3*(8630.255144032923_dp)*Hr5(1,1,0,1,1) + z*(8899.95061728395_dp)*Hr5(1,1,1,0,0) + &
                        z**3*(8899.95061728395_dp)*Hr5(1,1,1,0,0) + z*(8630.255144032923_dp)*Hr5(1,1,1,0,1) + &
                        z**3*(8630.255144032923_dp)*Hr5(1,1,1,0,1) + z*(8630.255144032923_dp)*Hr5(1,1,1,1,0) + &
                        z**3*(8630.255144032923_dp)*Hr5(1,1,1,1,0))/(z*((-1._dp) + z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-214656.29903213575_dp) + nf**2*(-5048.973537330817_dp) + nf**4*(-2.2978555612816147_dp) &
                            + nf**3*(182.16604600021253_dp) + nf*(56656.23246393325_dp)

                        else
                            Ibar_select = (-214656.29903213575_dp) + nf**2*(-5048.973537330817_dp) + nf**4*(-2.2978555612816147_dp) &
                            + nf**3*(182.16604600021253_dp) + nf*(56656.23246393325_dp)

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-11040.464477732992_dp) + nf**2*(-688.5635930976614_dp) + &
                            nf**4*(-1.755829903978052_dp) + nf**3*(85.25362055837081_dp) + (107585.57577873324_dp)

                        else
                            Ibar_select = nf*(-11040.464477732992_dp) + nf**2*(-688.5635930976614_dp) + &
                            nf**4*(-1.755829903978052_dp) + nf**3*(85.25362055837081_dp) + (107585.57577873324_dp)

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-23815.97447535026_dp) + nf**3*(-51.50434385002286_dp) + nf**2*(2141.152664766941_dp) &
                            + (69609.18235654884_dp)

                        else
                            Ibar_select = nf*(-23815.97447535026_dp) + nf**3*(-51.50434385002286_dp) + nf**2*(2141.152664766941_dp) &
                            + (69609.18235654884_dp)

                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-53513.496138672184_dp) + nf**2*(-337.119341563786_dp) + nf*(9019.271227709527_dp)
                        else
                            Ibar_select = (-53513.496138672184_dp) + nf**2*(-337.119341563786_dp) + nf*(9019.271227709527_dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-599.3232738911752_dp) + (6427.068364624576_dp)
                        else
                            Ibar_select = nf*(-599.3232738911752_dp) + (6427.068364624576_dp)
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
                        Ibar_select = (z**2*(-23359.431448872096_dp) + z**4*(-10847.762647276684_dp) + (-10166.453905323033_dp) + &
                        nf*z*(-5314.9142268490305_dp) + z**3*(-2596.555453024232_dp) + nf**2*z**2*(-139.88696844993143_dp) + &
                        nf**2*(-24.97960676726109_dp) + nf**2*z**4*(-24.97960676726109_dp) + nf**3*z**3*(-6.127846364883402_dp) + &
                        nf**4*z*(-0.10534979423868314_dp) + nf**4*z**3*(0.10534979423868314_dp) + nf**3*(0.21069958847736628_dp) + &
                        nf**3*z**4*(0.21069958847736628_dp) + nf**3*z*(2.5108367626886148_dp) + nf**3*z**2*(3.195610425240055_dp) + &
                        nf**2*z**3*(64.42130104552423_dp) + nf**2*z*(125.42488093892938_dp) + nf*z**3*(249.97017966148633_dp) + &
                        nf*(909.4256397571147_dp) + nf*z**4*(952.74900612678_dp) + nf*z**2*(2986.643197722121_dp) + &
                        z*(50374.403330723864_dp) + z**4*(-4610.9553510738615_dp)*Hr1(0) + z**2*(-4542.7253670158225_dp)*Hr1(0) + &
                        (-3750.342016911999_dp)*Hr1(0) + nf*z**3*(-2026.7521226898195_dp)*Hr1(0) + &
                        nf*z*(-290.1689848187169_dp)*Hr1(0) + nf**2*z**2*(-24.625514403292183_dp)*Hr1(0) + &
                        nf**2*(-3.9681755829903973_dp)*Hr1(0) + nf**3*z**3*(-2.949794238683128_dp)*Hr1(0) + &
                        nf**2*z**4*(-2.5283950617283955_dp)*Hr1(0) + nf**3*z*(-0.5618655692729768_dp)*Hr1(0) + &
                        nf**2*z*(14.819204389574757_dp)*Hr1(0) + nf**2*z**3*(145.8831275720165_dp)*Hr1(0) + &
                        nf*z**4*(215.73882030178328_dp)*Hr1(0) + nf*(249.73168724279836_dp)*Hr1(0) + &
                        nf*z**2*(587.5010773438325_dp)*Hr1(0) + z*(3408.133680559081_dp)*Hr1(0) + &
                        z**3*(11367.081406232726_dp)*Hr1(0) + z**2*(-26849.29405010864_dp)*Hr1(1) + (-8111.754777696587_dp)*Hr1(1) + &
                        z**4*(-8111.754777696587_dp)*Hr1(1) + nf*z**3*(-2551.1561083471934_dp)*Hr1(1) + &
                        nf*z*(-319.0112441905299_dp)*Hr1(1) + nf**2*z*(-98.45816186556927_dp)*Hr1(1) + &
                        nf**2*z**2*(-49.251028806584365_dp)*Hr1(1) + nf**2*(-6.496570644718793_dp)*Hr1(1) + &
                        nf**2*z**4*(-6.496570644718793_dp)*Hr1(1) + nf**3*z**3*(-3.5116598079561046_dp)*Hr1(1) + &
                        nf**3*z*(3.5116598079561046_dp)*Hr1(1) + nf**2*z**3*(160.70233196159123_dp)*Hr1(1) + &
                        nf*(465.4705075445816_dp)*Hr1(1) + nf*z**4*(465.4705075445816_dp)*Hr1(1) + &
                        nf*z**2*(1939.2263374485601_dp)*Hr1(1) + z**3*(21195.127228645328_dp)*Hr1(1) + &
                        z*(22742.181191182597_dp)*Hr1(1) + z**2*(-6156.0010325896155_dp)*Hr2(0,0) + &
                        nf*z**3*(-959.7975740190537_dp)*Hr2(0,0) + (-704._dp)*Hr2(0,0) + z**4*(-538.8641975308642_dp)*Hr2(0,0) + &
                        nf*z*(-202.6935096572152_dp)*Hr2(0,0) + nf**2*z**2*(-7.743209876543209_dp)*Hr2(0,0) + &
                        nf**2*z*(1.4573388203017839_dp)*Hr2(0,0) + nf*z**4*(11.799176954732513_dp)*Hr2(0,0) + &
                        nf*(24.9679012345679_dp)*Hr2(0,0) + nf**2*z**3*(25.95116598079561_dp)*Hr2(0,0) + &
                        nf*z**2*(611.6170096021947_dp)*Hr2(0,0) + z*(2764.515191070828_dp)*Hr2(0,0) + &
                        z**3*(7669.756351793233_dp)*Hr2(0,0) + z**2*(-6200.493827160493_dp)*Hr2(0,1) + &
                        (-2085.925925925926_dp)*Hr2(0,1) + z**4*(-1755.6543209876543_dp)*Hr2(0,1) + &
                        nf*z**3*(-1739.3251028806587_dp)*Hr2(0,1) + nf*z*(-100.34567901234561_dp)*Hr2(0,1) + &
                        nf**2*z*(5.004115226337454_dp)*Hr2(0,1) + nf*z**4*(42.139917695473244_dp)*Hr2(0,1) + &
                        nf**2*z**3*(53.99176954732511_dp)*Hr2(0,1) + nf*(68.47736625514405_dp)*Hr2(0,1) + &
                        nf*z**2*(186.73251028806575_dp)*Hr2(0,1) + z*(2736.1055604580056_dp)*Hr2(0,1) + &
                        z**3*(15524.924908596062_dp)*Hr2(0,1) + z**2*(-6200.493827160494_dp)*Hr2(1,0) + &
                        (-1920.7901234567903_dp)*Hr2(1,0) + z**4*(-1920.7901234567903_dp)*Hr2(1,0) + &
                        nf*z*(-919.8353909465021_dp)*Hr2(1,0) + nf*z**3*(-919.8353909465021_dp)*Hr2(1,0) + &
                        nf**2*z**3*(29.49794238683128_dp)*Hr2(1,0) + nf**2*z*(29.497942386831287_dp)*Hr2(1,0) + &
                        nf*(55.30864197530864_dp)*Hr2(1,0) + nf*z**4*(55.30864197530864_dp)*Hr2(1,0) + &
                        nf*z**2*(186.73251028806584_dp)*Hr2(1,0) + z*(9426.269415610615_dp)*Hr2(1,0) + &
                        z**3*(9426.269415610615_dp)*Hr2(1,0) + (z**2*(-12400.987654320988_dp) + (-3841.5802469135806_dp) + &
                        z**4*(-3841.5802469135806_dp) + nf**2*(z*(-58.99588477366256_dp) + z**3*(58.99588477366256_dp)) + &
                        nf*(z**3*(-1839.6707818930042_dp) + (110.61728395061728_dp) + z**4*(110.61728395061728_dp) + &
                        z**2*(373.4650205761317_dp) + z*(1244.971193415638_dp)) + z*(1823.1176790940806_dp) + &
                        z**3*(18261.030469054065_dp))*Hr2(1,1) + ((-51.2_dp) + z**4*(-17.979698216735255_dp) + &
                        z*((-292.2842249657066_dp) + nf**2*(-2.5283950617283955_dp) + nf*(55.50178326474622_dp)) + &
                        z**2*((-2318.705075445817_dp) + nf*(116.14814814814815_dp)) + z**3*(nf*(-216.5991769547325_dp) + &
                        nf**2*(2.5283950617283955_dp) + (3286.9838134430734_dp)))*Hr3(0,0,0) + &
                        z**2*(-4708.433470507545_dp)*Hr3(0,0,1) + nf*z**3*(-356.4334705075446_dp)*Hr3(0,0,1) + &
                        (-199.11111111111114_dp)*Hr3(0,0,1) + z**4*(-84.27983539094649_dp)*Hr3(0,0,1) + &
                        nf*z*(-55.659807956104274_dp)*Hr3(0,0,1) + nf*z**2*(232.29629629629633_dp)*Hr3(0,0,1) + &
                        z*(1516.1591220850476_dp)*Hr3(0,0,1) + z**3*(5902.924554183813_dp)*Hr3(0,0,1) + &
                        nf*z**3*(-285.23456790123464_dp)*Hr3(0,1,0) + (-274.962962962963_dp)*Hr3(0,1,0) + &
                        z**2*(-213.0699588477363_dp)*Hr3(0,1,0) + z**4*(-160.13168724279836_dp)*Hr3(0,1,0) + &
                        z*(-48.855967078189224_dp)*Hr3(0,1,0) + nf*z*(15.539094650205758_dp)*Hr3(0,1,0) + &
                        z**3*(4337.909465020576_dp)*Hr3(0,1,0) + nf*z**3*(-570.4691358024693_dp)*Hr3(0,1,1) + &
                        (-549.925925925926_dp)*Hr3(0,1,1) + z**2*(-426.1399176954726_dp)*Hr3(0,1,1) + &
                        z**4*(-320.2633744855967_dp)*Hr3(0,1,1) + z*(-97.71193415637845_dp)*Hr3(0,1,1) + &
                        nf*z*(31.078189300411516_dp)*Hr3(0,1,1) + z**3*(8675.818930041152_dp)*Hr3(0,1,1) + &
                        (-145.0315500685871_dp)*Hr3(1,0,0) + z**4*(-145.0315500685871_dp)*Hr3(1,0,0) + &
                        z**2*(-142.04663923182446_dp)*Hr3(1,0,0) + nf*z*(-89.8984910836763_dp)*Hr3(1,0,0) + &
                        nf*z**3*(-89.89849108367626_dp)*Hr3(1,0,0) + z*(1429.684499314129_dp)*Hr3(1,0,0) + &
                        z**3*(1429.684499314129_dp)*Hr3(1,0,0) + (-435.0946502057613_dp)*Hr3(1,0,1) + &
                        z**4*(-435.0946502057613_dp)*Hr3(1,0,1) + z**2*(-426.1399176954737_dp)*Hr3(1,0,1) + &
                        nf*z*(-269.6954732510289_dp)*Hr3(1,0,1) + nf*z**3*(-269.6954732510289_dp)*Hr3(1,0,1) + &
                        z*(4289.0534979423865_dp)*Hr3(1,0,1) + z**3*(4289.0534979423865_dp)*Hr3(1,0,1) + &
                        (-435.0946502057613_dp)*Hr3(1,1,0) + z**4*(-435.0946502057613_dp)*Hr3(1,1,0) + &
                        z**2*(-426.1399176954737_dp)*Hr3(1,1,0) + nf*z*(-269.6954732510289_dp)*Hr3(1,1,0) + &
                        nf*z**3*(-269.6954732510289_dp)*Hr3(1,1,0) + z*(4289.0534979423865_dp)*Hr3(1,1,0) + &
                        z**3*(4289.0534979423865_dp)*Hr3(1,1,0) + z*(-5985.448559670782_dp)*Hr3(1,1,1) + &
                        (-870.1893004115226_dp)*Hr3(1,1,1) + z**4*(-870.1893004115226_dp)*Hr3(1,1,1) + &
                        z**2*(-852.2798353909474_dp)*Hr3(1,1,1) + nf*z**3*(-539.3909465020578_dp)*Hr3(1,1,1) + &
                        nf*z*(539.3909465020578_dp)*Hr3(1,1,1) + z**3*(8578.106995884773_dp)*Hr3(1,1,1) + &
                        z**2*(-827.7333333333335_dp)*Hr4(0,0,0,0) + nf*z**3*(-26.548148148148147_dp)*Hr4(0,0,0,0) + &
                        nf*z*(-7.585185185185186_dp)*Hr4(0,0,0,0) + nf*z**2*(34.133333333333326_dp)*Hr4(0,0,0,0) + &
                        z*(51.79698216735254_dp)*Hr4(0,0,0,0) + z**3*(811.8957475994513_dp)*Hr4(0,0,0,0) + &
                        z**2*(-1237.3333333333335_dp)*Hr4(0,0,0,1) + z*(-205.2565157750343_dp)*Hr4(0,0,0,1) + &
                        nf*z**3*(-18.962962962962962_dp)*Hr4(0,0,0,1) + nf*z*(18.962962962962962_dp)*Hr4(0,0,0,1) + &
                        z**3*(1622.3868312757202_dp)*Hr4(0,0,0,1) + z**2*(-1237.3333333333335_dp)*Hr4(0,0,1,0) + &
                        z*(337.29492455418386_dp)*Hr4(0,0,1,0) + z**3*(1259.6323731138546_dp)*Hr4(0,0,1,0) + &
                        z**2*(-2474.666666666667_dp)*Hr4(0,0,1,1) + z*(674.5898491083677_dp)*Hr4(0,0,1,1) + &
                        z**3*(2519.264746227709_dp)*Hr4(0,0,1,1) + z*(-127.64883401920434_dp)*Hr4(0,1,0,0) + &
                        z**3*(487.24279835390945_dp)*Hr4(0,1,0,0) + z*(-382.9465020576131_dp)*Hr4(0,1,0,1) + &
                        z**3*(1461.7283950617286_dp)*Hr4(0,1,0,1) + z*(-382.9465020576131_dp)*Hr4(0,1,1,0) + &
                        z**3*(1461.7283950617286_dp)*Hr4(0,1,1,0) + z*(-765.8930041152262_dp)*Hr4(0,1,1,1) + &
                        z**3*(2923.4567901234573_dp)*Hr4(0,1,1,1) + z**3*(89.89849108367626_dp)*Hr4(1,0,0,0) + &
                        z*(89.8984910836763_dp)*Hr4(1,0,0,0) + z**3*(359.59396433470505_dp)*Hr4(1,0,0,1) + &
                        z*(359.5939643347052_dp)*Hr4(1,0,0,1) + z*(539.3909465020578_dp)*Hr4(1,0,1,0) + &
                        z**3*(539.3909465020578_dp)*Hr4(1,0,1,0) + z*(1078.7818930041155_dp)*Hr4(1,0,1,1) + &
                        z**3*(1078.7818930041155_dp)*Hr4(1,0,1,1) + z**3*(359.59396433470505_dp)*Hr4(1,1,0,0) + &
                        z*(359.5939643347052_dp)*Hr4(1,1,0,0) + z*(1078.7818930041155_dp)*Hr4(1,1,0,1) + &
                        z**3*(1078.7818930041155_dp)*Hr4(1,1,0,1) + z*(1078.7818930041155_dp)*Hr4(1,1,1,0) + &
                        z**3*(1078.7818930041155_dp)*Hr4(1,1,1,0) + z*(-2157.563786008231_dp)*Hr4(1,1,1,1) + &
                        z**3*(2157.563786008231_dp)*Hr4(1,1,1,1))/(z*((-1._dp) + z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-21369.375091003178_dp) + nf**2*(-429.0804214240642_dp) + nf**4*(-0.1580246913580247_dp) &
                            + nf**3*(14.230769836943017_dp) + nf*(5214.295952850753_dp)

                        else
                            Ibar_select = (-21369.375091003178_dp) + nf**2*(-429.0804214240642_dp) + nf**4*(-0.1580246913580247_dp) &
                            + nf**3*(14.230769836943017_dp) + nf*(5214.295952850753_dp)

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-2466.27614917806_dp) + nf**2*(-19.577881554172123_dp) + &
                            nf**4*(-0.21069958847736625_dp) + nf**3*(8.638683127572017_dp) + (22187.28351083629_dp)

                        else
                            Ibar_select = nf*(-2466.27614917806_dp) + nf**2*(-19.577881554172123_dp) + &
                            nf**4*(-0.21069958847736625_dp) + nf**3*(8.638683127572017_dp) + (22187.28351083629_dp)

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-2232.144864156664_dp) + (-1547.0539625372685_dp) + nf**3*(-7.023319615912208_dp) + &
                            nf**2*(259.1604938271605_dp)

                        else
                            Ibar_select = nf*(-2232.144864156664_dp) + (-1547.0539625372685_dp) + nf**3*(-7.023319615912208_dp) + &
                            nf**2*(259.1604938271605_dp)

                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-8218.956394979994_dp) + nf**2*(-58.99588477366255_dp) + nf*(1542.320987654321_dp)
                        else
                            Ibar_select = (-8218.956394979994_dp) + nf**2*(-58.99588477366255_dp) + nf*(1542.320987654321_dp)
                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-179.79698216735252_dp) + (2427.259259259259_dp)
                        else
                            Ibar_select = nf*(-179.79698216735252_dp) + (2427.259259259259_dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (-179.79698216735252_dp)
                        else
                            Ibar_select = (-179.79698216735252_dp)
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
                        Ibar_select = (z**2*(-282672.1442567002_dp) + z**4*(-150595.72949352118_dp) + (-148537.61617422395_dp) + &
                        nf*z*(-100540.0629007812_dp) + z**3*(-25117.96537503645_dp) + nf**2*z**2*(-3075.2965788097654_dp) + &
                        nf**2*z**4*(-936.3362909285732_dp) + nf**2*(-894.0526853517797_dp) + nf**2*z**3*(-129.47180888895036_dp) + &
                        nf**3*z**3*(-75.10411948270566_dp) + nf**3*z*(-68.55528203129116_dp) + nf**4*z*(-2.126505105928974_dp) + &
                        nf**4*z**2*(-1.9626276482243563_dp) + nf**4*(-0.11705532693187014_dp) + nf**4*z**4*(-0.11705532693187014_dp) &
                        + nf**5*z**3*(-0.05852766346593507_dp) + nf**5*z*(0.05852766346593507_dp) + &
                        nf**4*z**3*(4.3232434080170705_dp) + nf**3*(17.225341665396535_dp) + nf**3*z**4*(17.225341665396535_dp) + &
                        nf**3*z**2*(109.20871818320377_dp) + nf**2*z*(5239.276556250513_dp) + nf*z**3*(6198.40542802125_dp) + &
                        nf*(19416.34072939799_dp) + nf*z**4*(21022.417179890606_dp) + nf*z**2*(47436.770356626024_dp) + &
                        z*(656279.6044826803_dp) + z**4*(-93646.69234958958_dp)*Hr1(0) + (-61963.74924070922_dp)*Hr1(0) + &
                        z**2*(-47327.79624386877_dp)*Hr1(0) + nf*z**3*(-25517.853842391614_dp)*Hr1(0) + &
                        nf*z*(-5652.988641757365_dp)*Hr1(0) + nf**2*z**2*(-663.7042384176309_dp)*Hr1(0) + &
                        nf**2*(-231.72662703856113_dp)*Hr1(0) + nf**2*z**4*(-179.62139917695472_dp)*Hr1(0) + &
                        nf**3*z**3*(-118.65313214449016_dp)*Hr1(0) + nf**3*z*(-12.956073769242494_dp)*Hr1(0) + &
                        nf**4*z*(0.3589696692577351_dp)*Hr1(0) + nf**3*z**4*(1.4983081847279376_dp)*Hr1(0) + &
                        nf**4*z**3*(1.7792409693644262_dp)*Hr1(0) + nf**3*(2.5908245694253926_dp)*Hr1(0) + &
                        nf**3*z**2*(18.003109282121624_dp)*Hr1(0) + nf**2*z*(306.6036291137005_dp)*Hr1(0) + &
                        nf**2*z**3*(2560.723596183306_dp)*Hr1(0) + nf*(6684.803817022826_dp)*Hr1(0) + &
                        nf*z**4*(7280.491409874511_dp)*Hr1(0) + nf*z**2*(8222.294882233005_dp)*Hr1(0) + &
                        z*(65864.89524387354_dp)*Hr1(0) + z**3*(131161.168680524_dp)*Hr1(0) + z**2*(-392164.87008948444_dp)*Hr1(1) + &
                        (-144109.62722332927_dp)*Hr1(1) + z**4*(-144109.62722332927_dp)*Hr1(1) + &
                        nf*z**3*(-41187.941875229924_dp)*Hr1(1) + nf*z*(-33121.60023269498_dp)*Hr1(1) + &
                        nf**2*z**2*(-1827.1868312757201_dp)*Hr1(1) + nf**2*(-411.348026215516_dp)*Hr1(1) + &
                        nf**2*z**4*(-411.34802621551586_dp)*Hr1(1) + nf**2*z*(-327.67182218768244_dp)*Hr1(1) + &
                        nf**3*z**3*(-131.60920591373264_dp)*Hr1(1) + nf**4*z*(-2.138210638622161_dp)*Hr1(1) + &
                        nf**4*z**3*(2.138210638622161_dp)*Hr1(1) + nf**3*(4.08913275415333_dp)*Hr1(1) + &
                        nf**3*z**4*(4.08913275415333_dp)*Hr1(1) + nf**3*z**2*(36.006218564243255_dp)*Hr1(1) + &
                        nf**3*z*(87.42472184118273_dp)*Hr1(1) + nf**2*z**3*(2977.554705894434_dp)*Hr1(1) + &
                        nf*(13651.672564039956_dp)*Hr1(1) + nf*z**4*(13651.672564039956_dp)*Hr1(1) + &
                        nf*z**2*(45565.35562263481_dp)*Hr1(1) + z**3*(303476.58858085563_dp)*Hr1(1) + &
                        z*(399818.0451072959_dp)*Hr1(1) + z**2*(-82753.24905042429_dp)*Hr2(0,0) + &
                        z**4*(-16020.39216149099_dp)*Hr2(0,0) + nf*z**3*(-15988.980214393174_dp)*Hr2(0,0) + &
                        (-15266.82056315651_dp)*Hr2(0,0) + nf*z*(-5407.845661656312_dp)*Hr2(0,0) + &
                        nf**2*z**2*(-555.8352690138698_dp)*Hr2(0,0) + nf**2*(-21.53037646700198_dp)*Hr2(0,0) + &
                        nf**3*z**3*(-17.113488797439416_dp)*Hr2(0,0) + nf**2*z**4*(-8.67770156988264_dp)*Hr2(0,0) + &
                        nf**3*z*(-0.44481024234110755_dp)*Hr2(0,0) + nf**3*z**2*(3.5116598079561046_dp)*Hr2(0,0) + &
                        nf**2*z*(177.73672456902338_dp)*Hr2(0,0) + nf*z**4*(742.3921963115378_dp)*Hr2(0,0) + &
                        nf**2*z**3*(976.2590687551649_dp)*Hr2(0,0) + nf*(1161.1654320987655_dp)*Hr2(0,0) + &
                        nf*z**2*(12379.234406783968_dp)*Hr2(0,0) + z*(49995.862884520044_dp)*Hr2(0,0) + &
                        z**3*(90583.1496599093_dp)*Hr2(0,0) + z**2*(-136625.82619872817_dp)*Hr2(0,1) + &
                        z**4*(-47763.548966558556_dp)*Hr2(0,1) + (-45351.788206985_dp)*Hr2(0,1) + &
                        nf*z**3*(-30868.19943199648_dp)*Hr2(0,1) + nf*z*(-3781.228426425445_dp)*Hr2(0,1) + &
                        nf**2*z**2*(-188.51760402377664_dp)*Hr2(0,1) + nf**2*(-56.98253315043438_dp)*Hr2(0,1) + &
                        nf**3*z**3*(-37.73863740283493_dp)*Hr2(0,1) + nf**2*z**4*(-31.2771833561957_dp)*Hr2(0,1) + &
                        nf**3*z*(-4.401280292638319_dp)*Hr2(0,1) + nf**2*z*(105.84142661179678_dp)*Hr2(0,1) + &
                        nf**2*z**3*(1874.793232738912_dp)*Hr2(0,1) + nf*z**4*(2458.6027739673827_dp)*Hr2(0,1) + &
                        nf*(3296.149245541838_dp)*Hr2(0,1) + nf*z**2*(9031.345221764976_dp)*Hr2(0,1) + &
                        z*(63766.83305836311_dp)*Hr2(0,1) + z**3*(223125.57053551328_dp)*Hr2(0,1) + &
                        z**2*(-136673.19307929234_dp)*Hr2(1,0) + (-47577.66108802618_dp)*Hr2(1,0) + &
                        z**4*(-47577.66108802618_dp)*Hr2(1,0) + nf*z*(-17817.63756435026_dp)*Hr2(1,0) + &
                        nf*z**3*(-17817.63756435026_dp)*Hr2(1,0) + nf**2*z**2*(-188.5176040237768_dp)*Hr2(1,0) + &
                        nf**2*(-44.12985825331504_dp)*Hr2(1,0) + nf**2*z**4*(-44.12985825331504_dp)*Hr2(1,0) + &
                        nf**3*z*(-21.069958847736626_dp)*Hr2(1,0) + nf**3*z**3*(-21.069958847736626_dp)*Hr2(1,0) + &
                        nf**2*z*(990.3173296753544_dp)*Hr2(1,0) + nf**2*z**3*(990.3173296753544_dp)*Hr2(1,0) + &
                        nf*(2877.3760097546105_dp)*Hr2(1,0) + nf*z**4*(2877.3760097546105_dp)*Hr2(1,0) + &
                        nf*z**2*(9031.345221764974_dp)*Hr2(1,0) + z**3*(152016.43753682286_dp)*Hr2(1,0) + &
                        z*(152016.4375368229_dp)*Hr2(1,0) + (z**2*(-273251.65239745635_dp) + (-93115.33717354356_dp) + &
                        z**4*(-93115.33717354356_dp) + nf**3*(z**3*(-42.13991769547325_dp) + z*(42.13991769547325_dp)) + &
                        nf**2*(z*(-1427.0800182898947_dp) + z**2*(-377.0352080475536_dp) + (-88.25971650663008_dp) + &
                        z**4*(-88.25971650663008_dp) + z**3*(1980.6346593507087_dp)) + nf*(z**3*(-34649.42785842192_dp) + &
                        z*(5077.23337587353_dp) + (5754.752019509221_dp) + z**4*(5754.752019509221_dp) + &
                        z**2*(18062.690443529944_dp)) + z*(177200.61549373966_dp) + z**3*(286892.40359387634_dp))*Hr2(1,1) + &
                        (z**2*(-32085.793155363484_dp) + z*(-3941.881728801615_dp) + (-1971.2_dp) + z**4*(-1059.18683127572_dp) + &
                        nf**3*(z**3*(-2.03676268861454_dp) + z*(2.03676268861454_dp)) + nf**2*z*((-81.24029873494895_dp) + &
                        z*(-75.95720164609052_dp) + z**2*(199.64956561499767_dp)) + nf*(z**3*(-5902.306333257158_dp) + &
                        z**4*(21.975186709343085_dp) + (77.74814814814813_dp) + z*(1154.453783803025_dp) + &
                        z**2*(3471.9572895852107_dp)) + z**3*(46861.54806494618_dp))*Hr3(0,0,0) + &
                        z**2*(-79116.19431750168_dp)*Hr3(0,0,1) + nf*z**3*(-10853.135802469134_dp)*Hr3(0,0,1) + &
                        (-7770.074074074072_dp)*Hr3(0,0,1) + z**4*(-4589.0370370370365_dp)*Hr3(0,0,1) + &
                        nf*z*(-2316.5541838134423_dp)*Hr3(0,0,1) + nf**2*z**2*(-151.9144032921811_dp)*Hr3(0,0,1) + &
                        nf**2*z*(30.98844688309708_dp)*Hr3(0,0,1) + nf*z**4*(102.25953360768175_dp)*Hr3(0,0,1) + &
                        nf**2*z**3*(290.7342173449169_dp)*Hr3(0,0,1) + nf*(292.240329218107_dp)*Hr3(0,0,1) + &
                        nf*z**2*(8070.502423411064_dp)*Hr3(0,0,1) + z*(25198.51645873047_dp)*Hr3(0,0,1) + &
                        z**3*(95597.90760896885_dp)*Hr3(0,0,1) + z**2*(-22786.736168267027_dp)*Hr3(0,1,0) + &
                        (-11038.024691358025_dp)*Hr3(0,1,0) + nf*z**3*(-8319.443987197072_dp)*Hr3(0,1,0) + &
                        z**4*(-7856.987654320987_dp)*Hr3(0,1,0) + nf**2*z*(-2.516689529035233_dp)*Hr3(0,1,0) + &
                        nf*z**4*(197.5893918609968_dp)*Hr3(0,1,0) + nf*z*(217.13763145861867_dp)*Hr3(0,1,0) + &
                        nf**2*z**3*(257.22908093278465_dp)*Hr3(0,1,0) + nf*(387.57018747142206_dp)*Hr3(0,1,0) + &
                        nf*z**2*(460.115226337449_dp)*Hr3(0,1,0) + z*(6496.36781739856_dp)*Hr3(0,1,0) + &
                        z**3*(80744.414287623_dp)*Hr3(0,1,0) + z**2*(-45573.47233653408_dp)*Hr3(0,1,1) + &
                        (-22076.04938271605_dp)*Hr3(0,1,1) + nf*z**3*(-16638.887974394143_dp)*Hr3(0,1,1) + &
                        z**4*(-15713.975308641973_dp)*Hr3(0,1,1) + nf**2*z*(-5.033379058070466_dp)*Hr3(0,1,1) + &
                        nf*z**4*(395.1787837219936_dp)*Hr3(0,1,1) + nf*z*(434.27526291723734_dp)*Hr3(0,1,1) + &
                        nf**2*z**3*(514.4581618655693_dp)*Hr3(0,1,1) + nf*(775.1403749428441_dp)*Hr3(0,1,1) + &
                        nf*z**2*(920.230452674898_dp)*Hr3(0,1,1) + z*(15264.035322337455_dp)*Hr3(0,1,1) + &
                        z**3*(156062.8176228142_dp)*Hr3(0,1,1) + z**2*(-15191.157445511353_dp)*Hr3(1,0,0) + &
                        (-6298.3374485596705_dp)*Hr3(1,0,0) + z**4*(-6298.3374485596705_dp)*Hr3(1,0,0) + &
                        nf*z*(-2700.768785246151_dp)*Hr3(1,0,0) + nf*z**3*(-2700.768785246151_dp)*Hr3(1,0,0) + &
                        nf**2*z*(84.90413046791646_dp)*Hr3(1,0,0) + nf**2*z**3*(84.90413046791646_dp)*Hr3(1,0,0) + &
                        nf*(195.0531931108063_dp)*Hr3(1,0,0) + nf*z**4*(195.0531931108063_dp)*Hr3(1,0,0) + &
                        nf*z**2*(306.7434842249651_dp)*Hr3(1,0,0) + z*(29343.153307081488_dp)*Hr3(1,0,0) + &
                        z**3*(29343.15330708149_dp)*Hr3(1,0,0) + z**2*(-45573.472336534054_dp)*Hr3(1,0,1) + &
                        (-18895.012345679013_dp)*Hr3(1,0,1) + z**4*(-18895.012345679013_dp)*Hr3(1,0,1) + &
                        nf*z*(-8102.306355738456_dp)*Hr3(1,0,1) + nf*z**3*(-8102.306355738454_dp)*Hr3(1,0,1) + &
                        nf**2*z*(254.7123914037494_dp)*Hr3(1,0,1) + nf**2*z**3*(254.7123914037494_dp)*Hr3(1,0,1) + &
                        nf*(585.1595793324187_dp)*Hr3(1,0,1) + nf*z**4*(585.1595793324187_dp)*Hr3(1,0,1) + &
                        nf*z**2*(920.2304526748972_dp)*Hr3(1,0,1) + z**3*(85663.4264725758_dp)*Hr3(1,0,1) + &
                        z*(85663.42647257583_dp)*Hr3(1,0,1) + z**2*(-45573.472336534054_dp)*Hr3(1,1,0) + &
                        (-18895.012345679013_dp)*Hr3(1,1,0) + z**4*(-18895.012345679013_dp)*Hr3(1,1,0) + &
                        nf*z*(-8102.306355738456_dp)*Hr3(1,1,0) + nf*z**3*(-8102.306355738454_dp)*Hr3(1,1,0) + &
                        nf**2*z*(254.7123914037494_dp)*Hr3(1,1,0) + nf**2*z**3*(254.7123914037494_dp)*Hr3(1,1,0) + &
                        nf*(585.1595793324187_dp)*Hr3(1,1,0) + nf*z**4*(585.1595793324187_dp)*Hr3(1,1,0) + &
                        nf*z**2*(920.2304526748972_dp)*Hr3(1,1,0) + z*(87240.78210502157_dp)*Hr3(1,1,0) + &
                        z**3*(87240.78210502157_dp)*Hr3(1,1,0) + z**2*(-91146.94467306811_dp)*Hr3(1,1,1) + &
                        (-37790.02469135803_dp)*Hr3(1,1,1) + z**4*(-37790.02469135803_dp)*Hr3(1,1,1) + &
                        nf*z**3*(-16204.612711476908_dp)*Hr3(1,1,1) + z*(-4599.85888936747_dp)*Hr3(1,1,1) + &
                        nf**2*z*(-509.4247828074988_dp)*Hr3(1,1,1) + nf**2*z**3*(509.4247828074988_dp)*Hr3(1,1,1) + &
                        nf*(1170.3191586648375_dp)*Hr3(1,1,1) + nf*z**4*(1170.3191586648375_dp)*Hr3(1,1,1) + &
                        nf*z**2*(1840.4609053497943_dp)*Hr3(1,1,1) + nf*z*(12023.513488797436_dp)*Hr3(1,1,1) + &
                        z**3*(171326.8529451516_dp)*Hr3(1,1,1) + z**2*(-18229.195793324187_dp)*Hr4(0,0,0,0) + &
                        nf*z**3*(-1507.2278006401464_dp)*Hr4(0,0,0,0) + nf*z*(-221.33991769547322_dp)*Hr4(0,0,0,0) + &
                        (-102.39999999999998_dp)*Hr4(0,0,0,0) + nf**2*z**2*(-42.982716049382724_dp)*Hr4(0,0,0,0) + &
                        z**4*(-19.977442463039168_dp)*Hr4(0,0,0,0) + nf**2*z*(8.427983539094651_dp)*Hr4(0,0,0,0) + &
                        nf**2*z**3*(34.55473251028807_dp)*Hr4(0,0,0,0) + z*(1269.1431184270689_dp)*Hr4(0,0,0,0) + &
                        nf*z**2*(1668.635390946502_dp)*Hr4(0,0,0,0) + z**3*(17927.475933546717_dp)*Hr4(0,0,0,0) + &
                        z**2*(-26816.409693644266_dp)*Hr4(0,0,0,1) + z*(-4076.738545953361_dp)*Hr4(0,0,0,1) + &
                        nf*z**3*(-2373.039231824417_dp)*Hr4(0,0,0,1) + (-500.62222222222215_dp)*Hr4(0,0,0,1) + &
                        z**4*(-115.86916628562717_dp)*Hr4(0,0,0,1) + nf**2*z*(-26.12674897119341_dp)*Hr4(0,0,0,1) + &
                        nf**2*z**3*(26.12674897119341_dp)*Hr4(0,0,0,1) + nf*z*(641.1471422039322_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(1432.2304526748972_dp)*Hr4(0,0,0,1) + z**3*(35734.86870903826_dp)*Hr4(0,0,0,1) + &
                        z**2*(-26830.807498856884_dp)*Hr4(0,0,1,0) + nf*z**3*(-1679.3927754915408_dp)*Hr4(0,0,1,0) + &
                        (-948.1481481481482_dp)*Hr4(0,0,1,0) + nf*z*(-352.1609510745313_dp)*Hr4(0,0,1,0) + &
                        z**4*(-288.424325560128_dp)*Hr4(0,0,1,0) + nf*z**2*(1432.2304526748972_dp)*Hr4(0,0,1,0) + &
                        z*(8516.506630086875_dp)*Hr4(0,0,1,0) + z**3*(28001.331504343856_dp)*Hr4(0,0,1,0) + &
                        z**2*(-53661.61499771377_dp)*Hr4(0,0,1,1) + nf*z**3*(-3358.7855509830815_dp)*Hr4(0,0,1,1) + &
                        (-1896.2962962962963_dp)*Hr4(0,0,1,1) + nf*z*(-704.3219021490626_dp)*Hr4(0,0,1,1) + &
                        z**4*(-576.848651120256_dp)*Hr4(0,0,1,1) + nf*z**2*(2864.4609053497943_dp)*Hr4(0,0,1,1) + &
                        z*(17033.01326017375_dp)*Hr4(0,0,1,1) + z**3*(56002.66300868771_dp)*Hr4(0,0,1,1) + &
                        z*(-1635.233653406494_dp)*Hr4(0,1,0,0) + (-839.9890260631001_dp)*Hr4(0,1,0,0) + &
                        nf*z**3*(-742.0722450845908_dp)*Hr4(0,1,0,0) + z**4*(-400.1731443377534_dp)*Hr4(0,1,0,0) + &
                        z**2*(-28.795610425238454_dp)*Hr4(0,1,0,0) + nf*z*(142.74897119341563_dp)*Hr4(0,1,0,0) + &
                        z**3*(11354.649596098154_dp)*Hr4(0,1,0,0) + z*(-4905.700960219481_dp)*Hr4(0,1,0,1) + &
                        (-2519.9670781893005_dp)*Hr4(0,1,0,1) + nf*z**3*(-2226.216735253772_dp)*Hr4(0,1,0,1) + &
                        z**4*(-1200.51943301326_dp)*Hr4(0,1,0,1) + z**2*(-86.38683127571888_dp)*Hr4(0,1,0,1) + &
                        nf*z*(428.2469135802468_dp)*Hr4(0,1,0,1) + z**3*(34063.94878829446_dp)*Hr4(0,1,0,1) + &
                        z*(-4905.700960219481_dp)*Hr4(0,1,1,0) + (-2519.9670781893005_dp)*Hr4(0,1,1,0) + &
                        nf*z**3*(-2226.216735253772_dp)*Hr4(0,1,1,0) + z**4*(-1200.51943301326_dp)*Hr4(0,1,1,0) + &
                        z**2*(-86.38683127571888_dp)*Hr4(0,1,1,0) + nf*z*(428.2469135802468_dp)*Hr4(0,1,1,0) + &
                        z**3*(34063.94878829446_dp)*Hr4(0,1,1,0) + z*(-9811.401920438962_dp)*Hr4(0,1,1,1) + &
                        (-5039.934156378601_dp)*Hr4(0,1,1,1) + nf*z**3*(-4452.433470507544_dp)*Hr4(0,1,1,1) + &
                        z**4*(-2401.03886602652_dp)*Hr4(0,1,1,1) + z**2*(-172.77366255143775_dp)*Hr4(0,1,1,1) + &
                        nf*z*(856.4938271604937_dp)*Hr4(0,1,1,1) + z**3*(68127.89757658892_dp)*Hr4(0,1,1,1) + &
                        (-310.0405426002133_dp)*Hr4(1,0,0,0) + z**4*(-310.0405426002133_dp)*Hr4(1,0,0,0) + &
                        nf*z*(-149.8308184727938_dp)*Hr4(1,0,0,0) + nf*z**3*(-149.8308184727938_dp)*Hr4(1,0,0,0) + &
                        z**2*(-14.397805212619666_dp)*Hr4(1,0,0,0) + z*(2429.8539856729144_dp)*Hr4(1,0,0,0) + &
                        z**3*(2429.8539856729153_dp)*Hr4(1,0,0,0) + (-1240.1621704008533_dp)*Hr4(1,0,0,1) + &
                        z**4*(-1240.1621704008533_dp)*Hr4(1,0,0,1) + nf*z*(-599.3232738911752_dp)*Hr4(1,0,0,1) + &
                        nf*z**3*(-599.3232738911752_dp)*Hr4(1,0,0,1) + z**2*(-57.59122085047866_dp)*Hr4(1,0,0,1) + &
                        z*(9719.415942691658_dp)*Hr4(1,0,0,1) + z**3*(9719.415942691661_dp)*Hr4(1,0,0,1) + &
                        (-1860.2432556012802_dp)*Hr4(1,0,1,0) + z**4*(-1860.2432556012802_dp)*Hr4(1,0,1,0) + &
                        nf*z*(-898.9849108367627_dp)*Hr4(1,0,1,0) + nf*z**3*(-898.9849108367627_dp)*Hr4(1,0,1,0) + &
                        z**2*(-86.38683127572064_dp)*Hr4(1,0,1,0) + z**3*(14579.123914037493_dp)*Hr4(1,0,1,0) + &
                        z*(14579.1239140375_dp)*Hr4(1,0,1,0) + (-3720.4865112025605_dp)*Hr4(1,0,1,1) + &
                        z**4*(-3720.4865112025605_dp)*Hr4(1,0,1,1) + nf*z*(-1797.9698216735253_dp)*Hr4(1,0,1,1) + &
                        nf*z**3*(-1797.9698216735253_dp)*Hr4(1,0,1,1) + z**2*(-172.77366255144128_dp)*Hr4(1,0,1,1) + &
                        z**3*(29158.247828074986_dp)*Hr4(1,0,1,1) + z*(29158.247828075_dp)*Hr4(1,0,1,1) + &
                        (-1240.1621704008533_dp)*Hr4(1,1,0,0) + z**4*(-1240.1621704008533_dp)*Hr4(1,1,0,0) + &
                        nf*z*(-599.3232738911752_dp)*Hr4(1,1,0,0) + nf*z**3*(-599.3232738911752_dp)*Hr4(1,1,0,0) + &
                        z**2*(-57.59122085047866_dp)*Hr4(1,1,0,0) + z*(9719.415942691658_dp)*Hr4(1,1,0,0) + &
                        z**3*(9719.415942691661_dp)*Hr4(1,1,0,0) + (-3720.4865112025605_dp)*Hr4(1,1,0,1) + &
                        z**4*(-3720.4865112025605_dp)*Hr4(1,1,0,1) + nf*z*(-1797.9698216735253_dp)*Hr4(1,1,0,1) + &
                        nf*z**3*(-1797.9698216735253_dp)*Hr4(1,1,0,1) + z**2*(-172.77366255144128_dp)*Hr4(1,1,0,1) + &
                        z**3*(29158.247828074986_dp)*Hr4(1,1,0,1) + z*(29158.247828075_dp)*Hr4(1,1,0,1) + &
                        (-3720.4865112025605_dp)*Hr4(1,1,1,0) + z**4*(-3720.4865112025605_dp)*Hr4(1,1,1,0) + &
                        nf*z*(-1797.9698216735253_dp)*Hr4(1,1,1,0) + nf*z**3*(-1797.9698216735253_dp)*Hr4(1,1,1,0) + &
                        z**2*(-172.77366255144128_dp)*Hr4(1,1,1,0) + z**3*(29158.247828074986_dp)*Hr4(1,1,1,0) + &
                        z*(29158.247828075_dp)*Hr4(1,1,1,0) + z*(-43089.00228623685_dp)*Hr4(1,1,1,1) + &
                        (-7440.973022405121_dp)*Hr4(1,1,1,1) + z**4*(-7440.973022405121_dp)*Hr4(1,1,1,1) + &
                        nf*z**3*(-3595.9396433470506_dp)*Hr4(1,1,1,1) + z**2*(-345.54732510288255_dp)*Hr4(1,1,1,1) + &
                        nf*z*(3595.9396433470506_dp)*Hr4(1,1,1,1) + z**3*(58316.49565614997_dp)*Hr4(1,1,1,1) + &
                        z**2*(-3006.8938271604934_dp)*Hr5(0,0,0,0,0) + nf*z**3*(-170.24526748971192_dp)*Hr5(0,0,0,0,0) + &
                        z*(-102.71995122694709_dp)*Hr5(0,0,0,0,0) + nf**2*z*(-0.8427983539094649_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z**3*(0.8427983539094649_dp)*Hr5(0,0,0,0,0) + nf*z*(22.33415637860082_dp)*Hr5(0,0,0,0,0) + &
                        nf*z**2*(147.9111111111111_dp)*Hr5(0,0,0,0,0) + z**3*(3141.5776863283036_dp)*Hr5(0,0,0,0,0) + &
                        z**2*(-7903.762962962963_dp)*Hr5(0,0,0,0,1) + nf*z**3*(-225.86995884773663_dp)*Hr5(0,0,0,0,1) + &
                        nf*z*(-69.95226337448558_dp)*Hr5(0,0,0,0,1) + nf*z**2*(295.8222222222222_dp)*Hr5(0,0,0,0,1) + &
                        z*(525.3443072702331_dp)*Hr5(0,0,0,0,1) + z**3*(7570.202103337905_dp)*Hr5(0,0,0,0,1) + &
                        z**2*(-5669.925925925924_dp)*Hr5(0,0,0,1,0) + z*(-1032.6620941929582_dp)*Hr5(0,0,0,1,0) + &
                        nf*z**3*(-77.95884773662551_dp)*Hr5(0,0,0,1,0) + nf*z*(77.95884773662551_dp)*Hr5(0,0,0,1,0) + &
                        z**3*(7182.046639231823_dp)*Hr5(0,0,0,1,0) + z**2*(-11339.851851851849_dp)*Hr5(0,0,0,1,1) + &
                        z*(-2065.3241883859164_dp)*Hr5(0,0,0,1,1) + nf*z**3*(-155.91769547325103_dp)*Hr5(0,0,0,1,1) + &
                        nf*z*(155.91769547325103_dp)*Hr5(0,0,0,1,1) + z**3*(14364.093278463646_dp)*Hr5(0,0,0,1,1) + &
                        z**2*(-3779.9506172839497_dp)*Hr5(0,0,1,0,0) + z*(1039.7634506934917_dp)*Hr5(0,0,1,0,0) + &
                        z**3*(3379.465325407712_dp)*Hr5(0,0,1,0,0) + z**2*(-11339.85185185185_dp)*Hr5(0,0,1,0,1) + &
                        z*(3119.2903520804753_dp)*Hr5(0,0,1,0,1) + z**3*(10138.395976223137_dp)*Hr5(0,0,1,0,1) + &
                        z**2*(-11339.85185185185_dp)*Hr5(0,0,1,1,0) + z*(3119.2903520804753_dp)*Hr5(0,0,1,1,0) + &
                        z**3*(10138.395976223137_dp)*Hr5(0,0,1,1,0) + z**2*(-22679.7037037037_dp)*Hr5(0,0,1,1,1) + &
                        z*(6238.5807041609505_dp)*Hr5(0,0,1,1,1) + z**3*(20276.791952446274_dp)*Hr5(0,0,1,1,1) + &
                        z*(-345.19615912208496_dp)*Hr5(0,1,0,0,0) + z**3*(824.654778235025_dp)*Hr5(0,1,0,0,0) + &
                        z*(-1380.7846364883399_dp)*Hr5(0,1,0,0,1) + z**3*(3298.6191129401_dp)*Hr5(0,1,0,0,1) + &
                        z*(-2071.1769547325107_dp)*Hr5(0,1,0,1,0) + z**3*(4947.928669410151_dp)*Hr5(0,1,0,1,0) + &
                        z*(-4142.353909465021_dp)*Hr5(0,1,0,1,1) + z**3*(9895.857338820302_dp)*Hr5(0,1,0,1,1) + &
                        z*(-1380.7846364883399_dp)*Hr5(0,1,1,0,0) + z**3*(3298.6191129401_dp)*Hr5(0,1,1,0,0) + &
                        z*(-4142.353909465021_dp)*Hr5(0,1,1,0,1) + z**3*(9895.857338820302_dp)*Hr5(0,1,1,0,1) + &
                        z*(-4142.353909465021_dp)*Hr5(0,1,1,1,0) + z**3*(9895.857338820302_dp)*Hr5(0,1,1,1,0) + &
                        z*(-8284.707818930043_dp)*Hr5(0,1,1,1,1) + z**3*(19791.714677640604_dp)*Hr5(0,1,1,1,1) + &
                        z*(95.89172382258796_dp)*Hr5(1,0,0,0,0) + z**3*(95.891723822588_dp)*Hr5(1,0,0,0,0) + &
                        z**3*(479.4586191129401_dp)*Hr5(1,0,0,0,1) + z*(479.45861911294014_dp)*Hr5(1,0,0,0,1) + &
                        z**3*(958.9172382258802_dp)*Hr5(1,0,0,1,0) + z*(958.9172382258803_dp)*Hr5(1,0,0,1,0) + &
                        z**3*(1917.8344764517603_dp)*Hr5(1,0,0,1,1) + z*(1917.8344764517606_dp)*Hr5(1,0,0,1,1) + &
                        z**3*(958.9172382258802_dp)*Hr5(1,0,1,0,0) + z*(958.9172382258803_dp)*Hr5(1,0,1,0,0) + &
                        z*(2876.7517146776404_dp)*Hr5(1,0,1,0,1) + z**3*(2876.7517146776404_dp)*Hr5(1,0,1,0,1) + &
                        z*(2876.7517146776404_dp)*Hr5(1,0,1,1,0) + z**3*(2876.7517146776404_dp)*Hr5(1,0,1,1,0) + &
                        z*(5753.503429355281_dp)*Hr5(1,0,1,1,1) + z**3*(5753.503429355281_dp)*Hr5(1,0,1,1,1) + &
                        z**3*(479.4586191129401_dp)*Hr5(1,1,0,0,0) + z*(479.45861911294014_dp)*Hr5(1,1,0,0,0) + &
                        z**3*(1917.8344764517603_dp)*Hr5(1,1,0,0,1) + z*(1917.8344764517606_dp)*Hr5(1,1,0,0,1) + &
                        z*(2876.7517146776404_dp)*Hr5(1,1,0,1,0) + z**3*(2876.7517146776404_dp)*Hr5(1,1,0,1,0) + &
                        z*(5753.503429355281_dp)*Hr5(1,1,0,1,1) + z**3*(5753.503429355281_dp)*Hr5(1,1,0,1,1) + &
                        z**3*(1917.8344764517603_dp)*Hr5(1,1,1,0,0) + z*(1917.8344764517606_dp)*Hr5(1,1,1,0,0) + &
                        z*(5753.503429355281_dp)*Hr5(1,1,1,0,1) + z**3*(5753.503429355281_dp)*Hr5(1,1,1,0,1) + &
                        z*(5753.503429355281_dp)*Hr5(1,1,1,1,0) + z**3*(5753.503429355281_dp)*Hr5(1,1,1,1,0) + &
                        z*(-11507.006858710562_dp)*Hr5(1,1,1,1,1) + z**3*(11507.006858710562_dp)*Hr5(1,1,1,1,1))/(z*((-1._dp) + &
                        z*(1._dp)))

                    case (2)
                        if (z == 1._dp) then
                            Ibar_select = (-180685.19013962927_dp) + nf**2*(-6963.06819459863_dp) + nf**4*(-9.557270391251393_dp) + &
                            nf**5*(0.0877914951989026_dp) + nf**3*(379.4751057355212_dp) + nf*(58934.28136172893_dp)

                        else
                            Ibar_select = (-180685.19013962927_dp) + nf**2*(-6963.06819459863_dp) + nf**4*(-9.557270391251393_dp) + &
                            nf**5*(0.0877914951989026_dp) + nf**3*(379.4751057355212_dp) + nf*(58934.28136172893_dp)

                        endif
                    case (3)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-45522.08963549428_dp) + nf**4*(-6.449748513946044_dp) + &
                            nf**5*(0.11705532693187014_dp) + nf**3*(61.38659186066184_dp) + nf**2*(2090.1549053124127_dp) + &
                            (240877.0920262052_dp)

                        else
                            Ibar_select = nf*(-45522.08963549428_dp) + nf**4*(-6.449748513946044_dp) + &
                            nf**5*(0.11705532693187014_dp) + nf**3*(61.38659186066184_dp) + nf**2*(2090.1549053124127_dp) + &
                            (240877.0920262052_dp)

                        endif
                    case (4)
                        if (z == 1._dp) then
                            Ibar_select = (-96341.45652644028_dp) + nf*(-8066.341642534953_dp) + nf**3*(-219.0339277549154_dp) + &
                            nf**4*(4.276421277244323_dp) + nf**2*(3305.2265280821166_dp)

                        else
                            Ibar_select = (-96341.45652644028_dp) + nf*(-8066.341642534953_dp) + nf**3*(-219.0339277549154_dp) + &
                            nf**4*(4.276421277244323_dp) + nf**2*(3305.2265280821166_dp)

                        endif
                    case (5)
                        if (z == 1._dp) then
                            Ibar_select = (-54845.89405006837_dp) + nf**2*(-1703.8573388203017_dp) + nf**3*(42.13991769547325_dp) + &
                            nf*(19863.330617147727_dp)

                        else
                            Ibar_select = (-54845.89405006837_dp) + nf**2*(-1703.8573388203017_dp) + nf**3*(42.13991769547325_dp) + &
                            nf*(19863.330617147727_dp)

                        endif
                    case (6)
                        if (z == 1._dp) then
                            Ibar_select = nf*(-4704.687700045724_dp) + nf**2*(169.80826093583295_dp) + (29321.11863908652_dp)
                        else
                            Ibar_select = nf*(-4704.687700045724_dp) + nf**2*(169.80826093583295_dp) + (29321.11863908652_dp)
                        endif
                    case (7)
                        if (z == 1._dp) then
                            Ibar_select = (-4225.229080932784_dp) + nf*(299.6616369455876_dp)
                        else
                            Ibar_select = (-4225.229080932784_dp) + nf*(299.6616369455876_dp)
                        endif
                    case (8)
                        if (z == 1._dp) then
                            Ibar_select = (191.78344764517604_dp)
                        else
                            Ibar_select = (191.78344764517604_dp)
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
