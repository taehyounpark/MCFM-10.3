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
                        Ibar_select = z*(z*(-2._dp) + (2._dp))
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
                        Ibar_select = z**2*(-2._dp) + (-1._dp) + z*(2._dp)
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
                        Ibar_select = z*(-116.91531083870387_dp) + z**3*(-103.13731364576645_dp) + z**4*(-16.769336143308745_dp) + &
                        z**5*(4.082696945337621_dp) + (5.9516385763253_dp)/z + (6.51093768146994_dp) + z**2*(204.80522133400248_dp) &
                        + (z**2*(-24.73187968628803_dp) + z*(6.666666666666667_dp) + (11.333333333333334_dp) + &
                        z**3*(21.54247493927264_dp))*Log(z) + (z**2*(-20.465598027127005_dp) + (-1.1666666666666667_dp) + &
                        z*(11.666666666666666_dp) + z**3*(15.293497059840886_dp))*Log(z)**2 + (z**2*(-1.0086275130599969_dp) + &
                        z**3*(-0.5403004572175049_dp) + (0.7777777777777778_dp) + z*(2.4444444444444446_dp))*Log(z)**3 + &
                        (z**2*(-179.51737530784183_dp) + (-66.54875304640971_dp) + z**3*(54.8176444640494_dp) + &
                        z*(189.5818172235355_dp))*Log(z*(-1._dp) + (1._dp)) + (z*(-8.33947549500082_dp) + &
                        z**3*(-2.5205581395348835_dp) + (0.5761253443996348_dp) + z**2*(10.283908290136068_dp))*Log(z*(-1._dp) + &
                        (1._dp))**2 + (z**2*(-7.76637769548793_dp) + (-2.489245468238055_dp) + z**3*(2.916343891402715_dp) + &
                        z*(7.8948348278788245_dp))*Log(z*(-1._dp) + (1._dp))**3

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
                        Ibar_select = z*(-50.80691404374448_dp) + (-17.333333333333332_dp)/z + z**4*(-10.858290621969802_dp) + &
                        z**3*(-9.145891105850524_dp) + z**5*(1.6911060433295324_dp) + (21.632893830259437_dp) + &
                        z**2*(53.2503853487307_dp) + (z**2*(-72.65790509942902_dp) + (-7.333333333333333_dp) + &
                        z*(17.333333333333332_dp) + z**3*(18.23998653651969_dp))*Log(z) + (z**2*(-5.319453991947364_dp) + &
                        z**3*(1.9444742694991948_dp) + (4.666666666666667_dp) + z*(14.666666666666666_dp))*Log(z)**2 + &
                        (z*(-5.497916892354098_dp) + (-2.352382972419054_dp) + z**2*(2.0529827019653624_dp) + &
                        z**3*(5.797317162807791_dp))*Log(z*(-1._dp) + (1._dp)) + (z*(-4.526119553055392_dp) + &
                        z**3*(0.764066016504126_dp) + (2.645092784779759_dp) + z**2*(4.45029408510484_dp))*Log(z*(-1._dp) + &
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
                        Ibar_select = z**2*(-31.000000000000007_dp) + (2.333333333333332_dp) + (4._dp)/z + z*(26.66666666666667_dp) &
                        + (z**2*(-5.333333333333333_dp) + (4.666666666666667_dp) + z*(26.666666666666668_dp))*Log(z) + &
                        (z*(-17.333333333333332_dp) + (8.666666666666668_dp) + z**2*(17.333333333333332_dp))*Log(z*(-1._dp) + &
                        (1._dp))

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
                        Ibar_select = (z**2*(-694771.5000060453_dp) + z*(-1159.4130990003464_dp) + (-1122.3366645410142_dp) + &
                        Nf*(z**2*(-758.3733709715073_dp) + z**5*(-287.42787276181235_dp) + z**6*(-2.0060722275487377_dp) + &
                        (11.954578641471228_dp) + z*(256.9582959516055_dp) + z**3*(417.1572453922548_dp) + &
                        z**4*(442.25008552856656_dp)) + z**6*(1265.6419248539194_dp) + z**5*(4365.665280900331_dp) + &
                        z**3*(95793.58871805332_dp) + z**4*(594981.9678463419_dp))/z + (z**3*(-768074.2462285938_dp) + &
                        z**2*(-536598.7264135503_dp) + (-2355.8728839622704_dp) + z*(-1769.1986942544231_dp) + &
                        (-177.71551065150712_dp)/z + Nf*(z**2*(-561.2871173243674_dp) + z*(-71.98083276696133_dp) + &
                        (216.89773656629333_dp) + z**3*(453.907736333136_dp)))*Log(z) + (z**2*(-91729.69787691534_dp) + &
                        (-100.67928621642898_dp) + z*(-74.87166810714174_dp) + Nf*(z**2*(-202.06556363074068_dp) + &
                        z**3*(-119.87409778812572_dp) + z*(-75.78325472728147_dp) + (67.23388703965145_dp)) + &
                        z**3*(221190.33563270248_dp))*Log(z)**2 + (z**3*(-34083.821526722284_dp) + z**2*(-8962.498892970549_dp) + &
                        (-172.79383354273608_dp) + Nf*(z*(-21.637860082304528_dp) + z**2*(-20.142619680851062_dp) + &
                        (14.263374485596708_dp) + z**3*(120.4826609640504_dp)) + z*(170.13609352560715_dp))*Log(z)**3 + &
                        (z**2*(-463.9959882901442_dp) + Nf*(z**3*(-4.682166732206056_dp) + z*(-3.654320987654321_dp) + &
                        z**2*(-3.104622898300718_dp) + (1.271604938271605_dp)) + (5.080246913580247_dp) + z*(65.88271604938272_dp) + &
                        z**3*(3952.9617292300754_dp))*Log(z)**4 + (z**3*(-259.0994455049541_dp) + z**2*(-10.131536475145255_dp) + &
                        (-1.9851851851851852_dp) + Nf*(z*(-0.35555555555555557_dp) + z**2*(-0.05652110413319702_dp) + &
                        (0.17777777777777778_dp) + z**3*(3.236072206660442_dp)) + z*(10.264197530864198_dp))*Log(z)**5 + &
                        (z**2*(-1.5566149919806875e6_dp) + (-693546.9978914413_dp) + Nf*(z**2*(-1006.8647775213462_dp) + &
                        (-467.64705600611205_dp) + z**3*(267.94036979969184_dp) + z*(1228.8232996089628_dp)) + &
                        z**3*(431546.870521429_dp) + z*(1.8185152483274792e6_dp))*Log(z*(-1._dp) + (1._dp)) + &
                        (z*(-237028.73863117665_dp) + z**3*(-107803.26885563733_dp) + Nf*(z*(-276.84972079414723_dp) + &
                        z**3*(-116.9679720483331_dp) + (83.34828178031447_dp) + z**2*(311.9507223282516_dp)) + (64530.1134571651_dp) &
                        + z**2*(280267.68273180583_dp))*Log(z*(-1._dp) + (1._dp))**2 + (z**2*(-74362.79530918912_dp) + &
                        (-27325.900826678386_dp) + Nf*(z**2*(-34.8770972460037_dp) + (-13.356512894823108_dp) + &
                        z**3*(9.97839917147507_dp) + z*(37.43216570186203_dp)) + z**3*(23529.678926234235_dp) + &
                        z*(78171.94971539348_dp))*Log(z*(-1._dp) + (1._dp))**3 + (z*(-5013.1927508007975_dp) + &
                        z**3*(-1791.2076153131773_dp) + Nf*(z*(-7.3892448014815235_dp) + z**3*(-2.6122063956230543_dp) + &
                        (2.388519202929235_dp) + z**2*(7.458611006521022_dp)) + (1613.8938023117116_dp) + &
                        z**2*(5192.867674913375_dp))*Log(z*(-1._dp) + (1._dp))**4 + (z**2*(-799.8276852791904_dp) + &
                        (-269.2723548239249_dp) + z**3*(264.81470226466973_dp) + z*(803.3594119125196_dp))*Log(z*(-1._dp) + &
                        (1._dp))**5

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
                        Ibar_select = (z**3*(-129853.46299731074_dp) + z**5*(-6806.1606047527985_dp) + z*(-2983.878573188233_dp) + &
                        z**6*(857.6941500905778_dp) + (1177.71869631644_dp) + Nf*(z**3*(-2289.180446617168_dp) + &
                        z**2*(-237.23538046745549_dp) + z**5*(-105.80302534215178_dp) + (-0.7456847915865921_dp) + &
                        z**6*(0.8514062046970137_dp) + z*(225.71154970332464_dp) + z**4*(2455.2093695051494_dp)) + &
                        z**2*(68034.69546514687_dp) + z**4*(69266.04261275455_dp))/z + (z*(-406.09413789678945_dp) + &
                        (-142.03277970340008_dp) + Nf*(z**3*(-1318.2484107185967_dp) + z**2*(-789.4531928480204_dp) + &
                        z*(-41.05622138871932_dp) + (170.13855597230702_dp)) + (220.75299624923696_dp)/z + &
                        z**2*(2706.0171943820083_dp) + z**3*(4337.465476956942_dp))*Log(z) + (z**2*(-943.3366050665269_dp) + &
                        (-686.5753907919999_dp) + z*(368.46612175335054_dp) + Nf*(z**2*(-108.8747624649353_dp) + &
                        z*(-44.51851851851852_dp) + (58.03703703703704_dp) + z**3*(447.7548033189073_dp)) + &
                        z**3*(7877.884784996017_dp))*Log(z)**2 + (z**2*(-150.96362083673992_dp) + Nf*(z**3*(-50.083411045487516_dp) &
                        + z*(-16.74074074074074_dp) + z**2*(-14.79422398785877_dp) + (7.037037037037037_dp)) + &
                        (19.814814814814817_dp) + z*(210.96296296296296_dp) + z**3*(748.3480083894103_dp))*Log(z)**3 + &
                        ((-12.407407407407407_dp) + z**2*(-4.964967156709415_dp) + Nf*(z*(-2.2222222222222223_dp) + &
                        z**2*(-0.30381875741986547_dp) + (1.1111111111111112_dp) + z**3*(10.512828074314362_dp)) + &
                        z*(65.44444444444444_dp) + z**3*(245.79940434477922_dp))*Log(z)**4 + (z*(-164498.7812499988_dp) + &
                        z**3*(-36531.95653134886_dp) + Nf*(z**3*(-15.811826411840782_dp) + z*(-5.1282702059057215_dp) + &
                        (7.8804441192546975_dp) + z**2*(11.417166543662496_dp)) + (63642.33315537591_dp) + &
                        z**2*(137401.6990514342_dp))*Log(z*(-1._dp) + (1._dp)) + (z**2*(-41550.78233079086_dp) + &
                        (-11006.106814397284_dp) + Nf*(z**2*(-8.082547674742711_dp) + (-3.7232480914246526_dp) + &
                        z**3*(-0.5610909490817112_dp) + z*(8.663183011545371_dp)) + z**3*(15337.26041997336_dp) + &
                        z*(37323.251826545704_dp))*Log(z*(-1._dp) + (1._dp))**2 + (z*(-4693.423208402095_dp) + &
                        z**3*(-1371.582746152572_dp) + Nf*(z**2*(-1.0949968769720386_dp) + (-0.6014752054707988_dp) + &
                        z**3*(-0.12360953461975029_dp) + z*(1.0793408763218468_dp)) + (1679.2906014951318_dp) + &
                        z**2*(4395.715353059536_dp))*Log(z*(-1._dp) + (1._dp))**3 + (z**2*(-1102.370407606106_dp) + &
                        (-364.2159460270891_dp) + z**3*(364.8179715302491_dp) + z*(1093.2498635844274_dp))*Log(z*(-1._dp) + &
                        (1._dp))**4

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
                        Ibar_select = (z**2*(-1761.738355258621_dp) + z**4*(-1623.8628770433393_dp) + (-526.482813147809_dp) + &
                        z**5*(-511.5079224819828_dp) + z*(-29.82600718599596_dp) + Nf*(z**2*(-126.39616515865728_dp) + &
                        z**4*(-23.250231248743212_dp) + z**5*(-1.66601235914213_dp) + z**6*(-0.13370149687545416_dp) + &
                        (1.876543209876543_dp) + z*(72.35424987622878_dp) + z**3*(78.19436723004364_dp)) + &
                        z**6*(85.32015506632303_dp) + z**3*(4366.624363057692_dp))/z + (z**2*(-830.997242227052_dp) + &
                        (-704.081670491375_dp) + z*(-250.14255445319202_dp) + (-104._dp)/z + Nf*(z*(-34.51851851851852_dp) + &
                        z**3*(17.158043940467753_dp) + z**2*(38.95244150559512_dp) + (50.370370370370374_dp)) + &
                        z**3*(559.775732516814_dp))*Log(z) + (z**2*(-224.0978772125856_dp) + Nf*(z*(-20.444444444444443_dp) + &
                        z**3*(-3.468035118090763_dp) + z**2*(-2.232350732256853_dp) + (8.88888888888889_dp)) + &
                        (13.333333333333336_dp) + z**3*(324.3509263900062_dp) + z*(352.00000000000006_dp))*Log(z)**2 + &
                        ((-19.851851851851855_dp) + z**2*(-10.81862877519045_dp) + Nf*(z*(-3.555555555555556_dp) + &
                        z**2*(0.06560190073917635_dp) + z**3*(1.2826381569040808_dp) + (1.777777777777778_dp)) + &
                        z**3*(31.844111969111964_dp) + z*(114.96296296296295_dp))*Log(z)**3 + (z**2*(-3162.5688947847793_dp) + &
                        (-947.2665819032279_dp) + Nf*(z**2*(-19.936639321570453_dp) + (-8.676849507603912_dp) + &
                        z**3*(-1.8515865744982094_dp) + z*(20.83544577404294_dp)) + z**3*(1290.873378662998_dp) + &
                        z*(3122.295431358342_dp))*Log(z*(-1._dp) + (1._dp)) + (z*(-255.23219381134624_dp) + &
                        z**3*(-110.5392076214939_dp) + Nf*(z**2*(-2.4704520791863733_dp) + (-1.194809684709898_dp) + &
                        z**3*(0.08226564168268204_dp) + z*(2.471885011102477_dp)) + (80.79093753937059_dp) + &
                        z**2*(296.6471305601362_dp))*Log(z*(-1._dp) + (1._dp))**2 + (z**2*(-248.67192543336725_dp) + &
                        (-94.2694498179906_dp) + z**3*(62.75679336324389_dp) + z*(251.29569299922508_dp))*Log(z*(-1._dp) + &
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
                        Ibar_select = z**2*(-249.87674040142244_dp) + (-160.4130925141413_dp) + z**3*(-20.612096754936847_dp) + &
                        z**4*(-0.4942348008385712_dp) + z**5*(-0.010576844845821838_dp) + Nf*(z*(-19.851851851851844_dp) + &
                        (-1.6790123456790123_dp)/z + (7.6296296296296315_dp) + z**2*(13.012345679012347_dp)) + &
                        (90.66666666666666_dp)/z + z*(417.5626614919019_dp) + (z**2*(-72.97707069550948_dp) + &
                        z**3*(-39.692601948754955_dp) + Nf*(z*(-6.518518518518518_dp) + (1.9259259259259256_dp)) + &
                        (12.22222222222222_dp) + (16._dp)/z + z*(303.7037037037037_dp))*Log(z) + ((-9.925925925925926_dp) + &
                        z**2*(-4.904683840749415_dp) + Nf*(z*(-1.7777777777777775_dp) + (0.8888888888888887_dp)) + &
                        z**3*(4.875785001847063_dp) + z*(107.85185185185183_dp))*Log(z)**2 + (z*(-233.53429955542347_dp) + &
                        Nf*(z**2*(-5.037037037037037_dp) + (-2.518518518518518_dp) + z*(5.037037037037037_dp)) + &
                        z**3*(53.08120874462546_dp) + (64.49293933520963_dp) + z**2*(142.40459592003285_dp))*Log(z*(-1._dp) + &
                        (1._dp)) + (z**2*(-53.77538749685294_dp) + (-31.221091862450702_dp) + z**3*(-8.426555886502587_dp) + &
                        z*(54.01562783839883_dp))*Log(z*(-1._dp) + (1._dp))**2

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
                        Ibar_select = z**2*(-1.0081713193428655e6_dp) + (-112235.81636832503_dp) + nf*z*(-7493.229574468216_dp) + &
                        nf*(-3573.150030592713_dp) + (nf*(-881.2405095436995_dp))/z + nf**2*(-434.54706940171667_dp) + &
                        nf**2*z**2*(-309.753517266152_dp) + (nf**2*(0.3353562670682706_dp))/z + nf**2*z*(480.4808830410886_dp) + &
                        nf*z**2*(13329.638633372593_dp) + (53885.37423687335_dp)/z + z*(54435.55681434879_dp) + &
                        ((-8.772981689857207_dp)*(z*(-1133.1922932411449_dp) + nf*(z**2*(-491.17374727753617_dp) + &
                        z**3*(-423.84967320346215_dp) + z*(-227.46187363876808_dp) + (-20.37037037037037_dp)) + nf**2*z*((0.5_dp) + &
                        z*(1._dp) + z**2*(1._dp)) + (117.03012023224566_dp) + z**2*(3833.3801473046365_dp) + &
                        z**3*(6324.3125297809665_dp))*Hr1(-1))/z + nf*z**2*(-12645.906117325805_dp)*Hr1(0) + &
                        nf*(-6420.113992397535_dp)*Hr1(0) + nf*z*(-3096.728987784415_dp)*Hr1(0) + &
                        nf**2*(-260.47351805845256_dp)*Hr1(0) + (nf*(-32.381204042022915_dp)*Hr1(0))/z + &
                        nf**2*z*(52.52950247436462_dp)*Hr1(0) + nf**2*z**2*(116.49794238683127_dp)*Hr1(0) + &
                        ((12814.827186642691_dp)*Hr1(0))/z + z*(31687.383676041092_dp)*Hr1(0) + (167219.77026446795_dp)*Hr1(0) + &
                        z**2*(905349.1496050835_dp)*Hr1(0) + z**2*(-61682.13256298966_dp)*Hr1(1) + (-10019.85638845875_dp)*Hr1(1) + &
                        nf*z*(-2722.9733919969294_dp)*Hr1(1) + nf**2*z*(-87.44512530602756_dp)*Hr1(1) + &
                        nf**2*z**2*(-20.53841378862264_dp)*Hr1(1) + (nf*(-17.261114192176436_dp)*Hr1(1))/z + &
                        (nf**2*(4.576131687242798_dp)*Hr1(1))/z + nf**2*(91.13614289992736_dp)*Hr1(1) + &
                        nf*(822.0343743993235_dp)*Hr1(1) + nf*z**2*(2428.3060848246523_dp)*Hr1(1) + &
                        ((10325.415698690585_dp)*Hr1(1))/z + z*(70519.80071428683_dp)*Hr1(1) + &
                        nf*z*(-1421.2230337568676_dp)*Hr2(-1,-1) + nf*z**2*(-1421.2230337568676_dp)*Hr2(-1,-1) + &
                        nf*(-710.6115168784338_dp)*Hr2(-1,-1) + ((6819.531100249002_dp)*Hr2(-1,-1))/z + &
                        (33563.83362010985_dp)*Hr2(-1,-1) + z**2*(100673.59967074703_dp)*Hr2(-1,-1) + &
                        z*(105389.56471691692_dp)*Hr2(-1,-1) + z**2*(-252761.11419159643_dp)*Hr2(-1,0) + &
                        z*(-248881.90262270483_dp)*Hr2(-1,0) + (-48362.94232259875_dp)*Hr2(-1,0) + &
                        ((-38475.15436991376_dp)*Hr2(-1,0))/z + nf*z*(-3160.1185589253455_dp)*Hr2(-1,0) + &
                        nf*(-2448.1580448947716_dp)*Hr2(-1,0) + nf*z**2*(-268.66176880188846_dp)*Hr2(-1,0) + &
                        nf**2*(8.88888888888889_dp)*Hr2(-1,0) + nf**2*z*(17.77777777777778_dp)*Hr2(-1,0) + &
                        nf**2*z**2*(17.77777777777778_dp)*Hr2(-1,0) + (nf*(166.22222222222223_dp)*Hr2(-1,0))/z + &
                        z**2*(-174689.27661534474_dp)*Hr2(0,-1) + z*(-115806.4472341878_dp)*Hr2(0,-1) + &
                        (-9151.806723329963_dp)*Hr2(0,-1) + ((-4737.410112522892_dp)*Hr2(0,-1))/z + &
                        nf*(773.2408583865811_dp)*Hr2(0,-1) + nf*z*(1503.104196195535_dp)*Hr2(0,-1) + &
                        nf*z**2*(1516.7510566019794_dp)*Hr2(0,-1) + z**2*(-92135.48370632797_dp)*Hr2(0,0) + &
                        (-17930.002555808744_dp)*Hr2(0,0) + nf*z*(-3570.9683392227794_dp)*Hr2(0,0) + &
                        nf*(-1212.9470396454997_dp)*Hr2(0,0) + nf**2*(-150.0119017902365_dp)*Hr2(0,0) + &
                        nf**2*z**2*(-48.49382716049383_dp)*Hr2(0,0) + nf**2*z*(198.76454432121375_dp)*Hr2(0,0) + &
                        ((1324.5179774954217_dp)*Hr2(0,0))/z + nf*z**2*(12402.057466366434_dp)*Hr2(0,0) + &
                        z*(30977.06396369575_dp)*Hr2(0,0) + z*(-22439.45522924295_dp)*Hr2(0,1) + &
                        ((-4687.152684170228_dp)*Hr2(0,1))/z + nf*z*(-2204.566125679464_dp)*Hr2(0,1) + &
                        nf**2*z**2*(5.827160493827161_dp)*Hr2(0,1) + nf**2*z*(19.234567901234566_dp)*Hr2(0,1) + &
                        (nf*(38.51851851851852_dp)*Hr2(0,1))/z + nf**2*(68.82716049382717_dp)*Hr2(0,1) + &
                        nf*(1292.7647902590236_dp)*Hr2(0,1) + nf*z**2*(1510.794095362739_dp)*Hr2(0,1) + &
                        (6743.541590408113_dp)*Hr2(0,1) + z**2*(108871.6242428635_dp)*Hr2(0,1) + z*(-173522.6545439702_dp)*Hr2(1,0) &
                        + ((-32313.303350129027_dp)*Hr2(1,0))/z + nf*z**2*(-2707.6627942163072_dp)*Hr2(1,0) + &
                        nf**2*z*(-68.34567901234568_dp)*Hr2(1,0) + (nf**2*(-13.037037037037036_dp)*Hr2(1,0))/z + &
                        nf**2*(15.209876543209877_dp)*Hr2(1,0) + nf**2*z**2*(58.27160493827161_dp)*Hr2(1,0) + &
                        nf*(398.5965864309411_dp)*Hr2(1,0) + (nf*(771.522633744856_dp)*Hr2(1,0))/z + &
                        nf*z*(1661.7204073850316_dp)*Hr2(1,0) + (19591.64765046048_dp)*Hr2(1,0) + &
                        z**2*(198022.49689281458_dp)*Hr2(1,0) + ((9.876543209876543_dp)*(z*(-1066.5926005738806_dp) + &
                        nf**2*z*(z*(-1._dp) + (0.5_dp) + z**2*(1._dp)) + nf*(z*(-43.45696044010894_dp) + &
                        z**3*(-16.927254213551205_dp) + (-3.326666666666667_dp) + z**2*(49.253920880217876_dp)) + &
                        (366.52865098275817_dp) + z**3*(839.5923635426063_dp) + z**2*(1552.994371080597_dp))*Hr2(1,1))/z + &
                        z*(-27933.173700505347_dp)*Hr3(-1,-1,-1) + z**2*(-27933.173700505347_dp)*Hr3(-1,-1,-1) + &
                        (-13966.586850252674_dp)*Hr3(-1,-1,-1) + nf**2*z*(-10.666666666666666_dp)*Hr3(-1,-1,0) + &
                        nf**2*z**2*(-10.666666666666666_dp)*Hr3(-1,-1,0) + nf**2*(-5.333333333333333_dp)*Hr3(-1,-1,0) + &
                        (nf*(39.50617283950617_dp)*Hr3(-1,-1,0))/z + nf*(882.4691358024692_dp)*Hr3(-1,-1,0) + &
                        nf*z**2*(1511.111111111111_dp)*Hr3(-1,-1,0) + nf*z*(1794.8641975308642_dp)*Hr3(-1,-1,0) + &
                        ((6877.037037037037_dp)*Hr3(-1,-1,0))/z + (32426.187435523447_dp)*Hr3(-1,-1,0) + &
                        z**2*(60350.59709326911_dp)*Hr3(-1,-1,0) + z*(72876.67116734319_dp)*Hr3(-1,-1,0) + &
                        (16415.71090533781_dp)*Hr3(-1,0,-1) + z*(32831.42181067562_dp)*Hr3(-1,0,-1) + &
                        z**2*(32831.42181067562_dp)*Hr3(-1,0,-1) + nf*z*(-4796.098765432099_dp)*Hr3(-1,0,0) + &
                        nf*z**2*(-3732.8888888888887_dp)*Hr3(-1,0,0) + (-3469.8287521840994_dp)*Hr3(-1,0,0) + &
                        nf*(-1924.3456790123457_dp)*Hr3(-1,0,0) + (nf*(-337.9753086419753_dp)*Hr3(-1,0,0))/z + &
                        nf**2*(8._dp)*Hr3(-1,0,0) + nf**2*z*(16._dp)*Hr3(-1,0,0) + nf**2*z**2*(16._dp)*Hr3(-1,0,0) + &
                        ((10889.481481481482_dp)*Hr3(-1,0,0))/z + z*(68553.82397711328_dp)*Hr3(-1,0,0) + &
                        z**2*(77452.78694007624_dp)*Hr3(-1,0,0) + (-4024.980320290796_dp)*Hr3(-1,0,1) + &
                        nf*z*(-1252.7407407407406_dp)*Hr3(-1,0,1) + nf*z**2*(-1035.5555555555557_dp)*Hr3(-1,0,1) + &
                        nf*(-537.1851851851852_dp)*Hr3(-1,0,1) + (nf*(-88.88888888888889_dp)*Hr3(-1,0,1))/z + &
                        ((1530.6666666666667_dp)*Hr3(-1,0,1))/z + z*(18379.965285344333_dp)*Hr3(-1,0,1) + &
                        z**2*(30529.15047052952_dp)*Hr3(-1,0,1) + (8065.294500208726_dp)*Hr3(0,-1,-1) + &
                        z**2*(33407.514274976245_dp)*Hr3(0,-1,-1) + z*(53733.53807463207_dp)*Hr3(0,-1,-1) + &
                        z*(-94342.17542910855_dp)*Hr3(0,-1,0) + ((-4704._dp)*Hr3(0,-1,0))/z + &
                        nf*z**2*(-2653.3333333333335_dp)*Hr3(0,-1,0) + nf*(-2192.641975308642_dp)*Hr3(0,-1,0) + &
                        nf*z*(-1468.6913580246915_dp)*Hr3(0,-1,0) + (nf*(-3.5555555555555554_dp)*Hr3(0,-1,0))/z + &
                        nf**2*(5.333333333333333_dp)*Hr3(0,-1,0) + nf**2*z*(10.666666666666666_dp)*Hr3(0,-1,0) + &
                        nf**2*z**2*(10.666666666666666_dp)*Hr3(0,-1,0) + (20294.93034972544_dp)*Hr3(0,-1,0) + &
                        z**2*(81466.64364801004_dp)*Hr3(0,-1,0) + z*(-55226.407125522775_dp)*Hr3(0,0,-1) + &
                        z**2*(-7205.542294602719_dp)*Hr3(0,0,-1) + nf*(27.781108684547824_dp)*Hr3(0,0,-1) + &
                        nf*z*(55.56221736909565_dp)*Hr3(0,0,-1) + (27924.64441275132_dp)*Hr3(0,0,-1) + &
                        z**2*(-156236.68535606167_dp)*Hr3(0,0,0) + nf*z*(-5544.1465858326765_dp)*Hr3(0,0,0) + &
                        nf*(-2339.057583431956_dp)*Hr3(0,0,0) + nf**2*(-59.25925925925926_dp)*Hr3(0,0,0) + &
                        nf**2*z**2*(14.222222222222221_dp)*Hr3(0,0,0) + nf**2*z*(132.74074074074073_dp)*Hr3(0,0,0) + &
                        nf*z**2*(392.98765432098764_dp)*Hr3(0,0,0) + (18200.492920971323_dp)*Hr3(0,0,0) + &
                        z*(30058.275787768616_dp)*Hr3(0,0,0) + z**2*(-133545.75465842115_dp)*Hr3(0,0,1) + &
                        (-13255.98499739277_dp)*Hr3(0,0,1) + ((-2496._dp)*Hr3(0,0,1))/z + nf*z*(-1431.3975482019202_dp)*Hr3(0,0,1) + &
                        nf**2*z*(-39.111111111111114_dp)*Hr3(0,0,1) + nf**2*z**2*(-18.962962962962962_dp)*Hr3(0,0,1) + &
                        nf**2*(26.666666666666668_dp)*Hr3(0,0,1) + nf*(1151.8098852120713_dp)*Hr3(0,0,1) + &
                        nf*z**2*(2188.4444444444443_dp)*Hr3(0,0,1) + z*(13669.445067555514_dp)*Hr3(0,0,1) + &
                        z*(-75343.57489887223_dp)*Hr3(0,1,0) + z**2*(-41127.3105862651_dp)*Hr3(0,1,0) + &
                        (-25303.218848530516_dp)*Hr3(0,1,0) + ((-7493.333333333333_dp)*Hr3(0,1,0))/z + &
                        nf**2*z*(-39.111111111111114_dp)*Hr3(0,1,0) + nf**2*z**2*(-28.444444444444443_dp)*Hr3(0,1,0) + &
                        nf**2*(-3.5555555555555554_dp)*Hr3(0,1,0) + (nf*(124.44444444444444_dp)*Hr3(0,1,0))/z + &
                        nf*z*(897.283950617284_dp)*Hr3(0,1,0) + nf*(1203.0617283950617_dp)*Hr3(0,1,0) + &
                        nf*z**2*(1538.3703703703704_dp)*Hr3(0,1,0) + z*(-28333.63831545777_dp)*Hr3(0,1,1) + &
                        (-9024.741836714544_dp)*Hr3(0,1,1) + ((-4280._dp)*Hr3(0,1,1))/z + nf*z**2*(-370.3703703703704_dp)*Hr3(0,1,1) &
                        + nf*(-141.7283950617284_dp)*Hr3(0,1,1) + nf**2*z*(-4.7407407407407405_dp)*Hr3(0,1,1) + &
                        nf**2*(2.3703703703703702_dp)*Hr3(0,1,1) + nf**2*z**2*(4.7407407407407405_dp)*Hr3(0,1,1) + &
                        (nf*(32._dp)*Hr3(0,1,1))/z + nf*z*(226.5679012345679_dp)*Hr3(0,1,1) + &
                        z**2*(28430.514907055334_dp)*Hr3(0,1,1) + z**2*(-11869.844226376801_dp)*Hr3(1,0,-1) + &
                        (-5934.9221131884005_dp)*Hr3(1,0,-1) + z*(11869.844226376801_dp)*Hr3(1,0,-1) + &
                        z**2*(-17957.0864294809_dp)*Hr3(1,0,0) + nf*z*(-2664.320987654321_dp)*Hr3(1,0,0) + &
                        (-1085.3024739997072_dp)*Hr3(1,0,0) + (nf*(-96.5925925925926_dp)*Hr3(1,0,0))/z + &
                        nf**2*z**2*(-19.555555555555557_dp)*Hr3(1,0,0) + nf**2*(-9.777777777777779_dp)*Hr3(1,0,0) + &
                        nf**2*z*(19.555555555555557_dp)*Hr3(1,0,0) + ((638.8148148148148_dp)*Hr3(1,0,0))/z + &
                        nf*(1085.3086419753085_dp)*Hr3(1,0,0) + nf*z**2*(2137.432098765432_dp)*Hr3(1,0,0) + &
                        z*(21103.7901331846_dp)*Hr3(1,0,0) + z*(-46060.690131258176_dp)*Hr3(1,0,1) + &
                        ((-11262.222222222223_dp)*Hr3(1,0,1))/z + nf*z**2*(-1303.3086419753085_dp)*Hr3(1,0,1) + &
                        nf*(-163.25925925925927_dp)*Hr3(1,0,1) + nf**2*z*(-10.666666666666666_dp)*Hr3(1,0,1) + &
                        nf**2*(5.333333333333333_dp)*Hr3(1,0,1) + nf**2*z**2*(10.666666666666666_dp)*Hr3(1,0,1) + &
                        (nf*(215.30864197530863_dp)*Hr3(1,0,1))/z + nf*z*(1265.1851851851852_dp)*Hr3(1,0,1) + &
                        (1889.7524730364958_dp)*Hr3(1,0,1) + z**2*(48116.83827940632_dp)*Hr3(1,0,1) + &
                        z*(-68944.87385005968_dp)*Hr3(1,1,0) + ((-12819.555555555555_dp)*Hr3(1,1,0))/z + &
                        nf*z**2*(-971.7530864197531_dp)*Hr3(1,1,0) + nf*(36.148148148148145_dp)*Hr3(1,1,0) + &
                        (nf*(218.8641975308642_dp)*Hr3(1,1,0))/z + nf*z*(858.6666666666666_dp)*Hr3(1,1,0) + &
                        (10569.32581391873_dp)*Hr3(1,1,0) + z**2*(70847.3923685782_dp)*Hr3(1,1,0) + &
                        z*(-1895.6878785148544_dp)*Hr3(1,1,1) + nf*z**2*(-411.358024691358_dp)*Hr3(1,1,1) + &
                        nf*(-216.2962962962963_dp)*Hr3(1,1,1) + (nf*(-25.679012345679013_dp)*Hr3(1,1,1))/z + &
                        nf**2*z**2*(-5.925925925925926_dp)*Hr3(1,1,1) + nf**2*(-2.962962962962963_dp)*Hr3(1,1,1) + &
                        nf**2*z*(5.925925925925926_dp)*Hr3(1,1,1) + nf*z*(352.5925925925926_dp)*Hr3(1,1,1) + &
                        ((1945.7777777777778_dp)*Hr3(1,1,1))/z + z**2*(4218.058248885224_dp)*Hr3(1,1,1) + &
                        (4287.1031985166865_dp)*Hr3(1,1,1) + ((-668.4444444444445_dp)*Hr4(-1,-1,-1,0))/z + &
                        nf*z*(-237.03703703703704_dp)*Hr4(-1,-1,-1,0) + nf*z**2*(-237.03703703703704_dp)*Hr4(-1,-1,-1,0) + &
                        nf*(-118.51851851851852_dp)*Hr4(-1,-1,-1,0) + z**2*(850.9629629629629_dp)*Hr4(-1,-1,-1,0) + &
                        z*(5730.37037037037_dp)*Hr4(-1,-1,-1,0) + (6660.740740740741_dp)*Hr4(-1,-1,-1,0) + &
                        z**2*(-75110.22222222222_dp)*Hr4(-1,-1,0,0) + z*(-69742.22222222222_dp)*Hr4(-1,-1,0,0) + &
                        (-9617.333333333334_dp)*Hr4(-1,-1,0,0) + ((-7253.333333333333_dp)*Hr4(-1,-1,0,0))/z + &
                        nf*(540.1481481481482_dp)*Hr4(-1,-1,0,0) + nf*z*(1080.2962962962963_dp)*Hr4(-1,-1,0,0) + &
                        nf*z**2*(1080.2962962962963_dp)*Hr4(-1,-1,0,0) + z*(-47770.07407407407_dp)*Hr4(-1,-1,0,1) + &
                        z**2*(-47342.81481481482_dp)*Hr4(-1,-1,0,1) + (-10357.037037037036_dp)*Hr4(-1,-1,0,1) + &
                        ((-4480._dp)*Hr4(-1,-1,0,1))/z + nf*(372.74074074074076_dp)*Hr4(-1,-1,0,1) + &
                        nf*z*(745.4814814814815_dp)*Hr4(-1,-1,0,1) + nf*z**2*(745.4814814814815_dp)*Hr4(-1,-1,0,1) + &
                        z**2*(-26465.48148148148_dp)*Hr4(-1,0,-1,0) + z*(-24523.25925925926_dp)*Hr4(-1,0,-1,0) + &
                        (-3448._dp)*Hr4(-1,0,-1,0) + ((-2485.3333333333335_dp)*Hr4(-1,0,-1,0))/z + &
                        nf*(188.44444444444446_dp)*Hr4(-1,0,-1,0) + nf*z*(376.8888888888889_dp)*Hr4(-1,0,-1,0) + &
                        nf*z**2*(376.8888888888889_dp)*Hr4(-1,0,-1,0) + nf*z*(-713.4814814814815_dp)*Hr4(-1,0,0,0) + &
                        nf*z**2*(-713.4814814814815_dp)*Hr4(-1,0,0,0) + nf*(-356.74074074074076_dp)*Hr4(-1,0,0,0) + &
                        (4845.185185185185_dp)*Hr4(-1,0,0,0) + ((6561.777777777777_dp)*Hr4(-1,0,0,0))/z + &
                        z*(56070.81481481482_dp)*Hr4(-1,0,0,0) + z**2*(62542.07407407407_dp)*Hr4(-1,0,0,0) + &
                        nf*z*(-735.4074074074074_dp)*Hr4(-1,0,0,1) + nf*z**2*(-735.4074074074074_dp)*Hr4(-1,0,0,1) + &
                        nf*(-367.7037037037037_dp)*Hr4(-1,0,0,1) + ((6400._dp)*Hr4(-1,0,0,1))/z + &
                        (10274.814814814816_dp)*Hr4(-1,0,0,1) + z*(59804.444444444445_dp)*Hr4(-1,0,0,1) + &
                        z**2*(60894.51851851852_dp)*Hr4(-1,0,0,1) + nf*z*(-118.51851851851852_dp)*Hr4(-1,0,1,0) + &
                        nf*z**2*(-118.51851851851852_dp)*Hr4(-1,0,1,0) + nf*(-59.25925925925926_dp)*Hr4(-1,0,1,0) + &
                        ((2481.777777777778_dp)*Hr4(-1,0,1,0))/z + (4222.814814814815_dp)*Hr4(-1,0,1,0) + &
                        z**2*(20263.11111111111_dp)*Hr4(-1,0,1,0) + z*(21275.25925925926_dp)*Hr4(-1,0,1,0) + &
                        nf*(53.333333333333336_dp)*Hr4(-1,0,1,1) + nf*z*(106.66666666666667_dp)*Hr4(-1,0,1,1) + &
                        nf*z**2*(106.66666666666667_dp)*Hr4(-1,0,1,1) + ((1024._dp)*Hr4(-1,0,1,1))/z + &
                        z**2*(4472.888888888889_dp)*Hr4(-1,0,1,1) + (8445.037037037036_dp)*Hr4(-1,0,1,1) + &
                        z*(12613.925925925925_dp)*Hr4(-1,0,1,1) + z**2*(-50522.07407407407_dp)*Hr4(0,-1,-1,0) + &
                        z*(-20705.48148148148_dp)*Hr4(0,-1,-1,0) + ((-1152._dp)*Hr4(0,-1,-1,0))/z + &
                        nf*(295.4074074074074_dp)*Hr4(0,-1,-1,0) + nf*z*(495.4074074074074_dp)*Hr4(0,-1,-1,0) + &
                        nf*z**2*(502.51851851851853_dp)*Hr4(0,-1,-1,0) + (1034.0740740740741_dp)*Hr4(0,-1,-1,0) + &
                        (-4201.62962962963_dp)*Hr4(0,-1,0,0) + nf*z*(-1595.2592592592594_dp)*Hr4(0,-1,0,0) + &
                        nf*(-480.14814814814815_dp)*Hr4(0,-1,0,0) + nf*z**2*(-354.3703703703704_dp)*Hr4(0,-1,0,0) + &
                        ((4032._dp)*Hr4(0,-1,0,0))/z + z**2*(69601.77777777778_dp)*Hr4(0,-1,0,0) + &
                        z*(79678.51851851853_dp)*Hr4(0,-1,0,0) + nf*z**2*(-670.8148148148148_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(-666.074074074074_dp)*Hr4(0,-1,0,1) + (-370.0740740740741_dp)*Hr4(0,-1,0,1) + &
                        nf*(-322.3703703703704_dp)*Hr4(0,-1,0,1) + ((2304._dp)*Hr4(0,-1,0,1))/z + &
                        z*(39807.11111111111_dp)*Hr4(0,-1,0,1) + z**2*(63422.81481481482_dp)*Hr4(0,-1,0,1) + &
                        (-6693.925925925926_dp)*Hr4(0,0,-1,0) + nf*z*(-1540.148148148148_dp)*Hr4(0,0,-1,0) + &
                        nf*(-530.2222222222222_dp)*Hr4(0,0,-1,0) + nf*z**2*(525.0370370370371_dp)*Hr4(0,0,-1,0) + &
                        z**2*(15389.037037037036_dp)*Hr4(0,0,-1,0) + z*(25501.037037037036_dp)*Hr4(0,0,-1,0) + &
                        z**2*(-19877.037037037036_dp)*Hr4(0,0,0,0) + nf*z*(-5120.592592592592_dp)*Hr4(0,0,0,0) + &
                        (-3846.962962962963_dp)*Hr4(0,0,0,0) + nf*z**2*(-218.07407407407408_dp)*Hr4(0,0,0,0) + &
                        nf*(-30.51851851851852_dp)*Hr4(0,0,0,0) + nf**2*(-21.333333333333332_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(42.666666666666664_dp)*Hr4(0,0,0,0) + z*(52043.851851851854_dp)*Hr4(0,0,0,0) + &
                        z**2*(-47471.7037037037_dp)*Hr4(0,0,0,1) + nf*z*(-2112.8888888888887_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(-49.18518518518518_dp)*Hr4(0,0,0,1) + nf**2*z*(-28.444444444444443_dp)*Hr4(0,0,0,1) + &
                        nf**2*(14.222222222222221_dp)*Hr4(0,0,0,1) + nf*(551.0370370370371_dp)*Hr4(0,0,0,1) + &
                        (3290.074074074074_dp)*Hr4(0,0,0,1) + z*(31200._dp)*Hr4(0,0,0,1) + z**2*(-29600._dp)*Hr4(0,0,1,0) + &
                        (-65.03703703703704_dp)*Hr4(0,0,1,0) + nf**2*z*(-14.222222222222221_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(-9.481481481481481_dp)*Hr4(0,0,1,0) + nf**2*(7.111111111111111_dp)*Hr4(0,0,1,0) + &
                        nf*z*(146.96296296296296_dp)*Hr4(0,0,1,0) + nf*(197.33333333333334_dp)*Hr4(0,0,1,0) + &
                        ((576._dp)*Hr4(0,0,1,0))/z + z*(3603.1111111111113_dp)*Hr4(0,0,1,0) + &
                        z**2*(-25865.777777777777_dp)*Hr4(0,0,1,1) + (-5755.851851851852_dp)*Hr4(0,0,1,1) + &
                        z*(-4747.555555555556_dp)*Hr4(0,0,1,1) + nf*(-9.481481481481481_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(71.11111111111111_dp)*Hr4(0,0,1,1) + nf*z*(530.9629629629629_dp)*Hr4(0,0,1,1) + &
                        ((864._dp)*Hr4(0,0,1,1))/z + z**2*(-19727.40740740741_dp)*Hr4(0,1,0,0) + &
                        nf*z*(-656.2962962962963_dp)*Hr4(0,1,0,0) + nf*z**2*(21.925925925925927_dp)*Hr4(0,1,0,0) + &
                        (26.074074074074073_dp)*Hr4(0,1,0,0) + nf*(179.1851851851852_dp)*Hr4(0,1,0,0) + ((1408._dp)*Hr4(0,1,0,0))/z &
                        + z*(27118.814814814814_dp)*Hr4(0,1,0,0) + z*(-22421.925925925927_dp)*Hr4(0,1,0,1) + &
                        z**2*(-8734.222222222223_dp)*Hr4(0,1,0,1) + (-4218.37037037037_dp)*Hr4(0,1,0,1) + &
                        ((-768._dp)*Hr4(0,1,0,1))/z + nf*z**2*(130.37037037037038_dp)*Hr4(0,1,0,1) + &
                        nf*(137.4814814814815_dp)*Hr4(0,1,0,1) + nf*z*(862.8148148148148_dp)*Hr4(0,1,0,1) + &
                        z*(-25355.25925925926_dp)*Hr4(0,1,1,0) + z**2*(-11413.333333333334_dp)*Hr4(0,1,1,0) + &
                        (-4848.592592592592_dp)*Hr4(0,1,1,0) + ((-1088._dp)*Hr4(0,1,1,0))/z + nf*(177.1851851851852_dp)*Hr4(0,1,1,0) &
                        + nf*z**2*(267.85185185185185_dp)*Hr4(0,1,1,0) + nf*z*(815.4074074074074_dp)*Hr4(0,1,1,0) + &
                        z*(-15758.222222222223_dp)*Hr4(0,1,1,1) + ((-2208._dp)*Hr4(0,1,1,1))/z + &
                        nf*z**2*(-343.7037037037037_dp)*Hr4(0,1,1,1) + nf*(-168.88888888888889_dp)*Hr4(0,1,1,1) + &
                        nf*z*(231.11111111111111_dp)*Hr4(0,1,1,1) + (1520._dp)*Hr4(0,1,1,1) + &
                        z**2*(26761.777777777777_dp)*Hr4(0,1,1,1) + z*(-15633.777777777777_dp)*Hr4(1,0,-1,0) + &
                        ((-2321.777777777778_dp)*Hr4(1,0,-1,0))/z + (-944.8888888888889_dp)*Hr4(1,0,-1,0) + &
                        nf*z**2*(-101.92592592592592_dp)*Hr4(1,0,-1,0) + nf*(-50.96296296296296_dp)*Hr4(1,0,-1,0) + &
                        nf*z*(101.92592592592592_dp)*Hr4(1,0,-1,0) + z**2*(19476.444444444445_dp)*Hr4(1,0,-1,0) + &
                        z**2*(-37174.96296296296_dp)*Hr4(1,0,0,0) + nf*z*(-581.925925925926_dp)*Hr4(1,0,0,0) + &
                        nf*(290.962962962963_dp)*Hr4(1,0,0,0) + nf*z**2*(581.925925925926_dp)*Hr4(1,0,0,0) + &
                        (697.7777777777778_dp)*Hr4(1,0,0,0) + ((3489.777777777778_dp)*Hr4(1,0,0,0))/z + &
                        z*(29309.185185185186_dp)*Hr4(1,0,0,0) + z**2*(-13690.074074074075_dp)*Hr4(1,0,0,1) + &
                        nf*z*(-260.74074074074076_dp)*Hr4(1,0,0,1) + nf*(130.37037037037038_dp)*Hr4(1,0,0,1) + &
                        nf*z**2*(260.74074074074076_dp)*Hr4(1,0,0,1) + ((917.3333333333334_dp)*Hr4(1,0,0,1))/z + &
                        (2730.5185185185187_dp)*Hr4(1,0,0,1) + z*(8623.555555555555_dp)*Hr4(1,0,0,1) + &
                        z*(-11691.25925925926_dp)*Hr4(1,0,1,0) + ((-1322.6666666666667_dp)*Hr4(1,0,1,0))/z + &
                        nf*z**2*(-280.8888888888889_dp)*Hr4(1,0,1,0) + nf*(-140.44444444444446_dp)*Hr4(1,0,1,0) + &
                        nf*z*(280.8888888888889_dp)*Hr4(1,0,1,0) + (3816.5925925925926_dp)*Hr4(1,0,1,0) + &
                        z**2*(11632._dp)*Hr4(1,0,1,0) + z*(-27093.333333333332_dp)*Hr4(1,0,1,1) + &
                        ((-3018.6666666666665_dp)*Hr4(1,0,1,1))/z + nf*z**2*(-512._dp)*Hr4(1,0,1,1) + nf*(-256._dp)*Hr4(1,0,1,1) + &
                        nf*z*(512._dp)*Hr4(1,0,1,1) + (4756.444444444444_dp)*Hr4(1,0,1,1) + &
                        z**2*(28966.222222222223_dp)*Hr4(1,0,1,1) + z**2*(-6977.185185185185_dp)*Hr4(1,1,0,0) + &
                        nf*z*(-247.7037037037037_dp)*Hr4(1,1,0,0) + ((85.33333333333333_dp)*Hr4(1,1,0,0))/z + &
                        nf*(123.85185185185185_dp)*Hr4(1,1,0,0) + nf*z**2*(247.7037037037037_dp)*Hr4(1,1,0,0) + &
                        (2265.037037037037_dp)*Hr4(1,1,0,0) + z*(3262.3703703703704_dp)*Hr4(1,1,0,0) + &
                        z*(-28641.777777777777_dp)*Hr4(1,1,0,1) + ((-3658.6666666666665_dp)*Hr4(1,1,0,1))/z + &
                        nf*z**2*(-408.8888888888889_dp)*Hr4(1,1,0,1) + nf*(-204.44444444444446_dp)*Hr4(1,1,0,1) + &
                        nf*z*(408.8888888888889_dp)*Hr4(1,1,0,1) + (3255.1111111111113_dp)*Hr4(1,1,0,1) + &
                        z**2*(32132.444444444445_dp)*Hr4(1,1,0,1) + z*(-28850.962962962964_dp)*Hr4(1,1,1,0) + &
                        ((-4156.444444444444_dp)*Hr4(1,1,1,0))/z + nf*z**2*(-282.0740740740741_dp)*Hr4(1,1,1,0) + &
                        nf*(-141.03703703703704_dp)*Hr4(1,1,1,0) + nf*z*(282.0740740740741_dp)*Hr4(1,1,1,0) + &
                        (1641.1851851851852_dp)*Hr4(1,1,1,0) + z**2*(33805.333333333336_dp)*Hr4(1,1,1,0) + &
                        z**2*(-11795.555555555555_dp)*Hr4(1,1,1,1) + (-4183.7037037037035_dp)*Hr4(1,1,1,1) + &
                        nf*z*(-533.3333333333334_dp)*Hr4(1,1,1,1) + nf*(266.6666666666667_dp)*Hr4(1,1,1,1) + &
                        ((533.3333333333334_dp)*Hr4(1,1,1,1))/z + nf*z**2*(533.3333333333334_dp)*Hr4(1,1,1,1) + &
                        z*(11597.037037037036_dp)*Hr4(1,1,1,1) + (2609.777777777778_dp)*Hr5(-1,-1,-1,-1,0) + &
                        z*(5219.555555555556_dp)*Hr5(-1,-1,-1,-1,0) + z**2*(5219.555555555556_dp)*Hr5(-1,-1,-1,-1,0) + &
                        (12951.111111111111_dp)*Hr5(-1,-1,-1,0,0) + z*(25902.222222222223_dp)*Hr5(-1,-1,-1,0,0) + &
                        z**2*(25902.222222222223_dp)*Hr5(-1,-1,-1,0,0) + (9795.555555555555_dp)*Hr5(-1,-1,-1,0,1) + &
                        z*(19591.11111111111_dp)*Hr5(-1,-1,-1,0,1) + z**2*(19591.11111111111_dp)*Hr5(-1,-1,-1,0,1) + &
                        (3476.740740740741_dp)*Hr5(-1,-1,0,-1,0) + z*(6953.481481481482_dp)*Hr5(-1,-1,0,-1,0) + &
                        z**2*(6953.481481481482_dp)*Hr5(-1,-1,0,-1,0) + z*(-24517.333333333332_dp)*Hr5(-1,-1,0,0,0) + &
                        z**2*(-24517.333333333332_dp)*Hr5(-1,-1,0,0,0) + (-12258.666666666666_dp)*Hr5(-1,-1,0,0,0) + &
                        z*(-29161.48148148148_dp)*Hr5(-1,-1,0,0,1) + z**2*(-29161.48148148148_dp)*Hr5(-1,-1,0,0,1) + &
                        (-14580.74074074074_dp)*Hr5(-1,-1,0,0,1) + z*(-13089.185185185184_dp)*Hr5(-1,-1,0,1,0) + &
                        z**2*(-13089.185185185184_dp)*Hr5(-1,-1,0,1,0) + (-6544.592592592592_dp)*Hr5(-1,-1,0,1,0) + &
                        z*(-12804.74074074074_dp)*Hr5(-1,-1,0,1,1) + z**2*(-12804.74074074074_dp)*Hr5(-1,-1,0,1,1) + &
                        (-6402.37037037037_dp)*Hr5(-1,-1,0,1,1) + (3476.740740740741_dp)*Hr5(-1,0,-1,-1,0) + &
                        z*(6953.481481481482_dp)*Hr5(-1,0,-1,-1,0) + z**2*(6953.481481481482_dp)*Hr5(-1,0,-1,-1,0) + &
                        z*(-21076.14814814815_dp)*Hr5(-1,0,-1,0,0) + z**2*(-21076.14814814815_dp)*Hr5(-1,0,-1,0,0) + &
                        (-10538.074074074075_dp)*Hr5(-1,0,-1,0,0) + z*(-16482.37037037037_dp)*Hr5(-1,0,-1,0,1) + &
                        z**2*(-16482.37037037037_dp)*Hr5(-1,0,-1,0,1) + (-8241.185185185184_dp)*Hr5(-1,0,-1,0,1) + &
                        z*(-6205.62962962963_dp)*Hr5(-1,0,0,-1,0) + z**2*(-6205.62962962963_dp)*Hr5(-1,0,0,-1,0) + &
                        (-3102.814814814815_dp)*Hr5(-1,0,0,-1,0) + (4664.592592592592_dp)*Hr5(-1,0,0,0,0) + &
                        z*(9329.185185185184_dp)*Hr5(-1,0,0,0,0) + z**2*(9329.185185185184_dp)*Hr5(-1,0,0,0,0) + &
                        (8865.481481481482_dp)*Hr5(-1,0,0,0,1) + z*(17730.962962962964_dp)*Hr5(-1,0,0,0,1) + &
                        z**2*(17730.962962962964_dp)*Hr5(-1,0,0,0,1) + (5507.555555555556_dp)*Hr5(-1,0,0,1,0) + &
                        z*(11015.111111111111_dp)*Hr5(-1,0,0,1,0) + z**2*(11015.111111111111_dp)*Hr5(-1,0,0,1,0) + &
                        (5766.518518518518_dp)*Hr5(-1,0,0,1,1) + z*(11533.037037037036_dp)*Hr5(-1,0,0,1,1) + &
                        z**2*(11533.037037037036_dp)*Hr5(-1,0,0,1,1) + (1147.851851851852_dp)*Hr5(-1,0,1,0,0) + &
                        z*(2295.703703703704_dp)*Hr5(-1,0,1,0,0) + z**2*(2295.703703703704_dp)*Hr5(-1,0,1,0,0) + &
                        (865.1851851851852_dp)*Hr5(-1,0,1,0,1) + z*(1730.3703703703704_dp)*Hr5(-1,0,1,0,1) + &
                        z**2*(1730.3703703703704_dp)*Hr5(-1,0,1,0,1) + (865.1851851851852_dp)*Hr5(-1,0,1,1,0) + &
                        z*(1730.3703703703704_dp)*Hr5(-1,0,1,1,0) + z**2*(1730.3703703703704_dp)*Hr5(-1,0,1,1,0) + &
                        z*(-4160._dp)*Hr5(-1,0,1,1,1) + z**2*(-4160._dp)*Hr5(-1,0,1,1,1) + (-2080._dp)*Hr5(-1,0,1,1,1) + &
                        z*(2088.296296296296_dp)*Hr5(0,-1,-1,-1,0) + (6424.888888888889_dp)*Hr5(0,-1,-1,-1,0) + &
                        z**2*(10055.111111111111_dp)*Hr5(0,-1,-1,-1,0) + z*(-50277.333333333336_dp)*Hr5(0,-1,-1,0,0) + &
                        z**2*(-15342.222222222223_dp)*Hr5(0,-1,-1,0,0) + (3295.1111111111113_dp)*Hr5(0,-1,-1,0,0) + &
                        z*(-31621.925925925927_dp)*Hr5(0,-1,-1,0,1) + z**2*(-15281.777777777777_dp)*Hr5(0,-1,-1,0,1) + &
                        (-1690.6666666666667_dp)*Hr5(0,-1,-1,0,1) + z*(-18513.777777777777_dp)*Hr5(0,-1,0,-1,0) + &
                        z**2*(-6132.148148148148_dp)*Hr5(0,-1,0,-1,0) + (740.7407407407408_dp)*Hr5(0,-1,0,-1,0) + &
                        (-5198.222222222223_dp)*Hr5(0,-1,0,0,0) + z**2*(9397.333333333334_dp)*Hr5(0,-1,0,0,0) + &
                        z*(42619.851851851854_dp)*Hr5(0,-1,0,0,0) + (86.81481481481481_dp)*Hr5(0,-1,0,0,1) + &
                        z**2*(17022.814814814814_dp)*Hr5(0,-1,0,0,1) + z*(42027.25925925926_dp)*Hr5(0,-1,0,0,1) + &
                        (1448.2962962962963_dp)*Hr5(0,-1,0,1,0) + z**2*(8139.851851851852_dp)*Hr5(0,-1,0,1,0) + &
                        z*(15525.925925925925_dp)*Hr5(0,-1,0,1,0) + z*(6637.037037037037_dp)*Hr5(0,-1,0,1,1) + &
                        (8417.185185185184_dp)*Hr5(0,-1,0,1,1) + z**2*(14724.74074074074_dp)*Hr5(0,-1,0,1,1) + &
                        z*(-15723.25925925926_dp)*Hr5(0,0,-1,-1,0) + z**2*(-2863.4074074074074_dp)*Hr5(0,0,-1,-1,0) + &
                        nf*(33.77777777777778_dp)*Hr5(0,0,-1,-1,0) + nf*z*(67.55555555555556_dp)*Hr5(0,0,-1,-1,0) + &
                        (11743.703703703704_dp)*Hr5(0,0,-1,-1,0) + (-8494.222222222223_dp)*Hr5(0,0,-1,0,0) + &
                        nf*z*(-545.7777777777778_dp)*Hr5(0,0,-1,0,0) + nf*(-272.8888888888889_dp)*Hr5(0,0,-1,0,0) + &
                        z**2*(3133.6296296296296_dp)*Hr5(0,0,-1,0,0) + z*(50641.18518518518_dp)*Hr5(0,0,-1,0,0) + &
                        (-11104.296296296296_dp)*Hr5(0,0,-1,0,1) + z**2*(2948.740740740741_dp)*Hr5(0,0,-1,0,1) + &
                        z*(25712._dp)*Hr5(0,0,-1,0,1) + nf*z*(-592._dp)*Hr5(0,0,0,-1,0) + nf*(-569.7777777777778_dp)*Hr5(0,0,0,-1,0) &
                        + z**2*(493.037037037037_dp)*Hr5(0,0,0,-1,0) + (6752.2962962962965_dp)*Hr5(0,0,0,-1,0) + &
                        z*(18136.88888888889_dp)*Hr5(0,0,0,-1,0) + nf*z*(-2088.8888888888887_dp)*Hr5(0,0,0,0,0) + &
                        z**2*(-834.3703703703703_dp)*Hr5(0,0,0,0,0) + nf*(-547.5555555555555_dp)*Hr5(0,0,0,0,0) + &
                        (3878.5185185185187_dp)*Hr5(0,0,0,0,0) + z*(36983.7037037037_dp)*Hr5(0,0,0,0,0) + &
                        (-6442.222222222223_dp)*Hr5(0,0,0,0,1) + z**2*(-3555.5555555555557_dp)*Hr5(0,0,0,0,1) + &
                        nf*z*(-1003.5555555555555_dp)*Hr5(0,0,0,0,1) + nf*(501.77777777777777_dp)*Hr5(0,0,0,0,1) + &
                        z*(41321.48148148148_dp)*Hr5(0,0,0,0,1) + z**2*(-6109.62962962963_dp)*Hr5(0,0,0,1,0) + &
                        (-833.6296296296297_dp)*Hr5(0,0,0,1,0) + nf*z*(-355.55555555555554_dp)*Hr5(0,0,0,1,0) + &
                        nf*(241.77777777777777_dp)*Hr5(0,0,0,1,0) + z*(28397.333333333332_dp)*Hr5(0,0,0,1,0) + &
                        z**2*(-9943.703703703704_dp)*Hr5(0,0,0,1,1) + nf*(-85.33333333333333_dp)*Hr5(0,0,0,1,1) + &
                        nf*z*(170.66666666666666_dp)*Hr5(0,0,0,1,1) + (5950.518518518518_dp)*Hr5(0,0,0,1,1) + &
                        z*(17237.925925925927_dp)*Hr5(0,0,0,1,1) + z**2*(-2697.4814814814813_dp)*Hr5(0,0,1,0,0) + &
                        nf*z*(-370.6666666666667_dp)*Hr5(0,0,1,0,0) + nf*(185.33333333333334_dp)*Hr5(0,0,1,0,0) + &
                        (694.9629629629629_dp)*Hr5(0,0,1,0,0) + z*(29736.88888888889_dp)*Hr5(0,0,1,0,0) + &
                        z**2*(-3166.814814814815_dp)*Hr5(0,0,1,0,1) + nf*(-32._dp)*Hr5(0,0,1,0,1) + nf*z*(64._dp)*Hr5(0,0,1,0,1) + &
                        (2337.185185185185_dp)*Hr5(0,0,1,0,1) + z*(2619.259259259259_dp)*Hr5(0,0,1,0,1) + &
                        z**2*(-3546.074074074074_dp)*Hr5(0,0,1,1,0) + z*(-699.2592592592592_dp)*Hr5(0,0,1,1,0) + &
                        nf*(-32._dp)*Hr5(0,0,1,1,0) + nf*z*(64._dp)*Hr5(0,0,1,1,0) + (2876.4444444444443_dp)*Hr5(0,0,1,1,0) + &
                        z*(-18706.962962962964_dp)*Hr5(0,0,1,1,1) + (-6662.518518518518_dp)*Hr5(0,0,1,1,1) + &
                        nf*z*(-35.55555555555556_dp)*Hr5(0,0,1,1,1) + nf*(17.77777777777778_dp)*Hr5(0,0,1,1,1) + &
                        z**2*(388.74074074074076_dp)*Hr5(0,0,1,1,1) + z*(-14151.111111111111_dp)*Hr5(0,1,0,-1,0) + &
                        (-3372.4444444444443_dp)*Hr5(0,1,0,-1,0) + z**2*(440.8888888888889_dp)*Hr5(0,1,0,-1,0) + &
                        z**2*(-2460.4444444444443_dp)*Hr5(0,1,0,0,0) + (5416.148148148148_dp)*Hr5(0,1,0,0,0) + &
                        z*(25980.74074074074_dp)*Hr5(0,1,0,0,0) + z**2*(-3477.3333333333335_dp)*Hr5(0,1,0,0,1) + &
                        (1709.7777777777778_dp)*Hr5(0,1,0,0,1) + z*(13697.481481481482_dp)*Hr5(0,1,0,0,1) + &
                        z**2*(-4096._dp)*Hr5(0,1,0,1,0) + (-2194.962962962963_dp)*Hr5(0,1,0,1,0) + &
                        z*(-602.074074074074_dp)*Hr5(0,1,0,1,0) + z*(-13890.37037037037_dp)*Hr5(0,1,0,1,1) + &
                        (-3182.814814814815_dp)*Hr5(0,1,0,1,1) + z**2*(52.148148148148145_dp)*Hr5(0,1,0,1,1) + &
                        z**2*(-3268.740740740741_dp)*Hr5(0,1,1,0,0) + (513.925925925926_dp)*Hr5(0,1,1,0,0) + &
                        z*(8601.185185185184_dp)*Hr5(0,1,1,0,0) + z*(-17232.59259259259_dp)*Hr5(0,1,1,0,1) + &
                        (-4391.7037037037035_dp)*Hr5(0,1,1,0,1) + z**2*(-943.4074074074074_dp)*Hr5(0,1,1,0,1) + &
                        z*(-19626.666666666668_dp)*Hr5(0,1,1,1,0) + (-5434.666666666667_dp)*Hr5(0,1,1,1,0) + &
                        z**2*(-2128.5925925925926_dp)*Hr5(0,1,1,1,0) + z*(-4370.962962962963_dp)*Hr5(0,1,1,1,1) + &
                        (4585.481481481482_dp)*Hr5(0,1,1,1,1) + z**2*(6836.148148148148_dp)*Hr5(0,1,1,1,1) + &
                        z**2*(-6000.592592592592_dp)*Hr5(1,0,-1,-1,0) + (-3000.296296296296_dp)*Hr5(1,0,-1,-1,0) + &
                        z*(6000.592592592592_dp)*Hr5(1,0,-1,-1,0) + z*(-5108.148148148148_dp)*Hr5(1,0,-1,0,0) + &
                        (2554.074074074074_dp)*Hr5(1,0,-1,0,0) + z**2*(5108.148148148148_dp)*Hr5(1,0,-1,0,0) + &
                        z*(-4215.7037037037035_dp)*Hr5(1,0,-1,0,1) + (2107.8518518518517_dp)*Hr5(1,0,-1,0,1) + &
                        z**2*(4215.7037037037035_dp)*Hr5(1,0,-1,0,1) + z*(-6703.407407407408_dp)*Hr5(1,0,0,-1,0) + &
                        (3351.703703703704_dp)*Hr5(1,0,0,-1,0) + z**2*(6703.407407407408_dp)*Hr5(1,0,0,-1,0) + &
                        z**2*(-7110.518518518518_dp)*Hr5(1,0,0,0,0) + (-3555.259259259259_dp)*Hr5(1,0,0,0,0) + &
                        z*(7110.518518518518_dp)*Hr5(1,0,0,0,0) + z**2*(-9739.25925925926_dp)*Hr5(1,0,0,0,1) + &
                        (-4869.62962962963_dp)*Hr5(1,0,0,0,1) + z*(9739.25925925926_dp)*Hr5(1,0,0,0,1) + &
                        z**2*(-4097.185185185185_dp)*Hr5(1,0,0,1,0) + (-2048.5925925925926_dp)*Hr5(1,0,0,1,0) + &
                        z*(4097.185185185185_dp)*Hr5(1,0,0,1,0) + z**2*(-746.6666666666666_dp)*Hr5(1,0,0,1,1) + &
                        (-373.3333333333333_dp)*Hr5(1,0,0,1,1) + z*(746.6666666666666_dp)*Hr5(1,0,0,1,1) + &
                        z**2*(-2257.777777777778_dp)*Hr5(1,0,1,0,0) + (-1128.888888888889_dp)*Hr5(1,0,1,0,0) + &
                        z*(2257.777777777778_dp)*Hr5(1,0,1,0,0) + z*(-6840.888888888889_dp)*Hr5(1,0,1,0,1) + &
                        (3420.4444444444443_dp)*Hr5(1,0,1,0,1) + z**2*(6840.888888888889_dp)*Hr5(1,0,1,0,1) + &
                        z*(-6994.962962962963_dp)*Hr5(1,0,1,1,0) + (3497.4814814814813_dp)*Hr5(1,0,1,1,0) + &
                        z**2*(6994.962962962963_dp)*Hr5(1,0,1,1,0) + z*(-10788.74074074074_dp)*Hr5(1,0,1,1,1) + &
                        (5394.37037037037_dp)*Hr5(1,0,1,1,1) + z**2*(10788.74074074074_dp)*Hr5(1,0,1,1,1) + &
                        z*(-7406.222222222223_dp)*Hr5(1,1,0,-1,0) + (3703.1111111111113_dp)*Hr5(1,1,0,-1,0) + &
                        z**2*(7406.222222222223_dp)*Hr5(1,1,0,-1,0) + z**2*(-9495.111111111111_dp)*Hr5(1,1,0,0,0) + &
                        (-4747.555555555556_dp)*Hr5(1,1,0,0,0) + z*(9495.111111111111_dp)*Hr5(1,1,0,0,0) + &
                        z*(-1152._dp)*Hr5(1,1,0,0,1) + (576._dp)*Hr5(1,1,0,0,1) + z**2*(1152._dp)*Hr5(1,1,0,0,1) + &
                        z*(-7765.333333333333_dp)*Hr5(1,1,0,1,0) + (3882.6666666666665_dp)*Hr5(1,1,0,1,0) + &
                        z**2*(7765.333333333333_dp)*Hr5(1,1,0,1,0) + z*(-12329.481481481482_dp)*Hr5(1,1,0,1,1) + &
                        (6164.740740740741_dp)*Hr5(1,1,0,1,1) + z**2*(12329.481481481482_dp)*Hr5(1,1,0,1,1) + &
                        z*(-3856.5925925925926_dp)*Hr5(1,1,1,0,0) + (1928.2962962962963_dp)*Hr5(1,1,1,0,0) + &
                        z**2*(3856.5925925925926_dp)*Hr5(1,1,1,0,0) + z*(-13253.925925925925_dp)*Hr5(1,1,1,0,1) + &
                        (6626.962962962963_dp)*Hr5(1,1,1,0,1) + z**2*(13253.925925925925_dp)*Hr5(1,1,1,0,1) + &
                        z*(-13562.074074074075_dp)*Hr5(1,1,1,1,0) + (6781.037037037037_dp)*Hr5(1,1,1,1,0) + &
                        z**2*(13562.074074074075_dp)*Hr5(1,1,1,1,0) + z**2*(-13451.851851851852_dp)*Hr5(1,1,1,1,1) + &
                        (-6725.925925925926_dp)*Hr5(1,1,1,1,1) + z*(13451.851851851852_dp)*Hr5(1,1,1,1,1)

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
                        Ibar_select = z*(-25062.618087033014_dp) + (-12379.976043668541_dp)/z + nf*z**2*(-1702.9237784280365_dp) + &
                        nf*(-1234.8422064318481_dp) + nf*z*(-680.3242900146914_dp) + nf**2*(-66.57164594039419_dp) + &
                        nf**2*z**2*(-14.224807724575667_dp) + (nf**2*(-1.3607681755829903_dp))/z + nf**2*z*(90.35636581594531_dp) + &
                        (nf*(321.8485748241723_dp))/z + (34121.32903631323_dp) + z**2*(170449.15385993893_dp) + &
                        (z**2*(-33493.398499606374_dp) + z*(-29526.386149552054_dp) + (-4766.542901807625_dp) + &
                        (-2624.0963010106225_dp)/z + nf*((194.79268603960725_dp) + z*(389.5853720792145_dp) + &
                        z**2*(389.5853720792145_dp)))*Hr1(-1) + z**2*(-53590.06282238449_dp)*Hr1(0) + (-9947.657111773726_dp)*Hr1(0) &
                        + z*(-4048.001878535186_dp)*Hr1(0) + ((-3631.2645859245695_dp)*Hr1(0))/z + &
                        nf*z*(-1803.3764862671242_dp)*Hr1(0) + nf**2*(-43.14190132717933_dp)*Hr1(0) + &
                        nf**2*z**2*(-18.765432098765434_dp)*Hr1(0) + nf**2*z*(11.617135987692002_dp)*Hr1(0) + &
                        (nf*(53.925925925925924_dp)*Hr1(0))/z + nf*(111.08552037829526_dp)*Hr1(0) + &
                        nf*z**2*(3391.4742775240934_dp)*Hr1(0) + z*(-21175.711717849503_dp)*Hr1(1) + &
                        ((-6174.603786946309_dp)*Hr1(1))/z + nf*z*(-1535.4831094605377_dp)*Hr1(1) + &
                        nf**2*z*(-5.004115226337449_dp)*Hr1(1) + nf**2*z**2*(-4.609053497942387_dp)*Hr1(1) + &
                        (nf**2*(-1.0534979423868314_dp)*Hr1(1))/z + nf**2*(5.0699588477366255_dp)*Hr1(1) + &
                        (nf*(27.116598079561044_dp)*Hr1(1))/z + nf*(731.4452584339725_dp)*Hr1(1) + &
                        nf*z**2*(1073.7986101464087_dp)*Hr1(1) + (1664.1747010568315_dp)*Hr1(1) + z**2*(22156.96964892434_dp)*Hr1(1) &
                        + (4561.300628230202_dp)*Hr2(-1,-1) + z*(9122.601256460404_dp)*Hr2(-1,-1) + &
                        z**2*(9122.601256460404_dp)*Hr2(-1,-1) + nf*z*(-930.7654320987655_dp)*Hr2(-1,0) + &
                        nf*z**2*(-746.3703703703703_dp)*Hr2(-1,0) + nf*(-308.8395061728395_dp)*Hr2(-1,0) + &
                        (nf*(-81.77777777777777_dp)*Hr2(-1,0))/z + nf**2*(1.7777777777777777_dp)*Hr2(-1,0) + &
                        nf**2*z*(3.5555555555555554_dp)*Hr2(-1,0) + nf**2*z**2*(3.5555555555555554_dp)*Hr2(-1,0) + &
                        (1396.7872846928685_dp)*Hr2(-1,0) + ((3608.8888888888887_dp)*Hr2(-1,0))/z + &
                        z**2*(20343.84617432401_dp)*Hr2(-1,0) + z*(20938.26592741043_dp)*Hr2(-1,0) + &
                        z*(-16369.733982778005_dp)*Hr2(0,-1) + z**2*(-1250.3123534285382_dp)*Hr2(0,-1) + &
                        (3623.566363158799_dp)*Hr2(0,-1) + z**2*(-33015.56312630529_dp)*Hr2(0,0) + &
                        nf*z*(-1440.4571137848777_dp)*Hr2(0,0) + nf*(-434.1336431389341_dp)*Hr2(0,0) + ((-416._dp)*Hr2(0,0))/z + &
                        nf**2*(-15.802469135802468_dp)*Hr2(0,0) + nf**2*z**2*(4.7407407407407405_dp)*Hr2(0,0) + &
                        nf**2*z*(31.604938271604937_dp)*Hr2(0,0) + nf*z**2*(240.6255144032922_dp)*Hr2(0,0) + &
                        (4204.879973438621_dp)*Hr2(0,0) + z*(8609.005101326773_dp)*Hr2(0,0) + z**2*(-35283.17859149011_dp)*Hr2(0,1) &
                        + (-5031.1196328471_dp)*Hr2(0,1) + nf*z*(-910.156378600823_dp)*Hr2(0,1) + &
                        (nf*(-30.814814814814813_dp)*Hr2(0,1))/z + nf**2*z**2*(-4.7407407407407405_dp)*Hr2(0,1) + &
                        nf**2*z*(1.5802469135802468_dp)*Hr2(0,1) + nf**2*(2.765432098765432_dp)*Hr2(0,1) + &
                        ((247.11111111111111_dp)*Hr2(0,1))/z + nf*(422.2880658436214_dp)*Hr2(0,1) + &
                        nf*z**2*(813.925925925926_dp)*Hr2(0,1) + z*(10776.385430090992_dp)*Hr2(0,1) + &
                        z**2*(-20082.995325954176_dp)*Hr2(1,0) + (-3744.0655642116553_dp)*Hr2(1,0) + &
                        nf*z*(-922.9958847736625_dp)*Hr2(1,0) + (nf*(-58.205761316872426_dp)*Hr2(1,0))/z + &
                        nf**2*z**2*(-4.7407407407407405_dp)*Hr2(1,0) + nf**2*(-2.3703703703703702_dp)*Hr2(1,0) + &
                        nf**2*z*(4.7407407407407405_dp)*Hr2(1,0) + nf*(324.0164609053498_dp)*Hr2(1,0) + &
                        nf*z**2*(801.6460905349794_dp)*Hr2(1,0) + ((2211.1604938271603_dp)*Hr2(1,0))/z + &
                        z*(20533.06940002825_dp)*Hr2(1,0) + ((1.9753086419753085_dp)*(z**2*(-1683.0870458620461_dp) + (-1096.6_dp) + &
                        z*(-1076.381477068977_dp) + nf**2*z*(z*(-1._dp) + (0.5_dp) + z**2*(1._dp)) + &
                        nf*(z**2*(-58.53333333333334_dp) + (11.666666666666668_dp) + z*(50.41666666666667_dp) + &
                        z**3*(73.11666666666666_dp)) + z**3*(1376.8120458620465_dp))*Hr2(1,1))/z + &
                        z**2*(-13376.592592592593_dp)*Hr3(-1,-1,0) + z*(-11408.79012345679_dp)*Hr3(-1,-1,0) + &
                        ((-1313.1851851851852_dp)*Hr3(-1,-1,0))/z + (-697.3827160493827_dp)*Hr3(-1,-1,0) + &
                        nf*(98.17283950617283_dp)*Hr3(-1,-1,0) + nf*z*(196.34567901234567_dp)*Hr3(-1,-1,0) + &
                        nf*z**2*(196.34567901234567_dp)*Hr3(-1,-1,0) + nf*z*(-167.50617283950618_dp)*Hr3(-1,0,0) + &
                        nf*z**2*(-167.50617283950618_dp)*Hr3(-1,0,0) + nf*(-83.75308641975309_dp)*Hr3(-1,0,0) + &
                        (585.4814814814815_dp)*Hr3(-1,0,0) + ((1125.9259259259259_dp)*Hr3(-1,0,0))/z + &
                        z*(9751.703703703704_dp)*Hr3(-1,0,0) + z**2*(11449.481481481482_dp)*Hr3(-1,0,0) + &
                        nf*z*(-138.66666666666666_dp)*Hr3(-1,0,1) + nf*z**2*(-138.66666666666666_dp)*Hr3(-1,0,1) + &
                        nf*(-69.33333333333333_dp)*Hr3(-1,0,1) + (473.58024691358025_dp)*Hr3(-1,0,1) + &
                        ((938.6666666666666_dp)*Hr3(-1,0,1))/z + z*(8094.617283950617_dp)*Hr3(-1,0,1) + &
                        z**2*(9522.37037037037_dp)*Hr3(-1,0,1) + (-1164.148148148148_dp)*Hr3(0,-1,0) + &
                        nf*z*(-338.17283950617286_dp)*Hr3(0,-1,0) + nf*(-25.08641975308642_dp)*Hr3(0,-1,0) + &
                        nf*z**2*(48.98765432098765_dp)*Hr3(0,-1,0) + ((576._dp)*Hr3(0,-1,0))/z + &
                        z**2*(3125.925925925926_dp)*Hr3(0,-1,0) + z*(13699.555555555555_dp)*Hr3(0,-1,0) + &
                        z**2*(-3847.3086419753085_dp)*Hr3(0,0,0) + nf*z*(-1155.1604938271605_dp)*Hr3(0,0,0) + &
                        (-1148.148148148148_dp)*Hr3(0,0,0) + nf*z**2*(-18.962962962962962_dp)*Hr3(0,0,0) + &
                        nf**2*(-7.111111111111111_dp)*Hr3(0,0,0) + nf**2*z*(14.222222222222221_dp)*Hr3(0,0,0) + &
                        nf*(72.88888888888889_dp)*Hr3(0,0,0) + z*(15039.20987654321_dp)*Hr3(0,0,0) + &
                        z**2*(-9015.111111111111_dp)*Hr3(0,0,1) + nf*z*(-655.4074074074074_dp)*Hr3(0,0,1) + &
                        nf*z**2*(-15.802469135802468_dp)*Hr3(0,0,1) + nf**2*z*(-4.7407407407407405_dp)*Hr3(0,0,1) + &
                        nf**2*(2.3703703703703702_dp)*Hr3(0,0,1) + nf*(115.55555555555556_dp)*Hr3(0,0,1) + ((384._dp)*Hr3(0,0,1))/z &
                        + (1272._dp)*Hr3(0,0,1) + z*(16783.604938271605_dp)*Hr3(0,0,1) + z**2*(-9496.888888888889_dp)*Hr3(0,1,0) + &
                        nf*z*(-331.85185185185185_dp)*Hr3(0,1,0) + nf*(45.03703703703704_dp)*Hr3(0,1,0) + &
                        nf*z**2*(64.79012345679013_dp)*Hr3(0,1,0) + (485.9259259259259_dp)*Hr3(0,1,0) + &
                        ((938.6666666666666_dp)*Hr3(0,1,0))/z + z*(15373.432098765432_dp)*Hr3(0,1,0) + &
                        z**2*(-11507.555555555555_dp)*Hr3(0,1,1) + (-895.7037037037037_dp)*Hr3(0,1,1) + &
                        nf*z*(-4.7407407407407405_dp)*Hr3(0,1,1) + nf*(64.5925925925926_dp)*Hr3(0,1,1) + &
                        nf*z**2*(115.35802469135803_dp)*Hr3(0,1,1) + ((832._dp)*Hr3(0,1,1))/z + z*(4792.493827160494_dp)*Hr3(0,1,1) &
                        + z**2*(-9580.641975308641_dp)*Hr3(1,0,0) + (-628.2469135802469_dp)*Hr3(1,0,0) + &
                        nf*z*(-148.93827160493828_dp)*Hr3(1,0,0) + nf*(74.46913580246914_dp)*Hr3(1,0,0) + &
                        nf*z**2*(148.93827160493828_dp)*Hr3(1,0,0) + ((912.5925925925926_dp)*Hr3(1,0,0))/z + &
                        z*(8319.407407407407_dp)*Hr3(1,0,0) + z**2*(-8317.037037037036_dp)*Hr3(1,0,1) + &
                        (-266.5679012345679_dp)*Hr3(1,0,1) + nf*z*(-100.74074074074075_dp)*Hr3(1,0,1) + &
                        nf*(50.370370370370374_dp)*Hr3(1,0,1) + nf*z**2*(100.74074074074075_dp)*Hr3(1,0,1) + &
                        ((839.1111111111111_dp)*Hr3(1,0,1))/z + z*(7107.753086419753_dp)*Hr3(1,0,1) + &
                        z**2*(-8317.037037037036_dp)*Hr3(1,1,0) + (-266.5679012345679_dp)*Hr3(1,1,0) + &
                        nf*z*(-100.74074074074075_dp)*Hr3(1,1,0) + nf*(50.370370370370374_dp)*Hr3(1,1,0) + &
                        nf*z**2*(100.74074074074075_dp)*Hr3(1,1,0) + ((839.1111111111111_dp)*Hr3(1,1,0))/z + &
                        z*(7107.753086419753_dp)*Hr3(1,1,0) + z*(-8606.222222222223_dp)*Hr3(1,1,1) + &
                        ((-782.2222222222222_dp)*Hr3(1,1,1))/z + nf*z**2*(-231.11111111111111_dp)*Hr3(1,1,1) + &
                        nf*(-115.55555555555556_dp)*Hr3(1,1,1) + nf*z*(231.11111111111111_dp)*Hr3(1,1,1) + &
                        (1420.7407407407406_dp)*Hr3(1,1,1) + z**2*(9572.148148148148_dp)*Hr3(1,1,1) + &
                        (2300.0493827160494_dp)*Hr4(-1,-1,-1,0) + z*(4600.098765432099_dp)*Hr4(-1,-1,-1,0) + &
                        z**2*(4600.098765432099_dp)*Hr4(-1,-1,-1,0) + z*(-3922.962962962963_dp)*Hr4(-1,-1,0,0) + &
                        z**2*(-3922.962962962963_dp)*Hr4(-1,-1,0,0) + (-1961.4814814814815_dp)*Hr4(-1,-1,0,0) + &
                        z*(-3245.8271604938273_dp)*Hr4(-1,-1,0,1) + z**2*(-3245.8271604938273_dp)*Hr4(-1,-1,0,1) + &
                        (-1622.9135802469136_dp)*Hr4(-1,-1,0,1) + z*(-2300.0493827160494_dp)*Hr4(-1,0,-1,0) + &
                        z**2*(-2300.0493827160494_dp)*Hr4(-1,0,-1,0) + (-1150.0246913580247_dp)*Hr4(-1,0,-1,0) + &
                        (811.4567901234568_dp)*Hr4(-1,0,0,0) + z*(1622.9135802469136_dp)*Hr4(-1,0,0,0) + &
                        z**2*(1622.9135802469136_dp)*Hr4(-1,0,0,0) + (1284.3456790123457_dp)*Hr4(-1,0,0,1) + &
                        z*(2568.6913580246915_dp)*Hr4(-1,0,0,1) + z**2*(2568.6913580246915_dp)*Hr4(-1,0,0,1) + &
                        (472.8888888888889_dp)*Hr4(-1,0,1,0) + z*(945.7777777777778_dp)*Hr4(-1,0,1,0) + &
                        z**2*(945.7777777777778_dp)*Hr4(-1,0,1,0) + (945.7777777777778_dp)*Hr4(-1,0,1,1) + &
                        z*(1891.5555555555557_dp)*Hr4(-1,0,1,1) + z**2*(1891.5555555555557_dp)*Hr4(-1,0,1,1) + &
                        z*(-8209.382716049382_dp)*Hr4(0,-1,-1,0) + z**2*(-660.5432098765432_dp)*Hr4(0,-1,-1,0) + &
                        (1804.641975308642_dp)*Hr4(0,-1,-1,0) + (-1552.5925925925926_dp)*Hr4(0,-1,0,0) + &
                        z**2*(545.1851851851852_dp)*Hr4(0,-1,0,0) + z*(7028.148148148148_dp)*Hr4(0,-1,0,0) + &
                        (-1300.5432098765432_dp)*Hr4(0,-1,0,1) + z**2*(429.82716049382714_dp)*Hr4(0,-1,0,1) + &
                        z*(5846.913580246914_dp)*Hr4(0,-1,0,1) + nf*z*(-128._dp)*Hr4(0,0,-1,0) + nf*(-64._dp)*Hr4(0,0,-1,0) + &
                        z**2*(88.49382716049382_dp)*Hr4(0,0,-1,0) + (754.1728395061729_dp)*Hr4(0,0,-1,0) + &
                        z*(7417.679012345679_dp)*Hr4(0,0,-1,0) + nf*z*(-440.8888888888889_dp)*Hr4(0,0,0,0) + &
                        z**2*(-151.7037037037037_dp)*Hr4(0,0,0,0) + nf*(-99.55555555555556_dp)*Hr4(0,0,0,0) + &
                        (705.1851851851852_dp)*Hr4(0,0,0,0) + z*(8131.950617283951_dp)*Hr4(0,0,0,0) + &
                        (-1786.469135802469_dp)*Hr4(0,0,0,1) + z**2*(-581.5308641975308_dp)*Hr4(0,0,0,1) + &
                        nf*z*(-265.48148148148147_dp)*Hr4(0,0,0,1) + nf*(132.74074074074073_dp)*Hr4(0,0,0,1) + &
                        z*(12787.358024691359_dp)*Hr4(0,0,0,1) + z**2*(-1001.8765432098766_dp)*Hr4(0,0,1,0) + &
                        nf*z*(-80.5925925925926_dp)*Hr4(0,0,1,0) + nf*(40.2962962962963_dp)*Hr4(0,0,1,0) + &
                        (675.1604938271605_dp)*Hr4(0,0,1,0) + z*(11591.901234567902_dp)*Hr4(0,0,1,0) + &
                        z**2*(-1637.1358024691358_dp)*Hr4(0,0,1,1) + nf*(-11.851851851851851_dp)*Hr4(0,0,1,1) + &
                        nf*z*(23.703703703703702_dp)*Hr4(0,0,1,1) + (2686.814814814815_dp)*Hr4(0,0,1,1) + &
                        z*(9438.814814814816_dp)*Hr4(0,0,1,1) + z**2*(-1091.9506172839506_dp)*Hr4(0,1,0,0) + &
                        (990.8148148148148_dp)*Hr4(0,1,0,0) + z*(6231.7037037037035_dp)*Hr4(0,1,0,0) + &
                        z**2*(-1779.358024691358_dp)*Hr4(0,1,0,1) + (498.5679012345679_dp)*Hr4(0,1,0,1) + &
                        z*(6554.864197530864_dp)*Hr4(0,1,0,1) + z**2*(-1779.358024691358_dp)*Hr4(0,1,1,0) + &
                        (498.5679012345679_dp)*Hr4(0,1,1,0) + z*(6554.864197530864_dp)*Hr4(0,1,1,0) + &
                        z*(-2574.222222222222_dp)*Hr4(0,1,1,1) + (-2232.8888888888887_dp)*Hr4(0,1,1,1) + &
                        z**2*(-1716.148148148148_dp)*Hr4(0,1,1,1) + z**2*(-1366.1234567901236_dp)*Hr4(1,0,0,0) + &
                        (-683.0617283950618_dp)*Hr4(1,0,0,0) + z*(1366.1234567901236_dp)*Hr4(1,0,0,0) + &
                        z**2*(-2883.9506172839506_dp)*Hr4(1,0,0,1) + (-1441.9753086419753_dp)*Hr4(1,0,0,1) + &
                        z*(2883.9506172839506_dp)*Hr4(1,0,0,1) + z**2*(-3035.6543209876545_dp)*Hr4(1,0,1,0) + &
                        (-1517.8271604938273_dp)*Hr4(1,0,1,0) + z*(3035.6543209876545_dp)*Hr4(1,0,1,0) + &
                        z**2*(-1459.358024691358_dp)*Hr4(1,0,1,1) + (-729.679012345679_dp)*Hr4(1,0,1,1) + &
                        z*(1459.358024691358_dp)*Hr4(1,0,1,1) + z**2*(-2883.9506172839506_dp)*Hr4(1,1,0,0) + &
                        (-1441.9753086419753_dp)*Hr4(1,1,0,0) + z*(2883.9506172839506_dp)*Hr4(1,1,0,0) + &
                        z**2*(-1459.358024691358_dp)*Hr4(1,1,0,1) + (-729.679012345679_dp)*Hr4(1,1,0,1) + &
                        z*(1459.358024691358_dp)*Hr4(1,1,0,1) + z**2*(-1459.358024691358_dp)*Hr4(1,1,1,0) + &
                        (-729.679012345679_dp)*Hr4(1,1,1,0) + z*(1459.358024691358_dp)*Hr4(1,1,1,0) + &
                        z*(-6305.185185185185_dp)*Hr4(1,1,1,1) + (3152.5925925925926_dp)*Hr4(1,1,1,1) + &
                        z**2*(6305.185185185185_dp)*Hr4(1,1,1,1)

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
                        Ibar_select = z**2*(-10114.055128069549_dp) + (-2173.0639466267185_dp) + nf*z*(-341.66861008966373_dp) + &
                        (nf*(-78.98765432098766_dp))/z + nf**2*z**2*(-6.378600823045267_dp) + nf**2*(-5.2592592592592595_dp) + &
                        (nf**2*(0.823045267489712_dp))/z + nf**2*z*(11.25925925925926_dp) + nf*(76.55922239749904_dp) + &
                        nf*z**2*(564.7863041563032_dp) + (1546.3709009477598_dp)/z + z*(2107.2553250420533_dp) + &
                        z**2*(-3527.6028945262215_dp)*Hr1(0) + nf*z*(-359.5338835697139_dp)*Hr1(0) + &
                        nf*(-48.307132289217115_dp)*Hr1(0) + (nf*(-12.74074074074074_dp)*Hr1(0))/z + &
                        nf**2*(-1.6296296296296295_dp)*Hr1(0) + nf**2*z**2*(0.3950617283950617_dp)*Hr1(0) + &
                        nf**2*z*(2.3703703703703702_dp)*Hr1(0) + nf*z**2*(58.03292181069959_dp)*Hr1(0) + ((448._dp)*Hr1(0))/z + &
                        (550.7802232282645_dp)*Hr1(0) + z*(4114.483848768637_dp)*Hr1(0) + z**2*(-7414.585565454451_dp)*Hr1(1) + &
                        (-780.3298197642627_dp)*Hr1(1) + nf*z*(-323.4567901234568_dp)*Hr1(1) + &
                        (nf*(-27.588477366255145_dp)*Hr1(1))/z + nf**2*z**2*(-2.074074074074074_dp)*Hr1(1) + &
                        nf**2*(-1.037037037037037_dp)*Hr1(1) + nf**2*z*(2.074074074074074_dp)*Hr1(1) + &
                        nf*(91.75308641975309_dp)*Hr1(1) + nf*z**2*(282.2551440329218_dp)*Hr1(1) + &
                        ((1127.5555555555557_dp)*Hr1(1))/z + z*(7566.733713602599_dp)*Hr1(1) + &
                        z**2*(-359.85185185185185_dp)*Hr2(0,0) + (-185.01234567901236_dp)*Hr2(0,0) + &
                        nf*z*(-159.60493827160494_dp)*Hr2(0,0) + nf**2*(-1.1851851851851851_dp)*Hr2(0,0) + &
                        nf*z**2*(1.5802469135802468_dp)*Hr2(0,0) + nf**2*z*(2.3703703703703702_dp)*Hr2(0,0) + &
                        nf*(20.691358024691358_dp)*Hr2(0,0) + ((48._dp)*Hr2(0,0))/z + z*(3474.6666666666665_dp)*Hr2(0,0) + &
                        z**2*(-1345.037037037037_dp)*Hr2(0,1) + nf*z*(-139.45679012345678_dp)*Hr2(0,1) + &
                        nf*(-0.49382716049382713_dp)*Hr2(0,1) + nf*z**2*(9.481481481481481_dp)*Hr2(0,1) + &
                        ((165.33333333333334_dp)*Hr2(0,1))/z + (389.48148148148147_dp)*Hr2(0,1) + z*(4696.592592592592_dp)*Hr2(0,1) &
                        + z**2*(-2281.185185185185_dp)*Hr2(1,0) + (-253.25925925925927_dp)*Hr2(1,0) + &
                        nf*z*(-46.123456790123456_dp)*Hr2(1,0) + nf*(23.061728395061728_dp)*Hr2(1,0) + &
                        nf*z**2*(46.123456790123456_dp)*Hr2(1,0) + ((196.14814814814815_dp)*Hr2(1,0))/z + &
                        z*(2017.3333333333333_dp)*Hr2(1,0) + (z**2*(-4562.37037037037_dp) + (-506.51851851851853_dp) + &
                        nf*(z*(-92.24691358024691_dp) + (46.123456790123456_dp) + z**2*(92.24691358024691_dp)) + &
                        (392.2962962962963_dp)/z + z*(4034.6666666666665_dp))*Hr2(1,1) + (nf*z*(-47.407407407407405_dp) + &
                        z**2*(-12.641975308641975_dp) + nf*(-8.296296296296296_dp) + (58.76543209876543_dp) + &
                        z*(1295.8024691358025_dp))*Hr3(0,0,0) + (-210.07407407407408_dp)*Hr3(0,0,1) + &
                        z**2*(-66.37037037037037_dp)*Hr3(0,0,1) + nf*z*(-30.814814814814813_dp)*Hr3(0,0,1) + &
                        nf*(15.407407407407407_dp)*Hr3(0,0,1) + z*(2185.4814814814813_dp)*Hr3(0,0,1) + &
                        z**2*(-158.8148148148148_dp)*Hr3(0,1,0) + (254.5185185185185_dp)*Hr3(0,1,0) + &
                        z*(1256.2962962962963_dp)*Hr3(0,1,0) + z**2*(-317.6296296296296_dp)*Hr3(0,1,1) + &
                        (509.037037037037_dp)*Hr3(0,1,1) + z*(2512.5925925925926_dp)*Hr3(0,1,1) + &
                        z**2*(-249.08641975308643_dp)*Hr3(1,0,0) + (-124.54320987654322_dp)*Hr3(1,0,0) + &
                        z*(249.08641975308643_dp)*Hr3(1,0,0) + z**2*(-747.2592592592592_dp)*Hr3(1,0,1) + &
                        (-373.6296296296296_dp)*Hr3(1,0,1) + z*(747.2592592592592_dp)*Hr3(1,0,1) + &
                        z**2*(-747.2592592592592_dp)*Hr3(1,1,0) + (-373.6296296296296_dp)*Hr3(1,1,0) + &
                        z*(747.2592592592592_dp)*Hr3(1,1,0) + z**2*(-1494.5185185185185_dp)*Hr3(1,1,1) + &
                        (-747.2592592592592_dp)*Hr3(1,1,1) + z*(1494.5185185185185_dp)*Hr3(1,1,1)

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
                        Ibar_select = z*(-469412.71240597754_dp) + (-256503.1730595722_dp)/z + nf*z**2*(-179575.3084562035_dp) + &
                        nf*(-20740.03203515044_dp) + (nf**2*(-217.15352709606694_dp))/z + nf**2*(-69.4346605283962_dp) + &
                        nf**3*z*(-61.56266003147559_dp) + (nf**3*(0.943758573388203_dp))/z + nf**3*z**2*(4.7148256035288485_dp) + &
                        nf**3*(50.49292374132484_dp) + nf**2*z**2*(809.2344778215572_dp) + nf**2*z*(2950.6301927274108_dp) + &
                        nf*z*(13890.112928231356_dp) + (nf*(14799.25174218644_dp))/z + (714935.2023348046_dp) + &
                        z**2*(4.776781458354928e6_dp) + ((-231.50923903789854_dp)*(nf*(z**2*(-167.88341762560685_dp) + &
                        z**3*(-156.77448975036322_dp) + z*(-51.26100705841748_dp) + (-12.632826510721248_dp)) + nf**2*z*((0.5_dp) + &
                        z*(1._dp) + z**2*(1._dp)) + (612.121609527476_dp) + z*(653.194857131777_dp) + z**2*(4828.745512695094_dp) + &
                        z**3*(5031.034033311976_dp))*Hr1(-1))/z + z**2*(-511786.7774891476_dp)*Hr1(0) + &
                        (-270559.62785357656_dp)*Hr1(0) + z*(-211141.10465730372_dp)*Hr1(0) + ((-87375.1185450103_dp)*Hr1(0))/z + &
                        nf*z*(-17731.368410972005_dp)*Hr1(0) + nf**2*z**2*(-2557.44208628283_dp)*Hr1(0) + &
                        nf**2*(-766.2776968971109_dp)*Hr1(0) + (nf**2*(-33.39003200731596_dp)*Hr1(0))/z + &
                        nf**3*z*(-4.045119526786653_dp)*Hr1(0) + nf**3*z**2*(9.777777777777779_dp)*Hr1(0) + &
                        nf**3*(31.948485689319252_dp)*Hr1(0) + nf**2*z*(1579.7251588985093_dp)*Hr1(0) + &
                        (nf*(3626.1810660900674_dp)*Hr1(0))/z + nf*(27370.47345111955_dp)*Hr1(0) + &
                        nf*z**2*(116464.10138522225_dp)*Hr1(0) + z*(-1.094946585257935e6_dp)*Hr1(1) + &
                        ((-211175.38024484238_dp)*Hr1(1))/z + nf*z*(-8186.160505709948_dp)*Hr1(1) + &
                        nf**2*(-753.2769164884136_dp)*Hr1(1) + nf**2*z**2*(-634.3206368177051_dp)*Hr1(1) + &
                        (nf**2*(-47.57750342935528_dp)*Hr1(1))/z + nf**3*(-7.209876543209877_dp)*Hr1(1) + &
                        nf**3*z**2*(-0.4609053497942387_dp)*Hr1(1) + (nf**3*(1.0534979423868314_dp)*Hr1(1))/z + &
                        nf**3*z*(10.074074074074074_dp)*Hr1(1) + nf**2*z*(1198.6251635666763_dp)*Hr1(1) + &
                        nf*(1632.9931285868952_dp)*Hr1(1) + (nf*(4470.05782031292_dp)*Hr1(1))/z + &
                        nf*z**2*(5031.790422779148_dp)*Hr1(1) + (221101.81737432192_dp)*Hr1(1) + &
                        z**2*(1.0816866578797677e6_dp)*Hr1(1) + nf*z*(-7998.3598991675935_dp)*Hr2(-1,-1) + &
                        nf*z**2*(-7998.3598991675935_dp)*Hr2(-1,-1) + nf*(-3999.1799495837968_dp)*Hr2(-1,-1) + &
                        ((36014.7144631027_dp)*Hr2(-1,-1))/z + (88148.89934938354_dp)*Hr2(-1,-1) + &
                        z*(450177.2900461982_dp)*Hr2(-1,-1) + z**2*(505664.7746083061_dp)*Hr2(-1,-1) + &
                        nf*z*(-31416.013492202233_dp)*Hr2(-1,0) + nf*z**2*(-30930.888663670004_dp)*Hr2(-1,0) + &
                        (nf*(-4817.152263374485_dp)*Hr2(-1,0))/z + nf*(-2960.1494072808173_dp)*Hr2(-1,0) + &
                        nf**3*z*(-1.7777777777777777_dp)*Hr2(-1,0) + nf**3*z**2*(-1.7777777777777777_dp)*Hr2(-1,0) + &
                        nf**3*(-0.8888888888888888_dp)*Hr2(-1,0) + (nf**2*(58.666666666666664_dp)*Hr2(-1,0))/z + &
                        nf**2*(310.22222222222223_dp)*Hr2(-1,0) + nf**2*z**2*(539.1111111111111_dp)*Hr2(-1,0) + &
                        nf**2*z*(756.4444444444445_dp)*Hr2(-1,0) + (52649.388397839684_dp)*Hr2(-1,0) + &
                        ((72993.0242635347_dp)*Hr2(-1,0))/z + z**2*(354740.23517945776_dp)*Hr2(-1,0) + &
                        z*(391546.36371499597_dp)*Hr2(-1,0) + z*(-633637.4378100501_dp)*Hr2(0,-1) + &
                        z**2*(-96613.15031196545_dp)*Hr2(0,-1) + ((-16083.799764738213_dp)*Hr2(0,-1))/z + &
                        nf*(-909.1408432674247_dp)*Hr2(0,-1) + nf*z**2*(-775.4882580170074_dp)*Hr2(0,-1) + &
                        nf*z*(14222.627945497396_dp)*Hr2(0,-1) + (90318.27368035447_dp)*Hr2(0,-1) + &
                        z**2*(-1.0046028976533923e6_dp)*Hr2(0,0) + nf*z*(-36913.89474757732_dp)*Hr2(0,0) + &
                        ((-16613.793757773707_dp)*Hr2(0,0))/z + (-8683.15801075079_dp)*Hr2(0,0) + &
                        nf*(-3257.7771423546665_dp)*Hr2(0,0) + nf**2*(-36.88580785551619_dp)*Hr2(0,0) + &
                        nf**3*z*(-22.51851851851852_dp)*Hr2(0,0) + nf**2*z**2*(-18.88065843621399_dp)*Hr2(0,0) + &
                        nf**3*z**2*(-3.1604938271604937_dp)*Hr2(0,0) + nf**3*(11.851851851851851_dp)*Hr2(0,0) + &
                        (nf*(370.3703703703704_dp)*Hr2(0,0))/z + nf**2*z*(1602.2052865975018_dp)*Hr2(0,0) + &
                        nf*z**2*(37261.34170147207_dp)*Hr2(0,0) + z*(80391.06212019565_dp)*Hr2(0,0) + &
                        z**2*(-826168.08531365_dp)*Hr2(0,1) + (-111078.45997963127_dp)*Hr2(0,1) + &
                        nf*z*(-40736.04978021934_dp)*Hr2(0,1) + ((-29658.976945997823_dp)*Hr2(0,1))/z + &
                        (nf*(-780.9382716049382_dp)*Hr2(0,1))/z + nf**2*z**2*(-580.0438957475994_dp)*Hr2(0,1) + &
                        nf**2*(-490.40329218106996_dp)*Hr2(0,1) + nf**3*(-3.1604938271604937_dp)*Hr2(0,1) + &
                        nf**3*z*(-0.7901234567901234_dp)*Hr2(0,1) + nf**3*z**2*(3.950617283950617_dp)*Hr2(0,1) + &
                        (nf**2*(13.82716049382716_dp)*Hr2(0,1))/z + nf**2*z*(318.9794238683128_dp)*Hr2(0,1) + &
                        nf*(6142.380477688363_dp)*Hr2(0,1) + z*(25813.744239432213_dp)*Hr2(0,1) + &
                        nf*z**2*(49961.88638986642_dp)*Hr2(0,1) + z**2*(-372682.3332671119_dp)*Hr2(1,0) + &
                        (-68496.9616049669_dp)*Hr2(1,0) + nf*z*(-38898.89746833752_dp)*Hr2(1,0) + &
                        (nf*(-3612.1371742112483_dp)*Hr2(1,0))/z + nf**2*z**2*(-593.5912208504801_dp)*Hr2(1,0) + &
                        nf**2*(-245.33333333333334_dp)*Hr2(1,0) + nf**3*z*(-3.1604938271604937_dp)*Hr2(1,0) + &
                        nf**3*(1.5802469135802468_dp)*Hr2(1,0) + nf**3*z**2*(3.1604938271604937_dp)*Hr2(1,0) + &
                        (nf**2*(35.072702331961594_dp)*Hr2(1,0))/z + nf**2*z*(666.1399176954733_dp)*Hr2(1,0) + &
                        nf*(8400.104426898526_dp)*Hr2(1,0) + nf*z**2*(37027.418730340265_dp)*Hr2(1,0) + &
                        ((41377.644673635834_dp)*Hr2(1,0))/z + z*(384960.44086821366_dp)*Hr2(1,0) + &
                        ((-0.9876543209876543_dp)*(z**3*(-216283.6685705656_dp) + nf**3*z*(z*(-1._dp) + (0.5_dp) + z**2*(1._dp)) + &
                        nf**2*(z*(-103.86666666666667_dp) + z**3*(15.427777777777777_dp) + (34.82222222222223_dp) + &
                        z**2*(190.03333333333333_dp)) + z*(322.8573478640527_dp) + nf*(z**3*(-12159.942594447513_dp) + &
                        z*(-8552.907408334868_dp) + (-1418.5555555555557_dp) + z**2*(14465.781483336403_dp)) + &
                        (77480.28562270424_dp) + z**2*(226520.44619339827_dp))*Hr2(1,1))/z + z*(-123120.6748859516_dp)*Hr3(-1,-1,-1) &
                        + z**2*(-123120.6748859516_dp)*Hr3(-1,-1,-1) + (-61560.3374429758_dp)*Hr3(-1,-1,-1) + &
                        z**2*(-332323.7295988669_dp)*Hr3(-1,-1,0) + z*(-330134.19873466936_dp)*Hr3(-1,-1,0) + &
                        ((-54522.07407407407_dp)*Hr3(-1,-1,0))/z + (-19979.000601902597_dp)*Hr3(-1,-1,0) + &
                        nf**2*z*(-112._dp)*Hr3(-1,-1,0) + nf**2*z**2*(-112._dp)*Hr3(-1,-1,0) + nf**2*(-56._dp)*Hr3(-1,-1,0) + &
                        (nf*(1405.9807956104253_dp)*Hr3(-1,-1,0))/z + nf*(4113.283950617284_dp)*Hr3(-1,-1,0) + &
                        nf*z**2*(14543.758573388202_dp)*Hr3(-1,-1,0) + nf*z*(15652.543209876543_dp)*Hr3(-1,-1,0) + &
                        (30780.1687214879_dp)*Hr3(-1,0,-1) + z*(61560.3374429758_dp)*Hr3(-1,0,-1) + &
                        z**2*(61560.3374429758_dp)*Hr3(-1,0,-1) + nf*z*(-13870.222222222223_dp)*Hr3(-1,0,0) + &
                        nf*z**2*(-12811.286694101509_dp)*Hr3(-1,0,0) + nf*(-3707.1111111111113_dp)*Hr3(-1,0,0) + &
                        (nf*(-1240.471879286694_dp)*Hr3(-1,0,0))/z + nf**2*(49.18518518518518_dp)*Hr3(-1,0,0) + &
                        nf**2*z*(98.37037037037037_dp)*Hr3(-1,0,0) + nf**2*z**2*(98.37037037037037_dp)*Hr3(-1,0,0) + &
                        (39612.60754553488_dp)*Hr3(-1,0,0) + ((48200.88888888889_dp)*Hr3(-1,0,0))/z + &
                        z*(334984.37558489694_dp)*Hr3(-1,0,0) + z**2*(336619.3138565019_dp)*Hr3(-1,0,0) + &
                        nf*z*(-12087.901234567902_dp)*Hr3(-1,0,1) + nf*z**2*(-11078.814814814816_dp)*Hr3(-1,0,1) + &
                        nf*(-3300.938271604938_dp)*Hr3(-1,0,1) + (nf*(-1074.962962962963_dp)*Hr3(-1,0,1))/z + &
                        nf**2*(42.370370370370374_dp)*Hr3(-1,0,1) + nf**2*z*(84.74074074074075_dp)*Hr3(-1,0,1) + &
                        nf**2*z**2*(84.74074074074075_dp)*Hr3(-1,0,1) + ((41879.7037037037_dp)*Hr3(-1,0,1))/z + &
                        (42191.538084084736_dp)*Hr3(-1,0,1) + z*(305725.19962495955_dp)*Hr3(-1,0,1) + &
                        z**2*(306805.5453039719_dp)*Hr3(-1,0,1) + (-50252.93882049317_dp)*Hr3(0,-1,-1) + &
                        z**2*(15076.531496643496_dp)*Hr3(0,-1,-1) + z*(223626.5525269379_dp)*Hr3(0,-1,-1) + &
                        nf*z*(-23036.246913580246_dp)*Hr3(0,-1,0) + nf*z**2*(-799.8573388203018_dp)*Hr3(0,-1,0) + &
                        (nf*(-707.5555555555555_dp)*Hr3(0,-1,0))/z + nf**2*z**2*(-64._dp)*Hr3(0,-1,0) + nf**2*(56._dp)*Hr3(0,-1,0) + &
                        nf**2*z*(192._dp)*Hr3(0,-1,0) + nf*(2162.2716049382716_dp)*Hr3(0,-1,0) + (15714.566078368625_dp)*Hr3(0,-1,0) &
                        + ((23066.666666666668_dp)*Hr3(0,-1,0))/z + z**2*(136411.32347492952_dp)*Hr3(0,-1,0) + &
                        z*(245250.54805941368_dp)*Hr3(0,-1,0) + z*(-204988.84031470795_dp)*Hr3(0,0,-1) + &
                        (-21461.312615372906_dp)*Hr3(0,0,-1) + z**2*(-1767.593347882341_dp)*Hr3(0,0,-1) + &
                        nf*(1468.662119931651_dp)*Hr3(0,0,-1) + nf*z*(2937.324239863302_dp)*Hr3(0,0,-1) + &
                        z**2*(-207924.6090739719_dp)*Hr3(0,0,0) + (-57346.74021085067_dp)*Hr3(0,0,0) + &
                        nf*z*(-34400.74412519822_dp)*Hr3(0,0,0) + ((-1248._dp)*Hr3(0,0,0))/z + &
                        nf**2*(-227.81893004115227_dp)*Hr3(0,0,0) + nf**3*z*(-9.481481481481481_dp)*Hr3(0,0,0) + &
                        nf**3*(4.7407407407407405_dp)*Hr3(0,0,0) + nf**2*z**2*(21.069958847736626_dp)*Hr3(0,0,0) + &
                        nf**2*z*(1237.8600823045267_dp)*Hr3(0,0,0) + nf*z**2*(2245.9588477366256_dp)*Hr3(0,0,0) + &
                        nf*(6216.673849118416_dp)*Hr3(0,0,0) + z*(233334.9971298778_dp)*Hr3(0,0,0) + &
                        z**2*(-413137.79208649223_dp)*Hr3(0,0,1) + nf*z*(-35490.05725558168_dp)*Hr3(0,0,1) + &
                        nf*(-2439.0701376412567_dp)*Hr3(0,0,1) + (-1967.3629137217285_dp)*Hr3(0,0,1) + &
                        (nf*(-401.77777777777777_dp)*Hr3(0,0,1))/z + nf**2*(-166.97942386831275_dp)*Hr3(0,0,1) + &
                        nf**3*(-2.3703703703703702_dp)*Hr3(0,0,1) + nf**3*z*(4.7407407407407405_dp)*Hr3(0,0,1) + &
                        nf**2*z**2*(34.23868312757202_dp)*Hr3(0,0,1) + nf**2*z*(479.34156378600824_dp)*Hr3(0,0,1) + &
                        nf*z**2*(5673.6570644718795_dp)*Hr3(0,0,1) + ((8037.333333333333_dp)*Hr3(0,0,1))/z + &
                        z*(314841.42499234935_dp)*Hr3(0,0,1) + z**2*(-367647.50183608226_dp)*Hr3(0,1,0) + &
                        nf*z*(-27435.456790123455_dp)*Hr3(0,1,0) + (-26346.795132831525_dp)*Hr3(0,1,0) + &
                        (nf*(-993.1851851851852_dp)*Hr3(0,1,0))/z + nf**2*z**2*(-53.465020576131685_dp)*Hr3(0,1,0) + &
                        nf**2*(-27.45679012345679_dp)*Hr3(0,1,0) + nf**2*z*(210.17283950617283_dp)*Hr3(0,1,0) + &
                        nf*(482.0740740740741_dp)*Hr3(0,1,0) + nf*z**2*(9330.611796982168_dp)*Hr3(0,1,0) + &
                        ((29705.48148148148_dp)*Hr3(0,1,0))/z + z*(424503.8161966983_dp)*Hr3(0,1,0) + &
                        z**2*(-428451.08802174754_dp)*Hr3(0,1,1) + (-48349.72557334837_dp)*Hr3(0,1,1) + &
                        nf*z*(-12348.04938271605_dp)*Hr3(0,1,1) + (nf*(-744.2962962962963_dp)*Hr3(0,1,1))/z + &
                        nf**2*z**2*(-131.4238683127572_dp)*Hr3(0,1,1) + nf**2*z*(2.897119341563786_dp)*Hr3(0,1,1) + &
                        nf**2*(29.069958847736626_dp)*Hr3(0,1,1) + nf*(4390.222222222223_dp)*Hr3(0,1,1) + &
                        nf*z**2*(13331.753086419752_dp)*Hr3(0,1,1) + ((16632.296296296296_dp)*Hr3(0,1,1))/z + &
                        z*(223835.94108479907_dp)*Hr3(0,1,1) + z**2*(-307239.2288355537_dp)*Hr3(1,0,0) + &
                        (-47679.93129020484_dp)*Hr3(1,0,0) + nf*z*(-12958.485596707818_dp)*Hr3(1,0,0) + &
                        (nf*(-998.8696844993141_dp)*Hr3(1,0,0))/z + nf**2*z**2*(-107.39094650205762_dp)*Hr3(1,0,0) + &
                        nf**2*(-53.69547325102881_dp)*Hr3(1,0,0) + nf**2*z*(107.39094650205762_dp)*Hr3(1,0,0) + &
                        nf*(3674.1563786008232_dp)*Hr3(1,0,0) + nf*z**2*(12165.997256515775_dp)*Hr3(1,0,0) + &
                        ((38832.0658436214_dp)*Hr3(1,0,0))/z + z*(308175.549823208_dp)*Hr3(1,0,0) + &
                        z**2*(-279225.79831414483_dp)*Hr3(1,0,1) + (-50038.38063855389_dp)*Hr3(1,0,1) + &
                        nf*z*(-12479.802469135802_dp)*Hr3(1,0,1) + (nf*(-857.9862825788751_dp)*Hr3(1,0,1))/z + &
                        nf**2*z**2*(-73.94238683127573_dp)*Hr3(1,0,1) + nf**2*(-36.97119341563786_dp)*Hr3(1,0,1) + &
                        nf**2*z*(73.94238683127573_dp)*Hr3(1,0,1) + nf*(3833.58024691358_dp)*Hr3(1,0,1) + &
                        nf*z**2*(11732.06035665295_dp)*Hr3(1,0,1) + ((31297.975308641977_dp)*Hr3(1,0,1))/z + &
                        z*(279602.88473389787_dp)*Hr3(1,0,1) + z**2*(-276835.2149578667_dp)*Hr3(1,1,0) + &
                        (-48843.08896041483_dp)*Hr3(1,1,0) + nf*z*(-12479.802469135802_dp)*Hr3(1,1,0) + &
                        (nf*(-857.9862825788751_dp)*Hr3(1,1,0))/z + nf**2*z**2*(-73.94238683127573_dp)*Hr3(1,1,0) + &
                        nf**2*(-36.97119341563786_dp)*Hr3(1,1,0) + nf**2*z*(73.94238683127573_dp)*Hr3(1,1,0) + &
                        nf*(3833.58024691358_dp)*Hr3(1,1,0) + nf*z**2*(11732.06035665295_dp)*Hr3(1,1,0) + &
                        ((31297.975308641977_dp)*Hr3(1,1,0))/z + z*(277212.30137761973_dp)*Hr3(1,1,0) + &
                        z*(-168578.75429844903_dp)*Hr3(1,1,1) + ((-41565.82716049383_dp)*Hr3(1,1,1))/z + &
                        (-17270.52408534338_dp)*Hr3(1,1,1) + nf*z**2*(-6499.22633744856_dp)*Hr3(1,1,1) + &
                        nf*(-740.2469135802469_dp)*Hr3(1,1,1) + nf**2*z*(-141.23456790123456_dp)*Hr3(1,1,1) + &
                        nf**2*(70.61728395061728_dp)*Hr3(1,1,1) + nf**2*z**2*(141.23456790123456_dp)*Hr3(1,1,1) + &
                        (nf*(919.1769547325102_dp)*Hr3(1,1,1))/z + nf*z*(6783.0123456790125_dp)*Hr3(1,1,1) + &
                        z**2*(168564.8283725231_dp)*Hr3(1,1,1) + nf*z*(-3916.641975308642_dp)*Hr4(-1,-1,-1,0) + &
                        nf*z**2*(-3916.641975308642_dp)*Hr4(-1,-1,-1,0) + nf*(-1958.320987654321_dp)*Hr4(-1,-1,-1,0) + &
                        ((17186.765432098764_dp)*Hr4(-1,-1,-1,0))/z + (20197.333333333332_dp)*Hr4(-1,-1,-1,0) + &
                        z*(171163.25925925927_dp)*Hr4(-1,-1,-1,0) + z**2*(197504.3950617284_dp)*Hr4(-1,-1,-1,0) + &
                        z**2*(-174093.33333333334_dp)*Hr4(-1,-1,0,0) + z*(-150641.8765432099_dp)*Hr4(-1,-1,0,0) + &
                        (-17350.172839506173_dp)*Hr4(-1,-1,0,0) + ((-15243.851851851852_dp)*Hr4(-1,-1,0,0))/z + &
                        nf*(1705.1851851851852_dp)*Hr4(-1,-1,0,0) + nf*z*(3410.3703703703704_dp)*Hr4(-1,-1,0,0) + &
                        nf*z**2*(3410.3703703703704_dp)*Hr4(-1,-1,0,0) + z**2*(-150682.27160493826_dp)*Hr4(-1,-1,0,1) + &
                        z*(-130120.49382716049_dp)*Hr4(-1,-1,0,1) + (-14503.012345679012_dp)*Hr4(-1,-1,0,1) + &
                        ((-13300.938271604939_dp)*Hr4(-1,-1,0,1))/z + nf*(1452.0493827160494_dp)*Hr4(-1,-1,0,1) + &
                        nf*z*(2904.098765432099_dp)*Hr4(-1,-1,0,1) + nf*z**2*(2904.098765432099_dp)*Hr4(-1,-1,0,1) + &
                        z**2*(-98752.1975308642_dp)*Hr4(-1,0,-1,0) + z*(-85581.62962962964_dp)*Hr4(-1,0,-1,0) + &
                        (-10098.666666666666_dp)*Hr4(-1,0,-1,0) + ((-8593.382716049382_dp)*Hr4(-1,0,-1,0))/z + &
                        nf*(979.1604938271605_dp)*Hr4(-1,0,-1,0) + nf*z*(1958.320987654321_dp)*Hr4(-1,0,-1,0) + &
                        nf*z**2*(1958.320987654321_dp)*Hr4(-1,0,-1,0) + nf*z*(-1452.0493827160494_dp)*Hr4(-1,0,0,0) + &
                        nf*z**2*(-1452.0493827160494_dp)*Hr4(-1,0,0,0) + nf*(-726.0246913580247_dp)*Hr4(-1,0,0,0) + &
                        ((6650.469135802469_dp)*Hr4(-1,0,0,0))/z + (7251.506172839506_dp)*Hr4(-1,0,0,0) + &
                        z*(65060.246913580246_dp)*Hr4(-1,0,0,0) + z**2*(75341.13580246913_dp)*Hr4(-1,0,0,0) + &
                        nf*z*(-2397.8271604938273_dp)*Hr4(-1,0,0,1) + nf*z**2*(-2397.8271604938273_dp)*Hr4(-1,0,0,1) + &
                        nf*(-1198.9135802469136_dp)*Hr4(-1,0,0,1) + ((11358.024691358025_dp)*Hr4(-1,0,0,1))/z + &
                        (11655.851851851852_dp)*Hr4(-1,0,0,1) + z*(109599.11111111111_dp)*Hr4(-1,0,0,1) + &
                        z**2*(127271.20987654322_dp)*Hr4(-1,0,0,1) + nf*z*(-945.7777777777778_dp)*Hr4(-1,0,1,0) + &
                        nf*z**2*(-945.7777777777778_dp)*Hr4(-1,0,1,0) + nf*(-472.8888888888889_dp)*Hr4(-1,0,1,0) + &
                        (4404.3456790123455_dp)*Hr4(-1,0,1,0) + ((4707.555555555556_dp)*Hr4(-1,0,1,0))/z + &
                        z*(44538.864197530864_dp)*Hr4(-1,0,1,0) + z**2*(51930.07407407407_dp)*Hr4(-1,0,1,0) + &
                        nf*z*(-1891.5555555555557_dp)*Hr4(-1,0,1,1) + nf*z**2*(-1891.5555555555557_dp)*Hr4(-1,0,1,1) + &
                        nf*(-945.7777777777778_dp)*Hr4(-1,0,1,1) + (8808.691358024691_dp)*Hr4(-1,0,1,1) + &
                        ((9415.111111111111_dp)*Hr4(-1,0,1,1))/z + z*(89077.72839506173_dp)*Hr4(-1,0,1,1) + &
                        z**2*(103860.14814814815_dp)*Hr4(-1,0,1,1) + z*(-221826.37037037036_dp)*Hr4(0,-1,-1,0) + &
                        z**2*(-41589.72839506173_dp)*Hr4(0,-1,-1,0) + ((-7651.555555555556_dp)*Hr4(0,-1,-1,0))/z + &
                        nf*(-493.4320987654321_dp)*Hr4(0,-1,-1,0) + nf*z**2*(-301.30041152263374_dp)*Hr4(0,-1,-1,0) + &
                        nf*z*(6932.543209876543_dp)*Hr4(0,-1,-1,0) + (24821.925925925927_dp)*Hr4(0,-1,-1,0) + &
                        (-21769.62962962963_dp)*Hr4(0,-1,0,0) + nf*z*(-6056.2962962962965_dp)*Hr4(0,-1,0,0) + &
                        nf*z**2*(311.04526748971193_dp)*Hr4(0,-1,0,0) + nf*(399.7037037037037_dp)*Hr4(0,-1,0,0) + &
                        ((6801.777777777777_dp)*Hr4(0,-1,0,0))/z + z**2*(36292.444444444445_dp)*Hr4(0,-1,0,0) + &
                        z*(195294.22222222222_dp)*Hr4(0,-1,0,0) + (-18717.333333333332_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(-5180.049382716049_dp)*Hr4(0,-1,0,1) + nf*(305.9753086419753_dp)*Hr4(0,-1,0,1) + &
                        nf*z**2*(320.7901234567901_dp)*Hr4(0,-1,0,1) + ((5952._dp)*Hr4(0,-1,0,1))/z + &
                        z**2*(30995.16049382716_dp)*Hr4(0,-1,0,1) + z*(168762.07407407407_dp)*Hr4(0,-1,0,1) + &
                        nf*z*(-8749.037037037036_dp)*Hr4(0,0,-1,0) + nf*(-1242.6666666666667_dp)*Hr4(0,0,-1,0) + &
                        nf**2*(64._dp)*Hr4(0,0,-1,0) + nf**2*z*(128._dp)*Hr4(0,0,-1,0) + &
                        nf*z**2*(136.95473251028807_dp)*Hr4(0,0,-1,0) + ((2304._dp)*Hr4(0,0,-1,0))/z + &
                        z**2*(9733.530864197532_dp)*Hr4(0,0,-1,0) + (11713.283950617284_dp)*Hr4(0,0,-1,0) + &
                        z*(178422.5185185185_dp)*Hr4(0,0,-1,0) + nf*z*(-19525.53086419753_dp)*Hr4(0,0,0,0) + &
                        z**2*(-13610.930041152264_dp)*Hr4(0,0,0,0) + nf*(-1606.320987654321_dp)*Hr4(0,0,0,0) + &
                        nf**2*(44.24691358024691_dp)*Hr4(0,0,0,0) + nf**2*z*(442.4691358024691_dp)*Hr4(0,0,0,0) + &
                        (10795.98353909465_dp)*Hr4(0,0,0,0) + z*(199474.96296296295_dp)*Hr4(0,0,0,0) + &
                        z**2*(-44490.33744855967_dp)*Hr4(0,0,0,1) + (-26049.975308641977_dp)*Hr4(0,0,0,1) + &
                        nf*z*(-20137.876543209877_dp)*Hr4(0,0,0,1) + nf**2*(-139.85185185185185_dp)*Hr4(0,0,0,1) + &
                        nf*z**2*(23.17695473251029_dp)*Hr4(0,0,0,1) + nf**2*z*(194.37037037037038_dp)*Hr4(0,0,0,1) + &
                        ((1728._dp)*Hr4(0,0,0,1))/z + nf*(2326.1234567901233_dp)*Hr4(0,0,0,1) + &
                        z*(322460.83950617287_dp)*Hr4(0,0,0,1) + z**2*(-65503.60493827161_dp)*Hr4(0,0,1,0) + &
                        nf*z*(-12795.654320987655_dp)*Hr4(0,0,1,0) + nf**2*(-22.91358024691358_dp)*Hr4(0,0,1,0) + &
                        nf**2*z*(45.82716049382716_dp)*Hr4(0,0,1,0) + nf*z**2*(245.9917695473251_dp)*Hr4(0,0,1,0) + &
                        nf*(884.3456790123457_dp)*Hr4(0,0,1,0) + ((5952._dp)*Hr4(0,0,1,0))/z + (6837.761316872428_dp)*Hr4(0,0,1,0) + &
                        z*(309175.1111111111_dp)*Hr4(0,0,1,0) + z**2*(-105388.44444444444_dp)*Hr4(0,0,1,1) + &
                        nf*z*(-9833.876543209877_dp)*Hr4(0,0,1,1) + nf*(-295.5061728395062_dp)*Hr4(0,0,1,1) + &
                        nf**2*z*(-105.87654320987654_dp)*Hr4(0,0,1,1) + nf**2*(52.93827160493827_dp)*Hr4(0,0,1,1) + &
                        nf*z**2*(420.34567901234567_dp)*Hr4(0,0,1,1) + ((7616._dp)*Hr4(0,0,1,1))/z + &
                        (34470.84773662552_dp)*Hr4(0,0,1,1) + z*(279785.0864197531_dp)*Hr4(0,0,1,1) + &
                        z**2*(-66497.74485596707_dp)*Hr4(0,1,0,0) + nf*z*(-5641.481481481482_dp)*Hr4(0,1,0,0) + &
                        nf*(-128.2962962962963_dp)*Hr4(0,1,0,0) + nf*z**2*(642.8971193415638_dp)*Hr4(0,1,0,0) + &
                        ((8195.555555555555_dp)*Hr4(0,1,0,0))/z + (14775.325102880659_dp)*Hr4(0,1,0,0) + &
                        z*(191551.67078189302_dp)*Hr4(0,1,0,0) + z**2*(-114243.55555555556_dp)*Hr4(0,1,0,1) + &
                        nf*z*(-5659.6543209876545_dp)*Hr4(0,1,0,1) + nf*(368.1975308641975_dp)*Hr4(0,1,0,1) + &
                        nf*z**2*(1361.6460905349795_dp)*Hr4(0,1,0,1) + (9000.493827160493_dp)*Hr4(0,1,0,1) + &
                        ((12650.666666666666_dp)*Hr4(0,1,0,1))/z + z*(218116.34567901236_dp)*Hr4(0,1,0,1) + &
                        z**2*(-114243.55555555556_dp)*Hr4(0,1,1,0) + nf*z*(-5659.6543209876545_dp)*Hr4(0,1,1,0) + &
                        nf*(368.1975308641975_dp)*Hr4(0,1,1,0) + nf*z**2*(1361.6460905349795_dp)*Hr4(0,1,1,0) + &
                        (9000.493827160493_dp)*Hr4(0,1,1,0) + ((12650.666666666666_dp)*Hr4(0,1,1,0))/z + &
                        z*(218116.34567901236_dp)*Hr4(0,1,1,0) + z**2*(-127345.58024691358_dp)*Hr4(0,1,1,1) + &
                        (-26481.382716049382_dp)*Hr4(0,1,1,1) + z*(97.58024691358025_dp)*Hr4(0,1,1,1) + &
                        nf*(1770.6666666666667_dp)*Hr4(0,1,1,1) + nf*z**2*(2127.0123456790125_dp)*Hr4(0,1,1,1) + &
                        nf*z*(2574.222222222222_dp)*Hr4(0,1,1,1) + ((10816._dp)*Hr4(0,1,1,1))/z + &
                        z**2*(-67039.76954732511_dp)*Hr4(1,0,0,0) + (-7325.382716049383_dp)*Hr4(1,0,0,0) + &
                        nf*z*(-1356.4444444444443_dp)*Hr4(1,0,0,0) + nf*(678.2222222222222_dp)*Hr4(1,0,0,0) + &
                        nf*z**2*(1356.4444444444443_dp)*Hr4(1,0,0,0) + ((5753.679012345679_dp)*Hr4(1,0,0,0))/z + &
                        z*(58770.732510288064_dp)*Hr4(1,0,0,0) + z**2*(-138821.49794238684_dp)*Hr4(1,0,0,1) + &
                        (-14011.506172839507_dp)*Hr4(1,0,0,1) + nf*z*(-2693.7283950617284_dp)*Hr4(1,0,0,1) + &
                        nf*(1346.8641975308642_dp)*Hr4(1,0,0,1) + nf*z**2*(2693.7283950617284_dp)*Hr4(1,0,0,1) + &
                        ((12142.617283950618_dp)*Hr4(1,0,0,1))/z + z*(121152.75720164609_dp)*Hr4(1,0,0,1) + &
                        z**2*(-143563.45679012345_dp)*Hr4(1,0,1,0) + (-13372.246913580248_dp)*Hr4(1,0,1,0) + &
                        nf*z*(-2674.567901234568_dp)*Hr4(1,0,1,0) + nf*(1337.283950617284_dp)*Hr4(1,0,1,0) + &
                        nf*z**2*(2674.567901234568_dp)*Hr4(1,0,1,0) + ((12777.876543209877_dp)*Hr4(1,0,1,0))/z + &
                        z*(124764.04938271605_dp)*Hr4(1,0,1,0) + z**2*(-87556.14814814815_dp)*Hr4(1,0,1,1) + &
                        (-2912.9876543209875_dp)*Hr4(1,0,1,1) + nf*z*(-1098.2716049382716_dp)*Hr4(1,0,1,1) + &
                        nf*(549.1358024691358_dp)*Hr4(1,0,1,1) + nf*z**2*(1098.2716049382716_dp)*Hr4(1,0,1,1) + &
                        ((8854.913580246914_dp)*Hr4(1,0,1,1))/z + z*(74033.77777777778_dp)*Hr4(1,0,1,1) + &
                        z**2*(-138821.49794238684_dp)*Hr4(1,1,0,0) + (-14011.506172839507_dp)*Hr4(1,1,0,0) + &
                        nf*z*(-2693.7283950617284_dp)*Hr4(1,1,0,0) + nf*(1346.8641975308642_dp)*Hr4(1,1,0,0) + &
                        nf*z**2*(2693.7283950617284_dp)*Hr4(1,1,0,0) + ((12142.617283950618_dp)*Hr4(1,1,0,0))/z + &
                        z*(121152.75720164609_dp)*Hr4(1,1,0,0) + z**2*(-87556.14814814815_dp)*Hr4(1,1,0,1) + &
                        (-2912.9876543209875_dp)*Hr4(1,1,0,1) + nf*z*(-1098.2716049382716_dp)*Hr4(1,1,0,1) + &
                        nf*(549.1358024691358_dp)*Hr4(1,1,0,1) + nf*z**2*(1098.2716049382716_dp)*Hr4(1,1,0,1) + &
                        ((8854.913580246914_dp)*Hr4(1,1,0,1))/z + z*(74033.77777777778_dp)*Hr4(1,1,0,1) + &
                        z**2*(-87556.14814814815_dp)*Hr4(1,1,1,0) + (-2912.9876543209875_dp)*Hr4(1,1,1,0) + &
                        nf*z*(-1098.2716049382716_dp)*Hr4(1,1,1,0) + nf*(549.1358024691358_dp)*Hr4(1,1,1,0) + &
                        nf*z**2*(1098.2716049382716_dp)*Hr4(1,1,1,0) + ((8854.913580246914_dp)*Hr4(1,1,1,0))/z + &
                        z*(74033.77777777778_dp)*Hr4(1,1,1,0) + z*(-202921.0864197531_dp)*Hr4(1,1,1,1) + &
                        ((-15691.851851851852_dp)*Hr4(1,1,1,1))/z + nf*z**2*(-6305.185185185185_dp)*Hr4(1,1,1,1) + &
                        nf*(-3152.5925925925926_dp)*Hr4(1,1,1,1) + nf*z*(6305.185185185185_dp)*Hr4(1,1,1,1) + &
                        (41837.03703703704_dp)*Hr4(1,1,1,1) + z**2*(224029.23456790124_dp)*Hr4(1,1,1,1) + &
                        z*(-59031.7037037037_dp)*Hr5(-1,-1,-1,-1,0) + z**2*(-59031.7037037037_dp)*Hr5(-1,-1,-1,-1,0) + &
                        (-29515.85185185185_dp)*Hr5(-1,-1,-1,-1,0) + (26091.061728395063_dp)*Hr5(-1,-1,-1,0,0) + &
                        z*(52182.12345679013_dp)*Hr5(-1,-1,-1,0,0) + z**2*(52182.12345679013_dp)*Hr5(-1,-1,-1,0,0) + &
                        (22666.271604938273_dp)*Hr5(-1,-1,-1,0,1) + z*(45332.543209876545_dp)*Hr5(-1,-1,-1,0,1) + &
                        z**2*(45332.543209876545_dp)*Hr5(-1,-1,-1,0,1) + (14757.925925925925_dp)*Hr5(-1,-1,0,-1,0) + &
                        z*(29515.85185185185_dp)*Hr5(-1,-1,0,-1,0) + z**2*(29515.85185185185_dp)*Hr5(-1,-1,0,-1,0) + &
                        z*(-22666.271604938273_dp)*Hr5(-1,-1,0,0,0) + z**2*(-22666.271604938273_dp)*Hr5(-1,-1,0,0,0) + &
                        (-11333.135802469136_dp)*Hr5(-1,-1,0,0,0) + z*(-38482.96296296296_dp)*Hr5(-1,-1,0,0,1) + &
                        z**2*(-38482.96296296296_dp)*Hr5(-1,-1,0,0,1) + (-19241.48148148148_dp)*Hr5(-1,-1,0,0,1) + &
                        z*(-15816.691358024691_dp)*Hr5(-1,-1,0,1,0) + z**2*(-15816.691358024691_dp)*Hr5(-1,-1,0,1,0) + &
                        (-7908.3456790123455_dp)*Hr5(-1,-1,0,1,0) + z*(-31633.382716049382_dp)*Hr5(-1,-1,0,1,1) + &
                        z**2*(-31633.382716049382_dp)*Hr5(-1,-1,0,1,1) + (-15816.691358024691_dp)*Hr5(-1,-1,0,1,1) + &
                        (14757.925925925925_dp)*Hr5(-1,0,-1,-1,0) + z*(29515.85185185185_dp)*Hr5(-1,0,-1,-1,0) + &
                        z**2*(29515.85185185185_dp)*Hr5(-1,0,-1,-1,0) + z*(-26091.061728395063_dp)*Hr5(-1,0,-1,0,0) + &
                        z**2*(-26091.061728395063_dp)*Hr5(-1,0,-1,0,0) + (-13045.530864197532_dp)*Hr5(-1,0,-1,0,0) + &
                        z*(-22666.271604938273_dp)*Hr5(-1,0,-1,0,1) + z**2*(-22666.271604938273_dp)*Hr5(-1,0,-1,0,1) + &
                        (-11333.135802469136_dp)*Hr5(-1,0,-1,0,1) + z*(-9838.617283950618_dp)*Hr5(-1,0,0,-1,0) + &
                        z**2*(-9838.617283950618_dp)*Hr5(-1,0,0,-1,0) + (-4919.308641975309_dp)*Hr5(-1,0,0,-1,0) + &
                        (3206.9135802469136_dp)*Hr5(-1,0,0,0,0) + z*(6413.827160493827_dp)*Hr5(-1,0,0,0,0) + &
                        z**2*(6413.827160493827_dp)*Hr5(-1,0,0,0,0) + (7908.3456790123455_dp)*Hr5(-1,0,0,0,1) + &
                        z*(15816.691358024691_dp)*Hr5(-1,0,0,0,1) + z**2*(15816.691358024691_dp)*Hr5(-1,0,0,0,1) + &
                        (6195.950617283951_dp)*Hr5(-1,0,0,1,0) + z*(12391.901234567902_dp)*Hr5(-1,0,0,1,0) + &
                        z**2*(12391.901234567902_dp)*Hr5(-1,0,0,1,0) + (12391.901234567902_dp)*Hr5(-1,0,0,1,1) + &
                        z*(24783.802469135804_dp)*Hr5(-1,0,0,1,1) + z**2*(24783.802469135804_dp)*Hr5(-1,0,0,1,1) + &
                        (1494.5185185185185_dp)*Hr5(-1,0,1,0,0) + z*(2989.037037037037_dp)*Hr5(-1,0,1,0,0) + &
                        z**2*(2989.037037037037_dp)*Hr5(-1,0,1,0,0) + (4483.555555555556_dp)*Hr5(-1,0,1,0,1) + &
                        z*(8967.111111111111_dp)*Hr5(-1,0,1,0,1) + z**2*(8967.111111111111_dp)*Hr5(-1,0,1,0,1) + &
                        (4483.555555555556_dp)*Hr5(-1,0,1,1,0) + z*(8967.111111111111_dp)*Hr5(-1,0,1,1,0) + &
                        z**2*(8967.111111111111_dp)*Hr5(-1,0,1,1,0) + (8967.111111111111_dp)*Hr5(-1,0,1,1,1) + &
                        z*(17934.222222222223_dp)*Hr5(-1,0,1,1,1) + z**2*(17934.222222222223_dp)*Hr5(-1,0,1,1,1) + &
                        (-23912.296296296296_dp)*Hr5(0,-1,-1,-1,0) + z**2*(7471.407407407408_dp)*Hr5(0,-1,-1,-1,0) + &
                        z*(106856.29629629629_dp)*Hr5(0,-1,-1,-1,0) + z*(-94688.3950617284_dp)*Hr5(0,-1,-1,0,0) + &
                        z**2*(-6450.567901234568_dp)*Hr5(0,-1,-1,0,0) + (21253.135802469136_dp)*Hr5(0,-1,-1,0,0) + &
                        z*(-82520.49382716049_dp)*Hr5(0,-1,-1,0,1) + z**2*(-5429.728395061728_dp)*Hr5(0,-1,-1,0,1) + &
                        (18593.975308641977_dp)*Hr5(0,-1,-1,0,1) + z*(-53428.148148148146_dp)*Hr5(0,-1,0,-1,0) + &
                        z**2*(-3735.703703703704_dp)*Hr5(0,-1,0,-1,0) + (11956.148148148148_dp)*Hr5(0,-1,0,-1,0) + &
                        (-9296.987654320988_dp)*Hr5(0,-1,0,0,0) + z**2*(2714.864197530864_dp)*Hr5(0,-1,0,0,0) + &
                        z*(41260.246913580246_dp)*Hr5(0,-1,0,0,0) + (-15934.814814814816_dp)*Hr5(0,-1,0,0,1) + &
                        z**2*(4408.888888888889_dp)*Hr5(0,-1,0,0,1) + z*(70352.5925925926_dp)*Hr5(0,-1,0,0,1) + &
                        (-6637.827160493827_dp)*Hr5(0,-1,0,1,0) + z**2*(1694.0246913580247_dp)*Hr5(0,-1,0,1,0) + &
                        z*(29092.345679012345_dp)*Hr5(0,-1,0,1,0) + (-13275.654320987655_dp)*Hr5(0,-1,0,1,1) + &
                        z**2*(3388.0493827160494_dp)*Hr5(0,-1,0,1,1) + z*(58184.69135802469_dp)*Hr5(0,-1,0,1,1) + &
                        z*(-97697.18518518518_dp)*Hr5(0,0,-1,-1,0) + (-10178.37037037037_dp)*Hr5(0,0,-1,-1,0) + &
                        z**2*(-910.2222222222222_dp)*Hr5(0,0,-1,-1,0) + nf*(676.3456790123457_dp)*Hr5(0,0,-1,-1,0) + &
                        nf*z*(1352.6913580246915_dp)*Hr5(0,0,-1,-1,0) + nf*z*(-1231.0123456790122_dp)*Hr5(0,0,-1,0,0) + &
                        nf*(-615.5061728395061_dp)*Hr5(0,0,-1,0,0) + z**2*(764.8395061728395_dp)*Hr5(0,0,-1,0,0) + &
                        (9068.04938271605_dp)*Hr5(0,0,-1,0,0) + z*(86733.43209876544_dp)*Hr5(0,0,-1,0,0) + &
                        nf*z*(-1109.3333333333333_dp)*Hr5(0,0,-1,0,1) + nf*(-554.6666666666666_dp)*Hr5(0,0,-1,0,1) + &
                        z**2*(619.4567901234568_dp)*Hr5(0,0,-1,0,1) + (7957.728395061728_dp)*Hr5(0,0,-1,0,1) + &
                        z*(75769.67901234567_dp)*Hr5(0,0,-1,0,1) + (-2970.864197530864_dp)*Hr5(0,0,0,-1,0) + &
                        nf*z*(-2212.3456790123455_dp)*Hr5(0,0,0,-1,0) + z**2*(101.1358024691358_dp)*Hr5(0,0,0,-1,0) + &
                        nf*(429.82716049382714_dp)*Hr5(0,0,0,-1,0) + z*(60049.38271604938_dp)*Hr5(0,0,0,-1,0) + &
                        nf*z*(-4917.728395061728_dp)*Hr5(0,0,0,0,0) + (-2829.1687242798353_dp)*Hr5(0,0,0,0,0) + &
                        z**2*(-269.69547325102883_dp)*Hr5(0,0,0,0,0) + nf**2*(-18.962962962962962_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z*(37.925925925925924_dp)*Hr5(0,0,0,0,0) + nf*(578.3703703703703_dp)*Hr5(0,0,0,0,0) + &
                        z*(64364.510288065845_dp)*Hr5(0,0,0,0,0) + nf*z*(-5685.728395061728_dp)*Hr5(0,0,0,0,1) + &
                        z**2*(-1331.6213991769548_dp)*Hr5(0,0,0,0,1) + nf*(-1253.1358024691358_dp)*Hr5(0,0,0,0,1) + &
                        (9941.069958847736_dp)*Hr5(0,0,0,0,1) + z*(135245.95884773662_dp)*Hr5(0,0,0,0,1) + &
                        (-9867.720164609054_dp)*Hr5(0,0,0,1,0) + nf*z*(-3151.0123456790125_dp)*Hr5(0,0,0,1,0) + &
                        z**2*(-3025.6460905349795_dp)*Hr5(0,0,0,1,0) + nf*(466.17283950617286_dp)*Hr5(0,0,0,1,0) + &
                        z*(148004.87242798353_dp)*Hr5(0,0,0,1,0) + (-25135.27572016461_dp)*Hr5(0,0,0,1,1) + &
                        z**2*(-5427.621399176955_dp)*Hr5(0,0,0,1,1) + nf*z*(-2762.2716049382716_dp)*Hr5(0,0,0,1,1) + &
                        nf*(1594.469135802469_dp)*Hr5(0,0,0,1,1) + z*(182400.5267489712_dp)*Hr5(0,0,0,1,1) + &
                        z**2*(-4482.633744855967_dp)*Hr5(0,0,1,0,0) + (-1781.3991769547324_dp)*Hr5(0,0,1,0,0) + &
                        nf*z*(-883.358024691358_dp)*Hr5(0,0,1,0,0) + nf*(441.679012345679_dp)*Hr5(0,0,1,0,0) + &
                        z*(100666.20576131687_dp)*Hr5(0,0,1,0,0) + z**2*(-10230.518518518518_dp)*Hr5(0,0,1,0,1) + &
                        nf*z*(-799.604938271605_dp)*Hr5(0,0,1,0,1) + nf*(399.8024691358025_dp)*Hr5(0,0,1,0,1) + &
                        (10186.271604938273_dp)*Hr5(0,0,1,0,1) + z*(146859.45679012345_dp)*Hr5(0,0,1,0,1) + &
                        z**2*(-10230.518518518518_dp)*Hr5(0,0,1,1,0) + nf*z*(-799.604938271605_dp)*Hr5(0,0,1,1,0) + &
                        nf*(399.8024691358025_dp)*Hr5(0,0,1,1,0) + (10186.271604938273_dp)*Hr5(0,0,1,1,0) + &
                        z*(146859.45679012345_dp)*Hr5(0,0,1,1,0) + z**2*(-16807.506172839505_dp)*Hr5(0,0,1,1,1) + &
                        nf*(-308.14814814814815_dp)*Hr5(0,0,1,1,1) + nf*z*(616.2962962962963_dp)*Hr5(0,0,1,1,1) + &
                        (39233.97530864197_dp)*Hr5(0,0,1,1,1) + z*(105688.49382716049_dp)*Hr5(0,0,1,1,1) + &
                        z**2*(-4928.263374485597_dp)*Hr5(0,1,0,0,0) + (7066.3374485596705_dp)*Hr5(0,1,0,0,0) + &
                        z*(37650.43621399177_dp)*Hr5(0,1,0,0,0) + z**2*(-12475.522633744857_dp)*Hr5(0,1,0,0,1) + &
                        (13766.584362139918_dp)*Hr5(0,1,0,0,1) + z*(81750.38683127573_dp)*Hr5(0,1,0,0,1) + &
                        z**2*(-15094.518518518518_dp)*Hr5(0,1,0,1,0) + (13400.493827160493_dp)*Hr5(0,1,0,1,0) + &
                        z*(88199.90123456791_dp)*Hr5(0,1,0,1,0) + z**2*(-19888.987654320987_dp)*Hr5(0,1,0,1,1) + &
                        (4324.740740740741_dp)*Hr5(0,1,0,1,1) + z*(71044.74074074074_dp)*Hr5(0,1,0,1,1) + &
                        z**2*(-12475.522633744857_dp)*Hr5(0,1,1,0,0) + (13766.584362139918_dp)*Hr5(0,1,1,0,0) + &
                        z*(81750.38683127573_dp)*Hr5(0,1,1,0,0) + z**2*(-19888.987654320987_dp)*Hr5(0,1,1,0,1) + &
                        (4324.740740740741_dp)*Hr5(0,1,1,0,1) + z*(71044.74074074074_dp)*Hr5(0,1,1,0,1) + &
                        z**2*(-19888.987654320987_dp)*Hr5(0,1,1,1,0) + (4324.740740740741_dp)*Hr5(0,1,1,1,0) + &
                        z*(71044.74074074074_dp)*Hr5(0,1,1,1,0) + z*(-68620.64197530864_dp)*Hr5(0,1,1,1,1) + &
                        (-36303.01234567901_dp)*Hr5(0,1,1,1,1) + z**2*(-19177.876543209877_dp)*Hr5(0,1,1,1,1) + &
                        z**2*(-5812.01646090535_dp)*Hr5(1,0,0,0,0) + (-2906.008230452675_dp)*Hr5(1,0,0,0,0) + &
                        z*(5812.01646090535_dp)*Hr5(1,0,0,0,0) + z**2*(-19200.263374485596_dp)*Hr5(1,0,0,0,1) + &
                        (-9600.131687242798_dp)*Hr5(1,0,0,0,1) + z*(19200.263374485596_dp)*Hr5(1,0,0,0,1) + &
                        z**2*(-28540.70781893004_dp)*Hr5(1,0,0,1,0) + (-14270.35390946502_dp)*Hr5(1,0,0,1,0) + &
                        z*(28540.70781893004_dp)*Hr5(1,0,0,1,0) + z**2*(-36947.22633744856_dp)*Hr5(1,0,0,1,1) + &
                        (-18473.61316872428_dp)*Hr5(1,0,0,1,1) + z*(36947.22633744856_dp)*Hr5(1,0,0,1,1) + &
                        z**2*(-28540.70781893004_dp)*Hr5(1,0,1,0,0) + (-14270.35390946502_dp)*Hr5(1,0,1,0,0) + &
                        z*(28540.70781893004_dp)*Hr5(1,0,1,0,0) + z**2*(-35493.92592592593_dp)*Hr5(1,0,1,0,1) + &
                        (-17746.962962962964_dp)*Hr5(1,0,1,0,1) + z*(35493.92592592593_dp)*Hr5(1,0,1,0,1) + &
                        z**2*(-35493.92592592593_dp)*Hr5(1,0,1,1,0) + (-17746.962962962964_dp)*Hr5(1,0,1,1,0) + &
                        z*(35493.92592592593_dp)*Hr5(1,0,1,1,0) + z**2*(-10585.283950617284_dp)*Hr5(1,0,1,1,1) + &
                        (-5292.641975308642_dp)*Hr5(1,0,1,1,1) + z*(10585.283950617284_dp)*Hr5(1,0,1,1,1) + &
                        z**2*(-19200.263374485596_dp)*Hr5(1,1,0,0,0) + (-9600.131687242798_dp)*Hr5(1,1,0,0,0) + &
                        z*(19200.263374485596_dp)*Hr5(1,1,0,0,0) + z**2*(-36947.22633744856_dp)*Hr5(1,1,0,0,1) + &
                        (-18473.61316872428_dp)*Hr5(1,1,0,0,1) + z*(36947.22633744856_dp)*Hr5(1,1,0,0,1) + &
                        z**2*(-35493.92592592593_dp)*Hr5(1,1,0,1,0) + (-17746.962962962964_dp)*Hr5(1,1,0,1,0) + &
                        z*(35493.92592592593_dp)*Hr5(1,1,0,1,0) + z**2*(-10585.283950617284_dp)*Hr5(1,1,0,1,1) + &
                        (-5292.641975308642_dp)*Hr5(1,1,0,1,1) + z*(10585.283950617284_dp)*Hr5(1,1,0,1,1) + &
                        z**2*(-36947.22633744856_dp)*Hr5(1,1,1,0,0) + (-18473.61316872428_dp)*Hr5(1,1,1,0,0) + &
                        z*(36947.22633744856_dp)*Hr5(1,1,1,0,0) + z**2*(-10585.283950617284_dp)*Hr5(1,1,1,0,1) + &
                        (-5292.641975308642_dp)*Hr5(1,1,1,0,1) + z*(10585.283950617284_dp)*Hr5(1,1,1,0,1) + &
                        z**2*(-10585.283950617284_dp)*Hr5(1,1,1,1,0) + (-5292.641975308642_dp)*Hr5(1,1,1,1,0) + &
                        z*(10585.283950617284_dp)*Hr5(1,1,1,1,0) + z*(-99634.56790123456_dp)*Hr5(1,1,1,1,1) + &
                        (49817.28395061728_dp)*Hr5(1,1,1,1,1) + z**2*(99634.56790123456_dp)*Hr5(1,1,1,1,1)

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
                        Ibar_select = z**2*(-151051.76825481033_dp) + (-41963.678091156915_dp) + nf*z*(-4320.549201653714_dp) + &
                        (nf*(-2152.5925663806047_dp))/z + nf**2*z**2*(-429.41735895976643_dp) + nf**2*(-84.11791050852364_dp) + &
                        nf**3*z*(-6.360493827160494_dp) + (nf**3*(-0.4345679012345679_dp))/z + nf**3*(3.1901234567901233_dp) + &
                        nf**3*z**2*(3.3679012345679014_dp) + (nf**2*(57.225788751714674_dp))/z + nf**2*z*(345.02184234551737_dp) + &
                        nf*(3353.306020454973_dp) + z*(9343.664372666477_dp) + nf*z**2*(15714.886815792368_dp) + &
                        (23775.64942456555_dp)/z + z**2*(-81942.84798341214_dp)*Hr1(0) + nf*z*(-6553.636528491377_dp)*Hr1(0) + &
                        nf*(-794.6760608299096_dp)*Hr1(0) + (nf*(-567.2296296296296_dp)*Hr1(0))/z + &
                        nf**2*z**2*(-20.20411522633745_dp)*Hr1(0) + nf**3*z*(-1.0271604938271606_dp)*Hr1(0) + &
                        nf**3*z**2*(-0.3160493827160494_dp)*Hr1(0) + nf**3*(1.0469135802469136_dp)*Hr1(0) + &
                        (nf**2*(8.61673525377229_dp)*Hr1(0))/z + nf**2*(18.99555339861041_dp)*Hr1(0) + &
                        nf**2*z*(251.54963394351992_dp)*Hr1(0) + (2608.5036620392884_dp)*Hr1(0) + &
                        nf*z**2*(3617.7445477089896_dp)*Hr1(0) + ((8664.193757773708_dp)*Hr1(0))/z + z*(41063.17669416033_dp)*Hr1(0) &
                        + z**2*(-132489.89833336763_dp)*Hr1(1) + (-18210.525900760178_dp)*Hr1(1) + &
                        nf*z*(-9720.858089387748_dp)*Hr1(1) + (nf*(-1314.0411522633744_dp)*Hr1(1))/z + &
                        nf**2*z**2*(-183.2647462277092_dp)*Hr1(1) + nf**2*(-86.43621399176955_dp)*Hr1(1) + &
                        nf**3*z*(-0.9876543209876543_dp)*Hr1(1) + nf**3*(0.49382716049382713_dp)*Hr1(1) + &
                        nf**3*z**2*(0.9876543209876543_dp)*Hr1(1) + (nf**2*(17.338820301783265_dp)*Hr1(1))/z + &
                        nf**2*z*(235.22633744855966_dp)*Hr1(1) + nf*(1327.3261640354376_dp)*Hr1(1) + &
                        nf*z**2*(9347.12969432602_dp)*Hr1(1) + ((22016.227345333817_dp)*Hr1(1))/z + z*(137760.63058796513_dp)*Hr1(1) &
                        + z**2*(-12920.920783172625_dp)*Hr2(0,0) + nf*z*(-4859.849003535628_dp)*Hr2(0,0) + &
                        (-4454.7813849307995_dp)*Hr2(0,0) + (nf*(-55.46666666666667_dp)*Hr2(0,0))/z + &
                        nf**2*(-27.805761316872427_dp)*Hr2(0,0) + nf**3*z*(-1.4222222222222223_dp)*Hr2(0,0) + &
                        nf**2*z**2*(0.4213991769547325_dp)*Hr2(0,0) + nf**3*(0.7111111111111111_dp)*Hr2(0,0) + &
                        nf**2*z*(123.0485596707819_dp)*Hr2(0,0) + nf*z**2*(249.61316872427983_dp)*Hr2(0,0) + &
                        nf*(573.0633282516694_dp)*Hr2(0,0) + ((1603.2_dp)*Hr2(0,0))/z + z*(47239.56579301484_dp)*Hr2(0,0) + &
                        z**2*(-40508.0050738305_dp)*Hr2(0,1) + nf*z*(-6367.20987654321_dp)*Hr2(0,1) + &
                        nf*(-453.1358024691358_dp)*Hr2(0,1) + (nf*(-179.35802469135803_dp)*Hr2(0,1))/z + &
                        nf**2*(-11.390946502057613_dp)*Hr2(0,1) + nf**2*z**2*(0.5267489711934157_dp)*Hr2(0,1) + &
                        nf**2*z*(73.7448559670782_dp)*Hr2(0,1) + (528.7718483085318_dp)*Hr2(0,1) + &
                        nf*z**2*(969.5473251028807_dp)*Hr2(0,1) + ((5600._dp)*Hr2(0,1))/z + z*(85746.90570356837_dp)*Hr2(0,1) + &
                        z**2*(-52489.414221795836_dp)*Hr2(1,0) + (-6693.645382502857_dp)*Hr2(1,0) + &
                        nf*z*(-2438.9135802469136_dp)*Hr2(1,0) + (nf*(-195.6872427983539_dp)*Hr2(1,0))/z + &
                        nf**2*z**2*(-25.744855967078188_dp)*Hr2(1,0) + nf**2*(-12.872427983539094_dp)*Hr2(1,0) + &
                        nf**2*z*(25.744855967078188_dp)*Hr2(1,0) + nf*(631.0123456790124_dp)*Hr2(1,0) + &
                        nf*z**2*(2366.1563786008232_dp)*Hr2(1,0) + ((6845.62962962963_dp)*Hr2(1,0))/z + &
                        z*(51803.09323414152_dp)*Hr2(1,0) + ((-51.489711934156375_dp)*(z**2*(-1934.0571425582775_dp) + &
                        (-265.9028132992328_dp) + nf**2*z*(z*(-1._dp) + (0.5_dp) + z**2*(1._dp)) + nf*(z**3*(-91.90792838874681_dp) &
                        + z*(-24.510230179028135_dp) + (7.601023017902814_dp) + z**2*(94.73401534526855_dp)) + &
                        z*(220.94161475739952_dp) + z**3*(1960.7157103332133_dp))*Hr2(1,1))/z + (z**2*(-754.2913580246914_dp) + &
                        nf*(z*(-1937.0666666666666_dp) + (-153.6_dp) + z**2*(5.05679012345679_dp)) + nf**2*((5.530864197530864_dp) + &
                        z*(43.930864197530866_dp)) + (115.2_dp)/z + (892.3917695473251_dp) + z*(26513.093004115228_dp))*Hr3(0,0,0) + &
                        z**2*(-3439.4074074074074_dp)*Hr3(0,0,1) + (-3078.2551440329216_dp)*Hr3(0,0,1) + &
                        nf*z*(-2396.4444444444443_dp)*Hr3(0,0,1) + nf**2*(-15.012345679012345_dp)*Hr3(0,0,1) + &
                        nf*z**2*(25.28395061728395_dp)*Hr3(0,0,1) + nf**2*z*(30.02469135802469_dp)*Hr3(0,0,1) + &
                        nf*(334.22222222222223_dp)*Hr3(0,0,1) + ((512._dp)*Hr3(0,0,1))/z + z*(47178.27160493827_dp)*Hr3(0,0,1) + &
                        z**2*(-6820.740740740741_dp)*Hr3(0,1,0) + nf*z*(-1073.7777777777778_dp)*Hr3(0,1,0) + &
                        nf*(-103.11111111111111_dp)*Hr3(0,1,0) + nf*z**2*(74.27160493827161_dp)*Hr3(0,1,0) + &
                        ((867.5555555555555_dp)*Hr3(0,1,0))/z + (3402.6666666666665_dp)*Hr3(0,1,0) + &
                        z*(30138.074074074073_dp)*Hr3(0,1,0) + z**2*(-13641.481481481482_dp)*Hr3(0,1,1) + &
                        nf*z*(-2147.5555555555557_dp)*Hr3(0,1,1) + nf*(-206.22222222222223_dp)*Hr3(0,1,1) + &
                        nf*z**2*(148.54320987654322_dp)*Hr3(0,1,1) + ((1735.111111111111_dp)*Hr3(0,1,1))/z + &
                        (6805.333333333333_dp)*Hr3(0,1,1) + z*(60276.148148148146_dp)*Hr3(0,1,1) + &
                        z**2*(-8830.024691358025_dp)*Hr3(1,0,0) + (-1370.5349794238682_dp)*Hr3(1,0,0) + &
                        nf*z*(-220.44444444444446_dp)*Hr3(1,0,0) + nf*(110.22222222222223_dp)*Hr3(1,0,0) + &
                        nf*z**2*(220.44444444444446_dp)*Hr3(1,0,0) + ((670.0246913580247_dp)*Hr3(1,0,0))/z + &
                        z*(7875.028806584362_dp)*Hr3(1,0,0) + z**2*(-26490.074074074073_dp)*Hr3(1,0,1) + &
                        (-4111.604938271605_dp)*Hr3(1,0,1) + nf*z*(-661.3333333333334_dp)*Hr3(1,0,1) + &
                        nf*(330.6666666666667_dp)*Hr3(1,0,1) + nf*z**2*(661.3333333333334_dp)*Hr3(1,0,1) + &
                        ((2010.0740740740741_dp)*Hr3(1,0,1))/z + z*(23625.086419753086_dp)*Hr3(1,0,1) + &
                        z**2*(-26490.074074074073_dp)*Hr3(1,1,0) + (-4111.604938271605_dp)*Hr3(1,1,0) + &
                        nf*z*(-661.3333333333334_dp)*Hr3(1,1,0) + nf*(330.6666666666667_dp)*Hr3(1,1,0) + &
                        nf*z**2*(661.3333333333334_dp)*Hr3(1,1,0) + ((2010.0740740740741_dp)*Hr3(1,1,0))/z + &
                        z*(23625.086419753086_dp)*Hr3(1,1,0) + z**2*(-52980.148148148146_dp)*Hr3(1,1,1) + &
                        (-8223.20987654321_dp)*Hr3(1,1,1) + nf*z*(-1322.6666666666667_dp)*Hr3(1,1,1) + &
                        nf*(661.3333333333334_dp)*Hr3(1,1,1) + nf*z**2*(1322.6666666666667_dp)*Hr3(1,1,1) + &
                        ((4020.1481481481483_dp)*Hr3(1,1,1))/z + z*(47250.17283950617_dp)*Hr3(1,1,1) + &
                        nf*z*(-356.5037037037037_dp)*Hr4(0,0,0,0) + (-141.45843621399177_dp)*Hr4(0,0,0,0) + &
                        z**2*(-13.48477366255144_dp)*Hr4(0,0,0,0) + nf**2*(-0.9481481481481482_dp)*Hr4(0,0,0,0) + &
                        nf**2*z*(1.8962962962962964_dp)*Hr4(0,0,0,0) + nf*(28.918518518518518_dp)*Hr4(0,0,0,0) + &
                        z*(6220.694650205762_dp)*Hr4(0,0,0,0) + nf*z*(-436.14814814814815_dp)*Hr4(0,0,0,1) + &
                        z**2*(-84.2798353909465_dp)*Hr4(0,0,0,1) + nf*(-80.5925925925926_dp)*Hr4(0,0,0,1) + &
                        (641.843621399177_dp)*Hr4(0,0,0,1) + z*(13606.979423868313_dp)*Hr4(0,0,0,1) + &
                        (-1125.1358024691358_dp)*Hr4(0,0,1,0) + z**2*(-240.19753086419752_dp)*Hr4(0,0,1,0) + &
                        nf*z*(-137.4814814814815_dp)*Hr4(0,0,1,0) + nf*(68.74074074074075_dp)*Hr4(0,0,1,0) + &
                        z*(11295.604938271605_dp)*Hr4(0,0,1,0) + (-2250.2716049382716_dp)*Hr4(0,0,1,1) + &
                        z**2*(-480.39506172839504_dp)*Hr4(0,0,1,1) + nf*z*(-274.962962962963_dp)*Hr4(0,0,1,1) + &
                        nf*(137.4814814814815_dp)*Hr4(0,0,1,1) + z*(22591.20987654321_dp)*Hr4(0,0,1,1) + &
                        z**2*(-435.0946502057613_dp)*Hr4(0,1,0,0) + (896.2633744855967_dp)*Hr4(0,1,0,0) + &
                        z*(4237.6954732510285_dp)*Hr4(0,1,0,0) + z**2*(-1305.283950617284_dp)*Hr4(0,1,0,1) + &
                        (2688.7901234567903_dp)*Hr4(0,1,0,1) + z*(12713.086419753086_dp)*Hr4(0,1,0,1) + &
                        z**2*(-1305.283950617284_dp)*Hr4(0,1,1,0) + (2688.7901234567903_dp)*Hr4(0,1,1,0) + &
                        z*(12713.086419753086_dp)*Hr4(0,1,1,0) + z**2*(-2610.567901234568_dp)*Hr4(0,1,1,1) + &
                        (5377.580246913581_dp)*Hr4(0,1,1,1) + z*(25426.172839506173_dp)*Hr4(0,1,1,1) + &
                        z**2*(-611.2921810699588_dp)*Hr4(1,0,0,0) + (-305.6460905349794_dp)*Hr4(1,0,0,0) + &
                        z*(611.2921810699588_dp)*Hr4(1,0,0,0) + z**2*(-2445.1687242798353_dp)*Hr4(1,0,0,1) + &
                        (-1222.5843621399176_dp)*Hr4(1,0,0,1) + z*(2445.1687242798353_dp)*Hr4(1,0,0,1) + &
                        z**2*(-3667.753086419753_dp)*Hr4(1,0,1,0) + (-1833.8765432098764_dp)*Hr4(1,0,1,0) + &
                        z*(3667.753086419753_dp)*Hr4(1,0,1,0) + z**2*(-7335.506172839506_dp)*Hr4(1,0,1,1) + &
                        (-3667.753086419753_dp)*Hr4(1,0,1,1) + z*(7335.506172839506_dp)*Hr4(1,0,1,1) + &
                        z**2*(-2445.1687242798353_dp)*Hr4(1,1,0,0) + (-1222.5843621399176_dp)*Hr4(1,1,0,0) + &
                        z*(2445.1687242798353_dp)*Hr4(1,1,0,0) + z**2*(-7335.506172839506_dp)*Hr4(1,1,0,1) + &
                        (-3667.753086419753_dp)*Hr4(1,1,0,1) + z*(7335.506172839506_dp)*Hr4(1,1,0,1) + &
                        z**2*(-7335.506172839506_dp)*Hr4(1,1,1,0) + (-3667.753086419753_dp)*Hr4(1,1,1,0) + &
                        z*(7335.506172839506_dp)*Hr4(1,1,1,0) + z**2*(-14671.012345679012_dp)*Hr4(1,1,1,1) + &
                        (-7335.506172839506_dp)*Hr4(1,1,1,1) + z*(14671.012345679012_dp)*Hr4(1,1,1,1)

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
                        Ibar_select = z**2*(-1.854778471941995e6_dp) + (-610227.6213618444_dp) + (nf*(-46377.89882109445_dp))/z + &
                        z*(-45106.9988960973_dp) + nf*z*(-43829.68606236862_dp) + nf**2*z**2*(-16320.397813739768_dp) + &
                        nf**2*(-3293.0133742115695_dp) + nf**3*z*(-292.98343059244337_dp) + (nf**3*(-39.255235482395975_dp))/z + &
                        nf**4*(-1.9006858710562413_dp) + nf**4*z**2*(-1.859716506630087_dp) + (nf**4*(0.2399634202103338_dp))/z + &
                        nf**4*z*(3.6521262002743486_dp) + nf**3*(87.47777539989266_dp) + nf**3*z**2*(302.00810066536843_dp) + &
                        (nf**2*(2125.7992968182834_dp))/z + nf**2*z*(6172.061960907413_dp) + nf*(77528.05906104947_dp) + &
                        nf*z**2*(332208.53787173564_dp) + (346833.14130914246_dp)/z + z**2*(-1.6194187808168097e6_dp)*Hr1(0) + &
                        nf*z*(-93645.37412185513_dp)*Hr1(0) + (-17375.811735904193_dp)*Hr1(0) + &
                        (nf*(-15436.400096134023_dp)*Hr1(0))/z + nf*(-11878.89676682071_dp)*Hr1(0) + &
                        nf**2*z**2*(-2738.2360619868764_dp)*Hr1(0) + nf**3*z*(-160.29840655156485_dp)*Hr1(0) + &
                        (nf**3*(-5.580612711476909_dp)*Hr1(0))/z + nf**4*(-0.640877914951989_dp)*Hr1(0) + &
                        nf**4*z**2*(0.21069958847736625_dp)*Hr1(0) + nf**4*z*(0.4916323731138546_dp)*Hr1(0) + &
                        nf**3*(1.8677217943009374_dp)*Hr1(0) + nf**3*z**2*(4.171315348270081_dp)*Hr1(0) + &
                        (nf**2*(522.8444444444444_dp)*Hr1(0))/z + nf**2*(816.1608658254386_dp)*Hr1(0) + &
                        nf**2*z*(7086.0288632354805_dp)*Hr1(0) + nf*z**2*(129642.1286643231_dp)*Hr1(0) + &
                        ((143714.15968221176_dp)*Hr1(0))/z + z*(346417.06659947825_dp)*Hr1(0) + &
                        z**2*(-2.1501961489460045e6_dp)*Hr1(1) + (-321041.8260492884_dp)*Hr1(1) + &
                        nf*z*(-247716.22196973802_dp)*Hr1(1) + (nf*(-37337.86759928072_dp)*Hr1(1))/z + &
                        nf**2*z**2*(-8550.736960374461_dp)*Hr1(1) + nf**2*(-1930.3341866344188_dp)*Hr1(1) + &
                        nf**3*z*(-158.62752629172383_dp)*Hr1(1) + (nf**3*(-10.677396738302088_dp)*Hr1(1))/z + &
                        nf**4*z**2*(-0.5091906721536351_dp)*Hr1(1) + nf**4*(-0.25459533607681756_dp)*Hr1(1) + &
                        nf**4*z*(0.5091906721536351_dp)*Hr1(1) + nf**3*(65.15153177869227_dp)*Hr1(1) + &
                        nf**3*z**2*(116.3564090839811_dp)*Hr1(1) + (nf**2*(1142.855784179241_dp)*Hr1(1))/z + &
                        nf**2*z*(9516.91784743436_dp)*Hr1(1) + nf*(31711.350694543267_dp)*Hr1(1) + &
                        nf*z**2*(241360.23194346018_dp)*Hr1(1) + ((361183.50285641535_dp)*Hr1(1))/z + &
                        z*(2.2089474504057355e6_dp)*Hr1(1) + z**2*(-362001.53790032986_dp)*Hr2(0,0) + &
                        nf*z*(-92141.36205679426_dp)*Hr2(0,0) + (-85859.11855607109_dp)*Hr2(0,0) + &
                        (nf*(-2625.8962962962964_dp)*Hr2(0,0))/z + nf**2*(-753.842991293259_dp)*Hr2(0,0) + &
                        nf**2*z**2*(-132.53662551440328_dp)*Hr2(0,0) + nf**3*z*(-92.41518061271148_dp)*Hr2(0,0) + &
                        nf**3*z**2*(-0.7023319615912208_dp)*Hr2(0,0) + nf**4*(-0.4213991769547325_dp)*Hr2(0,0) + &
                        nf**4*z*(0.842798353909465_dp)*Hr2(0,0) + nf**3*(27.57384545038866_dp)*Hr2(0,0) + &
                        (nf**2*(47.43374485596708_dp)*Hr2(0,0))/z + nf**2*z*(4755.684183149737_dp)*Hr2(0,0) + &
                        nf*(11775.287872884634_dp)*Hr2(0,0) + nf*z**2*(16579.31028427252_dp)*Hr2(0,0) + &
                        ((35013.39470654559_dp)*Hr2(0,0))/z + z*(552333.9756581552_dp)*Hr2(0,0) + &
                        z**2*(-1.0153671697048068e6_dp)*Hr2(0,1) + nf*z*(-160850.432957618_dp)*Hr2(0,1) + &
                        (-75794.66579019974_dp)*Hr2(0,1) + (nf*(-8735.841975308642_dp)*Hr2(0,1))/z + &
                        nf*(-5508.901980600944_dp)*Hr2(0,1) + nf**2*z**2*(-504.59844535893916_dp)*Hr2(0,1) + &
                        nf**3*z*(-39.950983081847276_dp)*Hr2(0,1) + nf**3*z**2*(-2.8561499771376315_dp)*Hr2(0,1) + &
                        nf**3*(11.626520347508002_dp)*Hr2(0,1) + nf**2*(140.99753086419753_dp)*Hr2(0,1) + &
                        (nf**2*(146.84298125285778_dp)*Hr2(0,1))/z + nf**2*z*(5669.465020576132_dp)*Hr2(0,1) + &
                        nf*z**2*(53270.32999075175_dp)*Hr2(0,1) + ((117988.60233914478_dp)*Hr2(0,1))/z + &
                        z*(1.3755551335704024e6_dp)*Hr2(0,1) + z**2*(-1.0819356009340223e6_dp)*Hr2(1,0) + &
                        (-151896.04859824275_dp)*Hr2(1,0) + nf*z*(-81724.59725783492_dp)*Hr2(1,0) + &
                        (nf*(-10024.0768175583_dp)*Hr2(1,0))/z + nf**2*z**2*(-1893.1614083219022_dp)*Hr2(1,0) + &
                        nf**2*(-690.721536351166_dp)*Hr2(1,0) + nf**3*z*(-14.588020118884316_dp)*Hr2(1,0) + &
                        nf**3*(7.294010059442158_dp)*Hr2(1,0) + nf**3*z**2*(14.588020118884316_dp)*Hr2(1,0) + &
                        (nf**2*(151.10379515317786_dp)*Hr2(1,0))/z + nf**2*z*(2121.371742112483_dp)*Hr2(1,0) + &
                        nf*(13104.529996087096_dp)*Hr2(1,0) + nf*z**2*(80667.6132614929_dp)*Hr2(1,0) + &
                        ((159916.01820455655_dp)*Hr2(1,0))/z + z*(1.0940193688928154e6_dp)*Hr2(1,0) + &
                        ((29.176040237768632_dp)*(z**3*(-67895.74801074543_dp) + z*(-9073.783091835396_dp) + nf**3*z*(z*(-1._dp) + &
                        (0.5_dp) + z**2*(1._dp)) + nf**2*(z**3*(-129.77507522567703_dp) + z*(-47.348545636910735_dp) + &
                        (10.358074222668003_dp) + z**2*(145.41875626880642_dp)) + nf*(z**2*(-5438.094323731427_dp) + &
                        (-687.1444332998998_dp) + z*(816.2688268606985_dp) + z**3*(5365.638706880875_dp)) + (10568.035786109822_dp) &
                        + z**2*(69308.89751515933_dp))*Hr2(1,1))/z + nf*z*(-51339.23452619498_dp)*Hr3(0,0,0) + &
                        z**2*(-34691.19864786722_dp)*Hr3(0,0,0) + nf*(-3352.044684079132_dp)*Hr3(0,0,0) + &
                        (nf*(-173.51111111111112_dp)*Hr3(0,0,0))/z + nf**3*z*(-31.885871056241427_dp)*Hr3(0,0,0) + &
                        nf**3*(-2.7039780521262005_dp)*Hr3(0,0,0) + nf**2*z**2*(-0.9364426154549611_dp)*Hr3(0,0,0) + &
                        nf**2*(194.4273612721275_dp)*Hr3(0,0,0) + nf*z**2*(764.858039932937_dp)*Hr3(0,0,0) + &
                        nf**2*z*(2084.784783628584_dp)*Hr3(0,0,0) + ((4473.6_dp)*Hr3(0,0,0))/z + (14117.068045815977_dp)*Hr3(0,0,0) &
                        + z*(410447.2849956735_dp)*Hr3(0,0,0) + z**2*(-142713.49135107239_dp)*Hr3(0,0,1) + &
                        nf*z*(-83598.03841076273_dp)*Hr3(0,0,1) + (-34893.62581854445_dp)*Hr3(0,0,1) + &
                        (nf*(-741.4518518518519_dp)*Hr3(0,0,1))/z + nf**2*(-447.7322359396434_dp)*Hr3(0,0,1) + &
                        nf**3*z*(-22.685322359396434_dp)*Hr3(0,0,1) + nf**2*z**2*(-1.9665294924554184_dp)*Hr3(0,0,1) + &
                        nf**3*(11.342661179698217_dp)*Hr3(0,0,1) + nf**2*z*(2172.066941015089_dp)*Hr3(0,0,1) + &
                        nf*z**2*(3350.416095107453_dp)*Hr3(0,0,1) + nf*(6362.867399254252_dp)*Hr3(0,0,1) + &
                        ((20038.4_dp)*Hr3(0,0,1))/z + z*(809043.5451352277_dp)*Hr3(0,0,1) + z**2*(-239911.84746160943_dp)*Hr3(0,1,0) &
                        + nf*z*(-47774.55144032922_dp)*Hr3(0,1,0) + nf*(-4499.753086419753_dp)*Hr3(0,1,0) + &
                        (nf*(-1196.7736625514403_dp)*Hr3(0,1,0))/z + nf**2*z**2*(-14.748971193415638_dp)*Hr3(0,1,0) + &
                        nf**2*(2.392318244170096_dp)*Hr3(0,1,0) + nf**2*z*(737.5363511659808_dp)*Hr3(0,1,0) + &
                        nf*z**2*(6759.798811156836_dp)*Hr3(0,1,0) + (32036.3343028865_dp)*Hr3(0,1,0) + &
                        ((34434.37037037037_dp)*Hr3(0,1,0))/z + z*(652812.5231273256_dp)*Hr3(0,1,0) + &
                        z**2*(-473703.7399156924_dp)*Hr3(0,1,1) + nf*z*(-95549.10288065844_dp)*Hr3(0,1,1) + &
                        nf*(-8999.506172839507_dp)*Hr3(0,1,1) + (nf*(-2393.5473251028807_dp)*Hr3(0,1,1))/z + &
                        nf**2*z**2*(-29.497942386831276_dp)*Hr3(0,1,1) + nf**2*(4.784636488340192_dp)*Hr3(0,1,1) + &
                        nf**2*z*(1475.0727023319616_dp)*Hr3(0,1,1) + nf*z**2*(13519.597622313671_dp)*Hr3(0,1,1) + &
                        (48354.75873856348_dp)*Hr3(0,1,1) + ((68868.74074074074_dp)*Hr3(0,1,1))/z + &
                        z*(1.2335734742745233e6_dp)*Hr3(0,1,1) + z**2*(-231793.1565887786_dp)*Hr3(1,0,0) + &
                        (-32373.07395053928_dp)*Hr3(1,0,0) + nf*z*(-11666.363054412437_dp)*Hr3(1,0,0) + &
                        (nf*(-868.3554336229233_dp)*Hr3(1,0,0))/z + nf**2*z**2*(-156.02011888431642_dp)*Hr3(1,0,0) + &
                        nf**2*(-78.01005944215821_dp)*Hr3(1,0,0) + nf**2*z*(156.02011888431642_dp)*Hr3(1,0,0) + &
                        nf*(3011.197073616827_dp)*Hr3(1,0,0) + nf*z**2*(11878.418533760097_dp)*Hr3(1,0,0) + &
                        ((27394.897119341564_dp)*Hr3(1,0,0))/z + z*(224061.91561941418_dp)*Hr3(1,0,0) + &
                        z**2*(-664917.6555820544_dp)*Hr3(1,0,1) + (-81888.3147594771_dp)*Hr3(1,0,1) + &
                        nf*z*(-34999.08916323731_dp)*Hr3(1,0,1) + (nf*(-2605.06630086877_dp)*Hr3(1,0,1))/z + &
                        nf**2*z**2*(-468.0603566529492_dp)*Hr3(1,0,1) + nf**2*(-234.0301783264746_dp)*Hr3(1,0,1) + &
                        nf**2*z*(468.0603566529492_dp)*Hr3(1,0,1) + nf*(9033.59122085048_dp)*Hr3(1,0,1) + &
                        nf*z**2*(35635.255601280296_dp)*Hr3(1,0,1) + ((82184.69135802469_dp)*Hr3(1,0,1))/z + &
                        z*(641723.9326739612_dp)*Hr3(1,0,1) + z**2*(-685225.5317049087_dp)*Hr3(1,1,0) + &
                        (-92042.25282090426_dp)*Hr3(1,1,0) + nf*z*(-34999.08916323731_dp)*Hr3(1,1,0) + &
                        (nf*(-2605.06630086877_dp)*Hr3(1,1,0))/z + nf**2*z**2*(-468.0603566529492_dp)*Hr3(1,1,0) + &
                        nf**2*(-234.0301783264746_dp)*Hr3(1,1,0) + nf**2*z*(468.0603566529492_dp)*Hr3(1,1,0) + &
                        nf*(9033.59122085048_dp)*Hr3(1,1,0) + nf*z**2*(35635.255601280296_dp)*Hr3(1,1,0) + &
                        ((82184.69135802469_dp)*Hr3(1,1,0))/z + z*(662031.8087968155_dp)*Hr3(1,1,0) + &
                        z**2*(-1.3298353111641088e6_dp)*Hr3(1,1,1) + (-163776.6295189542_dp)*Hr3(1,1,1) + &
                        nf*z*(-69998.17832647462_dp)*Hr3(1,1,1) + (nf*(-5210.13260173754_dp)*Hr3(1,1,1))/z + &
                        nf**2*z**2*(-936.1207133058984_dp)*Hr3(1,1,1) + nf**2*(-468.0603566529492_dp)*Hr3(1,1,1) + &
                        nf**2*z*(936.1207133058984_dp)*Hr3(1,1,1) + nf*(18067.18244170096_dp)*Hr3(1,1,1) + &
                        nf*z**2*(71270.51120256059_dp)*Hr3(1,1,1) + ((164369.38271604938_dp)*Hr3(1,1,1))/z + &
                        z*(1.2834478653479223e6_dp)*Hr3(1,1,1) + nf*z*(-15487.730772748057_dp)*Hr4(0,0,0,0) + &
                        (-2971.0485596707817_dp)*Hr4(0,0,0,0) + z**2*(-1312.4916323731138_dp)*Hr4(0,0,0,0) + &
                        nf**2*(-48.724279835390945_dp)*Hr4(0,0,0,0) + nf**3*z*(-2.528395061728395_dp)*Hr4(0,0,0,0) + &
                        nf**3*(1.2641975308641975_dp)*Hr4(0,0,0,0) + nf*z**2*(8.989849108367627_dp)*Hr4(0,0,0,0) + &
                        ((230.4_dp)*Hr4(0,0,0,0))/z + nf**2*z*(475.12757201646093_dp)*Hr4(0,0,0,0) + &
                        nf*(723.6448102423411_dp)*Hr4(0,0,0,0) + z*(157300.8973022405_dp)*Hr4(0,0,0,0) + &
                        nf*z*(-25887.558116140834_dp)*Hr4(0,0,0,1) + z**2*(-7092.297393689986_dp)*Hr4(0,0,0,1) + &
                        nf*(-2179.499954275263_dp)*Hr4(0,0,0,1) + nf*z**2*(51.31705532693187_dp)*Hr4(0,0,0,1) + &
                        nf**2*(81.64609053497942_dp)*Hr4(0,0,0,1) + nf**2*z*(540.8658436213992_dp)*Hr4(0,0,0,1) + &
                        ((1254.4_dp)*Hr4(0,0,0,1))/z + (13230.058344764517_dp)*Hr4(0,0,0,1) + z*(352833.5568358482_dp)*Hr4(0,0,0,1) &
                        + (-22519.645176040238_dp)*Hr4(0,0,1,0) + z**2*(-16505.942386831277_dp)*Hr4(0,0,1,0) + &
                        nf*z*(-16286.463648834018_dp)*Hr4(0,0,1,0) + nf**2*(-88.23045267489712_dp)*Hr4(0,0,1,0) + &
                        nf*z**2*(140.46639231824417_dp)*Hr4(0,0,1,0) + nf**2*z*(176.46090534979425_dp)*Hr4(0,0,1,0) + &
                        nf*(2317.6515775034295_dp)*Hr4(0,0,1,0) + ((2759.1111111111113_dp)*Hr4(0,0,1,0))/z + &
                        z*(305786.9080932785_dp)*Hr4(0,0,1,0) + (-45039.290352080476_dp)*Hr4(0,0,1,1) + &
                        z**2*(-33011.884773662554_dp)*Hr4(0,0,1,1) + nf*z*(-32572.927297668037_dp)*Hr4(0,0,1,1) + &
                        nf**2*(-176.46090534979425_dp)*Hr4(0,0,1,1) + nf*z**2*(280.93278463648835_dp)*Hr4(0,0,1,1) + &
                        nf**2*z*(352.9218106995885_dp)*Hr4(0,0,1,1) + nf*(4635.303155006859_dp)*Hr4(0,0,1,1) + &
                        ((5518.222222222223_dp)*Hr4(0,0,1,1))/z + z*(611573.816186557_dp)*Hr4(0,0,1,1) + &
                        z**2*(-23084.15912208505_dp)*Hr4(0,1,0,0) + nf*z*(-4942.719707361683_dp)*Hr4(0,1,0,0) + &
                        nf*(-701.249199817101_dp)*Hr4(0,1,0,0) + nf*z**2*(305.16323731138544_dp)*Hr4(0,1,0,0) + &
                        ((3075.1604938271603_dp)*Hr4(0,1,0,0))/z + (16768.0438957476_dp)*Hr4(0,1,0,0) + &
                        z*(125184.6145404664_dp)*Hr4(0,1,0,0) + z**2*(-69252.47736625515_dp)*Hr4(0,1,0,1) + &
                        nf*z*(-14828.159122085048_dp)*Hr4(0,1,0,1) + nf*(-2103.7475994513034_dp)*Hr4(0,1,0,1) + &
                        nf*z**2*(915.4897119341564_dp)*Hr4(0,1,0,1) + ((9225.481481481482_dp)*Hr4(0,1,0,1))/z + &
                        (50304.1316872428_dp)*Hr4(0,1,0,1) + z*(375553.8436213992_dp)*Hr4(0,1,0,1) + &
                        z**2*(-69252.47736625515_dp)*Hr4(0,1,1,0) + nf*z*(-14828.159122085048_dp)*Hr4(0,1,1,0) + &
                        nf*(-2103.7475994513034_dp)*Hr4(0,1,1,0) + nf*z**2*(915.4897119341564_dp)*Hr4(0,1,1,0) + &
                        ((9225.481481481482_dp)*Hr4(0,1,1,0))/z + (50304.1316872428_dp)*Hr4(0,1,1,0) + &
                        z*(375553.8436213992_dp)*Hr4(0,1,1,0) + z**2*(-138504.9547325103_dp)*Hr4(0,1,1,1) + &
                        nf*z*(-29656.318244170096_dp)*Hr4(0,1,1,1) + nf*(-4207.495198902607_dp)*Hr4(0,1,1,1) + &
                        nf*z**2*(1830.9794238683128_dp)*Hr4(0,1,1,1) + ((18450.962962962964_dp)*Hr4(0,1,1,1))/z + &
                        (100608.2633744856_dp)*Hr4(0,1,1,1) + z*(751107.6872427984_dp)*Hr4(0,1,1,1) + &
                        z**2*(-25548.422496570645_dp)*Hr4(1,0,0,0) + (-4807.930498399634_dp)*Hr4(1,0,0,0) + &
                        nf*z*(-727.5573845450389_dp)*Hr4(1,0,0,0) + nf*(363.77869227251944_dp)*Hr4(1,0,0,0) + &
                        nf*z**2*(727.5573845450389_dp)*Hr4(1,0,0,0) + ((1747.5775034293554_dp)*Hr4(1,0,0,0))/z + &
                        z*(22955.222679469593_dp)*Hr4(1,0,0,0) + z**2*(-102193.68998628258_dp)*Hr4(1,0,0,1) + &
                        (-19231.721993598538_dp)*Hr4(1,0,0,1) + nf*z*(-2910.2295381801555_dp)*Hr4(1,0,0,1) + &
                        nf*(1455.1147690900777_dp)*Hr4(1,0,0,1) + nf*z**2*(2910.2295381801555_dp)*Hr4(1,0,0,1) + &
                        ((6990.3100137174215_dp)*Hr4(1,0,0,1))/z + z*(91820.89071787837_dp)*Hr4(1,0,0,1) + &
                        z**2*(-153290.53497942386_dp)*Hr4(1,0,1,0) + (-28847.582990397805_dp)*Hr4(1,0,1,0) + &
                        nf*z*(-4365.344307270233_dp)*Hr4(1,0,1,0) + nf*(2182.6721536351165_dp)*Hr4(1,0,1,0) + &
                        nf*z**2*(4365.344307270233_dp)*Hr4(1,0,1,0) + ((10485.465020576132_dp)*Hr4(1,0,1,0))/z + &
                        z*(137731.33607681756_dp)*Hr4(1,0,1,0) + z**2*(-306581.0699588477_dp)*Hr4(1,0,1,1) + &
                        (-57695.16598079561_dp)*Hr4(1,0,1,1) + nf*z*(-8730.688614540466_dp)*Hr4(1,0,1,1) + &
                        nf*(4365.344307270233_dp)*Hr4(1,0,1,1) + nf*z**2*(8730.688614540466_dp)*Hr4(1,0,1,1) + &
                        ((20970.930041152264_dp)*Hr4(1,0,1,1))/z + z*(275462.6721536351_dp)*Hr4(1,0,1,1) + &
                        z**2*(-102193.68998628258_dp)*Hr4(1,1,0,0) + (-19231.721993598538_dp)*Hr4(1,1,0,0) + &
                        nf*z*(-2910.2295381801555_dp)*Hr4(1,1,0,0) + nf*(1455.1147690900777_dp)*Hr4(1,1,0,0) + &
                        nf*z**2*(2910.2295381801555_dp)*Hr4(1,1,0,0) + ((6990.3100137174215_dp)*Hr4(1,1,0,0))/z + &
                        z*(91820.89071787837_dp)*Hr4(1,1,0,0) + z**2*(-306581.0699588477_dp)*Hr4(1,1,0,1) + &
                        (-57695.16598079561_dp)*Hr4(1,1,0,1) + nf*z*(-8730.688614540466_dp)*Hr4(1,1,0,1) + &
                        nf*(4365.344307270233_dp)*Hr4(1,1,0,1) + nf*z**2*(8730.688614540466_dp)*Hr4(1,1,0,1) + &
                        ((20970.930041152264_dp)*Hr4(1,1,0,1))/z + z*(275462.6721536351_dp)*Hr4(1,1,0,1) + &
                        z**2*(-306581.0699588477_dp)*Hr4(1,1,1,0) + (-57695.16598079561_dp)*Hr4(1,1,1,0) + &
                        nf*z*(-8730.688614540466_dp)*Hr4(1,1,1,0) + nf*(4365.344307270233_dp)*Hr4(1,1,1,0) + &
                        nf*z**2*(8730.688614540466_dp)*Hr4(1,1,1,0) + ((20970.930041152264_dp)*Hr4(1,1,1,0))/z + &
                        z*(275462.6721536351_dp)*Hr4(1,1,1,0) + z**2*(-613162.1399176954_dp)*Hr4(1,1,1,1) + &
                        (-115390.33196159122_dp)*Hr4(1,1,1,1) + nf*z*(-17461.377229080932_dp)*Hr4(1,1,1,1) + &
                        nf*(8730.688614540466_dp)*Hr4(1,1,1,1) + nf*z**2*(17461.377229080932_dp)*Hr4(1,1,1,1) + &
                        ((41941.86008230453_dp)*Hr4(1,1,1,1))/z + z*(550925.3443072703_dp)*Hr4(1,1,1,1) + &
                        nf*z*(-1944.0548696844994_dp)*Hr5(0,0,0,0,0) + nf*(-74.7281207133059_dp)*Hr5(0,0,0,0,0) + &
                        z**2*(-11.986465477823502_dp)*Hr5(0,0,0,0,0) + nf**2*(4.424691358024691_dp)*Hr5(0,0,0,0,0) + &
                        nf**2*z*(25.28395061728395_dp)*Hr5(0,0,0,0,0) + (282.72958390489254_dp)*Hr5(0,0,0,0,0) + &
                        z*(24883.153177869226_dp)*Hr5(0,0,0,0,0) + nf*z*(-3426.8181069958846_dp)*Hr5(0,0,0,0,1) + &
                        (-1568.4770004572474_dp)*Hr5(0,0,0,0,1) + z**2*(-86.9018747142204_dp)*Hr5(0,0,0,0,1) + &
                        nf**2*(-8.217283950617285_dp)*Hr5(0,0,0,0,1) + nf**2*z*(16.43456790123457_dp)*Hr5(0,0,0,0,1) + &
                        nf*(291.1868312757202_dp)*Hr5(0,0,0,0,1) + z*(66873.05276634659_dp)*Hr5(0,0,0,0,1) + &
                        nf*z*(-2052.213991769547_dp)*Hr5(0,0,0,1,0) + nf*(-396.1152263374486_dp)*Hr5(0,0,0,1,0) + &
                        z**2*(-288.424325560128_dp)*Hr5(0,0,0,1,0) + (3526.467306812986_dp)*Hr5(0,0,0,1,0) + &
                        z*(72411.36168267032_dp)*Hr5(0,0,0,1,0) + nf*z*(-4104.427983539094_dp)*Hr5(0,0,0,1,1) + &
                        nf*(-792.2304526748972_dp)*Hr5(0,0,0,1,1) + z**2*(-576.848651120256_dp)*Hr5(0,0,0,1,1) + &
                        (7052.934613625972_dp)*Hr5(0,0,0,1,1) + z*(144822.72336534064_dp)*Hr5(0,0,0,1,1) + &
                        (-4057.7814357567445_dp)*Hr5(0,0,1,0,0) + z**2*(-600.2597165066301_dp)*Hr5(0,0,1,0,0) + &
                        nf*z*(-419.99451303155007_dp)*Hr5(0,0,1,0,0) + nf*(209.99725651577504_dp)*Hr5(0,0,1,0,0) + &
                        z*(39571.95793324189_dp)*Hr5(0,0,1,0,0) + (-12173.344307270234_dp)*Hr5(0,0,1,0,1) + &
                        z**2*(-1800.7791495198903_dp)*Hr5(0,0,1,0,1) + nf*z*(-1259.9835390946503_dp)*Hr5(0,0,1,0,1) + &
                        nf*(629.9917695473251_dp)*Hr5(0,0,1,0,1) + z*(118715.87379972565_dp)*Hr5(0,0,1,0,1) + &
                        (-12173.344307270234_dp)*Hr5(0,0,1,1,0) + z**2*(-1800.7791495198903_dp)*Hr5(0,0,1,1,0) + &
                        nf*z*(-1259.9835390946503_dp)*Hr5(0,0,1,1,0) + nf*(629.9917695473251_dp)*Hr5(0,0,1,1,0) + &
                        z*(118715.87379972565_dp)*Hr5(0,0,1,1,0) + (-24346.688614540468_dp)*Hr5(0,0,1,1,1) + &
                        z**2*(-3601.5582990397806_dp)*Hr5(0,0,1,1,1) + nf*z*(-2519.9670781893005_dp)*Hr5(0,0,1,1,1) + &
                        nf*(1259.9835390946503_dp)*Hr5(0,0,1,1,1) + z*(237431.7475994513_dp)*Hr5(0,0,1,1,1) + &
                        z**2*(-930.1216278006401_dp)*Hr5(0,1,0,0,0) + (2388.835848193873_dp)*Hr5(0,1,0,0,0) + &
                        z*(10950.525834476452_dp)*Hr5(0,1,0,0,0) + z**2*(-3720.4865112025605_dp)*Hr5(0,1,0,0,1) + &
                        (9555.343392775492_dp)*Hr5(0,1,0,0,1) + z*(43802.103337905806_dp)*Hr5(0,1,0,0,1) + &
                        z**2*(-5580.729766803841_dp)*Hr5(0,1,0,1,0) + (14333.015089163237_dp)*Hr5(0,1,0,1,0) + &
                        z*(65703.15500685871_dp)*Hr5(0,1,0,1,0) + z**2*(-11161.459533607682_dp)*Hr5(0,1,0,1,1) + &
                        (28666.030178326473_dp)*Hr5(0,1,0,1,1) + z*(131406.31001371742_dp)*Hr5(0,1,0,1,1) + &
                        z**2*(-3720.4865112025605_dp)*Hr5(0,1,1,0,0) + (9555.343392775492_dp)*Hr5(0,1,1,0,0) + &
                        z*(43802.103337905806_dp)*Hr5(0,1,1,0,0) + z**2*(-11161.459533607682_dp)*Hr5(0,1,1,0,1) + &
                        (28666.030178326473_dp)*Hr5(0,1,1,0,1) + z*(131406.31001371742_dp)*Hr5(0,1,1,0,1) + &
                        z**2*(-11161.459533607682_dp)*Hr5(0,1,1,1,0) + (28666.030178326473_dp)*Hr5(0,1,1,1,0) + &
                        z*(131406.31001371742_dp)*Hr5(0,1,1,1,0) + z**2*(-22322.919067215364_dp)*Hr5(0,1,1,1,1) + &
                        (57332.060356652946_dp)*Hr5(0,1,1,1,1) + z*(262812.62002743484_dp)*Hr5(0,1,1,1,1) + &
                        z**2*(-1234.570827617741_dp)*Hr5(1,0,0,0,0) + (-617.2854138088705_dp)*Hr5(1,0,0,0,0) + &
                        z*(1234.570827617741_dp)*Hr5(1,0,0,0,0) + z**2*(-6172.854138088706_dp)*Hr5(1,0,0,0,1) + &
                        (-3086.427069044353_dp)*Hr5(1,0,0,0,1) + z*(6172.854138088706_dp)*Hr5(1,0,0,0,1) + &
                        z**2*(-12345.708276177413_dp)*Hr5(1,0,0,1,0) + (-6172.854138088706_dp)*Hr5(1,0,0,1,0) + &
                        z*(12345.708276177413_dp)*Hr5(1,0,0,1,0) + z**2*(-24691.416552354825_dp)*Hr5(1,0,0,1,1) + &
                        (-12345.708276177413_dp)*Hr5(1,0,0,1,1) + z*(24691.416552354825_dp)*Hr5(1,0,0,1,1) + &
                        z**2*(-12345.708276177413_dp)*Hr5(1,0,1,0,0) + (-6172.854138088706_dp)*Hr5(1,0,1,0,0) + &
                        z*(12345.708276177413_dp)*Hr5(1,0,1,0,0) + z**2*(-37037.124828532236_dp)*Hr5(1,0,1,0,1) + &
                        (-18518.562414266118_dp)*Hr5(1,0,1,0,1) + z*(37037.124828532236_dp)*Hr5(1,0,1,0,1) + &
                        z**2*(-37037.124828532236_dp)*Hr5(1,0,1,1,0) + (-18518.562414266118_dp)*Hr5(1,0,1,1,0) + &
                        z*(37037.124828532236_dp)*Hr5(1,0,1,1,0) + z**2*(-74074.24965706447_dp)*Hr5(1,0,1,1,1) + &
                        (-37037.124828532236_dp)*Hr5(1,0,1,1,1) + z*(74074.24965706447_dp)*Hr5(1,0,1,1,1) + &
                        z**2*(-6172.854138088706_dp)*Hr5(1,1,0,0,0) + (-3086.427069044353_dp)*Hr5(1,1,0,0,0) + &
                        z*(6172.854138088706_dp)*Hr5(1,1,0,0,0) + z**2*(-24691.416552354825_dp)*Hr5(1,1,0,0,1) + &
                        (-12345.708276177413_dp)*Hr5(1,1,0,0,1) + z*(24691.416552354825_dp)*Hr5(1,1,0,0,1) + &
                        z**2*(-37037.124828532236_dp)*Hr5(1,1,0,1,0) + (-18518.562414266118_dp)*Hr5(1,1,0,1,0) + &
                        z*(37037.124828532236_dp)*Hr5(1,1,0,1,0) + (z**2*(-74074.24965706447_dp) + (-37037.124828532236_dp) + &
                        z*(74074.24965706447_dp))*Hr5(1,1,0,1,1) + z**2*(-24691.416552354825_dp)*Hr5(1,1,1,0,0) + &
                        (-12345.708276177413_dp)*Hr5(1,1,1,0,0) + z*(24691.416552354825_dp)*Hr5(1,1,1,0,0) + &
                        z**2*(-74074.24965706447_dp)*Hr5(1,1,1,0,1) + (-37037.124828532236_dp)*Hr5(1,1,1,0,1) + &
                        z*(74074.24965706447_dp)*Hr5(1,1,1,0,1) + z**2*(-74074.24965706447_dp)*Hr5(1,1,1,1,0) + &
                        (-37037.124828532236_dp)*Hr5(1,1,1,1,0) + z*(74074.24965706447_dp)*Hr5(1,1,1,1,0) + &
                        z**2*(-148148.49931412895_dp)*Hr5(1,1,1,1,1) + (-74074.24965706447_dp)*Hr5(1,1,1,1,1) + &
                        z*(148148.49931412895_dp)*Hr5(1,1,1,1,1)

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
