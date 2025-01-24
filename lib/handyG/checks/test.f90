
! These tests assume that GPLInfinity = 30

PROGRAM TEST
  use globals
  use utils
  use shuffle
  use maths_functions
  use mpl_module
  use gpl_module
  use chenreftest, only: do_chen_test
  use ttools
  implicit none
  integer :: i
  character(len=32) :: arg


  i = 1
  do
    call get_command_argument(i, arg)
    if (len_trim(arg) == 0) exit

    ! parse verbosity
    select case(trim(arg))
      case('-verb')
#ifdef DEBUG
        verb = readint(trim(arg),i)
#else
        call iprint("Argument -verb is not available, compile with --debug", 2)
#endif
      case('-polylog-test')
        tol = zero
        call do_poly_tests
      case('-mpl-test')
        tol = zero * 1.e5_prec
        call do_MPL_tests
      case('-gpl-test')
        tol = zero * 1.e5_prec
        call do_GPL_tests
      case('-chen-test')
        tol = 8.0e-7
        call set_options(mpldel=1e-10_prec, lidel=1e-10_prec)
        tests_successful = tests_successful .and. do_chen_test()

#if defined(HAVE_GINAC)
      case('-ginac-tests')
        tol = 8.0e-7
        call do_ginac_tests
      case('-speed-tests')
        call do_timing_tests(readint(trim(arg),i))
      case('-hw-tests')
        tol = 8.0e-7
        call do_high_weight_tests
#else
      case('-ginac-tests', '-speed-tests', '-hw-tests')
        call iprint("Argument "//trim(arg)//" is not available, compile with --with-ginac", 2)
#endif

#ifdef HAVE_GINAC
      case('-long-test')
        tol = 8.0e-7
        call do_long_test(readint(trim(arg),i))
#else
      case('-long-test')
        call iprint("Argument "//trim(arg)//" is not available, compile with --with-ginac", 2)
#endif
#ifdef DEBUG
      case('-report')
        verb = 1000
        tol = zero * 1.e5_prec
        call do_MPL_tests
        call do_GPL_tests
        tol = 8.0e-7
        tests_successful = tests_successful .and. do_chen_test()
#if defined(HAVE_GINAC)
        call do_ginac_tests
#endif
#else
      case('-report')
        call iprint("Argument -report is not available, compile with --debug", 2)
#endif

      case('--help','--h','-h','-help')
        call printhelp
        stop 0
      case default
        call iprint("Unknown argument "//trim(arg)//". Try -h to get help",2)

    end select
    i = i+1
  end do

  ! call do_shuffle_tests() ! put this somewhere else


  if(tests_successful) then
    call iprint('All tests passed. ', 0)
  else 
    call iprint('Some tests failed. ',2)
  end if

CONTAINS
   

  subroutine printhelp
    character(len=32) :: arg
    call get_command_argument(0,arg)

#ifdef DEBUG
    print*,"Usage: "//trim(arg)//" [-verb <n>] [opts]"
#else
    print*,"Usage: "//trim(arg)//" [opts]"
#endif
    print*,"Runs a set of tests for handyG"
    print*,""
    print*,"Possible tests are:"
    print*,"    -polylog-test     compares a few complex polylogs"
    print*,"    -mpl-test         performs tests on the series expansion of"
    print*,"                      convergent MPLs"
    print*,"    -gpl-test         tests GPLs and their reduction. This includes"
    print*,"                      real, ieps, and complex arguments"
    print*,"    -chen-test        compares all GPLs needed in [1811.06461] to"
    print*,"                      reference values"
#ifdef HAVE_GINAC
    print*,"    -ginac-tests      compare all GPLs needed by [1801.01033],"
    print*,"                      [1709.07435], and [1806.08241] to GiNaC"
    print*,"    -speed-tests <n>  compare the evaluation speed of the"
    print*,"                      aforementioned GPLs to GiNaC by averaging over"
    print*,"                      <n> evaluations"
    print*,"    -hw-tests         compares `random' GPLs with high weight to GiNaC"
    print*,"    -long-test <seed> compares many `random' GPLs with weight up to"
    print*,"                      four to GiNaC. Use <seed> as the random seed"
#endif
#ifdef DEBUG
    print*,"    -verb <n>         sets the verbosity level to <n>"
    print*,"    -report           performs a coverage test"
#endif

  end subroutine

  subroutine do_poly_tests()
    implicit none
    complex(kind=prec) :: ref, z

    z = 0.9_prec * exp((0._prec,1._prec))
    ref = ( 0.31714670838088083420005254335964850081059000145548_prec,0.90620928887768788007657902653829197810196658407346_prec)
    call check( polylog(2, z), ref)
    ref = ( 0.47248076790053852605402172732597392822155906295748_prec,0.78010930007060002834519910647340859560884976749326_prec)
    call check( polylog(5, z), ref)
    ref = ( 0.48593029219713200876281747476986815689740323314048_prec,0.75804436457088007970152497542119942774917319813883_prec)
    call check( polylog(10, z), ref)

    z = 0.99_prec * exp((0._prec,1._prec))
    ref = ( 0.32369036161935845024104448122222955899340391857079_prec,1.00324349284774627379699444852820582088595888924728_prec)
    call check( polylog(2, z), ref)
    ref = ( 0.51780997596604204699322616691391794080129634794226_prec,0.86049504083221113347441620727610317572983128003060_prec)
    call check( polylog(5, z), ref)
    ref = ( 0.53448415482661295991549353367903775505364004919109_prec,0.83392812071435332971141463249398007205963468167980_prec)
    call check( polylog(10, z), ref)

    z = 0.999_prec * exp((0._prec,1._prec))
    ref = ( 0.32409544945505588451045078881514007688743827135543_prec,1.01288825835849929172881742864723514784379050415565_prec)
    call check( polylog(2, z), ref)
    ref = ( 0.52231995926604156766783810033494360274235745715174_prec,0.86855495448337982380397597048602597837498452189533_prec)
    call check( polylog(5, z), ref)
    ref = ( 0.53933913318980396043645098080232673436372398744714_prec,0.84151728931442886718138858740874188686013946305784_prec)
    call check( polylog(10, z), ref)

    z = 0.9999_prec * exp((0._prec,1._prec))
    ref = ( 0.32413353539256702649540641848025485576423369942320_prec,1.01385205195042738164387173721328870295049412205747_prec)
    call check( polylog(2, z), ref)
    ref = ( 0.52277072515537604456567059944364018719426776187131_prec,0.86936115589767703671512280052476965418527852729045_prec)
    call check( polylog(5, z), ref)
    ref = ( 0.53982462693059649557938673066095705131895675203284_prec,0.84227621410312363023331730853369343023972319744457_prec)
    call check( polylog(10, z), ref)

    z = 0.99999_prec * exp((0._prec,1._prec))
    ref = ( 0.32413731983117050553258243423708701423605015249363_prec,1.01394842438972303272195829575740012923040998957878_prec)
    call check( polylog(2, z), ref)
    ref = ( 0.52281579941693425566978147340712469715561075153075_prec,0.86944177813618755451719373878155049326662610449505_prec)
    call check( polylog(5, z), ref)
    ref = ( 0.53987317626370393558176605497595292698896795484936_prec,0.84235210666127870713197117713665972968191528367673_prec)
    call check( polylog(10, z), ref)

    z = 1.00001_prec * exp((0._prec,1._prec))
    ref = ( 0.32413816022128717836771254536463031191720755084579_prec,1.01396984031625872927469992897403420152117604112105_prec)
    call check( polylog(2, z), ref)
    ref = ( 0.52282581586202930584572041341446950240136278312399_prec,0.86945969424096113443707836590455779449482967107881_prec)
    call check( polylog(5, z), ref)
    ref = ( 0.53988396500338278677493363054405590339711449183908_prec,0.84236897167615972350483622438537793816611895748722_prec)
    call check( polylog(10, z), ref)

    z = 0.999999_prec * exp((0._prec,1._prec))
    ref = ( 0.32413769803355298203649114917074410549435060587269_prec,1.01395806156436393326567566256155241345384396358214_prec)
    call check( polylog(2, z), ref)
    ref = ( 0.52282030681981336566989095922161219600778245313269_prec,0.86944984038100599719483571253130901162490098836993_prec)
    call check( polylog(5, z), ref)
    ref = ( 0.53987803119660494489725285096895119223918936255085_prec,0.84235969591788706954609196967489789002490184573847_prec)
    call check( polylog(10, z), ref)

    z = 0.9_prec * exp((0._prec,2._prec))
    ref = (-0.44452028964890283007349990328300389929008760996086_prec,0.66878743447993535281130919066420953251950351295302_prec)
    call check( polylog(2, z), ref)
    ref = (-0.38839535801823932606323111787411257099116023643543_prec,0.79888940673536917702425494564670039829270801271179_prec)
    call check( polylog(5, z), ref)
    ref = (-0.37503747568892401565030939135002191537632464837044_prec,0.81776617491393185281903993184823060507449761681675_prec)
    call check( polylog(10, z), ref)

    z = 0.99_prec * exp((0._prec,2._prec))
    ref = (-0.49145214297542995287792244783916287476339763066814_prec,0.72142557030218510441094486726609727591395278662659_prec)
    call check( polylog(2, z), ref)
    ref = (-0.42847005378062827950272801697834832300896766072327_prec,0.87664867860304046775940447885176343497500132376663_prec)
    call check( polylog(5, z), ref)
    ref = (-0.41259541382394926678059535801164943627251464211660_prec,0.89947635154294082643326725592105204892079432293970_prec)
    call check( polylog(10, z), ref)

    z = 0.999999_prec * exp((0._prec,2._prec))
    ref = (-0.49665806619812223943709229668650064888627387195443_prec,0.72714548006682757749406286437255854455016840984651_prec)
    call check( polylog(2, z), ref)
    ref = (-0.43293336180419752186586468054516010134686757646849_prec,0.88526440139529870422603576510716149051864398952584_prec)
    call check( polylog(5, z), ref)
    ref = (-0.41676869338400351335064877272829801919405402540781_prec,0.90855360420316407253701176433346888385813040933393_prec)
    call check( polylog(10, z), ref)
  end subroutine

  subroutine do_MPL_tests()
    complex(kind=prec) :: ref
    print*, 'doing MPL tests...'
    
    ref = (0.0226974036143118403367048874190211851348_prec,0.)
    call test_one_MPL((/ 1,1 /),cmplx([ 119._prec/377._prec,  35._prec / 102._prec ],kind=prec),ref, '1.1')

    ref = (2.31346156303083313941538625891818239608127E-4_prec,0.)
    call test_one_MPL((/ 1,1 /),cmplx([0.03_prec, 100._prec-30._prec*sqrt(11._prec) ], kind=prec),ref, '1.2')
    
    ref = (2.34617486908178608355021488772249510275409004880846912378356e-5_prec,0.)
    call test_one_MPL((/ 2,1,2 /),cmplx([0.03_prec, 100._prec-30._prec*sqrt(11._prec), 250._prec/9._prec + sqrt(6875._prec)/3._prec ], kind=prec),ref, '1.3')
    
    ref = (-0.06565799418838388900030511939327151809905663800733251335678_prec,0.)
    call test_one_MPL((/1, 1/), (/(-0.25_prec,0.),(-2._prec,0.) /), ref, '1.4')

    ref = (-0.0319989639656484125354488353677231989818697987217232462808_prec,0.)
    call test_one_MPL((/2, 1/), (/(-0.25_prec,0.),(-2._prec,0.) /), ref, '1.5')
  end subroutine do_MPL_tests

  subroutine do_GPL_tests()
    complex(kind=prec) :: ref, res
    real(kind=prec) :: z, xchen
    print*, 'doing GPL tests...'
    
    ref = (0.0819393734128676670196050250748145117414_prec,0.)
    call test_one_condensed((/ 1,1 /),(/ (1.3_prec,0.), (1.1_prec,0.) /),(0.4_prec,0.),2,ref,'2.1')
    
    ref = (0.0159279595253714815409530431586662402173_prec,0.)
    call test_one_condensed((/ 3,2 /),(/ (1.3_prec,0.), (1.1_prec,0.) /),(0.4_prec,0.),2,ref,'2.2')

    ref = (0.0020332632172573960663410568998748341439_prec,0.)
    call test_one_condensed((/ 4 /),(/ (0._prec,0.) /),(1.6_prec,0.),1,ref,'2.3')

    ! requires making convergent
    ref = (0.0959304167763934268889451723647597269243_prec,-0.882935179519785044294098901083860819793_prec)
    call test_one_flat(cmplx([0,1,3,2],kind=prec),ref,'2.4')

    ref = (0.00994794778992862285848386211876302459484_prec,0.0)
    call test_one_flat([ (0._prec,0.), (0._prec,0.), cmplx(-10._prec/3._prec,kind=prec), cmplx(-10._prec/3._prec,kind=prec), (1._prec,0.) ],ref,'2.5')

    ! requires hoelder convolution
    ref = (-0.01270994282825128082040127705282437854604_prec,0.0)
    call test_one_flat([ (0._prec,0.), cmplx(10._prec/3._prec,kind=prec), (1._prec,0.), cmplx(10._prec/3._prec,kind=prec), (1._prec,0.) ],ref,'2.6')

    ! here the tests from the mathematica nb start
    ! --------------------------------------------------------------------------

    z = 1._prec/200; xchen = 0.3_prec;

    ref = (-0.00501254182354428204309373895836778138658_prec,0.0)
    call test_one_flat(cmplx([1./z,1.0_prec],kind=prec),ref,'3.1')

    ref = (-0.001501126126267145650881556104825717433567_prec,0.0)
    call test_one_flat(cmplx([1./(xchen*z),1.0_prec],kind=prec),ref,'3.2')

    ref = (-7.50286081781062398521044689345496255877E-4_prec,0.0)
    call test_one_flat(cmplx([(1+sqrt(1-z**2))/(xchen*z),1.0_prec],kind=prec),ref,'3.3')

    ref = (0.0074335969894769006767497160582533005912_prec,0.0)
    call test_one_flat(cmplx([-1._prec/xchen,-1._prec/xchen,1._prec,1._prec,1._prec],kind=prec),ref,'3.4')
    
    ref = (-8.40378597494480117619893114520705717404E-6_prec,0.0)
    call test_one_flat(cmplx([-1._prec/xchen,0._prec,-1._prec/xchen,1._prec/(xchen*z),1.0_prec],kind=prec),ref,'3.5')
    
    ref = (0.492575584745015188956099120118429131062_prec,2.63892140547433046711653044177665570864_prec)
    call test_one_flat(cmplx([-1._prec,-1._prec,z,z,1._prec],kind=prec),ref,'3.6')
    
    ref = (-0.001531771317885848316868878120481266022195_prec,-4.50039113668843745227527662921010599E-4_prec)
    call test_one_flat(cmplx([0._prec,0._prec,(1._prec-sqrt(1-z**2))/(xchen*z), 1._prec/(xchen*z),1._prec],kind=prec),ref,'3.7')


    ! here the chen integral tests start
    ! ----------------------------------------------------------------------

    z = 1/100._prec
    ref = (-1.203972804325935992622746217761838502954_prec,0.0)
    call test_one_flat((/ (0._prec,0.), (0.3_prec,0.) /),ref,'4.1')

    ref = (10.6037962209567960211233327771880353832_prec,0.0)
    call test_one_flat( (/ (0._prec,0.), (0._prec,0.), (0.01_prec,0.) /),ref,'4.2')

    ref = (5.0429900651807654463744449762857296353E-4_prec,0.0)
    call test_one_flat((/(100._prec,0.), (1._prec,0.), (0.3_prec,0.) /),ref,'4.3')

    ref = (0.05630861877185110944118569266547272128_prec,0.0)
    call test_one_flat( (/ (1._prec,0.), (0._prec,0.), (0.01_prec,0.) /) ,ref,'4.4')

    ref = (0.0403289515087272334386799760345612761274_prec,0.29227096475688417627233177065883250904_prec)
    call test_one_flat([cmplx(z,sqrt(1._prec-z**2),kind=prec), (0.3_prec,0._prec)],ref,'4.5')

    ref = (2.52106450983403521007132254962389653479E-5_prec,0.)
    call test_one_flat(cmplx([1._prec/z,  (1+sqrt(1-z**2))/z, 1._prec], kind=prec),ref,'4.6')

    ref = (0.0556470547634805907802456701296765199105_prec,0.)
    call test_one_flat((/ (-1._prec,0.), (0.01_prec,0.), (0._prec,0.), (0.01_prec,0.) /),ref,'4.7')

    ref = (3.79489584700899518570534417045601144619E-5_prec,0.)
    call test_one_flat((/ (100._prec,0.), (100._prec,0.), (1._prec,0.), (0._prec,0.), (1._prec,0.) /),ref,'4.8')

    ref = (7.57428443370789293292663421167820669347E-5_prec,0.)
    call test_one_flat((/ (100._prec,0.), (1._prec,0.), (100._prec,0.), (0._prec,0.), (1._prec,0.) /),ref,'4.9')

    ref = (1.88019115860102242751826218248903485831_prec,2.50943413859806280406885100241705929527_prec)
    call test_one_flat((/ (0.01_prec,0.), (-1.0_prec,0.), (0.01_prec,0.), (1._prec,0.) /),ref,'4.10')

    ref = (-0.0125391083150538283259376060156453327489_prec,-0.0154142501684376766660782277248386124052_prec)
    call test_one_flat(cmplx([z, (1+sqrt(1-z**2))/z, z, 1._prec],kind=prec),ref,'4.11')


    ! Increase test coverage
    ref = (-1.203972804325935992622746217761838502953611_prec, pi)
    call test_one_flat((/ (0._prec,0.), (-0.3_prec,0.) /),ref,'E.1')

    ref = (-0.539404830685964989184608085018712655769250_prec, 1.0303768265243124637877433270311515319634364_prec)
    call test_one_flat([(0._prec,0._prec), (0.3_prec, 0.5_prec)],ref,'E.2')

    ref = (3.79489584700899518570534417045601144619E-5_prec,0.)
    call iprint('  testing GPL E.3 ...',-1)
    res = G((/100._prec, 100._prec, 1._prec, 0._prec, 1._prec/))
    call check(res,ref)

    call iprint('  testing GPL E.4 ...',-1)
    res = G((/100,100,1,0,1/))
    call check(res,ref)

    call iprint('  testing GPL E.5 ...',-1)
    res = G_superflatn(toinum((/100,100,1,0,1/)), 5)
    call check(res,ref)

    call iprint('  testing GPL E.6 ...',-1)
    res = G((/100._prec,100._prec,1._prec,0._prec/), 1._prec)
    call check(res,ref)

    call iprint('  testing GPL E.7 ...',-1)
    res = G((/(100._prec,0.),(100._prec,0.),(1._prec,0.),(0._prec,0.)/), (1._prec,0.))
    call check(res,ref)

    call iprint('  testing GPL E.8 ...',-1)
    res = G( (/1,1,1,1/) , (/ 100._prec, 100._prec, 1._prec, 0._prec /), 1._prec)
    call check(res,ref)

    !ref = cmplx(0.01592795952537145)
    ref = (0.0485035531168477720320998551314479928648_prec,-2.82326885118931135859594797186107474877_prec)
    call iprint('  testing GPL E.9 ...',-1)
    res = G( (/ 3, 2 /), toinum((/ 1.3_prec, 1.1_prec /)), toinum(4._prec) )
    call check(res,ref)

    call iprint('  testing GPL E.10 ...',-1)
    res = G( (/ 3, 2 /), (/ 1.3_prec, 1.1_prec /), 4._prec )
    call check(res,ref)

    call iprint('  testing GPL E.11 ...',-1)
    res = G( (/ 3, 2 /), (/ (1.3_prec,0.), (1.1_prec,0.) /), (4._prec,0.) )
    call check(res,ref)

    call iprint('  testing GPL E.12 ...',-1)
    ref = (0.107961231635965271810236308594279564021882373573100109007684973130605145011141_prec, 0.)
    res = G((/ 0._prec, 1.3_prec, 0._prec, 0._prec, 4._prec, 1.1_prec /))
    call check(res,ref)

    call iprint('  testing GPL E.13 ...',-1)
    ref = (-1.229846637025798984989235266565630713100_prec,-1.05705870670612913651142424798975000765_prec)
    res = G([inum(2._prec, +1), inum(7._prec, +1)], inum(5._prec, +1))
    call check(res,ref)

    call iprint('  testing GPL E.14 ...',-1)
    ref = (-1.229846637025798984989235266565630713100_prec,+1.05705870670612913651142424798975000765_prec)
    res = G([inum(2._prec, -1), inum(7._prec, +1)], inum(5._prec, +1))
    call check(res,ref)

    call iprint('  testing GPL E.15 ...',-1)
    ref = (0.2982254208688088675254638762780704094718_prec,0.)
    res = G([inum(2._prec, -1), inum(7._prec, +1)], inum(-5._prec, -1))
    call check(res,ref)

    call iprint('  testing GPL E.17 ...',-1)
    ref = (0.190800137777535619036913153766083992418_prec, 0.)
    res = G((/ 6, 1, 1 /))
    call check(res,ref)

    call iprint('  testing GPL E.18 ...',-1)
    ref = (0.058192342415778512650117048874978455691_prec, 0.)
    res = G((/ 6, 1, -1 /))
    call check(res,ref)

    ref = cmplx(log(0.5_prec), pi,kind=prec)
    call iprint('  testing GPL E.19a ...',-1)
    res = G((/ inum(2._prec, +1) /), inum(3._prec, +1))
    call check(res,ref)
    call iprint('  testing GPL E.19b ...',-1)
    res = G((/ inum(2._prec, -1) /), inum(3._prec, +1))
    call check(res,conjg(ref))
    call iprint('  testing GPL E.19c ...',-1)
    res = G((/ inum(2._prec, +1) /), inum(3._prec, -1))
    call check(res,ref)
    call iprint('  testing GPL E.19d ...',-1)
    res = G((/ inum(2._prec, -1) /), inum(3._prec, -1))
    call check(res,conjg(ref))

    call iprint('  testing GPL E.19e ...',-1)
    res = G((/ inum(-2._prec, +1) /), inum(-3._prec, +1))
    call check(res,conjg(ref))
    call iprint('  testing GPL E.19f ...',-1)
    res = G((/ inum(-2._prec, -1) /), inum(-3._prec, +1))
    call check(res,ref)
    call iprint('  testing GPL E.19g ...',-1)
    res = G((/ inum(-2._prec, +1) /), inum(-3._prec, -1))
    call check(res,conjg(ref))
    call iprint('  testing GPL E.19h ...',-1)
    res = G((/ inum(-2._prec, -1) /), inum(-3._prec, -1))
    call check(res,ref)


    ref = -(2.374395270272480200677499763071638424_prec, - 1.273806204919600530933131685580471698_prec)
    call iprint('  testing GPL E.20a ...',-1)
    res = G((/ izero, inum(2._prec, +1) /), inum(3._prec, +1))
    call check(res,ref)
    call iprint('  testing GPL E.20b ...',-1)
    res = G((/ izero, inum(2._prec, -1) /), inum(3._prec, +1))
    call check(res,conjg(ref))
    call iprint('  testing GPL E.20c ...',-1)
    res = G((/ izero, inum(2._prec, +1) /), inum(3._prec, -1))
    call check(res,ref)
    call iprint('  testing GPL E.20d ...',-1)
    res = G((/ izero, inum(2._prec, -1) /), inum(3._prec, -1))
    call check(res,conjg(ref))

    call iprint('  testing GPL E.20e ...',-1)
    res = G((/ izero, inum(-2._prec, +1) /), inum(-3._prec, +1))
    call check(res,conjg(ref))
    call iprint('  testing GPL E.20f ...',-1)
    res = G((/ izero, inum(-2._prec, -1) /), inum(-3._prec, +1))
    call check(res,ref)
    call iprint('  testing GPL E.20g ...',-1)
    res = G((/ izero, inum(-2._prec, +1) /), inum(-3._prec, -1))
    call check(res,conjg(ref))
    call iprint('  testing GPL E.20h ...',-1)
    res = G((/ izero, inum(-2._prec, -1) /), inum(-3._prec, -1))
    call check(res,ref)


    ! Thanks to Roman Zwicky and Ben Pullin for these tests
    ref = (0.392707112217551702879328061598355445_prec, - 1.274969448494380061814943180491890080_prec)
    call test_one_flat((/ (0._prec,0.), (1._prec,0.), (0._prec, 1.5_prec) /),ref,'F.1')

    ref = (2.09624324167194065961839363660174566_prec, 0.62644605179087454819421905065313983_prec)
    call test_one_flat( (/ (1._prec,0.), (-1._prec,0.), (2.5_prec, -2.4_prec) /) , ref, 'F.2')

    ref = conjg(ref)
    call test_one_flat( (/ (1._prec,0.), (-1._prec,0.), (2.5_prec, +2.4_prec) /) , ref, 'F.3')

    ref = (0.6874619495289224183167486785286777066_prec, -1.8261934916106546308783576818077830287_prec)
    call test_one_flat( (/ (-1._prec,0.), (1._prec,0.), (2.5_prec, +2.4_prec) /) , ref, 'F.4')

    ref = (0.320016770069038023050391296154549_prec , 0.064263450879286017225353602844105_prec)
    call test_one_flat( (/ (0._prec,0.), (1._prec,0.), (-1._prec,0.), (0.3_prec, -1.2_prec) /) , ref, 'F.5')

    ref = conjg(ref)
    call test_one_flat( (/ (0._prec,0.), (1._prec,0.), (-1._prec,0.), (0.3_prec, +1.2_prec) /) , ref, 'F.6')

    ref = (0.8382358254435272068734922352_prec, - 0.3702062785327198487992149341_prec)
    call test_one_flat( (/ (0._prec,0.), (1._prec,0.), (-1._prec,0.), (1.1_prec,2._prec) /) , ref, 'F.7')

    ref = (0.185156872427485220072774923047908301422_prec,-0.249989197161429773744237773045427322847_prec)
    call test_one_flat( (/ (1._prec,0.), (-1._prec,0.), (2.3_prec,-1._prec), (1.1_prec,-0.1_prec) /) , ref, 'F.8')

    ref = (-1.11161035333074623447215317094270858897_prec,-0.875225967273157437459269290416981432294_prec)
    call test_one_flat( (/ (1._prec,0.), (1.1_prec,-0.1_prec), (1._prec,0.), (2.3_prec,-1.2_prec), (3.4_prec,9.7_prec) /) , ref, 'F.9')

    ref = (-0.1109203962012759970680116512959180346924_prec,0.114779444698808046532782358409267479097_prec)
    call test_one_flat( (/ (0._prec, 0.) , (-0.1_prec, - 0.3_prec), (1._prec,0.), (0.5_prec, 0._prec) /) , ref, 'F.10')

    ref = (-0.108930339927322068475622176344923224046_prec,0.1122893584726586134006078500824852031925_prec)
    call test_one_flat( (/ (-0.005_prec, 0.) , (-0.1_prec, - 0.3_prec), (1._prec,0.), (0.5_prec, 0._prec)/) , ref, 'F.11')

    ref = (-1.00642162671475190543135817048062920163_prec,0.75972801214517468120067717498723544972_prec)
    call test_one_flat( (/ (-0.005_prec, 0.) , (-0.1_prec, - 0.3_prec), (1._prec,0.), (1.5_prec, 0._prec)/) , ref, 'F.12')

    ref = (-0.764543307991553066235740993950606803628_prec, 0.54861507469010653128533842086292369104_prec)
    call test_one_flat( (/ (-0.1_prec, - 0.3_prec), (1._prec,0.), (-0.005_prec, 0.) , (0.5_prec, 0._prec)/) , ref, 'F.13')

    ref = (4.52004218277686921073832347986485146573_prec,-1.812384329980915889835338149845575221705_prec)
    call test_one_flat( (/ (1.005_prec,0.), (0._prec,0.), (1.1_prec,0.3_prec), (1._prec,0.) /) , ref, 'F.14')

    ref = (1.6706939608077788063504210821631008732_prec,-0.631820518538505729899318406678510198266_prec)
    call test_one_flat( (/ (1.1_prec,0.), (0._prec,0.), (1.1_prec,0.3_prec), (1._prec,0.) /) , ref, 'F.10')

    ref = (0.198693810038777868937967610639933071867_prec,+0.72807783717224553483870820314753060868_prec)
    call test_one_flat( (/ (0._prec,0.), (1.0_prec,0._prec), (1.1_prec,0.), (1.1_prec,+0.3_prec) /) , ref, 'F.11')


    ref = (1.6706939608077788063504210821631008732_prec,-0.631820518538505729899318406678510198266_prec)
    call test_one_flat( (/ (1.1_prec,0.), (0._prec,0._prec), (1.1_prec,+0.3_prec), (1._prec,0.) /) , ref, 'F.12')

    ref = (0.1810339553848393655117844582129810543006_prec,+0.82543851493400141056822586864697825094_prec)
    call test_one_flat( (/ (0._prec,0.), (1.0_prec,0._prec), (1._prec,0.), (1.1_prec,+0.3_prec) /) , ref, 'F.13')


    ref = (-3.028608056170828558746818459366566689807_prec,-0.52999686156911157999452896083600882192_prec)
    call test_one_flat( (/ (1._prec,0._prec), (0.0_prec,0._prec), (0._prec,1._prec), (1.1_prec,0.0_prec) /) , ref, 'F.14')

    ! Here the branch cut matters, in Mathematica this is entered as G(1_-, 0, -I, -1.1)
    ref = (0.0538179677556747874824671812943783830465_prec,0.340519960719077282936661573786600733497_prec)
    call iprint('  testing GPL F.15 ...',-1)
    res = G((/ inum(1._prec,-1), izero, inum((0._prec,-1._prec), -1) /), inum(-1.1_prec, +1))
    call check(res,ref)

    ref = (0.0075181720252360930529649927100295277234_prec,0.009279206321129807482628652573319115651_prec)
    call test_one_flat( (/ (0._prec,0.), (0._prec,0.), (1.0_prec,0._prec), (0.1_prec,-0.1_prec), (0.11_prec,0._prec) /) , ref, 'F.16')

    ref = (0.0135493735561310633082925053080300174678_prec,+0.01851696361979639851211857931403696163_prec)
    call test_one_flat( (/ (0._prec,0.), (0._prec,0.), (1.0_prec,0._prec), (0.1_prec,-0.1_prec), (0.15_prec,0._prec) /) , ref, 'F.17')

    ref = (-2.41176903997917647290181947143140323970805e-4_prec,0.00129690184788114442363777169655122_prec)
    call test_one_flat( (/ (0._prec,0.), (1._prec,0.), (0.5_prec,-0.1_prec), (-1.3_prec,10.2_prec), (0.4_prec,0._prec) /) , ref, 'F.18')

    ref = (-0.003113755848404291707649322614093090379815_prec,-4.6914973893242872720303973859400498154E-4_prec)
    call test_one_flat( (/ (0._prec,0.), (1._prec,0.), (0.1_prec,-0.1_prec), (-1.3_prec,10.2_prec), (0.4_prec,0._prec) /) , ref, 'F.19')

    call clearcache

    ref = (0.6492441149301267712849185177612894261440132_prec,2.76692869057070821067245201199781913943_prec)
    call iprint('  testing GPL G.1 ...',-1)
    res = G([ toinum( (-1641._prec,4682._prec)/5000. ),&
              toinum( (-1187._prec,4827._prec)/5000. ),&
              toinum( 1069._prec/5000., +1_1 ), toinum( 787._prec / 5000., -1_1)], toinum(1._prec) )
    call check(res,ref, ttol=5.e-9_prec)



  end subroutine do_GPL_tests
  



#ifdef HAVE_GINAC

  function evalt(arr, what)
#if KINDREAL==16
    use ieps, only: inum2inum
#endif
    implicit none
    complex(kind=prec) :: arr(:), evalt
    complex(kind=8) :: geval
    integer what, i, l

    evalt =0.

    do i=1,size(arr)
      if (abs(arr(i)) .gt. 1.e10) then ! isnan?
        l = i-1
        goto 123
      endif
    enddo
    l = size(arr)

123 continue
    if (l==0) return

    if (what .eq. 0) then
      evalt = G(arr(1:l))
    elseif (what.eq.1) then
#if KINDREAL==16
      evalt = cmplx(geval(inum2inum(toinum(arr(1:l))),l),kind=prec)
#else
      evalt = geval(toinum(arr(1:l)),l)
#endif
    endif
  end function



  subroutine perform_ginacv(n, args)
    use maths_functions, only:clearcache
    complex(kind=prec) :: args(:,:)
    integer i,n
    character(len=40) :: msg
    call clearcache

    do i=1,size(args,1)
      write(msg,900) n,i
      call iprint(msg, -1)
      call check(evalt(args(i,:),0),evalt(args(i,:),1))
    enddo

900 format("   testing GPL ",I1,".",I4.4," ...")

  end subroutine

  subroutine do_ginac_tests
    use maths_functions, only:clearcache
    use gtestchen   , only: inichen   =>args
    use gtestmuone  , only: inimuone  =>args
    use gtestmuonenp, only: inimuonenp=>args
    implicit none
    tol = 6.0e-9
    call perform_ginacv( 6, inichen   ((0.3_prec,0.), (0.1_prec,0.)) )
    call perform_ginacv( 7, inimuone  ((0.5_prec,0.), (0.6_prec,0.)) )
    call perform_ginacv( 8, inimuonenp((0.3_prec,0.), (0.6_prec,0.)) )

  end subroutine

  subroutine do_one_speed_test(args, u, msg)
    use maths_functions, only:clearcache
    implicit none
    complex(kind=prec) :: args(:,:,:),res
    integer(kind=8) cstart, cend, count_rate
    real(kind=prec) :: time(2), ttime(2)
    integer i,j, u
    character(len=*) msg
    character(len=100) msg2
    do j=1,size(args,1)
      ! try function a bunch of times
      call system_clock(cstart, count_rate=count_rate)
      do i=1,size(args,3)
        res=evalt(args(j,:,i),0)
      enddo
      call system_clock(cend, count_rate=count_rate)
      time(1) = real(cend-cstart)/real(count_rate,kind=prec)/size(args,3)
      if (time(1).lt.zero) print*,j, cend-cstart,count_rate
      ttime(1) = ttime(1) + time(1)

      call system_clock(cstart, count_rate=count_rate)
      do i=1,size(args,3)
        res=evalt(args(j,:,i),1)
      enddo
      call system_clock(cend, count_rate=count_rate)
      time(2) = real(cend-cstart)/real(count_rate,kind=prec)/size(args,3)
      if (time(2).lt.zero) print*,j
      ttime(2) = ttime(2) + time(2)

      write(u,*) time

      write(msg2, 900) j, size(args,1), msg
      call iprint(msg2,-1)
    enddo
    write(msg2,901) size(args,1)/ttime(2)/1000., size(args,1)/ttime(1)/1000., int(ttime(2)/ttime(1))
    call iprint(msg2,4)

    ! Lets to another, fair comparison
    call system_clock(cstart, count_rate=count_rate)
    do j=1,size(args,1)
      res=evalt(args(j,:,1),0)
    enddo
    call system_clock(cend, count_rate=count_rate)
    ttime(1) = real(cend-cstart)/real(count_rate,kind=prec)

    call system_clock(cstart, count_rate=count_rate)
    do j=1,size(args,1)
      res=evalt(args(j,:,1),1)
    enddo
    call system_clock(cend, count_rate=count_rate)
    ttime(2) = real(cend-cstart)/real(count_rate,kind=prec)

    write(msg2,902) msg, size(args,1)/ttime(2)/1000., size(args,1)/ttime(1)/1000., int(ttime(2)/ttime(1))
    call iprint(msg2,4)

900 FORMAT('Evaluating function ',i4,'/',i4,' for ',a)
901 format(' using GiNaC at ',F9.2,'kG/s and handyG at ',F9.2,'kG/s (',I3,'x)')
902 format('Evaluating ',A,' using GiNaC at ',F9.2,'kG/s and handyG at ',F9.2,'kG/s (',I3,'x)')
  end subroutine

  subroutine do_timing_tests(n)
    use gtestchen   , only: inichen   =>args
    use gtestchenff , only: inichenff =>args
    use gtestmuone  , only: inimuone  =>args
    use gtestmuonenp, only: inimuonenp=>args
    implicit none
    integer, intent(in) :: n
    integer i
    complex(kind=prec) :: cargs( 1399,5,n)
    complex(kind=prec) :: fargs(  540,5,n)
    complex(kind=prec) :: pargs(  198,5,n)
    complex(kind=prec) :: nargs( 1733,5,n)
    real(kind=prec) :: z, x, y, w
    integer ranseed
    ranseed = 233123

    do i=1,n
      z = ran2(ranseed) / 2.
      x = ran2(ranseed)*(1-z) + z
      cargs(:,:,i) = inichen  (cmplx(x,kind=prec), cmplx(z,kind=prec))
      fargs(:,:,i) = inichenff(cmplx(x,kind=prec), cmplx(z,kind=prec))

      w = ran2(ranseed) ! 0<w<1
      z = ran2(ranseed) * (sqrt(1-w+w**2)-sqrt(w)) + sqrt(w)
      nargs(:,:,i) = inimuonenp(cmplx(w,kind=prec), cmplx(z,kind=prec))

      x = ran2(ranseed)
      y = ran2(ranseed)
      pargs(:,:,i) = inimuone(cmplx(x,kind=prec), cmplx(y,kind=prec))
    enddo

    cargs(1181,:,:)=1.e15
    fargs(367,:,:)=1.e15


    open(unit=9, file="stats.txt")
    write(9,*) "Chen form factor"
    call do_one_speed_test(fargs,9,"Chen FF")
    write(9,*) "Chen"
    call do_one_speed_test(cargs,9,"Chen")
    write(9,*) "MUonE-planar"
    call do_one_speed_test(pargs,9,"Muone planar")
    write(9,*) "MUonE-non-planar"
    call do_one_speed_test(nargs,9,"Muone non planar")
    close(unit=9)

  end subroutine

  subroutine fyshuffle(list, seed)
    real(kind=prec) :: list(:), tmp
    integer j,i, seed

    do i=size(list),2,-1
      j = int(ran2(seed) * i) + 1
      tmp = list(j)
      list(j) = list(i)
      list(i) = tmp
    enddo
  end subroutine

  function test_one_high_weight(length, fzero, funiq, alphabet, seed, what) result(res)
    integer length,seed, i, what
    integer nzero, nuniq, nsmall, c
    real(kind=prec) :: fzero, funiq, fsmall
    real(kind=prec) :: arguments(length), alphabet(10)
    complex(kind=prec) :: Geval, res
    nzero  = nint(fzero*length)
    nuniq  = nint(       funiq*length)
    nsmall = nint(fsmall*funiq*length)

    if (nzero+nuniq > length) nuniq = length - nzero

    !do i=1,nsmall
    !  alphabet(i) = 0.8_prec * ran2(seed)
    !enddo
    !do i=nsmall+1,nuniq
    !  alphabet(i) = 1.5_prec/(0.15_prec + ran2(seed))
    !enddo

    arguments(1:nzero) = 0._prec
    do i=1,nuniq
      arguments(nzero+i) = alphabet(i)
    enddo
    do i=nzero+nuniq,length
      c = int(ran2(seed)*nuniq)+1
      arguments(i) = alphabet(c)
    enddo

    call fyshuffle(arguments, seed)
    do while(arguments(length)<zero)
      call fyshuffle(arguments, seed)
    enddo
    write(9,*)"args",arguments
    if (what.eq.1) then
      res = G(arguments, 1._prec)
    elseif (what.eq.2) then
      res = Geval(toinum(cmplx([arguments, 1._prec],kind=prec)),length+1)
    endif
  end function


  subroutine do_high_weight_tests
    integer, parameter :: ntests=10
    real(kind=prec), parameter :: fzero = 0.28
    real(kind=prec), parameter :: funiq = 0.63
    real(kind=prec), parameter :: fsmal = 0.42
    integer(kind=8) cstart, cend, count_rate
    complex(kind=prec) :: res(ntests,2)
    real(kind=prec) :: alphabet(10)

    integer i,j
    integer seed, seedold
    seed = 123123

    open(unit=9, file="stats-hi.txt")
    alphabet = (/ 0.3, 0.6, 1.8, 5.3, 3.6, 8.3, 0.1, 0.4, 0.2, 10.2 /)
    do i=2,7
      seedold=seed
      do j=1,ntests
        call system_clock(cstart, count_rate=count_rate)
        res(j,1) = test_one_high_weight(i, fzero, funiq, alphabet, seed,1)
        call system_clock(cend, count_rate=count_rate)
        write(9,*) "time",i,j,real(cend-cstart)/real(count_rate,kind=prec)!/ntests
      enddo

      seed=seedold
      do j=1,ntests
        call system_clock(cstart, count_rate=count_rate)
        res(j,2) = test_one_high_weight(i, fzero, funiq, alphabet, seed,2)
        call system_clock(cend, count_rate=count_rate)
        write(9,*)"ginac",i,j,real(cend-cstart)/real(count_rate,kind=prec)!/ntests
      enddo

      write(9,*)"del",abs(res(:,1)-res(:,2)) / tol
      write(9,*)
      write(9,*)
    enddo
    close(unit=9)
  end subroutine


  SUBROUTINE DO_LONG_TEST(seed)
#if KINDREAL==16
    use ieps,only:inum2inum
#endif
    implicit none
    integer,parameter :: nzero = 10
    integer,parameter :: nieps = 50
    integer,parameter :: ncmpl = 50

    integer :: perweight(5) = (/ 1000, 5000, 100000, 500000, 500 /)
    real(kind=prec), parameter :: rrange = 5
    type(inum), dimension(nzero+nieps+ncmpl) :: basis
    type(inum), dimension(size(perweight)) :: args
    type(inum), parameter :: ione = inum((1._prec,0._prec), di0)
    integer i, j, w, seed, oldseed
    integer(kind=1) i0
    real(kind=prec) :: v, maxd
    complex(kind=prec) :: ans(2)
    complex(kind=8) :: geval
    character(len=80) :: msg
    maxd=0._prec

    tol = 0.01

    write (msg, "(A10,I0.5,A4)") "long-test-", seed, ".bin"
    open(unit=9, action='write', form='unformatted', file=trim(msg))
    basis(1:nzero) = izero

    do i=nzero+1,nzero+nieps
      if (ran2(seed).gt.0.5) then
        i0 = +1_1
      else
        i0 = -1_1
      endif
      v = 2*rrange*(ran2(seed) - 0.5)
      basis(i) = inum(cmplx(v,kind=prec), i0)
    enddo
    do i=nzero+nieps+1,nzero+nieps+ncmpl
      v = 2*rrange*(ran2(seed) - 0.5)
      basis(i) = toinum(v*exp(i_*2*pi*ran2(seed)))
    enddo

    do w=1,size(perweight)
      do i=1,perweight(w)
        oldseed = seed
        write(msg, 900) i, perweight(w), w
        call iprint(msg,-1)
        args(1:w) = basis((/ (1+int(size(basis)*ran2(seed)),j=1,w) /))
#if KINDREAL==16
        ans(1) = cmplx(geval(inum2inum([args(1:w),ione]), w+1),kind=prec)
#else
        ans(1) = cmplx(geval([args(1:w),ione], w+1),kind=prec)
#endif
        ans(2) = G(args(1:w),ione)
        if ((abs(ans(1)) .gt. 1.e10).or.(abs(ans(2)) .gt. 1.e10)) then
          write(msg,903) oldseed
          call iprint(msg,3)
          print*,"Args are ",args(1:w)
          cycle
        endif
        write(9) w,i,oldseed, abs(ans(1)-ans(2))
        flush(9)
        if(abs(ans(1)-ans(2)) > maxd) maxd = abs(ans(1)-ans(2))
        if(abs(ans(1)-ans(2)) > tol) goto 123
      enddo
      write(msg,901) maxd
      if (maxd < 1.e-3_prec) then
        call iprint(msg,0)
      else
        call iprint(msg,1)
      endif
      maxd=0.
    enddo

    close(9)
    tests_successful = tests_successful .and. .true.

    return
123 continue
    write(msg,902) abs(ans(1)-ans(2)), oldseed
    call iprint(msg,1)
    print*,"Offending G was",args(1:w)
    close(9)
    tests_successful = tests_successful .and. .false.

900 FORMAT('Testing ',i7,'/',i7,' GPLs with w=',i1)
901 FORMAT(' done. Largest delta = ',ES10.3)
902 FORMAT(' failed with delta = ',ES10.3, ' (seed=',I10,')')
903 FORMAT(' GiNaC problems with seed=',I10)


  END SUBROUTINE
      FUNCTION RAN2(randy)

            ! This is the usual "random"

      implicit none
      real(kind=prec) :: MINV,RAN2
      integer m,a,Qran,r,hi,lo,randy
      PARAMETER(M=2147483647,A=16807,Qran=127773,R=2836)
      PARAMETER(MINV=0.46566128752458e-09)
      HI = RANDY/Qran
      LO = MOD(RANDY,Qran)
      RANDY = A*LO - R*HI
      IF(RANDY.LE.0) RANDY = RANDY + M
      RAN2 = RANDY*MINV
      END FUNCTION RAN2




#endif

  ! subroutine do_shuffle_tests() 
  !   complex(kind=prec) :: v(2) = cmplx((/1,2/))
  !   complex(kind=prec) :: w(2) = cmplx((/3,4/))

  !   call print_matrix(shuffle_product(v,w))
  ! end subroutine do_shuffle_tests

END PROGRAM TEST

! In terminal kann man den exit code bekommen via echo $? 
