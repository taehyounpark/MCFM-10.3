
MODULE maths_functions
  use globals
  use ieps
  use utils, only: factorial, binom
  implicit none
  interface polylog
    module procedure polylog1, polylog2
  end interface polylog

  real(kind=prec), parameter :: zeta(2:10) = (/ 1.6449340668482264364724151666460251892189499012068_prec,  &
      1.2020569031595942853997381615114499907649862923405_prec, 1.0823232337111381915160036965411679027747509519187_prec,  &
      1.0369277551433699263313654864570341680570809195019_prec, 1.0173430619844491397145179297909205279018174900329_prec,  &
      1.0083492773819228268397975498497967595998635605652_prec, 1.0040773561979443393786852385086524652589607906499_prec,  &
      1.0020083928260822144178527692324120604856058513949_prec, 1.0009945751278180853371459589003190170060195315645_prec /)

  real(kind=prec), parameter :: DirichletBeta(2:10) = (/ 0.91596559417721901505460351493238411077414937428167_prec, &
    0.96894614625936938048363484584691860006954026768391_prec, 0.98894455174110533610842263322837782131586088706273_prec, &
    0.99615782807708806400631936863097528151139552938826_prec, 0.99868522221843813544160078786020654967836454612651_prec, &
    0.99955450789053990949634654989905898300218848194998_prec, 0.99984999024682965633806705924046378147600743300743_prec, &
    0.99994968418722008982135887329384752737274799691796_prec, 0.99998316402619687740554072995833414145685781649717_prec /)
  type el
    complex(kind=prec) :: c
    complex(kind=prec) ans
  end type el

  type(el) :: cache(PolyLogCacheSize(1),PolyLogCacheSize(2))
  integer :: plcachesize(PolyLogCacheSize(1)) = 0
CONTAINS

  FUNCTION naive_polylog(m,x) result(res)
    ! Computes the classical polylogarithm Li_m(x) using series representation up to order n
    integer :: m
    complex(kind=prec) :: x, res, del
    integer(kind=ikin) :: i
    res=0._prec

    i = 1
    del = 1._prec
    do while (abs(del) > zero)
      if(i**m.lt.0) return ! roll over
      if(abs(x**i).lt.1.e-250_prec) return
      del = x**i/i**m
      if (abs(del) < LiDelta) return
      res = res+del
      i = i+1
    end do
  END FUNCTION naive_polylog

  FUNCTION bernoullinumber(n)
    ! This returns the n-th Bernoulli number by computing all Bernoulli numbers
    ! up to the n-th recursively using the relation
    !   Sum[ Binomial[m+1, k] BernoulliB[k], {k,0,m} ] = 0
    ! for m > 0 (https://mathworld.wolfram.com/BernoulliNumber.html).
    ! Solving this for BernoulliB[m] results in
    !   BernoulliB[m] = - Sum[Binomial[m, k] BernoulliB[k] / (m-k+1), {k,0,m-1}]
    ! for m > 0 and BernoulliB[0] = 1.
    ! Care is taken to avoid multiple computation.
    integer, intent(in) :: n
    ! the cache is dynamic, the first entry is the initial size, the
    ! second the increment
    integer, parameter :: cachecontr(2) = (/ 20, 10 /)
    real(kind=prec) :: bernoullinumber
    ! keep track of the bernoulli numbers we have already calculated
    real(kind=prec), save, allocatable :: bernoulli(:)
    real(kind=prec), allocatable :: buffer(:)
    integer, save :: m = 2
    integer k

    if (.not. allocated(bernoulli)) then
#ifdef DEBUG
      if(verb >= 150) print*, 'initialising Bernoulli number system for m=0,', cachecontr(1)
#endif
      allocate(bernoulli(0:cachecontr(1)))
      bernoulli = 0._prec
      bernoulli(0) = 1.
      bernoulli(1) = -0.5
    endif

    if (n > size(bernoulli)) then
      ! We need more buffer
      k = size(bernoulli) - 1
      allocate(buffer(0:k))
      buffer = bernoulli
      deallocate(bernoulli)

      allocate(bernoulli(0:n + cachecontr(2)))
      bernoulli = 0._prec
      bernoulli(0:k) = buffer(0:k)

#ifdef DEBUG
      if(verb >= 150) then
        print*, 'increasing Bernoulli number system from m=0,', k, &
                ' to m=0,',size(bernoulli)-1
      endif
#endif
      deallocate(buffer)
    endif



    do m=m,n,2
      bernoulli(m) = 0.
      do k=0,m-1
        bernoulli(m) = bernoulli(m) - binom(m, k) * bernoulli(k) / (m - k + 1)
      enddo
#ifdef DEBUG
      if(verb >= 120) then
        write(*,'(A,I3,A,E25.10)') "calculating B(",m,") = ",bernoulli(m)
      endif
#endif
    enddo

    bernoullinumber = bernoulli(n)
  END FUNCTION bernoullinumber

  FUNCTION harmonicnumber(n)
    integer, intent(in) :: n
    real(kind=prec) :: harmonicnumber
    integer, parameter :: nmax = 40
    real(kind=prec), save :: Harmonic(0:nmax) = 0
    integer, save :: m = 0

    if (n > nmax) then
      print*,"Harmonic numbers bigger ",nmax," are not supported."
      stop 1
    endif
    do m=m, n+1
      harmonic(m+1) = harmonic(m) + 1._prec / real(m+1)
#ifdef DEBUG
      if(verb >= 120) then
        write(*,'(A,I3,A,E25.10)') "calculating H(",m,") = ",harmonic(m)
      endif
#endif
    enddo
    harmonicnumber = harmonic(n)

  END FUNCTION

  FUNCTION logz_polylog(n, z) result(res)
    ! Computes the classical polylogarithm Li_m(z) using series
    ! representation in log(z). valid for log z < 2pi.
    !
    ! The algorithm works by using (1.4) or [Crandall 2006]
    !
    !  PolyLog[n, z] = Ssum[Zeta[n-m] Log[z]^m / m!, {m,0,Infinity}]
    !                + Log[z]^(n-1) (HarmonicNumber[n-1] -  Log[-Log[z]]) / (n-1)!
    !
    ! where Ssum[..] excludes the singular Zeta[1] term at m = n-1. In
    ! Fortran, we split the Ssum in a sum from 0..n-2 with positive
    ! arguments in the Zeta function. The next term m=n we do manually
    ! to not have to implement Zeta[0] = -1/2 and then we use
    ! Zeta[n-m] = (-1)**(m-n) * bernoullinumber(1+m-n) / (1+m-n)
    ! for the remaining terms.
    !
    ! References
    ! R. E. Crandall, Note on fast polylogarithm computation,
    ! www.reed.edu/~crandall/papers/Polylog.pdf, January 2006.
    ! or http://functions.wolfram.com/10.08.06.0024.01
    real(kind=prec) :: fac, zetamn
    complex(kind=prec), intent(in) :: z
    complex(kind=prec) :: res, logz, del
    integer, intent(in) :: n
    integer m

    logz = log(z)

    ! The factorial will become a problem later. We do the first few
    ! terms and iterate later.
    fac = real(factorial(n-1),kind=prec)

    ! Non-sum term
    res = logz**(n-1) / fac * (harmonicnumber(n-1) - log(-logz))

    ! positive arguments in the Zeta function, 0..,n-2
    do m=0,n-2
      res = res + zeta(n-m) / factorial(m) * logz**m
    enddo

    ! Zeta[0], i.e. m=n case
    res = res - 0.5_prec * logz**n / fac / n

    ! All remaining terms
    m = n + 1
    del = 1._prec
    do while (abs(del) > LiDelta)
      zetamn = (-1)**(m-n) * bernoullinumber(1+m-n) / (1+m-n)
      fac = fac * m * (m-1)
      del = zetamn / fac * logz**m
      res = res + del
      m = m + 2
    enddo
  END FUNCTION logz_polylog


  FUNCTION Li2(x)

   !! Dilogarithm for arguments x < = 1.0

   real (kind=prec):: X,Y,T,S,A,ZERO,ONE,HALF,MALF,MONE,MTWO
   real (kind=prec):: C(0:42),H,ALFA,B0,B1,B2,LI2_OLD
   real (kind=prec):: Li2
   integer :: i, maxi

   DATA ZERO /0.0_prec/, ONE /1.0_prec/
   DATA HALF /0.5_prec/, MALF /-0.5_prec/ 
   DATA MONE /-1.0_prec/, MTWO /-2.0_prec/

   DATA C( 0) / 0.42996693560813697203703367869938799_prec/
   DATA C( 1) / 0.40975987533077105846826371092525528_prec/
   DATA C( 2) /-0.01858843665014591964764164914021227_prec/
   DATA C( 3) / 0.00145751084062267855367392841645949_prec/
   DATA C( 4) /-0.00014304184442340048774883012009088_prec/
   DATA C( 5) / 0.00001588415541879553236190550471677_prec/
   DATA C( 6) /-0.00000190784959386582722719632114209_prec/
   DATA C( 7) / 0.00000024195180854164749499461464343_prec/
   DATA C( 8) /-0.00000003193341274251783460496014143_prec/
   DATA C( 9) / 0.00000000434545062676912298795717848_prec/
   DATA C(10) /-0.00000000060578480118407444429705331_prec/
   DATA C(11) / 0.00000000008612097799359498244283685_prec/
   DATA C(12) /-0.00000000001244331659938867989642421_prec/
   DATA C(13) / 0.00000000000182255696235736330065548_prec/
   DATA C(14) /-0.00000000000027006766049114651808572_prec/
   DATA C(15) / 0.00000000000004042209263152664648833_prec/
   DATA C(16) /-0.00000000000000610325145269187950378_prec/
   DATA C(17) / 0.00000000000000092862975330195758613_prec/
   DATA C(18) /-0.00000000000000014226020855112446840_prec/
   DATA C(19) / 0.00000000000000002192631718153957354_prec/
   DATA C(20) /-0.00000000000000000339797324215897863_prec/
   DATA C(21) / 0.00000000000000000052919542448331471_prec/
   DATA C(22) /-0.00000000000000000008278580814278998_prec/
   DATA C(23) / 0.00000000000000000001300371734545560_prec/
   DATA C(24) /-0.00000000000000000000205022224255282_prec/
   DATA C(25) / 0.00000000000000000000032435785491489_prec/
   DATA C(26) /-0.00000000000000000000005147799903343_prec/
   DATA C(27) / 0.00000000000000000000000819387747717_prec/
   DATA C(28) /-0.00000000000000000000000130778354057_prec/
   DATA C(29) / 0.00000000000000000000000020925629306_prec/
   DATA C(30) /-0.00000000000000000000000003356166151_prec/
   DATA C(31) / 0.00000000000000000000000000539465777_prec/
   DATA C(32) /-0.00000000000000000000000000086891932_prec/
   DATA C(33) / 0.00000000000000000000000000014022817_prec/
   DATA C(34) /-0.00000000000000000000000000002267156_prec/
   DATA C(35) / 0.00000000000000000000000000000367174_prec/
   DATA C(36) /-0.00000000000000000000000000000059562_prec/
   DATA C(37) / 0.00000000000000000000000000000009677_prec/
   DATA C(38) /-0.00000000000000000000000000000001574_prec/
   DATA C(39) / 0.00000000000000000000000000000000257_prec/
   DATA C(40) /-0.00000000000000000000000000000000042_prec/
   DATA C(41) / 0.00000000000000000000000000000000007_prec/
   DATA C(42) /-0.00000000000000000000000000000000001_prec/

   if(X > 1.00000000001_prec) then
     print*, 'crashes because Li called with bad arguments'
   elseif(X > 1.0_prec) then
     X = 1._prec
   endif    

   IF(X > 0.999999_prec) THEN
    LI2_OLD=zeta(2)
    Li2 = Real(LI2_OLD,prec)
    RETURN
   ELSE IF(abs(x-MONE) < zero) THEN
    LI2_OLD=MALF*zeta(2)
    RETURN
   END IF
   T=-X
   IF(T .LE. MTWO) THEN
    Y=MONE/(ONE+T)
    S=ONE
    A=-2*zeta(2)+HALF*(LOG(-T)**2-LOG(ONE+ONE/T)**2)
   ELSE IF(T .LT. MONE) THEN
    Y=MONE-T
    S=MONE
    A=LOG(-T)
    A=-zeta(2)+A*(A+LOG(ONE+ONE/T))
   ELSE IF(T .LE. MALF) THEN
    Y=(MONE-T)/T
    S=ONE
    A=LOG(-T)
    A=-zeta(2)+A*(MALF*A+LOG(ONE+T))
   ELSE IF(T .LT. ZERO) THEN
    Y=-T/(ONE+T)
    S=MONE
    A=HALF*LOG(ONE+T)**2
   ELSE IF(T .LE. ONE) THEN
    Y=T
    S=ONE
    A=ZERO
   ELSE
    Y=ONE/T
    S=MONE
    A=zeta(2)+HALF*LOG(T)**2
   END IF

   H=Y+Y-ONE
   ALFA=H+H
   B1=ZERO
   B2=ZERO
   if (precision(1._prec) < 20) then
     maxi = 18
   else
     maxi = 42
   endif
   DO  I = maxi,0,-1
     B0=C(I)+ALFA*B1-B2
     B2=B1
     B1=B0
   ENDDO
   LI2_OLD=-(S*(B0-H*B2)+A)
         ! Artificial conversion           
   Li2 = Real(LI2_OLD,prec)
  END FUNCTION Li2

  RECURSIVE FUNCTION dilog(x) result(res)
    ! evaluates dilog for any argument |x|<1
    complex(kind=prec) :: res
    complex(kind=prec) :: x

    if(abs(aimag(x)) < zero ) then
      res = Li2(real(x,kind=prec))
    else if ( (0.5_prec .lt. abs(x)) .and. (abs(x) .lt. 2._prec) ) then
      res = logz_polylog(2,x)
    else
      res = naive_polylog(2,x)
    endif
  END FUNCTION dilog

  FUNCTION Li3(x)
    ! Trilogarithm for arguments x < = 1.0
    ! This was hacked from LI2 to also follow C332
    ! In theory this could also produce Re[Li [x]] for x>1

    real (kind=prec):: X,S,A
    real (kind=prec):: CA(0:52),HA,ALFAA,BA0,BA1,BA2, YA
    real (kind=prec):: CB(0:52),HB,ALFAB,BB0,BB1,BB2, YB

    DATA CA( 0) / 0.46172939286012092817954516381760016_prec/
    DATA CA( 1) / 0.45017399588550288560580364647352070_prec/
    DATA CA( 2) /-0.01091284195229295374494914320402658_prec/
    DATA CA( 3) / 0.00059324547127243642952756961712713_prec/
    DATA CA( 4) /-0.00004479593219280756178998757870776_prec/
    DATA CA( 5) / 0.00000405154578580684540984293800468_prec/
    DATA CA( 6) /-0.00000041095398606214457668547736075_prec/
    DATA CA( 7) / 0.00000004513178770934181970313262557_prec/
    DATA CA( 8) /-0.00000000525466158515342604029419927_prec/
    DATA CA( 9) / 0.00000000063982547910449452549291936_prec/
    DATA CA(10) /-0.00000000008071938839872532971820424_prec/
    DATA CA(11) / 0.00000000001048073680087126094928657_prec/
    DATA CA(12) /-0.00000000000139365085138335067524094_prec/
    DATA CA(13) / 0.00000000000018907205037339730044704_prec/
    DATA CA(14) /-0.00000000000002609371657183250621931_prec/
    DATA CA(15) / 0.00000000000000365481859219879483309_prec/
    DATA CA(16) /-0.00000000000000051855842492271228151_prec/
    DATA CA(17) / 0.00000000000000007441491722173878908_prec/
    DATA CA(18) /-0.00000000000000001078686838424874221_prec/
    DATA CA(19) / 0.00000000000000000157774237809543778_prec/
    DATA CA(20) /-0.00000000000000000023264073800573828_prec/
    DATA CA(21) / 0.00000000000000000003455457587964154_prec/
    DATA CA(22) /-0.00000000000000000000516658458392580_prec/
    DATA CA(23) / 0.00000000000000000000077718849383139_prec/
    DATA CA(24) /-0.00000000000000000000011755815708807_prec/
    DATA CA(25) / 0.00000000000000000000001787262690583_prec/
    DATA CA(26) /-0.00000000000000000000000272999302683_prec/
    DATA CA(27) / 0.00000000000000000000000041881267359_prec/
    DATA CA(28) /-0.00000000000000000000000006451004176_prec/
    DATA CA(29) / 0.00000000000000000000000000997383916_prec/
    DATA CA(30) /-0.00000000000000000000000000154744603_prec/
    DATA CA(31) / 0.00000000000000000000000000024087296_prec/
    DATA CA(32) /-0.00000000000000000000000000003760889_prec/
    DATA CA(33) / 0.00000000000000000000000000000588900_prec/
    DATA CA(34) /-0.00000000000000000000000000000092463_prec/
    DATA CA(35) / 0.00000000000000000000000000000014555_prec/
    DATA CA(36) /-0.00000000000000000000000000000002297_prec/
    DATA CA(37) / 0.00000000000000000000000000000000363_prec/
    DATA CA(38) /-0.00000000000000000000000000000000058_prec/
    DATA CA(39) / 0.00000000000000000000000000000000009_prec/
    DATA CA(40) /-0.00000000000000000000000000000000001_prec/
    DATA CB( 0) /-0.01601618044919582873670691984756338_prec/
    DATA CB( 1) /-0.50364244007530129181209541016960792_prec/
    DATA CB( 2) /-0.01615099243050025888745446951929454_prec/
    DATA CB( 3) /-0.00124402421042449361265610524413112_prec/
    DATA CB( 4) /-0.00013757218124461673921971996409271_prec/
    DATA CB( 5) /-0.00001856381852603773316486795183129_prec/
    DATA CB( 6) /-0.00000284173534515440415934505790039_prec/
    DATA CB( 7) /-0.00000047459967905789937221638951390_prec/
    DATA CB( 8) /-0.00000008448038543781200676091819474_prec/
    DATA CB( 9) /-0.00000001578767124043400543246870475_prec/
    DATA CB(10) /-0.00000000306576207139903798128889004_prec/
    DATA CB(11) /-0.00000000061407921728125845808062189_prec/
    DATA CB(12) /-0.00000000012618830243156719690872484_prec/
    DATA CB(13) /-0.00000000002649314179819609957126783_prec/
    DATA CB(14) /-0.00000000000566470854636425926158812_prec/
    DATA CB(15) /-0.00000000000123041115779581117517467_prec/
    DATA CB(16) /-0.00000000000027093457836786768143960_prec/
    DATA CB(17) /-0.00000000000006038026463383701279197_prec/
    DATA CB(18) /-0.00000000000001360008993995749682352_prec/
    DATA CB(19) /-0.00000000000000309244740631856875855_prec/
    DATA CB(20) /-0.00000000000000070917249609207158220_prec/
    DATA CB(21) /-0.00000000000000016388083639226002471_prec/
    DATA CB(22) /-0.00000000000000003813464350168994613_prec/
    DATA CB(23) /-0.00000000000000000893010739611811656_prec/
    DATA CB(24) /-0.00000000000000000210331341599359416_prec/
    DATA CB(25) /-0.00000000000000000049802988416537866_prec/
    DATA CB(26) /-0.00000000000000000011850292695597351_prec/
    DATA CB(27) /-0.00000000000000000002832460494402074_prec/
    DATA CB(28) /-0.00000000000000000000679854955943073_prec/
    DATA CB(29) /-0.00000000000000000000163816629435900_prec/
    DATA CB(30) /-0.00000000000000000000039616291258646_prec/
    DATA CB(31) /-0.00000000000000000000009613022139972_prec/
    DATA CB(32) /-0.00000000000000000000002340035706102_prec/
    DATA CB(33) /-0.00000000000000000000000571315840877_prec/
    DATA CB(34) /-0.00000000000000000000000139876183805_prec/
    DATA CB(35) /-0.00000000000000000000000034336361321_prec/
    DATA CB(36) /-0.00000000000000000000000008449733573_prec/
    DATA CB(37) /-0.00000000000000000000000002084253881_prec/
    DATA CB(38) /-0.00000000000000000000000000515255292_prec/
    DATA CB(39) /-0.00000000000000000000000000127646290_prec/
    DATA CB(40) /-0.00000000000000000000000000031685555_prec/
    DATA CB(41) /-0.00000000000000000000000000007880228_prec/
    DATA CB(42) /-0.00000000000000000000000000001963363_prec/
    DATA CB(43) /-0.00000000000000000000000000000490016_prec/
    DATA CB(44) /-0.00000000000000000000000000000122499_prec/
    DATA CB(45) /-0.00000000000000000000000000000030671_prec/
    DATA CB(46) /-0.00000000000000000000000000000007691_prec/
    DATA CB(47) /-0.00000000000000000000000000000001931_prec/
    DATA CB(48) /-0.00000000000000000000000000000000486_prec/
    DATA CB(49) /-0.00000000000000000000000000000000122_prec/
    DATA CB(50) /-0.00000000000000000000000000000000031_prec/
    DATA CB(51) /-0.00000000000000000000000000000000008_prec/
    DATA CB(52) /-0.00000000000000000000000000000000002_prec/
    real (kind=prec):: Li3
    integer :: i, maxi


    if(x > 1.00000000001_prec) then
      print*, 'need to crash Li3, since not convergent'
    elseif(x > 1.0_prec) then
      x = 1._prec
    endif

    IF(X > 0.999999_prec) THEN
      LI3=zeta(3)
    RETURN
    ELSE IF( abs(x+1) < zero) THEN
      LI3=-0.75_prec*zeta(3)
    RETURN
    END IF
    IF(X .LE. -1._prec) THEN
      YA=1._prec/x ; YB=0._prec
      S=-1._prec
      A=-LOG(-X)*(zeta(2)+LOG(-x)**2/6._prec)
    ELSE IF(X .LE. 0._prec) THEN
      YA=x ; YB=0._prec
      S=-1._prec
      A=0._prec
    ELSE IF(X .LE. 0.5_prec) THEN
      YA=0._prec ; YB=x
      S=-1._prec
      A=0._prec
    ELSE IF(X .LE. 1._prec) THEN
      YA=(x-1._prec)/x ; YB=1._prec-x
      S=1._prec
      A=zeta(3) + zeta(2)*Log(x) - (Log(1._prec - X)*Log(X)**2)/2._prec + Log(X)**3/6._prec
    ELSE IF(X .LE. 2._prec) THEN
      YA=1._prec - X ; YB=(X-1._prec)/X
      S=1._prec
      A=zeta(3) + zeta(2)*Log(x) - (Log(X - 1._prec)*Log(X)**2)/2._prec + Log(X)**3/6._prec
    ELSE
      YA=0._prec ; YB=1._prec/X
      S=-1._prec
      A=2*zeta(2)*Log(x)-Log(x)**3/6._prec
    END IF


    HA=-2._prec*YA-1._prec ; HB= 2._prec*YB
    ALFAA=HA+HA ; ALFAB = HB+HB

    BA0 = 0. ; BA1=0. ; BA2=0.
    BB0 = 0. ; BB1=0. ; BB2=0.
    if (precision(1._prec) < 20) then
      maxi = 18
    else
      maxi = 42
    endif
    DO  I = maxi,0,-1
       BA0=CA(I)+ALFAA*BA1-BA2 ; BA2=BA1 ; BA1=BA0
       BB0=CB(I)+ALFAB*BB1-BB2 ; BB2=BB1 ; BB1=BB0
    ENDDO
    Li3 = A + S * (  (BA0 - HA*BA2) + (BB0 - HB*BB2) )
  END FUNCTION Li3

  FUNCTION trilog(x) result(res)
    ! evaluates trilog for any argument |x|<1
    complex(kind=prec) :: res
    complex(kind=prec) :: x
    if(abs(aimag(x)) < zero ) then
      res = Li3(real(x,kind=prec))
    else if ( (0.5_prec .lt. abs(x)) .and. (abs(x) .lt. 2._prec) ) then
      res = logz_polylog(3,x)
    else
      res = naive_polylog(3,x)
    endif
  END FUNCTION trilog

  FUNCTION BERNOULLI_POLYNOMIAL(n, x) result(res)
    integer, parameter :: maxn = 15
    integer n
    complex(kind=prec) :: x, res
    complex(kind=prec) :: xpow(maxn+1)
    integer, parameter :: coeffN(maxn+1, maxn) = reshape((/ &
        -   1, +   1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0, &
        +   1, -   1, +   1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0, &
            0, +   1, -   3, +   1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0, &
        -   1,     0, +   1, -   2, +   1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0, &
            0, -   1,     0, +   5, -   5, +   1,     0,     0,     0,     0,     0,     0,     0,     0,     0,     0, &
        +   1,     0, -   1,     0, +   5, -   3, +   1,     0,     0,     0,     0,     0,     0,     0,     0,     0, &
            0, +   1,     0, -   7,     0, +   7, -   7, +   1,     0,     0,     0,     0,     0,     0,     0,     0, &
        -   1,     0, +   2,     0, -   7,     0, +  14, -   4, +   1,     0,     0,     0,     0,     0,     0,     0, &
            0, -   3,     0, +   2,     0, -  21,     0, +   6, -   9, +   1,     0,     0,     0,     0,     0,     0, &
        +   5,     0, -   3,     0, +   5,     0, -   7,     0, +  15, -   5, +   1,     0,     0,     0,     0,     0, &
            0, +   5,     0, -  11,     0, +  11,     0, -  11,     0, +  55, -  11, +   1,     0,     0,     0,     0, &
        - 691,     0, +   5,     0, -  33,     0, +  22,     0, -  33,     0, +  11, -   6, +   1,     0,     0,     0, &
            0, - 691,     0, +  65,     0, - 429,     0, + 286,     0, - 143,     0, +  13, -  13, +   1,     0,     0, &
        +   7,     0, - 691,     0, + 455,     0, -1001,     0, + 143,     0, -1001,     0, +  91, -   7, +   1,     0, &
            0, +  35,     0, - 691,     0, + 455,     0, - 429,     0, + 715,     0, -  91,     0, +  35, -  15, +   1 /), &
            (/maxn+1, maxn/))
    integer, parameter :: coeffD(maxn+1, maxn) = reshape((/ &
        +   2, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +   6, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +   1, +   2, +   2, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +  30, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +   1, +   6, +   1, +   3, +   2, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +  42, +   1, +   2, +   1, +   2, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +   1, +   6, +   1, +   6, +   1, +   2, +   2, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +  30, +   1, +   3, +   1, +   3, +   1, +   3, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +   1, +  10, +   1, +   1, +   1, +   5, +   1, +   1, +   2, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +  66, +   1, +   2, +   1, +   1, +   1, +   1, +   1, +   2, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +   1, +   6, +   1, +   2, +   1, +   1, +   1, +   1, +   1, +   6, +   2, +   1, +   1, +   1, +   1, +   1, &
        +2730, +   1, +   1, +   1, +   2, +   1, +   1, +   1, +   2, +   1, +   1, +   1, +   1, +   1, +   1, +   1, &
        +   1, + 210, +   1, +   3, +   1, +  10, +   1, +   7, +   1, +   6, +   1, +   1, +   2, +   1, +   1, +   1, &
        +   6, +   1, +  30, +   1, +   6, +   1, +  10, +   1, +   2, +   1, +  30, +   1, +   6, +   1, +   1, +   1, &
        +   1, +   2, +   1, +   6, +   1, +   2, +   1, +   2, +   1, +   6, +   1, +   2, +   1, +   2, +   2, +   1 /), &
            (/maxn+1, maxn/))

    real(kind=prec), parameter :: coeff(maxn+1,maxn) = coeffN/real(coeffD,kind=prec)
    integer i

    if (n>maxn) then
      print*,"Bernoulli beyond 15 is not implemented"
      stop
    endif

    xpow(1:n+1) = (/ ( x**i, i = 0, n ) /)
    res = sum( xpow(1:n+1) * coeff(1:n+1,n) )

  END FUNCTION

  FUNCTION MYLOG(x)
    complex(kind=prec) :: x, mylog
    if (abs(aimag(x)) < zero) then
      if (real(x) > 0) then
        mylog = cmplx(log(real(+x,kind=prec)),     kind=prec)
      else
        mylog = cmplx(log(real(-x,kind=prec)), pi, kind=prec)
      endif
    else
      mylog = log(x)
    endif
  END FUNCTION

  RECURSIVE FUNCTION polylog1(m,x) result(res)
    ! computes the polylog
    
    integer :: m
    complex(kind=prec) :: x, inv
    complex(kind=prec) :: res
    integer i

    
#ifdef DEBUG
    if(verb >= 70) print*, 'called polylog(',m,',',x,')'
#endif
#ifndef NOCACHE
    if (m.le.5) then
      do i=1,plcachesize(m)
        if( abs(cache(m,i)%c-x).lt.zero ) then
          res = cache(m,i)%ans
          return
        endif
      enddo
    endif
#endif
    if ((m.le.9).and.(abs(x-1.).lt.zero)) then
      res = zeta(m)
    else if ((m.le.9).and.(abs(x+1._prec).lt.zero)) then
      res = -(1._prec - 2._prec**(1-m))*zeta(m)
    else if ((m.le.9).and.(abs(x-i_).lt.zero)) then
      res = -(0.5_prec**m - 0.5_prec**(2*m-1)) * zeta(m) + i_*dirichletbeta(m)
    else if ((m.le.9).and.(abs(x+i_).lt.zero)) then
      res = -(0.5_prec**m - 0.5_prec**(2*m-1)) * zeta(m) - i_*dirichletbeta(m)
    else if (abs(x) .gt. 1) then
      inv = 1._prec/x
      res = (-1)**(m-1)*polylog(m,inv) &
          - (2._prec*pi*i_)**m * bernoulli_polynomial(m, 0.5_prec-i_*mylog(-x)/2._prec/pi) / factorial(m)
    else if(m == 2) then
      res = dilog(x)
    else if(m == 3) then
      res = trilog(x)
    else if ( (0.5_prec .lt. abs(x)) .and. (abs(x) .lt. 2._prec) ) then
      res = logz_polylog(m,x)
    else
      res = naive_polylog(m,x)
    end if

#ifndef NOCACHE
    if (m.le.PolyLogCacheSize(1)) then
      if (plcachesize(m).lt.PolyLogCacheSize(2)) then
        plcachesize(m) = plcachesize(m) + 1
        cache(m,plcachesize(m)) = el(x,res)
      endif
    endif
#endif
  END FUNCTION polylog1




  RECURSIVE FUNCTION polylog2(m,x,y) result(res)
    type(inum) :: x, y
    integer m
    complex(kind=prec) :: res
    res=polylog1(m,x%c/y%c)
    if ( (abs(aimag(x)).lt.zero).and.(abs(aimag(y)).lt.zero) ) then
      ! Both arguments are real, only here does the ieps matter
      ! FIXME this is rather ugly..
      if (real(x).gt.real(y) .and. real(y).gt. 0) then
        ! Force the sign to be -b%i0
        res = cmplx(real(res), sign(aimag(res), -y%i0*1._prec), kind=prec)
      elseif(real(x).lt.real(y) .and. real(y).lt. 0) then
        res = cmplx(real(res), sign(aimag(res), +y%i0*1._prec), kind=prec)
      endif
    endif
  END FUNCTION POLYLOG2


  FUNCTION PLOG1(a,b)
  ! calculates log(1-a/b)
  implicit none
  type(inum) :: a,b
  complex(kind=prec) plog1

  if ( (abs(aimag(a)).lt.zero).and.(abs(aimag(b)).lt.zero) ) then
    ! Both arguments are real, only here does the ieps matter
    plog1 = log(abs(1.-a%c/b%c))
    ! this does not depend on the sign of a
    if (real(a).gt.real(b) .and. real(b).gt. 0) then
      plog1 = plog1 + b%i0*i_*pi
    elseif(real(a).lt.real(b) .and. real(b).lt. 0) then
      plog1 = plog1 - b%i0*i_*pi
    endif
  else
    plog1 = mylog(1.-a%c/b%c)
  endif
  END FUNCTION

#ifndef NOCACHE
  SUBROUTINE CLEARCACHE
  plcachesize=0
  END SUBROUTINE
#endif

END MODULE maths_functions

! PROGRAM test
!   use maths_functions
!   implicit none
!   complex(kind=prec) :: res
!   res = Li3(0.4d0)
!   print*, res
! END PROGRAM test

