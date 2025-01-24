      module mod_vvamp 
          use mod_vvamp_eval1
          use mod_vvamp_eval2
          use mod_vvamp_eval3
          use mod_vvamp_eval4
          use mod_vvamp_eval5
          use mod_vvamp_eval6
          use mod_vvamp_eval7
          use mod_vvamp_eval8
          use mod_vvamp_eval9
          use mod_vvamp_eval10
          use mod_vvamp_eval11
          use mod_vvamp_eval12

          public :: qqbAjinit
          public :: qqbAjfill

          private

      contains

      subroutine qqbAjinit(s,t,u,p3sq,p4sq,xqqbr)
      implicit none
      include 'types.f'
      include 'constants.f'
      real(dp):: s,t,u,p3sq,p4sq
      real(dp):: ma2,mb2,smin,t0,t1,bb
      real(dp):: x,y,z,xqqbr(7)

      ma2=p3sq
      mb2=p4sq

! check whether phase space point in physical region

      if ((ma2.le.zip) .or. (mb2.le.zip)) then
         write(6,*) 'Error: illegal ma2 or mb2 value'
         write(6,*) 'ma2: ',ma2
         write(6,*) 'mb2: ',mb2
         stop
      endif

      smin=(sqrt(ma2) + sqrt(mb2))**2
      if (s.lt.smin) then
         write(6,*) 'Error: illegal s value'
         write(6,*) 'smin ',smin
         write(6,*) 's: ',s
         stop
      endif

!      root=dsqrt((s-ma2-mb2)**2-4._dp*ma2*mb2)
!      t0=(ma2 + mb2 - s - root)/2._dp
!      t1=(ma2 + mb2 - s + root)/2._dp
! accurate determination of the roots of a quadratic
      bb=s-ma2-mb2
      t0=-half*(bb+sign(one,bb)*sqrt(bb**2-4._dp*ma2*mb2))
      t1=ma2*mb2/t0
      if ((t.lt.t0) .or. (t.gt.t1)) then
         write(6,*) 'Error: illegal t value'
         write(6,*) 't0: ',t0
         write(6,*) 't1: ',t1
         write(6,*) 't: ',t
         stop
      endif

! prepare x,y,z
      x =  -t0/ma2
      y =  t1*(u+t)/(ma2*mb2)-1._dp
      z =  t1*t/(ma2*mb2)

!
      xqqbr(1) = ma2
      xqqbr(2) = t/ma2
      xqqbr(3) = u/ma2
      xqqbr(4) = mb2/ma2
      xqqbr(5) = x
      xqqbr(6) = y
      xqqbr(7) = z

      return
      end

      subroutine qqbAjfill(s12,t,u,p3sq,p4sq,order,nflav,qqbAj)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'rcpp.f'
      real(dp):: s12,t,u,p3sq,p4sq,xqqbr(7)
      complex(dp)::qqbAj(0:2,4,10)
      integer:: order,nflav
      integer,parameter:: N=3
      integer,parameter::A=1,B=2,C=3,F=4
      real(dp)::qqbAi0A1(0:1),qqbAi0A2(0:1),qqbAi0A3(0:1),qqbAi0A4(0:1),qqbAi0A5(0:1),
     &          qqbAi0A6(0:1),qqbAi0A7(0:1),qqbAi0A8(0:1),qqbAi0A9(0:1),qqbAi0A10(0:1)
      real(dp)::qqbAi0B1(0:1),qqbAi0B2(0:1),qqbAi0B3(0:1),qqbAi0B4(0:1),qqbAi0B5(0:1),
     &          qqbAi0B6(0:1),qqbAi0B7(0:1),qqbAi0B8(0:1),qqbAi0B9(0:1),qqbAi0B10(0:1)
      real(dp)::qqbAi0C1(0:1),qqbAi0C2(0:1),qqbAi0C3(0:1),qqbAi0C4(0:1),qqbAi0C5(0:1),
     &          qqbAi0C6(0:1),qqbAi0C7(0:1),qqbAi0C8(0:1),qqbAi0C9(0:1),qqbAi0C10(0:1)
      real(dp)::qqbAi1A1(0:1),qqbAi1A2(0:1),qqbAi1A3(0:1),qqbAi1A4(0:1),qqbAi1A5(0:1),
     &          qqbAi1A6(0:1),qqbAi1A7(0:1),qqbAi1A8(0:1),qqbAi1A9(0:1),qqbAi1A10(0:1)
      real(dp)::qqbAi1B1(0:1),qqbAi1B2(0:1),qqbAi1B3(0:1),qqbAi1B4(0:1),qqbAi1B5(0:1),
     &          qqbAi1B6(0:1),qqbAi1B7(0:1),qqbAi1B8(0:1),qqbAi1B9(0:1),qqbAi1B10(0:1)
      real(dp)::qqbAi1C1(0:1),qqbAi1C2(0:1),qqbAi1C3(0:1),qqbAi1C4(0:1),qqbAi1C5(0:1),
     &          qqbAi1C6(0:1),qqbAi1C7(0:1),qqbAi1C8(0:1),qqbAi1C9(0:1),qqbAi1C10(0:1)
      real(dp)::qqbAi2A1(0:1),qqbAi2A2(0:1),qqbAi2A3(0:1),qqbAi2A4(0:1),qqbAi2A5(0:1),
     &          qqbAi2A6(0:1),qqbAi2A7(0:1),qqbAi2A8(0:1),qqbAi2A9(0:1),qqbAi2A10(0:1)
      real(dp)::qqbAi2B1(0:1),qqbAi2B2(0:1),qqbAi2B3(0:1),qqbAi2B4(0:1),qqbAi2B5(0:1),
     &          qqbAi2B6(0:1),qqbAi2B7(0:1),qqbAi2B8(0:1),qqbAi2B9(0:1),qqbAi2B10(0:1)
      real(dp)::qqbAi2C1(0:1),qqbAi2C2(0:1),qqbAi2C3(0:1),qqbAi2C4(0:1),qqbAi2C5(0:1),
     &          qqbAi2C6(0:1),qqbAi2C7(0:1),qqbAi2C8(0:1),qqbAi2C9(0:1),qqbAi2C10(0:1)
!      complex(dp)::qqbEj2A1,qqbEj2A2,qqbEj2A3,qqbEj2A4,qqbEj2A5,qqbEj2A6,qqbEj2A7,qqbEj2A8,qqbEj2A9

      allocate(qqbr(1:611275), source=0._dp)

      call qqbAjinit(s12,t,u,p3sq,p4sq,xqqbr)
      qqbr(1:7)=xqqbr(1:7)

      if (order.eq.0) then
          call qqb0276(qqbr); call qqb3050(qqbr);call qqb3051(qqbr)
      elseif (order.eq.1) then
          call qqb0002(qqbr);call qqb0003(qqbr);call qqb0005(qqbr);call qqb0006(qqbr);call qqb0010(qqbr);call qqb0011(qqbr);call
     & qqb0012(qqbr);call qqb0015(qqbr);call qqb0276(qqbr);call qqb0277(qqbr)
          call qqb0278(qqbr);call qqb0279(qqbr);call qqb0280(qqbr);call qqb0281(qqbr);call qqb0282(qqbr);call qqb0283(qqbr);call
     & qqb0284(qqbr);call qqb0285(qqbr);call qqb0286(qqbr);call qqb0287(qqbr)
          call qqb0288(qqbr);call qqb0291(qqbr);call qqb0292(qqbr);call qqb0293(qqbr);call qqb0294(qqbr);call qqb0295(qqbr);call
     & qqb0297(qqbr);call qqb0303(qqbr);call qqb0304(qqbr);call qqb0305(qqbr)
          call qqb0310(qqbr);call qqb0311(qqbr);call qqb0312(qqbr);call qqb0313(qqbr);call qqb0318(qqbr);call qqb0319(qqbr);call
     & qqb0320(qqbr);call qqb0321(qqbr);call qqb0326(qqbr);call qqb0327(qqbr)
          call qqb0333(qqbr);call qqb0334(qqbr);call qqb0335(qqbr);call qqb0340(qqbr);call qqb0341(qqbr);call qqb0351(qqbr);call
     & qqb0352(qqbr);call qqb0353(qqbr);call qqb0354(qqbr);call qqb0355(qqbr)
          call qqb0357(qqbr);call qqb0358(qqbr);call qqb0361(qqbr);call qqb0366(qqbr);call qqb0367(qqbr);call qqb0368(qqbr);call
     & qqb0369(qqbr);call qqb0370(qqbr);call qqb0373(qqbr);call qqb0378(qqbr)
          call qqb0379(qqbr);call qqb0380(qqbr);call qqb0382(qqbr);call qqb0383(qqbr);call qqb0410(qqbr);call qqb0411(qqbr);call
     & qqb0414(qqbr);call qqb0425(qqbr);call qqb0426(qqbr);call qqb0428(qqbr)
          call qqb0444(qqbr);call qqb0445(qqbr);call qqb0446(qqbr);call qqb0447(qqbr);call qqb0448(qqbr);call qqb0461(qqbr);call
     & qqb0463(qqbr);call qqb0468(qqbr);call qqb0482(qqbr);call qqb0483(qqbr)
          call qqb0484(qqbr);call qqb0485(qqbr);call qqb0488(qqbr);call qqb0499(qqbr);call qqb0500(qqbr);call qqb0501(qqbr);call
     & qqb0502(qqbr);call qqb0505(qqbr);call qqb0518(qqbr);call qqb0521(qqbr)
          call qqb0523(qqbr);call qqb0539(qqbr);call qqb0540(qqbr);call qqb0542(qqbr);call qqb0545(qqbr);call qqb0564(qqbr);call
     & qqb0565(qqbr);call qqb0568(qqbr);call qqb0571(qqbr);call qqb0591(qqbr)
          call qqb0594(qqbr);call qqb0614(qqbr);call qqb0615(qqbr);call qqb0617(qqbr);call qqb0621(qqbr);call qqb0649(qqbr);call
     & qqb0664(qqbr);call qqb0666(qqbr);call qqb0692(qqbr);call qqb0712(qqbr)
          call qqb0715(qqbr);call qqb0737(qqbr);call qqb0755(qqbr);call qqb0904(qqbr);call qqb0906(qqbr);call qqb0927(qqbr);call
     & qqb0929(qqbr);call qqb0930(qqbr);call qqb0949(qqbr);call qqb0950(qqbr)
          call qqb0970(qqbr);call qqb1068(qqbr);call qqb1069(qqbr);call qqb1091(qqbr);call qqb1092(qqbr);call qqb1093(qqbr);call
     & qqb1095(qqbr);call qqb1097(qqbr);call qqb1098(qqbr);call qqb1118(qqbr)
          call qqb1121(qqbr);call qqb1122(qqbr);call qqb1125(qqbr);call qqb1145(qqbr);call qqb1148(qqbr);call qqb1196(qqbr);call
     & qqb1217(qqbr);call qqb1301(qqbr);call qqb1303(qqbr);call qqb1306(qqbr)
          call qqb1327(qqbr);call qqb1328(qqbr);call qqb1330(qqbr);call qqb1331(qqbr);call qqb1354(qqbr);call qqb1357(qqbr);call
     & qqb1362(qqbr);call qqb1404(qqbr);call qqb1405(qqbr);call qqb1513(qqbr)
          call qqb1515(qqbr);call qqb1516(qqbr);call qqb1521(qqbr);call qqb1543(qqbr);call qqb1548(qqbr);call qqb1552(qqbr);call
     & qqb1577(qqbr);call qqb1609(qqbr);call qqb1610(qqbr);call qqb1803(qqbr)
          call qqb1804(qqbr);call qqb2539(qqbr);call qqb2540(qqbr);call qqb2541(qqbr);call qqb2549(qqbr);call qqb2550(qqbr);call
     & qqb2555(qqbr);call qqb2556(qqbr);call qqb2557(qqbr);call qqb2562(qqbr)
          call qqb2563(qqbr);call qqb2564(qqbr);call qqb2565(qqbr);call qqb2569(qqbr);call qqb2575(qqbr);call qqb2576(qqbr);call
     & qqb2584(qqbr);call qqb2590(qqbr);call qqb2596(qqbr);call qqb2600(qqbr)
          call qqb2605(qqbr);call qqb2606(qqbr);call qqb2607(qqbr);call qqb2608(qqbr);call qqb2843(qqbr);call qqb2846(qqbr);call
     & qqb2848(qqbr);call qqb2912(qqbr);call qqb2913(qqbr);call qqb2914(qqbr)
          call qqb2918(qqbr);call qqb2948(qqbr);call qqb2949(qqbr);call qqb2950(qqbr);call qqb2951(qqbr);call qqb2952(qqbr);call
     & qqb2953(qqbr);call qqb2954(qqbr);call qqb2957(qqbr);call qqb2960(qqbr)
          call qqb2962(qqbr);call qqb2965(qqbr);call qqb2966(qqbr);call qqb2968(qqbr);call qqb2969(qqbr);call qqb2973(qqbr);call
     & qqb3016(qqbr);call qqb3017(qqbr);call qqb3018(qqbr);call qqb3019(qqbr)
          call qqb3020(qqbr);call qqb3021(qqbr);call qqb3022(qqbr);call qqb3023(qqbr);call qqb3024(qqbr);call qqb3025(qqbr);call
     & qqb3038(qqbr);call qqb3039(qqbr);call qqb3050(qqbr);call qqb3051(qqbr)
      elseif (order.eq.2) then
                         call qqb0002(qqbr);call qqb0003(qqbr);call qqb0004(qqbr);call qqb0005(qqbr);call qqb0006(qqbr);call
     & qqb0007(qqbr);call qqb0008(qqbr);call qqb0009(qqbr);call qqb0010(qqbr)
          call qqb0011(qqbr);call qqb0012(qqbr);call qqb0013(qqbr);call qqb0014(qqbr);call qqb0015(qqbr);call qqb0016(qqbr);call
     & qqb0017(qqbr);call qqb0018(qqbr);call qqb0019(qqbr);call qqb0020(qqbr)
          call qqb0021(qqbr);call qqb0022(qqbr);call qqb0023(qqbr);call qqb0024(qqbr);call qqb0025(qqbr);call qqb0026(qqbr);call
     & qqb0027(qqbr);call qqb0028(qqbr);call qqb0029(qqbr);call qqb0030(qqbr)
          call qqb0031(qqbr);call qqb0032(qqbr);call qqb0033(qqbr);call qqb0034(qqbr);call qqb0035(qqbr);call qqb0036(qqbr);call
     & qqb0037(qqbr);call qqb0038(qqbr);call qqb0039(qqbr);call qqb0040(qqbr)
          call qqb0041(qqbr);call qqb0042(qqbr);call qqb0043(qqbr);call qqb0044(qqbr);call qqb0045(qqbr);call qqb0046(qqbr);call
     & qqb0047(qqbr);call qqb0048(qqbr);call qqb0049(qqbr);call qqb0050(qqbr)
          call qqb0051(qqbr);call qqb0052(qqbr);call qqb0053(qqbr);call qqb0054(qqbr);call qqb0055(qqbr);call qqb0056(qqbr);call
     & qqb0057(qqbr);call qqb0058(qqbr);call qqb0059(qqbr);call qqb0060(qqbr)
          call qqb0061(qqbr);call qqb0062(qqbr);call qqb0063(qqbr);call qqb0064(qqbr);call qqb0065(qqbr);call qqb0066(qqbr);call
     & qqb0067(qqbr);call qqb0068(qqbr);call qqb0069(qqbr);call qqb0070(qqbr)
          call qqb0071(qqbr);call qqb0072(qqbr);call qqb0073(qqbr);call qqb0074(qqbr);call qqb0075(qqbr);call qqb0076(qqbr);call
     & qqb0077(qqbr);call qqb0078(qqbr);call qqb0079(qqbr);call qqb0080(qqbr)
          call qqb0081(qqbr);call qqb0082(qqbr);call qqb0083(qqbr);call qqb0084(qqbr);call qqb0085(qqbr);call qqb0086(qqbr);call
     & qqb0087(qqbr);call qqb0088(qqbr);call qqb0089(qqbr);call qqb0090(qqbr)
          call qqb0091(qqbr);call qqb0092(qqbr);call qqb0093(qqbr);call qqb0094(qqbr);call qqb0095(qqbr);call qqb0096(qqbr);call
     & qqb0097(qqbr);call qqb0098(qqbr);call qqb0099(qqbr);call qqb0100(qqbr)
          call qqb0101(qqbr);call qqb0102(qqbr);call qqb0103(qqbr);call qqb0104(qqbr);call qqb0105(qqbr);call qqb0106(qqbr);call
     & qqb0107(qqbr);call qqb0108(qqbr);call qqb0109(qqbr);call qqb0110(qqbr)
          call qqb0111(qqbr);call qqb0112(qqbr);call qqb0113(qqbr);call qqb0114(qqbr);call qqb0115(qqbr);call qqb0116(qqbr);call
     & qqb0117(qqbr);call qqb0118(qqbr);call qqb0119(qqbr);call qqb0120(qqbr)
          call qqb0121(qqbr);call qqb0122(qqbr);call qqb0123(qqbr);call qqb0124(qqbr);call qqb0125(qqbr);call qqb0126(qqbr);call
     & qqb0127(qqbr);call qqb0128(qqbr);call qqb0129(qqbr);call qqb0130(qqbr)
          call qqb0131(qqbr);call qqb0132(qqbr);call qqb0133(qqbr);call qqb0134(qqbr);call qqb0135(qqbr);call qqb0136(qqbr);call
     & qqb0137(qqbr);call qqb0138(qqbr);call qqb0139(qqbr);call qqb0140(qqbr)
          call qqb0141(qqbr);call qqb0142(qqbr);call qqb0143(qqbr);call qqb0144(qqbr);call qqb0145(qqbr);call qqb0146(qqbr);call
     & qqb0147(qqbr);call qqb0148(qqbr);call qqb0149(qqbr);call qqb0150(qqbr)
          call qqb0151(qqbr);call qqb0152(qqbr);call qqb0153(qqbr);call qqb0154(qqbr);call qqb0155(qqbr);call qqb0156(qqbr);call
     & qqb0157(qqbr);call qqb0158(qqbr);call qqb0159(qqbr);call qqb0160(qqbr)
          call qqb0161(qqbr);call qqb0162(qqbr);call qqb0163(qqbr);call qqb0164(qqbr);call qqb0165(qqbr);call qqb0166(qqbr);call
     & qqb0167(qqbr);call qqb0168(qqbr);call qqb0169(qqbr);call qqb0170(qqbr)
          call qqb0171(qqbr);call qqb0172(qqbr);call qqb0173(qqbr);call qqb0174(qqbr);call qqb0175(qqbr);call qqb0176(qqbr);call
     & qqb0177(qqbr);call qqb0178(qqbr);call qqb0179(qqbr);call qqb0180(qqbr)
          call qqb0181(qqbr);call qqb0182(qqbr);call qqb0183(qqbr);call qqb0184(qqbr);call qqb0185(qqbr);call qqb0186(qqbr);call
     & qqb0187(qqbr);call qqb0188(qqbr);call qqb0189(qqbr);call qqb0190(qqbr)
          call qqb0191(qqbr);call qqb0192(qqbr);call qqb0193(qqbr);call qqb0194(qqbr);call qqb0195(qqbr);call qqb0196(qqbr);call
     & qqb0197(qqbr);call qqb0198(qqbr);call qqb0199(qqbr);call qqb0200(qqbr)
          call qqb0201(qqbr);call qqb0202(qqbr);call qqb0203(qqbr);call qqb0204(qqbr);call qqb0205(qqbr);call qqb0206(qqbr);call
     & qqb0207(qqbr);call qqb0208(qqbr);call qqb0209(qqbr);call qqb0210(qqbr)
          call qqb0211(qqbr);call qqb0212(qqbr);call qqb0213(qqbr);call qqb0214(qqbr);call qqb0215(qqbr);call qqb0216(qqbr);call
     & qqb0217(qqbr);call qqb0218(qqbr);call qqb0219(qqbr);call qqb0220(qqbr)
          call qqb0221(qqbr);call qqb0222(qqbr);call qqb0223(qqbr);call qqb0224(qqbr);call qqb0225(qqbr);call qqb0226(qqbr);call
     & qqb0227(qqbr);call qqb0228(qqbr);call qqb0229(qqbr);call qqb0230(qqbr)
          call qqb0231(qqbr);call qqb0232(qqbr);call qqb0233(qqbr);call qqb0234(qqbr);call qqb0235(qqbr);call qqb0236(qqbr);call
     & qqb0237(qqbr);call qqb0238(qqbr);call qqb0239(qqbr);call qqb0240(qqbr)
          call qqb0241(qqbr);call qqb0242(qqbr);call qqb0243(qqbr);call qqb0244(qqbr);call qqb0245(qqbr);call qqb0246(qqbr);call
     & qqb0247(qqbr);call qqb0248(qqbr);call qqb0249(qqbr);call qqb0250(qqbr)
          call qqb0251(qqbr);call qqb0252(qqbr);call qqb0253(qqbr);call qqb0254(qqbr);call qqb0255(qqbr);call qqb0256(qqbr);call
     & qqb0257(qqbr);call qqb0258(qqbr);call qqb0259(qqbr);call qqb0260(qqbr)
          call qqb0261(qqbr);call qqb0262(qqbr);call qqb0263(qqbr);call qqb0264(qqbr);call qqb0265(qqbr);call qqb0266(qqbr);call
     & qqb0267(qqbr);call qqb0268(qqbr);call qqb0269(qqbr);call qqb0270(qqbr)
          call qqb0271(qqbr);call qqb0272(qqbr);call qqb0273(qqbr);call qqb0274(qqbr);call qqb0275(qqbr);call qqb0276(qqbr);call
     & qqb0277(qqbr);call qqb0278(qqbr);call qqb0279(qqbr);call qqb0280(qqbr)
          call qqb0281(qqbr);call qqb0282(qqbr);call qqb0283(qqbr);call qqb0284(qqbr);call qqb0285(qqbr);call qqb0286(qqbr);call
     & qqb0287(qqbr);call qqb0288(qqbr);call qqb0289(qqbr);call qqb0290(qqbr)
          call qqb0291(qqbr);call qqb0292(qqbr);call qqb0293(qqbr);call qqb0294(qqbr);call qqb0295(qqbr);call qqb0296(qqbr);call
     & qqb0297(qqbr);call qqb0298(qqbr);call qqb0299(qqbr);call qqb0300(qqbr)
          call qqb0301(qqbr);call qqb0302(qqbr);call qqb0303(qqbr);call qqb0304(qqbr);call qqb0305(qqbr);call qqb0306(qqbr);call
     & qqb0307(qqbr);call qqb0308(qqbr);call qqb0309(qqbr);call qqb0310(qqbr)
          call qqb0311(qqbr);call qqb0312(qqbr);call qqb0313(qqbr);call qqb0314(qqbr);call qqb0315(qqbr);call qqb0316(qqbr);call
     & qqb0317(qqbr);call qqb0318(qqbr);call qqb0319(qqbr);call qqb0320(qqbr)
          call qqb0321(qqbr);call qqb0322(qqbr);call qqb0323(qqbr);call qqb0324(qqbr);call qqb0325(qqbr);call qqb0326(qqbr);call
     & qqb0327(qqbr);call qqb0328(qqbr);call qqb0329(qqbr);call qqb0330(qqbr)
          call qqb0331(qqbr);call qqb0332(qqbr);call qqb0333(qqbr);call qqb0334(qqbr);call qqb0335(qqbr);call qqb0336(qqbr);call
     & qqb0337(qqbr);call qqb0338(qqbr);call qqb0339(qqbr);call qqb0340(qqbr)
          call qqb0341(qqbr);call qqb0342(qqbr);call qqb0343(qqbr);call qqb0344(qqbr);call qqb0345(qqbr);call qqb0346(qqbr);call
     & qqb0347(qqbr);call qqb0348(qqbr);call qqb0349(qqbr);call qqb0350(qqbr)
          call qqb0351(qqbr);call qqb0352(qqbr);call qqb0353(qqbr);call qqb0354(qqbr);call qqb0355(qqbr);call qqb0356(qqbr);call
     & qqb0357(qqbr);call qqb0358(qqbr);call qqb0359(qqbr);call qqb0360(qqbr)
          call qqb0361(qqbr);call qqb0362(qqbr);call qqb0363(qqbr);call qqb0364(qqbr);call qqb0365(qqbr);call qqb0366(qqbr);call
     & qqb0367(qqbr);call qqb0368(qqbr);call qqb0369(qqbr);call qqb0370(qqbr)
          call qqb0371(qqbr);call qqb0372(qqbr);call qqb0373(qqbr);call qqb0374(qqbr);call qqb0375(qqbr);call qqb0376(qqbr);call
     & qqb0377(qqbr);call qqb0378(qqbr);call qqb0379(qqbr);call qqb0380(qqbr)
          call qqb0381(qqbr);call qqb0382(qqbr);call qqb0383(qqbr);call qqb0384(qqbr);call qqb0385(qqbr);call qqb0386(qqbr);call
     & qqb0387(qqbr);call qqb0388(qqbr);call qqb0389(qqbr);call qqb0390(qqbr)
          call qqb0391(qqbr);call qqb0392(qqbr);call qqb0393(qqbr);call qqb0394(qqbr);call qqb0395(qqbr);call qqb0396(qqbr);call
     & qqb0397(qqbr);call qqb0398(qqbr);call qqb0399(qqbr);call qqb0400(qqbr)
          call qqb0401(qqbr);call qqb0402(qqbr);call qqb0403(qqbr);call qqb0404(qqbr);call qqb0405(qqbr);call qqb0406(qqbr);call
     & qqb0407(qqbr);call qqb0408(qqbr);call qqb0409(qqbr);call qqb0410(qqbr)
          call qqb0411(qqbr);call qqb0412(qqbr);call qqb0413(qqbr);call qqb0414(qqbr);call qqb0415(qqbr);call qqb0416(qqbr);call
     & qqb0417(qqbr);call qqb0418(qqbr);call qqb0419(qqbr);call qqb0420(qqbr)
          call qqb0421(qqbr);call qqb0422(qqbr);call qqb0423(qqbr);call qqb0424(qqbr);call qqb0425(qqbr);call qqb0426(qqbr);call
     & qqb0427(qqbr);call qqb0428(qqbr);call qqb0429(qqbr);call qqb0430(qqbr)
          call qqb0431(qqbr);call qqb0432(qqbr);call qqb0433(qqbr);call qqb0434(qqbr);call qqb0435(qqbr);call qqb0436(qqbr);call
     & qqb0437(qqbr);call qqb0438(qqbr);call qqb0439(qqbr);call qqb0440(qqbr)
          call qqb0441(qqbr);call qqb0442(qqbr);call qqb0443(qqbr);call qqb0444(qqbr);call qqb0445(qqbr);call qqb0446(qqbr);call
     & qqb0447(qqbr);call qqb0448(qqbr);call qqb0449(qqbr);call qqb0450(qqbr)
          call qqb0451(qqbr);call qqb0452(qqbr);call qqb0453(qqbr);call qqb0454(qqbr);call qqb0455(qqbr);call qqb0456(qqbr);call
     & qqb0457(qqbr);call qqb0458(qqbr);call qqb0459(qqbr);call qqb0460(qqbr)
          call qqb0461(qqbr);call qqb0462(qqbr);call qqb0463(qqbr);call qqb0464(qqbr);call qqb0465(qqbr);call qqb0466(qqbr);call
     & qqb0467(qqbr);call qqb0468(qqbr);call qqb0469(qqbr);call qqb0470(qqbr)
          call qqb0471(qqbr);call qqb0472(qqbr);call qqb0473(qqbr);call qqb0474(qqbr);call qqb0475(qqbr);call qqb0476(qqbr);call
     & qqb0477(qqbr);call qqb0478(qqbr);call qqb0479(qqbr);call qqb0480(qqbr)
          call qqb0481(qqbr);call qqb0482(qqbr);call qqb0483(qqbr);call qqb0484(qqbr);call qqb0485(qqbr);call qqb0486(qqbr);call
     & qqb0487(qqbr);call qqb0488(qqbr);call qqb0489(qqbr);call qqb0490(qqbr)
          call qqb0491(qqbr);call qqb0492(qqbr);call qqb0493(qqbr);call qqb0494(qqbr);call qqb0495(qqbr);call qqb0496(qqbr);call
     & qqb0497(qqbr);call qqb0498(qqbr);call qqb0499(qqbr);call qqb0500(qqbr)
          call qqb0501(qqbr);call qqb0502(qqbr);call qqb0503(qqbr);call qqb0504(qqbr);call qqb0505(qqbr);call qqb0506(qqbr);call
     & qqb0507(qqbr);call qqb0508(qqbr);call qqb0509(qqbr);call qqb0510(qqbr)
          call qqb0511(qqbr);call qqb0512(qqbr);call qqb0513(qqbr);call qqb0514(qqbr);call qqb0515(qqbr);call qqb0516(qqbr);call
     & qqb0517(qqbr);call qqb0518(qqbr);call qqb0519(qqbr);call qqb0520(qqbr)
          call qqb0521(qqbr);call qqb0522(qqbr);call qqb0523(qqbr);call qqb0524(qqbr);call qqb0525(qqbr);call qqb0526(qqbr);call
     & qqb0527(qqbr);call qqb0528(qqbr);call qqb0529(qqbr);call qqb0530(qqbr)
          call qqb0531(qqbr);call qqb0532(qqbr);call qqb0533(qqbr);call qqb0534(qqbr);call qqb0535(qqbr);call qqb0536(qqbr);call
     & qqb0537(qqbr);call qqb0538(qqbr);call qqb0539(qqbr);call qqb0540(qqbr)
          call qqb0541(qqbr);call qqb0542(qqbr);call qqb0543(qqbr);call qqb0544(qqbr);call qqb0545(qqbr);call qqb0546(qqbr);call
     & qqb0547(qqbr);call qqb0548(qqbr);call qqb0549(qqbr);call qqb0550(qqbr)
          call qqb0551(qqbr);call qqb0552(qqbr);call qqb0553(qqbr);call qqb0554(qqbr);call qqb0555(qqbr);call qqb0556(qqbr);call
     & qqb0557(qqbr);call qqb0558(qqbr);call qqb0559(qqbr);call qqb0560(qqbr)
          call qqb0561(qqbr);call qqb0562(qqbr);call qqb0563(qqbr);call qqb0564(qqbr);call qqb0565(qqbr);call qqb0566(qqbr);call
     & qqb0567(qqbr);call qqb0568(qqbr);call qqb0569(qqbr);call qqb0570(qqbr)
          call qqb0571(qqbr);call qqb0572(qqbr);call qqb0573(qqbr);call qqb0574(qqbr);call qqb0575(qqbr);call qqb0576(qqbr);call
     & qqb0577(qqbr);call qqb0578(qqbr);call qqb0579(qqbr);call qqb0580(qqbr)
          call qqb0581(qqbr);call qqb0582(qqbr);call qqb0583(qqbr);call qqb0584(qqbr);call qqb0585(qqbr);call qqb0586(qqbr);call
     & qqb0587(qqbr);call qqb0588(qqbr);call qqb0589(qqbr);call qqb0590(qqbr)
          call qqb0591(qqbr);call qqb0592(qqbr);call qqb0593(qqbr);call qqb0594(qqbr);call qqb0595(qqbr);call qqb0596(qqbr);call
     & qqb0597(qqbr);call qqb0598(qqbr);call qqb0599(qqbr);call qqb0600(qqbr)
          call qqb0601(qqbr);call qqb0602(qqbr);call qqb0603(qqbr);call qqb0604(qqbr);call qqb0605(qqbr);call qqb0606(qqbr);call
     & qqb0607(qqbr);call qqb0608(qqbr);call qqb0609(qqbr);call qqb0610(qqbr)
          call qqb0611(qqbr);call qqb0612(qqbr);call qqb0613(qqbr);call qqb0614(qqbr);call qqb0615(qqbr);call qqb0616(qqbr);call
     & qqb0617(qqbr);call qqb0618(qqbr);call qqb0619(qqbr);call qqb0620(qqbr)
          call qqb0621(qqbr);call qqb0622(qqbr);call qqb0623(qqbr);call qqb0624(qqbr);call qqb0625(qqbr);call qqb0626(qqbr);call
     & qqb0627(qqbr);call qqb0628(qqbr);call qqb0629(qqbr);call qqb0630(qqbr)
          call qqb0631(qqbr);call qqb0632(qqbr);call qqb0633(qqbr);call qqb0634(qqbr);call qqb0635(qqbr);call qqb0636(qqbr);call
     & qqb0637(qqbr);call qqb0638(qqbr);call qqb0639(qqbr);call qqb0640(qqbr)
          call qqb0641(qqbr);call qqb0642(qqbr);call qqb0643(qqbr);call qqb0644(qqbr);call qqb0645(qqbr);call qqb0646(qqbr);call
     & qqb0647(qqbr);call qqb0648(qqbr);call qqb0649(qqbr);call qqb0650(qqbr)
          call qqb0651(qqbr);call qqb0652(qqbr);call qqb0653(qqbr);call qqb0654(qqbr);call qqb0655(qqbr);call qqb0656(qqbr);call
     & qqb0657(qqbr);call qqb0658(qqbr);call qqb0659(qqbr);call qqb0660(qqbr)
          call qqb0661(qqbr);call qqb0662(qqbr);call qqb0663(qqbr);call qqb0664(qqbr);call qqb0665(qqbr);call qqb0666(qqbr);call
     & qqb0667(qqbr);call qqb0668(qqbr);call qqb0669(qqbr);call qqb0670(qqbr)
          call qqb0671(qqbr);call qqb0672(qqbr);call qqb0673(qqbr);call qqb0674(qqbr);call qqb0675(qqbr);call qqb0676(qqbr);call
     & qqb0677(qqbr);call qqb0678(qqbr);call qqb0679(qqbr);call qqb0680(qqbr)
          call qqb0681(qqbr);call qqb0682(qqbr);call qqb0683(qqbr);call qqb0684(qqbr);call qqb0685(qqbr);call qqb0686(qqbr);call
     & qqb0687(qqbr);call qqb0688(qqbr);call qqb0689(qqbr);call qqb0690(qqbr)
          call qqb0691(qqbr);call qqb0692(qqbr);call qqb0693(qqbr);call qqb0694(qqbr);call qqb0695(qqbr);call qqb0696(qqbr);call
     & qqb0697(qqbr);call qqb0698(qqbr);call qqb0699(qqbr);call qqb0700(qqbr)
          call qqb0701(qqbr);call qqb0702(qqbr);call qqb0703(qqbr);call qqb0704(qqbr);call qqb0705(qqbr);call qqb0706(qqbr);call
     & qqb0707(qqbr);call qqb0708(qqbr);call qqb0709(qqbr);call qqb0710(qqbr)
          call qqb0711(qqbr);call qqb0712(qqbr);call qqb0713(qqbr);call qqb0714(qqbr);call qqb0715(qqbr);call qqb0716(qqbr);call
     & qqb0717(qqbr);call qqb0718(qqbr);call qqb0719(qqbr);call qqb0720(qqbr)
          call qqb0721(qqbr);call qqb0722(qqbr);call qqb0723(qqbr);call qqb0724(qqbr);call qqb0725(qqbr);call qqb0726(qqbr);call
     & qqb0727(qqbr);call qqb0728(qqbr);call qqb0729(qqbr);call qqb0730(qqbr)
          call qqb0731(qqbr);call qqb0732(qqbr);call qqb0733(qqbr);call qqb0734(qqbr);call qqb0735(qqbr);call qqb0736(qqbr);call
     & qqb0737(qqbr);call qqb0738(qqbr);call qqb0739(qqbr);call qqb0740(qqbr)
          call qqb0741(qqbr);call qqb0742(qqbr);call qqb0743(qqbr);call qqb0744(qqbr);call qqb0745(qqbr);call qqb0746(qqbr);call
     & qqb0747(qqbr);call qqb0748(qqbr);call qqb0749(qqbr);call qqb0750(qqbr)
          call qqb0751(qqbr);call qqb0752(qqbr);call qqb0753(qqbr);call qqb0754(qqbr);call qqb0755(qqbr);call qqb0756(qqbr);call
     & qqb0757(qqbr);call qqb0758(qqbr);call qqb0759(qqbr);call qqb0760(qqbr)
          call qqb0761(qqbr);call qqb0762(qqbr);call qqb0763(qqbr);call qqb0764(qqbr);call qqb0765(qqbr);call qqb0766(qqbr);call
     & qqb0767(qqbr);call qqb0768(qqbr);call qqb0769(qqbr);call qqb0770(qqbr)
          call qqb0771(qqbr);call qqb0772(qqbr);call qqb0773(qqbr);call qqb0774(qqbr);call qqb0775(qqbr);call qqb0776(qqbr);call
     & qqb0777(qqbr);call qqb0778(qqbr);call qqb0779(qqbr);call qqb0780(qqbr)
          call qqb0781(qqbr);call qqb0782(qqbr);call qqb0783(qqbr);call qqb0784(qqbr);call qqb0785(qqbr);call qqb0786(qqbr);call
     & qqb0787(qqbr);call qqb0788(qqbr);call qqb0789(qqbr);call qqb0790(qqbr)
          call qqb0791(qqbr);call qqb0792(qqbr);call qqb0793(qqbr);call qqb0794(qqbr);call qqb0795(qqbr);call qqb0796(qqbr);call
     & qqb0797(qqbr);call qqb0798(qqbr);call qqb0799(qqbr);call qqb0800(qqbr)
          call qqb0801(qqbr);call qqb0802(qqbr);call qqb0803(qqbr);call qqb0804(qqbr);call qqb0805(qqbr);call qqb0806(qqbr);call
     & qqb0807(qqbr);call qqb0808(qqbr);call qqb0809(qqbr);call qqb0810(qqbr)
          call qqb0811(qqbr);call qqb0812(qqbr);call qqb0813(qqbr);call qqb0814(qqbr);call qqb0815(qqbr);call qqb0816(qqbr);call
     & qqb0817(qqbr);call qqb0818(qqbr);call qqb0819(qqbr);call qqb0820(qqbr)
          call qqb0821(qqbr);call qqb0822(qqbr);call qqb0823(qqbr);call qqb0824(qqbr);call qqb0825(qqbr);call qqb0826(qqbr);call
     & qqb0827(qqbr);call qqb0828(qqbr);call qqb0829(qqbr);call qqb0830(qqbr)
          call qqb0831(qqbr);call qqb0832(qqbr);call qqb0833(qqbr);call qqb0834(qqbr);call qqb0835(qqbr);call qqb0836(qqbr);call
     & qqb0837(qqbr);call qqb0838(qqbr);call qqb0839(qqbr);call qqb0840(qqbr)
          call qqb0841(qqbr);call qqb0842(qqbr);call qqb0843(qqbr);call qqb0844(qqbr);call qqb0845(qqbr);call qqb0846(qqbr);call
     & qqb0847(qqbr);call qqb0848(qqbr);call qqb0849(qqbr);call qqb0850(qqbr)
          call qqb0851(qqbr);call qqb0852(qqbr);call qqb0853(qqbr);call qqb0854(qqbr);call qqb0855(qqbr);call qqb0856(qqbr);call
     & qqb0857(qqbr);call qqb0858(qqbr);call qqb0859(qqbr);call qqb0860(qqbr)
          call qqb0861(qqbr);call qqb0862(qqbr);call qqb0863(qqbr);call qqb0864(qqbr);call qqb0865(qqbr);call qqb0866(qqbr);call
     & qqb0867(qqbr);call qqb0868(qqbr);call qqb0869(qqbr);call qqb0870(qqbr)
          call qqb0871(qqbr);call qqb0872(qqbr);call qqb0873(qqbr);call qqb0874(qqbr);call qqb0875(qqbr);call qqb0876(qqbr);call
     & qqb0877(qqbr);call qqb0878(qqbr);call qqb0879(qqbr);call qqb0880(qqbr)
          call qqb0881(qqbr);call qqb0882(qqbr);call qqb0883(qqbr);call qqb0884(qqbr);call qqb0885(qqbr);call qqb0886(qqbr);call
     & qqb0887(qqbr);call qqb0888(qqbr);call qqb0889(qqbr);call qqb0890(qqbr)
          call qqb0891(qqbr);call qqb0892(qqbr);call qqb0893(qqbr);call qqb0894(qqbr);call qqb0895(qqbr);call qqb0896(qqbr);call
     & qqb0897(qqbr);call qqb0898(qqbr);call qqb0899(qqbr);call qqb0900(qqbr)
          call qqb0901(qqbr);call qqb0902(qqbr);call qqb0903(qqbr);call qqb0904(qqbr);call qqb0905(qqbr);call qqb0906(qqbr);call
     & qqb0907(qqbr);call qqb0908(qqbr);call qqb0909(qqbr);call qqb0910(qqbr)
          call qqb0911(qqbr);call qqb0912(qqbr);call qqb0913(qqbr);call qqb0914(qqbr);call qqb0915(qqbr);call qqb0916(qqbr);call
     & qqb0917(qqbr);call qqb0918(qqbr);call qqb0919(qqbr);call qqb0920(qqbr)
          call qqb0921(qqbr);call qqb0922(qqbr);call qqb0923(qqbr);call qqb0924(qqbr);call qqb0925(qqbr);call qqb0926(qqbr);call
     & qqb0927(qqbr);call qqb0928(qqbr);call qqb0929(qqbr);call qqb0930(qqbr)
          call qqb0931(qqbr);call qqb0932(qqbr);call qqb0933(qqbr);call qqb0934(qqbr);call qqb0935(qqbr);call qqb0936(qqbr);call
     & qqb0937(qqbr);call qqb0938(qqbr);call qqb0939(qqbr);call qqb0940(qqbr)
          call qqb0941(qqbr);call qqb0942(qqbr);call qqb0943(qqbr);call qqb0944(qqbr);call qqb0945(qqbr);call qqb0946(qqbr);call
     & qqb0947(qqbr);call qqb0948(qqbr);call qqb0949(qqbr);call qqb0950(qqbr)
          call qqb0951(qqbr);call qqb0952(qqbr);call qqb0953(qqbr);call qqb0954(qqbr);call qqb0955(qqbr);call qqb0956(qqbr);call
     & qqb0957(qqbr);call qqb0958(qqbr);call qqb0959(qqbr);call qqb0960(qqbr)
          call qqb0961(qqbr);call qqb0962(qqbr);call qqb0963(qqbr);call qqb0964(qqbr);call qqb0965(qqbr);call qqb0966(qqbr);call
     & qqb0967(qqbr);call qqb0968(qqbr);call qqb0969(qqbr);call qqb0970(qqbr)
          call qqb0971(qqbr);call qqb0972(qqbr);call qqb0973(qqbr);call qqb0974(qqbr);call qqb0975(qqbr);call qqb0976(qqbr);call
     & qqb0977(qqbr);call qqb0978(qqbr);call qqb0979(qqbr);call qqb0980(qqbr)
          call qqb0981(qqbr);call qqb0982(qqbr);call qqb0983(qqbr);call qqb0984(qqbr);call qqb0985(qqbr);call qqb0986(qqbr);call
     & qqb0987(qqbr);call qqb0988(qqbr);call qqb0989(qqbr);call qqb0990(qqbr)
          call qqb0991(qqbr);call qqb0992(qqbr);call qqb0993(qqbr);call qqb0994(qqbr);call qqb0995(qqbr);call qqb0996(qqbr);call
     & qqb0997(qqbr);call qqb0998(qqbr);call qqb0999(qqbr);call qqb1000(qqbr)
          call qqb1001(qqbr);call qqb1002(qqbr);call qqb1003(qqbr);call qqb1004(qqbr);call qqb1005(qqbr);call qqb1006(qqbr);call
     & qqb1007(qqbr);call qqb1008(qqbr);call qqb1009(qqbr);call qqb1010(qqbr)
          call qqb1011(qqbr);call qqb1012(qqbr);call qqb1013(qqbr);call qqb1014(qqbr);call qqb1015(qqbr);call qqb1016(qqbr);call
     & qqb1017(qqbr);call qqb1018(qqbr);call qqb1019(qqbr);call qqb1020(qqbr)
          call qqb1021(qqbr);call qqb1022(qqbr);call qqb1023(qqbr);call qqb1024(qqbr);call qqb1025(qqbr);call qqb1026(qqbr);call
     & qqb1027(qqbr);call qqb1028(qqbr);call qqb1029(qqbr);call qqb1030(qqbr)
          call qqb1031(qqbr);call qqb1032(qqbr);call qqb1033(qqbr);call qqb1034(qqbr);call qqb1035(qqbr);call qqb1036(qqbr);call
     & qqb1037(qqbr);call qqb1038(qqbr);call qqb1039(qqbr);call qqb1040(qqbr)
          call qqb1041(qqbr);call qqb1042(qqbr);call qqb1043(qqbr);call qqb1044(qqbr);call qqb1045(qqbr);call qqb1046(qqbr);call
     & qqb1047(qqbr);call qqb1048(qqbr);call qqb1049(qqbr);call qqb1050(qqbr)
          call qqb1051(qqbr);call qqb1052(qqbr);call qqb1053(qqbr);call qqb1054(qqbr);call qqb1055(qqbr);call qqb1056(qqbr);call
     & qqb1057(qqbr);call qqb1058(qqbr);call qqb1059(qqbr);call qqb1060(qqbr)
          call qqb1061(qqbr);call qqb1062(qqbr);call qqb1063(qqbr);call qqb1064(qqbr);call qqb1065(qqbr);call qqb1066(qqbr);call
     & qqb1067(qqbr);call qqb1068(qqbr);call qqb1069(qqbr);call qqb1070(qqbr)
          call qqb1071(qqbr);call qqb1072(qqbr);call qqb1073(qqbr);call qqb1074(qqbr);call qqb1075(qqbr);call qqb1076(qqbr);call
     & qqb1077(qqbr);call qqb1078(qqbr);call qqb1079(qqbr);call qqb1080(qqbr)
          call qqb1081(qqbr);call qqb1082(qqbr);call qqb1083(qqbr);call qqb1084(qqbr);call qqb1085(qqbr);call qqb1086(qqbr);call
     & qqb1087(qqbr);call qqb1088(qqbr);call qqb1089(qqbr);call qqb1090(qqbr)
          call qqb1091(qqbr);call qqb1092(qqbr);call qqb1093(qqbr);call qqb1094(qqbr);call qqb1095(qqbr);call qqb1096(qqbr);call
     & qqb1097(qqbr);call qqb1098(qqbr);call qqb1099(qqbr);call qqb1100(qqbr)
          call qqb1101(qqbr);call qqb1102(qqbr);call qqb1103(qqbr);call qqb1104(qqbr);call qqb1105(qqbr);call qqb1106(qqbr);call
     & qqb1107(qqbr);call qqb1108(qqbr);call qqb1109(qqbr);call qqb1110(qqbr)
          call qqb1111(qqbr);call qqb1112(qqbr);call qqb1113(qqbr);call qqb1114(qqbr);call qqb1115(qqbr);call qqb1116(qqbr);call
     & qqb1117(qqbr);call qqb1118(qqbr);call qqb1119(qqbr);call qqb1120(qqbr)
          call qqb1121(qqbr);call qqb1122(qqbr);call qqb1123(qqbr);call qqb1124(qqbr);call qqb1125(qqbr);call qqb1126(qqbr);call
     & qqb1127(qqbr);call qqb1128(qqbr);call qqb1129(qqbr);call qqb1130(qqbr)
          call qqb1131(qqbr);call qqb1132(qqbr);call qqb1133(qqbr);call qqb1134(qqbr);call qqb1135(qqbr);call qqb1136(qqbr);call
     & qqb1137(qqbr);call qqb1138(qqbr);call qqb1139(qqbr);call qqb1140(qqbr)
          call qqb1141(qqbr);call qqb1142(qqbr);call qqb1143(qqbr);call qqb1144(qqbr);call qqb1145(qqbr);call qqb1146(qqbr);call
     & qqb1147(qqbr);call qqb1148(qqbr);call qqb1149(qqbr);call qqb1150(qqbr)
          call qqb1151(qqbr);call qqb1152(qqbr);call qqb1153(qqbr);call qqb1154(qqbr);call qqb1155(qqbr);call qqb1156(qqbr);call
     & qqb1157(qqbr);call qqb1158(qqbr);call qqb1159(qqbr);call qqb1160(qqbr)
          call qqb1161(qqbr);call qqb1162(qqbr);call qqb1163(qqbr);call qqb1164(qqbr);call qqb1165(qqbr);call qqb1166(qqbr);call
     & qqb1167(qqbr);call qqb1168(qqbr);call qqb1169(qqbr);call qqb1170(qqbr)
          call qqb1171(qqbr);call qqb1172(qqbr);call qqb1173(qqbr);call qqb1174(qqbr);call qqb1175(qqbr);call qqb1176(qqbr);call
     & qqb1177(qqbr);call qqb1178(qqbr);call qqb1179(qqbr);call qqb1180(qqbr)
          call qqb1181(qqbr);call qqb1182(qqbr);call qqb1183(qqbr);call qqb1184(qqbr);call qqb1185(qqbr);call qqb1186(qqbr);call
     & qqb1187(qqbr);call qqb1188(qqbr);call qqb1189(qqbr);call qqb1190(qqbr)
          call qqb1191(qqbr);call qqb1192(qqbr);call qqb1193(qqbr);call qqb1194(qqbr);call qqb1195(qqbr);call qqb1196(qqbr);call
     & qqb1197(qqbr);call qqb1198(qqbr);call qqb1199(qqbr);call qqb1200(qqbr)
          call qqb1201(qqbr);call qqb1202(qqbr);call qqb1203(qqbr);call qqb1204(qqbr);call qqb1205(qqbr);call qqb1206(qqbr);call
     & qqb1207(qqbr);call qqb1208(qqbr);call qqb1209(qqbr);call qqb1210(qqbr)
          call qqb1211(qqbr);call qqb1212(qqbr);call qqb1213(qqbr);call qqb1214(qqbr);call qqb1215(qqbr);call qqb1216(qqbr);call
     & qqb1217(qqbr);call qqb1218(qqbr);call qqb1219(qqbr);call qqb1220(qqbr)
          call qqb1221(qqbr);call qqb1222(qqbr);call qqb1223(qqbr);call qqb1224(qqbr);call qqb1225(qqbr);call qqb1226(qqbr);call
     & qqb1227(qqbr);call qqb1228(qqbr);call qqb1229(qqbr);call qqb1230(qqbr)
          call qqb1231(qqbr);call qqb1232(qqbr);call qqb1233(qqbr);call qqb1234(qqbr);call qqb1235(qqbr);call qqb1236(qqbr);call
     & qqb1237(qqbr);call qqb1238(qqbr);call qqb1239(qqbr);call qqb1240(qqbr)
          call qqb1241(qqbr);call qqb1242(qqbr);call qqb1243(qqbr);call qqb1244(qqbr);call qqb1245(qqbr);call qqb1246(qqbr);call
     & qqb1247(qqbr);call qqb1248(qqbr);call qqb1249(qqbr);call qqb1250(qqbr)
          call qqb1251(qqbr);call qqb1252(qqbr);call qqb1253(qqbr);call qqb1254(qqbr);call qqb1255(qqbr);call qqb1256(qqbr);call
     & qqb1257(qqbr);call qqb1258(qqbr);call qqb1259(qqbr);call qqb1260(qqbr)
          call qqb1261(qqbr);call qqb1262(qqbr);call qqb1263(qqbr);call qqb1264(qqbr);call qqb1265(qqbr);call qqb1266(qqbr);call
     & qqb1267(qqbr);call qqb1268(qqbr);call qqb1269(qqbr);call qqb1270(qqbr)
          call qqb1271(qqbr);call qqb1272(qqbr);call qqb1273(qqbr);call qqb1274(qqbr);call qqb1275(qqbr);call qqb1276(qqbr);call
     & qqb1277(qqbr);call qqb1278(qqbr);call qqb1279(qqbr);call qqb1280(qqbr)
          call qqb1281(qqbr);call qqb1282(qqbr);call qqb1283(qqbr);call qqb1284(qqbr);call qqb1285(qqbr);call qqb1286(qqbr);call
     & qqb1287(qqbr);call qqb1288(qqbr);call qqb1289(qqbr);call qqb1290(qqbr)
          call qqb1291(qqbr);call qqb1292(qqbr);call qqb1293(qqbr);call qqb1294(qqbr);call qqb1295(qqbr);call qqb1296(qqbr);call
     & qqb1297(qqbr);call qqb1298(qqbr);call qqb1299(qqbr);call qqb1300(qqbr)
          call qqb1301(qqbr);call qqb1302(qqbr);call qqb1303(qqbr);call qqb1304(qqbr);call qqb1305(qqbr);call qqb1306(qqbr);call
     & qqb1307(qqbr);call qqb1308(qqbr);call qqb1309(qqbr);call qqb1310(qqbr)
          call qqb1311(qqbr);call qqb1312(qqbr);call qqb1313(qqbr);call qqb1314(qqbr);call qqb1315(qqbr);call qqb1316(qqbr);call
     & qqb1317(qqbr);call qqb1318(qqbr);call qqb1319(qqbr);call qqb1320(qqbr)
          call qqb1321(qqbr);call qqb1322(qqbr);call qqb1323(qqbr);call qqb1324(qqbr);call qqb1325(qqbr);call qqb1326(qqbr);call
     & qqb1327(qqbr);call qqb1328(qqbr);call qqb1329(qqbr);call qqb1330(qqbr)
          call qqb1331(qqbr);call qqb1332(qqbr);call qqb1333(qqbr);call qqb1334(qqbr);call qqb1335(qqbr);call qqb1336(qqbr);call
     & qqb1337(qqbr);call qqb1338(qqbr);call qqb1339(qqbr);call qqb1340(qqbr)
          call qqb1341(qqbr);call qqb1342(qqbr);call qqb1343(qqbr);call qqb1344(qqbr);call qqb1345(qqbr);call qqb1346(qqbr);call
     & qqb1347(qqbr);call qqb1348(qqbr);call qqb1349(qqbr);call qqb1350(qqbr)
          call qqb1351(qqbr);call qqb1352(qqbr);call qqb1353(qqbr);call qqb1354(qqbr);call qqb1355(qqbr);call qqb1356(qqbr);call
     & qqb1357(qqbr);call qqb1358(qqbr);call qqb1359(qqbr);call qqb1360(qqbr)
          call qqb1361(qqbr);call qqb1362(qqbr);call qqb1363(qqbr);call qqb1364(qqbr);call qqb1365(qqbr);call qqb1366(qqbr);call
     & qqb1367(qqbr);call qqb1368(qqbr);call qqb1369(qqbr);call qqb1370(qqbr)
          call qqb1371(qqbr);call qqb1372(qqbr);call qqb1373(qqbr);call qqb1374(qqbr);call qqb1375(qqbr);call qqb1376(qqbr);call
     & qqb1377(qqbr);call qqb1378(qqbr);call qqb1379(qqbr);call qqb1380(qqbr)
          call qqb1381(qqbr);call qqb1382(qqbr);call qqb1383(qqbr);call qqb1384(qqbr);call qqb1385(qqbr);call qqb1386(qqbr);call
     & qqb1387(qqbr);call qqb1388(qqbr);call qqb1389(qqbr);call qqb1390(qqbr)
          call qqb1391(qqbr);call qqb1392(qqbr);call qqb1393(qqbr);call qqb1394(qqbr);call qqb1395(qqbr);call qqb1396(qqbr);call
     & qqb1397(qqbr);call qqb1398(qqbr);call qqb1399(qqbr);call qqb1400(qqbr)
          call qqb1401(qqbr);call qqb1402(qqbr);call qqb1403(qqbr);call qqb1404(qqbr);call qqb1405(qqbr);call qqb1406(qqbr);call
     & qqb1407(qqbr);call qqb1408(qqbr);call qqb1409(qqbr);call qqb1410(qqbr)
          call qqb1411(qqbr);call qqb1412(qqbr);call qqb1413(qqbr);call qqb1414(qqbr);call qqb1415(qqbr);call qqb1416(qqbr);call
     & qqb1417(qqbr);call qqb1418(qqbr);call qqb1419(qqbr);call qqb1420(qqbr)
          call qqb1421(qqbr);call qqb1422(qqbr);call qqb1423(qqbr);call qqb1424(qqbr);call qqb1425(qqbr);call qqb1426(qqbr);call
     & qqb1427(qqbr);call qqb1428(qqbr);call qqb1429(qqbr);call qqb1430(qqbr)
          call qqb1431(qqbr);call qqb1432(qqbr);call qqb1433(qqbr);call qqb1434(qqbr);call qqb1435(qqbr);call qqb1436(qqbr);call
     & qqb1437(qqbr);call qqb1438(qqbr);call qqb1439(qqbr);call qqb1440(qqbr)
          call qqb1441(qqbr);call qqb1442(qqbr);call qqb1443(qqbr);call qqb1444(qqbr);call qqb1445(qqbr);call qqb1446(qqbr);call
     & qqb1447(qqbr);call qqb1448(qqbr);call qqb1449(qqbr);call qqb1450(qqbr)
          call qqb1451(qqbr);call qqb1452(qqbr);call qqb1453(qqbr);call qqb1454(qqbr);call qqb1455(qqbr);call qqb1456(qqbr);call
     & qqb1457(qqbr);call qqb1458(qqbr);call qqb1459(qqbr);call qqb1460(qqbr)
          call qqb1461(qqbr);call qqb1462(qqbr);call qqb1463(qqbr);call qqb1464(qqbr);call qqb1465(qqbr);call qqb1466(qqbr);call
     & qqb1467(qqbr);call qqb1468(qqbr);call qqb1469(qqbr);call qqb1470(qqbr)
          call qqb1471(qqbr);call qqb1472(qqbr);call qqb1473(qqbr);call qqb1474(qqbr);call qqb1475(qqbr);call qqb1476(qqbr);call
     & qqb1477(qqbr);call qqb1478(qqbr);call qqb1479(qqbr);call qqb1480(qqbr)
          call qqb1481(qqbr);call qqb1482(qqbr);call qqb1483(qqbr);call qqb1484(qqbr);call qqb1485(qqbr);call qqb1486(qqbr);call
     & qqb1487(qqbr);call qqb1488(qqbr);call qqb1489(qqbr);call qqb1490(qqbr)
          call qqb1491(qqbr);call qqb1492(qqbr);call qqb1493(qqbr);call qqb1494(qqbr);call qqb1495(qqbr);call qqb1496(qqbr);call
     & qqb1497(qqbr);call qqb1498(qqbr);call qqb1499(qqbr);call qqb1500(qqbr)
          call qqb1501(qqbr);call qqb1502(qqbr);call qqb1503(qqbr);call qqb1504(qqbr);call qqb1505(qqbr);call qqb1506(qqbr);call
     & qqb1507(qqbr);call qqb1508(qqbr);call qqb1509(qqbr);call qqb1510(qqbr)
          call qqb1511(qqbr);call qqb1512(qqbr);call qqb1513(qqbr);call qqb1514(qqbr);call qqb1515(qqbr);call qqb1516(qqbr);call
     & qqb1517(qqbr);call qqb1518(qqbr);call qqb1519(qqbr);call qqb1520(qqbr)
          call qqb1521(qqbr);call qqb1522(qqbr);call qqb1523(qqbr);call qqb1524(qqbr);call qqb1525(qqbr);call qqb1526(qqbr);call
     & qqb1527(qqbr);call qqb1528(qqbr);call qqb1529(qqbr);call qqb1530(qqbr)
          call qqb1531(qqbr);call qqb1532(qqbr);call qqb1533(qqbr);call qqb1534(qqbr);call qqb1535(qqbr);call qqb1536(qqbr);call
     & qqb1537(qqbr);call qqb1538(qqbr);call qqb1539(qqbr);call qqb1540(qqbr)
          call qqb1541(qqbr);call qqb1542(qqbr);call qqb1543(qqbr);call qqb1544(qqbr);call qqb1545(qqbr);call qqb1546(qqbr);call
     & qqb1547(qqbr);call qqb1548(qqbr);call qqb1549(qqbr);call qqb1550(qqbr)
          call qqb1551(qqbr);call qqb1552(qqbr);call qqb1553(qqbr);call qqb1554(qqbr);call qqb1555(qqbr);call qqb1556(qqbr);call
     & qqb1557(qqbr);call qqb1558(qqbr);call qqb1559(qqbr);call qqb1560(qqbr)
          call qqb1561(qqbr);call qqb1562(qqbr);call qqb1563(qqbr);call qqb1564(qqbr);call qqb1565(qqbr);call qqb1566(qqbr);call
     & qqb1567(qqbr);call qqb1568(qqbr);call qqb1569(qqbr);call qqb1570(qqbr)
          call qqb1571(qqbr);call qqb1572(qqbr);call qqb1573(qqbr);call qqb1574(qqbr);call qqb1575(qqbr);call qqb1576(qqbr);call
     & qqb1577(qqbr);call qqb1578(qqbr);call qqb1579(qqbr);call qqb1580(qqbr)
          call qqb1581(qqbr);call qqb1582(qqbr);call qqb1583(qqbr);call qqb1584(qqbr);call qqb1585(qqbr);call qqb1586(qqbr);call
     & qqb1587(qqbr);call qqb1588(qqbr);call qqb1589(qqbr);call qqb1590(qqbr)
          call qqb1591(qqbr);call qqb1592(qqbr);call qqb1593(qqbr);call qqb1594(qqbr);call qqb1595(qqbr);call qqb1596(qqbr);call
     & qqb1597(qqbr);call qqb1598(qqbr);call qqb1599(qqbr);call qqb1600(qqbr)
          call qqb1601(qqbr);call qqb1602(qqbr);call qqb1603(qqbr);call qqb1604(qqbr);call qqb1605(qqbr);call qqb1606(qqbr);call
     & qqb1607(qqbr);call qqb1608(qqbr);call qqb1609(qqbr);call qqb1610(qqbr)
          call qqb1611(qqbr);call qqb1612(qqbr);call qqb1613(qqbr);call qqb1614(qqbr);call qqb1615(qqbr);call qqb1616(qqbr);call
     & qqb1617(qqbr);call qqb1618(qqbr);call qqb1619(qqbr);call qqb1620(qqbr)
          call qqb1621(qqbr);call qqb1622(qqbr);call qqb1623(qqbr);call qqb1624(qqbr);call qqb1625(qqbr);call qqb1626(qqbr);call
     & qqb1627(qqbr);call qqb1628(qqbr);call qqb1629(qqbr);call qqb1630(qqbr)
          call qqb1631(qqbr);call qqb1632(qqbr);call qqb1633(qqbr);call qqb1634(qqbr);call qqb1635(qqbr);call qqb1636(qqbr);call
     & qqb1637(qqbr);call qqb1638(qqbr);call qqb1639(qqbr);call qqb1640(qqbr)
          call qqb1641(qqbr);call qqb1642(qqbr);call qqb1643(qqbr);call qqb1644(qqbr);call qqb1645(qqbr);call qqb1646(qqbr);call
     & qqb1647(qqbr);call qqb1648(qqbr);call qqb1649(qqbr);call qqb1650(qqbr)
          call qqb1651(qqbr);call qqb1652(qqbr);call qqb1653(qqbr);call qqb1654(qqbr);call qqb1655(qqbr);call qqb1656(qqbr);call
     & qqb1657(qqbr);call qqb1658(qqbr);call qqb1659(qqbr);call qqb1660(qqbr)
          call qqb1661(qqbr);call qqb1662(qqbr);call qqb1663(qqbr);call qqb1664(qqbr);call qqb1665(qqbr);call qqb1666(qqbr);call
     & qqb1667(qqbr);call qqb1668(qqbr);call qqb1669(qqbr);call qqb1670(qqbr)
          call qqb1671(qqbr);call qqb1672(qqbr);call qqb1673(qqbr);call qqb1674(qqbr);call qqb1675(qqbr);call qqb1676(qqbr);call
     & qqb1677(qqbr);call qqb1678(qqbr);call qqb1679(qqbr);call qqb1680(qqbr)
          call qqb1681(qqbr);call qqb1682(qqbr);call qqb1683(qqbr);call qqb1684(qqbr);call qqb1685(qqbr);call qqb1686(qqbr);call
     & qqb1687(qqbr);call qqb1688(qqbr);call qqb1689(qqbr);call qqb1690(qqbr)
          call qqb1691(qqbr);call qqb1692(qqbr);call qqb1693(qqbr);call qqb1694(qqbr);call qqb1695(qqbr);call qqb1696(qqbr);call
     & qqb1697(qqbr);call qqb1698(qqbr);call qqb1699(qqbr);call qqb1700(qqbr)
          call qqb1701(qqbr);call qqb1702(qqbr);call qqb1703(qqbr);call qqb1704(qqbr);call qqb1705(qqbr);call qqb1706(qqbr);call
     & qqb1707(qqbr);call qqb1708(qqbr);call qqb1709(qqbr);call qqb1710(qqbr)
          call qqb1711(qqbr);call qqb1712(qqbr);call qqb1713(qqbr);call qqb1714(qqbr);call qqb1715(qqbr);call qqb1716(qqbr);call
     & qqb1717(qqbr);call qqb1718(qqbr);call qqb1719(qqbr);call qqb1720(qqbr)
          call qqb1721(qqbr);call qqb1722(qqbr);call qqb1723(qqbr);call qqb1724(qqbr);call qqb1725(qqbr);call qqb1726(qqbr);call
     & qqb1727(qqbr);call qqb1728(qqbr);call qqb1729(qqbr);call qqb1730(qqbr)
          call qqb1731(qqbr);call qqb1732(qqbr);call qqb1733(qqbr);call qqb1734(qqbr);call qqb1735(qqbr);call qqb1736(qqbr);call
     & qqb1737(qqbr);call qqb1738(qqbr);call qqb1739(qqbr);call qqb1740(qqbr)
          call qqb1741(qqbr);call qqb1742(qqbr);call qqb1743(qqbr);call qqb1744(qqbr);call qqb1745(qqbr);call qqb1746(qqbr);call
     & qqb1747(qqbr);call qqb1748(qqbr);call qqb1749(qqbr);call qqb1750(qqbr)
          call qqb1751(qqbr);call qqb1752(qqbr);call qqb1753(qqbr);call qqb1754(qqbr);call qqb1755(qqbr);call qqb1756(qqbr);call
     & qqb1757(qqbr);call qqb1758(qqbr);call qqb1759(qqbr);call qqb1760(qqbr)
          call qqb1761(qqbr);call qqb1762(qqbr);call qqb1763(qqbr);call qqb1764(qqbr);call qqb1765(qqbr);call qqb1766(qqbr);call
     & qqb1767(qqbr);call qqb1768(qqbr);call qqb1769(qqbr);call qqb1770(qqbr)
          call qqb1771(qqbr);call qqb1772(qqbr);call qqb1773(qqbr);call qqb1774(qqbr);call qqb1775(qqbr);call qqb1776(qqbr);call
     & qqb1777(qqbr);call qqb1778(qqbr);call qqb1779(qqbr);call qqb1780(qqbr)
          call qqb1781(qqbr);call qqb1782(qqbr);call qqb1783(qqbr);call qqb1784(qqbr);call qqb1785(qqbr);call qqb1786(qqbr);call
     & qqb1787(qqbr);call qqb1788(qqbr);call qqb1789(qqbr);call qqb1790(qqbr)
          call qqb1791(qqbr);call qqb1792(qqbr);call qqb1793(qqbr);call qqb1794(qqbr);call qqb1795(qqbr);call qqb1796(qqbr);call
     & qqb1797(qqbr);call qqb1798(qqbr);call qqb1799(qqbr);call qqb1800(qqbr)
          call qqb1801(qqbr);call qqb1802(qqbr);call qqb1803(qqbr);call qqb1804(qqbr);call qqb1805(qqbr);call qqb1806(qqbr);call
     & qqb1807(qqbr);call qqb1808(qqbr);call qqb1809(qqbr);call qqb1810(qqbr)
          call qqb1811(qqbr);call qqb1812(qqbr);call qqb1813(qqbr);call qqb1814(qqbr);call qqb1815(qqbr);call qqb1816(qqbr);call
     & qqb1817(qqbr);call qqb1818(qqbr);call qqb1819(qqbr);call qqb1820(qqbr)
          call qqb1821(qqbr);call qqb1822(qqbr);call qqb1823(qqbr);call qqb1824(qqbr);call qqb1825(qqbr);call qqb1826(qqbr);call
     & qqb1827(qqbr);call qqb1828(qqbr);call qqb1829(qqbr);call qqb1830(qqbr)
          call qqb1831(qqbr);call qqb1832(qqbr);call qqb1833(qqbr);call qqb1834(qqbr);call qqb1835(qqbr);call qqb1836(qqbr);call
     & qqb1837(qqbr);call qqb1838(qqbr);call qqb1839(qqbr);call qqb1840(qqbr)
          call qqb1841(qqbr);call qqb1842(qqbr);call qqb1843(qqbr);call qqb1844(qqbr);call qqb1845(qqbr);call qqb1846(qqbr);call
     & qqb1847(qqbr);call qqb1848(qqbr);call qqb1849(qqbr);call qqb1850(qqbr)
          call qqb1851(qqbr);call qqb1852(qqbr);call qqb1853(qqbr);call qqb1854(qqbr);call qqb1855(qqbr);call qqb1856(qqbr);call
     & qqb1857(qqbr);call qqb1858(qqbr);call qqb1859(qqbr);call qqb1860(qqbr)
          call qqb1861(qqbr);call qqb1862(qqbr);call qqb1863(qqbr);call qqb1864(qqbr);call qqb1865(qqbr);call qqb1866(qqbr);call
     & qqb1867(qqbr);call qqb1868(qqbr);call qqb1869(qqbr);call qqb1870(qqbr)
          call qqb1871(qqbr);call qqb1872(qqbr);call qqb1873(qqbr);call qqb1874(qqbr);call qqb1875(qqbr);call qqb1876(qqbr);call
     & qqb1877(qqbr);call qqb1878(qqbr);call qqb1879(qqbr);call qqb1880(qqbr)
          call qqb1881(qqbr);call qqb1882(qqbr);call qqb1883(qqbr);call qqb1884(qqbr);call qqb1885(qqbr);call qqb1886(qqbr);call
     & qqb1887(qqbr);call qqb1888(qqbr);call qqb1889(qqbr);call qqb1890(qqbr)
          call qqb1891(qqbr);call qqb1892(qqbr);call qqb1893(qqbr);call qqb1894(qqbr);call qqb1895(qqbr);call qqb1896(qqbr);call
     & qqb1897(qqbr);call qqb1898(qqbr);call qqb1899(qqbr);call qqb1900(qqbr)
          call qqb1901(qqbr);call qqb1902(qqbr);call qqb1903(qqbr);call qqb1904(qqbr);call qqb1905(qqbr);call qqb1906(qqbr);call
     & qqb1907(qqbr);call qqb1908(qqbr);call qqb1909(qqbr);call qqb1910(qqbr)
          call qqb1911(qqbr);call qqb1912(qqbr);call qqb1913(qqbr);call qqb1914(qqbr);call qqb1915(qqbr);call qqb1916(qqbr);call
     & qqb1917(qqbr);call qqb1918(qqbr);call qqb1919(qqbr);call qqb1920(qqbr)
          call qqb1921(qqbr);call qqb1922(qqbr);call qqb1923(qqbr);call qqb1924(qqbr);call qqb1925(qqbr);call qqb1926(qqbr);call
     & qqb1927(qqbr);call qqb1928(qqbr);call qqb1929(qqbr);call qqb1930(qqbr)
          call qqb1931(qqbr);call qqb1932(qqbr);call qqb1933(qqbr);call qqb1934(qqbr);call qqb1935(qqbr);call qqb1936(qqbr);call
     & qqb1937(qqbr);call qqb1938(qqbr);call qqb1939(qqbr);call qqb1940(qqbr)
          call qqb1941(qqbr);call qqb1942(qqbr);call qqb1943(qqbr);call qqb1944(qqbr);call qqb1945(qqbr);call qqb1946(qqbr);call
     & qqb1947(qqbr);call qqb1948(qqbr);call qqb1949(qqbr);call qqb1950(qqbr)
          call qqb1951(qqbr);call qqb1952(qqbr);call qqb1953(qqbr);call qqb1954(qqbr);call qqb1955(qqbr);call qqb1956(qqbr);call
     & qqb1957(qqbr);call qqb1958(qqbr);call qqb1959(qqbr);call qqb1960(qqbr)
          call qqb1961(qqbr);call qqb1962(qqbr);call qqb1963(qqbr);call qqb1964(qqbr);call qqb1965(qqbr);call qqb1966(qqbr);call
     & qqb1967(qqbr);call qqb1968(qqbr);call qqb1969(qqbr);call qqb1970(qqbr)
          call qqb1971(qqbr);call qqb1972(qqbr);call qqb1973(qqbr);call qqb1974(qqbr);call qqb1975(qqbr);call qqb1976(qqbr);call
     & qqb1977(qqbr);call qqb1978(qqbr);call qqb1979(qqbr);call qqb1980(qqbr)
          call qqb1981(qqbr);call qqb1982(qqbr);call qqb1983(qqbr);call qqb1984(qqbr);call qqb1985(qqbr);call qqb1986(qqbr);call
     & qqb1987(qqbr);call qqb1988(qqbr);call qqb1989(qqbr);call qqb1990(qqbr)
          call qqb1991(qqbr);call qqb1992(qqbr);call qqb1993(qqbr);call qqb1994(qqbr);call qqb1995(qqbr);call qqb1996(qqbr);call
     & qqb1997(qqbr);call qqb1998(qqbr);call qqb1999(qqbr);call qqb2000(qqbr)
          call qqb2001(qqbr);call qqb2002(qqbr);call qqb2003(qqbr);call qqb2004(qqbr);call qqb2005(qqbr);call qqb2006(qqbr);call
     & qqb2007(qqbr);call qqb2008(qqbr);call qqb2009(qqbr);call qqb2010(qqbr)
          call qqb2011(qqbr);call qqb2012(qqbr);call qqb2013(qqbr);call qqb2014(qqbr);call qqb2015(qqbr);call qqb2016(qqbr);call
     & qqb2017(qqbr);call qqb2018(qqbr);call qqb2019(qqbr);call qqb2020(qqbr)
          call qqb2021(qqbr);call qqb2022(qqbr);call qqb2023(qqbr);call qqb2024(qqbr);call qqb2025(qqbr);call qqb2026(qqbr);call
     & qqb2027(qqbr);call qqb2028(qqbr);call qqb2029(qqbr);call qqb2030(qqbr)
          call qqb2031(qqbr);call qqb2032(qqbr);call qqb2033(qqbr);call qqb2034(qqbr);call qqb2035(qqbr);call qqb2036(qqbr);call
     & qqb2037(qqbr);call qqb2038(qqbr);call qqb2039(qqbr);call qqb2040(qqbr)
          call qqb2041(qqbr);call qqb2042(qqbr);call qqb2043(qqbr);call qqb2044(qqbr);call qqb2045(qqbr);call qqb2046(qqbr);call
     & qqb2047(qqbr);call qqb2048(qqbr);call qqb2049(qqbr);call qqb2050(qqbr)
          call qqb2051(qqbr);call qqb2052(qqbr);call qqb2053(qqbr);call qqb2054(qqbr);call qqb2055(qqbr);call qqb2056(qqbr);call
     & qqb2057(qqbr);call qqb2058(qqbr);call qqb2059(qqbr);call qqb2060(qqbr)
          call qqb2061(qqbr);call qqb2062(qqbr);call qqb2063(qqbr);call qqb2064(qqbr);call qqb2065(qqbr);call qqb2066(qqbr);call
     & qqb2067(qqbr);call qqb2068(qqbr);call qqb2069(qqbr);call qqb2070(qqbr)
          call qqb2071(qqbr);call qqb2072(qqbr);call qqb2073(qqbr);call qqb2074(qqbr);call qqb2075(qqbr);call qqb2076(qqbr);call
     & qqb2077(qqbr);call qqb2078(qqbr);call qqb2079(qqbr);call qqb2080(qqbr)
          call qqb2081(qqbr);call qqb2082(qqbr);call qqb2083(qqbr);call qqb2084(qqbr);call qqb2085(qqbr);call qqb2086(qqbr);call
     & qqb2087(qqbr);call qqb2088(qqbr);call qqb2089(qqbr);call qqb2090(qqbr)
          call qqb2091(qqbr);call qqb2092(qqbr);call qqb2093(qqbr);call qqb2094(qqbr);call qqb2095(qqbr);call qqb2096(qqbr);call
     & qqb2097(qqbr);call qqb2098(qqbr);call qqb2099(qqbr);call qqb2100(qqbr)
          call qqb2101(qqbr);call qqb2102(qqbr);call qqb2103(qqbr);call qqb2104(qqbr);call qqb2105(qqbr);call qqb2106(qqbr);call
     & qqb2107(qqbr);call qqb2108(qqbr);call qqb2109(qqbr);call qqb2110(qqbr)
          call qqb2111(qqbr);call qqb2112(qqbr);call qqb2113(qqbr);call qqb2114(qqbr);call qqb2115(qqbr);call qqb2116(qqbr);call
     & qqb2117(qqbr);call qqb2118(qqbr);call qqb2119(qqbr);call qqb2120(qqbr)
          call qqb2121(qqbr);call qqb2122(qqbr);call qqb2123(qqbr);call qqb2124(qqbr);call qqb2125(qqbr);call qqb2126(qqbr);call
     & qqb2127(qqbr);call qqb2128(qqbr);call qqb2129(qqbr);call qqb2130(qqbr)
          call qqb2131(qqbr);call qqb2132(qqbr);call qqb2133(qqbr);call qqb2134(qqbr);call qqb2135(qqbr);call qqb2136(qqbr);call
     & qqb2137(qqbr);call qqb2138(qqbr);call qqb2139(qqbr);call qqb2140(qqbr)
          call qqb2141(qqbr);call qqb2142(qqbr);call qqb2143(qqbr);call qqb2144(qqbr);call qqb2145(qqbr);call qqb2146(qqbr);call
     & qqb2147(qqbr);call qqb2148(qqbr);call qqb2149(qqbr);call qqb2150(qqbr)
          call qqb2151(qqbr);call qqb2152(qqbr);call qqb2153(qqbr);call qqb2154(qqbr);call qqb2155(qqbr);call qqb2156(qqbr);call
     & qqb2157(qqbr);call qqb2158(qqbr);call qqb2159(qqbr);call qqb2160(qqbr)
          call qqb2161(qqbr);call qqb2162(qqbr);call qqb2163(qqbr);call qqb2164(qqbr);call qqb2165(qqbr);call qqb2166(qqbr);call
     & qqb2167(qqbr);call qqb2168(qqbr);call qqb2169(qqbr);call qqb2170(qqbr)
          call qqb2171(qqbr);call qqb2172(qqbr);call qqb2173(qqbr);call qqb2174(qqbr);call qqb2175(qqbr);call qqb2176(qqbr);call
     & qqb2177(qqbr);call qqb2178(qqbr);call qqb2179(qqbr);call qqb2180(qqbr)
          call qqb2181(qqbr);call qqb2182(qqbr);call qqb2183(qqbr);call qqb2184(qqbr);call qqb2185(qqbr);call qqb2186(qqbr);call
     & qqb2187(qqbr);call qqb2188(qqbr);call qqb2189(qqbr);call qqb2190(qqbr)
          call qqb2191(qqbr);call qqb2192(qqbr);call qqb2193(qqbr);call qqb2194(qqbr);call qqb2195(qqbr);call qqb2196(qqbr);call
     & qqb2197(qqbr);call qqb2198(qqbr);call qqb2199(qqbr);call qqb2200(qqbr)
          call qqb2201(qqbr);call qqb2202(qqbr);call qqb2203(qqbr);call qqb2204(qqbr);call qqb2205(qqbr);call qqb2206(qqbr);call
     & qqb2207(qqbr);call qqb2208(qqbr);call qqb2209(qqbr);call qqb2210(qqbr)
          call qqb2211(qqbr);call qqb2212(qqbr);call qqb2213(qqbr);call qqb2214(qqbr);call qqb2215(qqbr);call qqb2216(qqbr);call
     & qqb2217(qqbr);call qqb2218(qqbr);call qqb2219(qqbr);call qqb2220(qqbr)
          call qqb2221(qqbr);call qqb2222(qqbr);call qqb2223(qqbr);call qqb2224(qqbr);call qqb2225(qqbr);call qqb2226(qqbr);call
     & qqb2227(qqbr);call qqb2228(qqbr);call qqb2229(qqbr);call qqb2230(qqbr)
          call qqb2231(qqbr);call qqb2232(qqbr);call qqb2233(qqbr);call qqb2234(qqbr);call qqb2235(qqbr);call qqb2236(qqbr);call
     & qqb2237(qqbr);call qqb2238(qqbr);call qqb2239(qqbr);call qqb2240(qqbr)
          call qqb2241(qqbr);call qqb2242(qqbr);call qqb2243(qqbr);call qqb2244(qqbr);call qqb2245(qqbr);call qqb2246(qqbr);call
     & qqb2247(qqbr);call qqb2248(qqbr);call qqb2249(qqbr);call qqb2250(qqbr)
          call qqb2251(qqbr);call qqb2252(qqbr);call qqb2253(qqbr);call qqb2254(qqbr);call qqb2255(qqbr);call qqb2256(qqbr);call
     & qqb2257(qqbr);call qqb2258(qqbr);call qqb2259(qqbr);call qqb2260(qqbr)
          call qqb2261(qqbr);call qqb2262(qqbr);call qqb2263(qqbr);call qqb2264(qqbr);call qqb2265(qqbr);call qqb2266(qqbr);call
     & qqb2267(qqbr);call qqb2268(qqbr);call qqb2269(qqbr);call qqb2270(qqbr)
          call qqb2271(qqbr);call qqb2272(qqbr);call qqb2273(qqbr);call qqb2274(qqbr);call qqb2275(qqbr);call qqb2276(qqbr);call
     & qqb2277(qqbr);call qqb2278(qqbr);call qqb2279(qqbr);call qqb2280(qqbr)
          call qqb2281(qqbr);call qqb2282(qqbr);call qqb2283(qqbr);call qqb2284(qqbr);call qqb2285(qqbr);call qqb2286(qqbr);call
     & qqb2287(qqbr);call qqb2288(qqbr);call qqb2289(qqbr);call qqb2290(qqbr)
          call qqb2291(qqbr);call qqb2292(qqbr);call qqb2293(qqbr);call qqb2294(qqbr);call qqb2295(qqbr);call qqb2296(qqbr);call
     & qqb2297(qqbr);call qqb2298(qqbr);call qqb2299(qqbr);call qqb2300(qqbr)
          call qqb2301(qqbr);call qqb2302(qqbr);call qqb2303(qqbr);call qqb2304(qqbr);call qqb2305(qqbr);call qqb2306(qqbr);call
     & qqb2307(qqbr);call qqb2308(qqbr);call qqb2309(qqbr);call qqb2310(qqbr)
          call qqb2311(qqbr);call qqb2312(qqbr);call qqb2313(qqbr);call qqb2314(qqbr);call qqb2315(qqbr);call qqb2316(qqbr);call
     & qqb2317(qqbr);call qqb2318(qqbr);call qqb2319(qqbr);call qqb2320(qqbr)
          call qqb2321(qqbr);call qqb2322(qqbr);call qqb2323(qqbr);call qqb2324(qqbr);call qqb2325(qqbr);call qqb2326(qqbr);call
     & qqb2327(qqbr);call qqb2328(qqbr);call qqb2329(qqbr);call qqb2330(qqbr)
          call qqb2331(qqbr);call qqb2332(qqbr);call qqb2333(qqbr);call qqb2334(qqbr);call qqb2335(qqbr);call qqb2336(qqbr);call
     & qqb2337(qqbr);call qqb2338(qqbr);call qqb2339(qqbr);call qqb2340(qqbr)
          call qqb2341(qqbr);call qqb2342(qqbr);call qqb2343(qqbr);call qqb2344(qqbr);call qqb2345(qqbr);call qqb2346(qqbr);call
     & qqb2347(qqbr);call qqb2348(qqbr);call qqb2349(qqbr);call qqb2350(qqbr)
          call qqb2351(qqbr);call qqb2352(qqbr);call qqb2353(qqbr);call qqb2354(qqbr);call qqb2355(qqbr);call qqb2356(qqbr);call
     & qqb2357(qqbr);call qqb2358(qqbr);call qqb2359(qqbr);call qqb2360(qqbr)
          call qqb2361(qqbr);call qqb2362(qqbr);call qqb2363(qqbr);call qqb2364(qqbr);call qqb2365(qqbr);call qqb2366(qqbr);call
     & qqb2367(qqbr);call qqb2368(qqbr);call qqb2369(qqbr);call qqb2370(qqbr)
          call qqb2371(qqbr);call qqb2372(qqbr);call qqb2373(qqbr);call qqb2374(qqbr);call qqb2375(qqbr);call qqb2376(qqbr);call
     & qqb2377(qqbr);call qqb2378(qqbr);call qqb2379(qqbr);call qqb2380(qqbr)
          call qqb2381(qqbr);call qqb2382(qqbr);call qqb2383(qqbr);call qqb2384(qqbr);call qqb2385(qqbr);call qqb2386(qqbr);call
     & qqb2387(qqbr);call qqb2388(qqbr);call qqb2389(qqbr);call qqb2390(qqbr)
          call qqb2391(qqbr);call qqb2392(qqbr);call qqb2393(qqbr);call qqb2394(qqbr);call qqb2395(qqbr);call qqb2396(qqbr);call
     & qqb2397(qqbr);call qqb2398(qqbr);call qqb2399(qqbr);call qqb2400(qqbr)
          call qqb2401(qqbr);call qqb2402(qqbr);call qqb2403(qqbr);call qqb2404(qqbr);call qqb2405(qqbr);call qqb2406(qqbr);call
     & qqb2407(qqbr);call qqb2408(qqbr);call qqb2409(qqbr);call qqb2410(qqbr)
          call qqb2411(qqbr);call qqb2412(qqbr);call qqb2413(qqbr);call qqb2414(qqbr);call qqb2415(qqbr);call qqb2416(qqbr);call
     & qqb2417(qqbr);call qqb2418(qqbr);call qqb2419(qqbr);call qqb2420(qqbr)
          call qqb2421(qqbr);call qqb2422(qqbr);call qqb2423(qqbr);call qqb2424(qqbr);call qqb2425(qqbr);call qqb2426(qqbr);call
     & qqb2427(qqbr);call qqb2428(qqbr);call qqb2429(qqbr);call qqb2430(qqbr)
          call qqb2431(qqbr);call qqb2432(qqbr);call qqb2433(qqbr);call qqb2434(qqbr);call qqb2435(qqbr);call qqb2436(qqbr);call
     & qqb2437(qqbr);call qqb2438(qqbr);call qqb2439(qqbr);call qqb2440(qqbr)
          call qqb2441(qqbr);call qqb2442(qqbr);call qqb2443(qqbr);call qqb2444(qqbr);call qqb2445(qqbr);call qqb2446(qqbr);call
     & qqb2447(qqbr);call qqb2448(qqbr);call qqb2449(qqbr);call qqb2450(qqbr)
          call qqb2451(qqbr);call qqb2452(qqbr);call qqb2453(qqbr);call qqb2454(qqbr);call qqb2455(qqbr);call qqb2456(qqbr);call
     & qqb2457(qqbr);call qqb2458(qqbr);call qqb2459(qqbr);call qqb2460(qqbr)
          call qqb2461(qqbr);call qqb2462(qqbr);call qqb2463(qqbr);call qqb2464(qqbr);call qqb2465(qqbr);call qqb2466(qqbr);call
     & qqb2467(qqbr);call qqb2468(qqbr);call qqb2469(qqbr);call qqb2470(qqbr)
          call qqb2471(qqbr);call qqb2472(qqbr);call qqb2473(qqbr);call qqb2474(qqbr);call qqb2475(qqbr);call qqb2476(qqbr);call
     & qqb2477(qqbr);call qqb2478(qqbr);call qqb2479(qqbr);call qqb2480(qqbr)
          call qqb2481(qqbr);call qqb2482(qqbr);call qqb2483(qqbr);call qqb2484(qqbr);call qqb2485(qqbr);call qqb2486(qqbr);call
     & qqb2487(qqbr);call qqb2488(qqbr);call qqb2489(qqbr);call qqb2490(qqbr)
          call qqb2491(qqbr);call qqb2492(qqbr);call qqb2493(qqbr);call qqb2494(qqbr);call qqb2495(qqbr);call qqb2496(qqbr);call
     & qqb2497(qqbr);call qqb2498(qqbr);call qqb2499(qqbr);call qqb2500(qqbr)
          call qqb2501(qqbr);call qqb2502(qqbr);call qqb2503(qqbr);call qqb2504(qqbr);call qqb2505(qqbr);call qqb2506(qqbr);call
     & qqb2507(qqbr);call qqb2508(qqbr);call qqb2509(qqbr);call qqb2510(qqbr)
          call qqb2511(qqbr);call qqb2512(qqbr);call qqb2513(qqbr);call qqb2514(qqbr);call qqb2515(qqbr);call qqb2516(qqbr);call
     & qqb2517(qqbr);call qqb2518(qqbr);call qqb2519(qqbr);call qqb2520(qqbr)
          call qqb2521(qqbr);call qqb2522(qqbr);call qqb2523(qqbr);call qqb2524(qqbr);call qqb2525(qqbr);call qqb2526(qqbr);call
     & qqb2527(qqbr);call qqb2528(qqbr);call qqb2529(qqbr);call qqb2530(qqbr)
          call qqb2531(qqbr);call qqb2532(qqbr);call qqb2533(qqbr);call qqb2534(qqbr);call qqb2535(qqbr);call qqb2536(qqbr);call
     & qqb2537(qqbr);call qqb2538(qqbr);call qqb2539(qqbr);call qqb2540(qqbr)
          call qqb2541(qqbr);call qqb2542(qqbr);call qqb2543(qqbr);call qqb2544(qqbr);call qqb2545(qqbr);call qqb2546(qqbr);call
     & qqb2547(qqbr);call qqb2548(qqbr);call qqb2549(qqbr);call qqb2550(qqbr)
          call qqb2551(qqbr);call qqb2552(qqbr);call qqb2553(qqbr);call qqb2554(qqbr);call qqb2555(qqbr);call qqb2556(qqbr);call
     & qqb2557(qqbr);call qqb2558(qqbr);call qqb2559(qqbr);call qqb2560(qqbr)
          call qqb2561(qqbr);call qqb2562(qqbr);call qqb2563(qqbr);call qqb2564(qqbr);call qqb2565(qqbr);call qqb2566(qqbr);call
     & qqb2567(qqbr);call qqb2568(qqbr);call qqb2569(qqbr);call qqb2570(qqbr)
          call qqb2571(qqbr);call qqb2572(qqbr);call qqb2573(qqbr);call qqb2574(qqbr);call qqb2575(qqbr);call qqb2576(qqbr);call
     & qqb2577(qqbr);call qqb2578(qqbr);call qqb2579(qqbr);call qqb2580(qqbr)
          call qqb2581(qqbr);call qqb2582(qqbr);call qqb2583(qqbr);call qqb2584(qqbr);call qqb2585(qqbr);call qqb2586(qqbr);call
     & qqb2587(qqbr);call qqb2588(qqbr);call qqb2589(qqbr);call qqb2590(qqbr)
          call qqb2591(qqbr);call qqb2592(qqbr);call qqb2593(qqbr);call qqb2594(qqbr);call qqb2595(qqbr);call qqb2596(qqbr);call
     & qqb2597(qqbr);call qqb2598(qqbr);call qqb2599(qqbr);call qqb2600(qqbr)
          call qqb2601(qqbr);call qqb2602(qqbr);call qqb2603(qqbr);call qqb2604(qqbr);call qqb2605(qqbr);call qqb2606(qqbr);call
     & qqb2607(qqbr);call qqb2608(qqbr);call qqb2609(qqbr);call qqb2610(qqbr)
          call qqb2611(qqbr);call qqb2612(qqbr);call qqb2613(qqbr);call qqb2614(qqbr);call qqb2615(qqbr);call qqb2616(qqbr);call
     & qqb2617(qqbr);call qqb2618(qqbr);call qqb2619(qqbr);call qqb2620(qqbr)
          call qqb2621(qqbr);call qqb2622(qqbr);call qqb2623(qqbr);call qqb2624(qqbr);call qqb2625(qqbr);call qqb2626(qqbr);call
     & qqb2627(qqbr);call qqb2628(qqbr);call qqb2629(qqbr);call qqb2630(qqbr)
          call qqb2631(qqbr);call qqb2632(qqbr);call qqb2633(qqbr);call qqb2634(qqbr);call qqb2635(qqbr);call qqb2636(qqbr);call
     & qqb2637(qqbr);call qqb2638(qqbr);call qqb2639(qqbr);call qqb2640(qqbr)
          call qqb2641(qqbr);call qqb2642(qqbr);call qqb2643(qqbr);call qqb2644(qqbr);call qqb2645(qqbr);call qqb2646(qqbr);call
     & qqb2647(qqbr);call qqb2648(qqbr);call qqb2649(qqbr);call qqb2650(qqbr)
          call qqb2651(qqbr);call qqb2652(qqbr);call qqb2653(qqbr);call qqb2654(qqbr);call qqb2655(qqbr);call qqb2656(qqbr);call
     & qqb2657(qqbr);call qqb2658(qqbr);call qqb2659(qqbr);call qqb2660(qqbr)
          call qqb2661(qqbr);call qqb2662(qqbr);call qqb2663(qqbr);call qqb2664(qqbr);call qqb2665(qqbr);call qqb2666(qqbr);call
     & qqb2667(qqbr);call qqb2668(qqbr);call qqb2669(qqbr);call qqb2670(qqbr)
          call qqb2671(qqbr);call qqb2672(qqbr);call qqb2673(qqbr);call qqb2674(qqbr);call qqb2675(qqbr);call qqb2676(qqbr);call
     & qqb2677(qqbr);call qqb2678(qqbr);call qqb2679(qqbr);call qqb2680(qqbr)
          call qqb2681(qqbr);call qqb2682(qqbr);call qqb2683(qqbr);call qqb2684(qqbr);call qqb2685(qqbr);call qqb2686(qqbr);call
     & qqb2687(qqbr);call qqb2688(qqbr);call qqb2689(qqbr);call qqb2690(qqbr)
          call qqb2691(qqbr);call qqb2692(qqbr);call qqb2693(qqbr);call qqb2694(qqbr);call qqb2695(qqbr);call qqb2696(qqbr);call
     & qqb2697(qqbr);call qqb2698(qqbr);call qqb2699(qqbr);call qqb2700(qqbr)
          call qqb2701(qqbr);call qqb2702(qqbr);call qqb2703(qqbr);call qqb2704(qqbr);call qqb2705(qqbr);call qqb2706(qqbr);call
     & qqb2707(qqbr);call qqb2708(qqbr);call qqb2709(qqbr);call qqb2710(qqbr)
          call qqb2711(qqbr);call qqb2712(qqbr);call qqb2713(qqbr);call qqb2714(qqbr);call qqb2715(qqbr);call qqb2716(qqbr);call
     & qqb2717(qqbr);call qqb2718(qqbr);call qqb2719(qqbr);call qqb2720(qqbr)
          call qqb2721(qqbr);call qqb2722(qqbr);call qqb2723(qqbr);call qqb2724(qqbr);call qqb2725(qqbr);call qqb2726(qqbr);call
     & qqb2727(qqbr);call qqb2728(qqbr);call qqb2729(qqbr);call qqb2730(qqbr)
          call qqb2731(qqbr);call qqb2732(qqbr);call qqb2733(qqbr);call qqb2734(qqbr);call qqb2735(qqbr);call qqb2736(qqbr);call
     & qqb2737(qqbr);call qqb2738(qqbr);call qqb2739(qqbr);call qqb2740(qqbr)
          call qqb2741(qqbr);call qqb2742(qqbr);call qqb2743(qqbr);call qqb2744(qqbr);call qqb2745(qqbr);call qqb2746(qqbr);call
     & qqb2747(qqbr);call qqb2748(qqbr);call qqb2749(qqbr);call qqb2750(qqbr)
          call qqb2751(qqbr);call qqb2752(qqbr);call qqb2753(qqbr);call qqb2754(qqbr);call qqb2755(qqbr);call qqb2756(qqbr);call
     & qqb2757(qqbr);call qqb2758(qqbr);call qqb2759(qqbr);call qqb2760(qqbr)
          call qqb2761(qqbr);call qqb2762(qqbr);call qqb2763(qqbr);call qqb2764(qqbr);call qqb2765(qqbr);call qqb2766(qqbr);call
     & qqb2767(qqbr);call qqb2768(qqbr);call qqb2769(qqbr);call qqb2770(qqbr)
          call qqb2771(qqbr);call qqb2772(qqbr);call qqb2773(qqbr);call qqb2774(qqbr);call qqb2775(qqbr);call qqb2776(qqbr);call
     & qqb2777(qqbr);call qqb2778(qqbr);call qqb2779(qqbr);call qqb2780(qqbr)
          call qqb2781(qqbr);call qqb2782(qqbr);call qqb2783(qqbr);call qqb2784(qqbr);call qqb2785(qqbr);call qqb2786(qqbr);call
     & qqb2787(qqbr);call qqb2788(qqbr);call qqb2789(qqbr);call qqb2790(qqbr)
          call qqb2791(qqbr);call qqb2792(qqbr);call qqb2793(qqbr);call qqb2794(qqbr);call qqb2795(qqbr);call qqb2796(qqbr);call
     & qqb2797(qqbr);call qqb2798(qqbr);call qqb2799(qqbr);call qqb2800(qqbr)
          call qqb2801(qqbr);call qqb2802(qqbr);call qqb2803(qqbr);call qqb2804(qqbr);call qqb2805(qqbr);call qqb2806(qqbr);call
     & qqb2807(qqbr);call qqb2808(qqbr);call qqb2809(qqbr);call qqb2810(qqbr)
          call qqb2811(qqbr);call qqb2812(qqbr);call qqb2813(qqbr);call qqb2814(qqbr);call qqb2815(qqbr);call qqb2816(qqbr);call
     & qqb2817(qqbr);call qqb2818(qqbr);call qqb2819(qqbr);call qqb2820(qqbr)
          call qqb2821(qqbr);call qqb2822(qqbr);call qqb2823(qqbr);call qqb2824(qqbr);call qqb2825(qqbr);call qqb2826(qqbr);call
     & qqb2827(qqbr);call qqb2828(qqbr);call qqb2829(qqbr);call qqb2830(qqbr)
          call qqb2831(qqbr);call qqb2832(qqbr);call qqb2833(qqbr);call qqb2834(qqbr);call qqb2835(qqbr);call qqb2836(qqbr);call
     & qqb2837(qqbr);call qqb2838(qqbr);call qqb2839(qqbr);call qqb2840(qqbr)
          call qqb2841(qqbr);call qqb2842(qqbr);call qqb2843(qqbr);call qqb2844(qqbr);call qqb2845(qqbr);call qqb2846(qqbr);call
     & qqb2847(qqbr);call qqb2848(qqbr);call qqb2849(qqbr);call qqb2850(qqbr)
          call qqb2851(qqbr);call qqb2852(qqbr);call qqb2853(qqbr);call qqb2854(qqbr);call qqb2855(qqbr);call qqb2856(qqbr);call
     & qqb2857(qqbr);call qqb2858(qqbr);call qqb2859(qqbr);call qqb2860(qqbr)
          call qqb2861(qqbr);call qqb2862(qqbr);call qqb2863(qqbr);call qqb2864(qqbr);call qqb2865(qqbr);call qqb2866(qqbr);call
     & qqb2867(qqbr);call qqb2868(qqbr);call qqb2869(qqbr);call qqb2870(qqbr)
          call qqb2871(qqbr);call qqb2872(qqbr);call qqb2873(qqbr);call qqb2874(qqbr);call qqb2875(qqbr);call qqb2876(qqbr);call
     & qqb2877(qqbr);call qqb2878(qqbr);call qqb2879(qqbr);call qqb2880(qqbr)
          call qqb2881(qqbr);call qqb2882(qqbr);call qqb2883(qqbr);call qqb2884(qqbr);call qqb2885(qqbr);call qqb2886(qqbr);call
     & qqb2887(qqbr);call qqb2888(qqbr);call qqb2889(qqbr);call qqb2890(qqbr)
          call qqb2891(qqbr);call qqb2892(qqbr);call qqb2893(qqbr);call qqb2894(qqbr);call qqb2895(qqbr);call qqb2896(qqbr);call
     & qqb2897(qqbr);call qqb2898(qqbr);call qqb2899(qqbr);call qqb2900(qqbr)
          call qqb2901(qqbr);call qqb2902(qqbr);call qqb2903(qqbr);call qqb2904(qqbr);call qqb2905(qqbr);call qqb2906(qqbr);call
     & qqb2907(qqbr);call qqb2908(qqbr);call qqb2909(qqbr);call qqb2910(qqbr)
          call qqb2911(qqbr);call qqb2912(qqbr);call qqb2913(qqbr);call qqb2914(qqbr);call qqb2915(qqbr);call qqb2916(qqbr);call
     & qqb2917(qqbr);call qqb2918(qqbr);call qqb2919(qqbr);call qqb2920(qqbr)
          call qqb2921(qqbr);call qqb2922(qqbr);call qqb2923(qqbr);call qqb2924(qqbr);call qqb2925(qqbr);call qqb2926(qqbr);call
     & qqb2927(qqbr);call qqb2928(qqbr);call qqb2929(qqbr);call qqb2930(qqbr)
          call qqb2931(qqbr);call qqb2932(qqbr);call qqb2933(qqbr);call qqb2934(qqbr);call qqb2935(qqbr);call qqb2936(qqbr);call
     & qqb2937(qqbr);call qqb2938(qqbr);call qqb2939(qqbr);call qqb2940(qqbr)
          call qqb2941(qqbr);call qqb2942(qqbr);call qqb2943(qqbr);call qqb2944(qqbr);call qqb2945(qqbr);call qqb2946(qqbr);call
     & qqb2947(qqbr);call qqb2948(qqbr);call qqb2949(qqbr);call qqb2950(qqbr)
          call qqb2951(qqbr);call qqb2952(qqbr);call qqb2953(qqbr);call qqb2954(qqbr);call qqb2955(qqbr);call qqb2956(qqbr);call
     & qqb2957(qqbr);call qqb2958(qqbr);call qqb2959(qqbr);call qqb2960(qqbr)
          call qqb2961(qqbr);call qqb2962(qqbr);call qqb2963(qqbr);call qqb2964(qqbr);call qqb2965(qqbr);call qqb2966(qqbr);call
     & qqb2967(qqbr);call qqb2968(qqbr);call qqb2969(qqbr);call qqb2970(qqbr)
          call qqb2971(qqbr);call qqb2972(qqbr);call qqb2973(qqbr);call qqb2974(qqbr);call qqb2975(qqbr);call qqb2976(qqbr);call
     & qqb2977(qqbr);call qqb2978(qqbr);call qqb2979(qqbr);call qqb2980(qqbr)
          call qqb2981(qqbr);call qqb2982(qqbr);call qqb2983(qqbr);call qqb2984(qqbr);call qqb2985(qqbr);call qqb2986(qqbr);call
     & qqb2987(qqbr);call qqb2988(qqbr);call qqb2989(qqbr);call qqb2990(qqbr)
          call qqb2991(qqbr);call qqb2992(qqbr);call qqb2993(qqbr);call qqb2994(qqbr);call qqb2995(qqbr);call qqb2996(qqbr);call
     & qqb2997(qqbr);call qqb2998(qqbr);call qqb2999(qqbr);call qqb3000(qqbr)
          call qqb3001(qqbr);call qqb3002(qqbr);call qqb3003(qqbr);call qqb3004(qqbr);call qqb3005(qqbr);call qqb3006(qqbr);call
     & qqb3007(qqbr);call qqb3008(qqbr);call qqb3009(qqbr);call qqb3010(qqbr)
          call qqb3011(qqbr);call qqb3012(qqbr);call qqb3013(qqbr);call qqb3014(qqbr);call qqb3015(qqbr);call qqb3016(qqbr);call
     & qqb3017(qqbr);call qqb3018(qqbr);call qqb3019(qqbr);call qqb3020(qqbr)
          call qqb3021(qqbr);call qqb3022(qqbr);call qqb3023(qqbr);call qqb3024(qqbr);call qqb3025(qqbr);call qqb3026(qqbr);call
     & qqb3027(qqbr);call qqb3028(qqbr);call qqb3029(qqbr);call qqb3030(qqbr)
          call qqb3031(qqbr);call qqb3032(qqbr);call qqb3033(qqbr);call qqb3034(qqbr);call qqb3035(qqbr);call qqb3036(qqbr);call
     & qqb3037(qqbr);call qqb3038(qqbr);call qqb3039(qqbr);call qqb3040(qqbr)
          call qqb3041(qqbr);call qqb3042(qqbr);call qqb3043(qqbr);call qqb3044(qqbr);call qqb3045(qqbr);call qqb3046(qqbr);call
     & qqb3047(qqbr);call qqb3048(qqbr);call qqb3049(qqbr);call qqb3050(qqbr)
          call qqb3051(qqbr)
      endif

      qqbAi0A1(0) = (0._dp)
      qqbAi0A1(1) = (0._dp)
      qqbAi0B1(0) = (0._dp)
      qqbAi0B1(1) = (0._dp)
      qqbAi0C1(0) = (0._dp)
      qqbAi0C1(1) = (0._dp)
      qqbAi0A2(0) = (0._dp)
      qqbAi0A2(1) = (0._dp)
      qqbAi0B2(0) = (0._dp)
      qqbAi0B2(1) = (0._dp)
      qqbAi0C2(0) = (0._dp)
      qqbAi0C2(1) = (0._dp)
      qqbAi0A3(0) = (0._dp)
      qqbAi0A3(1) = (0._dp)
      qqbAi0B3(0) = (0._dp)
      qqbAi0B3(1) = (0._dp)
      qqbAi0C3(0) = (0._dp)
      qqbAi0C3(1) = (0._dp)
      qqbAi0A4(0) = (0._dp)
      qqbAi0A4(1) = (0._dp)
      qqbAi0B4(0) = (0._dp)
      qqbAi0B4(1) = (0._dp)
      qqbAi0C4(0) = (0._dp)
      qqbAi0C4(1) = (0._dp)
      qqbAi0A5(0) = (0._dp)
      qqbAi0A5(1) = (0._dp)
      qqbAi0B5(0) = (0._dp)
      qqbAi0B5(1) = (0._dp)
      qqbAi0C5(0) = (0._dp)
      qqbAi0C5(1) = (0._dp)
      qqbAi0A6(0) = (0._dp)
      qqbAi0A6(1) = (0._dp)
      qqbAi0B6(0) = (0._dp)
      qqbAi0B6(1) = (0._dp)
      qqbAi0C6(0) = (0._dp)
      qqbAi0C6(1) = (0._dp)
      qqbAi0A7(0) = qqbr(611273)
      qqbAi0A7(1) = (0._dp)
      qqbAi0B7(0) = (0._dp)
      qqbAi0B7(1) = (0._dp)
      qqbAi0C7(0) = (0._dp)
      qqbAi0C7(1) = (0._dp)
      qqbAi0A8(0) = (0._dp)
      qqbAi0A8(1) = (0._dp)
      qqbAi0B8(0) = qqbr(611274)
      qqbAi0B8(1) = (0._dp)
      qqbAi0C8(0) = (0._dp)
      qqbAi0C8(1) = (0._dp)
      qqbAi0A9(0) = (0._dp)
      qqbAi0A9(1) = (0._dp)
      qqbAi0B9(0) = qqbr(611275)
      qqbAi0B9(1) = (0._dp)
      qqbAi0C9(0) = (0._dp)
      qqbAi0C9(1) = (0._dp)
      qqbAi0A10(0) = qqbr(611091)
      qqbAi0A10(1) = (0._dp)
      qqbAi0B10(0) = (0._dp)
      qqbAi0B10(1) = (0._dp)
      qqbAi0C10(0) = (0._dp)
      qqbAi0C10(1) = (0._dp)

      qqbAi1A1(0) = CF*qqbr(611126)
      qqbAi1A1(1) = CF*pi*qqbr(611099)
      qqbAi1B1(0) = CF*qqbr(611127)
      qqbAi1B1(1) = CF*pi*qqbr(611100)
      qqbAi1C1(0) = (0._dp)
      qqbAi1C1(1) = (0._dp)
      qqbAi1A2(0) = CF*qqbr(611124)
      qqbAi1A2(1) = CF*pi*qqbr(611110)
      qqbAi1B2(0) = CF*qqbr(611123)
      qqbAi1B2(1) = CF*pi*qqbr(611101)
      qqbAi1C2(0) = (0._dp)
      qqbAi1C2(1) = (0._dp)
      qqbAi1A3(0) = CF*qqbr(611121)
      qqbAi1A3(1) = CF*pi*qqbr(611096)
      qqbAi1B3(0) = CF*qqbr(611122)
      qqbAi1B3(1) = CF*pi*qqbr(611109)
      qqbAi1C3(0) = (0._dp)
      qqbAi1C3(1) = (0._dp)
      qqbAi1A4(0) = CF*qqbr(611125)
      qqbAi1A4(1) = CF*pi*qqbr(611095)
      qqbAi1B4(0) = CF*qqbr(611128)
      qqbAi1B4(1) = CF*pi*qqbr(611104)
      qqbAi1C4(0) = (0._dp)
      qqbAi1C4(1) = (0._dp)
      qqbAi1A5(0) = CF*qqbr(611117)
      qqbAi1A5(1) = CF*pi*qqbr(611097)
      qqbAi1B5(0) = CF*qqbr(611119)
      qqbAi1B5(1) = CF*pi*qqbr(611107)
      qqbAi1C5(0) = (0._dp)
      qqbAi1C5(1) = (0._dp)
      qqbAi1A6(0) = CF*qqbr(611120)
      qqbAi1A6(1) = CF*pi*qqbr(611108)
      qqbAi1B6(0) = CF*qqbr(611115)
      qqbAi1B6(1) = CF*pi*qqbr(611102)
      qqbAi1C6(0) = (0._dp)
      qqbAi1C6(1) = (0._dp)
      qqbAi1A7(0) = CF*qqbr(611133)
      qqbAi1A7(1) = CF*pi*qqbr(611111)
      qqbAi1B7(0) = CF*qqbr(611116)
      qqbAi1B7(1) = CF*pi*qqbr(611103)
      qqbAi1C7(0) = (0._dp)
      qqbAi1C7(1) = (0._dp)
      qqbAi1A8(0) = CF*qqbr(611118)
      qqbAi1A8(1) = CF*pi*qqbr(611098)
      qqbAi1B8(0) = CF*qqbr(611134)
      qqbAi1B8(1) = CF*pi*qqbr(611112)
      qqbAi1C8(0) = (0._dp)
      qqbAi1C8(1) = (0._dp)
      qqbAi1A9(0) = CF*qqbr(611113)
      qqbAi1A9(1) = CF*pi*qqbr(611093)
      qqbAi1B9(0) = CF*qqbr(611130)
      qqbAi1B9(1) = CF*pi*qqbr(611106)
      qqbAi1C9(0) = (0._dp)
      qqbAi1C9(1) = (0._dp)
      qqbAi1A10(0) = CF*qqbr(611129)
      qqbAi1A10(1) = CF*pi*qqbr(611105)
      qqbAi1B10(0) = CF*qqbr(611114)
      qqbAi1B10(1) = CF*pi*qqbr(611094)
      qqbAi1C10(0) = (0._dp)
      qqbAi1C10(1) = (0._dp)

      qqbAi2A1(0) = CF*(nflav)*qqbr(611166)+CF*CF*qqbr(611215)+CF*((1._dp)/N)*qqbr(611246)
      qqbAi2A1(1) = pi*(CF*(nflav)*qqbr(611138)+CF*CF*qqbr(611174)+CF*((1._dp)/N)*qqbr(611216))
      qqbAi2B1(0) = CF*(nflav)*qqbr(611167)+CF*CF*qqbr(611218)+CF*((1._dp)/N)*qqbr(611247)
      qqbAi2B1(1) = pi*(CF*(nflav)*qqbr(611139)+CF*CF*qqbr(611175)+CF*((1._dp)/N)*qqbr(611217))
      qqbAi2C1(0) = CF*qqbr(611256)
      qqbAi2C1(1) = CF*pi*qqbr(611238)
      qqbAi2A2(0) = CF*(nflav)*qqbr(611164)+CF*CF*qqbr(611201)+CF*((1._dp)/N)*qqbr(611268)
      qqbAi2A2(1) = pi*(CF*(nflav)*qqbr(611148)+CF*CF*qqbr(611179)+CF*((1._dp)/N)*qqbr(611230))
      qqbAi2B2(0) = CF*(nflav)*qqbr(611163)+CF*CF*qqbr(611202)+CF*((1._dp)/N)*qqbr(611267)
      qqbAi2B2(1) = pi*(CF*(nflav)*qqbr(611144)+CF*CF*qqbr(611180)+CF*((1._dp)/N)*qqbr(611224))
      qqbAi2C2(0) = CF*qqbr(611254)
      qqbAi2C2(1) = CF*pi*qqbr(611236)
      qqbAi2A3(0) = CF*(nflav)*qqbr(611161)+CF*CF*qqbr(611203)+CF*((1._dp)/N)*qqbr(611265)
      qqbAi2A3(1) = pi*(CF*(nflav)*qqbr(611143)+CF*CF*qqbr(611181)+CF*((1._dp)/N)*qqbr(611223))
      qqbAi2B3(0) = CF*(nflav)*qqbr(611162)+CF*CF*qqbr(611204)+CF*((1._dp)/N)*qqbr(611266)
      qqbAi2B3(1) = pi*(CF*(nflav)*qqbr(611147)+CF*CF*qqbr(611182)+CF*((1._dp)/N)*qqbr(611229))
      qqbAi2C3(0) = CF*qqbr(611253)
      qqbAi2C3(1) = CF*pi*qqbr(611235)
      qqbAi2A4(0) = CF*(nflav)*qqbr(611165)+CF*CF*qqbr(611214)+CF*((1._dp)/N)*qqbr(611245)
      qqbAi2A4(1) = pi*(CF*(nflav)*qqbr(611135)+CF*CF*qqbr(611173)+CF*((1._dp)/N)*qqbr(611213))
      qqbAi2B4(0) = CF*(nflav)*qqbr(611168)+CF*CF*qqbr(611219)+CF*((1._dp)/N)*qqbr(611248)
      qqbAi2B4(1) = pi*(CF*(nflav)*qqbr(611142)+CF*CF*qqbr(611176)+CF*((1._dp)/N)*qqbr(611220))
      qqbAi2C4(0) = CF*qqbr(611255)
      qqbAi2C4(1) = CF*pi*qqbr(611237)
      qqbAi2A5(0) = CF*(nflav)*qqbr(611171)+CF*CF*qqbr(611207)+CF*((1._dp)/N)*qqbr(611272)
      qqbAi2A5(1) = pi*(CF*(nflav)*qqbr(611136)+CF*CF*qqbr(611184)+CF*((1._dp)/N)*qqbr(611228))
      qqbAi2B5(0) = CF*(nflav)*qqbr(611157)+CF*CF*qqbr(611189)+CF*((1._dp)/N)*qqbr(611257)
      qqbAi2B5(1) = pi*(CF*(nflav)*qqbr(611145)+CF*CF*qqbr(611185)+CF*((1._dp)/N)*qqbr(611221))
      qqbAi2C5(0) = CF*qqbr(611250)
      qqbAi2C5(1) = CF*pi*qqbr(611232)
      qqbAi2A6(0) = CF*(nflav)*qqbr(611158)+CF*CF*qqbr(611190)+CF*((1._dp)/N)*qqbr(611258)
      qqbAi2A6(1) = pi*(CF*(nflav)*qqbr(611146)+CF*CF*qqbr(611186)+CF*((1._dp)/N)*qqbr(611222))
      qqbAi2B6(0) = CF*(nflav)*qqbr(611169)+CF*CF*qqbr(611209)+CF*((1._dp)/N)*qqbr(611270)
      qqbAi2B6(1) = pi*(CF*(nflav)*qqbr(611140)+CF*CF*qqbr(611187)+CF*((1._dp)/N)*qqbr(611226))
      qqbAi2C6(0) = CF*qqbr(611252)
      qqbAi2C6(1) = CF*pi*qqbr(611234)
      qqbAi2A7(0) = CF*(nflav)*qqbr(611159)+CF*CF*qqbr(611239)+CF*((1._dp)/N)*qqbr(611261)
      qqbAi2A7(1) = pi*(CF*(nflav)*qqbr(611153)+CF*((1._dp)/N)*qqbr(611195)+CF*CF*qqbr(611205))
      qqbAi2B7(0) = CF*(nflav)*qqbr(611170)+CF*CF*qqbr(611206)+CF*((1._dp)/N)*qqbr(611271)
      qqbAi2B7(1) = pi*(CF*(nflav)*qqbr(611141)+CF*CF*qqbr(611183)+CF*((1._dp)/N)*qqbr(611227))
      qqbAi2C7(0) = CF*qqbr(611251)
      qqbAi2C7(1) = CF*pi*qqbr(611233)
      qqbAi2A8(0) = CF*(nflav)*qqbr(611172)+CF*CF*qqbr(611210)+CF*((1._dp)/N)*qqbr(611269)
      qqbAi2A8(1) = pi*(CF*(nflav)*qqbr(611137)+CF*CF*qqbr(611188)+CF*((1._dp)/N)*qqbr(611225))
      qqbAi2B8(0) = CF*(nflav)*qqbr(611160)+CF*CF*qqbr(611241)+CF*((1._dp)/N)*qqbr(611263)
      qqbAi2B8(1) = pi*(CF*(nflav)*qqbr(611154)+CF*((1._dp)/N)*qqbr(611198)+CF*CF*qqbr(611208))
      qqbAi2C8(0) = CF*qqbr(611249)
      qqbAi2C8(1) = CF*pi*qqbr(611231)
      qqbAi2A9(0) = CF*(nflav)*qqbr(611149)+CF*CF*qqbr(611191)+CF*((1._dp)/N)*qqbr(611260)
      qqbAi2A9(1) = pi*(CF*(nflav)*qqbr(611131)+CF*CF*qqbr(611177)+CF*((1._dp)/N)*qqbr(611211))
      qqbAi2B9(0) = CF*(nflav)*qqbr(611156)+CF*CF*qqbr(611242)+CF*((1._dp)/N)*qqbr(611264)
      qqbAi2B9(1) = pi*(CF*(nflav)*qqbr(611152)+CF*((1._dp)/N)*qqbr(611197)+CF*CF*qqbr(611200))
      qqbAi2C9(0) = CF*qqbr(611243)
      qqbAi2C9(1) = CF*pi*qqbr(611193)
      qqbAi2A10(0) = CF*(nflav)*qqbr(611155)+CF*CF*qqbr(611240)+CF*((1._dp)/N)*qqbr(611262)
      qqbAi2A10(1) = pi*(CF*(nflav)*qqbr(611151)+CF*((1._dp)/N)*qqbr(611196)+CF*CF*qqbr(611199))
      qqbAi2B10(0) = CF*(nflav)*qqbr(611150)+CF*CF*qqbr(611192)+CF*((1._dp)/N)*qqbr(611259)
      qqbAi2B10(1) = pi*(CF*(nflav)*qqbr(611132)+CF*CF*qqbr(611178)+CF*((1._dp)/N)*qqbr(611212))
      qqbAi2C10(0) = CF*qqbr(611244)
      qqbAi2C10(1) = CF*pi*qqbr(611194)


      qqbAj(0,A,1) = cmplx(qqbAi0A1(0),qqbAi0A1(1),kind=dp)
      qqbAj(0,A,2) = cmplx(qqbAi0A2(0),qqbAi0A2(1),kind=dp)
      qqbAj(0,A,3) = cmplx(qqbAi0A3(0),qqbAi0A3(1),kind=dp)
      qqbAj(0,A,4) = cmplx(qqbAi0A4(0),qqbAi0A4(1),kind=dp)
      qqbAj(0,A,5) = cmplx(qqbAi0A5(0),qqbAi0A5(1),kind=dp)
      qqbAj(0,A,6) = cmplx(qqbAi0A6(0),qqbAi0A6(1),kind=dp)
      qqbAj(0,A,7) = cmplx(qqbAi0A7(0),qqbAi0A7(1),kind=dp)
      qqbAj(0,A,8) = cmplx(qqbAi0A8(0),qqbAi0A8(1),kind=dp)
      qqbAj(0,A,9) = cmplx(qqbAi0A9(0),qqbAi0A9(1),kind=dp)
      qqbAj(0,A,10) = cmplx(qqbAi0A10(0),qqbAi0A10(1),kind=dp)

      qqbAj(0,B,1) = cmplx(qqbAi0B1(0),qqbAi0B1(1),kind=dp)
      qqbAj(0,B,2) = cmplx(qqbAi0B2(0),qqbAi0B2(1),kind=dp)
      qqbAj(0,B,3) = cmplx(qqbAi0B3(0),qqbAi0B3(1),kind=dp)
      qqbAj(0,B,4) = cmplx(qqbAi0B4(0),qqbAi0B4(1),kind=dp)
      qqbAj(0,B,5) = cmplx(qqbAi0B5(0),qqbAi0B5(1),kind=dp)
      qqbAj(0,B,6) = cmplx(qqbAi0B6(0),qqbAi0B6(1),kind=dp)
      qqbAj(0,B,7) = cmplx(qqbAi0B7(0),qqbAi0B7(1),kind=dp)
      qqbAj(0,B,8) = cmplx(qqbAi0B8(0),qqbAi0B8(1),kind=dp)
      qqbAj(0,B,9) = cmplx(qqbAi0B9(0),qqbAi0B9(1),kind=dp)
      qqbAj(0,B,10) = cmplx(qqbAi0B10(0),qqbAi0B10(1),kind=dp)

      qqbAj(0,C,1) = cmplx(qqbAi0C1(0),qqbAi0C1(1),kind=dp)
      qqbAj(0,C,2) = cmplx(qqbAi0C2(0),qqbAi0C2(1),kind=dp)
      qqbAj(0,C,3) = cmplx(qqbAi0C3(0),qqbAi0C3(1),kind=dp)
      qqbAj(0,C,4) = cmplx(qqbAi0C4(0),qqbAi0C4(1),kind=dp)
      qqbAj(0,C,5) = cmplx(qqbAi0C5(0),qqbAi0C5(1),kind=dp)
      qqbAj(0,C,6) = cmplx(qqbAi0C6(0),qqbAi0C6(1),kind=dp)
      qqbAj(0,C,7) = cmplx(qqbAi0C7(0),qqbAi0C7(1),kind=dp)
      qqbAj(0,C,8) = cmplx(qqbAi0C8(0),qqbAi0C8(1),kind=dp)
      qqbAj(0,C,9) = cmplx(qqbAi0C9(0),qqbAi0C9(1),kind=dp)
      qqbAj(0,C,10) = cmplx(qqbAi0C10(0),qqbAi0C10(1),kind=dp)

      qqbAj(1,A,1) = cmplx(qqbAi1A1(0),qqbAi1A1(1),kind=dp)
      qqbAj(1,A,2) = cmplx(qqbAi1A2(0),qqbAi1A2(1),kind=dp)
      qqbAj(1,A,3) = cmplx(qqbAi1A3(0),qqbAi1A3(1),kind=dp)
      qqbAj(1,A,4) = cmplx(qqbAi1A4(0),qqbAi1A4(1),kind=dp)
      qqbAj(1,A,5) = cmplx(qqbAi1A5(0),qqbAi1A5(1),kind=dp)
      qqbAj(1,A,6) = cmplx(qqbAi1A6(0),qqbAi1A6(1),kind=dp)
      qqbAj(1,A,7) = cmplx(qqbAi1A7(0),qqbAi1A7(1),kind=dp)
      qqbAj(1,A,8) = cmplx(qqbAi1A8(0),qqbAi1A8(1),kind=dp)
      qqbAj(1,A,9) = cmplx(qqbAi1A9(0),qqbAi1A9(1),kind=dp)
      qqbAj(1,A,10) = cmplx(qqbAi1A10(0),qqbAi1A10(1),kind=dp)

      qqbAj(1,B,1) = cmplx(qqbAi1B1(0),qqbAi1B1(1),kind=dp)
      qqbAj(1,B,2) = cmplx(qqbAi1B2(0),qqbAi1B2(1),kind=dp)
      qqbAj(1,B,3) = cmplx(qqbAi1B3(0),qqbAi1B3(1),kind=dp)
      qqbAj(1,B,4) = cmplx(qqbAi1B4(0),qqbAi1B4(1),kind=dp)
      qqbAj(1,B,5) = cmplx(qqbAi1B5(0),qqbAi1B5(1),kind=dp)
      qqbAj(1,B,6) = cmplx(qqbAi1B6(0),qqbAi1B6(1),kind=dp)
      qqbAj(1,B,7) = cmplx(qqbAi1B7(0),qqbAi1B7(1),kind=dp)
      qqbAj(1,B,8) = cmplx(qqbAi1B8(0),qqbAi1B8(1),kind=dp)
      qqbAj(1,B,9) = cmplx(qqbAi1B9(0),qqbAi1B9(1),kind=dp)
      qqbAj(1,B,10) = cmplx(qqbAi1B10(0),qqbAi1B10(1),kind=dp)

      qqbAj(1,C,1) = cmplx(qqbAi1C1(0),qqbAi1C1(1),kind=dp)
      qqbAj(1,C,2) = cmplx(qqbAi1C2(0),qqbAi1C2(1),kind=dp)
      qqbAj(1,C,3) = cmplx(qqbAi1C3(0),qqbAi1C3(1),kind=dp)
      qqbAj(1,C,4) = cmplx(qqbAi1C4(0),qqbAi1C4(1),kind=dp)
      qqbAj(1,C,5) = cmplx(qqbAi1C5(0),qqbAi1C5(1),kind=dp)
      qqbAj(1,C,6) = cmplx(qqbAi1C6(0),qqbAi1C6(1),kind=dp)
      qqbAj(1,C,7) = cmplx(qqbAi1C7(0),qqbAi1C7(1),kind=dp)
      qqbAj(1,C,8) = cmplx(qqbAi1C8(0),qqbAi1C8(1),kind=dp)
      qqbAj(1,C,9) = cmplx(qqbAi1C9(0),qqbAi1C9(1),kind=dp)
      qqbAj(1,C,10) = cmplx(qqbAi1C10(0),qqbAi1C10(1),kind=dp)

      qqbAj(2,A,1) = cmplx(qqbAi2A1(0),qqbAi2A1(1),kind=dp)
      qqbAj(2,A,2) = cmplx(qqbAi2A2(0),qqbAi2A2(1),kind=dp)
      qqbAj(2,A,3) = cmplx(qqbAi2A3(0),qqbAi2A3(1),kind=dp)
      qqbAj(2,A,4) = cmplx(qqbAi2A4(0),qqbAi2A4(1),kind=dp)
      qqbAj(2,A,5) = cmplx(qqbAi2A5(0),qqbAi2A5(1),kind=dp)
      qqbAj(2,A,6) = cmplx(qqbAi2A6(0),qqbAi2A6(1),kind=dp)
      qqbAj(2,A,7) = cmplx(qqbAi2A7(0),qqbAi2A7(1),kind=dp)
      qqbAj(2,A,8) = cmplx(qqbAi2A8(0),qqbAi2A8(1),kind=dp)
      qqbAj(2,A,9) = cmplx(qqbAi2A9(0),qqbAi2A9(1),kind=dp)
      qqbAj(2,A,10) = cmplx(qqbAi2A10(0),qqbAi2A10(1),kind=dp)

      qqbAj(2,B,1) = cmplx(qqbAi2B1(0),qqbAi2B1(1),kind=dp)
      qqbAj(2,B,2) = cmplx(qqbAi2B2(0),qqbAi2B2(1),kind=dp)
      qqbAj(2,B,3) = cmplx(qqbAi2B3(0),qqbAi2B3(1),kind=dp)
      qqbAj(2,B,4) = cmplx(qqbAi2B4(0),qqbAi2B4(1),kind=dp)
      qqbAj(2,B,5) = cmplx(qqbAi2B5(0),qqbAi2B5(1),kind=dp)
      qqbAj(2,B,6) = cmplx(qqbAi2B6(0),qqbAi2B6(1),kind=dp)
      qqbAj(2,B,7) = cmplx(qqbAi2B7(0),qqbAi2B7(1),kind=dp)
      qqbAj(2,B,8) = cmplx(qqbAi2B8(0),qqbAi2B8(1),kind=dp)
      qqbAj(2,B,9) = cmplx(qqbAi2B9(0),qqbAi2B9(1),kind=dp)
      qqbAj(2,B,10) = cmplx(qqbAi2B10(0),qqbAi2B10(1),kind=dp)

      qqbAj(2,C,1) = cmplx(qqbAi2C1(0),qqbAi2C1(1),kind=dp)
      qqbAj(2,C,2) = cmplx(qqbAi2C2(0),qqbAi2C2(1),kind=dp)
      qqbAj(2,C,3) = cmplx(qqbAi2C3(0),qqbAi2C3(1),kind=dp)
      qqbAj(2,C,4) = cmplx(qqbAi2C4(0),qqbAi2C4(1),kind=dp)
      qqbAj(2,C,5) = cmplx(qqbAi2C5(0),qqbAi2C5(1),kind=dp)
      qqbAj(2,C,6) = cmplx(qqbAi2C6(0),qqbAi2C6(1),kind=dp)
      qqbAj(2,C,7) = cmplx(qqbAi2C7(0),qqbAi2C7(1),kind=dp)
      qqbAj(2,C,8) = cmplx(qqbAi2C8(0),qqbAi2C8(1),kind=dp)
      qqbAj(2,C,9) = cmplx(qqbAi2C9(0),qqbAi2C9(1),kind=dp)
      qqbAj(2,C,10) = cmplx(qqbAi2C10(0),qqbAi2C10(1),kind=dp)

      return
      end

      end module
