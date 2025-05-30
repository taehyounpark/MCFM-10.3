!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine vertices_tt2(mtsq,ep,facuLl,facuRl,
     &  vert25x1,vert25x2,vert25x3,vert25x4,vert25x5,vert25x6,
     &  vert25x7,vert25x8,vert25x9,vert25x10,vert25x11,
     &  vert25x12)
      use loopI2_generic
      use loopI3_generic
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'poles.f'
      include 'scale.f'
      include 'masses.f'
      include 'decl_kininv.f'
      include 'scalarselect.f'

      real(dp):: mtsq
      integer:: ep
      complex(dp):: eploopI3,ep2loopI3
      complex(dp):: facuLl,facuRl,facuRLdiff
      complex(dp):: vert25x1,vert25x2,vert25x3,vert25x4
      complex(dp):: vert25x5,vert25x6,vert25x7,vert25x8,vert25x9
      complex(dp):: vert25x10, vert25x11,vert25x12
      real(dp):: p5Dp34,den5x34


      mtsq=mt**2
      p5Dp34=(s345-s34-mtsq)/2d0
      den5x34=1d0/(s34*mtsq-p5Dp34**2)
      facuRLdiff=facuRl-facuLl

      vert25x1= + loopI2(s34,mtsq,mtsq,musq,ep) * ( mt*s34*den5x34*
     &    facuRLdiff + 3.D0*mt*s34*p5Dp34**2*den5x34**2*facuRLdiff + 3.D
     &    0/2.D0*mt*s34**2*p5Dp34*den5x34**2*facuRLdiff - facuLl*mt*s34
     &    *den5x34 + 3.D0*facuLl*mt*s34*p5Dp34**2*den5x34**2 + 3.D0/2.D0
     &    *facuLl*mt*s34**2*p5Dp34*den5x34**2 )
      vert25x1 = vert25x1 + loopI2(mtsq,zip,mtsq,musq,ep) * ( mt*p5Dp34*
     &    den5x34*facuRLdiff + 2.D0*mt*p5Dp34**3*den5x34**2*facuRLdiff
     &     + mt*s34*p5Dp34**2*den5x34**2*facuRLdiff + mt**3*p5Dp34*
     &    den5x34*facuRLdiff*s345**(-1) + mt**3*s34*den5x34*facuRLdiff*
     &    s345**(-1) + mt**3*s34*p5Dp34*den5x34**2*facuRLdiff + 1.D0/2.D
     &    0*mt**3*s34**2*den5x34**2*facuRLdiff - facuLl*mt*p5Dp34*
     &    den5x34 + 2.D0*facuLl*mt*p5Dp34**3*den5x34**2 + facuLl*mt*s34
     &    *p5Dp34**2*den5x34**2 + facuLl*mt**3*p5Dp34*den5x34*
     &    s345**(-1) + facuLl*mt**3*s34*den5x34*s345**(-1) + facuLl*
     &    mt**3*s34*p5Dp34*den5x34**2 + 1.D0/2.D0*facuLl*mt**3*s34**2*
     &    den5x34**2 )
      vert25x1 = vert25x1 + loopI2(s345,zip,mtsq,musq,ep) * (  - mt*
     &    p5Dp34*den5x34*facuRLdiff - 2.D0*mt*p5Dp34**3*den5x34**2*
     &    facuRLdiff - mt*s34*den5x34*facuRLdiff - 4.D0*mt*s34*
     &    p5Dp34**2*den5x34**2*facuRLdiff - 3.D0/2.D0*mt*s34**2*p5Dp34*
     &    den5x34**2*facuRLdiff - mt**3*p5Dp34*den5x34*facuRLdiff*
     &    s345**(-1) - mt**3*s34*den5x34*facuRLdiff*s345**(-1) - mt**3*
     &    s34*p5Dp34*den5x34**2*facuRLdiff - 1.D0/2.D0*mt**3*s34**2*
     &    den5x34**2*facuRLdiff + facuLl*mt*p5Dp34*den5x34 - 2.D0*
     &    facuLl*mt*p5Dp34**3*den5x34**2 + facuLl*mt*s34*den5x34 - 4.D0
     &    *facuLl*mt*s34*p5Dp34**2*den5x34**2 - 3.D0/2.D0*facuLl*mt*
     &    s34**2*p5Dp34*den5x34**2 - facuLl*mt**3*p5Dp34*den5x34*
     &    s345**(-1) - facuLl*mt**3*s34*den5x34*s345**(-1) - facuLl*
     &    mt**3*s34*p5Dp34*den5x34**2 - 1.D0/2.D0*facuLl*mt**3*s34**2*
     &    den5x34**2 )
      vert25x1 = vert25x1 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * (  - 4.D0*mt*p5Dp34**2*den5x34*facuRLdiff - 4.D0*mt*p5Dp34**4
     &    *den5x34**2*facuRLdiff - 2.D0*mt*s34*p5Dp34*den5x34*
     &    facuRLdiff - 4.D0*mt*s34*p5Dp34**3*den5x34**2*facuRLdiff - mt
     &    *s34**2*p5Dp34**2*den5x34**2*facuRLdiff - 2.D0*mt**3*s34*
     &    p5Dp34**2*den5x34**2*facuRLdiff - 2.D0*mt**3*s34**2*p5Dp34*
     &    den5x34**2*facuRLdiff - 1.D0/2.D0*mt**3*s34**3*den5x34**2*
     &    facuRLdiff + 4.D0*facuLl*mt - 4.D0*facuLl*mt*p5Dp34**4*
     &    den5x34**2 - 4.D0*facuLl*mt*s34*p5Dp34**3*den5x34**2 - facuLl
     &    *mt*s34**2*p5Dp34**2*den5x34**2 - 2.D0*facuLl*mt**3*s34*
     &    p5Dp34**2*den5x34**2 - 2.D0*facuLl*mt**3*s34**2*p5Dp34*
     &    den5x34**2 - 1.D0/2.D0*facuLl*mt**3*s34**3*den5x34**2 )
      vert25x1 = vert25x1 + fp(ep)*mt*p5Dp34*den5x34*facuRLdiff + fp(ep
     &    )*mt*s34*den5x34*facuRLdiff - fp(ep)*mt**3*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) - fp(ep)*mt**3*s34*den5x34*facuRLdiff*
     &    s345**(-1) + fp(ep)*facuLl*mt*p5Dp34*den5x34 + fp(ep)*facuLl*
     &    mt*s34*den5x34 - fp(ep)*facuLl*mt**3*p5Dp34*den5x34*
     &    s345**(-1) - fp(ep)*facuLl*mt**3*s34*den5x34*s345**(-1) - 2.D0
     &    *eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*s34*
     &    p5Dp34**2*den5x34**2*facuRLdiff - 2.D0*eploopI3(mtsq,s34,s345,
     &    zip,mtsq,mtsq,musq,ep)*mt**3*s34**2*p5Dp34*den5x34**2*
     &    facuRLdiff - 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     &    musq,ep)*mt**3*s34**3*den5x34**2*facuRLdiff - 2.D0*eploopI3(
     &    mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**3*s34*
     &    p5Dp34**2*den5x34**2 - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*facuLl*mt**3*s34**2*p5Dp34*den5x34**2 - 1.D0/2.D
     &    0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**3*
     &    s34**3*den5x34**2
      vert25x1 = vert25x1 - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     & musq,ep)*mt**3*s34*p5Dp34**2*den5x34**2*facuRLdiff - 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*s34**2*
     &    p5Dp34*den5x34**2*facuRLdiff - 1.D0/2.D0*ep2loopI3(mtsq,s34,
     &    s345,zip,mtsq,mtsq,musq,ep)*mt**3*s34**3*den5x34**2*
     &    facuRLdiff - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep
     &    )*facuLl*mt**3*s34*p5Dp34**2*den5x34**2 - 2.D0*ep2loopI3(mtsq,
     &    s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**3*s34**2*p5Dp34*
     &    den5x34**2 - 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     &    musq,ep)*facuLl*mt**3*s34**3*den5x34**2

      vert25x2= + loopI2(s34,mtsq,mtsq,musq,ep) * ( mt*p5Dp34*den5x34*
     &    facuRLdiff - mt*p5Dp34**3*den5x34**2*facuRLdiff - 1.D0/2.D0*
     &    mt*s34*p5Dp34**2*den5x34**2*facuRLdiff - 2.D0*mt**3*s34*
     &    p5Dp34*den5x34**2*facuRLdiff - mt**3*s34**2*den5x34**2*
     &    facuRLdiff + 3.D0*facuLl*mt*p5Dp34*den5x34 - facuLl*mt*
     &    p5Dp34**3*den5x34**2 - 1.D0/2.D0*facuLl*mt*s34*p5Dp34**2*
     &    den5x34**2 - 2.D0*facuLl*mt**3*s34*p5Dp34*den5x34**2 - facuLl
     &    *mt**3*s34**2*den5x34**2 )
      vert25x2 = vert25x2 + loopI2(mtsq,zip,mtsq,musq,ep) * ( mt**3*
     &    den5x34*facuRLdiff - mt**3*p5Dp34*den5x34*facuRLdiff*
     &    s345**(-1) - 3.D0*mt**3*p5Dp34**2*den5x34**2*facuRLdiff - 3.D0
     &    /2.D0*mt**3*s34*p5Dp34*den5x34**2*facuRLdiff - mt**5*den5x34*
     &    facuRLdiff*s345**(-1) + 3.D0*facuLl*mt**3*den5x34 - facuLl*
     &    mt**3*p5Dp34*den5x34*s345**(-1) - 3.D0*facuLl*mt**3*p5Dp34**2
     &    *den5x34**2 - 3.D0/2.D0*facuLl*mt**3*s34*p5Dp34*den5x34**2 -
     &    facuLl*mt**5*den5x34*s345**(-1) )
      vert25x2 = vert25x2 + loopI2(s345,zip,mtsq,musq,ep) * (  - mt*
     &    p5Dp34*den5x34*facuRLdiff + mt*p5Dp34**3*den5x34**2*
     &    facuRLdiff + 1.D0/2.D0*mt*s34*p5Dp34**2*den5x34**2*facuRLdiff
     &     - mt**3*den5x34*facuRLdiff + mt**3*p5Dp34*den5x34*facuRLdiff
     &    *s345**(-1) + 3.D0*mt**3*p5Dp34**2*den5x34**2*facuRLdiff + 7.D
     &    0/2.D0*mt**3*s34*p5Dp34*den5x34**2*facuRLdiff + mt**3*s34**2*
     &    den5x34**2*facuRLdiff + mt**5*den5x34*facuRLdiff*s345**(-1)
     &     - 3.D0*facuLl*mt*p5Dp34*den5x34 + facuLl*mt*p5Dp34**3*
     &    den5x34**2 + 1.D0/2.D0*facuLl*mt*s34*p5Dp34**2*den5x34**2 - 3.
     &    D0*facuLl*mt**3*den5x34 + facuLl*mt**3*p5Dp34*den5x34*
     &    s345**(-1) + 3.D0*facuLl*mt**3*p5Dp34**2*den5x34**2 + 7.D0/2.D
     &    0*facuLl*mt**3*s34*p5Dp34*den5x34**2 + facuLl*mt**3*s34**2*
     &    den5x34**2 + facuLl*mt**5*den5x34*s345**(-1) )
      vert25x2 = vert25x2 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * ( 4.D0*mt*facuRLdiff + 6.D0*mt**3*p5Dp34**3*den5x34**2*
     &    facuRLdiff + 6.D0*mt**3*s34*p5Dp34**2*den5x34**2*facuRLdiff
     &     + 3.D0/2.D0*mt**3*s34**2*p5Dp34*den5x34**2*facuRLdiff + 4.D0
     &    *facuLl*mt - 4.D0*facuLl*mt**3*p5Dp34*den5x34 + 6.D0*facuLl*
     &    mt**3*p5Dp34**3*den5x34**2 - 2.D0*facuLl*mt**3*s34*den5x34 +
     &    6.D0*facuLl*mt**3*s34*p5Dp34**2*den5x34**2 + 3.D0/2.D0*facuLl
     &    *mt**3*s34**2*p5Dp34*den5x34**2 )
      vert25x2 = vert25x2 - fp(ep)*mt*p5Dp34*den5x34*facuRLdiff - fp(ep
     &    )*mt**3*den5x34*facuRLdiff + fp(ep)*mt**3*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) + fp(ep)*mt**5*den5x34*facuRLdiff*
     &    s345**(-1) - fp(ep)*facuLl*mt*p5Dp34*den5x34 - fp(ep)*facuLl*
     &    mt**3*den5x34 + fp(ep)*facuLl*mt**3*p5Dp34*den5x34*s345**(-1)
     &     + fp(ep)*facuLl*mt**5*den5x34*s345**(-1) + 2.D0*eploopI3(mtsq,
     &    s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*p5Dp34**3*den5x34**2*
     &    facuRLdiff + 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &    *mt**3*s34*p5Dp34**2*den5x34**2*facuRLdiff + 1.D0/2.D0*
     &    eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*s34**2*
     &    p5Dp34*den5x34**2*facuRLdiff + 2.D0*eploopI3(mtsq,s34,s345,zip,
     &    mtsq,mtsq,musq,ep)*facuLl*mt**3*p5Dp34**3*den5x34**2 + 2.D0*
     &    eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**3*s34*
     &    p5Dp34**2*den5x34**2 + 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,
     &    mtsq,mtsq,musq,ep)*facuLl*mt**3*s34**2*p5Dp34*den5x34**2 + 2.D
     &    0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*
     &    p5Dp34**3*den5x34**2*facuRLdiff
      vert25x2 = vert25x2 + 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     & musq,ep)*mt**3*s34*p5Dp34**2*den5x34**2*facuRLdiff + 1.D0/2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*s34**2*
     &    p5Dp34*den5x34**2*facuRLdiff + 2.D0*ep2loopI3(mtsq,s34,s345,zip
     &    ,mtsq,mtsq,musq,ep)*facuLl*mt**3*p5Dp34**3*den5x34**2 + 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**3*s34
     &    *p5Dp34**2*den5x34**2 + 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,
     &    mtsq,mtsq,musq,ep)*facuLl*mt**3*s34**2*p5Dp34*den5x34**2

      vert25x3= + loopI2(s34,mtsq,mtsq,musq,ep) * ( mt*s34*den5x34*
     &    facuRLdiff + facuLl*mt*s34*den5x34 )
      vert25x3 = vert25x3 + loopI2(mtsq,zip,mtsq,musq,ep) * ( mt*p5Dp34*
     &    den5x34*facuRLdiff + facuLl*mt*p5Dp34*den5x34 )
      vert25x3 = vert25x3 + loopI2(s345,zip,mtsq,musq,ep) * (  - mt*
     &    p5Dp34*den5x34*facuRLdiff - mt*s34*den5x34*facuRLdiff -
     &    facuLl*mt*p5Dp34*den5x34 - facuLl*mt*s34*den5x34 )
      vert25x3 = vert25x3 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * (  - 2.D0*mt*facuRLdiff - 2.D0*mt*p5Dp34**2*den5x34*
     &    facuRLdiff - mt*s34*p5Dp34*den5x34*facuRLdiff - 2.D0*facuLl*
     &    mt - 2.D0*facuLl*mt*p5Dp34**2*den5x34 - facuLl*mt*s34*p5Dp34*
     &    den5x34 )

      vert25x4= + loopI2(s34,mtsq,mtsq,musq,ep) * (  - mt*s34*den5x34*
     &    facuRLdiff )
      vert25x4 = vert25x4 + loopI2(mtsq,zip,mtsq,musq,ep) * (  - mt*
     &    p5Dp34*den5x34*facuRLdiff )
      vert25x4 = vert25x4 + loopI2(s345,zip,mtsq,musq,ep) * ( mt*p5Dp34*
     &    den5x34*facuRLdiff + mt*s34*den5x34*facuRLdiff )
      vert25x4 = vert25x4 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * ( 2.D0*mt*facuRLdiff + 2.D0*mt*p5Dp34**2*den5x34*facuRLdiff
     &     + mt*s34*p5Dp34*den5x34*facuRLdiff )

      vert25x5= + loopI2(s34,mtsq,mtsq,musq,ep) * (  - p5Dp34*den5x34*
     &    facuRLdiff - p5Dp34**3*den5x34**2*facuRLdiff + 2.D0*s34*
     &    den5x34*facuRLdiff - 1.D0/2.D0*s34*p5Dp34**2*den5x34**2*
     &    facuRLdiff - 2.D0*mt**2*s34*p5Dp34*den5x34**2*facuRLdiff -
     &    mt**2*s34**2*den5x34**2*facuRLdiff - facuLl*p5Dp34*den5x34 -
     &    facuLl*p5Dp34**3*den5x34**2 + 2.D0*facuLl*s34*den5x34 - 1.D0/
     &    2.D0*facuLl*s34*p5Dp34**2*den5x34**2 - 2.D0*facuLl*mt**2*s34*
     &    p5Dp34*den5x34**2 - facuLl*mt**2*s34**2*den5x34**2 )
      vert25x5 = vert25x5 + loopI2(mtsq,zip,mtsq,musq,ep) * ( 2.D0*p5Dp34
     &    *den5x34*facuRLdiff - mt**2*den5x34*facuRLdiff - mt**2*p5Dp34
     &    *den5x34*facuRLdiff*s345**(-1) - 3.D0*mt**2*p5Dp34**2*
     &    den5x34**2*facuRLdiff - 3.D0/2.D0*mt**2*s34*p5Dp34*den5x34**2
     &    *facuRLdiff - mt**4*den5x34*facuRLdiff*s345**(-1) + 2.D0*
     &    facuLl*p5Dp34*den5x34 - facuLl*mt**2*den5x34 - facuLl*mt**2*
     &    p5Dp34*den5x34*s345**(-1) - 3.D0*facuLl*mt**2*p5Dp34**2*
     &    den5x34**2 - 3.D0/2.D0*facuLl*mt**2*s34*p5Dp34*den5x34**2 -
     &    facuLl*mt**4*den5x34*s345**(-1) )
      vert25x5 = vert25x5 + loopI2(s345,zip,mtsq,musq,ep) * (  - p5Dp34*
     &    den5x34*facuRLdiff + p5Dp34**3*den5x34**2*facuRLdiff - 2.D0*
     &    s34*den5x34*facuRLdiff + 1.D0/2.D0*s34*p5Dp34**2*den5x34**2*
     &    facuRLdiff + mt**2*den5x34*facuRLdiff + mt**2*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) + 3.D0*mt**2*p5Dp34**2*den5x34**2*
     &    facuRLdiff + 7.D0/2.D0*mt**2*s34*p5Dp34*den5x34**2*facuRLdiff
     &     + mt**2*s34**2*den5x34**2*facuRLdiff + mt**4*den5x34*
     &    facuRLdiff*s345**(-1) - facuLl*p5Dp34*den5x34 + facuLl*
     &    p5Dp34**3*den5x34**2 - 2.D0*facuLl*s34*den5x34 + 1.D0/2.D0*
     &    facuLl*s34*p5Dp34**2*den5x34**2 + facuLl*mt**2*den5x34 +
     &    facuLl*mt**2*p5Dp34*den5x34*s345**(-1) + 3.D0*facuLl*mt**2*
     &    p5Dp34**2*den5x34**2 + 7.D0/2.D0*facuLl*mt**2*s34*p5Dp34*
     &    den5x34**2 + facuLl*mt**2*s34**2*den5x34**2 + facuLl*mt**4*
     &    den5x34*s345**(-1) )
      vert25x5 = vert25x5 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * (  - 4.D0*facuRLdiff - 4.D0*p5Dp34**2*den5x34*facuRLdiff - 2.D
     &    0*s34*p5Dp34*den5x34*facuRLdiff + 4.D0*mt**2*p5Dp34*den5x34*
     &    facuRLdiff + 6.D0*mt**2*p5Dp34**3*den5x34**2*facuRLdiff + 2.D0
     &    *mt**2*s34*den5x34*facuRLdiff + 6.D0*mt**2*s34*p5Dp34**2*
     &    den5x34**2*facuRLdiff + 3.D0/2.D0*mt**2*s34**2*p5Dp34*
     &    den5x34**2*facuRLdiff - 4.D0*facuLl - 4.D0*facuLl*p5Dp34**2*
     &    den5x34 - 2.D0*facuLl*s34*p5Dp34*den5x34 + 4.D0*facuLl*mt**2*
     &    p5Dp34*den5x34 + 6.D0*facuLl*mt**2*p5Dp34**3*den5x34**2 + 2.D0
     &    *facuLl*mt**2*s34*den5x34 + 6.D0*facuLl*mt**2*s34*p5Dp34**2*
     &    den5x34**2 + 3.D0/2.D0*facuLl*mt**2*s34**2*p5Dp34*den5x34**2
     &     )
      vert25x5 = vert25x5 - fp(ep)*p5Dp34*den5x34*facuRLdiff - fp(ep)*
     &    mt**2*den5x34*facuRLdiff + fp(ep)*mt**2*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) + fp(ep)*mt**4*den5x34*facuRLdiff*
     &    s345**(-1) - fp(ep)*facuLl*p5Dp34*den5x34 - fp(ep)*facuLl*
     &    mt**2*den5x34 + fp(ep)*facuLl*mt**2*p5Dp34*den5x34*s345**(-1)
     &     + fp(ep)*facuLl*mt**4*den5x34*s345**(-1) + 2.D0*eploopI3(mtsq,
     &    s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*p5Dp34**3*den5x34**2*
     &    facuRLdiff + 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &    *mt**2*s34*p5Dp34**2*den5x34**2*facuRLdiff + 1.D0/2.D0*
     &    eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34**2*
     &    p5Dp34*den5x34**2*facuRLdiff + 2.D0*eploopI3(mtsq,s34,s345,zip,
     &    mtsq,mtsq,musq,ep)*facuLl*mt**2*p5Dp34**3*den5x34**2 + 2.D0*
     &    eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**2*s34*
     &    p5Dp34**2*den5x34**2 + 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,
     &    mtsq,mtsq,musq,ep)*facuLl*mt**2*s34**2*p5Dp34*den5x34**2 + 2.D
     &    0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*
     &    p5Dp34**3*den5x34**2*facuRLdiff
      vert25x5 = vert25x5 + 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     & musq,ep)*mt**2*s34*p5Dp34**2*den5x34**2*facuRLdiff + 1.D0/2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34**2*
     &    p5Dp34*den5x34**2*facuRLdiff + 2.D0*ep2loopI3(mtsq,s34,s345,zip
     &    ,mtsq,mtsq,musq,ep)*facuLl*mt**2*p5Dp34**3*den5x34**2 + 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**2*s34
     &    *p5Dp34**2*den5x34**2 + 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,
     &    mtsq,mtsq,musq,ep)*facuLl*mt**2*s34**2*p5Dp34*den5x34**2

      vert25x6= + loopI2(s34,mtsq,mtsq,musq,ep) * (  - p5Dp34*den5x34*
     &    facuRLdiff + mt**2*den5x34*facuRLdiff + 3.D0*mt**2*p5Dp34**2*
     &    den5x34**2*facuRLdiff + 3.D0/2.D0*mt**2*s34*p5Dp34*den5x34**2
     &    *facuRLdiff - facuLl*p5Dp34*den5x34 + facuLl*mt**2*den5x34 +
     &    3.D0*facuLl*mt**2*p5Dp34**2*den5x34**2 + 3.D0/2.D0*facuLl*
     &    mt**2*s34*p5Dp34*den5x34**2 )
      vert25x6 = vert25x6 + loopI2(mtsq,zip,mtsq,musq,ep) * (  - 2.D0*
     &    mt**2*den5x34*facuRLdiff - mt**2*p5Dp34*den5x34*facuRLdiff*
     &    s345**(-1) - mt**4*den5x34*facuRLdiff*s345**(-1) + 3.D0*mt**4
     &    *p5Dp34*den5x34**2*facuRLdiff + 3.D0/2.D0*mt**4*s34*
     &    den5x34**2*facuRLdiff - 2.D0*facuLl*mt**2*den5x34 - facuLl*
     &    mt**2*p5Dp34*den5x34*s345**(-1) - facuLl*mt**4*den5x34*
     &    s345**(-1) + 3.D0*facuLl*mt**4*p5Dp34*den5x34**2 + 3.D0/2.D0*
     &    facuLl*mt**4*s34*den5x34**2 )
      vert25x6 = vert25x6 + loopI2(s345,zip,mtsq,musq,ep) * ( p5Dp34*
     &    den5x34*facuRLdiff + mt**2*den5x34*facuRLdiff + mt**2*p5Dp34*
     &    den5x34*facuRLdiff*s345**(-1) - 3.D0*mt**2*p5Dp34**2*
     &    den5x34**2*facuRLdiff - 3.D0/2.D0*mt**2*s34*p5Dp34*den5x34**2
     &    *facuRLdiff + mt**4*den5x34*facuRLdiff*s345**(-1) - 3.D0*
     &    mt**4*p5Dp34*den5x34**2*facuRLdiff - 3.D0/2.D0*mt**4*s34*
     &    den5x34**2*facuRLdiff + facuLl*p5Dp34*den5x34 + facuLl*mt**2*
     &    den5x34 + facuLl*mt**2*p5Dp34*den5x34*s345**(-1) - 3.D0*
     &    facuLl*mt**2*p5Dp34**2*den5x34**2 - 3.D0/2.D0*facuLl*mt**2*
     &    s34*p5Dp34*den5x34**2 + facuLl*mt**4*den5x34*s345**(-1) - 3.D0
     &    *facuLl*mt**4*p5Dp34*den5x34**2 - 3.D0/2.D0*facuLl*mt**4*s34*
     &    den5x34**2 )
      vert25x6 = vert25x6 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * ( 4.D0*mt**2*p5Dp34*den5x34*facuRLdiff + 2.D0*mt**2*s34*
     &    den5x34*facuRLdiff - 6.D0*mt**4*p5Dp34**2*den5x34**2*
     &    facuRLdiff - 6.D0*mt**4*s34*p5Dp34*den5x34**2*facuRLdiff - 3.D
     &    0/2.D0*mt**4*s34**2*den5x34**2*facuRLdiff + 4.D0*facuLl*mt**2
     &    *p5Dp34*den5x34 + 2.D0*facuLl*mt**2*s34*den5x34 - 6.D0*facuLl
     &    *mt**4*p5Dp34**2*den5x34**2 - 6.D0*facuLl*mt**4*s34*p5Dp34*
     &    den5x34**2 - 3.D0/2.D0*facuLl*mt**4*s34**2*den5x34**2 )
      vert25x6 = vert25x6 + fp(ep)*mt**2*den5x34*facuRLdiff + fp(ep)*
     &    mt**2*p5Dp34*den5x34*facuRLdiff*s345**(-1) + fp(ep)*mt**4*
     &    den5x34*facuRLdiff*s345**(-1) + fp(ep)*facuLl*mt**2*den5x34
     &     + fp(ep)*facuLl*mt**2*p5Dp34*den5x34*s345**(-1) + fp(ep)*
     &    facuLl*mt**4*den5x34*s345**(-1) - 2.D0*eploopI3(mtsq,s34,s345,
     &    zip,mtsq,mtsq,musq,ep)*mt**4*p5Dp34**2*den5x34**2*facuRLdiff
     &     - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34
     &    *p5Dp34*den5x34**2*facuRLdiff - 1.D0/2.D0*eploopI3(mtsq,s34,
     &    s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34**2*den5x34**2*
     &    facuRLdiff - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &    *facuLl*mt**4*p5Dp34**2*den5x34**2 - 2.D0*eploopI3(mtsq,s34,
     &    s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**4*s34*p5Dp34*
     &    den5x34**2 - 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     &    musq,ep)*facuLl*mt**4*s34**2*den5x34**2 - 2.D0*ep2loopI3(mtsq,
     &    s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*p5Dp34**2*den5x34**2*
     &    facuRLdiff
      vert25x6 = vert25x6 - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     & musq,ep)*mt**4*s34*p5Dp34*den5x34**2*facuRLdiff - 1.D0/2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34**2*
     &    den5x34**2*facuRLdiff - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*facuLl*mt**4*p5Dp34**2*den5x34**2 - 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**4*s34
     &    *p5Dp34*den5x34**2 - 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq
     &    ,mtsq,musq,ep)*facuLl*mt**4*s34**2*den5x34**2

      vert25x7= + loopI2(s34,mtsq,mtsq,musq,ep) * ( p5Dp34*den5x34*
     &    facuRLdiff + p5Dp34**3*den5x34**2*facuRLdiff - 2.D0*s34*
     &    den5x34*facuRLdiff + 1.D0/2.D0*s34*p5Dp34**2*den5x34**2*
     &    facuRLdiff + 2.D0*mt**2*s34*p5Dp34*den5x34**2*facuRLdiff +
     &    mt**2*s34**2*den5x34**2*facuRLdiff )
      vert25x7 = vert25x7 + loopI2(mtsq,zip,mtsq,musq,ep) * (  - 2.D0*
     &    p5Dp34*den5x34*facuRLdiff + mt**2*den5x34*facuRLdiff + mt**2*
     &    p5Dp34*den5x34*facuRLdiff*s345**(-1) + 3.D0*mt**2*p5Dp34**2*
     &    den5x34**2*facuRLdiff + 3.D0/2.D0*mt**2*s34*p5Dp34*den5x34**2
     &    *facuRLdiff + mt**4*den5x34*facuRLdiff*s345**(-1) )
      vert25x7 = vert25x7 + loopI2(s345,zip,mtsq,musq,ep) * ( p5Dp34*
     &    den5x34*facuRLdiff - p5Dp34**3*den5x34**2*facuRLdiff + 2.D0*
     &    s34*den5x34*facuRLdiff - 1.D0/2.D0*s34*p5Dp34**2*den5x34**2*
     &    facuRLdiff - mt**2*den5x34*facuRLdiff - mt**2*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) - 3.D0*mt**2*p5Dp34**2*den5x34**2*
     &    facuRLdiff - 7.D0/2.D0*mt**2*s34*p5Dp34*den5x34**2*facuRLdiff
     &     - mt**2*s34**2*den5x34**2*facuRLdiff - mt**4*den5x34*
     &    facuRLdiff*s345**(-1) )
      vert25x7 = vert25x7 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * ( 4.D0*facuRLdiff + 4.D0*p5Dp34**2*den5x34*facuRLdiff + 2.D0*
     &    s34*p5Dp34*den5x34*facuRLdiff - 4.D0*mt**2*p5Dp34*den5x34*
     &    facuRLdiff - 6.D0*mt**2*p5Dp34**3*den5x34**2*facuRLdiff - 2.D0
     &    *mt**2*s34*den5x34*facuRLdiff - 6.D0*mt**2*s34*p5Dp34**2*
     &    den5x34**2*facuRLdiff - 3.D0/2.D0*mt**2*s34**2*p5Dp34*
     &    den5x34**2*facuRLdiff )
      vert25x7 = vert25x7 + fp(ep)*p5Dp34*den5x34*facuRLdiff + fp(ep)*
     &    mt**2*den5x34*facuRLdiff - fp(ep)*mt**2*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) - fp(ep)*mt**4*den5x34*facuRLdiff*
     &    s345**(-1) - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &    *mt**2*p5Dp34**3*den5x34**2*facuRLdiff - 2.D0*eploopI3(mtsq,s34
     &    ,s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34*p5Dp34**2*den5x34**2*
     &    facuRLdiff - 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     &    musq,ep)*mt**2*s34**2*p5Dp34*den5x34**2*facuRLdiff - 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*p5Dp34**3*
     &    den5x34**2*facuRLdiff - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*mt**2*s34*p5Dp34**2*den5x34**2*facuRLdiff - 1.D0
     &    /2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*
     &    s34**2*p5Dp34*den5x34**2*facuRLdiff

      vert25x8= + loopI2(s34,mtsq,mtsq,musq,ep) * ( p5Dp34*den5x34*
     &    facuRLdiff - mt**2*den5x34*facuRLdiff - 3.D0*mt**2*p5Dp34**2*
     &    den5x34**2*facuRLdiff - 3.D0/2.D0*mt**2*s34*p5Dp34*den5x34**2
     &    *facuRLdiff )
      vert25x8 = vert25x8 + loopI2(mtsq,zip,mtsq,musq,ep) * ( 2.D0*mt**2*
     &    den5x34*facuRLdiff + mt**2*p5Dp34*den5x34*facuRLdiff*
     &    s345**(-1) + mt**4*den5x34*facuRLdiff*s345**(-1) - 3.D0*mt**4
     &    *p5Dp34*den5x34**2*facuRLdiff - 3.D0/2.D0*mt**4*s34*
     &    den5x34**2*facuRLdiff )
      vert25x8 = vert25x8 + loopI2(s345,zip,mtsq,musq,ep) * (  - p5Dp34*
     &    den5x34*facuRLdiff - mt**2*den5x34*facuRLdiff - mt**2*p5Dp34*
     &    den5x34*facuRLdiff*s345**(-1) + 3.D0*mt**2*p5Dp34**2*
     &    den5x34**2*facuRLdiff + 3.D0/2.D0*mt**2*s34*p5Dp34*den5x34**2
     &    *facuRLdiff - mt**4*den5x34*facuRLdiff*s345**(-1) + 3.D0*
     &    mt**4*p5Dp34*den5x34**2*facuRLdiff + 3.D0/2.D0*mt**4*s34*
     &    den5x34**2*facuRLdiff )
      vert25x8 = vert25x8 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * (  - 4.D0*mt**2*p5Dp34*den5x34*facuRLdiff - 2.D0*mt**2*s34*
     &    den5x34*facuRLdiff + 6.D0*mt**4*p5Dp34**2*den5x34**2*
     &    facuRLdiff + 6.D0*mt**4*s34*p5Dp34*den5x34**2*facuRLdiff + 3.D
     &    0/2.D0*mt**4*s34**2*den5x34**2*facuRLdiff )
      vert25x8 = vert25x8 - fp(ep)*mt**2*den5x34*facuRLdiff - fp(ep)*
     &    mt**2*p5Dp34*den5x34*facuRLdiff*s345**(-1) - fp(ep)*mt**4*
     &    den5x34*facuRLdiff*s345**(-1) + 2.D0*eploopI3(mtsq,s34,s345,zip
     &    ,mtsq,mtsq,musq,ep)*mt**4*p5Dp34**2*den5x34**2*facuRLdiff + 2.
     &    D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34*
     &    p5Dp34*den5x34**2*facuRLdiff + 1.D0/2.D0*eploopI3(mtsq,s34,s345
     &    ,zip,mtsq,mtsq,musq,ep)*mt**4*s34**2*den5x34**2*facuRLdiff +
     &    2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*
     &    p5Dp34**2*den5x34**2*facuRLdiff + 2.D0*ep2loopI3(mtsq,s34,s345,
     &    zip,mtsq,mtsq,musq,ep)*mt**4*s34*p5Dp34*den5x34**2*facuRLdiff
     &     + 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*
     &    mt**4*s34**2*den5x34**2*facuRLdiff

      vert25x9= + loopI2(s34,mtsq,mtsq,musq,ep) * ( facuRLdiff + 2.D0*
     &    p5Dp34**2*den5x34*facuRLdiff + p5Dp34**4*den5x34**2*
     &    facuRLdiff - s34*p5Dp34*den5x34*facuRLdiff + 1.D0/2.D0*s34*
     &    p5Dp34**3*den5x34**2*facuRLdiff - 2.D0*mt**2*s34*den5x34*
     &    facuRLdiff - mt**2*s34*p5Dp34**2*den5x34**2*facuRLdiff - 1.D0/
     &    2.D0*mt**2*s34**2*p5Dp34*den5x34**2*facuRLdiff + facuLl + 2.D0
     &    *facuLl*p5Dp34**2*den5x34 + facuLl*p5Dp34**4*den5x34**2 -
     &    facuLl*s34*p5Dp34*den5x34 + 1.D0/2.D0*facuLl*s34*p5Dp34**3*
     &    den5x34**2 - 2.D0*facuLl*mt**2*s34*den5x34 - facuLl*mt**2*s34
     &    *p5Dp34**2*den5x34**2 - 1.D0/2.D0*facuLl*mt**2*s34**2*p5Dp34*
     &    den5x34**2 )
      vert25x9 = vert25x9 + loopI2(mtsq,zip,mtsq,musq,ep) * (  - 2.D0*
     &    p5Dp34**2*den5x34*facuRLdiff + 1.D0/2.D0*mt**2*p5Dp34*den5x34
     &    *facuRLdiff + mt**2*p5Dp34**2*den5x34*facuRLdiff*s345**(-1)
     &     + 2.D0*mt**2*p5Dp34**3*den5x34**2*facuRLdiff + 3.D0/2.D0*
     &    mt**2*s34*den5x34*facuRLdiff + 1.D0/2.D0*mt**2*s34*p5Dp34*
     &    den5x34*facuRLdiff*s345**(-1) + mt**2*s34*p5Dp34**2*
     &    den5x34**2*facuRLdiff + 1.D0/2.D0*mt**4*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) - 2.D0*mt**4*s34*p5Dp34*den5x34**2*
     &    facuRLdiff - mt**4*s34**2*den5x34**2*facuRLdiff - 2.D0*facuLl
     &    *p5Dp34**2*den5x34 + 1.D0/2.D0*facuLl*mt**2*p5Dp34*den5x34 +
     &    facuLl*mt**2*p5Dp34**2*den5x34*s345**(-1) + 2.D0*facuLl*mt**2
     &    *p5Dp34**3*den5x34**2 + 3.D0/2.D0*facuLl*mt**2*s34*den5x34 +
     &    1.D0/2.D0*facuLl*mt**2*s34*p5Dp34*den5x34*s345**(-1) + facuLl
     &    *mt**2*s34*p5Dp34**2*den5x34**2 + 1.D0/2.D0*facuLl*mt**4*
     &    p5Dp34*den5x34*s345**(-1) - 2.D0*facuLl*mt**4*s34*p5Dp34*
     &    den5x34**2 )
      vert25x9 = vert25x9 + loopI2(mtsq,zip,mtsq,musq,ep) * (  - facuLl*
     &    mt**4*s34**2*den5x34**2 )
      vert25x9 = vert25x9 + loopI2(s345,zip,mtsq,musq,ep) * (  -
     &    p5Dp34**4*den5x34**2*facuRLdiff + s34*p5Dp34*den5x34*
     &    facuRLdiff - 1.D0/2.D0*s34*p5Dp34**3*den5x34**2*facuRLdiff -
     &    1.D0/2.D0*mt**2*p5Dp34*den5x34*facuRLdiff - mt**2*p5Dp34**2*
     &    den5x34*facuRLdiff*s345**(-1) - 2.D0*mt**2*p5Dp34**3*
     &    den5x34**2*facuRLdiff + 1.D0/2.D0*mt**2*s34*den5x34*
     &    facuRLdiff - 1.D0/2.D0*mt**2*s34*p5Dp34*den5x34*facuRLdiff*
     &    s345**(-1) + 1.D0/2.D0*mt**2*s34**2*p5Dp34*den5x34**2*
     &    facuRLdiff - 1.D0/2.D0*mt**4*p5Dp34*den5x34*facuRLdiff*
     &    s345**(-1) + 2.D0*mt**4*s34*p5Dp34*den5x34**2*facuRLdiff +
     &    mt**4*s34**2*den5x34**2*facuRLdiff - facuLl*p5Dp34**4*
     &    den5x34**2 + facuLl*s34*p5Dp34*den5x34 - 1.D0/2.D0*facuLl*s34
     &    *p5Dp34**3*den5x34**2 - 1.D0/2.D0*facuLl*mt**2*p5Dp34*den5x34
     &     - facuLl*mt**2*p5Dp34**2*den5x34*s345**(-1) - 2.D0*facuLl*
     &    mt**2*p5Dp34**3*den5x34**2 + 1.D0/2.D0*facuLl*mt**2*s34*
     &    den5x34 )
      vert25x9 = vert25x9 + loopI2(s345,zip,mtsq,musq,ep) * (  - 1.D0/2.D0
     &    *facuLl*mt**2*s34*p5Dp34*den5x34*s345**(-1) + 1.D0/2.D0*
     &    facuLl*mt**2*s34**2*p5Dp34*den5x34**2 - 1.D0/2.D0*facuLl*
     &    mt**4*p5Dp34*den5x34*s345**(-1) + 2.D0*facuLl*mt**4*s34*
     &    p5Dp34*den5x34**2 + facuLl*mt**4*s34**2*den5x34**2 )
      vert25x9 = vert25x9 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * ( 4.D0*p5Dp34*facuRLdiff + 4.D0*p5Dp34**3*den5x34*facuRLdiff
     &     + 2.D0*s34*p5Dp34**2*den5x34*facuRLdiff + 2.D0*mt**2*
     &    facuRLdiff - 2.D0*mt**2*p5Dp34**2*den5x34*facuRLdiff - 4.D0*
     &    mt**2*p5Dp34**4*den5x34**2*facuRLdiff - 4.D0*mt**2*s34*p5Dp34
     &    *den5x34*facuRLdiff - 4.D0*mt**2*s34*p5Dp34**3*den5x34**2*
     &    facuRLdiff - 3.D0/2.D0*mt**2*s34**2*den5x34*facuRLdiff -
     &    mt**2*s34**2*p5Dp34**2*den5x34**2*facuRLdiff + 4.D0*mt**4*s34
     &    *p5Dp34**2*den5x34**2*facuRLdiff + 4.D0*mt**4*s34**2*p5Dp34*
     &    den5x34**2*facuRLdiff + mt**4*s34**3*den5x34**2*facuRLdiff +
     &    4.D0*facuLl*p5Dp34 + 4.D0*facuLl*p5Dp34**3*den5x34 + 2.D0*
     &    facuLl*s34*p5Dp34**2*den5x34 - 2.D0*facuLl*mt**2*p5Dp34**2*
     &    den5x34 - 4.D0*facuLl*mt**2*p5Dp34**4*den5x34**2 - 4.D0*
     &    facuLl*mt**2*s34*p5Dp34*den5x34 - 4.D0*facuLl*mt**2*s34*
     &    p5Dp34**3*den5x34**2 - 3.D0/2.D0*facuLl*mt**2*s34**2*den5x34
     &     - facuLl*mt**2*s34**2*p5Dp34**2*den5x34**2 )
      vert25x9 = vert25x9 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * ( 4.D0*facuLl*mt**4*s34*p5Dp34**2*den5x34**2 + 4.D0*facuLl*
     &    mt**4*s34**2*p5Dp34*den5x34**2 + facuLl*mt**4*s34**3*
     &    den5x34**2 )
      vert25x9 = vert25x9 + fp(ep)*p5Dp34**2*den5x34*facuRLdiff + 1.D0/
     &    2.D0*fp(ep)*mt**2*p5Dp34*den5x34*facuRLdiff - fp(ep)*mt**2*
     &    p5Dp34**2*den5x34*facuRLdiff*s345**(-1) - fp(ep)*mt**2*s34*
     &    den5x34*facuRLdiff - 1.D0/2.D0*fp(ep)*mt**2*s34*p5Dp34*
     &    den5x34*facuRLdiff*s345**(-1) - 1.D0/2.D0*fp(ep)*mt**4*p5Dp34
     &    *den5x34*facuRLdiff*s345**(-1) + fp(ep)*facuLl*p5Dp34**2*
     &    den5x34 + 1.D0/2.D0*fp(ep)*facuLl*mt**2*p5Dp34*den5x34 - fp(
     &    ep)*facuLl*mt**2*p5Dp34**2*den5x34*s345**(-1) - fp(ep)*facuLl
     &    *mt**2*s34*den5x34 - 1.D0/2.D0*fp(ep)*facuLl*mt**2*s34*p5Dp34
     &    *den5x34*s345**(-1) - 1.D0/2.D0*fp(ep)*facuLl*mt**4*p5Dp34*
     &    den5x34*s345**(-1) - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     &    musq,ep)*mt**2*p5Dp34**2*den5x34*facuRLdiff - 2.D0*eploopI3(
     &    mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*p5Dp34**4*
     &    den5x34**2*facuRLdiff - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*mt**2*s34*p5Dp34*den5x34*facuRLdiff - 2.D0*
     &    eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34*
     &    p5Dp34**3*den5x34**2*facuRLdiff
      vert25x9 = vert25x9 - 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,
     & mtsq,musq,ep)*mt**2*s34**2*den5x34*facuRLdiff - 1.D0/2.D0*
     &    eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34**2*
     &    p5Dp34**2*den5x34**2*facuRLdiff + 2.D0*eploopI3(mtsq,s34,s345,
     &    zip,mtsq,mtsq,musq,ep)*mt**4*s34*p5Dp34**2*den5x34**2*
     &    facuRLdiff + 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &    *mt**4*s34**2*p5Dp34*den5x34**2*facuRLdiff + 1.D0/2.D0*
     &    eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34**3*
     &    den5x34**2*facuRLdiff - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*facuLl*mt**2*p5Dp34**2*den5x34 - 2.D0*eploopI3(
     &    mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**2*p5Dp34**4*
     &    den5x34**2 - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &    *facuLl*mt**2*s34*p5Dp34*den5x34 - 2.D0*eploopI3(mtsq,s34,s345,
     &    zip,mtsq,mtsq,musq,ep)*facuLl*mt**2*s34*p5Dp34**3*den5x34**2
     &     - 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*
     &    facuLl*mt**2*s34**2*den5x34
      vert25x9 = vert25x9 - 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,
     & mtsq,musq,ep)*facuLl*mt**2*s34**2*p5Dp34**2*den5x34**2 + 2.D0*
     &    eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**4*s34*
     &    p5Dp34**2*den5x34**2 + 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*facuLl*mt**4*s34**2*p5Dp34*den5x34**2 + 1.D0/2.D
     &    0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**4*
     &    s34**3*den5x34**2 - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     &    musq,ep)*mt**2*p5Dp34**2*den5x34*facuRLdiff - 2.D0*ep2loopI3(
     &    mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*p5Dp34**4*
     &    den5x34**2*facuRLdiff - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*mt**2*s34*p5Dp34*den5x34*facuRLdiff - 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34*
     &    p5Dp34**3*den5x34**2*facuRLdiff - 1.D0/2.D0*ep2loopI3(mtsq,s34,
     &    s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34**2*den5x34*facuRLdiff
     &     - 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*
     &    mt**2*s34**2*p5Dp34**2*den5x34**2*facuRLdiff
      vert25x9 = vert25x9 + 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     & musq,ep)*mt**4*s34*p5Dp34**2*den5x34**2*facuRLdiff + 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34**2*
     &    p5Dp34*den5x34**2*facuRLdiff + 1.D0/2.D0*ep2loopI3(mtsq,s34,
     &    s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34**3*den5x34**2*
     &    facuRLdiff - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep
     &    )*facuLl*mt**2*p5Dp34**2*den5x34 - 2.D0*ep2loopI3(mtsq,s34,s345
     &    ,zip,mtsq,mtsq,musq,ep)*facuLl*mt**2*p5Dp34**4*den5x34**2 - 2.
     &    D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**2*
     &    s34*p5Dp34*den5x34 - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq
     &    ,musq,ep)*facuLl*mt**2*s34*p5Dp34**3*den5x34**2 - 1.D0/2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**2*
     &    s34**2*den5x34 - 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*facuLl*mt**2*s34**2*p5Dp34**2*den5x34**2 + 2.D0
     &    *ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**4*
     &    s34*p5Dp34**2*den5x34**2
      vert25x9 = vert25x9 + 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     & musq,ep)*facuLl*mt**4*s34**2*p5Dp34*den5x34**2 + 1.D0/2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*facuLl*mt**4*
     &    s34**3*den5x34**2

      vert25x10= + loopI2(s34,mtsq,mtsq,musq,ep) * (  - facuRLdiff - 2.D0
     &    *p5Dp34**2*den5x34*facuRLdiff - p5Dp34**4*den5x34**2*
     &    facuRLdiff + s34*p5Dp34*den5x34*facuRLdiff - 1.D0/2.D0*s34*
     &    p5Dp34**3*den5x34**2*facuRLdiff + 2.D0*mt**2*s34*den5x34*
     &    facuRLdiff + mt**2*s34*p5Dp34**2*den5x34**2*facuRLdiff + 1.D0/
     &    2.D0*mt**2*s34**2*p5Dp34*den5x34**2*facuRLdiff )
      vert25x10 = vert25x10 + loopI2(mtsq,zip,mtsq,musq,ep) * ( 2.D0*
     &    p5Dp34**2*den5x34*facuRLdiff - 1.D0/2.D0*mt**2*p5Dp34*den5x34
     &    *facuRLdiff - mt**2*p5Dp34**2*den5x34*facuRLdiff*s345**(-1)
     &     - 2.D0*mt**2*p5Dp34**3*den5x34**2*facuRLdiff - 3.D0/2.D0*
     &    mt**2*s34*den5x34*facuRLdiff - 1.D0/2.D0*mt**2*s34*p5Dp34*
     &    den5x34*facuRLdiff*s345**(-1) - mt**2*s34*p5Dp34**2*
     &    den5x34**2*facuRLdiff - 1.D0/2.D0*mt**4*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) + 2.D0*mt**4*s34*p5Dp34*den5x34**2*
     &    facuRLdiff + mt**4*s34**2*den5x34**2*facuRLdiff )
      vert25x10 = vert25x10 + loopI2(s345,zip,mtsq,musq,ep) * ( p5Dp34**4
     &    *den5x34**2*facuRLdiff - s34*p5Dp34*den5x34*facuRLdiff + 1.D0/
     &    2.D0*s34*p5Dp34**3*den5x34**2*facuRLdiff + 1.D0/2.D0*mt**2*
     &    p5Dp34*den5x34*facuRLdiff + mt**2*p5Dp34**2*den5x34*
     &    facuRLdiff*s345**(-1) + 2.D0*mt**2*p5Dp34**3*den5x34**2*
     &    facuRLdiff - 1.D0/2.D0*mt**2*s34*den5x34*facuRLdiff + 1.D0/2.D
     &    0*mt**2*s34*p5Dp34*den5x34*facuRLdiff*s345**(-1) - 1.D0/2.D0*
     &    mt**2*s34**2*p5Dp34*den5x34**2*facuRLdiff + 1.D0/2.D0*mt**4*
     &    p5Dp34*den5x34*facuRLdiff*s345**(-1) - 2.D0*mt**4*s34*p5Dp34*
     &    den5x34**2*facuRLdiff - mt**4*s34**2*den5x34**2*facuRLdiff )
      vert25x10 = vert25x10 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * (  - 4.D0*p5Dp34*facuRLdiff - 4.D0*p5Dp34**3*den5x34*
     &    facuRLdiff - 2.D0*s34*p5Dp34**2*den5x34*facuRLdiff - 4.D0*
     &    mt**2*facuRLdiff + 2.D0*mt**2*p5Dp34**2*den5x34*facuRLdiff +
     &    4.D0*mt**2*p5Dp34**4*den5x34**2*facuRLdiff + 4.D0*mt**2*s34*
     &    p5Dp34*den5x34*facuRLdiff + 4.D0*mt**2*s34*p5Dp34**3*
     &    den5x34**2*facuRLdiff + 3.D0/2.D0*mt**2*s34**2*den5x34*
     &    facuRLdiff + mt**2*s34**2*p5Dp34**2*den5x34**2*facuRLdiff - 4.
     &    D0*mt**4*s34*p5Dp34**2*den5x34**2*facuRLdiff - 4.D0*mt**4*
     &    s34**2*p5Dp34*den5x34**2*facuRLdiff - mt**4*s34**3*den5x34**2
     &    *facuRLdiff )
      vert25x10 = vert25x10 - fp(ep)*p5Dp34**2*den5x34*facuRLdiff - 1.D0
     &    /2.D0*fp(ep)*mt**2*p5Dp34*den5x34*facuRLdiff + fp(ep)*mt**2*
     &    p5Dp34**2*den5x34*facuRLdiff*s345**(-1) + fp(ep)*mt**2*s34*
     &    den5x34*facuRLdiff + 1.D0/2.D0*fp(ep)*mt**2*s34*p5Dp34*
     &    den5x34*facuRLdiff*s345**(-1) + 1.D0/2.D0*fp(ep)*mt**4*p5Dp34
     &    *den5x34*facuRLdiff*s345**(-1) + 2.D0*eploopI3(mtsq,s34,s345,
     &    zip,mtsq,mtsq,musq,ep)*mt**2*p5Dp34**2*den5x34*facuRLdiff + 2.
     &    D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*
     &    p5Dp34**4*den5x34**2*facuRLdiff + 2.D0*eploopI3(mtsq,s34,s345,
     &    zip,mtsq,mtsq,musq,ep)*mt**2*s34*p5Dp34*den5x34*facuRLdiff +
     &    2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34*
     &    p5Dp34**3*den5x34**2*facuRLdiff + 1.D0/2.D0*eploopI3(mtsq,s34,
     &    s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34**2*den5x34*facuRLdiff
     &     + 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*
     &    mt**2*s34**2*p5Dp34**2*den5x34**2*facuRLdiff - 2.D0*eploopI3(
     &    mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34*p5Dp34**2*
     &    den5x34**2*facuRLdiff
      vert25x10 = vert25x10 - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     & musq,ep)*mt**4*s34**2*p5Dp34*den5x34**2*facuRLdiff - 1.D0/2.D0*
     &    eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34**3*
     &    den5x34**2*facuRLdiff + 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*mt**2*p5Dp34**2*den5x34*facuRLdiff + 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*p5Dp34**4*
     &    den5x34**2*facuRLdiff + 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*mt**2*s34*p5Dp34*den5x34*facuRLdiff + 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34*
     &    p5Dp34**3*den5x34**2*facuRLdiff + 1.D0/2.D0*ep2loopI3(mtsq,s34,
     &    s345,zip,mtsq,mtsq,musq,ep)*mt**2*s34**2*den5x34*facuRLdiff
     &     + 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*
     &    mt**2*s34**2*p5Dp34**2*den5x34**2*facuRLdiff - 2.D0*ep2loopI3(
     &    mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**4*s34*p5Dp34**2*
     &    den5x34**2*facuRLdiff - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*mt**4*s34**2*p5Dp34*den5x34**2*facuRLdiff
      vert25x10 = vert25x10 - 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,
     & mtsq,musq,ep)*mt**4*s34**3*den5x34**2*facuRLdiff

      vert25x11= + loopI2(s34,mtsq,mtsq,musq,ep) * (  - 3.D0*mt*s34*
     &    den5x34*facuRLdiff - 3.D0*mt*s34*p5Dp34**2*den5x34**2*
     &    facuRLdiff - 3.D0/2.D0*mt*s34**2*p5Dp34*den5x34**2*facuRLdiff
     &     )
      vert25x11 = vert25x11 + loopI2(mtsq,zip,mtsq,musq,ep) * (  - 3.D0*
     &    mt*p5Dp34*den5x34*facuRLdiff - 2.D0*mt*p5Dp34**3*den5x34**2*
     &    facuRLdiff - mt*s34*p5Dp34**2*den5x34**2*facuRLdiff - mt**3*
     &    p5Dp34*den5x34*facuRLdiff*s345**(-1) - mt**3*s34*den5x34*
     &    facuRLdiff*s345**(-1) - mt**3*s34*p5Dp34*den5x34**2*
     &    facuRLdiff - 1.D0/2.D0*mt**3*s34**2*den5x34**2*facuRLdiff )
      vert25x11 = vert25x11 + loopI2(s345,zip,mtsq,musq,ep) * ( 3.D0*mt*
     &    p5Dp34*den5x34*facuRLdiff + 2.D0*mt*p5Dp34**3*den5x34**2*
     &    facuRLdiff + 3.D0*mt*s34*den5x34*facuRLdiff + 4.D0*mt*s34*
     &    p5Dp34**2*den5x34**2*facuRLdiff + 3.D0/2.D0*mt*s34**2*p5Dp34*
     &    den5x34**2*facuRLdiff + mt**3*p5Dp34*den5x34*facuRLdiff*
     &    s345**(-1) + mt**3*s34*den5x34*facuRLdiff*s345**(-1) + mt**3*
     &    s34*p5Dp34*den5x34**2*facuRLdiff + 1.D0/2.D0*mt**3*s34**2*
     &    den5x34**2*facuRLdiff )
      vert25x11 = vert25x11 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * ( 4.D0*mt*facuRLdiff + 8.D0*mt*p5Dp34**2*den5x34*facuRLdiff
     &     + 4.D0*mt*p5Dp34**4*den5x34**2*facuRLdiff + 4.D0*mt*s34*
     &    p5Dp34*den5x34*facuRLdiff + 4.D0*mt*s34*p5Dp34**3*den5x34**2*
     &    facuRLdiff + mt*s34**2*p5Dp34**2*den5x34**2*facuRLdiff + 2.D0
     &    *mt**3*s34*p5Dp34**2*den5x34**2*facuRLdiff + 2.D0*mt**3*
     &    s34**2*p5Dp34*den5x34**2*facuRLdiff + 1.D0/2.D0*mt**3*s34**3*
     &    den5x34**2*facuRLdiff )
      vert25x11 = vert25x11 - fp(ep)*mt*p5Dp34*den5x34*facuRLdiff - fp(
     &    ep)*mt*s34*den5x34*facuRLdiff + fp(ep)*mt**3*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) + fp(ep)*mt**3*s34*den5x34*facuRLdiff*
     &    s345**(-1) + 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &    *mt**3*s34*p5Dp34**2*den5x34**2*facuRLdiff + 2.D0*eploopI3(mtsq
     &    ,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*s34**2*p5Dp34*
     &    den5x34**2*facuRLdiff + 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,
     &    mtsq,mtsq,musq,ep)*mt**3*s34**3*den5x34**2*facuRLdiff + 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*s34*
     &    p5Dp34**2*den5x34**2*facuRLdiff + 2.D0*ep2loopI3(mtsq,s34,s345,
     &    zip,mtsq,mtsq,musq,ep)*mt**3*s34**2*p5Dp34*den5x34**2*
     &    facuRLdiff + 1.D0/2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     &    musq,ep)*mt**3*s34**3*den5x34**2*facuRLdiff

      vert25x12= + loopI2(s34,mtsq,mtsq,musq,ep) * ( mt*p5Dp34*den5x34*
     &    facuRLdiff + mt*p5Dp34**3*den5x34**2*facuRLdiff + 1.D0/2.D0*
     &    mt*s34*p5Dp34**2*den5x34**2*facuRLdiff + 2.D0*mt**3*s34*
     &    p5Dp34*den5x34**2*facuRLdiff + mt**3*s34**2*den5x34**2*
     &    facuRLdiff )
      vert25x12 = vert25x12 + loopI2(mtsq,zip,mtsq,musq,ep) * ( mt**3*
     &    den5x34*facuRLdiff + mt**3*p5Dp34*den5x34*facuRLdiff*
     &    s345**(-1) + 3.D0*mt**3*p5Dp34**2*den5x34**2*facuRLdiff + 3.D0
     &    /2.D0*mt**3*s34*p5Dp34*den5x34**2*facuRLdiff + mt**5*den5x34*
     &    facuRLdiff*s345**(-1) )
      vert25x12 = vert25x12 + loopI2(s345,zip,mtsq,musq,ep) * (  - mt*
     &    p5Dp34*den5x34*facuRLdiff - mt*p5Dp34**3*den5x34**2*
     &    facuRLdiff - 1.D0/2.D0*mt*s34*p5Dp34**2*den5x34**2*facuRLdiff
     &     - mt**3*den5x34*facuRLdiff - mt**3*p5Dp34*den5x34*facuRLdiff
     &    *s345**(-1) - 3.D0*mt**3*p5Dp34**2*den5x34**2*facuRLdiff - 7.D
     &    0/2.D0*mt**3*s34*p5Dp34*den5x34**2*facuRLdiff - mt**3*s34**2*
     &    den5x34**2*facuRLdiff - mt**5*den5x34*facuRLdiff*s345**(-1) )
      vert25x12 = vert25x12 + loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &  * (  - 4.D0*mt*facuRLdiff - 4.D0*mt**3*p5Dp34*den5x34*
     &    facuRLdiff - 6.D0*mt**3*p5Dp34**3*den5x34**2*facuRLdiff - 2.D0
     &    *mt**3*s34*den5x34*facuRLdiff - 6.D0*mt**3*s34*p5Dp34**2*
     &    den5x34**2*facuRLdiff - 3.D0/2.D0*mt**3*s34**2*p5Dp34*
     &    den5x34**2*facuRLdiff )
      vert25x12 = vert25x12 + fp(ep)*mt*p5Dp34*den5x34*facuRLdiff + fp(
     &    ep)*mt**3*den5x34*facuRLdiff - fp(ep)*mt**3*p5Dp34*den5x34*
     &    facuRLdiff*s345**(-1) - fp(ep)*mt**5*den5x34*facuRLdiff*
     &    s345**(-1) - 2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)
     &    *mt**3*p5Dp34**3*den5x34**2*facuRLdiff - 2.D0*eploopI3(mtsq,s34
     &    ,s345,zip,mtsq,mtsq,musq,ep)*mt**3*s34*p5Dp34**2*den5x34**2*
     &    facuRLdiff - 1.D0/2.D0*eploopI3(mtsq,s34,s345,zip,mtsq,mtsq,
     &    musq,ep)*mt**3*s34**2*p5Dp34*den5x34**2*facuRLdiff - 2.D0*
     &    ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*p5Dp34**3*
     &    den5x34**2*facuRLdiff - 2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,
     &    mtsq,musq,ep)*mt**3*s34*p5Dp34**2*den5x34**2*facuRLdiff - 1.D0
     &    /2.D0*ep2loopI3(mtsq,s34,s345,zip,mtsq,mtsq,musq,ep)*mt**3*
     &    s34**2*p5Dp34*den5x34**2*facuRLdiff


      return
      end


