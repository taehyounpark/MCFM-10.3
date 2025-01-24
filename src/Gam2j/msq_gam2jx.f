!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      subroutine msq_gam2jqq(ja,ka,amp1_a,amp1_b,amp2_a,amp2_b,
     &                  msq0,msq1,msq2)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'ewcharge.f'
      integer:: ja,ka
      complex(dp):: a111,a112,a121,a211,a122,a212,a221,a222
      complex(dp):: b111,b112,b121,b211,b122,b212,b221,b222
      complex(dp):: amp1_a(2,2,2),amp1_b(2,2,2)
      complex(dp):: amp2_a(2,2,2),amp2_b(2,2,2)
      real(dp):: msq0,msq1,msq2

      if (ja  /=  ka) then
         a111=(Q(ja))*amp1_a(1,1,1)
     &       +(Q(ka))*amp1_b(1,1,1)
         a121=(Q(ja))*amp1_a(1,2,1)
     &       +(Q(ka))*amp1_b(1,2,1)
         a112=(Q(ja))*amp1_a(1,1,2)
     &       +(Q(ka))*amp1_b(1,1,2)
         a122=(Q(ja))*amp1_a(1,2,2)
     &       +(Q(ka))*amp1_b(1,2,2)
         a211=(Q(ja))*amp1_a(2,1,1)
     &       +(Q(ka))*amp1_b(2,1,1)
         a221=(Q(ja))*amp1_a(2,2,1)
     &       +(Q(ka))*amp1_b(2,2,1)
         a212=(Q(ja))*amp1_a(2,1,2)
     &       +(Q(ka))*amp1_b(2,1,2)
         a222=(Q(ja))*amp1_a(2,2,2)
     &       +(Q(ka))*amp1_b(2,2,2)

         msq0=zip
         msq1=
     &     real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &    +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &    +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &    +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp)
         msq2=zip

      elseif (ja == ka) then
         a111=(Q(ja))*(amp1_a(1,1,1)+amp1_b(1,1,1))
         b111=(Q(ja))*(amp2_a(1,1,1)+amp2_b(1,1,1))
         a112=(Q(ja))*(amp1_a(1,1,2)+amp1_b(1,1,2))
         b112=(Q(ja))*(amp2_a(1,1,2)+amp2_b(1,1,2))
         a221=(Q(ja))*(amp1_a(2,2,1)+amp1_b(2,2,1))
         b221=(Q(ja))*(amp2_a(2,2,1)+amp2_b(2,2,1))
         a222=(Q(ja))*(amp1_a(2,2,2)+amp1_b(2,2,2))
         b222=(Q(ja))*(amp2_a(2,2,2)+amp2_b(2,2,2))

         a121=(Q(ja))*amp1_a(1,2,1)
     &       +(Q(ka))*amp1_b(1,2,1)
         b121=(Q(ja))*amp2_a(1,2,1)
     &       +(Q(ka))*amp2_b(1,2,1)
         a122=(Q(ja))*amp1_a(1,2,2)
     &       +(Q(ka))*amp1_b(1,2,2)
         b122=(Q(ja))*amp2_a(1,2,2)
     &       +(Q(ka))*amp2_b(1,2,2)
         a211=(Q(ja))*amp1_a(2,1,1)
     &       +(Q(ka))*amp1_b(2,1,1)
         b211=(Q(ja))*amp2_a(2,1,1)
     &       +(Q(ka))*amp2_b(2,1,1)
         a212=(Q(ja))*amp1_a(2,1,2)
     &       +(Q(ka))*amp1_b(2,1,2)
         b212=(Q(ja))*amp2_a(2,1,2)
     &       +(Q(ka))*amp2_b(2,1,2)

         msq0=half*(
     &    +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &    +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
         msq1=half*(
     &     real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &    +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &    +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &    +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))
         msq2=half*(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))
      endif

      return
      end

      subroutine msq_gam2jqqb(ja,ka,amp1_a,amp1_b,amp2_a,amp2_b,
     &                  msq0,msq1,msq2,msq_up,msq_down)
      implicit none
      include 'types.f'
      include 'constants.f'
      include 'nf.f'
      include 'ewcharge.f'
      integer:: ja,ka
      complex(dp):: a111,a112,a121,a211,a122,a212,a221,a222
      complex(dp):: b111,b112,b121,b211,b122,b212,b221,b222
      complex(dp):: amp1_a(2,2,2),amp1_b(2,2,2)
      complex(dp):: amp2_a(2,2,2),amp2_b(2,2,2)
      real(dp):: msq0,msq1,msq2,msq_up,msq_down

      if (ja  /=  -ka) then
         a111=(Q(+ja))*amp1_a(1,1,1)
     &       +(Q(-ka))*amp1_b(1,1,1)
         a112=(Q(+ja))*amp1_a(1,1,2)
     &       +(Q(-ka))*amp1_b(1,1,2)
         a221=(Q(+ja))*amp1_a(2,2,1)
     &       +(Q(-ka))*amp1_b(2,2,1)
         a222=(Q(+ja))*amp1_a(2,2,2)
     &       +(Q(-ka))*amp1_b(2,2,2)

         a121=(Q(+ja))*amp1_a(1,2,1)
     &       +(Q(-ka))*amp1_b(1,2,1)
         a122=(Q(+ja))*amp1_a(1,2,2)
     &       +(Q(-ka))*amp1_b(1,2,2)
         a211=(Q(+ja))*amp1_a(2,1,1)
     &       +(Q(-ka))*amp1_b(2,1,1)
         a212=(Q(+ja))*amp1_a(2,1,2)
     &       +(Q(-ka))*amp1_b(2,1,2)

         msq0=zip
         msq1=zip
         msq2=
     &     real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &    +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &    +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &    +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp)

         msq_up=0._dp
         msq_down=0._dp

      elseif (ja == -ka) then
c--case where final state from annihilation diagrams is the same quark
         a111=(Q(ja))*(amp1_a(1,1,1)+amp1_b(1,1,1))
         b111=(Q(ja))*(amp2_a(1,1,1)+amp2_b(1,1,1))

         a112=(Q(ja))*(amp1_a(1,1,2)+amp1_b(1,1,2))
         b112=(Q(ja))*(amp2_a(1,1,2)+amp2_b(1,1,2))

         a221=(Q(ja))*(amp1_a(2,2,1)+amp1_b(2,2,1))
         b221=(Q(ja))*(amp2_a(2,2,1)+amp2_b(2,2,1))

         a222=(Q(ja))*(amp1_a(2,2,2)+amp1_b(2,2,2))
         b222=(Q(ja))*(amp2_a(2,2,2)+amp2_b(2,2,2))

         a121=(Q(+ja))*amp1_a(1,2,1)
     &       +(Q(-ka))*amp1_b(1,2,1)
         a122=(Q(+ja))*amp1_a(1,2,2)
     &       +(Q(-ka))*amp1_b(1,2,2)
         a211=(Q(+ja))*amp1_a(2,1,1)
     &       +(Q(-ka))*amp1_b(2,1,1)
         a212=(Q(+ja))*amp1_a(2,1,2)
     &       +(Q(-ka))*amp1_b(2,1,2)

         b121=(Q(+ja))*amp2_a(1,2,1)
     &       +(Q(-ka))*amp2_b(1,2,1)
         b122=(Q(+ja))*amp2_a(1,2,2)
     &       +(Q(-ka))*amp2_b(1,2,2)
         b211=(Q(+ja))*amp2_a(2,1,1)
     &       +(Q(-ka))*amp2_b(2,1,1)
         b212=(Q(+ja))*amp2_a(2,1,2)
     &       +(Q(-ka))*amp2_b(2,1,2)

         msq0=(
     &    +real(a111*conjg(b111),dp)+real(a112*conjg(b112),dp)
     &    +real(a221*conjg(b221),dp)+real(a222*conjg(b222),dp))*two/xn
         msq1=(
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp))
         msq2=(
     &     real(a111*conjg(a111),dp)+real(a112*conjg(a112),dp)
     &    +real(a221*conjg(a221),dp)+real(a222*conjg(a222),dp)
     &    +real(a122*conjg(a122),dp)+real(a212*conjg(a212),dp)
     &    +real(a121*conjg(a121),dp)+real(a211*conjg(a211),dp))

         b111=(Q(+ja))*amp2_a(1,1,1)
     &       +(Q(+1) )*amp2_b(1,1,1)
         b112=(Q(+ja))*amp2_a(1,1,2)
     &       +(Q(+1) )*amp2_b(1,1,2)
         b221=(Q(+ja))*amp2_a(2,2,1)
     &       +(Q(+1) )*amp2_b(2,2,1)
         b222=(Q(+ja))*amp2_a(2,2,2)
     &       +(Q(+1) )*amp2_b(2,2,2)
         b121=(Q(+ja))*amp2_a(1,2,1)
     &       +(Q(+1) )*amp2_b(1,2,1)
         b122=(Q(+ja))*amp2_a(1,2,2)
     &       +(Q(+1) )*amp2_b(1,2,2)
         b211=(Q(+ja))*amp2_a(2,1,1)
     &       +(Q(+1) )*amp2_b(2,1,1)
         b212=(Q(+ja))*amp2_a(2,1,2)
     &       +(Q(+1) )*amp2_b(2,1,2)

         msq_down=
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp)

         b111=(Q(+ja))*amp2_a(1,1,1)
     &       +(Q(+2) )*amp2_b(1,1,1)
         b112=(Q(+ja))*amp2_a(1,1,2)
     &       +(Q(+2) )*amp2_b(1,1,2)
         b221=(Q(+ja))*amp2_a(2,2,1)
     &       +(Q(+2) )*amp2_b(2,2,1)
         b222=(Q(+ja))*amp2_a(2,2,2)
     &       +(Q(+2) )*amp2_b(2,2,2)
         b121=(Q(+ja))*amp2_a(1,2,1)
     &       +(Q(+2) )*amp2_b(1,2,1)
         b122=(Q(+ja))*amp2_a(1,2,2)
     &       +(Q(+2) )*amp2_b(1,2,2)
         b211=(Q(+ja))*amp2_a(2,1,1)
     &       +(Q(+2) )*amp2_b(2,1,1)
         b212=(Q(+ja))*amp2_a(2,1,2)
     &       +(Q(+2) )*amp2_b(2,1,2)

         msq_up=
     &     real(b111*conjg(b111),dp)+real(b112*conjg(b112),dp)
     &    +real(b221*conjg(b221),dp)+real(b222*conjg(b222),dp)
     &    +real(b122*conjg(b122),dp)+real(b212*conjg(b212),dp)
     &    +real(b121*conjg(b121),dp)+real(b211*conjg(b211),dp)

      endif

      return
      end
