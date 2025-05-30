!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
      function q4ghpmp1(p1,p2,p3,p4,p5,za,zb)

c     This is the reduced matrix element squared
c     for the process
c     q(p1+)+qbar(p2) --> H((p5+p6)+Q(p3-)+qbar(p4)+g(p5+)
      include 'types.f'
      complex(dp)::q4ghpmp1

      include 'constants.f'
      include 'mxpart.f'
      include 'sprods_com.f'
      include 'zprods_decl.f'
      integer::p1,p2,p3,p4,p5
      real(dp)::s5h,t123,t234,t345,t125
      integer::i1,i2,i3,i4
      complex(dp)::zab,zaba

      zab(i1,i2,i3)=+za(i1,i2)*zb(i2,i3)
      zaba(i1,i2,i3,i4)=+za(i1,i2)*zb(i2,i3)*za(i3,i4)


      t123=s(p1,p2)+s(p2,p3)+s(p3,p1)
      t234=s(p2,p3)+s(p3,p4)+s(p4,p2)
      t125=s(p1,p2)+s(p2,p5)+s(p5,p1)
      t345=s(p3,p4)+s(p4,p5)+s(p5,p3)
      s5h=s(p1,p2)+s(p1,p3)+s(p1,p4)
     &            +s(p2,p3)+s(p2,p4)
     &                     +s(p3,p4)


c %\cite{DelDuca:2004wt}
c \bibitem{DelDuca:2004wt}
c V.~Del Duca, A.~Frizzo and F.~Maltoni,
c %``Higgs boson production in association with three jets,''
c JHEP {\bf 0405}, 064 (2004)
c [arXiv:hep-ph/0404013].
c %%CITATION = HEP-PH 0404013;%%
c---B13
      q4ghpmp1=
     & +one/(za(p1,p5)*s(p1,p2)*s(p3,p4)*t345)
     & *(za(p2,p3)*zb(p1,p2)*(-za(p2,p3)*zb(p4,p5)
     & *(zab(p1,p4,p3)+zab(p1,p5,p3))
     & +za(p1,p3)*zb(p3,p5)*(zab(p2,p3,p4)+zab(p2,p5,p4)))
     & +za(p1,p2)*(zab(p3,p4,p1)+zab(p3,p5,p1))
     & *(-za(p1,p3)*zb(p1,p4)*zb(p3,p5)+zb(p4,p5)*(s(p1,p4)+s(p1,p5)))
     & -s(p1,p2)*s(p4,p5)*za(p2,p3)*zb(p4,p5))
      q4ghpmp1=q4ghpmp1
     & -one/(s(p1,p2)*s(p3,p4)*s5h)
     & *(-za(p1,p2)*zb(p1,p5)**2*(zab(p3,p1,p4)+zab(p3,p2,p4))
     & -zb(p1,p5)*zb(p4,p5)*(two*zaba(p2,p1,p4,p3)
     & -za(p2,p3)*(s(p1,p3)+s(p1,p4)+s(p2,p3)+s(p2,p4)))
     & -za(p3,p4)*zb(p4,p5)**2*(zab(p2,p3,p1)+zab(p2,p4,p1)))
      q4ghpmp1=q4ghpmp1
     & +one/(za(p1,p5)*s(p1,p2)*s(p3,p4)*t125)
     & *(-za(p1,p2)**2*za(p3,p4)*zb(p1,p4)**2*zb(p2,p5)
     & +za(p1,p2)*(-za(p2,p3)**2*zb(p1,p2)*zb(p2,p5)*zb(p3,p4)
     & +za(p3,p4)*zb(p1,p4)*zb(p4,p5)*(s(p2,p5)-s(p1,p5)))
     & +zb(p1,p5)*(-za(p1,p3)*za(p2,p3)*zb(p3,p4)*(s(p1,p5)+s(p2,p5))
     & +za(p1,p5)*za(p2,p5)*za(p3,p4)*zb(p4,p5)**2))
      q4ghpmp1=q4ghpmp1
     & +one/(za(p1,p5)*s(p1,p2)*s(p3,p4))
     & *(za(p2,p3)*zb(p1,p5)*(zab(p1,p3,p4)-zab(p1,p5,p4))
     & +za(p1,p2)*(za(p1,p3)*zb(p1,p4)*zb(p1,p5)-zb(p4,p5)
     & *(zab(p3,p2,p1)-zab(p3,p4,p1))))
       q4ghpmp1=q4ghpmp1
     & -one/(za(p1,p5)*za(p4,p5)*s(p1,p2)*t345)
     & *(-za(p2,p3)**2*zb(p1,p2)*(zab(p1,p4,p3)+zab(p1,p5,p3))
     & +za(p1,p2)*(s(p1,p4)+s(p1,p5))*(zab(p3,p4,p1)+zab(p3,p5,p1))
     & -s(p1,p2)*s(p4,p5)*za(p2,p3))
       q4ghpmp1=q4ghpmp1
     & +za(p2,p3)*zb(p4,p5)/(s(p1,p2)*t123*s5h)
     & *(-zb(p1,p5)*t123-zb(p4,p5)*(zab(p4,p2,p1)+zab(p4,p3,p1)))
     & +za(p2,p3)*zb(p1,p5)/(s(p3,p4)*t234*s5h)
     & *(-zb(p1,p5)*(zab(p1,p2,p4)+zab(p1,p3,p4))-zb(p4,p5)*t234)
       q4ghpmp1=q4ghpmp1
     & +one/(za(p1,p5)*s(p3,p4)*t125)
     & *(za(p2,p3)*zb(p3,p4)*(zab(p3,p1,p5)+zab(p3,p2,p5))
     & +za(p3,p4)*zb(p4,p5)*(zab(p2,p1,p4)+zab(p2,p5,p4)))


      return
      end
