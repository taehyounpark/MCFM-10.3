!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!

      function Finite1x1qqgagau(s,t)
c     Results taken from hep-ph/0201274
c \bibitem{Anastasiou:2002zn}
c  C.~Anastasiou, E.~W.~N.~Glover and M.~E.~Tejeda-Yeomans,
c  %``Two loop QED and QCD corrections to massless fermion boson scattering,''
c  Nucl.\ Phys.\ B {\bf 629}, 255 (2002)
c  [hep-ph/0201274].
      implicit none
      include 'types.f'
      real(dp)::Finite1x1qqgagau
      include 'constants.f'
      real(dp)::s,t,u,x,y,Lx,Ly,AGTYG1u
      u=-s-t
      x=-t/s
      y=-u/s

      Lx=log(x)
      Ly=log(y)

      Finite1x1qqgagau=xn*CF**2*AGTYG1u(s,t,Lx,Ly)
      return
      end





