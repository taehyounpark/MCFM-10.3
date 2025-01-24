!
!  SPDX-License-Identifier: GPL-3.0-or-later
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
c@article{DelDuca:2001fn,
c      author         = "Del Duca, V. and Kilgore, W. and Oleari, C. and Schmidt,
c                        C. and Zeppenfeld, D.",
c      title          = "{Gluon fusion contributions to H + 2 jet production}",
c      journal        = "Nucl. Phys.",
c      volume         = "B616",
c      year           = "2001",
c      pages          = "367-399",
c      doi            = "10.1016/S0550-3213(01)00446-1",
c      eprint         = "hep-ph/0108030",
c      archivePrefix  = "arXiv",
c      primaryClass   = "hep-ph",
c      reportNumber   = "MADPH-01-1235, BNL-HET-01-28, MSUHEP-10709, DFTT-19-2001",
c      SLACcitation   = "%%CITATION = HEP-PH/0108030;%%"
c}
c  Formula taken from DelDuca:2001fn Eqs.C7-C10 with modified normalization.
c     Implementation of Eq.~(11.4) and (11.5) from arXiv:2002.04018 v2
      use mod_qcdloop_c
c      use iso_c_binding
c      use iso_fortran_env
      implicit none
      complex(dp)::FL,FT,B0k1,B0k2,B0k12,C0k1xk2
      integer p1,p2,p3,p4,e
      real(dp)::mtsq
      real(dp)::k1Dk1,k2Dk2,k1Dk2,k1Dk12,k2Dk12,k12Dk12,Gram,msq,mu2
c      real(real128)::msq_128,mu2_128
      e=0
      msq=mtsq
      mu2=mtsq
      k1Dk1=s(p1,p2)
      k2Dk2=s(p3,p4)
      k1Dk2=0.5_dp*(s(p1,p3)+s(p1,p4)+s(p2,p3)+s(p2,p4))
      k1Dk12=k1Dk1+k1Dk2
      k2Dk12=k2Dk2+k1Dk2
      k12Dk12=k1Dk12+k2Dk12
      C0k1xk2=qlI3(k1Dk1,k12Dk12,k2Dk2,msq,msq,msq,mu2,0)
      B0k1=qlI2(k1Dk1,msq,msq,mu2,e)
      B0k2=qlI2(k2Dk2,msq,msq,mu2,e)
      B0k12=qlI2(k12Dk12,msq,msq,mu2,e)
c catch results that need to be computed in quad precision (?!)
c        msq_128 = real(msq,kind=real128); mu2_128=real(mu2,kind=real128)
c        B0k1=qlI2q(real(k1Dk1,kind=real128),msq_128,msq_128,mu2_128,0)
c        msq_128 = real(msq,kind=real128); mu2_128=real(mu2,kind=real128)
c        B0k2=qlI2q(real(k2Dk2,kind=real128),msq_128,msq_128,mu2_128,0)
c        msq_128 = real(msq,kind=real128); mu2_128=real(mu2,kind=real128)
c        B0k12=qlI2q(real(k12Dk12,kind=real128),msq_128,msq_128,mu2_128,0)

      Gram=k1Dk1*k2Dk2-k1Dk2**2
      FL=-((2._dp-3._dp*k1Dk1*k2Dk12/Gram)*(B0k1-B0k12)
     &    +(2._dp-3._dp*k2Dk2*k1Dk12/Gram)*(B0k2-B0k12)
     &  -(4._dp*mtsq+k1Dk1+k2Dk2+k12Dk12-3._dp*k1Dk1*k2Dk2*k12Dk12/Gram)
     &  *C0k1xk2-2._dp)/Gram;
      FT=-(k12Dk12*(B0k1+B0k2-2._dp*B0k12-2._dp*k1Dk2*C0k1xk2)
     & +(k1Dk1-k2Dk2)*(B0k1-B0k2))/Gram-k1Dk2*FL;
      return
