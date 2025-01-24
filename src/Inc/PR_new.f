
!--- The variables R and P provide the Regular and Plus pieces associated
!--- with radiation from leg 1 (Q1(a,b,c,is) and leg 2 (Q2(a,b,c,is)
!--- In each case the parton labelling is Using the normal QM notation of putting
!--- everything backward
!---       emitted line after emission =   a
!---       emitter before emission     =    b
!---       spectator                   =    c
!--- There is no label for he or she who is emitted.
!--- Note that in general each piece will be composed of many different
!--- dipole contributions

      real(dp) :: Q1(-1:1,-1:1,-1:1,3),Q2(-1:1,-1:1,-1:1,3)
      common/RP_new/Q1,Q2
!$omp threadprivate(/RP_new/)
