      real(dp):: taucut,qtcut
      integer:: ntau
      logical:: usescet,abovecut,tauboost,dynamictau,incpowcorr,onlypowcorr,useQT_nnlo,useGLY
      common/mcfmtaucut/taucut
!$omp threadprivate(/mcfmtaucut/)
      common/mcfmqtcut/qtcut
!$omp threadprivate(/mcfmqtcut/)
      common/mcfmusescet/usescet,abovecut,tauboost,dynamictau,incpowcorr,onlypowcorr,useQT_nnlo,useGLY
      common/mcfmntau/ntau

