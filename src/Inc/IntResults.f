! Storing the integral results from Integralfill.f in the common block
      double complex Dint(dmax),Cint(cmax),Bint(bmax)
      common/IntResults/Dint,Cint,Bint
!$omp threadprivate(/IntResults/)
