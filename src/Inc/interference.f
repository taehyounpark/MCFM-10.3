
      logical:: interference,bw34_56,closestZ
      common/interference/interference,bw34_56,closestZ
      real(dp):: vsymfact,moppmin
      common/vsymfact/vsymfact,moppmin
!$omp threadprivate(/interference/,/vsymfact/)
