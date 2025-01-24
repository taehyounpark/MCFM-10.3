#include <ginac/ginac.h>
using namespace GiNaC;
#include <cln/cln.h>
#include <stdio.h>
#include <iostream>

typedef struct {double r,i;} complex_t;
typedef struct {complex_t c; signed char i0;} inum_t;
extern "C"{
complex_t geval_(inum_t * z, int* n);
};

complex_t geval_(inum_t * z, int* n) {
  cln::cl_inhibit_floating_point_underflow = true;
  lst w,s;
  for(long i=0;i<(*n)-1;i++)
  {
    ex zz;
    if (abs(z->c.i) < 1e-15)
      w.append((z->c.r));
    else
      w.append((z->c.r)+(z->c.i)*I);
    s.append(z->i0);
    z++;
  }
  try {
    ex ans = G(w,s,z->c.r).evalf();
    return {
      .r = ex_to<numeric>(evalf(real_part(ans))).to_double(),
      .i = ex_to<numeric>(evalf(imag_part(ans))).to_double()
    };
  } catch (...) {
    std::cout << "Caught!!\n";
    return {.r = 1.e11, .i = 1.e11};
  }

}
