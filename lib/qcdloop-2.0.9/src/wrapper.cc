//
// QCDLoop 2016
//
// Authors: Stefano Carrazza: stefano.carrazza@cern.ch
//          Keith Ellis: keith.ellis@durham.ac.uk
//          Giulia Zanderighi: giulia.zanderighi@cern.ch
//
// OpenMP support: Tobias Neumann

#include "qcdloop/wrapper.h"
#include "qcdloop/qcdloop.h"
#include <stdexcept>
#include <iostream>
#include <omp.h>

extern"C" {

// result container
extern vector<complex> r;
#ifdef __INTEL_COMPILER
vector<complex> r(3);
#pragma omp threadprivate(r)
#else
#pragma omp threadprivate(r)
vector<complex> r(3);
#endif

extern vector<qcomplex> rq;
#ifdef __INTEL_COMPILER
vector<qcomplex> rq(3);
#pragma omp threadprivate(rq)
#else
#pragma omp threadprivate(rq)
vector<qcomplex> rq(3);
#endif

// topologies
extern TadPole<complex,double,double> td;
#ifdef __INTEL_COMPILER
TadPole<complex,double,double> td;
#pragma omp threadprivate(td)
#else
#pragma omp threadprivate(td)
TadPole<complex,double,double> td;
#endif

extern TadPole<complex,complex,double> tdc;
#ifdef __INTEL_COMPILER
TadPole<complex,complex,double> tdc;
#pragma omp threadprivate(tdc)
#else
#pragma omp threadprivate(tdc)
TadPole<complex,complex,double> tdc;
#endif


extern TadPole<qcomplex,qdouble,qdouble> tdq;
#ifdef __INTEL_COMPILER
TadPole<qcomplex,qdouble,qdouble> tdq;
#pragma omp threadprivate(tdq)
#else
#pragma omp threadprivate(tdq)
TadPole<qcomplex,qdouble,qdouble> tdq;
#endif

extern TadPole<qcomplex,qcomplex,qdouble> tdcq;
#ifdef __INTEL_COMPILER
TadPole<qcomplex,qcomplex,qdouble> tdcq;
#pragma omp threadprivate(tdcq)
#else
#pragma omp threadprivate(tdcq)
TadPole<qcomplex,qcomplex,qdouble> tdcq;
#endif

extern Bubble<complex,double,double> bb;
#ifdef __INTEL_COMPILER
Bubble<complex,double,double> bb;
#pragma omp threadprivate(bb)
#else
#pragma omp threadprivate(bb)
Bubble<complex,double,double> bb;
#endif

extern Bubble<complex,complex,double> bbc;
#ifdef __INTEL_COMPILER
Bubble<complex,complex,double> bbc;
#pragma omp threadprivate(bbc)
#else
#pragma omp threadprivate(bbc)
Bubble<complex,complex,double> bbc;
#endif

extern Bubble<qcomplex,qdouble,qdouble> bbq;
#ifdef __INTEL_COMPILER
Bubble<qcomplex,qdouble,qdouble> bbq;
#pragma omp threadprivate(bbq)
#else
#pragma omp threadprivate(bbq)
Bubble<qcomplex,qdouble,qdouble> bbq;
#endif

extern Bubble<qcomplex,qcomplex,qdouble> bbcq;
#ifdef __INTEL_COMPILER
Bubble<qcomplex,qcomplex,qdouble> bbcq;
#pragma omp threadprivate(bbcq)
#else
#pragma omp threadprivate(bbcq)
Bubble<qcomplex,qcomplex,qdouble> bbcq;
#endif

extern Triangle<complex,double,double> tr;
#ifdef __INTEL_COMPILER
Triangle<complex,double,double> tr;
#pragma omp threadprivate(tr)
#else
#pragma omp threadprivate(tr)
Triangle<complex,double,double> tr;
#endif

extern Triangle<complex,complex,double> trc;
#ifdef __INTEL_COMPILER
Triangle<complex,complex,double> trc;
#pragma omp threadprivate(trc)
#else
#pragma omp threadprivate(trc)
Triangle<complex,complex,double> trc;
#endif

extern Triangle<qcomplex,qdouble,qdouble> trq;
#ifdef __INTEL_COMPILER
Triangle<qcomplex,qdouble,qdouble> trq;
#pragma omp threadprivate(trq)
#else
#pragma omp threadprivate(trq)
Triangle<qcomplex,qdouble,qdouble> trq;
#endif

extern Triangle<qcomplex,qcomplex,qdouble> trcq;
#ifdef __INTEL_COMPILER
Triangle<qcomplex,qcomplex,qdouble> trcq;
#pragma omp threadprivate(trcq)
#else
#pragma omp threadprivate(trcq)
Triangle<qcomplex,qcomplex,qdouble> trcq;
#endif

extern Box<complex,double,double> bo;
#ifdef __INTEL_COMPILER
Box<complex,double,double> bo;
#pragma omp threadprivate(bo)
#else
#pragma omp threadprivate(bo)
Box<complex,double,double> bo;
#endif

extern Box<complex,complex,double> boc;
#ifdef __INTEL_COMPILER
Box<complex,complex,double> boc;
#pragma omp threadprivate(boc)
#else
#pragma omp threadprivate(boc)
Box<complex,complex,double> boc;
#endif

extern Box<qcomplex,qdouble,qdouble> boq;
#ifdef __INTEL_COMPILER
Box<qcomplex,qdouble,qdouble> boq;
#pragma omp threadprivate(boq)
#else
#pragma omp threadprivate(boq)
Box<qcomplex,qdouble,qdouble> boq;
#endif

extern Box<qcomplex,qcomplex,qdouble> bocq;
#ifdef __INTEL_COMPILER
Box<qcomplex,qcomplex,qdouble> bocq;
#pragma omp threadprivate(bocq)
#else
#pragma omp threadprivate(bocq)
Box<qcomplex,qcomplex,qdouble> bocq;
#endif

// global vectors for improving parsing speed
extern vector<double>   mI1;
#ifdef __INTEL_COMPILER
vector<double>   mI1(1);
#pragma omp threadprivate(mI1)
#else
#pragma omp threadprivate(mI1)
vector<double>   mI1(1);
#endif

extern vector<complex>  mI1c;
#ifdef __INTEL_COMPILER
vector<complex>  mI1c(1);
#pragma omp threadprivate(mI1c)
#else
#pragma omp threadprivate(mI1c)
vector<complex>  mI1c(1);
#endif

extern vector<qdouble>  mI1q;
#ifdef __INTEL_COMPILER
vector<qdouble>  mI1q(1);
#pragma omp threadprivate(mI1q)
#else
#pragma omp threadprivate(mI1q)
vector<qdouble>  mI1q(1);
#endif

extern vector<qcomplex> mI1cq;
#ifdef __INTEL_COMPILER
vector<qcomplex> mI1cq(1);
#pragma omp threadprivate(mI1cq)
#else
#pragma omp threadprivate(mI1cq)
vector<qcomplex> mI1cq(1);
#endif

extern vector<double>   mI2;
#ifdef __INTEL_COMPILER
vector<double>   mI2(2);
#pragma omp threadprivate(mI2)
#else
#pragma omp threadprivate(mI2)
vector<double>   mI2(2);
#endif

extern vector<complex>  mI2c;
#ifdef __INTEL_COMPILER
vector<complex>  mI2c(2);
#pragma omp threadprivate(mI2c)
#else
#pragma omp threadprivate(mI2c)
vector<complex>  mI2c(2);
#endif

extern vector<qdouble>  mI2q;
#ifdef __INTEL_COMPILER
vector<qdouble>  mI2q(2);
#pragma omp threadprivate(mI2q)
#else
#pragma omp threadprivate(mI2q)
vector<qdouble>  mI2q(2);
#endif

extern vector<qcomplex> mI2cq;
#ifdef __INTEL_COMPILER
vector<qcomplex> mI2cq(2);
#pragma omp threadprivate(mI2cq)
#else
#pragma omp threadprivate(mI2cq)
vector<qcomplex> mI2cq(2);
#endif

extern vector<double>   mI3;
#ifdef __INTEL_COMPILER
vector<double>   mI3(3);
#pragma omp threadprivate(mI3)
#else
#pragma omp threadprivate(mI3)
vector<double>   mI3(3);
#endif

extern vector<complex>  mI3c;
#ifdef __INTEL_COMPILER
vector<complex>  mI3c(3);
#pragma omp threadprivate(mI3c)
#else
#pragma omp threadprivate(mI3c)
vector<complex>  mI3c(3);
#endif

extern vector<qdouble>  mI3q;
#ifdef __INTEL_COMPILER
vector<qdouble>  mI3q(3);
#pragma omp threadprivate(mI3q)
#else
#pragma omp threadprivate(mI3q)
vector<qdouble>  mI3q(3);
#endif

extern vector<qcomplex> mI3cq;
#ifdef __INTEL_COMPILER
vector<qcomplex> mI3cq(3);
#pragma omp threadprivate(mI3cq)
#else
#pragma omp threadprivate(mI3cq)
vector<qcomplex> mI3cq(3);
#endif

extern vector<double>   mI4;
#ifdef __INTEL_COMPILER
vector<double>   mI4(4);
#pragma omp threadprivate(mI4)
#else
#pragma omp threadprivate(mI4)
vector<double>   mI4(4);
#endif

extern vector<complex>  mI4c;
#ifdef __INTEL_COMPILER
vector<complex>  mI4c(4);
#pragma omp threadprivate(mI4c)
#else
#pragma omp threadprivate(mI4c)
vector<complex>  mI4c(4);
#endif

extern vector<qdouble>  mI4q;
#ifdef __INTEL_COMPILER
vector<qdouble>  mI4q(4);
#pragma omp threadprivate(mI4q)
#else
#pragma omp threadprivate(mI4q)
vector<qdouble>  mI4q(4);
#endif

extern vector<qcomplex> mI4cq;
#ifdef __INTEL_COMPILER
vector<qcomplex> mI4cq(4);
#pragma omp threadprivate(mI4cq)
#else
#pragma omp threadprivate(mI4cq)
vector<qcomplex> mI4cq(4);
#endif

extern vector<double>   pI2;
#ifdef __INTEL_COMPILER
vector<double>   pI2(1);
#pragma omp threadprivate(pI2)
#else
#pragma omp threadprivate(pI2)
vector<double>   pI2(1);
#endif

extern vector<qdouble>  pI2q;
#ifdef __INTEL_COMPILER
vector<qdouble>  pI2q(1);
#pragma omp threadprivate(pI2q)
#else
#pragma omp threadprivate(pI2q)
vector<qdouble>  pI2q(1);
#endif

extern vector<double>   pI3;
#ifdef __INTEL_COMPILER
vector<double>   pI3(3);
#pragma omp threadprivate(pI3)
#else
#pragma omp threadprivate(pI3)
vector<double>   pI3(3);
#endif

extern vector<qdouble>  pI3q;
#ifdef __INTEL_COMPILER
vector<qdouble>  pI3q(3);
#pragma omp threadprivate(pI3q)
#else
#pragma omp threadprivate(pI3q)
vector<qdouble>  pI3q(3);
#endif

extern vector<double>   pI4;
#ifdef __INTEL_COMPILER
vector<double>   pI4(6);
#pragma omp threadprivate(pI4)
#else
#pragma omp threadprivate(pI4)
vector<double>   pI4(6);
#endif

extern vector<qdouble>  pI4q;
#ifdef __INTEL_COMPILER
vector<qdouble>  pI4q(6);
#pragma omp threadprivate(pI4q)
#else
#pragma omp threadprivate(pI4q)
vector<qdouble>  pI4q(6);
#endif

void qlcachesize(const int &size)
{
  td.setCacheSize(size);
  tdc.setCacheSize(size);
  tdq.setCacheSize(size);
  tdcq.setCacheSize(size);

  bb.setCacheSize(size);
  bbc.setCacheSize(size);
  bbq.setCacheSize(size);
  bbcq.setCacheSize(size);

  tr.setCacheSize(size);
  trc.setCacheSize(size);
  trq.setCacheSize(size);
  trcq.setCacheSize(size);

  bo.setCacheSize(size);
  boc.setCacheSize(size);
  boq.setCacheSize(size);
  bocq.setCacheSize(size);
}

bool qlzero(double const& x)
{
  return td.iszero(x);
}

bool qlzeroq(qdouble const& x)
{
  return tdq.iszero(x);
}

bool qlnonzero(double const& x)
{
  return !td.iszero(x);
}

bool qlnonzeroq(qdouble const& x)
{
  return !tdq.iszero(x);
}

complex cln(complex const& x, double const& isig)
{
  return td.cLn(x, isig);
}

qcomplex clnq(qcomplex const& x, qdouble const& isig)
{
  return tdq.cLn(x, isig);
}

void qltadpole(complex (&out)[3], double const& mu2, double const& m1)
{
  try
  {
    mI1[0] = m1;
    td.integral(r, mu2, mI1);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltadpolec(complex (&out)[3], double const& mu2, complex const& m1)
{
  try
  {
    mI1c[0] = m1;
    tdc.integral(r, mu2, mI1c);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltadpoleq(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1)
{
  try
  {
    mI1q[0] = m1;
    tdq.integral(rq, mu2, mI1q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltadpolecq(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1)
{
  try
  {
    mI1cq[0] = m1;
    tdcq.integral(rq, mu2, mI1cq);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubble(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& p1)
{
  try
  {
    mI2[0] = m1;
    mI2[1] = m2;
    pI2[0] = p1;
    bb.integral(r, mu2, mI2, pI2);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubblec(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, double const& p1)
{
  try
  {
    mI2c[0] = m1;
    mI2c[1] = m2;
    pI2[0]  = p1;
    bbc.integral(r, mu2, mI2c, pI2);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubbleq(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& p1)
{
  try
  {
    mI2q[0] = m1;
    mI2q[1] = m2;
    pI2q[0]  = p1;
    bbq.integral(rq, mu2, mI2q, pI2q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbubblecq(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qdouble const& p1)
{
  try
  {
    mI2cq[0] = m1;
    mI2cq[1] = m2;
    pI2q[0]  = p1;
    bbcq.integral(rq, mu2, mI2cq, pI2q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltriangle(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& m3, double const& p1, double const& p2, double const& p3)
{
  try
  {
    mI3[0] = m1;
    mI3[1] = m2;
    mI3[2] = m3;
    pI3[0] = p1;
    pI3[1] = p2;
    pI3[2] = p3;
    tr.integral(r, mu2, mI3, pI3);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltrianglec(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, complex const& m3, double const& p1, double const& p2, double const& p3)
{
  try
  {
    mI3c[0] = m1;
    mI3c[1] = m2;
    mI3c[2] = m3;
    pI3[0] = p1;
    pI3[1] = p2;
    pI3[2] = p3;
    trc.integral(r, mu2, mI3c, pI3);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qltriangleq(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& p1, qdouble const& p2, qdouble const& p3)
{
  try
  {
    mI3q[0] = m1;
    mI3q[1] = m2;
    mI3q[2] = m3;
    pI3q[0] = p1;
    pI3q[1] = p2;
    pI3q[2] = p3;
    trq.integral(rq, mu2, mI3q, pI3q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}
void qltrianglecq(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qdouble const& p1, qdouble const& p2, qdouble const& p3)
{
  try
  {
    mI3cq[0] = m1;
    mI3cq[1] = m2;
    mI3cq[2] = m3;
    pI3q[0] = p1;
    pI3q[1] = p2;
    pI3q[2] = p3;
    trcq.integral(rq, mu2, mI3cq, pI3q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlbox(complex (&out)[3], double const& mu2, double const& m1, double const& m2, double const& m3, double const& m4, double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23)
{
  try
  {
    mI4[0] = m1;
    mI4[1] = m2;
    mI4[2] = m3;
    mI4[3] = m4;
    pI4[0] = p1;
    pI4[1] = p2;
    pI4[2] = p3;
    pI4[3] = p4;
    pI4[4] = s12;
    pI4[5] = s23;
    bo.integral(r, mu2, mI4, pI4);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlboxc(complex (&out)[3], double const& mu2, complex const& m1, complex const& m2, complex const& m3, complex const& m4, double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23)
{
  try
  {
    mI4c[0] = m1;
    mI4c[1] = m2;
    mI4c[2] = m3;
    mI4c[3] = m4;
    pI4[0] = p1;
    pI4[1] = p2;
    pI4[2] = p3;
    pI4[3] = p4;
    pI4[4] = s12;
    pI4[5] = s23;
    boc.integral(r, mu2, mI4c, pI4);
    out[0] = r[0];
    out[1] = r[1];
    out[2] = r[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlboxq(qcomplex (&out)[3], qdouble const& mu2, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& m4, qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23)
{
  try
  {
    mI4q[0] = m1;
    mI4q[1] = m2;
    mI4q[2] = m3;
    mI4q[3] = m4;
    pI4q[0] = p1;
    pI4q[1] = p2;
    pI4q[2] = p3;
    pI4q[3] = p4;
    pI4q[4] = s12;
    pI4q[5] = s23;
    boq.integral(rq, mu2, mI4q, pI4q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

void qlboxcq(qcomplex (&out)[3], qdouble const& mu2, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qcomplex const& m4, qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23)
{
  try
  {
    mI4cq[0] = m1;
    mI4cq[1] = m2;
    mI4cq[2] = m3;
    mI4cq[3] = m4;
    pI4q[0] = p1;
    pI4q[1] = p2;
    pI4q[2] = p3;
    pI4q[3] = p4;
    pI4q[4] = s12;
    pI4q[5] = s23;
    bocq.integral(rq, mu2, mI4cq, pI4q);
    out[0] = rq[0];
    out[1] = rq[1];
    out[2] = rq[2];
  }
  catch (std::exception &e)
  {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
}

#ifdef QL_NAMES

  // for backward compatibility.
  void qlinit()
  {
    std::cout << ql::yellow << "[QCDLoop warning]: this wrapper is not thread-safe." << std::endl;
    std::cout << "[QCDLoop suggestion]: consider developing object-oriented code." << ql::def << std::endl;
  }

  complex qli1(double const& m1, double const& mu2, int const& ep)
  {
    try
    {
      mI1[0] = m1;
      td.integral(r, mu2, mI1);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  complex qli1c(complex const& m1, double const& mu2, int const& ep)
  {
    try
    {
      mI1c[0] = m1;
      tdc.integral(r, mu2, mI1c);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli1q(qdouble const& m1, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI1q[0] = m1;
      tdq.integral(rq, mu2, mI1q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli1qc(qcomplex const& m1, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI1cq[0] = m1;
      tdcq.integral(rq, mu2, mI1cq);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli2(double const& p1, double const& m1, double const& m2, double const& mu2, int const& ep)
  {
    try
    {
      mI2[0] = m1;
      mI2[1] = m2;
      pI2[0] = p1;
      bb.integral(r, mu2, mI2, pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }    
    return r[abs(ep)];
  }

  complex qli2c(double const& p1, complex const& m1, complex const& m2, double const& mu2, int const& ep)
  {
    try
    {
      mI2c[0] = m1;
      mI2c[1] = m2;
      pI2[0] = p1;
      bbc.integral(r, mu2, mI2c, pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli2q(qdouble const& p1, qdouble const& m1, qdouble const& m2, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI2q[0] = m1;
      mI2q[1] = m2;
      pI2q[0] = p1;
      bbq.integral(rq, mu2, mI2q, pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli2qc(qdouble const& p1, qcomplex const& m1, qcomplex const& m2, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI2cq[0] = m1;
      mI2cq[1] = m2;
      pI2q[0] = p1;
      bbcq.integral(rq, mu2, mI2cq, pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli2p(double const& p1, double const& m1, double const& m2, double const& mu2, int const& ep)
  {
    try
    {
      mI2[0] = m1;
      mI2[1] = m2;
      pI2[0] = p1;
      bb.derivative(r, mu2, mI2, pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  complex qli2pc(double const& p1, complex const& m1, complex const& m2, double const& mu2, int const& ep)
  {
    try
    {
      mI2c[0] = m1;
      mI2c[1] = m2;
      pI2[0] = p1;
      bbc.derivative(r, mu2, mI2c, pI2);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli2pq(qdouble const& p1, qdouble const& m1, qdouble const& m2, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI2q[0] = m1;
      mI2q[1] = m2;
      pI2q[0] = p1;
      bbq.derivative(rq, mu2, mI2q, pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli2pqc(qdouble const& p1, qcomplex const& m1, qcomplex const& m2, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI2cq[0] = m1;
      mI2cq[1] = m2;
      pI2q[0] = p1;
      bbcq.derivative(rq, mu2, mI2cq, pI2q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli3(double const& p1, double const& p2, double const& p3, double const& m1, double const& m2, double const& m3, double const& mu2, int const& ep)
  {
    try
    {
      mI3[0] = m1;
      mI3[1] = m2;
      mI3[2] = m3;
      pI3[0] = p1;
      pI3[1] = p2;
      pI3[2] = p3;
      tr.integral(r, mu2, mI3, pI3);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  complex qli3c(double const& p1, double const& p2, double const& p3, complex const& m1, complex const& m2, complex const& m3, double const& mu2, int const& ep)
  {
    try
    {
      mI3c[0] = m1;
      mI3c[1] = m2;
      mI3c[2] = m3;
      pI3[0] = p1;
      pI3[1] = p2;
      pI3[2] = p3;
      trc.integral(r, mu2, mI3c, pI3);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli3q(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI3q[0] = m1;
      mI3q[1] = m2;
      mI3q[2] = m3;
      pI3q[0] = p1;
      pI3q[1] = p2;
      pI3q[2] = p3;
      trq.integral(rq, mu2, mI3q, pI3q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli3qc(qdouble const& p1, qdouble const& p2, qdouble const& p3, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI3cq[0] = m1;
      mI3cq[1] = m2;
      mI3cq[2] = m3;
      pI3q[0] = p1;
      pI3q[1] = p2;
      pI3q[2] = p3;
      trcq.integral(rq, mu2, mI3cq, pI3q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  complex qli4(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, double const& m1, double const& m2, double const& m3, double const& m4, double const& mu2, int const& ep)
  {
    try
    {
      mI4[0] = m1;
      mI4[1] = m2;
      mI4[2] = m3;
      mI4[3] = m4;
      pI4[0] = p1;
      pI4[1] = p2;
      pI4[2] = p3;
      pI4[3] = p4;
      pI4[4] = s12;
      pI4[5] = s23;
      bo.integral(r, mu2, mI4, pI4);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }

    return r[abs(ep)];
  }

  complex qli4c(double const& p1, double const& p2, double const& p3, double const& p4, double const& s12, double const& s23, complex const& m1, complex const& m2, complex const& m3, complex const& m4, double const& mu2, int const& ep)
  {
    try
    {
      mI4c[0] = m1;
      mI4c[1] = m2;
      mI4c[2] = m3;
      mI4c[3] = m4;
      pI4[0] = p1;
      pI4[1] = p2;
      pI4[2] = p3;
      pI4[3] = p4;
      pI4[4] = s12;
      pI4[5] = s23;
      boc.integral(r, mu2, mI4c, pI4);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return r[abs(ep)];
  }

  qcomplex qli4q(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qdouble const& m1, qdouble const& m2, qdouble const& m3, qdouble const& m4, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI4q[0] = m1;
      mI4q[1] = m2;
      mI4q[2] = m3;
      mI4q[3] = m4;
      pI4q[0] = p1;
      pI4q[1] = p2;
      pI4q[2] = p3;
      pI4q[3] = p4;
      pI4q[4] = s12;
      pI4q[5] = s23;
      boq.integral(rq, mu2, mI4q, pI4q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

  qcomplex qli4qc(qdouble const& p1, qdouble const& p2, qdouble const& p3, qdouble const& p4, qdouble const& s12, qdouble const& s23, qcomplex const& m1, qcomplex const& m2, qcomplex const& m3, qcomplex const& m4, qdouble const& mu2, int const& ep)
  {
    try
    {
      mI4cq[0] = m1;
      mI4cq[1] = m2;
      mI4cq[2] = m3;
      mI4cq[3] = m4;
      pI4q[0] = p1;
      pI4q[1] = p2;
      pI4q[2] = p3;
      pI4q[3] = p4;
      pI4q[4] = s12;
      pI4q[5] = s23;
      bocq.integral(rq, mu2, mI4cq, pI4q);
    }
    catch (std::exception &e)
    {
      std::cout << e.what() << std::endl;
      exit(-1);
    }
    return rq[abs(ep)];
  }

#endif

}
