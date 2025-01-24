// Copyright (C) 2021 John Campbell, Stefan Hoeche, Christian T Preuss.

#ifndef MCFM__CXX_Wrapper_H
#define MCFM__CXX_Wrapper_H

#define MCFM_NF 5
#define MCFM_NMX 14

#include <complex>
#include <vector>
#include <iostream>
#include <iomanip>

extern "C" {

  // Common blocks
  extern struct {
    bool verbose;
  } verbose_;
  extern struct {
    double epinv;
  } epinv_;
  extern struct {
    double epinv2;
  } epinv2_;
  extern struct {
    char scheme[4];
  } scheme_;
  extern struct {
    double scale,musq;
  } mcfmscale_;
  extern struct {
    double gsq,as,ason2pi,ason4pi;
  } qcdcouple_;
  extern struct {
    double gf_inp,aemmz_inp,xw_inp,wmass_inp,zmass_inp;
  } ewinput_;
  extern struct {
    int ewscheme;
  } ewscheme_;
  extern struct {
    double Gf,gw,xw,gwsq,esq,vevsq;
  } ewcouple_;
  extern struct {
    int nflav;
  } nflav_;
  extern struct {
    double md,mu,ms,mc,mb,mt,
      mel,mmu,mtau,
      hmass,hwidth,
      wmass,wwidth,
      zmass,zwidth,
      twidth,
      tauwidth,
      mtausq,mcsq,mbsq;
  } masses_;
  extern struct {
    int n2, n3;
    double mass2, width2, mass3, width3;
  } breit_;
  extern struct {
    double Vud,Vus,Vub,Vcd,Vcs,Vcb;
  } cabib_;
  extern struct {
    int nlooprun;
  } nlooprun_;
  extern struct {
    char pdlabel[255];
  } pdlabel_;
  extern struct {
    double amz;
  } couple_;
  extern struct {
    int nproc,nprocbelow;
  } nproc_;
  extern struct {
    char hdecaymode[4], hdecaymode2[4];
  } hdecaymode_;
  extern struct {
    bool removebr;
  } removebr_;
  extern struct {
    bool zerowidth;
  } zerowidth_;
  extern struct {
    int qflag,gflag,qandgflag;
  } flags_;
  extern struct {
    double bbsqmin,bbsqmax,wsqmin,wsqmax,m3456min,m3456max;
  } limits_;
  extern struct {
    int nwz;
  } nwz_;
  extern struct {
    struct {double r,i;} zaemmz,zesq,rq1,rq2,
      zxw,zl,zr,zle,zln,zre,zrn,zsin2w,zl1,zr1,zl2,zr2;
  } zcouple_cms_;
  extern struct{
    double s[MCFM_NMX*MCFM_NMX];
  } sprods_;
  extern struct{
    int useblha, blhatype, blhafl[MCFM_NMX];
  } blha_;
  extern struct {
    char codeversion[6];
  } versionnumber_;
  extern struct {
    bool interference,bw34_56;
  } interference_;
  extern struct {
    double vsymfact;
  } vsymfact_;
  extern struct {
    double mt_yuk,mb_yuk,mc_yuk;
  } yukawas_;
  extern struct {
    int scalarselect;
  } scalarselect_;

#pragma omp threadprivate \
  (nflav_, masses_, breit_, epinv_, epinv2_,   \
   scheme_, mcfmscale_, qcdcouple_, ewcouple_, \
   zerowidth_, flags_, sprods_, interference_, vsymfact_)

  inline int mp(const int id, const int i) {return i*MCFM_NMX+id;}

  inline int mr(const int i, const int j) {
    return (j+MCFM_NF)*(2*MCFM_NF+1)+(i+MCFM_NF);}

}// end of extern "C"

namespace MCFM {

  // Simple four-vector struct.
  struct FourVec {
    double m_p[4];
    FourVec() {m_p[0]=0.;m_p[1]=0.;m_p[2]=0.;m_p[3]=0.;}
    FourVec(double e, double px, double py, double pz) {
      m_p[0]=e;m_p[1]=px;m_p[2]=py;m_p[3]=pz;}
    inline const double &operator [](int i) const {return m_p[i];}
    inline double& operator [](int i) {return m_p[i];}
    friend std::ostream& operator<<(std::ostream& os, const FourVec& p)
    { return os<<'('<<p[0]<<','<<p[1]<<','<<p[2]<<','<<p[3]<<')'; }
    inline bool operator==(const FourVec &p) const
    { return m_p[0]==p[0] && m_p[1]==p[1] && m_p[2]==p[2] && m_p[3]==p[3]; }
    inline FourVec operator-() const
    { return FourVec(-m_p[0],-m_p[1],-m_p[2],-m_p[3]); }
    inline void e(const double& in)  {m_p[0]=in;}
    inline void px(const double& in) {m_p[1]=in;}
    inline void py(const double& in) {m_p[2]=in;}
    inline void pz(const double& in) {m_p[3]=in;}
    inline double e()  const {return m_p[0];}
    inline double px() const {return m_p[1];}
    inline double py() const {return m_p[2];}
    inline double pz() const {return m_p[3];}
  };// end of struct FourVec

  // Calculate mass squared.
  inline double m2(const FourVec& p) {
    return p.e()*p.e()-p.px()*p.px()-p.py()*p.py()-p.pz()*p.pz();}

  // Order momenta as (px, py, pz, E).
  inline void GetMom(double* m, const int n, const FourVec& p) {
    m[mp(n,3)]=p[0]; for (int i(1);i<4;++i) m[mp(n,i-1)]=p[i];}

  // Translate PDG ID 21 to MCFM-used ID 0.
  inline int MCFMId(const int id) {return id==21?0:id;}

}// end of namespace MCFM

#endif
