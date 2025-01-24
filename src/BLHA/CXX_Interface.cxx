// Copyright (C) 2021 John Campbell, Stefan Hoeche, Christian T Preuss.
//
//  This program is free software: you can redistribute it and/or modify it under
//  the terms of the GNU General Public License as published by the Free Software
//  Foundation, either version 3 of the License, or (at your option) any later
//  version.
//
//  This program is distributed in the hope that it will be useful, but WITHOUT ANY
//  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
//  PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along with
//  this program. If not, see <http://www.gnu.org/licenses/>
 

#include "MCFM/CXX_Interface.h"
#include "MCFM/CXX_Wrapper.h"

// std headers
#include <algorithm>
#include <cstdlib>

#define PRINT_VAR(VAR) \
  std::cout<<__LINE__<<":"<<#VAR<<"="<<VAR<<std::endl

#define DEBUG_VAR(VAR) \
  if (CXX_Interface::Verbose()) std::cout<<__LINE__<<":"<<#VAR<<"="<<VAR<<std::endl

#define DEBUG_MSG() \
  if (CXX_Interface::Verbose()) std::cout

extern "C" { void chooser_(); }

using namespace MCFM;

std::vector<int> MCFM::ID(size_t id)
{
  std::vector<int> ids;
  for (size_t n(0);id>0;++n) {
    if (id&(1<<n)) {
      ids.push_back(n);
      id-=1<<n;
    }
  }
  return ids;
}

size_t MCFM::ID(const std::vector<int> &id)
{
  size_t ids(0);
  for (size_t n(0);n<id.size();++n) ids|=1<<n;
  return ids;
}

template <typename __Tp>
std::ostream &MCFM::operator<<
(std::ostream &str,const std::vector<__Tp> &v)
{
  str<<"(";
  if (v.size()>0) str<<v[0];
  else str<<"<no entry>";
  for (size_t i=1;i<v.size();++i) str<<","<<v[i];
    return str<<")";
}

template std::ostream &MCFM::operator<<
(std::ostream &str,const std::vector<int> &v);
template std::ostream &MCFM::operator<<
(std::ostream &str,const std::vector<double> &v);

int CXX_Interface::s_verbose(0);

Process::Process(const std::vector<Leg> &legs,int is1,int is2):
  m_legs(legs) {
  m_is[0]=is1;
  m_is[1]=is2;
  m_sf=ISSymmetryFactor(m_legs,0);
  m_sf*=FSSymmetryFactor(m_legs,0);
  m_cfac=ISSymmetryFactor(m_legs,1);
  m_cfac*=FSSymmetryFactor(m_legs,1);
  m_bfac=0.;
  for (size_t i(0);i<m_legs.size();++i)
    if (m_legs[i].Mass()==0.0) {
      if (abs(m_legs[i].Strong())==3) m_bfac+=4./3.;
      else if (m_legs[i].Strong()==8) m_bfac+=3.;
    }
  p_p = new double[4*MCFM_NMX];
  p_msqv = new double[(2*MCFM_NF+1)*(2*MCFM_NF+1)];
}

Process::~Process() {
  delete [] p_p;
  delete [] p_msqv;
}

double Process::GetSchemeFactor(int scheme,int mode) {
  int myscheme(GetScheme());
  if (scheme==myscheme) return 0.;
  double factor(0.0);
  for (size_t i(0);i<m_legs.size();++i)
    if (m_legs[i].Mass()==0.0) {
      if (abs(m_legs[i].Strong())==3) factor+=.5*4./3.;
      else if (m_legs[i].Strong()==8) factor+=3./6.;
    }
  if (mode==1) factor*=qcdcouple_.ason2pi;
  if (scheme==0 && myscheme==1) return factor;
  if (scheme==1 && myscheme==0) return -factor;
  return sqrt(-1.);
}

double Process::ISSymmetryFactor(const std::vector<Leg> &legs,
				 const int mode) const {
  double issf(1.);
  if (mode==0) {
    for (size_t i(0);i<legs.size();++i)
      if (legs[i].m_is) {
	if (legs[i].IsFermion()) issf*=2.;
	if (legs[i].IsVector()) issf*=legs[i].Mass()?3.:2.;
	if (legs[i].Strong()) issf*=abs(legs[i].Strong());
      }
  }
  else {
    for (size_t i(0);i<2;++i) {
      if (legs[m_is[i]].IsFermion()) issf*=2.;
      if (legs[m_is[i]].IsVector()) issf*=legs[m_is[i]].Mass()?3.:2.;
      if (legs[m_is[i]].Strong()) issf*=abs(legs[m_is[i]].Strong());
    }
  }
  return issf;
}

inline int factorial(int n) { return n<=1?1:factorial(n-1)*n; }

double Process::FSSymmetryFactor(const std::vector<Leg> &legs,
				 const int mode) const {
  std::vector<int> count(50,0);
  if (mode==0) {
    for (size_t i(0);i<m_legs.size();++i)
      if (!legs[i].m_is) ++count[legs[i].m_fl+24];
  }
  else {
    for (size_t i(0);i<m_legs.size();++i)
      if (i!=m_is[0] && i!=m_is[1]) ++count[legs[i].m_fl+24];
  }
  double fssf(1.);
  for (size_t i(0);i<count.size();++i) fssf*=factorial(count[i]);
  return fssf;
}

std::string Process::GetReferences() const
{
  return "";
}

CXX_Interface::CXX_Interface(const int banner): m_banner(banner) {
  // m_vno = std::string(versionnumber_.codeversion);
  m_vno = std::string("10.1  ");
  m_ref_list.push_back("J.M. Campbell, R.K. Ellis              hep-ph/9905386  ");
  m_ref_list.push_back("J.M. Campbell, R.K. Ellis, C. Williams arXiv:1105.0020 ");
  //m_ref_list.push_back("J.M. Campbell, R.K. Ellis, W. Giele    arXiv:1503.06182");
  m_ref_list.push_back("J.M. Campbell, T. Neumann              arXiv:1909.09117");
  m_ref_list.push_back("J.M. Campbell, S. Hoeche, C.T. Preuss  arXiv:2107.04472");
  m_cite_list.push_back("Virtual corrections computed by MCFM ");
  m_cite_list.push_back("\\cite{Campbell:1999ah,Campbell:2011bn,Campbell:2019dru,Campbell:2021vlt}.");
  if (m_banner) std::cout<<GetStartupMessage()<<std::endl;
}

CXX_Interface::~CXX_Interface() {
  if (m_banner) std::cout<<GetFinishMessage()<<std::endl;
}

std::string CXX_Interface::GetStartupMessage() const {
  std::string msg("*************************************************************\n");
  msg+="*                                                           *\n";
  msg+="*  MCFM Interface running MCFM v"+m_vno+"                      *\n";
  msg+="*                                                           *\n";
  msg+="*  Please cite:                                             *\n";
  for (size_t i(0);i<m_ref_list.size();++i) msg+="*  "+m_ref_list[i]+"  *\n";
  msg+="*                                                           *\n";
  msg+="*************************************************************\n";
  return msg;
}

std::string CXX_Interface::GetFinishMessage(const int aggregate) const {
  if (aggregate) return "";
  std::string msg("*************************************************************\n");
  msg+="*  Thanks for using MCFM.                                   *\n";
  msg+="*                                                           *\n";
  msg+="*  Please cite:                                             *\n";
  for (size_t i(0);i<m_ref_list.size();++i) msg+="*  "+m_ref_list[i]+"  *\n";
  msg+="*                                                           *\n";
  msg+="*************************************************************\n";
  return msg;
}

std::string CXX_Interface::GetReferences() const {
  std::string refs;
  for (size_t i(0);i<m_cite_list.size();++i) refs+=m_cite_list[i];
  return refs;
}

// Get parameter from list.
inline std::string CXX_Interface::ReadParam(std::string key,
  std::map<std::string,std::string>& params) {
  if (params.find(key) != params.end()) return params[key];
  std::cout << "Warning: parameter " << key << " not found!\n";
  return "0";
}

bool CXX_Interface::Initialize(std::map<std::string,std::string>& parameters) {
  // Check consistency of input parameters.
  if (!CheckInput(parameters)) {
    std::cout << " Error: inconsistent set of parameters." << std::endl;
    return false;
  }

  // General stuff.
  verbose_.verbose = s_verbose;
  nproc_.nproc     = -1;
  
  // Masses and widths.
  nflav_.nflav     = ToType<int>(ReadParam("n_flav", parameters));
  masses_.md       = ToType<double>(ReadParam("down_mass", parameters));
  masses_.mu       = ToType<double>(ReadParam("up_mass", parameters));
  masses_.ms       = ToType<double>(ReadParam("strange_mass", parameters));
  masses_.mc       = ToType<double>(ReadParam("charm_mass", parameters));
  masses_.mb       = ToType<double>(ReadParam("bottom_mass", parameters));
  masses_.mt       = ToType<double>(ReadParam("top_mass", parameters));
  masses_.twidth   = ToType<double>(ReadParam("top_width", parameters));
  masses_.mel      = ToType<double>(ReadParam("electron_mass", parameters));
  masses_.mmu      = ToType<double>(ReadParam("muon_mass", parameters));
  masses_.mtau     = ToType<double>(ReadParam("tau_mass", parameters));
  masses_.tauwidth = ToType<double>(ReadParam("tau_width", parameters));
  masses_.hmass    = ToType<double>(ReadParam("H_mass", parameters));
  masses_.hwidth   = ToType<double>(ReadParam("H_width", parameters));
  masses_.wmass    = ToType<double>(ReadParam("W_mass", parameters));
  masses_.wwidth   = ToType<double>(ReadParam("W_width", parameters));
  masses_.zmass    = ToType<double>(ReadParam("Z_mass", parameters));
  masses_.zwidth   = ToType<double>(ReadParam("Z_width", parameters));
  masses_.mtausq   = ToType<double>(ReadParam("tau_yukawa", parameters));
  yukawas_.mc_yuk  = ToType<double>(ReadParam("charm_yukawa", parameters));
  yukawas_.mb_yuk  = ToType<double>(ReadParam("bottom_yukawa", parameters));
  yukawas_.mt_yuk  = ToType<double>(ReadParam("top_yukawa", parameters));
  // masses_.mcsq     = ToType<double>(ReadParam("charm_yukawa", parameters));
  // masses_.mbsq     = ToType<double>(ReadParam("bottom_yukawa", parameters));
  //  masses_.mtsq     = ToType<double>(ReadParam("top_yukawa", parameters));
  
  s_flavors.m_mass[-1]  = s_flavors.m_mass[1]  = masses_.md;
  s_flavors.m_mass[-2]  = s_flavors.m_mass[2]  = masses_.mu;
  s_flavors.m_mass[-3]  = s_flavors.m_mass[3]  = masses_.ms;
  s_flavors.m_mass[-4]  = s_flavors.m_mass[4]  = masses_.mc;
  s_flavors.m_mass[-5]  = s_flavors.m_mass[5]  = masses_.mb;
  s_flavors.m_mass[-6]  = s_flavors.m_mass[6]  = masses_.mt;
  s_flavors.m_mass[-11] = s_flavors.m_mass[11] = masses_.mel;
  s_flavors.m_mass[-13] = s_flavors.m_mass[13] = masses_.mmu;
  s_flavors.m_mass[-15] = s_flavors.m_mass[15] = masses_.mtau;
  s_flavors.m_mass[-24] = s_flavors.m_mass[24] = masses_.wmass;
  s_flavors.m_mass[22]  = s_flavors.m_mass[21] = 0;
  s_flavors.m_mass[23]  = masses_.zmass;
  s_flavors.m_mass[25]  = masses_.hmass;

  s_flavors.m_width[-1]  = s_flavors.m_width[1]  = 0.;
  s_flavors.m_width[-2]  = s_flavors.m_width[2]  = 0.;
  s_flavors.m_width[-3]  = s_flavors.m_width[3]  = 0.;
  s_flavors.m_width[-4]  = s_flavors.m_width[4]  = 0.;
  s_flavors.m_width[-5]  = s_flavors.m_width[5]  = 0.;
  s_flavors.m_width[-6]  = s_flavors.m_width[6]  = masses_.twidth;
  s_flavors.m_width[-11] = s_flavors.m_width[11] = 0.;
  s_flavors.m_width[-13] = s_flavors.m_width[13] = 0.;
  s_flavors.m_width[-15] = s_flavors.m_width[15] = masses_.tauwidth;
  s_flavors.m_width[-24] = s_flavors.m_width[24] = masses_.wwidth;
  s_flavors.m_width[22]  = s_flavors.m_width[21] = 0;
  s_flavors.m_width[23]  = masses_.zwidth;
  s_flavors.m_width[25]  = masses_.hwidth;

  // EW parameters.
  ewscheme_.ewscheme = ToType<int>(ReadParam("ew_scheme", parameters));
  ewinput_.aemmz_inp = ToType<double>(ReadParam("alpha_EM", parameters));
  ewinput_.gf_inp    = ToType<double>(ReadParam("Gf", parameters));
  ewinput_.xw_inp    = ToType<double>(ReadParam("sin2_thetaW", parameters));
  ewinput_.wmass_inp = ToType<double>(ReadParam("W_mass", parameters));
  ewinput_.zmass_inp = ToType<double>(ReadParam("Z_mass", parameters));
  
  // CKM elements.
  cabib_.Vud = std::abs(ToType<std::complex<double> >(ReadParam("CKM_u_d", parameters)));
  cabib_.Vus = std::abs(ToType<std::complex<double> >(ReadParam("CKM_u_s", parameters)));
  cabib_.Vub = std::abs(ToType<std::complex<double> >(ReadParam("CKM_u_b", parameters)));
  cabib_.Vcd = std::abs(ToType<std::complex<double> >(ReadParam("CKM_c_d", parameters)));
  cabib_.Vcs = std::abs(ToType<std::complex<double> >(ReadParam("CKM_c_s", parameters)));
  cabib_.Vcb = std::abs(ToType<std::complex<double> >(ReadParam("CKM_c_b", parameters)));
  
  // Scales and strong coupling.
  mcfmscale_.scale   = ToType<double>(ReadParam("scale", parameters));
  mcfmscale_.musq    = (mcfmscale_.scale)*(mcfmscale_.scale);
  nlooprun_.nlooprun = ToType<int>(ReadParam("order_alpha_S", parameters));
  couple_.amz        = ToType<double>(ReadParam("alpha_S", parameters));
  if (!zerowidth_.zerowidth) limits_.bbsqmin = 1.;
  qcdcouple_.as      = ToType<double>(ReadParam("alpha_S", parameters));
  qcdcouple_.gsq     = 4.*M_PI*qcdcouple_.as;
  qcdcouple_.ason2pi = qcdcouple_.as/(2.*M_PI);
  qcdcouple_.ason4pi = qcdcouple_.as/(4.*M_PI);
  
  // Dummy PDF set.
  std::string dummy("mstw8lo");
  dummy.copy(pdlabel_.pdlabel,7);
  limits_.wsqmin = 1.e-6;
  limits_.wsqmax = 1.e99;
  
  // One-loop integral provider (QCDLoop = 1, OneLOop = 2)
  scalarselect_.scalarselect = 1;
  return true;
}

// Check consistency of parameters.

bool CXX_Interface::CheckInput(std::map<std::string,std::string>& parameters) {

  // All lepton masses and widths should be zero. (Except tau?)
  if (ToType<double>(ReadParam("electron_mass", parameters)) > 0.) return false;
  if (ToType<double>(ReadParam("muon_mass", parameters)) > 0.) return false;
  //if (ToType<double>(ReadParam("tau_mass", parameters)) > 0.) return false;
  //if (ToType<double>(ReadParam("tau_width", parameters)) > 0.) return false;

  // All quarks (except top and bottom) massless.
  if (ToType<double>(ReadParam("down_mass", parameters)) > 0.) return false;
  if (ToType<double>(ReadParam("up_mass", parameters)) > 0.) return false;
  if (ToType<double>(ReadParam("strange_mass", parameters)) > 0.) return false;
  if (ToType<double>(ReadParam("charm_mass", parameters)) > 0.) return false;
  // if (ToType<double>(ReadParam("bottom_mass", parameters)) > 0.) return false;

  // TODO: could check consistency of EW parameters here.

  return true;
}

// Print all settings.
void CXX_Interface::PrintSettings() {
  int nwdth = 29;//48 ws
  std::stringstream ss;
  ss << " **************************************************\n"
     << " *                                                *\n"
     << " *  Settings:                                     *\n"
     << " *                                                *\n";
  // General stuff.
  ss << " *   nproc       = " << std::setw(nwdth) << std::left << nproc_.nproc << "  *\n"
     << " *   verbose     = " << std::setw(nwdth) << std::left << verbose_.verbose << "  *\n"
     << " *                                                *\n";
  
  // Masses and widths.
  ss << " *  Masses and widths:                            *\n"
     << " *   nflav       = " << std::setw(nwdth) << std::left << nflav_.nflav << "  *\n"
     << " *   md          = " << std::setw(nwdth) << std::left << masses_.md << "  *\n"
     << " *   mu          = " << std::setw(nwdth) << std::left << masses_.mu << "  *\n"
     << " *   ms          = " << std::setw(nwdth) << std::left << masses_.ms << "  *\n"
     << " *   mc          = " << std::setw(nwdth) << std::left << masses_.mc << "  *\n"
     << " *   mb          = " << std::setw(nwdth) << std::left << masses_.mb << "  *\n"
     << " *   mt          = " << std::setw(nwdth) << std::left << masses_.mt << "  *\n"
     << " *   twidth      = " << std::setw(nwdth) << std::left << masses_.twidth << "  *\n"
     << " *   mel         = " << std::setw(nwdth) << std::left << masses_.mel << "  *\n"
     << " *   mmu         = " << std::setw(nwdth) << std::left << masses_.mmu << "  *\n"
     << " *   mtau        = " << std::setw(nwdth) << std::left << masses_.mtau << "  *\n"
     << " *   tauwidth    = " << std::setw(nwdth) << std::left << masses_.tauwidth << "  *\n"
     << " *   hmass       = " << std::setw(nwdth) << std::left << masses_.hmass << "  *\n"
     << " *   hwidth      = " << std::setw(nwdth) << std::left << masses_.hwidth << "  *\n"
     << " *   wmass       = " << std::setw(nwdth) << std::left << masses_.wmass << "  *\n"
     << " *   wwidth      = " << std::setw(nwdth) << std::left << masses_.wwidth << "  *\n"
     << " *   zmass       = " << std::setw(nwdth) << std::left << masses_.zmass << "  *\n"
     << " *   zwidth      = " << std::setw(nwdth) << std::left << masses_.zwidth << "  *\n"
     << " *                                                *\n";

  // EW parameters.
  ss << " *  EW parameters:                                *\n"
     << " *   ewscheme    = " << std::setw(nwdth) << std::left << ewscheme_.ewscheme << "  *\n"
     << " *   gf          = " << std::setw(nwdth) << std::left << ewinput_.gf_inp << "  *\n"
     << " *   xw_inp      = " << std::setw(nwdth) << std::left << ewinput_.xw_inp << "  *\n"
     << " *   zmass       = " << std::setw(nwdth) << std::left << ewinput_.zmass_inp << "  *\n"
     << " *   wmass       = " << std::setw(nwdth) << std::left << ewinput_.wmass_inp << "  *\n"
     << " *                                                *\n";

  // CKM elements.
  ss << " *  CKM elements:                                 *\n"
     << " *   Vud         = " << std::setw(nwdth) << std::left << cabib_.Vud << "  *\n"
     << " *   Vus         = " << std::setw(nwdth) << std::left << cabib_.Vus << "  *\n"
     << " *   Vub         = " << std::setw(nwdth) << std::left << cabib_.Vub << "  *\n"
     << " *   Vcs         = " << std::setw(nwdth) << std::left << cabib_.Vcs << "  *\n"
     << " *   Vcb         = " << std::setw(nwdth) << std::left << cabib_.Vcb << "  *\n"
     << " *                                                *\n";
  
  // Scales and strong coupling.
  ss << " *  Scales and couplings:                         *\n"
     << " *   alpha_S     = " << std::setw(nwdth) << std::left << qcdcouple_.as << "  *\n"
     << " *   scale       = " << std::setw(nwdth) << std::left << mcfmscale_.scale << "  *\n"
     << " *   musq        = " << std::setw(nwdth) << std::left << mcfmscale_.musq << "  *\n"
     << " *   nlooprun    = " << std::setw(nwdth) << std::left << nlooprun_.nlooprun << "  *\n"
     << " *                                                *\n";
  ss << " **************************************************\n";

  // Print.
  std::cout << ss.str();
}

// process specific code
#include "qqb_w.cxx"
#include "qqb_w1jet.cxx"
#include "qqb_w2jet.cxx"
#include "qqb_wh.cxx"
#include "qqb_wh1jet.cxx"
#include "qqb_z.cxx"
#include "qqb_z1jet.cxx"
#include "qqb_z2jet.cxx"
#include "qqb_zh.cxx"
#include "qqb_zh1jet.cxx"
#include "qqb_dirgam.cxx"
#include "qqb_gam2j.cxx"
#include "gg_2gam.cxx"
#include "qqb_gamgam.cxx"
#include "qqb_gmgmjt.cxx"
#include "qqb_trigam.cxx"
#include "qqb_fourgam.cxx"
#include "qqb_QQb.cxx"
#include "gg_h.cxx"
#include "gg_hg.cxx"
#include "gg_hgg.cxx"
#include "gg_h_mass.cxx"
#include "gg_hg_mass.cxx"
#include "gg_hgg_mass.cxx"
#include "qqb_ww.cxx"
#include "qqb_wz.cxx"
#include "qqb_zz.cxx"
#include "qqb_wgam.cxx"
#include "qqb_zgam.cxx"
#include "qqb_twojet.cxx"
#include "gg_hh.cxx"
#include "qqb_zzg.cxx"
#include "qqb_wwg.cxx"
#include "qqb_wzg.cxx"
#include "qqb_wgamg.cxx"
#include "qqb_zgamg.cxx"
// end process specific code

int CXX_Interface::InitializeProcess(const Process_Info &pi) {
  FlavorMulti_Map fmm;
  std::vector<Leg> legs(pi.m_ids.size());
  for (size_t i(0);i<pi.m_ids.size();++i) {
    legs[i]=i<pi.m_nin?Leg(pi.m_ids[i],i,0).Cross():Leg(pi.m_ids[i],i,0);
    if (fmm.find(abs(pi.m_ids[i]))==fmm.end()) fmm[abs(pi.m_ids[i])]=0;
    fmm[abs(pi.m_ids[i])]+=1;
  }
  std::sort(legs.begin(),legs.end(),Order_Flavor(&fmm));
  if (CXX_Interface::Verbose()) {
    std::cout<<__PRETTY_FUNCTION__<<": Requested process {\n";
    std::cout<<"  O(\\alpha_s)="<<pi.m_oqcd<<", O(\\alpha)="<<pi.m_oew<<"\n";
    if (pi.m_decids.size()) {
      std::cout<<"  Decays";
      for (size_t i(0);i<pi.m_decids.size();++i)
	std::cout<<" "<<ID(pi.m_decids[i])<<","<<pi.m_decfls[i];
      std::cout<<"\n";
    }
    for (size_t i(0);i<legs.size();++i)
      std::cout<<"  "<<i<<": "<<legs[i]<<"\n";
    std::cout<<"}\n";
  }
  for (size_t i(0);i<legs.size();++i)
    if (legs[i].Width()) {
      std::cout<<"External particle "<<legs[i].m_fl
	       <<" has width "<<legs[i].Width()<<"!"<<std::endl;
      return m_procmap[pi.UID()][pi.m_ids]=-1;
    }
  int procID(GetProcessID(pi));
  if (procID>=0) return procID;
  // process specific code
  if (pi.m_model=="sm") {
  if (qqb_w::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_w1jet::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_w2jet::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_wh::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_wh1jet::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_z::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_z1jet::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_z2jet::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_zh::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_zh1jet::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_dirgam::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_gam2j::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (gg_2gam::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_gamgam::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_gmgmjt::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_trigam::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_fourgam::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_QQb::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (gg_h_mass::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (gg_hg_mass::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (gg_hgg_mass::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_ww::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_wz::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_zz::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_wgam::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_zgam::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_twojet::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (gg_hh::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_zzg::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_wwg::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_wzg::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_wgamg::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (qqb_zgamg::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  }
  if (pi.m_model=="heft") {
  if (gg_h::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (gg_hg::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  if (gg_hgg::InitializeProcess(this,pi,legs)) return GetProcessID(pi);
  }
  // end process specific code
  return m_procmap[pi.UID()][pi.m_ids]=-1;
}

extern "C" {

  // This is needed on Linux and for the test programs only.
  // Make up for definition of qli2 in C++ version of QCDLoop.
  std::complex<double> qli2
  (double const& p1, double const& m1, double const& m2,
  double const& mu2, int const& ep);
  std::complex<double>
  qli2_(double const& p1, double const& m1, double const& m2,
  	double const& mu2, int const& ep)
  { return qli2(p1,m1,m2,mu2,ep); }

}

