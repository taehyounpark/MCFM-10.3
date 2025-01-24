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
 

// MCFM routines for tree-level and loop MEs.
extern "C" {
  void gggaga_mt_(double *p,double *msqv);
  void gg_2gam_(double *p,double *msqv);
  void gg_2gam_v_(double *p,double *msqv);
}

namespace MCFM {

  class gg_2gam: public Process {
  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=2 || pi.m_oew!=2) return false;
      if (legs.size()!=4) return false;
      if (pi.m_decids.size()) return false;
      for (size_t i(0);i<legs.size();++i)
	if (legs[i].Mass()) return false;
      if (legs[0].m_fl!=22 || legs[1].m_fl!=22) return false;
      if (legs[2].m_fl==21 && legs[3].m_fl==21)
	return interface->AddProcess(pi,new gg_2gam(legs))>-1;
      return false;
    }

    gg_2gam(const std::vector<Leg> &legs):
      Process(legs,3,2) {
      static int first(true);
      if (first) {
	first=false;
	nproc_.nproc=285;
	blha_.useblha=1;
	chooser_();
      }
      m_res.resize(4);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[3],0);// IS g
      SetMom(p_p,p,m_legs[2],1);// IS g
      SetMom(p_p,p,m_legs[0],2);// FS y
      SetMom(p_p,p,m_legs[1],3);// FS y
      // Calculate result.
      gggaga_mt_(p_p,&p_msqv[mr(0,0)]);
      // gg_2gam_(p_p,p_msqv);
      m_res[3]=p_msqv[mr(0,0)]*m_cfac;
    }

    int GetScheme() const { return 0; }

  };// end of class gg_2gam

  class gg_2gam_v: public Process {
  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=2 || pi.m_oew!=2) return false;
      if (legs.size()!=4) return false;
      for (size_t i(0);i<legs.size();++i)
	if (legs[i].Mass()) return false;
      if (legs[0].m_fl!=22 || legs[1].m_fl!=22) return false;
      if (legs[2].m_fl==21 && legs[3].m_fl==21)
	return interface->AddProcess(pi,new gg_2gam_v(legs))>-1;
      return false;
    }

    gg_2gam_v(const std::vector<Leg> &legs):
      Process(legs,3,2) {
      static int first(true);
      if (first) {
	first=false;
	nproc_.nproc=285;
	blha_.useblha=1;
	chooser_();
      }
      m_res.resize(4);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[3],0);// IS g
      SetMom(p_p,p,m_legs[2],1);// IS g
      SetMom(p_p,p,m_legs[0],2);// FS y
      SetMom(p_p,p,m_legs[1],3);// FS y
      // Calculate result.
      epinv_.epinv=epinv2_.epinv2=0.0;
      gg_2gam_v_(p_p,p_msqv);
      double res = p_msqv[mr(0,0)];
      m_res[0] = res*m_cfac;
      if (!m_polecheck) return;
      epinv_.epinv=1.0;
      gg_2gam_v_(p_p,p_msqv);
      double res1 = p_msqv[mr(0,0)];
      epinv2_.epinv2=1.0;
      gg_2gam_v_(p_p,p_msqv);
      double res2 = p_msqv[mr(0,0)];
      m_res[1]=(res1-res)*m_cfac;
      m_res[2]=(res2-res1)*m_cfac;
      m_res[3]=m_res[2]/(-qcdcouple_.ason2pi*m_bfac);
    }

    int GetScheme() const { return 0; }

  };// end of class gg_2gam_v

}// end of namespace MCFM
