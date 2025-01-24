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
 

// MCFM routine for loop MEs.
extern "C" {
  void qqb_w_v_(double *p,double *msqv);
  void qqb_w_(double *p,double *msqv);}

namespace MCFM {

  class qqb_w: public Process {
  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=1 || pi.m_oew!=2) return false;
      if (legs.size()!=4) return false;
      if (pi.m_decids.size()) return false;
      for (size_t i(0);i<legs.size();++i)
	if (legs[i].Mass()) return false;
      if (legs[0].m_fl>10 && legs[0].m_fl<17 &&
	  legs[1].m_fl==-legs[0].m_fl-1 &&
	  legs[2].m_fl>0 && legs[2].m_fl<6 && !(legs[2].m_fl%2) &&
	  legs[3].m_fl<0 && legs[3].m_fl>-6 && (legs[3].m_fl%2))
	return interface->AddProcess(pi,new qqb_w(legs))>-1;
      if (legs[0].m_fl>10 && legs[0].m_fl<17 &&
	  legs[1].m_fl==-legs[0].m_fl+1 &&
	  legs[2].m_fl>0 && legs[2].m_fl<6 && (legs[2].m_fl%2) &&
	  legs[3].m_fl<0 && legs[3].m_fl>-6 && !(legs[3].m_fl%2))
	return interface->AddProcess(pi,new qqb_w(legs))>-1;
      return false;
    }

    qqb_w(const std::vector<Leg> &legs):
      Process(legs,3,2) {
      static int first(true);
      if (first) {
	first=false;
	nproc_.nproc=(legs[0].m_fl%2)?6:1;
	blha_.useblha=1;
	chooser_();
      }
      m_res.resize(4);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[3],0);// IS q
      SetMom(p_p,p,m_legs[2],1);// IS qb'
      SetMom(p_p,p,m_legs[0],2);// FS e- / v
      SetMom(p_p,p,m_legs[1],3);// FS vb / e+
      // Calculate result.
      epinv_.epinv=epinv2_.epinv2=0.0;
      qqb_w_v_(p_p,p_msqv);
      double res = p_msqv[MSQId(m_legs[3],m_legs[2])];
      m_res[0] = res*m_cfac;
      if (!m_polecheck) return;
      epinv_.epinv=1.0;
      qqb_w_v_(p_p,p_msqv);
      double res1 = p_msqv[MSQId(m_legs[3],m_legs[2])];
      epinv2_.epinv2=1.0;
      qqb_w_v_(p_p,p_msqv);
      double res2 = p_msqv[MSQId(m_legs[3],m_legs[2])];
      m_res[1] = (res1-res)*m_cfac;
      m_res[2] = (res2-res1)*m_cfac;
      m_res[3] = m_res[2]/(-qcdcouple_.ason2pi*m_bfac);
    }

    int GetScheme() const { return 1; }

  };// end of class qqb_w

}// end of namespace MCFM
