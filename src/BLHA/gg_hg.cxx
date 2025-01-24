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
extern "C" { void gg_hg_v_nodecay_(double *p,int *iglue,double *msqv); }

namespace MCFM {

  class gg_hg: public Process {
  private:

    std::string m_hdecaymode;
    int m_id1, m_id2, m_id3, m_iglue, m_type;
    
  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=4 || pi.m_oew!=1) return false;
      if (legs.size()!=4 || legs[0].m_fl!=25) return false;
      if (pi.m_decids.size()) return false;
      for (size_t i(1);i<legs.size();++i)
	if (legs[i].Mass()) return false;
      if (legs[1].m_fl==21 && legs[2].m_fl==21 && legs[3].m_fl==21)
	return interface->AddProcess(pi,new gg_hg(legs,0))>-1;
      if (legs[2].m_fl>0 && legs[2].m_fl<6 &&
	  legs[3].m_fl==-legs[2].m_fl && legs[1].m_fl==21)
	return interface->AddProcess(pi,new gg_hg(legs,1))>-1;
      return false;
    }

    gg_hg(const std::vector<Leg> &legs,int type):
      Process(legs,3,2), m_iglue(4), m_type(type) {
      static int first(true);
      if (first) {
	first=false;
	removebr_.removebr=true;
	nproc_.nproc=203;
	blha_.useblha=1;
	chooser_();
        m_hdecaymode="none";
	m_hdecaymode.copy(hdecaymode_.hdecaymode,m_hdecaymode.size());
      }
      m_res.resize(4);
      //if (m_legs[2].m_is && m_legs[3].m_is) { m_id1=3; m_id2=2; m_id3=1; }
      //if (m_legs[2].m_is && m_legs[1].m_is) { m_id1=1; m_id2=2; m_id3=3; }
      //if (m_legs[3].m_is && m_legs[1].m_is) { m_id1=1; m_id2=3; m_id3=2; }
      m_id1=2; m_id2=3; m_id3=1;
//      m_cfac=ISSymmetryFactor(m_legs,0);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[m_id1],0);// IS 1
      SetMom(p_p,p,m_legs[m_id2],1);// IS 2
      SetMom(p_p,p,m_legs[m_id3],3);// FS
      SetMom(p_p,p,m_legs[0],2);// FS h
      // Calculate result.
      blha_.blhatype=m_type;
      epinv_.epinv=epinv2_.epinv2=0.0;
      gg_hg_v_nodecay_(p_p,&m_iglue,p_msqv);
      double res = p_msqv[MSQId(m_legs[m_id1],m_legs[m_id2])];
      m_res[0] = res*m_cfac;
      if (!m_polecheck) return;
      epinv_.epinv=1.0;
      gg_hg_v_nodecay_(p_p,&m_iglue,p_msqv);
      double res1 = p_msqv[MSQId(m_legs[m_id1],m_legs[m_id2])];
      epinv2_.epinv2=1.0;
      gg_hg_v_nodecay_(p_p,&m_iglue,p_msqv);
      double res2 = p_msqv[MSQId(m_legs[m_id1],m_legs[m_id2])];
      m_res[1] = (res1-res)*m_cfac;
      m_res[2] = (res2-res1)*m_cfac;
      m_res[3] = m_res[2]/(-qcdcouple_.ason2pi*m_bfac);
    }

    int GetScheme() const { return 0; }

  };// end of class gg_hg

}// end of namespace MCFM
