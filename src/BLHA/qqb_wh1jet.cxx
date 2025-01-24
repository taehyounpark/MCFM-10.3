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
 

// MCFM routine for loop and tree MEs.
extern "C" {
  void qqb_wh1jet_v_(double *p,double *msqv);
}

namespace MCFM {

  class qqb_wh1jet: public Process {
  private:

    int m_id[6];

  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=2 || pi.m_oew!=3) return false;
      if (legs.size()!=6) return false;
      if (pi.m_decids.size()) return false;
      for (size_t i(1);i<legs.size();++i)
	if (legs[i].Mass()) return false;
      if (legs[0].m_fl==25 && legs[3].m_fl==21 &&
	  legs[1].m_fl>10 && legs[1].m_fl<17 &&
	  legs[2].m_fl==-legs[1].m_fl-1 &&
	  legs[4].m_fl>0 && legs[4].m_fl<6 &&
	  legs[5].m_fl==-legs[4].m_fl+1)
        return interface->AddProcess(pi,new qqb_wh1jet(legs))>-1;
      if (legs[0].m_fl==25 && legs[3].m_fl==21 &&
	  legs[1].m_fl>10 && legs[1].m_fl<17 &&
	  legs[2].m_fl==-legs[1].m_fl+1 &&
	  legs[4].m_fl>0 && legs[4].m_fl<6 &&
	  legs[5].m_fl==-legs[4].m_fl-1)
	return interface->AddProcess(pi,new qqb_wh1jet(legs))>-1;
      return false;
    }

    inline void SwapIDs(int i1,int i2) {
      std::swap<int>(m_id[i1],m_id[i2]);
      m_cfac=ISSymmetryFactor(m_legs,0);
    }

    qqb_wh1jet(const std::vector<Leg> &legs):
      Process(legs,5,4), m_id{5,4,1,2,0,3} {
      static int first(true);
      if (first) {
	first=false;
	nproc_.nproc=(legs[1].m_fl%2)?617:612;
	blha_.useblha=1;
	chooser_();
	std::string hdecaymode="none";
	hdecaymode.copy(hdecaymode_.hdecaymode,4);
      }
      m_res.resize(4);
      if (m_legs[4].m_is && m_legs[3].m_is) SwapIDs(0,5);
      if (m_legs[5].m_is && m_legs[3].m_is) SwapIDs(1,5);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[m_id[0]],0);// IS q
      SetMom(p_p,p,m_legs[m_id[1]],1);// IS qb'
      SetMom(p_p,p,m_legs[m_id[2]],2);// FS e- / v
      SetMom(p_p,p,m_legs[m_id[3]],3);// FS vb / e+
      SetMom(p_p,p,m_legs[m_id[4]],4);// FS h
      SetMom(p_p,p,m_legs[m_id[5]],6);// FS g
      GetMom(p_p,5,FourVec(0.,0.,0.,0.));// dummy
      // Calculate result.
      blha_.blhafl[0]=-MCFMId(m_legs[m_id[0]].m_fl);
      blha_.blhafl[1]=-MCFMId(m_legs[m_id[1]].m_fl);
      epinv_.epinv=epinv2_.epinv2=0.0;
      qqb_wh1jet_v_(p_p,p_msqv);
      double res = p_msqv[MSQId(m_legs[m_id[0]],m_legs[m_id[1]])];
      m_res[0] = res*m_cfac;
      if (!m_polecheck) return;
      epinv_.epinv=1.0;
      qqb_wh1jet_v_(p_p,p_msqv);
      double res1 = p_msqv[MSQId(m_legs[m_id[0]],m_legs[m_id[1]])];
      epinv2_.epinv2=1.0;
      qqb_wh1jet_v_(p_p,p_msqv);
      double res2 = p_msqv[MSQId(m_legs[m_id[0]],m_legs[m_id[1]])];
      m_res[1] = (res1-res)*m_cfac;
      m_res[2] = (res2-res1)*m_cfac;
      m_res[3] = m_res[2]/(-qcdcouple_.ason2pi*m_bfac);
    }

    int GetScheme() const { return 1; }

  };// end of class qqb_wh1jet

}// end of namespace MCFM
