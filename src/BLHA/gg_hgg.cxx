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
extern "C" { void gg_hgg_v_(double *p,double *msqv); }

namespace MCFM {

  class gg_hgg: public Process {
  private:

    int m_type;
    
  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=5 || pi.m_oew!=1) return false;
      if (legs.size()!=5 || legs[0].m_fl!=25) return false;
      if (pi.m_decids.size()) return false;
      for (size_t i(1);i<legs.size();++i)
	if (legs[i].Mass()) return false;
      if (legs[1].m_fl==21 && legs[2].m_fl==21) {
	if (legs[3].m_fl==21 && legs[4].m_fl==21)// gggg
	  return interface->AddProcess(pi,new gg_hgg(legs,1))>-1;
	if (legs[3].m_fl>0 && legs[3].m_fl<6 && legs[4].m_fl==-legs[3].m_fl)// ggqa
	  return interface->AddProcess(pi,new gg_hgg(legs,2))>-1;
      }
      if (legs[1].m_fl>0 && legs[1].m_fl<6 && legs[3].m_fl==-legs[1].m_fl) {
       	if (legs[2].m_fl==legs[1].m_fl && legs[4].m_fl==-legs[2].m_fl)// aaaa
	  return interface->AddProcess(pi,new gg_hgg(legs,8))>-1;
	if (legs[2].m_fl>0 && legs[2].m_fl<6 && legs[4].m_fl==-legs[2].m_fl)// qrqr
	  return interface->AddProcess(pi,new gg_hgg(legs,10))>-1;
      }
      return false;
    }

    gg_hgg(const std::vector<Leg> &legs,int type):
      Process(legs,1,2), m_type(type) {
      static int first(true);
      if (first) {
	first=false;
	nproc_.nproc=271;
	blha_.useblha=1;
	chooser_();
	std::string dummy("none");
	dummy.copy(hdecaymode_.hdecaymode,4);
      }
      m_res.resize(4);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[1],0);// IS
      SetMom(p_p,p,m_legs[2],1);// IS
      SetMom(p_p,p,m_legs[3],4);// FS
      SetMom(p_p,p,m_legs[4],5);// FS
      GetMom(p_p,2,FourVec(0.,0.,0.,0.));// dummy
      SetMom(p_p,p,m_legs[0],3);// FS h
      // Calculate result.
      blha_.blhatype=m_type;
      blha_.blhafl[0]=MCFMId(m_legs[1].m_fl);
      blha_.blhafl[1]=MCFMId(m_legs[2].m_fl);
      blha_.blhafl[4]=MCFMId(m_legs[3].m_fl);
      blha_.blhafl[5]=MCFMId(m_legs[4].m_fl);
      epinv_.epinv=epinv2_.epinv2=0.0;
      gg_hgg_v_(p_p,p_msqv);
      double res = p_msqv[MSQId(m_legs[1],m_legs[2])];
      m_res[0] = res*m_cfac;
      if (!m_polecheck) return;
      epinv_.epinv=1.0;
      gg_hgg_v_(p_p,p_msqv);
      double res1 = p_msqv[MSQId(m_legs[1],m_legs[2])];
      epinv2_.epinv2=1.0;
      gg_hgg_v_(p_p,p_msqv);
      double res2 = p_msqv[MSQId(m_legs[1],m_legs[2])];
      m_res[1] = (res1-res)*m_cfac;
      m_res[2] = (res2-res1)*m_cfac;
      m_res[3] = m_res[2]/(-qcdcouple_.ason2pi*m_bfac);
    }

    int GetScheme() const { return 0; }

  };// end of class gg_hgg

}// end of namespace MCFM
