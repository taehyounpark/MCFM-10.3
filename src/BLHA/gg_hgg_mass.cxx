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
  void gg_hgg_mass_nodecay_(double *p,int *ig1,int *ig2,double *msqv);
  void gg_hgg_mass_tb_nodecay_(double *p,int *ig1,int *ig2,double *msqv);
}

namespace MCFM {

  class gg_hgg_mass: public Process {
  private:

    int m_ig1, m_ig2, m_type;
    
  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (pi.m_oqcd!=4 || pi.m_oew!=1) return false;
      if (legs.size()!=5 || legs[0].m_fl!=25) return false;
      if (pi.m_decids.size()) return false;
      if (legs[1].m_fl==21 && legs[2].m_fl==21) {
	if (legs[3].m_fl==21 && legs[4].m_fl==21)// gggg
	  return interface->AddProcess(pi,new gg_hgg_mass(legs,1))>-1;
	if (legs[3].m_fl>0 && legs[3].m_fl<6 && legs[4].m_fl==-legs[3].m_fl)// ggqa
	  return interface->AddProcess(pi,new gg_hgg_mass(legs,2))>-1;
      }
      if (legs[1].m_fl>0 && legs[1].m_fl<6 && legs[3].m_fl==-legs[1].m_fl) {
       	if (legs[2].m_fl==legs[1].m_fl && legs[4].m_fl==-legs[2].m_fl)// aaaa
	  return interface->AddProcess(pi,new gg_hgg_mass(legs,8))>-1;
	if (legs[2].m_fl>0 && legs[2].m_fl<6 && legs[4].m_fl==-legs[2].m_fl)// qrqr
	  return interface->AddProcess(pi,new gg_hgg_mass(legs,10))>-1;
      }
      return false;
    }

    gg_hgg_mass(const std::vector<Leg> &legs,int type):
      Process(legs,1,2), m_ig1(5), m_ig2(6), m_type(type) {
      static int first(true);
      if (first) {
	first=false;
	nproc_.nproc=269;
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
      if (masses_.mb==0.) gg_hgg_mass_nodecay_(p_p,&m_ig1,&m_ig2,p_msqv);
      else gg_hgg_mass_tb_nodecay_(p_p,&m_ig1,&m_ig2,p_msqv);
      double res = p_msqv[MSQId(m_legs[1],m_legs[2])];
      m_res[3] = res*m_cfac;
    }

    int GetScheme() const { return 0; }

  };// end of class gg_hgg_mass

}// end of namespace MCFM
