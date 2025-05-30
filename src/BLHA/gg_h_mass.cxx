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
 

// MCFM routine for loop ME.
extern "C" { void gg_h_mass_(double *p,double *msqv); }

namespace MCFM {

  class gg_h_mass: public Process {
  private:

    std::string m_hdecaymode;
    bool m_removebr;
    
  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (legs.size()<3 || legs.size()>4) return false;
      if (pi.m_decids.size()) return false;
      if (legs.size()==3) {
	if (pi.m_oqcd!=2 || pi.m_oew!=1) return false;
        if (legs[0].m_fl==25 && legs[1].m_fl==21 && legs[2].m_fl==21)
          return interface->AddProcess(pi,new gg_h_mass(legs))>-1;
      } else {
	// coupling order check needed
        if (legs[2].m_fl==21 && legs[3].m_fl==21) {
          if ( (legs[0].m_fl==22 && legs[1].m_fl==22)
            || (legs[0].m_fl==15 && legs[1].m_fl==-15)
            || (legs[0].m_fl==5 && legs[1].m_fl==-5) ) {
            return interface->AddProcess(pi,new gg_h_mass(legs))>-1;
          }
        }
      }
      return false;
    }

    gg_h_mass(const std::vector<Leg> &legs):
      Process(legs,legs.size()==3?1:2,legs.size()==3?2:3) {
      static int first(true);
      if (first) {
        first=false;
        if (legs.size()==3) {
          m_removebr=true;
          // g g -> H.
          if (legs[0].m_fl==25 && legs[1].m_fl==21 && legs[2].m_fl==21) {
            m_hdecaymode="none";
            nproc_.nproc=111;
          }
        } else {
          m_removebr=false;
          // g g -> H -> b bbar.
          if (legs[0].m_fl==5 && legs[1].m_fl==-5) {
            m_hdecaymode="bqba";
            nproc_.nproc=111;
          }
          // g g -> H -> tau- tau+.
          else if (legs[0].m_fl==15 && legs[1].m_fl==-15) {
            m_hdecaymode="tlta";
            nproc_.nproc=112;
          }
          // g g -> H -> gamma gamma.
          else if (legs[0].m_fl==22 && legs[1].m_fl==22) {
            m_hdecaymode="gaga";
            nproc_.nproc=119;
          }
          // others here...
        }
	removebr_.removebr=m_removebr;
	blha_.useblha=1;
	chooser_();
	m_hdecaymode.copy(hdecaymode_.hdecaymode,m_hdecaymode.size());
      }
      m_res.resize(4);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      if (m_removebr) {
        SetMom(p_p,p,m_legs[1],0);// IS g
        SetMom(p_p,p,m_legs[2],1);// IS g
        SetMom(p_p,p,m_legs[0],2);// FS h
      } else {
        SetMom(p_p,p,m_legs[2],0);// IS g
        SetMom(p_p,p,m_legs[3],1);// IS g
        SetMom(p_p,p,m_legs[0],2);// FS ptcl
        SetMom(p_p,p,m_legs[1],3);// FS antiptcl
      }
      // Calculate result.
      gg_h_mass_(p_p,p_msqv);
      m_res[3] = p_msqv[mr(0,0)]*m_cfac;
    }

    int GetScheme() const { return 1; }

  };// end of class gg_h_mass

}// end of namespace MCFM
