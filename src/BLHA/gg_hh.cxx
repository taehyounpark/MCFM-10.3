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
extern "C" { void gg_hh_(double *p,double *msqv); }

namespace MCFM {

  class gg_hh: public Process {
  private:

    std::string m_hdecaymode;
    
  public:

    static bool InitializeProcess(CXX_Interface *const interface,
				  const Process_Info &pi,
				  const std::vector<Leg> &legs) {
      if (legs.size()<4) return false;
      if (pi.m_decids.size()) return false;
      if (pi.m_oqcd!=2 || pi.m_oew!=2) return false;
      if (legs[0].m_fl==25 && legs[1].m_fl==25 && legs[2].m_fl==21 && legs[3].m_fl==21)
        return interface->AddProcess(pi,new gg_hh(legs))>-1;
      return false;
    }

    gg_hh(const std::vector<Leg> &legs):
      Process(legs,2,3) {
      static int first(true);
      if (first) {
        first=false;
        nproc_.nproc=601;
	blha_.useblha=1;
	chooser_();
        m_hdecaymode="none";
	m_hdecaymode.copy(hdecaymode_.hdecaymode,m_hdecaymode.size());
        masses_.hwidth   = 0.;
      }
      m_res.resize(4);
    }

    void Calc(const std::vector<FourVec> &p,int oqcd) {
      // Convert momenta.
      SetMom(p_p,p,m_legs[2],0);// IS g
      SetMom(p_p,p,m_legs[3],1);// IS g
      SetMom(p_p,p,m_legs[0],2);// FS h
      SetMom(p_p,p,m_legs[1],3);// FS h
      // Calculate result.
      gg_hh_(p_p,p_msqv);
      m_res[3] = p_msqv[mr(0,0)]*m_cfac;
    }

    int GetScheme() const { return 1; }

  };// end of class gg_hh

}// end of namespace MCFM
