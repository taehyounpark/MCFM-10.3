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
 
// Program to test MCFM interface against OpenLoops.

// Standard includes.
#include <ctime>
#include <fstream>
#include <algorithm>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

// MCFM includes
#include "MCFM/CXX_Interface.h"
#include "params.cxx"
#include "Rambo.cxx"

// Recola includes.
#include "recola.hpp"

void set_sm_parameters_rcl(std::map<std::string,std::string> &params,
			   const double lightFermionThresh=5.0);

using namespace MCFM;

std::string convert_pid_to_rcl(int pid) {
  switch (pid) {
  case -1:  return "d~";
  case 1:   return "d";
  case -2:  return "u~";
  case 2:   return "u";
  case -3:  return "s~";
  case 3:   return "s";
  case -4:  return "c~";
  case 4:   return "c";
  case -5:  return "b~";
  case 5:   return "b";
  case -6:  return "t~";
  case 6:   return "t";
  case -11: return "e+";
  case 11:  return "e-";
  case -13: return "mu+";
  case 13:  return "mu-";
  case -15: return "tau+";
  case 15:  return "tau-";
  case -12: return "nu_e~";
  case 12:  return "nu_e";
  case -14: return "nu_mu~";
  case 14:  return "nu_mu";
  case -16: return "nu_tau~";
  case 16:  return "nu_tau";
  case 25:  return "H";
  case 21:  return "g";
  case 22:  return "A";
  case 23:  return "Z";
  case 24:  return "W+";
  case -24: return "W-";
  default: break;
  }
  return "X";
}

using namespace MCFM;

// Main program.
int main(int argc, char** argv) {

  // Create instance of MCFM interface.
  CXX_Interface mcfm;

  // Read settings.
  std::map<std::string,std::string> params;
  read_parameters("params.lh",params);

  // Command-line options.
  int foEW = -1, foQCD = -1, diagonalCKM = 1, nIn = 2;
  int nPts = 1000000, checkPoles = 1;
  int printMom = 0, mode = 0, ampType = 11, fampType = -1;
  std::vector<int> decids, decfls;
  int opt, verbosity = 0;
  while ((opt = getopt(argc, argv, "vM:m:A:n:i:o:S:W:C:P:h")) != -1) {
    switch (opt) {
    case 'v': verbosity += 1; break;
    case 'M': printMom = atoi(optarg); break;
    case 'm': mode = atoi(optarg); break;
    case 'n': nPts = atoi(optarg); break;
    case 'i': nIn = atoi(optarg); break;
    case 'A': fampType = atoi(optarg); break;
    case 'W': foEW = atoi(optarg); break;
    case 'S': foQCD = atoi(optarg); break;
    case 'C': diagonalCKM = 0; break;
    case 'P': {
      size_t pos(0);
      std::string arg(optarg);
      for (;pos<arg.length();++pos) if (arg[pos]=='=') break;
      params[arg.substr(0,pos)]=arg.substr(pos+1,std::string::npos);
      break;
    }
    case 'h':
      std::cout<<argv[0]<<": Command line options\n"
	       <<"  -h        -- print this help message\n"
	       <<"  -v        -- increase verbosity\n"
	       <<"  -C        -- enable non-diagonal CKM\n"
	       <<"  -m mode   -- set stability check mode to 'mode'\n"
	       <<"  -n nin    -- set number of incoming particles to 'nin'\n"
	       <<"  -n points -- set number of phase space points to 'points'\n"
	       <<"  -S oQCD   -- set QCD order to 'oQCD'\n"
	       <<"  -W oEW    -- set EW order to 'oEW'\n"
	       <<"  -M print  -- set print mode to 'print'\n"
	       <<"  -P p=val  -- set parameter 'p' to 'val'\n"<<std::endl;
      exit(0);							
    default:
      std::cerr<<"Usage: "<<argv[0]<<" [-vmC] [-n #pts] [-S oQCD] [-W oEW] [-P p,val]"<<std::endl;
      exit(EXIT_FAILURE);
    }
  }
  if (optind >= argc) {
    std::cerr<<" Please provide an input file."<<std::endl;
    exit(EXIT_FAILURE);
  }

  //---------------------------------------------------------------------------
  // Read phase-space points.

  std::string infilename(argv[optind]);
  std::cout << "Input file: " << infilename << std::endl;
  std::ifstream infile(infilename);
  if (!infile.good())
    std::cout << "Input file " << infilename << " not found.\n";
  int nPt(0), nans(0), pid_mcfm(0), pid_rcl(1);
  std::vector<int> ids;
  std::vector<FourVec> pn;
  std::string line;
  std::vector<double> accs;
  while (std::getline(infile,line)) {
    line = reduce(line);
    line = chop(line,"!#");
    // Read configuration.
    if (line != "" && line != "eof") {
      int id(0);
      double e(0.0), px(0.0), py(0.0), pz(0.0);
      std::stringstream linestream(line);
      linestream.precision(16);
      linestream >> id >> e >> px >> py >> pz;
      ids.push_back(id);
      pn.push_back(FourVec(e,px,py,pz));
    }

    // At end of phase-space point do calculations.
    if (line == "" || line == "eof") {
      if (ids.size() > 0 && pn.size() > 0 && ids.size() == pn.size()) {
	// Count number of points and print this point.
	++nPt;
	if (printMom&2) {
	  std::cout << "Phase-space point " << nPt << std::endl;
	  FourVec psum(0.,0.,0.,0.);
	  for (size_t iPtcl(0); iPtcl<ids.size(); ++iPtcl) {
	    std::cout << ids[iPtcl] << " " << pn[iPtcl] << " "
		      << (pn[iPtcl][0]*pn[iPtcl][0]-pn[iPtcl][1]*pn[iPtcl][1]
			  -pn[iPtcl][2]*pn[iPtcl][2]-pn[iPtcl][3]*pn[iPtcl][3]) << std::endl;
	    psum[0]+=iPtcl<nIn?-pn[iPtcl][0]:pn[iPtcl][0];
	    psum[1]+=iPtcl<nIn?-pn[iPtcl][1]:pn[iPtcl][1];
	    psum[2]+=iPtcl<nIn?-pn[iPtcl][2]:pn[iPtcl][2];
	    psum[3]+=iPtcl<nIn?-pn[iPtcl][3]:pn[iPtcl][3];
	  }
	  std::cout<<"Sum "<<psum<<"\n";
	}
	int nOut = ids.size()-nIn;
	if (nPt==1) {
	  // Set up process parameters.
	  int oEW = 0, oQCD = 1, nh = 0;
	  for (size_t iPtcl(0); iPtcl<ids.size(); ++iPtcl) {
	    if (ids[iPtcl]==21 || (abs(ids[iPtcl])>0 && abs(ids[iPtcl])<7)) ++oQCD;
	    else {
	      if (ids[iPtcl]==25) ++nh;
	      ++oEW;
	    }
	  }
	  oQCD-=2;
	  if (nh && oEW == nh) {
	    params["H_width"] = "0";
	    if (params["model"] != "heft") {
	      ampType = 12;
	      oQCD -= 1;
	    }
	    oQCD += 2;
	  }
	  if (foEW>=0) oEW=foEW;
	  if (foQCD>=0) oQCD=foQCD;
	  if (fampType>0) ampType=fampType;
	  
	  // Initialize MCFM.
	  mcfm.Initialize(params);
	  Process_Info pi(ids,nIn,oQCD,oEW);
	  pi.m_decids=decids;
	  pi.m_decfls=decfls;
	  pi.m_model=params["model"];
	  pid_mcfm=mcfm.InitializeProcess(pi);
	  if (pid_mcfm < 0) {
	    std::cout << " Process not available in MCFM." << std::endl;
	    exit(EXIT_FAILURE);
	  }
	  mcfm.SetPoleCheck(checkPoles);
	
	  // Set Recola parameters.
	  set_sm_parameters_rcl(params);

	  // Verbosity.
	  if (verbosity) Recola::set_output_file_rcl("*");
	  Recola::set_print_level_parameters_rcl(verbosity);
	  Recola::set_print_level_amplitude_rcl(verbosity);
	  Recola::set_print_level_squared_amplitude_rcl(verbosity);
	
	  // Initialize process.
	  std::string processRcl;
	  for (int i(0);i<nIn;++i) processRcl += convert_pid_to_rcl(ids[i])+" ";
	  processRcl += "->";
	  for (int i(nIn);i<nIn+nOut;++i) processRcl += " "+convert_pid_to_rcl(ids[i]);
	  std::cout << " Recola process: " << processRcl << std::endl;
	  Recola::define_process_rcl(pid_rcl, processRcl, "NLO");

	  Recola::unselect_all_powers_BornAmpl_rcl(pid_rcl);
	  Recola::unselect_all_powers_LoopAmpl_rcl(pid_rcl);
	  // For loop-induced processes coupling specification is a little different 
	  if (ampType == 12) {
	    Recola::select_power_LoopAmpl_rcl(pid_rcl, "QED", oEW);
	    Recola::select_power_LoopAmpl_rcl(pid_rcl, "QCD", oQCD);
	  }
	  else {
	    Recola::select_power_BornAmpl_rcl(pid_rcl, "QED", oEW);
	    Recola::select_power_BornAmpl_rcl(pid_rcl, "QCD", oQCD-1);
	    Recola::select_power_LoopAmpl_rcl(pid_rcl, "QED", oEW);
	    Recola::select_power_LoopAmpl_rcl(pid_rcl, "QCD", oQCD+1);
	  }
  
	  // Generate process.
	  Recola::generate_processes_rcl();
	}
	std::vector<double> res_ol, res_mcfm;
	if (mode&1) {
	  double pn_rcl[nIn+nOut][4], res_rcl[2] = {0.,0.};
	  for (size_t i(0); i<nIn+nOut; ++i) {
	    pn_rcl[i][0] = pn[i][0];
	    pn_rcl[i][1] = pn[i][1];
	    pn_rcl[i][2] = pn[i][2];
	    pn_rcl[i][3] = pn[i][3];
	  }
	  Recola::compute_process_rcl(1, pn_rcl, "NLO", res_rcl);
	  res_ol={res_rcl[1],res_rcl[0]};
	}
	else {
	  mcfm.SetAlphaS(std::stod(params["alpha_S"]));
	  mcfm.Calc(pid_mcfm, pn, 1);
	  res_ol = mcfm.GetResult(pid_mcfm);
	  for (int i(0);i<4;++i)
	    res_ol[i] /= mcfm.GetProcess(pid_mcfm)->GetSymmetryFactor();
	  res_ol[0] -= res_ol[3]*mcfm.GetProcess(pid_mcfm)->GetSchemeFactor(0);
	}
	if (mode==0 || mode==3) {
	  for (size_t i(0);i<pn.size();++i)
	    pn[i] = FourVec(pn[i][0],pn[i][2],pn[i][3],pn[i][1]);
	}
	if (mode&2) {
	  double pn_rcl[nIn+nOut][4], res_rcl[2] = {0.,0.};
	  for (size_t i(0); i<nIn+nOut; ++i) {
	    pn_rcl[i][0] = pn[i][0];
	    pn_rcl[i][1] = pn[i][1];
	    pn_rcl[i][2] = pn[i][2];
	    pn_rcl[i][3] = pn[i][3];
	  }
	  Recola::compute_process_rcl(1, pn_rcl, "NLO", res_rcl);
	  res_mcfm={res_rcl[1],res_rcl[0]};
	}
	else {
	  mcfm.SetAlphaS(std::stod(params["alpha_S"]));
	  mcfm.Calc(pid_mcfm, pn, 1);
	  res_mcfm = mcfm.GetResult(pid_mcfm);
	  for (int i(0);i<4;++i)
	    res_mcfm[i] /= mcfm.GetProcess(pid_mcfm)->GetSymmetryFactor();
	  res_mcfm[0] -= res_mcfm[3]*mcfm.GetProcess(pid_mcfm)->GetSchemeFactor(0);
	}
	// Calculate accuracy.
	double acc = -log10(fabs((res_ol[0]-res_mcfm[0])/res_ol[0]));
	if (printMom&1) {
	  std::cout << std::setprecision(17) 
		    << ((mode&2)?"RCL:  ":"MCFM: ") << res_mcfm[0] << std::endl
		    << ((mode&1)?"RCL:  ":"MCFM: ") << res_ol[0] << std::endl
		    << "Acc:  " << acc << std::endl << std::endl;
	}
	if (std::isnan(acc) || std::isinf(acc)) ++nans;
	else accs.push_back(acc);
      }
      ids.clear();
      pn.clear();
      if (accs.size()>=nPts) break;
    }
  }
  double sum(0.), sum2(0.), n(accs.size());
  for (size_t i(0);i<n;++i) {
    sum+=accs[i];
    sum2+=accs[i]*accs[i];
  }
  std::cout << "Checked "
	    << ((mode&2)?"RCL":"MCFM") << " against "
	    << ((mode&1)?"RCL":"MCFM") << std::endl;
  std::cout<<"Average Acc: "<<sum/n<<" +- "
	   <<sqrt((sum2/n-pow(sum/n,2))/(n-1.0))<<"\n";
  std::sort(accs.begin(),accs.end());
  std::cout<<"Median  Acc: "<<accs[n/2]
	   <<" - "<<accs[n/2]-accs[n/4]
	   <<" + "<<accs[3*n/4]-accs[n/2]<<"\n";
  std::cout<<"Worst   Acc: "<<accs[0]<<" ("<<nans<<" nans)\n";
  std::cout<<"Best    Acc: "<<accs[n-1]<<"\n\n";
  // Done.
  Recola::reset_recola_rcl();
  return 0;
}

void set_sm_parameters_rcl(std::map<std::string,std::string> &params,
			   const double lightFermionThresh) {
  // Adapt the conventions from COLLIER to Catani-Seymour.
  Recola::set_delta_ir_rcl(0.0,M_PI*M_PI/6.0);
  // EW parameters.
  if (params["model"]!="heft") {
    Recola::use_alphaz_scheme_rcl(std::stod(params["alpha_EM"]));
    Recola::set_alphas_rcl(std::stod(params["alpha_S"]),std::stod(params["scale"]),-1);
    Recola::set_pole_mass_up_rcl(std::stod(params["up_mass"]));
    Recola::set_pole_mass_down_rcl(std::stod(params["down_mass"]));
    Recola::set_pole_mass_strange_rcl(std::stod(params["strange_mass"]));
    Recola::set_pole_mass_charm_rcl(std::stod(params["charm_mass"]),0.);
    Recola::set_pole_mass_bottom_rcl(std::stod(params["bottom_mass"]),0.);
    Recola::set_pole_mass_top_rcl(std::stod(params["top_mass"]),std::stod(params["top_width"]));
    Recola::set_pole_mass_electron_rcl(std::stod(params["electron_mass"]));
    Recola::set_pole_mass_muon_rcl(std::stod(params["muon_mass"]),0.);
    Recola::set_pole_mass_tau_rcl(std::stod(params["tau_mass"]),std::stod(params["tau_width"]));
    Recola::set_light_fermions_rcl(lightFermionThresh);
    Recola::set_pole_mass_z_rcl(std::stod(params["Z_mass"]),std::stod(params["Z_width"]));
    Recola::set_pole_mass_w_rcl(std::stod(params["W_mass"]),std::stod(params["W_width"]));
    Recola::set_pole_mass_h_rcl(std::stod(params["H_mass"]),std::stod(params["H_width"]));
  }
  else {
    Recola::set_parameter_rcl("aEW",std::stod(params["alpha_EM"]));
    Recola::set_parameter_rcl("aS",std::stod(params["alpha_S"]));
    Recola::set_parameter_rcl("MW",std::stod(params["W_mass"]));
    Recola::set_parameter_rcl("WW",std::stod(params["W_width"]));
    Recola::set_parameter_rcl("MZ",std::stod(params["Z_mass"]));
    Recola::set_parameter_rcl("WZ",std::stod(params["Z_width"]));
    Recola::set_parameter_rcl("MT",std::stod(params["top_mass"]));
    Recola::set_parameter_rcl("WT",std::stod(params["top_width"]));
  }
  Recola::set_mu_uv_rcl(std::stod(params["scale"]));
  Recola::set_mu_ir_rcl(std::stod(params["scale"]));
}
