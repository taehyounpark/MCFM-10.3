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

int convert_pid(const char* in) {
  if (strcmp("d~",in)==0) return -1;
  if (strcmp("d",in)==0) return 1;
  if (strcmp("u~",in)==0) return -2;
  if (strcmp("u",in)==0) return 2;
  if (strcmp("s~",in)==0) return -3;
  if (strcmp("s",in)==0) return 3;
  if (strcmp("c~",in)==0) return -4;
  if (strcmp("c",in)==0) return 4;
  if (strcmp("b~",in)==0) return -5;
  if (strcmp("b",in)==0) return 5;
  if (strcmp("t~",in)==0) return -6;
  if (strcmp("t",in)==0) return 6;
  if (strcmp("e+",in)==0) return -11;
  if (strcmp("e-",in)==0) return 11;
  if (strcmp("mu+",in)==0) return -13;
  if (strcmp("mu-",in)==0) return 13;
  if (strcmp("tau+",in)==0) return -15;
  if (strcmp("tau-",in)==0) return 15;
  if (strcmp("ve~",in)==0) return -12;
  if (strcmp("ve",in)==0) return 12;
  if (strcmp("vmu~",in)==0) return -14;
  if (strcmp("vmu",in)==0) return 14;
  if (strcmp("vtau~",in)==0) return -16;
  if (strcmp("vtau",in)==0) return 16;
  if (strcmp("h",in)==0) return 25;
  if (strcmp("g",in)==0) return 21;
  if (strcmp("a",in)==0) return 22;
  if (strcmp("Z",in)==0) return 23;
  if (strcmp("W+",in)==0) return 24;
  if (strcmp("W-",in)==0) return -24;
  return 0;
}

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

// Main program.
int main(int argc, char** argv) {

  double Ecm = 350., mcfm_nptsfac=1.;
  int foEW = -1, foQCD = -1, diagonalCKM = 1;
  int nPts = 10, checkPoles = 0;
  int printMoms = 0, ampType = 11, fampType = -1;
  std::vector<int> decids, decfls;

  // Collect settings.
  std::map<std::string,std::string> params;
  read_parameters("params.lh",params);

  // Optional arguments.
  int opt, verbosity = 0;
  while ((opt = getopt(argc, argv, "pvmA:n:e:o:S:W:CD:P:h")) != -1) {
    switch (opt) {
    case 'p': checkPoles = 1; break;
    case 'v': verbosity += 1; break;
    case 'm': printMoms = 1; break;
    case 'n': nPts = atoi(optarg); break;
    case 'e': mcfm_nptsfac = atoi(optarg); break;
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
    case 'D': {
      size_t pos(0);
      std::string arg(optarg);
      for (;pos<arg.length();++pos)
	if (arg[pos]==',') break;
      decids.push_back(atoi(arg.substr(0,pos).c_str()));
      decfls.push_back(convert_pid(arg.substr(pos+1,std::string::npos).c_str()));
      break;
    }
    case 'h':
      std::cout<<argv[0]<<": Command line options\n"
	       <<"  -h        -- print this help message\n"
	       <<"  -v        -- increase verbosity\n"
	       <<"  -p        -- print ME values and enable pole check\n"
	       <<"  -m        -- print momentum configuration\n"
	       <<"  -C        -- enable non-diagonal CKM\n"
	       <<"  -n points -- set number of phase space points to 'points'\n"
	       <<"  -e factor -- set factor for MCFM extra points to 'points'\n"
	       <<"  -S oQCD   -- set QCD order to 'oQCD'\n"
	       <<"  -W oEW    -- set EW order to 'oEW'\n"
	       <<"  -A type   -- set OpenLoops amplitude type to 'type'\n"
	       <<"  -P p=val  -- set parameter 'p' to 'val'\n"
	       <<"  -D id,fl  -- add decay id and flavor 'id' and 'fl'\n";
      exit(0);							
    default:
      std::cerr<<"Usage: "<<argv[0]
	       <<" [-vpmC] [-n #pts] [-S oQCD] [-W oEW] [-A type] [-P p,val] [-D id,fl] [process]"
	       <<std::endl;
      exit(EXIT_FAILURE);
    }
  }
  if (optind >= argc) {
    std::cerr<<"You did not specify a process"<<std::endl;
    exit(EXIT_FAILURE);
  }

  //---------------------------------------------------------------------------
  // Process and initialization.

  // Process.
  std::vector<int> ids;
  int oEW = 0, oQCD = 1, nh = 0; // Needed for OpenLoops.
  for (int i(optind);i<argc;++i) {
    ids.push_back(convert_pid(argv[i]));
    if (ids.back()==21 || (abs(ids.back()))>0 && abs(ids.back())<7) ++oQCD;
    else {
      if (ids.back()==25) ++nh;
      ++oEW;
    }
  }
  // Note: in general, this does not capture all cases.
  oQCD-=2;
  oQCD = std::max(oQCD,0);

  int nIn = 2, nOut = ids.size()-nIn;

  // Fix coupling orders and amplitude types for special cases.
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

  //---------------------------------------------------------------------------
  // MCFM.

  // Create instance of MCFM interface.
  CXX_Interface mcfm;

  // Initialize MCFM.
  mcfm.Initialize(params);
  mcfm.SetVerbose(verbosity);

  // Print settings.
  if (verbosity) mcfm.PrintSettings();

  // Initialize process.
  Process_Info pi(ids,2,oQCD,oEW);
  pi.m_decids=decids;
  pi.m_decfls=decfls;
  pi.m_model=params["model"];
  int pid_mcfm(mcfm.InitializeProcess(pi));
  if (pid_mcfm < 0) {
    std::cout << " Process not available in MCFM." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Enable calculation of pole tems.
  mcfm.SetPoleCheck(checkPoles);

  //---------------------------------------------------------------------------
  // Recola.

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
  int pid_rcl = 1;
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
  
  //---------------------------------------------------------------------------
  // Test on phase space points generated with RAMBO.

  // Compare CPU time.
  clock_t t;
  double  cput_mcfm;
  double  cput_ol;
  double  cput_rcl;
  double  cput_ps;

  std::cout << "\n Generating " << nPts << " flat phase space points\n\n";

  if (checkPoles) {
    for (int iPt(0);iPt<nPts;++iPt) {
      // Generate phase space point with OpenLoops.
      std::vector<FourVec> pn;
      double pn_rcl[nIn+nOut][4];
      GenPoint(Ecm, pn, ids, nIn, params);
      for (size_t i(0); i<nIn+nOut; ++i) {
	pn_rcl[i][0] = pn[i][0];
	pn_rcl[i][1] = pn[i][1];
	pn_rcl[i][2] = pn[i][2];
	pn_rcl[i][3] = pn[i][3];
      }
      
      // Calculate result with Recola.
      double res_rcl[2] = {0., 0.};
      t = clock();
      Recola::compute_process_rcl(1, pn_rcl, "NLO", res_rcl);
      // For loop-induced Higgs process, put result in Born entry
      if (ids[2] == 25) res_rcl[0] = res_rcl[1];
      t = clock() - t;
      cput_rcl += double(t)/CLOCKS_PER_SEC;
    
      // Calculate result with MCFM.
      mcfm.SetAlphaS(std::stod(params["alpha_S"]));
      t = clock();
      mcfm.Calc(pid_mcfm, pn, 1);
      t = clock() - t;
      cput_mcfm += double(t)/CLOCKS_PER_SEC;
      std::vector<double> res_mcfm = mcfm.GetResult(pid_mcfm);
      // Translate MCFM result to same normalization and scheme conventions as Recola.
      for (int i(0);i<4;++i)
	res_mcfm[i] /= mcfm.GetProcess(pid_mcfm)->GetSymmetryFactor();
      res_mcfm[0] -= res_mcfm[3]*mcfm.GetProcess(pid_mcfm)->GetSchemeFactor(0);
    
      //  Print result and optionally phase space point.
      if (ampType == 11) {
        std::cout << " Finite: MCFM = " << std::setprecision(17) << res_mcfm[0] << " Recola = " << res_rcl[1]
                  << " Ratio = " << res_mcfm[0]/res_rcl[1] << std::endl;
                  // << " IR:     MCFM = " << res_mcfm[1] << " Recola = " << res_rcl[1]
                  // << " Ratio = " << res_mcfm[1]/res_rcl[1] << std::endl
                  // << " IR2:    MCFM = " << res_mcfm[2] << " Recola = " << res_rcl[2]
                  // << " Ratio = " << res_mcfm[2]/res_rcl[2] << std::endl;
      }
      std::cout << " Born:   MCFM = " << std::setprecision(17) << res_mcfm[3] << " Recola = " << res_rcl[0]
                << " Ratio = " << res_mcfm[3]/res_rcl[0]
                << std::endl << std::endl;
      if (printMoms) {
	for (size_t i(0);i<pn.size();++i)
	  std::cout << " p[" << i << "] = " << std::setprecision(17)
		    << pn[i] << " ( p^2 = " << m2(pn[i]) << " )\n";
	std::cout << std::endl;
      }
    }
  } else {
    // Phase space.
    t = clock();
    for (int iPt(0);iPt<nPts;++iPt) {
      // Generate phase space point with Rambo.
      std::vector<FourVec> pn;
      double pn_rcl[nIn+nOut][4];
      GenPoint(Ecm, pn, ids, nIn, params);
      for (size_t i(0); i<nIn+nOut; ++i) {
	pn_rcl[i][0] = pn[i][0];
	pn_rcl[i][1] = pn[i][1];
	pn_rcl[i][2] = pn[i][2];
	pn_rcl[i][3] = pn[i][3];
      }
    }
    t = clock() - t;
    cput_ps = double(t)/CLOCKS_PER_SEC;
    
    // Recola.
    t = clock();
    for (int iPt(0);iPt<nPts;++iPt) {
      // Generate phase space point with Rambo.
      std::vector<FourVec> pn;
      double pn_rcl[nIn+nOut][4];
      GenPoint(Ecm, pn, ids, nIn, params);
      for (size_t i(0); i<nIn+nOut; ++i) {
	pn_rcl[i][0] = pn[i][0];
	pn_rcl[i][1] = pn[i][1];
	pn_rcl[i][2] = pn[i][2];
	pn_rcl[i][3] = pn[i][3];
      }
      
      // Calculate result with Recola.
      double res_rcl[2] = {0., 0.};
      Recola::compute_process_rcl(1, pn_rcl, "NLO", res_rcl);
    }
    cput_rcl=(double(clock()-t)-cput_ps)/CLOCKS_PER_SEC;
    
    // MCFM.
    t = clock();
    for (int iPt(0);iPt<nPts*mcfm_nptsfac;++iPt) {
      // Generate phase space point with Rambo.
      std::vector<FourVec> pn;
      double pn_rcl[nIn+nOut][4];
      GenPoint(Ecm, pn, ids, nIn, params);
      for (size_t i(0); i<nIn+nOut; ++i) {
	pn_rcl[i][0] = pn[i][0];
	pn_rcl[i][1] = pn[i][1];
	pn_rcl[i][2] = pn[i][2];
	pn_rcl[i][3] = pn[i][3];
      }

      // Calculate result with MCFM.
//      for (int i(0);i<mcfm_nptsfac;++i) {
	mcfm.SetAlphaS(std::stod(params["alpha_S"]));
	mcfm.Calc(pid_mcfm, pn, 1);
//      }
      std::vector<double> res_mcfm = mcfm.GetResult(pid_mcfm);
      // Translate MCFM result to same normalization and scheme conventions as OpenLoops.
      for (int i(0);i<4;++i)
	res_mcfm[i] /= mcfm.GetProcess(pid_mcfm)->GetSymmetryFactor();
      res_mcfm[0] -= res_mcfm[3]*mcfm.GetProcess(pid_mcfm)->GetSchemeFactor(0);
    }
    cput_mcfm=(double(clock()-t)-cput_ps*mcfm_nptsfac)/CLOCKS_PER_SEC/mcfm_nptsfac;
  }

  // Print average CPU time.
  std::cout << " ----------------------------------------------------\n"
	    << "  Average CPU time per event:\n"
	    << "    MCFM:      " << std::setprecision(17) << std::scientific << cput_mcfm/double(nPts) << "s\n"
	    << "    Recola:    " << std::setprecision(17) << std::scientific << cput_rcl/double(nPts) << "s\n"
	    << "    Ratio:     " << cput_rcl/cput_mcfm << "\n"
	    << " ----------------------------------------------------\n";

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
