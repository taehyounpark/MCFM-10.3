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

// OpenLoops includes.
#include "openloops.h"
void set_sm_parameters_ol(std::map<std::string,std::string> &params,
			  const int diagonalCKM=1);

using namespace MCFM;

// Convert FourVec struct to array used by OpenLoops.
void GetMomOL(const std::vector<FourVec>& pn, double* pn_ol);

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
      std::cerr<<"Usage: "<<argv[0]<<" [-vpmC] [-n #pts] [-S oQCD] [-W oEW] [-A type] [-P p,val] [-D id,fl] [process]"<<std::endl;
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

  int nIn = 2, nOut = ids.size()-nIn;

  //---------------------------------------------------------------------------
  // MCFM.

  // Create instance of MCFM interface.
  CXX_Interface mcfm;

  if (nh && oEW == nh) {
    params["H_width"] = "0";
    ol_setparameter_string("model", params["model"].c_str());
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

  // Set OpenLoops parameters.
  set_sm_parameters_ol(params,diagonalCKM);

  // Increase verbosity level to list loaded libraries.
  ol_setparameter_int("verbose", verbosity);

  // Initialize process.
  std::string process;
  for (int i(0);i<nIn;++i) process += std::to_string(ids[i])+" ";
  process += "->";
  for (int i(nIn);i<nIn+nOut;++i) process += " "+std::to_string(ids[i]);
  ol_setparameter_int("order_ew", oEW);
  int pid_ol = ol_register_process(process.c_str(), ampType);
  if (pid_ol < 0) {
    std::cout << " Process not available in OpenLoops." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Initialize OpenLoops.
  ol_start();
  
  //---------------------------------------------------------------------------
  // Test on phase space points generated with RAMBO.

  // Compare CPU time.
  clock_t t;
  double  cput_mcfm;
  double  cput_ol;
  double  cput_ps;

  std::cout << "\n Generating " << nPts << " flat phase space points\n\n";
  if (checkPoles) {
    for (int iPt(0);iPt<nPts;++iPt) {
      // Generate phase space point with OpenLoops.
      std::vector<FourVec> pn;
      double pn_ol[5*(nIn+nOut)];
      GenPoint(Ecm, pn, ids, nIn, params);
      GetMomOL(pn, pn_ol);

      // Calculate result with OpenLoops.
      double res_ol_tree, acc_ol;
      double res_ol_loop[3] = {0., 0., 0.};
      t = clock();
      if (ampType == 12) ol_evaluate_loop2(pid_ol, pn_ol, &res_ol_tree, &acc_ol);
      else ol_evaluate_loop(pid_ol, pn_ol, &res_ol_tree, res_ol_loop, &acc_ol);
      t = clock() - t;
      cput_ol += double(t)/CLOCKS_PER_SEC;
      std::vector<double> res_ol = {res_ol_loop[0], res_ol_loop[1],
                                    res_ol_loop[2], res_ol_tree};

      // Calculate result with MCFM.
      mcfm.SetAlphaS(std::stod(params["alpha_S"]));
      t = clock();
      mcfm.Calc(pid_mcfm, pn, 1);
      t = clock() - t;
      cput_mcfm += double(t)/CLOCKS_PER_SEC;
      std::vector<double> res_mcfm = mcfm.GetResult(pid_mcfm);
      // Translate MCFM result to same normalization and scheme conventions as OpenLoops.
      for (int i(0);i<4;++i)
        res_mcfm[i] /= mcfm.GetProcess(pid_mcfm)->GetSymmetryFactor();
      res_mcfm[0] -= res_mcfm[3]*mcfm.GetProcess(pid_mcfm)->GetSchemeFactor(0);
      
      // Print result and optionally phase space point.
      if (ampType == 11) {
        std::cout << " Finite: MCFM = " << std::setprecision(17) << res_mcfm[0] << " OpenLoops = " << res_ol[0]
                  << " Ratio = " << res_mcfm[0]/res_ol[0] << std::endl
                  << " IR:     MCFM = " << res_mcfm[1] << " OpenLoops = " << res_ol[1]
                  << " Ratio = " << res_mcfm[1]/res_ol[1] << std::endl
                  << " IR2:    MCFM = " << res_mcfm[2] << " OpenLoops = " << res_ol[2]
                  << " Ratio = " << res_mcfm[2]/res_ol[2] << std::endl;
      }
      std::cout << " Born:   MCFM = " << std::setprecision(17) << res_mcfm[3] << " OpenLoops = " << res_ol[3]
                << " Ratio = " << res_mcfm[3]/res_ol[3]
                << std::endl << std::endl;

      if (printMoms) {
        for (size_t i(0);i<pn.size();++i)
          std::cout << " p[" << i << "] = " << std::setprecision(17)
                    << pn[i] << " ( p^2 = " << m2(pn[i]) << " )\n";
        std::cout << std::endl;
      }
    }
  }
  else {
    t=clock();
    for (int iPt(0);iPt<nPts*100;++iPt) {
      // Generate phase space point with OpenLoops.
      std::vector<FourVec> pn;
      double pn_ol[5*(nIn+nOut)];
      GenPoint(Ecm, pn, ids, nIn, params);
      GetMomOL(pn, pn_ol);
    }
    cput_ps=double(clock()-t)/100.;

    t=clock();
    for (int iPt(0);iPt<nPts;++iPt) {
      // Generate phase space point with OpenLoops.
      std::vector<FourVec> pn;
      double pn_ol[5*(nIn+nOut)];
      GenPoint(Ecm, pn, ids, nIn, params);
      GetMomOL(pn, pn_ol);
      // Calculate result with OpenLoops.
      double res_ol_tree, res_ol_loop[3], acc_ol;
      ol_evaluate_loop(pid_ol, pn_ol, &res_ol_tree, res_ol_loop, &acc_ol);
    }
    cput_ol=(double(clock()-t)-cput_ps)/CLOCKS_PER_SEC;

    mcfm.SetPoleCheck(false);
    t=clock();
    for (int iPt(0);iPt<nPts*mcfm_nptsfac;++iPt) {
      // Generate phase space point with OpenLoops.
      std::vector<FourVec> pn;
      double pn_ol[5*(nIn+nOut)];
      GenPoint(Ecm, pn, ids, nIn, params);
      GetMomOL(pn, pn_ol);
      // Calculate result with MCFM.
//      for (int i(0);i<mcfm_nptsfac;++i) {
	mcfm.SetAlphaS(std::stod(params["alpha_S"]));
	mcfm.Calc(pid_mcfm, pn, 1);
//      }
    }
    cput_mcfm=(double(clock()-t)-cput_ps*mcfm_nptsfac)/CLOCKS_PER_SEC/mcfm_nptsfac;
  }

  // Print average CPU time.
  std::cout << " ----------------------------------------------------\n"
	    << "  Average CPU time per event:\n"
	    << "    MCFM:      " << std::setprecision(17) << std::scientific << cput_mcfm/double(nPts) << "s\n"
	    << "    OpenLoops: " << std::setprecision(17) << std::scientific << cput_ol/double(nPts) << "s\n"
	    << "    Ratio:     " << cput_ol/cput_mcfm << "\n"
	    << " ----------------------------------------------------\n";

  // Done.
  ol_finish();
  return 0;
}

// Get momenta as needed for OpenLoops.
void GetMomOL(const std::vector<FourVec>& pn, double* pn_ol) {
  for (size_t i(0); i<pn.size(); ++i) {
    pn_ol[5*i]   = pn[i].e();
    pn_ol[5*i+1] = pn[i].px();
    pn_ol[5*i+2] = pn[i].py();
    pn_ol[5*i+3] = pn[i].pz();
    pn_ol[5*i+4] = sqrt(m2(pn[i]));
  }
}

void set_sm_parameters_ol(std::map<std::string,std::string> &params,
			  const int diagonalCKM) {
  ol_setparameter_string("model", params["model"].c_str());
  ol_setparameter_double("mass(1)", std::stod(params["down_mass"]));
  ol_setparameter_double("mass(2)", std::stod(params["up_mass"]));
  ol_setparameter_double("mass(3)", std::stod(params["strange_mass"]));
  ol_setparameter_double("mass(4)", std::stod(params["charm_mass"]));
  ol_setparameter_double("yuk(4)", std::stod(params["charm_yukawa"]));
  ol_setparameter_double("mass(5)", std::stod(params["bottom_mass"]));
  ol_setparameter_double("yuk(5)", std::stod(params["bottom_yukawa"]));
  ol_setparameter_double("mass(6)", std::stod(params["top_mass"]));
  ol_setparameter_double("yuk(6)", std::stod(params["top_yukawa"]));
  ol_setparameter_double("width(6)", std::stod(params["top_width"]));
  ol_setparameter_double("mass(11)", std::stod(params["electron_mass"]));
  ol_setparameter_double("mass(13)", std::stod(params["muon_mass"]));
  ol_setparameter_double("mass(15)", std::stod(params["tau_mass"]));
  ol_setparameter_double("yuk(15)", std::stod(params["tau_yukawa"]));
  ol_setparameter_double("width(15)", std::stod(params["tau_width"]));
  ol_setparameter_double("mass(23)", std::stod(params["Z_mass"]));
  ol_setparameter_double("width(23)", std::stod(params["Z_width"]));
  ol_setparameter_double("mass(24)", std::stod(params["W_mass"]));
  ol_setparameter_double("width(24)", std::stod(params["W_width"]));
  ol_setparameter_double("mass(25)", std::stod(params["H_mass"]));
  ol_setparameter_double("width(25)", std::stod(params["H_width"]));
  // EW parameters.
  ol_setparameter_int("ew_scheme",2);
  ol_setparameter_double("alpha_qed_mz", std::stod(params["alpha_EM"]));
  // CKM elements.
  if (diagonalCKM) ol_setparameter_int("ckmorder", 0);
  else {
    ol_setparameter_double("VCKMdu", std::stod(params["CKM_u_d"]));
    ol_setparameter_double("VCKMsu", std::stod(params["CKM_u_s"]));
    ol_setparameter_double("VCKMbu", std::stod(params["CKM_u_b"]));
    ol_setparameter_double("VCKMdc", std::stod(params["CKM_c_d"]));
    ol_setparameter_double("VCKMsc", std::stod(params["CKM_c_s"]));
    ol_setparameter_double("VCKMbc", std::stod(params["CKM_c_b"]));
  }
  // Scales and couplings.
  ol_setparameter_double("alpha_s", std::stod(params["alpha_S"]));
  ol_setparameter_double("mu", std::stod(params["scale"]));
  ol_setparameter_double("muren", std::stod(params["scale"]));
}
