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
 
// Simple test for the MCFM interface.

#include "MCFM/CXX_Interface.h"
#include "params.cxx"

using namespace MCFM;

// Main program.
int main() {

  // Collect settings.
  std::map<std::string,std::string> params;
  get_sm_parameters(params);

  // Create instance of MCFM interface and initialize.
  CXX_Interface mcfm;
  mcfm.Initialize(params);

  // Initialize process.
  std::vector<int> ids = {1, -1, 11, -11};
  Process_Info pi(ids,2,1,2);
  int pid(mcfm.InitializeProcess(pi));
  if (pid < 0) {
    std::cout << " Process not available." << std::endl;
    exit(EXIT_FAILURE);
  }

  // Test on a single phase space point
  std::vector<FourVec> pn(4);
  pn[0] = FourVec(1.750000e+02,0.000000e+00,0.000000e+00,1.750000e+02);
  pn[1] = FourVec(1.750000e+02,0.000000e+00,0.000000e+00,-1.750000e+02);
  pn[2] = FourVec(1.750000e+02,-9.824269e+01,5.893446e+01,1.322880e+02);
  pn[3] = FourVec(1.750000e+02,9.824269e+01,-5.893446e+01,-1.322880e+02);
  mcfm.Calc(pid, pn, 1);
    
  // Print result and phase space point.
  std::cout <<" Finite remainder = " <<mcfm.GetResult(pid)[0]<< std::endl;

  return 0;
}
