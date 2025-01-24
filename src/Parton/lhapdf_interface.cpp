/*
!  Copyright (C) 2019-2022, respective authors of MCFM.
!
!  This program is free software: you can redistribute it and/or modify it under
!  the terms of the GNU General Public License as published by the Free Software
!  Foundation, either version 3 of the License, or (at your option) any later
!  version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY
!  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
!  PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License along with
!  this program. If not, see <http://www.gnu.org/licenses/>
*/
 
#include "LHAPDF/PDF.h"
#include "LHAPDF/PDFSet.h"
#include "LHAPDF/PDFIndex.h"
#include "LHAPDF/Factories.h"
#include "LHAPDF/Utils.h"
#include "LHAPDF/Paths.h"
#include "LHAPDF/Version.h"
#include "LHAPDF/GridPDF.h"
#include <cstring>

using namespace std;

extern "C" {

void lhapdf_info()
{
    std::vector<std::string> sets = LHAPDF::availablePDFSets();
    std::cout << "Available PDF sets in LHAPDF:" << std::endl;
    for (int i=0; i<sets.size(); ++i) {
        std::cout << sets[i] << std::endl;
    }
}

void lhapdf_pathsPrepend(char path[])
{
    std::string str(path);
    LHAPDF::pathsPrepend(str);
}

LHAPDF::PDF* lhapdf_loadMember(char name[], int member)
{
    std::string str(name);
    return LHAPDF::mkPDF(str,member);
}

double lhapdf_evolve(LHAPDF::PDF* pdf, double x, double q2, int flav)
{
    return pdf->xfxQ2(flav, x, q2);
}

int lhapdf_number(char name[])
{
    std::string str(name);
    return LHAPDF::PDFSet(str).size();
}

double lhapdf_alphas(LHAPDF::PDF* pdf, double q2)
{
    return pdf->alphasQ2(q2);
}

int lhapdf_numFlavors(LHAPDF::PDF* pdf, double q)
{
    return pdf->alphaS().numFlavorsQ(q);
}

double lhapdf_quarkThreshold(LHAPDF::PDF* pdf, int i)
{
    return pdf->alphaS().quarkThreshold(i);
}

double lhapdf_quarkMass(LHAPDF::PDF* pdf, int i)
{
    return pdf->alphaS().quarkMass(i);
}

int lhapdf_orderqcd(LHAPDF::PDF* pdf)
{
    return pdf->orderQCD(); 
}

LHAPDF::PDFUncertainty lhapdf_computeUncertainty(LHAPDF::PDF *pdf, double *values, int nval)
{
    LHAPDF::PDFSet set = pdf->set();
    std::vector<double> vals(&values[0], &values[nval]);
    return set.uncertainty(vals);
}

void lhapdf_getconfig(LHAPDF::PDF *pdf, char key[], char* value, int vlen)
{
    LHAPDF::PDFInfo info = pdf->info();
    std::string skey(key);
    strncpy(value, info.get_entry(skey).c_str(), vlen);
}

int lhapdf_numxKnots(LHAPDF::PDF *pdf)
{
    const LHAPDF::GridPDF &gpdf = *dynamic_cast<const LHAPDF::GridPDF*>(pdf);
    return gpdf.xKnots().size();
}

int lhapdf_numq2Knots(LHAPDF::PDF *pdf)
{
    const LHAPDF::GridPDF &gpdf = *dynamic_cast<const LHAPDF::GridPDF*>(pdf);
    return gpdf.q2Knots().size();
}

void lhapdf_getxKnots(LHAPDF::PDF *pdf, double *xknots)
{
    const LHAPDF::GridPDF &gpdf = *dynamic_cast<const LHAPDF::GridPDF*>(pdf);
    for (int i=0; i<gpdf.xKnots().size(); ++i) {
        xknots[i] = gpdf.xKnots()[i];
    }
}

void lhapdf_getq2Knots(LHAPDF::PDF *pdf, double *q2knots)
{
    const LHAPDF::GridPDF &gpdf = *dynamic_cast<const LHAPDF::GridPDF*>(pdf);
    for (int i=0; i<gpdf.q2Knots().size(); ++i) {
        q2knots[i] = gpdf.q2Knots()[i];
    }
}

}

