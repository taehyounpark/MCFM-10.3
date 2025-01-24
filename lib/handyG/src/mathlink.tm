:Evaluate: BeginPackage["handyG`"];
:Evaluate:
  Print["handyG-__VERSION__ by L. Naterop, Y. Ulrich, A. Signer"];
  G::usage = "G function";
  ClearCache::usage = "Clears handyG cache";
  SetHandyOptions::usage = "Sets evaluation options";
  MPLdel::usage = "difference between two successive terms at which the series expansion of MPLs is truncated.";
  LIdel::usage = "difference between two succesive terms at which the series expansion of PolyLogs is truncated.";
  hCircle::usage = "the size of the H\[ODoubleDot]lder circle \[Lambda]";
:Evaluate: Begin["`Private`"];
:Evaluate:
  args2r[a_]:=Re[N[a/.SubPlus|SubMinus->Identity]];
  args2i[a_]:=Im[N[a/.SubPlus|SubMinus->Identity]];
  args2e[a_]:=Switch[Head[#],
    SubPlus, 1,
    SubMinus, -1,
    Complex, Sign[Im[#]],
    _, 1]& /@ a;
:Begin:
:Function: gpl
:Pattern: G[a__]/;And @@ (NumberQ /@ ({a} /. SubPlus | SubMinus -> Identity))
:Arguments: {args2r[{a}],args2i[{a}],args2e[{a}] }
:ArgumentTypes: {RealList,RealList,IntegerList}
:ReturnType: Manual
:End:

:Begin:
:Function: wipecache
:Pattern: ClearCache[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Manual
:End:

:Begin:
:Function: setopts
:Pattern: SetHandyOptions[opts___Rule]
:Arguments: List@@( {MPLdel, LIdel, hCircle } /. {opts} /. s_Symbol->-1)
:ArgumentTypes: {Real, Real, Real}
:ReturnType: Manual
:End:

:Evaluate:
  End[];
  EndPackage[];

#include "mathlink.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#ifdef HAVE_QUAD
typedef __float128 real;
#define PUTREAL MLPutReal128
#else
typedef double real;
#define PUTREAL MLPutReal
#endif

typedef struct {real r,i;} complex;
typedef struct {complex c; signed char i0;} inum;

extern complex __gpl_module_MOD_g_superflatn(inum*,long*);

#ifdef HAVE_QUAD
void gpl(mlextended_double * re, int nr, mlextended_double * im, int ni, int*ieps, long ne)
#else
void gpl(double * re, long nr, double * im, long ni, int*ieps, long ne)
#endif
{
    assert(nr==ni);
    assert(nr==ne);
    inum input[nr];
    complex ans;
    long nnr = nr;
    for(long i=0;i<nr;i++)
    {
        input[i].c.r = *(re+i);
        input[i].c.i = *(im+i);
        input[i].i0  = *(ieps+i);
    }
    ans = __gpl_module_MOD_g_superflatn(&input[0],&nnr);
    
    
    if(ans.i == 0)
        PUTREAL(stdlink, ans.r);
    else
    {
        MLPutFunction(stdlink, "Complex", 2);
        PUTREAL(stdlink, ans.r);
        PUTREAL(stdlink, ans.i);
    }
}


void __maths_functions_MOD_clearcache();
void wipecache(void)
{
    __maths_functions_MOD_clearcache();
    MLPutSymbol(stdlink, "Null");
}

extern real __globals_MOD_hoeldercircle;
extern real __globals_MOD_lidelta;
extern real __globals_MOD_mpldelta;


#ifdef HAVE_QUAD
void setopts(mlextended_double MPLdel, mlextended_double LIdel, mlextended_double circ)
#else
void setopts(double MPLdel, double LIdel, double circ)
#endif
{
    if(MPLdel > 0)
        __globals_MOD_mpldelta = MPLdel;
    if (LIdel > 0)
        __globals_MOD_lidelta = LIdel;
    if (circ > 0)
        __globals_MOD_hoeldercircle = circ;

    MLPutFunction(stdlink, "List", 3);
      MLPutFunction(stdlink, "Rule", 2);
        MLPutSymbol(stdlink, "MPLdel");
        PUTREAL(stdlink, __globals_MOD_mpldelta);
      MLPutFunction(stdlink, "Rule", 2);
        MLPutSymbol(stdlink, "LIdel");
        PUTREAL(stdlink, __globals_MOD_lidelta);
      MLPutFunction(stdlink, "Rule", 2);
        MLPutSymbol(stdlink, "hCircle");
        PUTREAL(stdlink, __globals_MOD_hoeldercircle);
}

int main(int argc, char **argv)
{
    return MLMain(argc, argv);
}


