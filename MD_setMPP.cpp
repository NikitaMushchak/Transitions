#include "MD.h"

void MD::setMPP()
{
    P_M = 1.0;
    P_1d_M = 1.0/P_M;
    P_a = 1.0;
    P_a2 = MC_2ds3;
    P_alfa = 13.0;
    P_alfa2 =350.0;
    P_D = 1.0;
    P_D2= 0.4;
    P_a_cut = 0.5*(CSDFCC[9]+CSDFCC[10]); // было между 3 и 4 координационной сферой
    P_aa_cut = P_a_cut*P_a_cut;
    P_F = 2.0*P_D*P_alfa/P_a;
    P_F2 = 2.0*P_D2*P_alfa2;
    P_C = P_F*P_alfa/P_a;
    P_C2= P_F2;

}
