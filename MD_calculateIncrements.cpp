#include "MD.h"
#include "phys_constants.h"

void MD::calculateIncrements()
{
    int_fast32_t i;
    Ek=0;
    for(i=0; i<N; ++i)
    {
        V[i] += F[i]*P_dtM;
        Vav += V[i];
        R[i] += V[i]*P_dt;
        Ek += boost::qvm::mag_sqr(V[i]);
    }
    Ek *= 0.5*P_M;
    T = Ek*_1d_N*_1d_k_PhysConst;
    t += P_dt;
}
