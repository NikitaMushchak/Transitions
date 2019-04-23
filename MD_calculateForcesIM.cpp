#include "MD.h"
#include <iomanip>
#include "phys_constants.h"
void MD::calculateForcesIM()
{
    boost::qvm::vec<double,3> dr, ps, f;
    Ep=0;
    Stress = MC_0T33;
    memset(F,0,sizeof(boost::qvm::vec<double,3>)*N);

    int_fast32_t i, j, iN=0;
    int_fast32_t m;
    boost::qvm::vec<int_fast32_t,3> Rj, vp;
    boost::qvm::mat<double,3,3> S;
    double dr_m, exp_a_r, _1d_r, dr_mm;
	double exp_a2_r;
    ///double C_coeff;
    //Vector2Di LayerNumber = {SampleCreation.LayerNumber.x,SampleCreation.LayerNumber.y}, Layer;
    //Vector2D PeriodicShift;
    //uint_fast32_t x, y, j=0;
    //InteractionsNumber=0;
    //memset(MacroParameters.C, 0, 5*sizeof(double));
    //std::cerr<<"W1\n";
    S = MC_0T33;
    for(i=0; i<N; ++i)
    {
        //if(VP[i])continue;
        for(m=0; m<NMd2; ++m)
        {
            Rj = Ri[i] + M[m];
            vp.a[0] = int_fast32_t(Rj.a[0]<0)-int_fast32_t(Rj.a[0]>=n.a[0]);
            vp.a[1] = int_fast32_t(Rj.a[1]<0)-int_fast32_t(Rj.a[1]>=n.a[1]);
            vp.a[2] = int_fast32_t(Rj.a[2]<0)-int_fast32_t(Rj.a[2]>=n.a[2]);
            j = I[Rj.a[0]+vp.a[0]*n.a[0]+(Rj.a[1]+vp.a[1]*n.a[1])*n.a[0]+(Rj.a[2]+vp.a[2]*n.a[2])*n.a[1]*n.a[0]];
            ps = double(vp.a[0])*PS[0]+double(vp.a[1])*PS[1]+double(vp.a[2])*PS[2];
            //std::cerr<<"F0 "<<i<<" "<<j<<" "<<vp.a[0]<<" "<<vp.a[1]<<" "<<Rj.a[0]<<" "<<Rj.a[1]<<" "<<n.a[0]<<" "<<n.a[1]
            //<<" "<<ps.a[0]<<" "<<ps.a[1]<<"\n";
            dr = R[j]-R[i]-ps;
            //std::cerr<<boost::qvm::mag_sqr(dr)<<" "<<Rj.a[0]<<" "<<Rj.a[1]<<" "<<Rj.a[2]<<" "<<i<<" "<<j<<"\n";
            //dr_mm = boost::qvm::mag_sqr(dr);
            //if(fabs(boost::qvm::mag_sqr(dr)-1.0)>1e-5)
            //std::cerr<<"ERROR! "<<boost::qvm::mag_sqr(dr)-1.0<<" "<<i<<" "<<j<<" "<<dr.a[0]<<" "<<dr.a[1]<<" "<<dr.a[2]<<" "<<ps.a[0]<<" "<<ps.a[1]<<" "<<ps.a[2]<<"\n";
            //if(dr_mm>P_aa_cut){std::cerr<<"ERROR! "<<i<<" "<<j<<" "<<dr.a[0]<<" "<<dr.a[1]<<" "<<dr_mm<<" "<<ps.a[0]<<" "<<ps.a[1]<<"\n";std::cin.get();}
            {
                dr_m = boost::qvm::mag(dr);
                _1d_r = 1.0/dr_m;
                exp_a_r = exp(P_alfa*(P_a-dr_m));
				exp_a2_r = exp(P_alfa2 * (P_a2 -dr_m) *(dr_m - P_a2));
                ///C_coeff = _1d_r*_1d_r*Potential_F*exp_a_r*(Potential_alfa*(2.0*exp_a_r-1.0) + _1d_r*(exp_a_r-1.0));
                f = ((P_F*exp_a_r*(exp_a_r-1.0) + P_F2 * exp_a2_r *( P_a2 - dr_m))*_1d_r)*dr;
                F[i]-=f;
                F[j]+=f;

                Ep += P_D*exp_a_r*(exp_a_r-2.0) + P_D2 * exp_a2_r ;
                Stress -= tens(dr,f);
                //S = tens(f,dr);
                S = tens(dr,dr);
                //std::cerr<<"S "<<S.a[0][0]<<" "<<S.a[1][1]<<" "<<f.a[0]*dr.a[0]<<" "<<f.a[1]*dr.a[1]<<" "<<dr_mm<<"\n";
                //std::cerr<<"S "<<S.a[0][0]<<" "<<S.a[1][1]<<" "<<f.a[0]<<" "<<dr.a[0]<<" "<<f.a[1]<<" "<<dr.a[1]<<" "<<dr_mm<<"\n";
                //std::cerr.precision(2);
                //std::cerr<<"F "<<std::setw(5)<<dr_mm-1<<" "<<(P_F*exp_a_r*(exp_a_r-1.0)*_1d_r)<<" "<<i<<" "<<j<<" "<<dr.a[0]<<" "<<dr.a[1]<<" "<<dr_m<<" "<<ps.a[0]<<" "<<ps.a[1]<<"\n";
                //std::cin.get();
            }
            ++iN;
        }
        //std::cin.get();
    }
    Stress *= _1d_Vol;

    //std::cerr<<"Stress "<<Stress.a[0][0]<<" "<<Stress.a[0][1]<<" "<<Stress.a[1][0]<<" "<<Stress.a[1][1]<<" "<<Vol-1<<" "<<iN<<" "<<N<<" "<<iN/N<<"\n";
    //std::cin.get();
}

void MD::calculateStateIM()
{
    boost::qvm::vec<double,3> dr, ps, f;
    Ep=0;
    Ek=0;
    Stress = MC_0T33;
    int_fast32_t i, j, iN=0;
    int_fast32_t m;
    boost::qvm::vec<int_fast32_t,3> Rj, vp;
    boost::qvm::mat<double,3,3> S;
    double dr_m, exp_a_r, _1d_r, dr_mm, c;
	double exp_a2_r;
    S = MC_0T33;
    P_Cmax = P_C;
    for(i=0; i<N; ++i)
    {
        //std::cerr<<i<<" ";
        for(m=0; m<NMd2; ++m)
        {
            Rj = Ri[i] + M[m];
            vp.a[0] = int_fast32_t(Rj.a[0]<0)-int_fast32_t(Rj.a[0]>=n.a[0]);
            vp.a[1] = int_fast32_t(Rj.a[1]<0)-int_fast32_t(Rj.a[1]>=n.a[1]);
            vp.a[2] = int_fast32_t(Rj.a[2]<0)-int_fast32_t(Rj.a[2]>=n.a[2]);
            //if(i==0)
            //if(abs(Ri[i].a[0])>N || abs(Ri[i].a[1])>N || abs(Ri[i].a[2])>N || abs(M[m].a[0])>N || abs(M[m].a[1])>N || abs(M[m].a[2])>N){
                //std::cerr<<i<<" "<<Rj.a[0]<<" "<<Rj.a[1]<<" "<<Rj.a[2]<<" "<<vp.a[0]<<" "<<vp.a[1]<<" "<<vp.a[2]<<" "<<Rj.a[0]+vp.a[0]*n.a[0]+(Rj.a[1]+vp.a[1]*n.a[1])*n.a[0]+(Rj.a[2]+vp.a[2]*n.a[2])*n.a[1]*n.a[0]<<" "<<N<<"\n";
            //    std::cerr<<i<<" "<<m<<" "<<Ri[i].a[0]<<" "<<Ri[i].a[1]<<" "<<Ri[i].a[2]<<" "<<M[m].a[0]<<" "<<M[m].a[1]<<" "<<M[m].a[2]<<"\n";
            //    std::cin.get();            }
            j = I[Rj.a[0]+vp.a[0]*n.a[0]+(Rj.a[1]+vp.a[1]*n.a[1])*n.a[0]+(Rj.a[2]+vp.a[2]*n.a[2])*n.a[1]*n.a[0]];
            ps = double(vp.a[0])*PS[0]+double(vp.a[1])*PS[1]+double(vp.a[2])*PS[2];
            dr = R[j]-R[i]-ps;


            {
                dr_m = boost::qvm::mag(dr);
                _1d_r = 1.0/dr_m;
                exp_a_r = exp(P_alfa*(P_a-dr_m));
				exp_a2_r = exp(P_alfa2 * (P_a2 -dr_m) *(dr_m - P_a2));
                ///C_coeff = _1d_r*_1d_r*Potential_F*exp_a_r*(Potential_alfa*(2.0*exp_a_r-1.0) + _1d_r*(exp_a_r-1.0));
                f = ((P_F*exp_a_r*(exp_a_r-1.0) + P_F2 * exp_a2_r *( P_a2 - dr_m))*_1d_r)*dr;
                c = P_C*exp_a_r*(2.0*exp_a_r-1.0) + P_C2 *exp_a2_r*(1.0-2.0 * P_alfa2*(dr_m - P_alfa2)*(dr_m - P_alfa2));
                Ep += P_D*exp_a_r*(exp_a_r-2.0) + P_D2 * exp_a2_r ;
                Stress -= tens(f,dr);
                //S = tens(f,dr);
                S = tens(dr,dr);
                //std::cerr<<"S "<<S.a[0][0]<<" "<<S.a[1][1]<<" "<<f.a[0]*dr.a[0]<<" "<<f.a[1]*dr.a[1]<<" "<<dr_mm<<"\n";
                //std::cerr<<"S "<<S.a[0][0]<<" "<<S.a[1][1]<<" "<<f.a[0]<<" "<<dr.a[0]<<" "<<f.a[1]<<" "<<dr.a[1]<<" "<<dr_mm<<"\n";
                //std::cerr.precision(2);
                //std::cerr<<"F "<<std::setw(5)<<dr_mm-1<<" "<<(P_F*exp_a_r*(exp_a_r-1.0)*_1d_r)<<" "<<i<<" "<<j<<" "<<dr.a[0]<<" "<<dr.a[1]<<" "<<dr_m<<" "<<ps.a[0]<<" "<<ps.a[1]<<"\n";
                //std::cin.get();
            }
            if(P_Cmax<c)P_Cmax=c;
            ++iN;
        }

        Ek += boost::qvm::mag_sqr(V[i]);

    }
    P_dt = 0.01*MC_s2d2*pi*sqrt(P_M/P_Cmax);
    P_dtM = P_dt*P_1d_M;
    //std::cerr<<"P_dt "<<P_dt<<" "<<P_Cmax<<"\n";
    Stress *= _1d_Vol;
    Ek *= 0.5*P_M;
    T = Ek*_1d_N*_1d_k_PhysConst;
    //std::cin.get();
}
