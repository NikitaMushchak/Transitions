#include "MD.h"
#include "phys_constants.h"

void MD::reachSpecifiedStress(const boost::qvm::mat<double,3,3> &Svar)
{
    boost::qvm::mat<double,3,3> dS=Stress-Svar, eps_tmp, eps_tmpsum;
    boost::qvm::set_zero(eps_tmpsum);
    boost::qvm::mat<double,6,6>  TEMP;
    uint_fast32_t steps=0;
    for(;;)
    {
        calculateForcesIM();
        dS=Stress-Svar;
        if(fabs(dS.a[0][0])+fabs(dS.a[1][1])+fabs(dS.a[2][2])<1e-12)
            break;

        calculateActualStiffness();
        DuctilityVoigt = boost::qvm::inverse(StiffnessVoigt);
        /*std::cerr<<"Ductility\n";
        for(uint_fast16_t i=0; i<6; ++i)
        {
            for(uint_fast16_t j=0; j<6; ++j)
            {
                if(fabs(DuctilityVoigt.a[i][j])<1e-9)
                    DuctilityVoigt.a[i][j]=0;
                std::cerr.width(7);
                std::cerr.precision(3);
                std::cerr<<DuctilityVoigt.a[i][j]*stress_const*1e-9<<" ";
            }
            std::cerr<<"\n";
        }
        std::cerr<<"\n";/**/


        fromVoigtNotationS3D(DuctilityVoigt,Ductility);

        /*for(uint_fast16_t i=0; i<9; ++i)
        {
            for(uint_fast16_t j=0; j<9; ++j)
            {
                if(fabs(Ductility.a[i][j])<1e-9)
                    Ductility.a[i][j]=0;
                std::cerr.width(7);
                std::cerr.precision(3);
                std::cerr<<Ductility.a[i][j]*stress_const*1e-9<<" ";
            }
            std::cerr<<"\n";
        }
        std::cerr<<"\n";/**/

        eps_tmp = -dotdot(Ductility,dS);

        eps_tmpsum += eps_tmp*(MC_1T33+eps_tmpsum);
        //std::cerr<<"Stress "<<"\n"<<dS.a[0][0]<<" "<<dS.a[0][1]<<" "<<dS.a[0][2]<<"\n"
        //         <<dS.a[1][0]<<" "<<dS.a[1][1]<<" "<<dS.a[1][2]<<"\n"
        //         <<dS.a[2][0]<<" "<<dS.a[2][1]<<" "<<dS.a[2][2]<<"\n";
//std::cerr<<"Stress "<<" "<<dS.a[0][0]<<" "<<dS.a[0][1]<<" "<<dS.a[1][0]<<" "<<dS.a[1][1]<<"\n";
        //std::cerr<<"E "<<"\n"<<eps_tmp.a[0][0]<<" "<<eps_tmp.a[0][1]<<" "<<eps_tmp.a[0][2]<<"\n"
        //         <<eps_tmp.a[1][0]<<" "<<eps_tmp.a[1][1]<<" "<<eps_tmp.a[1][2]<<"\n"
        //         <<eps_tmp.a[2][0]<<" "<<eps_tmp.a[2][1]<<" "<<eps_tmp.a[2][2]<<"\n";
//<<" "<<eps_tmpsum.a[0][0]<<" "<<eps_tmpsum.a[1][1]<<std::endl<<std::endl;

        eps_tmp += MC_1T33;
        deformLattice(eps_tmp);
        ++steps;


        //std::cin.get();
    }
    dG0 = eps_tmpsum;
    //std::cerr<<"E "<<eps_tmpsum.a[0][0]<<" "<<eps_tmpsum.a[1][1]<<" "<<eps_tmpsum.a[2][2]<<" "<<steps<<"\n";
    //std::cin.get();
}
