#include "MD.h"
#include "phys_constants.h"

void MD::calculateActualStiffness(const boost::qvm::mat<double,3,3> &dG)
{
    double EpsDeform_const = 1e-6;
    double _1d_Eps2_const = 1.0/(EpsDeform_const*EpsDeform_const);
    boost::qvm::mat<double,3,3> defGS[18], StressS[18], dStress, dGtemp = dG-MC_1T33;
    for(uint_fast32_t i=0; i<3; ++i)
    {
        for(uint_fast32_t j=0; j<3; ++j)
        {
            boost::qvm::set_identity(defGS[3*i+j]);
            defGS[3*i+j].a[i][j]   += EpsDeform_const;
            boost::qvm::set_identity(defGS[3*i+j+9]);
            defGS[3*i+j+9].a[i][j] -= EpsDeform_const;
        }
    }
    //defGS[0].a[0][0]   += 0.1;//EpsDeform_const;
    //defGS[0+4].a[0][0] -= 0.1;//EpsDeform_const;
    //defGS[2*1+1].a[1][1]   += 0.1;//EpsDeform_const;
    //defGS[2*1+1+4].a[1][1] -= 0.1;//EpsDeform_const;

    saveLattice();
    for(uint_fast32_t i=0; i<18; ++i)
    {
        defGS[i] += dGtemp;
        deformLattice(defGS[i]);
        calculateForcesIM();
        StressS[i] = Stress;
        restoreLattice();
        //std::cerr<<"StressS "<<" "<<Stress.a[0][0]<<" "<<Stress.a[0][1]<<" "<<Stress.a[1][0]<<" "<<Stress.a[1][1]<<"\n";

    }
    boost::qvm::set_zero(Stiffness);
    //std::cerr<<"dStress\n";
    for(uint_fast32_t i=0; i<9; ++i)
    {
        dStress = 0.5*(StressS[i]-StressS[i+9]);
        for(uint_fast8_t j=0;j<3;++j)
            for(uint_fast8_t k=0;k<3;++k)
                if(fabs(dStress.a[j][k])<1e-8)dStress.a[j][k]=0;
        defGS[i]-=MC_1T33+dGtemp;
        //std::cerr<<dStress.a[0][0]<<" "<<dStress.a[0][1]<<" "<<dStress.a[0][2]<<"\n"
        //<<dStress.a[1][0]<<" "<<dStress.a[1][1]<<" "<<dStress.a[1][2]<<"\n"
        //<<dStress.a[2][0]<<" "<<dStress.a[2][1]<<" "<<dStress.a[2][2]<<"\n"<<"\n";
        //std::cerr<<"Stress "<<" "<<defGS[i].a[0][0]<<" "<<defGS[i].a[0][1]<<" "<<defGS[i].a[1][0]<<" "<<defGS[i].a[1][1]<<"\n";
        Stiffness += _1d_Eps2_const*tens(dStress,defGS[i]);
    }
    //std::cerr<<"E1\n";
    toVoigtNotationC3D(Stiffness,StiffnessVoigt);
    //std::cerr<<Vol<<" "<<_1d_Vol<<" "<<N<<" "<<Vol*_1d_N<<"\n";
    /*for(uint_fast16_t i=0; i<9; ++i)
    {
        for(uint_fast16_t j=0; j<9; ++j)
        {
            if(fabs(Stiffness.a[i][j])<1e-9)
                Stiffness.a[i][j]=0;
            std::cerr.width(7);
            std::cerr.precision(3);
            std::cerr<<Stiffness.a[i][j]*stress_const*1e-9<<" ";
        }
        std::cerr<<"\n";
    }
    std::cerr<<"\n";/**/
    /*for(uint_fast16_t i=0; i<6; ++i)
    {
        for(uint_fast16_t j=0; j<6; ++j)
        {
            if(fabs(StiffnessVoigt.a[i][j])<1e-9)
                StiffnessVoigt.a[i][j]=0;
            std::cerr.width(7);
            std::cerr.precision(3);
            std::cerr<<StiffnessVoigt.a[i][j]*stress_const*1e-9<<" ";
        }
        std::cerr<<"\n";
    }
    std::cerr<<"\n";
    //std::cin.get();
    /**/
}
