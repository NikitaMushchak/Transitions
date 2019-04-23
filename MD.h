#ifndef MD_H
#define MD_H
#include <boost/qvm/vec.hpp>
#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/mat.hpp>
#include <boost/qvm/vec_mat_operations.hpp>
#include <boost/qvm/vec_operations.hpp>
#include <stdint.h>
#include "math_constants.h"
#include "DataTypes.h"
#include "lattice_constans.h"
class MD
{
    public:
        MD();
        virtual ~MD();

        void calculateForcesIM();
        void calculateStateIM();
        void calculateIncrements();

        void createLattice3DT(const boost::qvm::mat<double,3,3> &dG=MC_1T33);
        void createIMatrix(const bool &recreate=1);
        void setMPP();
        void calculateActualStiffness(const boost::qvm::mat<double,3,3> &dG=MC_1T33);
        void reachSpecifiedStress(const boost::qvm::mat<double,3,3> &Svar=MC_0T33);
        void deformLattice(const boost::qvm::mat<double,3,3> &dG);
        void saveLattice();
        void restoreLattice();
        void setTemperatureNormal();


        boost::qvm::vec<double,3> *R, *F, *V, *R0;
        boost::qvm::vec<int_fast32_t,3> *Ri, *M;
        uint_fast32_t *I;
        //uint_fast8_t *VP;
        boost::qvm::vec<double,3> LV[3], PS[3], _1d_PSd2[3], SC, PS0[3], Vav;
        boost::qvm::vec<int_fast32_t,3> n, n0, nCenter;
        int_fast32_t N, iCenter, NM, NMd2;
        double Vol, _1d_Vol, eXYMax, Vol0;
        boost::qvm::mat<double,6,6> StiffnessVoigt, DuctilityVoigt;
        boost::qvm::mat<double,9,9> Stiffness, Ductility;
        boost::qvm::mat<double,3,3> Stress, dG0;
        InitialConditionType IC;
        double P_a, P_alfa, P_a_cut, P_aa_cut, P_F, P_D, P_M, P_1d_M, P_dt, P_dtM, P_C, P_Cmax;
		double P_a2,P_alfa2,                   P_F2,P_D2,                          P_C2;
        double Ek, Ek0, Ep, Ep0, _1d_N, t, dt, T;
    protected:
    private:
};

#endif // MD_H
