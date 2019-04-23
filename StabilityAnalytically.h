#ifndef STABILITYANALYTICALLY_H
#define STABILITYANALYTICALLY_H
#include "math_constants.h"
#include "CalcBorder.h"
#include <boost/random.hpp>
class CalcBorder;
class StabilityAnalytically
{
public:
    StabilityAnalytically();
    virtual ~StabilityAnalytically();
    void checkStability2D(boost::qvm::vec<double,3> *Pvar, StabilityPointType *Pdvar, boost::qvm::mat<double,3,3> *Svar, uint_fast32_t &Nvar);
    void checkStability3D(boost::qvm::vec<double,3> *Pvar, StabilityPointType *Pdvar, boost::qvm::mat<double,3,3> *Svar, uint_fast32_t &Nvar);
    void startTask();



    void deformVectors(boost::qvm::mat<double,3,3> &dG);
    void checkStability(StabilityPointType &Pd);
    void calculateForces();
    void createTriangularLattice();
    void createFCCLattice();
    void findCoordinationalSpheres();
    void decomposeVectors();
    void calculateActualStiffness(const boost::qvm::mat<double,3,3> &dG=MC_1T33);
    void reachSpecifiedStress(const boost::qvm::mat<double,3,3> &Svar=MC_0T33);



    boost::qvm::vec<double,3> *e, *ein, *ea, *R, *RS;
    double *A, *Ain, *_1d_A, *P1, *P2, *CSD;
    uint_fast32_t *CSN;
    //CalcBorder *CB;
    boost::random::mt19937 generator;
    boost::random::uniform_on_sphere<double> distribution;
    boost::random::variate_generator<boost::random::mt19937&, boost::random::uniform_on_sphere<double> > distr;
    boost::qvm::vec<double,3> aStep;
    boost::qvm::vec<uint_fast32_t,3> n, nCenter;
    boost::qvm::mat<double,3,3> Stress;
    boost::qvm::mat<double,6,6> StiffnessVoigt, DuctilityVoigt;
    boost::qvm::mat<double,9,9> Stiffness, Ductility;

    double V, _1d_V, E, NearZero, NearSphere,
    P_aCut, P_P1, P_P2, P_D, P_alfa, P_a,
    EpsDeform_const, eXYMax;
    uint_fast32_t N, iCenter, CSDN, MaxCSDN, RSN;
protected:
private:
};

#endif // STABILITYANALYTICALLY_H
