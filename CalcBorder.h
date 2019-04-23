#ifndef CALCBORDER_H
#define CALCBORDER_H
#include <boost/qvm/vec.hpp>
#include <boost/qvm/vec_access.hpp>
#include <boost/qvm/mat.hpp>
#include <boost/qvm/vec_mat_operations.hpp>
#include <boost/qvm/vec_operations.hpp>
#include <stdint.h>
#include "math_constants.h"
#include "DataTypes.h"
#include "StabilityAnalytically.h"
#include "MD.h"
void taskSRMPIi(MD *MDC);

class StabilityAnalytically;
class CalcBorder
{
public:
    CalcBorder();
    virtual ~CalcBorder();
    void createPoints1D();
    void createPoints2D();
    void createPoints3D();
    void createPoints3DLineBCC();
    void createIMatrix();
    void checkStability();
    void checkStabilityMD();
    void checkStability_MPI();
    void findBorder();
    void CheckBorderMD();
    void preciseBorderCreate();
    void preciseBorderCalculate();
    void checkStabilityMPIMD(boost::qvm::vec<double,3> *Pvar, StabilityPointType *Pdvar, boost::qvm::mat<double,3,3> *Svar, const uint_fast32_t &Nvar);
    void taskSRMPI0();

    //void preciseBorderCalculateMPI();
    void preciseBorderReCreate();
    //TaskType getNextBorderTask();

    TaskType3D getNextTask();
    void stopMPI();


    StabilityAnalytically *SA;
    MD *MDC;
    boost::qvm::vec<double,3> *P, *Pb, *Pbc, *PbcLn, *PS, *PbMD;
    StabilityPointType *Pdata, *Pbdata, *Pbcdata, *PdataR, *PbMDdata;
    uint_fast32_t *I, *IMD;
    boost::qvm::mat<double,3,3> *S, *Sbc, *SbMD;
    boost::qvm::vec<int_fast32_t,2> *PbcL;
    boost::qvm::vec<int_fast32_t,3> *M, *Pi;
    boost::qvm::mat<double,3,3> dG0;
    boost::qvm::vec<double,3> e[2], eStep, rCenter, _1d_n;
    boost::qvm::vec<int_fast32_t,3> n, nCenter;

    uint_fast32_t N, NM, iCenter, NbL, NbLAll, NbcL, NSR, NStep, NStepA, NStepB;
    uint_fast32_t NSended, NReceived;
    double AA, BorderPrecision, BPBP;
    double PotEnergy;

protected:

private:
};

#endif // CALCBORDER_H
