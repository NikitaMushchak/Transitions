#ifndef DATATYPES_H_INCLUDED
#define DATATYPES_H_INCLUDED
#include <boost/qvm/vec.hpp>
#include <boost/qvm/mat.hpp>
struct StabilityPointType
{
    uint_fast32_t StabilitySteps;
    double StabilityTime;
    int_fast32_t PointType, i;
    uint_fast8_t Stability;
    double PotEnergy;
};

struct ElasticModules2D
{
    double Ex,Ey,nuxy,nuyx,G,K;
};

struct ElasticModules3D
{
    double Ex,Ey,Ez,nuxy,nuxz,nuyx,nuyz,nuzx,nuzy,Gyz,Gzx,Gxy,K;
};

struct TaskType2D
{
    boost::qvm::mat<double,2,2> D;
    double Ek;
    int_fast32_t Step, i;
};

struct TaskType3D
{
    boost::qvm::mat<double,3,3> D;
    double Ek;
    int_fast32_t Step, i;
};

struct MacroParametrsType
{
    double Ek, Ep;
};

struct InitialConditionType
{
    double Ek;
};

#endif // DATATYPES_H_INCLUDED
