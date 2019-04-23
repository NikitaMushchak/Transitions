#include "math_constants.h"
#include "phys_constants.h"

boost::qvm::mat<double,2,2> tens(const boost::qvm::vec<double,2> &a, const boost::qvm::vec<double,2> &b)
{
    boost::qvm::mat<double,2,2> temp;
    temp.a[0][0] = a.a[0]*b.a[0];
    temp.a[0][1] = a.a[0]*b.a[1];
    temp.a[1][0] = a.a[1]*b.a[0];
    temp.a[1][1] = a.a[1]*b.a[1];
    //std::cerr<<"T "<<temp.a[0][0]<<" "<<temp.a[0][1]<<" "<<temp.a[1][0]<<" "<<temp.a[1][1]<<"\n";
    return temp;
}

boost::qvm::mat<double,4,4> tens(const boost::qvm::mat<double,2,2> &a, const boost::qvm::mat<double,2,2> &b)
{
    boost::qvm::mat<double,4,4> temp;
    temp.a[0][0] = a.a[0][0]*b.a[0][0];
    temp.a[0][1] = a.a[0][0]*b.a[0][1];
    temp.a[0][2] = a.a[0][0]*b.a[1][0];
    temp.a[0][3] = a.a[0][0]*b.a[1][1];

    temp.a[1][0] = a.a[0][1]*b.a[0][0];
    temp.a[1][1] = a.a[0][1]*b.a[0][1];
    temp.a[1][2] = a.a[0][1]*b.a[1][0];
    temp.a[1][3] = a.a[0][1]*b.a[1][1];

    temp.a[2][0] = a.a[1][0]*b.a[0][0];
    temp.a[2][1] = a.a[1][0]*b.a[0][1];
    temp.a[2][2] = a.a[1][0]*b.a[1][0];
    temp.a[2][3] = a.a[1][0]*b.a[1][1];

    temp.a[3][0] = a.a[1][1]*b.a[0][0];
    temp.a[3][1] = a.a[1][1]*b.a[0][1];
    temp.a[3][2] = a.a[1][1]*b.a[1][0];
    temp.a[3][3] = a.a[1][1]*b.a[1][1];
    return temp;
}

boost::qvm::mat<double,2,2> dotdot(const boost::qvm::mat<double,4,4> &a, const boost::qvm::mat<double,2,2> &b)
{
    boost::qvm::mat<double,2,2> temp;
    temp.a[0][0] = a.a[0][0]*b.a[0][0] + a.a[0][1]*b.a[0][1] + a.a[0][2]*b.a[1][0] + a.a[0][3]*b.a[1][1];
    temp.a[0][1] = a.a[1][0]*b.a[0][0] + a.a[1][1]*b.a[0][1] + a.a[1][2]*b.a[1][0] + a.a[1][3]*b.a[1][1];
    temp.a[1][0] = a.a[2][0]*b.a[0][0] + a.a[2][1]*b.a[0][1] + a.a[2][2]*b.a[1][0] + a.a[2][3]*b.a[1][1];
    temp.a[1][1] = a.a[3][0]*b.a[0][0] + a.a[3][1]*b.a[0][1] + a.a[3][2]*b.a[1][0] + a.a[3][3]*b.a[1][1];
    return temp;
}

void toVoigtNotationC2D(const boost::qvm::mat<double,4,4> &a, boost::qvm::mat<double,3,3> &b)
{
    boost::qvm::set_zero(b);
    for(uint_fast8_t i=0;i<4;++i)
    {
        for(uint_fast8_t j=0;j<4;++j)
        {
            b.a[VoigtNotation2D[i][j].i][VoigtNotation2D[i][j].j] += VoigtNotation2D[i][j].k*a.a[i][j];
        }
    }
}

void fromVoigtNotationC2D(const boost::qvm::mat<double,3,3> &a, boost::qvm::mat<double,4,4> &b)
{
    boost::qvm::set_zero(b);
    for(uint_fast8_t i=0;i<4;++i)
    {
        for(uint_fast8_t j=0;j<4;++j)
        {
            b.a[i][j] = a.a[VoigtNotation2D[i][j].i][VoigtNotation2D[i][j].j];
        }
    }
}

void toVoigtNotationS2D(const boost::qvm::mat<double,4,4> &a, boost::qvm::mat<double,3,3> &b)
{
    boost::qvm::set_zero(b);
    for(uint_fast8_t i=0;i<4;++i)
    {
        for(uint_fast8_t j=0;j<4;++j)
        {
            b.a[VoigtNotation2D[i][j].i][VoigtNotation2D[i][j].j] += a.a[i][j];
        }
    }
}

void fromVoigtNotationS2D(const boost::qvm::mat<double,3,3> &a, boost::qvm::mat<double,4,4> &b)
{
    boost::qvm::set_zero(b);
    for(uint_fast8_t i=0;i<4;++i)
    {
        for(uint_fast8_t j=0;j<4;++j)
        {
            b.a[i][j] = (1.0-0.5*uint_fast8_t(i%3>0))*(1.0-0.5*uint_fast8_t(j%3>0))*a.a[VoigtNotation2D[i][j].i][VoigtNotation2D[i][j].j];
        }
    }
}

void toElastic2D(const boost::qvm::mat<double,3,3> &a, ElasticModules2D &b)
{
    boost::qvm::mat<double,3,3> D = boost::qvm::inverse(a);
    b.Ex = 1.0/(D.a[0][0]);
    b.Ey = 1.0/D.a[1][1];
    b.nuxy = -D.a[1][0]/D.a[0][0];
    b.nuyx = -D.a[0][1]/D.a[1][1];
    b.G = 1.0/D.a[2][2];
    b.K = 0.25*(a.a[0][0]+a.a[1][1]+a.a[0][1]+a.a[1][0]);
}

