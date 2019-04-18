#include "math_constants.h"
#include "phys_constants.h"

boost::qvm::mat<double,3,3> tens(const boost::qvm::vec<double,3> &a, const boost::qvm::vec<double,3> &b)
{
    boost::qvm::mat<double,3,3> temp;
    temp.a[0][0] = a.a[0]*b.a[0];
    temp.a[0][1] = a.a[0]*b.a[1];
    temp.a[0][2] = a.a[0]*b.a[2];
    temp.a[1][0] = a.a[1]*b.a[0];
    temp.a[1][1] = a.a[1]*b.a[1];
    temp.a[1][2] = a.a[1]*b.a[2];
    temp.a[2][0] = a.a[2]*b.a[0];
    temp.a[2][1] = a.a[2]*b.a[1];
    temp.a[2][2] = a.a[2]*b.a[2];
    //std::cerr<<"T "<<temp.a[0][0]<<" "<<temp.a[0][1]<<" "<<temp.a[1][0]<<" "<<temp.a[1][1]<<"\n";
    return temp;
}

boost::qvm::mat<double,9,9> tens(const boost::qvm::mat<double,3,3> &a, const boost::qvm::mat<double,3,3> &b)
{
    boost::qvm::mat<double,9,9> temp;
    uint_fast16_t i,j,k,l;
    for(i=0; i<3; ++i)
        for(j=0; j<3; ++j)
            for(k=0; k<3; ++k)
                for(l=0; l<3; ++l)
                    temp.a[3*i+j][3*k+l] = a.a[i][j]*b.a[k][l];

    return temp;
}

boost::qvm::mat<double,3,3> dotdot(const boost::qvm::mat<double,9,9> &a, const boost::qvm::mat<double,3,3> &b)
{

    boost::qvm::mat<double,3,3> temp;
    boost::qvm::set_zero(temp);
    uint_fast16_t i,j,k,l;
    for(i=0; i<3; ++i)
        for(j=0; j<3; ++j)
            for(k=0; k<3; ++k)
                for(l=0; l<3; ++l)
                {
                    temp.a[i][j] += a.a[3*i+j][3*k+l]*b.a[l][k];
                    //std::cerr<<temp.a[i][j]<<" "<<a.a[3*i+j][3*k+l]<<" "<<b.a[l][k]<<"\n";
                }

    return temp;
}

void toVoigtNotationC3D(const boost::qvm::mat<double,9,9> &a, boost::qvm::mat<double,6,6> &b)
{
    boost::qvm::set_zero(b);
    for(uint_fast8_t i=0;i<9;++i)
    {
        for(uint_fast8_t j=0;j<9;++j)
        {
            b.a[VoigtNotation3D[i][j].i][VoigtNotation3D[i][j].j] += VoigtNotation3D[i][j].k*a.a[i][j];
        }
    }
}

void fromVoigtNotationC3D(const boost::qvm::mat<double,6,6> &a, boost::qvm::mat<double,9,9> &b)
{
    boost::qvm::set_zero(b);
    for(uint_fast8_t i=0;i<9;++i)
    {
        for(uint_fast8_t j=0;j<9;++j)
        {
            b.a[i][j] = a.a[VoigtNotation3D[i][j].i][VoigtNotation3D[i][j].j];
        }
    }
}

void toVoigtNotationS3D(const boost::qvm::mat<double,9,9> &a, boost::qvm::mat<double,6,6> &b)
{
    boost::qvm::set_zero(b);
    for(uint_fast8_t i=0;i<9;++i)
    {
        for(uint_fast8_t j=0;j<9;++j)
        {
            b.a[VoigtNotation3D[i][j].i][VoigtNotation3D[i][j].j] += a.a[i][j];
        }
    }
}

void fromVoigtNotationS3D(const boost::qvm::mat<double,6,6> &a, boost::qvm::mat<double,9,9> &b)
{
    boost::qvm::set_zero(b);
    for(uint_fast8_t i=0;i<9;++i)
    {
        for(uint_fast8_t j=0;j<9;++j)
        {
            b.a[i][j] = (1.0-0.5*uint_fast8_t(i%4>0))*(1.0-0.5*uint_fast8_t(j%4>0))*a.a[VoigtNotation3D[i][j].i][VoigtNotation3D[i][j].j];
        }
    }
}

void toElastic3D(const boost::qvm::mat<double,6,6> &a, ElasticModules3D &b)
{
    boost::qvm::mat<double,6,6> D = boost::qvm::inverse(a);
    b.Ex = 1.0/D.a[0][0];
    b.Ey = 1.0/D.a[1][1];
    b.Ez = 1.0/D.a[2][2];

    b.nuxy = -D.a[1][0]/D.a[0][0];
    b.nuxz = -D.a[2][0]/D.a[0][0];
    b.nuyx = -D.a[0][1]/D.a[1][1];
    b.nuyz = -D.a[2][1]/D.a[1][1];
    b.nuzx = -D.a[0][2]/D.a[2][2];
    b.nuzy = -D.a[1][2]/D.a[2][2];

    b.Gyz = 1.0/D.a[3][3];
    b.Gzx = 1.0/D.a[4][4];
    b.Gxy = 1.0/D.a[5][5];
    b.K = MC_1d9*(a.a[0][0]+a.a[0][1]+a.a[0][2]+a.a[1][0]+a.a[1][1]+a.a[1][2]+a.a[2][0]+a.a[2][1]+a.a[2][2]);
}

