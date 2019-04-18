#include "MD.h"
#include "phys_constants.h"

void MD::deformLattice(const boost::qvm::mat<double,3,3> &dG)
{
    Vol *= boost::qvm::determinant(dG);
    _1d_Vol = 1.0/Vol;
    //std::cerr<<"dG "<<dG.a[0][0]-1<<" "<<dG.a[0][1]<<" "<<dG.a[1][0]<<" "<<dG.a[1][1]-1<<" "<<(dG.a[0][0]*dG.a[1][1]-dG.a[0][1]*dG.a[1][0])-1<<"\n";
    double Norm, _1d_Norm;
    PS[0] = dG*PS[0];
    PS[1] = dG*PS[1];
    PS[2] = dG*PS[2];
    for(int_fast32_t i=0; i<N; ++i)
    {
        R[i]=dG*R[i];
        //std::cerr<<R[i].a[0]<<" "<<R[i].a[1]<<"\n";
    }
}

void MD::saveLattice()
{
    Vol0 = Vol;
    memcpy(R0,R,N*sizeof(boost::qvm::vec<double,3>));
    memcpy(PS0,PS,3*sizeof(boost::qvm::vec<double,3>));
}

void MD::restoreLattice()
{
    memcpy(R,R0,N*sizeof(boost::qvm::vec<double,3>));
    memcpy(PS,PS0,3*sizeof(boost::qvm::vec<double,3>));
    Vol = Vol0;
    _1d_Vol = 1.0/Vol;
}
