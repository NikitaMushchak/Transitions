#include "MD.h"
#include "phys_constants.h"

MD::MD()
{
    //ctor
    double a0 = MC_s2;
    LV[0].a[0] = 0.5;
    LV[0].a[1] = 0;
    LV[0].a[2] = 0;
    LV[0]*=a0;

    LV[1].a[0] = 0;
    LV[1].a[1] = 0.5;
    LV[1].a[2] = 0;
    LV[1]*=a0;

    LV[2].a[0] = 0;
    LV[2].a[1] = 0;
    LV[2].a[2] = 0.5;
    LV[2]*=a0;

    n.a[0] = 24;
    IC.Ek = 1e-6;
    boost::qvm::set_zero(n0);


    R = nullptr;
    R0 = nullptr;
    Ri = nullptr;
    F = nullptr;
    V = nullptr;
    I = nullptr;
    //VP = nullptr;
    M = nullptr;

}

MD::~MD()
{
    delete []R;
    R = nullptr;
    delete []R0;
    R0 = nullptr;
    delete []Ri;
    Ri = nullptr;
    delete []F;
    F = nullptr;
    delete []V;
    V = nullptr;
    delete []I;
    I = nullptr;
    //delete[] VP;
    //VP = nullptr;
    delete []M;
    M=nullptr;
    //dtor
}



