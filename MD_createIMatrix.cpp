#include "MD.h"
//#include <algorithm>

void MD::createIMatrix(const bool &recreate)
{
    boost::qvm::vec<double,3> ri,rj;
    boost::qvm::vec<int_fast32_t,3> tempMi, rI, rJ;
    double rmm, MMij, MMii, MMjj;
    NM=0;
    uint_fast32_t i=0, j=0, tempNM, calc;
    //for(uint8_t i=0;i<4;++i)
    if(M==nullptr)M = new boost::qvm::vec<int_fast32_t,3>[N];
    else if(recreate){delete []M;M = new boost::qvm::vec<int_fast32_t,3>[N];}
    //std::cerr<<"M0 "<<R[iCenter].a[0]<<" "<<R[iCenter].a[1]<<" "<<iCenter<<" "<<N<<" "<<n.a[0]<<" "<<n.a[1]<<"\n";
    for(uint_fast32_t i=0; i<N; ++i)
    {
        //if(i==iCenter)continue;
        ri = R[i]-R[iCenter];
        rmm = boost::qvm::mag_sqr(ri);
        //std::cerr<<"M1 "<<i<<" "<<r.a[0]<<" "<<r.a[1]<<" "<<rmm<<" "<<AA<<"\n";
        if(rmm<P_aa_cut && i!=iCenter)// && rmm>1e-10)
        {
            //M[NM].a[2] = int_fast32_t(i*_1d_n.a[0]*_1d_n.a[1])-int_fast32_t(nCenter.a[2]);
            //M[NM].a[1] = int_fast32_t((i-M[NM].a[2]*n.a[1])*_1d_n.a[0])-int_fast32_t(nCenter.a[1]);
            //M[NM].a[0] = int_fast32_t(i-M[NM].a[2]*n.a[1]*n.a[0]-M[NM].a[1]*n.a[0])-int_fast32_t(nCenter.a[0]);
            M[NM] = Ri[i]-Ri[iCenter];
            //std::cerr<<"M "<<NM<<" "<<M[NM].a[0]<<" "<<M[NM].a[1]<<" "<<M[NM].a[2]<<" "<<ri.a[0]<<" "<<ri.a[1]<<" "<<ri.a[2]<<" "<<rmm<<" "<<iCenter<<" "<<P_aa_cut<<"\n";
            ++NM;
        }
    }
    //std::cerr<<"M "<<NM<<" "<<P_aa_cut<<"\n";
    NMd2 = (NM)/2;

    boost::qvm::vec<int_fast32_t,3> *tempM = new boost::qvm::vec<int_fast32_t,3>[NM];
//std::cerr<<"Q2 "<<NM<<"\n";tempd
    tempNM = 0;
    calc = 0;
    for(i=0; i<NM; ++i)
    {
        for(j=0; j<tempNM; ++j)
        {
            //if(i==j)continue;
            rI = Ri[iCenter]+M[i];
            rJ = Ri[iCenter]+tempM[j];
            ri = R[I[rI.a[0]+rI.a[1]*n.a[0]+rI.a[2]*n.a[0]*n.a[1]]]-R[iCenter];
            rj = R[I[rJ.a[0]+rJ.a[1]*n.a[0]+rJ.a[2]*n.a[0]*n.a[1]]]-R[iCenter];
            MMij = boost::qvm::dot(ri,rj);
            MMii = boost::qvm::mag_sqr(ri);
            MMjj = boost::qvm::mag_sqr(rj);
            //std::cerr<<"M0 "<<i<<" "<<tempNM<<" "<<ri.a[0]<<" "<<ri.a[1]<<" "<<rj.a[0]<<" "<<rj.a[1]<<"\n";
            //std::cerr<<"M0 "<<i<<" "<<tempNM<<" "<<ri.a[0]<<" "<<ri.a[1]<<" "<<rj.a[0]<<" "<<rj.a[1]<<"\n";
            if(fabs(MMij*MMij-MMii*MMjj)<1e-10 && fabs(MMii-MMjj)<1e-10)
            {
                //std::cerr<<"M0 "<<i<<" "<<tempNM<<" "<<ri.a[0]<<" "<<ri.a[1]<<" "<<rj.a[0]<<" "<<rj.a[1]<<"\n";
                ++calc;
                goto MD_createIMatrix_Next;
            }
        }
        tempM[tempNM] = M[i];
        tempM[tempNM+NMd2] = -M[i];
        //std::cerr<<"M "<<i<<" "<<tempNM<<" "<<tempM[tempNM].a[0]<<" "<<tempM[tempNM].a[1]<<" "<<tempM[tempNM].a[2]<<"\n";
        ++tempNM;
MD_createIMatrix_Next:
        ;
    }
    delete []M;
    M = tempM;
    if(calc!=NMd2 || NM%2==1)
    {
        std::cerr<<"Error!!! M "<<NM<<" "<<calc<<" "<<tempNM<<" "<<P_a_cut<<" "<<P_aa_cut<<"\n";
        //std::cin.get();
    }
    //std::sort(&M[0],&M[NM-1]);
    //memcpy(&M[NM],&M[0],(NMd2+1)*sizeof(int_fast32_t));
    //memcpy(&M[0],&M[NMd2+1],NMd2*sizeof(int_fast32_t));
    //memcpy(&M[NMd2],&M[NM],(NMd2+1)*sizeof(int_fast32_t));
    //for(uint_fast32_t i=0; i<NM; ++i)
    //    std::cerr<<"M "<<i<<" "<<M[i].a[0]<<" "<<M[i].a[1]<<" "<<M[i].a[2]<<"\n";
    //for(uint_fast32_t i=0; i<NM; ++i)
    //    if(abs(M[i].a[0])>N || abs(M[i].a[1])>N || abs(M[i].a[2])>N)std::cerr<<"M "<<i<<" "<<M[i].a[0]<<" "<<M[i].a[1]<<" "<<M[i].a[2]<<"\n";

    //std::cin.get();
}
