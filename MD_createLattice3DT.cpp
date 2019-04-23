#include "MD.h"
#include "lattice_constans.h"

void MD::createLattice3DT(const boost::qvm::mat<double,3,3> &dG)
{
    n.a[0] = uint_fast32_t(3.0*P_a_cut/(LV[0].a[0]*dG.a[0][0]));
    n.a[1] = uint_fast32_t(3.0*P_a_cut/(LV[1].a[1]*dG.a[1][1]));
    n.a[2] = uint_fast32_t(3.0*P_a_cut/(LV[2].a[2]*dG.a[2][2]));
//std::cerr<<3.0*P_a_cut<<" "<<LV[0].a[0]<<" "<<dG.a[0][0]<<" "<<n.a[2]<<"\n";
    n.a[0] += n.a[0]%2;
    n.a[1] += n.a[1]%2;
    n.a[2] += n.a[2]%2;

    if(n.a[0]<10)
        n.a[0]=8;
    if(n.a[1]<10)
        n.a[1]=8;
    if(n.a[2]<10)
        n.a[2]=8;
    //ns.a[0] = n.a[0]+2*uint_fast32_t(P_a_cut/LV[0].a[0])+1;
    //ns.a[1] = n.a[1]+2*uint_fast32_t(P_a_cut/LV[1].a[1])+1;
    PS[0] = LV[0]*double(n.a[0]);
    PS[1] = LV[1]*double(n.a[1]);
    PS[2] = LV[2]*double(n.a[2]);
    SC = MC_1d3*(PS[0]+PS[1]+PS[2]);
    //LVd2[2] =
    //_1d_LVd2[2]
    N = n.a[0]*n.a[1]*n.a[2]/2;
    //std::cerr<<"Q3 "<<N<<"\n";
    _1d_N = 1.0/double(N);
    if(R==nullptr || (n0.a[0]!=n.a[0]) || (n0.a[1]!=n.a[1]) || (n0.a[2]!=n.a[2]) )
    {
        delete []R;R = new boost::qvm::vec<double,3>[N];
        delete []R0;R0 = new boost::qvm::vec<double,3>[N];
        delete []Ri;Ri = new boost::qvm::vec<int_fast32_t,3>[N];
        delete []V;V = new boost::qvm::vec<double,3>[N];
        delete []F;F = new boost::qvm::vec<double,3>[N];
        delete []I;I = new uint_fast32_t[2*N];
        n0 = n;
        //VP = new uint_fast8_t[N];
        //std::cerr<<"reCreate\n";
    }
    //memset(VP,0,sizeof(uint_fast8_t)*N);

    int_fast32_t i, j, k, iN = 0;
    boost::qvm::vec<double,3> rCenter = 0.5*(PS[0]+PS[1]+PS[2]);
    double drrCenterMin = 1e99, drrCenter;

    for(k=0; k<n.a[2]; ++k)
    for(j=0; j<n.a[1]; ++j)
    {
        for(i=0; i<n.a[0]; ++i)
        {
            if((i%2+j%2+k%2)%2)continue;
            R[iN] = double(i)*LV[0] + double(j)*LV[1] + double(k)*LV[2];
            Ri[iN].a[0] = i;
            Ri[iN].a[1] = j;
            Ri[iN].a[2] = k;
            //if(iN==31)std::cerr<<iN<<" "<<Ri[i].a[0]<<" "<<Ri[i].a[1]<<" "<<Ri[i].a[2]<<"\n";


            I[Ri[iN].a[0]+Ri[iN].a[1]*n.a[0]+Ri[iN].a[2]*n.a[0]*n.a[1]] = iN;
            drrCenter = boost::qvm::mag_sqr(rCenter-R[iN]);
            if(drrCenter<drrCenterMin)
            {
                drrCenterMin = drrCenter;
                nCenter = Ri[iN];

            }
            //std::cerr<<"PR "<<iN<<" "<<i<<" "<<j<<" "<<k<<" "<<R[iN].a[0]<<" "<<R[iN].a[1]<<" "<<R[iN].a[2]<<" "<<Ri[iN].a[0]<<" "<<Ri[iN].a[1]<<" "<<Ri[iN].a[2]<<"\n";

            ++iN;
            //VP[iN] = ((R[iN].a[0]-SC.a[0])*);
            //std::cin.get();
        }
    }
    iCenter = I[nCenter.a[0]+nCenter.a[1]*n.a[0]+nCenter.a[2]*n.a[0]*n.a[1]];
    //iCenter = n.a[0]*n.a[1]*(n.a[2]/2)+n.a[0]*(n.a[1]/2)+n.a[0]/2;
    //boost::qvm::vec<double,2> R1=R[n.a[0]-1]-R[0], R2=R[n.a[0]*(n.a[1]-1)]-R[0];
    Vol = boost::qvm::dot(PS[0],boost::qvm::cross(PS[1],PS[2]));
    //V = R1.a[0]*R2.a[1]-R1.a[1]*R2.a[0];
    //V=boost::qvm::mag(boost::qvm::cross(R[n.a[0]*(n.a[1]-1)]-R[0],));
    _1d_Vol = 1.0/Vol;
    ///eXYMax = 2.0*LV[1].a[1]/LV[0].a[0];
    //std::cerr<<"PS "<<PS[0].a[0]<<" "<<PS[0].a[1]<<" "<<PS[0].a[2]<<" "<<PS[1].a[0]<<" "<<PS[1].a[1]<<" "<<PS[1].a[2]<<" "<<PS[2].a[0]<<" "<<PS[2].a[1]<<" "<<PS[2].a[2]<<"\n";
    //std::cerr<<"NumberLatticeMD "<<N<<" "<<iN<<" "<<Vol<<" "<<iCenter<<" "<<R[iCenter].a[0]<<" "<<R[iCenter].a[1]<<" "<<R[iCenter].a[2]<<" "<<PS[0].a[0]<<" "<<PS[1].a[1]<<" "<<PS[2].a[2]<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<n.a[2]<<"\n";
    //std::cin.get();
    //std::cin.get();
}
