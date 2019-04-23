#include "CalcBorder.h"
#include <iomanip>
CalcBorder::CalcBorder()
{
    //ctor
    SA = nullptr;
    MDC = nullptr;
    boost::qvm::X(e[0])= -0.1;
    boost::qvm::Y(e[0])= boost::qvm::X(e[0]);//-0.6;
    boost::qvm::Z(e[0])= boost::qvm::X(e[0])*(MC_s6 - 3.)/(MC_2s3 - 3.);//-0.6; // -0.6
    boost::qvm::X(e[1])= 0.5;
    boost::qvm::Y(e[1])= boost::qvm::X(e[1]);//-0.1;
    boost::qvm::Z(e[1])= boost::qvm::X(e[1])*(MC_s6 - 3.)/(MC_2s3 - 3.);// -0.1;
    boost::qvm::X(n) = 150;
    boost::qvm::Y(n) = 150;//2*uint_fast32_t(X(n)*MC_1ds3);
    boost::qvm::Z(n) = 150;//2*uint_fast32_t(X(n)*MC_1ds3);
    BorderPrecision = 1e-2;
    NStepA = 100;
	NStepB = 1000;
    BPBP = BorderPrecision*BorderPrecision;
    ///SA = new StabilityAnalytically();
    MDC = new MD();
    P = nullptr;
    I = nullptr;
    M = nullptr;
    Pdata = nullptr;
    //SA = NULL;
    Pb = nullptr;
    Pbdata = nullptr;
    Pbc = nullptr;
    Pbcdata = nullptr;
    PbcL = nullptr;
    PbcLn = nullptr;
    S = nullptr;
    Sbc = nullptr;
}

CalcBorder::~CalcBorder()
{
    //dtor
    //std::cerr<<"W1\n";
    delete[] P;
    P = NULL;
    //std::cerr<<"W2\n";
    delete[] I;
    I = NULL;
    //std::cerr<<"W3\n";
    delete[] M;
    M=NULL;
    //std::cerr<<"W4\n";
    delete []Pdata;
    Pdata = NULL;
    //std::cerr<<"W5\n";
    delete SA;
    SA = NULL;
    //std::cerr<<"W6\n";
    delete[] Pb;
    Pb = NULL;
    //std::cerr<<"W7\n";
    delete[] Pbdata;
    Pbdata = NULL;
    //std::cerr<<"W8\n";
    delete[] Pbc;
    Pbc = NULL;
    //std::cerr<<"W9\n";
    delete[] Pbcdata;
    Pbcdata = NULL;
    //std::cerr<<"W10\n";
    delete SA;
    SA = NULL;
    //std::cerr<<"W11\n";
    delete[] PbcL;
    PbcL = NULL;
    //std::cerr<<"W12\n";
    delete[] PbcLn;
    PbcLn = NULL;
    //std::cerr<<"W13\n";
    delete[] S;
    S = NULL;
    //std::cerr<<"W14\n";
    delete[] Sbc;
    Sbc = NULL;
    //std::cerr<<"W15\n";
    delete MDC;
    MDC = NULL;
    //std::cerr<<"W16\n";
}

void CalcBorder::createPoints3D()
{
    boost::qvm::X(eStep) = (boost::qvm::X(e[1])-boost::qvm::X(e[0]))/double(boost::qvm::X(n));
    boost::qvm::Y(eStep) = boost::qvm::X(eStep);
    //boost::qvm::Z(eStep) = boost::qvm::X(eStep);
    boost::qvm::Z(eStep) = (boost::qvm::Z(e[1])-boost::qvm::Z(e[0]))/double(boost::qvm::Z(n));
    boost::qvm::Y(n) = (boost::qvm::Y(e[1])-boost::qvm::Y(e[0]))/boost::qvm::Y(eStep);
    //boost::qvm::Z(n) = (boost::qvm::Z(e[1])-boost::qvm::Z(e[0]))/boost::qvm::Z(eStep)+1;
    //N = 4*boost::qvm::X(n)*boost::qvm::Y(n)*boost::qvm::Z(n);
    n*=2;
    n+=MC_i1XYZV3;
    N = boost::qvm::X(n)*boost::qvm::Y(n)*boost::qvm::Z(n)/2+1;
    // N = boost::qvm::X(n);
    //std::cerr<<"C "<<N<<" "<<eStep.a[0]<<" "<<eStep.a[1]<<" "<<eStep.a[2]<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<n.a[2]<<"\n";
    //std::cin.get();
    _1d_n.a[0] = 1.0/double(n.a[0]);
    _1d_n.a[1] = 1.0/double(n.a[1]);
    _1d_n.a[2] = 1.0/double(n.a[2]);
    //AA = boost::qvm::mag_sqr(eStep);
    AA = boost::qvm::X(eStep)*boost::qvm::X(eStep)*0.5*2.1*2.1;
    //NL = 2*boost::qvm::Y(n)*(boost::qvm::X(n)-1)+2*(boost::qvm::Y(n)-1)*(2*boost::qvm::X(n)-1)+10;
    std::cerr<<"Number "<<N<<"\n";
    P = new boost::qvm::vec<double,3>[N];
    Pi = new boost::qvm::vec<int_fast32_t,3>[N];
    S = new boost::qvm::mat<double,3,3>[N];
    Pdata = new StabilityPointType[N];
    I = new uint_fast32_t[2*N];
    memset(Pdata,0,N*sizeof(StabilityPointType));
    //I = new uint_fast32_t[NL];
    boost::qvm::vec<double,3> rCenter = 0.5*(e[0]+e[1]);
    double drrCenterMin = 1e99, drrCenter;
    uint_fast32_t i, j, k, iN = 0;
    for(k=0; k<n.a[2]; ++k)
        for(j=0; j<n.a[1]; ++j)
            for(i=0; i<n.a[0]; ++i)
            {
                if((i%2+j%2+k%2)%2)
                    continue;

                P[iN].a[0] = e[0].a[0]+i*0.5*eStep.a[0];
                P[iN].a[1] = e[0].a[1]+j*0.5*eStep.a[1];
                P[iN].a[2] = e[0].a[2]+k*0.5*eStep.a[2];
                Pi[iN].a[0] = i;
                Pi[iN].a[1] = j;
                Pi[iN].a[2] = k;
                I[Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]] = iN;
                drrCenter = boost::qvm::mag_sqr(rCenter-P[iN]);
                if(drrCenter<drrCenterMin)
                {
                    drrCenterMin = drrCenter;
                    nCenter = Pi[iN];

                }
                //std::cerr<<"C "<<iN<<" "<<i<<" "<<j<<" "<<k<<" "<<P[iN].a[0]<<" "<<P[iN].a[1]<<" "<<P[iN].a[2]<<" "<<Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]<<" "<<I[Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]]<<"\n";
                //<<" "<<ParticleCoordinate[Particles.ParticlesNumber].y
                //<<" "<<ParticleCoordinate[Particles.ParticlesNumber].z<<"\n";
                //std::cin.get();
                ++iN;
            }

    nCenter.a[0] -= nCenter.a[0]%2;
    nCenter.a[1] -= nCenter.a[1]%2;
    nCenter.a[2] -= nCenter.a[2]%2;
    iCenter = I[nCenter.a[0]+nCenter.a[1]*n.a[0]+nCenter.a[2]*n.a[0]*n.a[1]];
    //std::cerr<<"C "<<N<<" "<<eStep.a[0]<<" "<<eStep.a[1]<<" "<<eStep.a[2]<<"\n";
    std::cerr<<"C "<<N<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<n.a[2]<<" "<<iCenter<<" "<<P[iCenter].a[0]<<" "<<P[iCenter].a[1]<<" "<<P[iCenter].a[2]<<"\n";
    //std::cin.get();
}
void CalcBorder::createPoints3DLineBCC()
{
    boost::qvm::X(eStep) = (boost::qvm::X(e[1])-boost::qvm::X(e[0]))/double(boost::qvm::X(n));
    boost::qvm::Y(eStep) = (boost::qvm::Y(e[1])-boost::qvm::Y(e[0]))/double(boost::qvm::Y(n));
    boost::qvm::Z(eStep) = (boost::qvm::Z(e[1])-boost::qvm::Z(e[0]))/double(boost::qvm::Z(n));//boost::qvm::X(eStep);
    // boost::qvm::Z(eStep) = (boost::qvm::Z(e[1])-boost::qvm::Z(e[0]))/double(boost::qvm::Z(n));
    boost::qvm::Y(n) = (boost::qvm::Y(e[1])-boost::qvm::Y(e[0]))/boost::qvm::Y(eStep);
    //boost::qvm::Z(n) = (boost::qvm::Z(e[1])-boost::qvm::Z(e[0]))/boost::qvm::Z(eStep)+1;
    //N = 4*boost::qvm::X(n)*boost::qvm::Y(n)*boost::qvm::Z(n);
    n*=2;
    n+=MC_i1XYZV3;
    // N = boost::qvm::X(n)*boost::qvm::Y(n)*boost::qvm::Z(n)/2+1;
    N = boost::qvm::X(n);
    //std::cerr<<"C "<<N<<" "<<eStep.a[0]<<" "<<eStep.a[1]<<" "<<eStep.a[2]<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<n.a[2]<<"\n";
    //std::cin.get();
    _1d_n.a[0] = 1.0/double(n.a[0]);
    _1d_n.a[1] = 1.0/double(n.a[1]);
    _1d_n.a[2] = 1.0/double(n.a[2]);
    //AA = boost::qvm::mag_sqr(eStep);
    AA = boost::qvm::X(eStep)*boost::qvm::X(eStep)*0.5*2.1*2.1;
    //NL = 2*boost::qvm::Y(n)*(boost::qvm::X(n)-1)+2*(boost::qvm::Y(n)-1)*(2*boost::qvm::X(n)-1)+10;
    std::cerr<<"Number "<<N<<"\n";
    P = new boost::qvm::vec<double,3>[N];
    Pi = new boost::qvm::vec<int_fast32_t,3>[N];
    S = new boost::qvm::mat<double,3,3>[N];
    Pdata = new StabilityPointType[N];
    I = new uint_fast32_t[2*N];
    memset(Pdata,0,N*sizeof(StabilityPointType));
    //I = new uint_fast32_t[NL];
    boost::qvm::vec<double,3> rCenter = 0.5*(e[0]+e[1]);
    double drrCenterMin = 1e99, drrCenter;
    uint_fast32_t i, j, k, iN = 0;
    // for(k=0; k<n.a[2]; ++k)
    //     for(j=0; j<n.a[1]; ++j)
    //         for(i=0; i<n.a[0]; ++i)

  std::ofstream myfile;
  myfile.open("line.txt");
  // myfile << "Writing this to a file.\n";
  // myfile.close();

  for(size_t i = 0 ; i< N; ++i){
                // if((i%2+j%2+k%2)%2)
                //     continue;
        P[iN].a[0] = e[0].a[0]+i*0.5*eStep.a[0];
        P[iN].a[1] = e[0].a[1]+i*0.5*eStep.a[1];
        P[iN].a[2] = e[0].a[2]+i*0.5*eStep.a[2];
        Pi[iN].a[0] = i;
        Pi[iN].a[1] = 0;//i;
        Pi[iN].a[2] = 0;//i;
        I[Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]] = iN;
        drrCenter = boost::qvm::mag_sqr(rCenter-P[iN]);
        if(drrCenter<drrCenterMin)
        {
            drrCenterMin = drrCenter;
            nCenter = Pi[iN];

        }
                //std::cerr<<"C "<<iN<<" "<<i<<" "<<j<<" "<<k<<" "<<P[iN].a[0]<<" "<<P[iN].a[1]<<" "<<P[iN].a[2]<<" "<<Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]<<" "<<I[Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]]<<"\n";
        myfile<<P[iN].a[0]<<" "<<P[iN].a[1]<<" "<<P[iN].a[2]<<"\n";//std::endl;
                //<<" "<<ParticleCoordinate[Particles.ParticlesNumber].y
                //<<" "<<ParticleCoordinate[Particles.ParticlesNumber].z<<"\n";
                //std::cin.get();
        ++iN;
    }

    nCenter.a[0] -= nCenter.a[0]%2;
    nCenter.a[1] -= nCenter.a[1]%2;
    nCenter.a[2] -= nCenter.a[2]%2;
    iCenter = I[nCenter.a[0]+nCenter.a[1]*n.a[0]+nCenter.a[2]*n.a[0]*n.a[1]];
    //std::cerr<<"C "<<N<<" "<<eStep.a[0]<<" "<<eStep.a[1]<<" "<<eStep.a[2]<<"\n";
    std::cerr<<"C "<<N<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<n.a[2]<<" "<<iCenter<<" "<<P[iCenter].a[0]<<" "<<P[iCenter].a[1]<<" "<<P[iCenter].a[2]<<"\n";
    //std::cin.get();
    myfile.close();
}

void CalcBorder::createPoints2D()
{
    boost::qvm::X(eStep) = (boost::qvm::X(e[1])-boost::qvm::X(e[0]))/double(boost::qvm::X(n)-1);
    boost::qvm::Y(eStep) = boost::qvm::X(eStep)*MC_s3d2;
    uint_fast32_t Yn = (boost::qvm::Y(e[1])-boost::qvm::Y(e[0]))/boost::qvm::Y(eStep);
    boost::qvm::Y(n) = (Yn%2)?Yn+1:Yn;
    //boost::qvm::Z(eStep) = boost::qvm::X(eStep);
    boost::qvm::Z(eStep) = 0;
    boost::qvm::Z(n) = 1;
    //N = 4*boost::qvm::X(n)*boost::qvm::Y(n)*boost::qvm::Z(n);
    N = boost::qvm::X(n)*boost::qvm::Y(n);
    //std::cerr<<"C "<<N<<" "<<eStep.a[0]<<" "<<eStep.a[1]<<" "<<eStep.a[2]<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<n.a[2]<<"\n";
    //std::cin.get();
    _1d_n.a[0] = 1.0/double(n.a[0]);
    _1d_n.a[1] = 1.0/double(n.a[1]);
    _1d_n.a[2] = 0;
    //AA = boost::qvm::mag_sqr(eStep);
    AA = boost::qvm::X(eStep)*boost::qvm::X(eStep)*1.44;
    //NL = 2*boost::qvm::Y(n)*(boost::qvm::X(n)-1)+2*(boost::qvm::Y(n)-1)*(2*boost::qvm::X(n)-1)+10;
    std::cerr<<"Number "<<N<<"\n";
    P = new boost::qvm::vec<double,3>[N];
    Pi = new boost::qvm::vec<int_fast32_t,3>[N];
    S = new boost::qvm::mat<double,3,3>[N];
    Pdata = new StabilityPointType[N];
    I = new uint_fast32_t[2*N];
    memset(Pdata,0,N*sizeof(StabilityPointType));
    //I = new uint_fast32_t[NL];
    boost::qvm::vec<double,3> rCenter = 0.5*(e[0]+e[1]);
    double drrCenterMin = 1e99, drrCenter;
    boost::qvm::vec<double,2> tmp= {0,0};
    uint_fast32_t i, j, iN = 0;
    for(j=0; j<boost::qvm::Y(n); ++j)
    {
        boost::qvm::X(tmp) = 0.5*boost::qvm::X(eStep)*(j%2);
        for(i=0; i<boost::qvm::X(n); ++i)
        {
            P[iN].a[0] = e[0].a[0]+i*eStep.a[0]+tmp.a[0];
            P[iN].a[1] = e[0].a[1]+j*eStep.a[1]+tmp.a[1];
            P[iN].a[2] = 0;
            Pi[iN].a[0] = i;
            Pi[iN].a[1] = j;
            Pi[iN].a[2] = 0;
            I[Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]] = iN;

            //std::cerr<<"C "<<iN<<" "<<N<<" "<<P[iN].a[0]<<" "<<P[iN].a[1]<<" "<<P[iN].a[2]<<"\n";
            drrCenter = boost::qvm::mag_sqr(rCenter-P[iN]);
            if(drrCenter<drrCenterMin)
            {
                drrCenterMin = drrCenter;
                nCenter = Pi[iN];

            }
            ++iN;
            //std::cerr<<"P "<<iN<<" "<<P[iN].a[0]<<" "<<P[iN].a[1]<<"\n";
        }
    }

    nCenter.a[0] -= nCenter.a[0]%2;
    nCenter.a[1] -= nCenter.a[1]%2;
    nCenter.a[2] -= nCenter.a[2]%2;
    iCenter = I[nCenter.a[0]+nCenter.a[1]*n.a[0]+nCenter.a[2]*n.a[0]*n.a[1]];
    //std::cerr<<"C "<<N<<" "<<eStep.a[0]<<" "<<eStep.a[1]<<" "<<eStep.a[2]<<"\n";
    std::cerr<<"C "<<N<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<n.a[2]<<" "<<iCenter<<" "<<P[iCenter].a[0]<<" "<<P[iCenter].a[1]<<" "<<P[iCenter].a[2]<<"\n";
    //std::cerr<<"C "<<nCenter.a[0]<<" "<<nCenter.a[1]<<" "<<nCenter.a[2]<<" "<<rCenter.a[0]<<" "<<rCenter.a[1]<<" "<<rCenter.a[2]<<"\n";
    //std::cin.get();
}

void CalcBorder::createPoints1D()
{
    boost::qvm::X(eStep) = (boost::qvm::X(e[1])-boost::qvm::X(e[0]))/double(boost::qvm::X(n)-1);
    // boost::qvm::Y(eStep) = boost::qvm::X(eStep)*MC_s3d2;
    // uint_fast32_t Yn = (boost::qvm::Y(e[1])-boost::qvm::Y(e[0]))/boost::qvm::Y(eStep);
    // boost::qvm::Y(n) = (Yn%2)?Yn+1:Yn;
    //boost::qvm::Z(eStep) = boost::qvm::X(eStep);
    boost::qvm::Z(eStep) = 0;
    boost::qvm::Y(eStep) = 0;
    boost::qvm::Z(n) = 1;
    boost::qvm::Y(n) = 1;
    //N = 4*boost::qvm::X(n)*boost::qvm::Y(n)*boost::qvm::Z(n);
    N = boost::qvm::X(n);
    // N = boost::qvm::X(n)*boost::qvm::Y(n);
    //std::cerr<<"C "<<N<<" "<<eStep.a[0]<<" "<<eStep.a[1]<<" "<<eStep.a[2]<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<n.a[2]<<"\n";
    //std::cin.get();
    _1d_n.a[0] = 1.0/double(n.a[0]);
    _1d_n.a[1] = 0;//1.0/double(n.a[1]);
    _1d_n.a[2] = 0;
    //AA = boost::qvm::mag_sqr(eStep);
    AA = boost::qvm::X(eStep)*boost::qvm::X(eStep)*1.44;
    //NL = 2*boost::qvm::Y(n)*(boost::qvm::X(n)-1)+2*(boost::qvm::Y(n)-1)*(2*boost::qvm::X(n)-1)+10;
    std::cerr<<"Number "<<N<<"\n";
    P = new boost::qvm::vec<double,3>[N];
    Pi = new boost::qvm::vec<int_fast32_t,3>[N];
    S = new boost::qvm::mat<double,3,3>[N];
    Pdata = new StabilityPointType[N];
    I = new uint_fast32_t[2*N];
    memset(Pdata,0,N*sizeof(StabilityPointType));
    //I = new uint_fast32_t[NL];
    boost::qvm::vec<double,3> rCenter = 0.5*(e[0]+e[1]);
    double drrCenterMin = 1e99, drrCenter;
    boost::qvm::vec<double,2> tmp= {0,0};
    uint_fast32_t i, j, iN = 0;
    // for(j=0; j<boost::qvm::Y(n); ++j)
    // {
        boost::qvm::X(tmp) = 0.5*boost::qvm::X(eStep)*(j%2);
        for(i=0; i<boost::qvm::X(n); ++i)
        {
            P[iN].a[0] = e[0].a[0]+i*eStep.a[0]+tmp.a[0];
            P[iN].a[1] = 0;//e[0].a[1]+j*eStep.a[1]+tmp.a[1];
            P[iN].a[2] = 0;
            Pi[iN].a[0] = i;
            Pi[iN].a[1] = 0;//j;
            Pi[iN].a[2] = 0;
            I[Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]] = iN;

            //std::cerr<<"C "<<iN<<" "<<N<<" "<<P[iN].a[0]<<" "<<P[iN].a[1]<<" "<<P[iN].a[2]<<"\n";
            drrCenter = boost::qvm::mag_sqr(rCenter-P[iN]);
            if(drrCenter<drrCenterMin)
            {
                drrCenterMin = drrCenter;
                nCenter = Pi[iN];

            }
            ++iN;
            //std::cerr<<"P "<<iN<<" "<<P[iN].a[0]<<" "<<P[iN].a[1]<<"\n";
        }
    // }

    nCenter.a[0] -= nCenter.a[0]%2;
    nCenter.a[1] -= nCenter.a[1]%2;
    nCenter.a[2] -= nCenter.a[2]%2;
    iCenter = I[nCenter.a[0]+nCenter.a[1]*n.a[0]+nCenter.a[2]*n.a[0]*n.a[1]];
    //std::cerr<<"C "<<N<<" "<<eStep.a[0]<<" "<<eStep.a[1]<<" "<<eStep.a[2]<<"\n";
    std::cerr<<"C "<<N<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<n.a[2]<<" "<<iCenter<<" "<<P[iCenter].a[0]<<" "<<P[iCenter].a[1]<<" "<<P[iCenter].a[2]<<"\n";
    //std::cerr<<"C "<<nCenter.a[0]<<" "<<nCenter.a[1]<<" "<<nCenter.a[2]<<" "<<rCenter.a[0]<<" "<<rCenter.a[1]<<" "<<rCenter.a[2]<<"\n";
    //std::cin.get();
}

void CalcBorder::createIMatrix()
{
    boost::qvm::vec<double,3> r;
    M = new boost::qvm::vec<int_fast32_t,3>[200];
    double rmm;
    NM=0;
    //std::cerr<<"M0 "<<P[iCenter].a[0]<<" "<<P[iCenter].a[1]<<" "<<iCenter<<" "<<N<<" "<<n.a[0]<<" "<<n.a[1]<<"\n";
    for(uint_fast32_t i=0; i<N; ++i)
    {
        //if(i==iCenter)continue;
        r = P[i]-P[iCenter];
        rmm = boost::qvm::mag_sqr(r);
        //std::cerr<<"M1 "<<i<<" "<<r.a[0]<<" "<<r.a[1]<<" "<<rmm<<" "<<AA<<"\n";
        if(rmm<AA)// && rmm>1e-10)
        {
            //M[NM].a[2] = int_fast32_t(i*_1d_n.a[0]*_1d_n.a[1])-int_fast32_t(nCenter.a[2]);
            //M[NM].a[1] = int_fast32_t((i-M[NM].a[2]*n.a[1])*_1d_n.a[0])-int_fast32_t(nCenter.a[1]);
            //M[NM].a[0] = int_fast32_t(i-M[NM].a[2]*n.a[1]*n.a[0]-M[NM].a[1]*n.a[0])-int_fast32_t(nCenter.a[0]);
            M[NM] = Pi[i]-nCenter;


            //std::cerr<<"M "<<Pi[i].a[0]<<" "<<Pi[i].a[1]<<" "<<M[NM].a[1]<<" "<<Pi[i].a[2]<<"\n";
            //std::cerr<<"M "<<NM<<" "<<M[NM].a[0]<<" "<<M[NM].a[1]<<" "<<M[NM].a[2]<<" "<<iCenter<<" "<<AA<<"\n";
            ++NM;
            //std::cin.get();
        }
    }
    //std::cin.get();
}

void CalcBorder::checkStability()
{
    SA->startTask();
    char taskdir[256]=".", filename[256]="";
    strcpy(filename, taskdir);
    strcat(filename, "/FCCSpheres.txt");
    std::ofstream Result_file(filename);
    Result_file<<"[";
    for(uint_fast32_t i=0; i<SA->MaxCSDN; ++i)
    {
        Result_file<<" "<<std::setprecision(14)<<SA->CSD[i]<<", ";
    }
    Result_file<<"];";
    Result_file.close();/**/
    //std::cerr<<"W0\n";
    //std::cin.get();
    SA->checkStability3D(P,Pdata,S,N);
    strcpy(filename, "./");
    strcat(filename, "ResultC.txt");
    std::ofstream ResultC_file(filename);
    for(uint_fast32_t i=0; i<N; ++i)
    {
        if(Pdata[i].Stability==2)
        ResultC_file<<i<<" "<<P[i].a[0]<<" "<<P[i].a[1]<<" "<<P[i].a[2]<<" "<<int_fast32_t(Pdata[i].Stability)<<" "<<Pdata[i].StabilitySteps<<" "<<Pdata[i].StabilityTime<<"\n";
    }
    ResultC_file.close();
    findBorder();
    preciseBorderCreate();
    strcpy(filename, "./");
    strcat(filename, "/ResultBc.txt");
    std::ofstream ResultBc_file(filename);
    for(uint_fast32_t i=0;i<NbL;++i)
    {

        ///toElastic2D(CB.Sbc[i],EM);
        //if(EM.Ex<0)
        ResultBc_file<<Pbc[i].a[0]<<" "<<Pbc[i].a[1]<<" "<<Pbc[i].a[2]<<"\n";
    }
    ResultBc_file.close();
    preciseBorderCalculate();
    strcpy(filename, "./");
    strcat(filename, "/ResultBcV.txt");
    std::ofstream ResultBcL_file(filename);
    for(uint_fast32_t i=0;i<NbL;++i)
    {
        ResultBcL_file<<Pbc[i].a[0]<<" "<<Pbc[i].a[1]<<" "<<Pbc[i].a[2]<<"\n";
    }

    ResultBcL_file.close();
}

void CalcBorder::checkStabilityMD()
{
    MDC->createLattice3DT();
    MDC->setMPP();
    MDC->createIMatrix();
    MDC->reachSpecifiedStress();
    MDC->setTemperatureNormal();
    MDC->calculateStateIM();
    int_fast32_t step;
    for(step=0; step<10000; ++step)
    {
        MDC->calculateForcesIM();
        //std::cerr<<"Q10"<<"\n";
        MDC->calculateIncrements();
        std::cerr<<"Q "<<MDC->t<<" "<<MDC->Ek*MDC->_1d_N<<" "<<MDC->Ep<<" "<<MDC->P_Cmax<<" "<<MDC->P_C<<"\n";
        //std::cin.get();
    }
    //MDC->calculateForcesIM();
    //MDC->calculateActualStiffness();
    //MDC->calculateForcesIM();

    //MD->startTask();
    //std::cerr<<"W0\n";
    //MD->checkStability2D(P,Pdata,S,N);
    //std::cerr<<"W1\n";
}

void CalcBorder::findBorder()
{
    boost::qvm::vec<double,3> r;
    boost::qvm::vec<int_fast32_t,3> nj;
    uint_fast32_t j;//,c=0;
    //uint_fast8_t s;
    NbL = 0;
    for(int_fast32_t i=0; i<N; ++i)
    {
        //s=0;
        //std::cerr<<"C0 "<<i<<"\n";
        if(Pdata[i].Stability==2)
            for(uint_fast32_t k=0; k<NM; ++k)
            {
                //ni.a[2] = int_fast32_t(i*_1d_n.a[0]*_1d_n.a[1]);
                //ni.a[1] = int_fast32_t((i-ni.a[2]*n.a[1])*_1d_n.a[0]);
                //ni.a[0] = i-n.a[0]*ni.a[1]-n.a[0]*n.a[1]*ni.a[2];
                nj = Pi[i]+M[k];

                if(nj.a[0]<0 || nj.a[0]>n.a[0]-1)
                    continue;
                if(nj.a[1]<0 || nj.a[1]>n.a[1]-1)
                    continue;
                if(nj.a[2]<0 || nj.a[2]>n.a[2]-1)
                    continue;
                j = I[nj.a[0] + nj.a[1]*n.a[0] + nj.a[2]*n.a[0]*n.a[1]];
                //std::cerr<<"Border1 "<<i<<" "<<nj.a[0]<<" "<<nj.a[1]<<" "<<nj.a[2]<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<j<<" "<<int(Pdata[j].Stability)<<"\n";

                //std::cerr<<"C1 "<<j<<" "<<nj.a[0]<<" "<<n.a[0]<<" "<<nj.a[1]<<" "<<n.a[1]<<"\n";
                //std::cerr<<"Border1 "<<i<<" "<<j<<"\n";
                if(Pdata[j].Stability==1)
                {
                    r = P[j]-P[iCenter];
                    Pdata[i].PointType = 1;
                    Pdata[j].PointType = -1;
                    ++NbL;
                    //if(i==1271 || j==1271)
                    //std::cerr<<"Border "<<N<<" "<<NbL<<" "<<i<<" "<<j<<" "<<Pdata[i].PointType<<" "<<M[k].a[0]<<" "<<M[k].a[1]<<" "<<nj.a[0]<<" "<<nj.a[1]<<" "<<nCenter.a[0]<<" "<<nCenter.a[1]<<" "<<boost::qvm::mag(r)<<"\n";

                    //std::cin.get();
                }
                /*if(I[0]>0)
                {
                    std::cerr<<i<<" "<<j<<" "<<I[0]<<"\n";
                    std::cin.get();
                }/**/
            }
//        if(Pdata[i].PointType==1)++c;
//SkipPoint:;
    }
std::cerr<<"FB "<<NbL<<"\n";
}

void CalcBorder::CheckBorderMD()
{
    PbMD = new boost::qvm::vec<double,3>[N];
    SbMD = new boost::qvm::mat<double,3,3>[N];
    PbMDdata = new StabilityPointType[N];
    IMD = new uint_fast32_t[N];
    boost::qvm::vec<int_fast32_t,3> nj;
    uint_fast32_t i,j;//,c=0;
    //uint_fast8_t s;

    do
    {
        NbL = 0;
        for(i=0; i<N; ++i)
        {
            if(Pdata[i].Stability==2 && Pdata[i].PointType!=2)
                for(uint_fast32_t k=0; k<NM; ++k)
                {
                    //ni.a[2] = int_fast32_t(i*_1d_n.a[0]*_1d_n.a[1]);
                    //ni.a[1] = int_fast32_t((i-ni.a[2]*n.a[1])*_1d_n.a[0]);
                    //ni.a[0] = i-n.a[0]*ni.a[1]-n.a[0]*n.a[1]*ni.a[2];
                    nj = Pi[i]+M[k];

                    if(nj.a[0]<0 || nj.a[0]>n.a[0]-1)
                        continue;
                    if(nj.a[1]<0 || nj.a[1]>n.a[1]-1)
                        continue;
                    if(nj.a[2]<0 || nj.a[2]>n.a[2]-1)
                        continue;
                    j = I[nj.a[0] + nj.a[1]*n.a[0] + nj.a[2]*n.a[0]*n.a[1]];
                    //std::cerr<<"Border1 "<<i<<" "<<nj.a[0]<<" "<<nj.a[1]<<" "<<nj.a[2]<<" "<<n.a[0]<<" "<<n.a[1]<<" "<<j<<" "<<int(Pdata[j].Stability)<<"\n";

                    //std::cerr<<"C1 "<<j<<" "<<nj.a[0]<<" "<<n.a[0]<<" "<<nj.a[1]<<" "<<n.a[1]<<"\n";
                    //std::cerr<<"Border1 "<<i<<" "<<j<<"\n";
                    if(Pdata[j].Stability==1)
                    {

                        PbMD[NbL] = P[i];
                        IMD[NbL] = i;
                        ++NbL;
                        break;
                        //if(i==1271 || j==1271)
                        //std::cerr<<"Border "<<N<<" "<<NbL<<" "<<i<<" "<<j<<" "<<Pdata[i].PointType<<" "<<M[k].a[0]<<" "<<M[k].a[1]<<" "<<nj.a[0]<<" "<<nj.a[1]<<" "<<nCenter.a[0]<<" "<<nCenter.a[1]<<" "<<boost::qvm::mag(r)<<"\n";

                        //std::cin.get();
                    }
                }
//        if(Pdata[i].PointType==1)++c;
//SkipPoint:;
        }
        std::cerr<<"FBMD "<<NbL<<"\n";
        NStepA = 500;
		NStepB = 10000;
        checkStabilityMPIMD(PbMD,PbMDdata,SbMD,NbL);
        for(i=0; i<NbL; ++i)
        {
            Pdata[IMD[i]]=PbMDdata[i];
            S[IMD[i]]=SbMD[i];
            if(PbMDdata[i].Stability==2)
            {
                Pdata[IMD[i]].PointType = 2;
                //IMD[i]=IMD[NbL-1];
                //PbMD[i]=PbMD[NbL-1];
                //PbMDdata[i]=PbMDdata[NbL-1];
                //SbMD[i]=SbMD[NbL-1];
                //--i;
                --NbL;
            }
        }
    }
    while(NbL>0);
    delete[] PbMD;
    delete[] SbMD;
    delete[] PbMDdata;
    delete[] IMD;
    PbMD = nullptr;
    SbMD = nullptr;
    PbMDdata = nullptr;
    IMD = nullptr;
}

void CalcBorder::preciseBorderCreate()
{
    boost::qvm::vec<double,3> r;
    boost::qvm::vec<int_fast32_t,3> ni,nj;
    uint_fast32_t j,l=0;//c=0;
    Pb = new boost::qvm::vec<double,3>[2*NbL];
    Pbdata = new StabilityPointType[2*NbL];
    Pbc = new boost::qvm::vec<double,3>[NbL];
    Pbcdata = new StabilityPointType[NbL];
    memset(Pbcdata,0,NbL*sizeof(StabilityPointType));
    Sbc = new boost::qvm::mat<double,3,3>[NbL];
    std::cerr<<"E1\n";
    for(int_fast32_t i=0; i<N; ++i)
    {
        if(Pdata[i].PointType==1)
            for(uint_fast32_t k=0; k<NM; ++k)
            {
                nj = Pi[i]+M[k];
                //std::cerr<<"Border1 "<<i<<" "<<nj.a[0]<<" "<<nj.a[1]<<" "<<n.a[0]<<" "<<n.a[1]<<"\n";
                if(nj.a[0]<0 || nj.a[0]>n.a[0]-1)
                    continue;
                if(nj.a[1]<0 || nj.a[1]>n.a[1]-1)
                    continue;
                if(nj.a[2]<0 || nj.a[2]>n.a[2]-1)
                    continue;
                j = I[nj.a[0] + nj.a[1]*n.a[0] + nj.a[2]*n.a[0]*n.a[1]];
                //std::cerr<<"Border1 "<<i<<" "<<j<<" "<<nj.a[0] + nj.a[1]*n.a[0] + nj.a[2]*n.a[0]*n.a[1]<<" "<<i<<" "<<Pi[i].a[0]<<" "<<Pi[i].a[1]<<" "<<Pi[i].a[2]<<"\n"
                //<<" "<<i<<" "<<nj.a[0]<<" "<<nj.a[1]<<" "<<nj.a[2]<<" "<<M[k].a[0]<<" "<<M[k].a[1]<<" "<<M[k].a[2]<<" "<<I[0]<<"\n";
                if(Pdata[j].PointType == -1)
                {
                    Pb[2*l] = P[i];
                    Pb[2*l+1] = P[j];
                    Pbc[l] = 0.5*(P[i]+P[j]);

                    Pbdata[2*l] = Pdata[i];
                    Pbdata[2*l+1] = Pdata[j];


                    //if(i==1271 || j==1271)
                    //std::cerr<<"Border2 "<<i<<" "<<j<<" "<<l<<" "<<Pbc[l].a[0]<<" "<<Pbc[l].a[1]<<"\n";
                    //std::cin.get();
                    ++l;
                }
            }
//SkipPoint:;
    }
//std::cerr<<"C "<<c<<"\n";
}

void CalcBorder::preciseBorderCalculate()
{
    NbLAll = NbL;
    boost::qvm::vec<double,3> tempP;
    StabilityPointType tempPdata;
    for(;NbL>0;)
    {
        //std::cerr<<"Q-1 "<<NbL<<"\n";
        SA->checkStability3D(Pbc,Pbcdata,Sbc,NbL);
        //std::cerr<<"Q0 "<<NbL<<"\n";
        for(int_fast32_t i=0; i<NbL; ++i)
        {
            if(boost::qvm::mag_sqr(Pb[2*i]-Pb[2*i+1])<BPBP)
            {
                if(Pbcdata[i].Stability==1){Pbc[i]=Pb[2*i];Pbcdata[i]=Pbdata[2*i];}
                tempP = Pbc[i];
                Pbc[i] = Pbc[NbL-1];
                Pbc[NbL-1] = tempP;
                tempPdata = Pbcdata[i];
                Pbcdata[i] = Pbcdata[NbL-1];
                Pbcdata[NbL-1] = tempPdata;

                tempP = Pb[2*i];
                Pb[2*i] = Pb[2*(NbL-1)];
                Pb[2*(NbL-1)] = tempP;
                tempPdata = Pbdata[2*i];
                Pbdata[2*i] = Pbdata[2*(NbL-1)];
                Pbdata[2*(NbL-1)] = tempPdata;

                tempP = Pb[2*i+1];
                Pb[2*i+1] = Pb[2*(NbL-1)+1];
                Pb[2*(NbL-1)+1] = tempP;
                tempPdata = Pbdata[2*i+1];
                Pbdata[2*i+1] = Pbdata[2*(NbL-1)+1];
                Pbdata[2*(NbL-1)+1] = tempPdata;
                //std::cerr<<"BorderP "<<i<<" "<<NbL<<" "<<Pbc[NbL-1].a[0]<<" "<<Pbc[NbL-1].a[1]<<"\n";
                --NbL;
                --i;
            } else if(Pbcdata[i].Stability==1)
            {
                Pb[2*i+1] = Pbc[i];
                Pbc[i] = 0.5*(Pb[2*i]+Pb[2*i+1]);
                Pbdata[2*i+1] = Pbcdata[i];
            } else if(Pbcdata[i].Stability==2)
            {
                Pb[2*i] = Pbc[i];
                Pbc[i] = 0.5*(Pb[2*i]+Pb[2*i+1]);
                Pbdata[2*i] = Pbcdata[i];
            }
        }
    }
    //std::cerr<<"E1\n";
    NbL = NbLAll;
    /**NbcL = 0;
    PbcL = new boost::qvm::vec<int_fast32_t,2>[NbL*NM];
    //std::cerr<<"E2\n";
    for(int_fast32_t i=0; i<NbL; ++i)
    {
        for(int_fast32_t j=0; j<NbL; ++j)
        {
            if(i==j)continue;
            if(boost::qvm::mag_sqr(Pbc[i]-Pbc[j])<AA)
            {
                PbcL[NbcL].a[0] = i;
                PbcL[NbcL].a[1] = j;

                std::cerr<<"PbcL "<<NbL*NM<<" "<<NbcL<<"\n";
                ++NbcL;
            }
        }
    }

    std::cerr<<"Q00 "<<NbL<<" "<<NbL*6<<" "<<NbcL<<"\n";
    int_fast32_t minj, tempNbcL=0, i, j;
    double rrmin, rr;
    boost::qvm::vec<double,3> r1, r2;
    boost::qvm::vec<int_fast32_t,2> *tempPbcL = new boost::qvm::vec<int_fast32_t,2>[NbL*2];
    //std::cerr<<"Q0 "<<NbL<<"\n";
    //std::cin.get();
    for(i=0; i<NbcL; ++i)
    {
        //std::cerr<<"Q1\n";
        rrmin = 1e99;
        minj = 0;
        for(j=0; PbcL[i+j].a[0]==PbcL[i].a[0]; ++j)
        {
            rr = boost::qvm::mag_sqr(Pbc[PbcL[i+j].a[0]]-Pbc[PbcL[i+j].a[1]]);
            if(rr<rrmin)
            {
                rrmin = rr;
                minj = j;
            }
        }
        //std::cerr<<"Q2\n";
        tempPbcL[tempNbcL].a[0] = PbcL[i+minj].a[0];
        tempPbcL[tempNbcL].a[1] = PbcL[i+minj].a[1];
        r1 = Pbc[tempPbcL[tempNbcL].a[0]]-Pbc[tempPbcL[tempNbcL].a[1]];
        ++tempNbcL;
        rrmin = 1e99;
        minj = 0;
        //std::cerr<<"Q3\n";

        for(j=0; PbcL[i+j].a[0]==PbcL[i].a[0]; ++j)
        {
            r2 = Pbc[PbcL[i+j].a[0]]-Pbc[PbcL[i+j].a[1]];

            rr = boost::qvm::mag_sqr(r2);
            if( (rr<rrmin) && (boost::qvm::dot( r1, r2 )<0) )
            {
                rrmin = rr;
                minj = j;
            }
        }
        //std::cerr<<"Q4\n";
        tempPbcL[tempNbcL].a[0] = PbcL[i+minj].a[0];
        tempPbcL[tempNbcL].a[1] = PbcL[i+minj].a[1];
        ++tempNbcL;
        i+=j;
        //std::cerr<<"Q "<<i<<" "<<tempNbcL<<" "<<NbL*2<<" "<<tempPbcL[tempNbcL-2].a[0]<<" "<<tempPbcL[tempNbcL-2].a[1]<<" "<<tempPbcL[tempNbcL-1].a[0]<<" "<<tempPbcL[tempNbcL-1].a[1]<<"\n";
        //std::cin.get();
    }
    //std::cerr<<"QQ\n";
    delete []PbcL;
    PbcL = tempPbcL;
    NbcL = tempNbcL;/**/

}
