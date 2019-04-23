#include "StabilityAnalytically.h"
#include "math_constants.h"
#include "phys_constants.h"
#include <boost/qvm/vec_mat_operations2.hpp>
#include <boost/qvm/mat_operations2.hpp>
#include <time.h>

StabilityAnalytically::StabilityAnalytically():generator(time(0)), distribution(3), distr(generator, distribution)
{
    //ctor
    NearZero = 1e-15;
    NearSphere = 1e-5;
    e = nullptr;
    ein = nullptr;
    ea = nullptr;
    R=nullptr;
    A = nullptr;
    Ain = nullptr;
    _1d_A = nullptr;
    P1 = nullptr;
    P2 = nullptr;
    CSD = nullptr;
    CSN = nullptr;

    n.a[0] = 40;
    aStep.a[0]=1.0*MC_s2;
    aStep.a[1]=1.0*MC_s2;
    aStep.a[2]=1.0*MC_s2;
    MaxCSDN = 134;
    EpsDeform_const = 1e-6;

    P_D = 1.0;
    P_alfa = 6.0;
    P_a = 1.0;
    P_aCut = 0;
    P_P1 = 2.0*P_D*P_alfa/P_a;
    P_P2 = 2.0*P_D*P_alfa*P_alfa/(P_a*P_a);
//std::cerr<<"R-2 "<<n.a[0]<<"\n";
}

StabilityAnalytically::~StabilityAnalytically()
{
    //dtor
    delete []e;
    e = NULL;
    delete []ein;
    ein = NULL;
    delete []ea;
    ea = NULL;
    delete []R;
    R = NULL;
    delete []A;
    A = NULL;
    delete []Ain;
    Ain = NULL;
    delete []_1d_A;
    _1d_A = NULL;
    delete []P1;
    P1 = NULL;
    delete []P2;
    P2 = NULL;
    delete []CSD;
    CSD = NULL;
    delete []CSN;
    CSN = NULL;
}

void StabilityAnalytically::startTask()
{
    //std::cerr<<"\n\n\n";
    //std::cerr<<"E0\n";
    //createTriangularLattice();
    createFCCLattice();
    //std::cerr<<"E1\n";
    //std::cin.get();
    findCoordinationalSpheres();
    P_aCut = 0.5*(CSD[3]+CSD[4]);
    std::cerr<<P_aCut<<"\n";
    //std::cerr<<"E2\n";
    //std::cin.get();
    decomposeVectors();
    //std::cerr<<"E3\n";
    //std::cin.get();
    calculateActualStiffness();
    reachSpecifiedStress();
    //std::cerr<<"E4\n";
    //std::cin.get();
}

void StabilityAnalytically::checkStability2D(boost::qvm::vec<double,3> *Pvar, StabilityPointType *Pdvar, boost::qvm::mat<double,3,3> *Svar, uint_fast32_t &Nvar)
{

    boost::qvm::mat<double,3,3> defGrad=MC_0T33, dG;
    double _1d_N = 1.0/double(Nvar);
    for(uint_fast32_t i=0; i<Nvar; ++i)
    {
        boost::qvm::set_identity(defGrad);
        defGrad.a[0][0] += Pvar[i].a[0];
        defGrad.a[0][1] += (eXYMax<Pvar[i].a[2])?Pvar[i].a[2]:(Pvar[i].a[2]-eXYMax*int_fast32_t(Pvar[i].a[2]/eXYMax));
        defGrad.a[1][1] += Pvar[i].a[1];


        deformVectors(defGrad);
        checkStability(Pdvar[i]);

        ///calculateActualStiffness(defGrad);
        ///Svar[i] = StiffnessVoigt;
        /*std::cerr<<Pvar[i].a[0]<<" "<<Pvar[i].a[1]<<"\n";
        for(uint_fast16_t i=0; i<3; ++i)
        {
            for(uint_fast16_t j=0; j<3; ++j)
            {
                if(fabs(StiffnessVoigt.a[i][j])<1e-10)
                    StiffnessVoigt.a[i][j]=0;
                std::cerr.width(7);
                std::cerr.precision(3);
                std::cerr<<StiffnessVoigt.a[i][j]*stress_const*1e-9<<" ";
            }
            std::cerr<<"\n";
        }
        std::cin.get();/**/
        //if(i%100)std::cout<<100.0*i*_1d_N<<"%\n";
    }
}

void StabilityAnalytically::checkStability3D(boost::qvm::vec<double,3> *Pvar, StabilityPointType *Pdvar, boost::qvm::mat<double,3,3> *Svar, uint_fast32_t &Nvar)
{

    boost::qvm::mat<double,3,3> defGrad=MC_0T33, dG;
    double _1d_N = 1.0/double(Nvar);
    for(uint_fast32_t i=0; i<Nvar; ++i)
    {
        boost::qvm::set_identity(defGrad);
        defGrad.a[0][0]+=Pvar[i].a[0];
        defGrad.a[1][1]+=Pvar[i].a[1];
        defGrad.a[2][2]+=Pvar[i].a[2];

        deformVectors(defGrad);
        checkStability(Pdvar[i]);
        if(i%10000==0)
        std::cerr<<"Done "<<100.0*i/double(Nvar)<<"%\n";

        ///calculateActualStiffness(defGrad);
        ///Svar[i] = StiffnessVoigt;
        /*std::cerr<<Pvar[i].a[0]<<" "<<Pvar[i].a[1]<<"\n";
        for(uint_fast16_t i=0; i<3; ++i)
        {
            for(uint_fast16_t j=0; j<3; ++j)
            {
                if(fabs(StiffnessVoigt.a[i][j])<1e-10)
                    StiffnessVoigt.a[i][j]=0;
                std::cerr.width(7);
                std::cerr.precision(3);
                std::cerr<<StiffnessVoigt.a[i][j]*stress_const*1e-9<<" ";
            }
            std::cerr<<"\n";
        }
        std::cin.get();/**/
        //if(i%100)std::cout<<100.0*i*_1d_N<<"%\n";
    }
}

void StabilityAnalytically::deformVectors(boost::qvm::mat<double,3,3> &dG)
{
    V *= boost::qvm::determinant(dG);
    _1d_V = 1.0/V;
    //std::cerr<<sizeof(boost::qvm::vec<double,2>)<<" "<<2*sizeof(double)<<"\n";
    memset(ea, 0, sizeof(boost::qvm::vec<double,3>)*N);
    memset(e, 0, sizeof(boost::qvm::vec<double,3>)*N);

    double Norm, _1d_Norm;
    //int periodicTg;

    //LatticeStep.x = Afina[0+0]*LatticeStep0.x;

    for (uint_fast32_t i=0; i<RSN; ++i)
    {
        //periodicTg = int(Afina[0+1]*e_in[2*k+1]/LatticeStep.x);
        ea[i] = dG*ein[i];
        Norm = boost::qvm::mag(ea[i]);
        _1d_Norm = 1.0/Norm;
        A[i] = Ain[i]*Norm;
        _1d_A[i] = 1.0/A[i];
        e[i] = ea[i]*_1d_Norm;;
    }
}

void StabilityAnalytically::checkStability(StabilityPointType &Pd)
{
    calculateForces();
    //std::cin.get();
    boost::qvm::mat<double,9,9> Q;
    boost::qvm::set_zero(Q);
    boost::qvm::vec<double,3> W;
    boost::qvm::mat<double,3,3> AA, WW, D;
    //boost::random::uniform_on_sphere<double>::result_of tmp;
    double Inv[3], Delta;
    std::vector<double> res;
    for(uint_fast32_t i=1; i<RSN; ++i)
    {
        AA = tens(e[i],e[i]);
        //std::cerr<<"Q1 "<<i<<" "<<e[i].a[0]<<" "<<e[i].a[1]<<"\n";
        Q += A[i]*A[i]*(P1[i]*_1d_A[i]*tens(MC_1T33,AA)+(P2[i]-P1[i]*_1d_A[i])*tens(AA,AA));
        //std::cerr<<"Q "<<i<<" "<<A[i]<<" "<<P1[i]<<" "<<_1d_A[i]<<" "<<P2[i]<<"\n";
    }

    for(uint_fast8_t i=1; i<9; ++i)
        for(uint_fast8_t j=1; j<9; ++j)
        {
            if (fabs(Q.a[i][j]) < near_zero)
            {
                Q.a[i][j]=0;
            }
        }


    Pd.StabilitySteps = 1;
    Pd.Stability = 2;
    for (int m=0; m<10000; m++)
    {
        res= distr();
        W.a[0] = res[0];
        W.a[1] = res[1];
        W.a[2] = res[2];

        WW = tens(W,W);
        D = dotdot(Q,WW);
        if(fabs((boost::qvm::mag_sqr(W))-1.0)>1e-14)
        {
            std::cerr<<W.a[0]<<" "<<W.a[1]<<" "<<W.a[2]<<" "<<(boost::qvm::mag_sqr(W))-1.0<<"\n";
            std::cin.get();
        }


        //Inv[0] = D.a[0][0]+D.a[1][1]+D.a[2][2];
        //Inv[1] = D.a[0][0]*D.a[1][1]+D.a[0][0]*D.a[2][2]+D.a[1][1]*D.a[2][2]-D.a[0][1]*D.a[1][0]-D.a[0][2]*D.a[2][0]-D.a[1][2]*D.a[2][1];

        Inv[0] = D.a[0][0];
        Inv[1] = D.a[0][0]*D.a[1][1]-D.a[0][1]*D.a[1][0];
        //D[3*0+0]*D[3*1+1]+D[3*0+0]*D[3*2+2]+D[3*1+1]*D[3*2+2]-D[3*0+1]*D[3*1+0]-D[3*0+2]*D[3*2+0]-D[3*1+2]*D[3*2+1];
        Inv[2] = boost::qvm::determinant(D);

        if(Inv[0]<near_zero || Inv[1]<near_zero || Inv[2]<near_zero)
        {
            Pd.StabilitySteps = 0;
            Pd.Stability = 1;
            break;
        }
    }


    /*for (uint_fast16_t i=0; i<4; ++i)
    {
        for (uint_fast16_t j=0; j<4; ++j)
        {
            if(fabs(Q.a[i][j])<NearZero)
            {
                Q.a[i][j] = 0;
            }
            std::cerr<<Q.a[i][j]<<" ";
        }
        std::cerr<<"\n";
    }
    std::cin.get();/**/

    //Pd.StabilitySteps = 1;
    //Pd.Stability = 1;

    /*double A=Q.a[0][0]*Q.a[3][0], B=0.5*(Q.a[0][0]*Q.a[3][3]+Q.a[0][3]*Q.a[3][0]-4.0*Q.a[1][1]*Q.a[1][1]), C=Q.a[0][3]*Q.a[3][3];

    if((A<=0) || (C<=0) || (B<-sqrt(A*C)))
    {
            //std::cerr<<"Q1\n";
            Pd.StabilitySteps = 0;
            Pd.Stability = 0;
    }/**/



    /*for (uint_fast32_t i=0; i<10000; ++i)
    {

        W.a[0] = rand()*_1d_RAND_MAX_double*2.0-1.0;
        W.a[1] = rand()*_1d_RAND_MAX_double*2.0-1.0;

        WW = tens(W,W);
        D = dotdot(Q,WW);

        Inv[0] = D.a[0][0];
        Inv[1] = boost::qvm::determinant(D);
        //Inv[2] = -D[2*2+0]*D[2*1+1]*D[2*0+2];


        if(fabs(Inv[0]) < NearZero)
            Inv[0] = 0;

        if(fabs(Inv[1]) < NearZero)
            Inv[1] = 0;

        //std::cerr<<"Stable "<<Inv[0]<<" "<<Inv[1]<<"\n";
        if((Inv[0]<=0) || (Inv[1]<=0))
        {
            Pd.StabilitySteps = 0;
            Pd.Stability = 0;
            break;
        }
    }/**/
}

void StabilityAnalytically::calculateForces()
{

    double exp_b_ra;
    boost::qvm::set_zero(Stress);
    E = 0;
    //std::cerr<<"F "<<RSN<<"\n";
    P1[0]=0;
    P2[0]=0;
    for (uint_fast32_t i=1; i<RSN; ++i)
    {
        //std::cerr<<"F "<<i<<" "<<A[i]<<"\n";
        //std::cin.get();
        if (A[i]<P_aCut)
        {
            exp_b_ra = exp( P_alfa*(P_a-A[i]) );
            P1[i] = -P_P1*( exp_b_ra - 1.0 )*exp_b_ra; //силы взаимодействий
            P2[i] = P_P2*( 2.0*exp_b_ra - 1.0 )*exp_b_ra; //жесткости связей ЦИКЛ ОТ 1 ДО Nvect
            E += P_D*( exp_b_ra - 2.0 )*exp_b_ra;

            //std::cerr<<"F "<<i<<" "<<A[i]<<" "<<P1[i]<<" "<<P2[i]<<"\n";

        }
        else
        {
            //break;
            P1[i] = 0;
            P2[i] = 0;
        }

        Stress += (_1d_V*A[i]*P1[i])*tens(e[i],e[i]);
        //std::cerr<<"Stress "<<_1d_V*A[i]*P1[i]<<" "<<Stress.a[0][0]<<" "<<Stress.a[0][1]<<" "<<Stress.a[1][0]<<" "<<Stress.a[1][1]<<"\n";
    }
    //std::cin.get();
}

void StabilityAnalytically::findCoordinationalSpheres()
{
    CSD = new double[N];
    CSN = new uint_fast32_t[N];
    boost::qvm::vec<double,3> dr;
    double dr_mm;
    bool find;
    //std::cerr<<"Q3"<<CoordSphereParticlesNumber<<"\n";
    CSDN = 0;
    memset(CSD, 0, N*sizeof(double));
    memset(CSN, 0, N*sizeof(uint_fast32_t));
    for(uint_fast32_t i=0; i<N; ++i)
    {
        dr_mm = boost::qvm::mag_sqr(R[i]-R[iCenter]);

        //std::cerr<<"R "<<boost::qvm::X(R[i]-R[iCenter])<<" "<<boost::qvm::Y(R[i]-R[iCenter])<<" "<<boost::qvm::Z(R[i]-R[iCenter])<<" "<<dr_mm<<"\n";
        //std::cin.get();
        for(uint_fast32_t j=0; j<CSDN; ++j)
        {
            if(fabs(dr_mm-CSD[j])<NearSphere)
            {

                ++CSN[j];
                //std::cerr<<"W1 "<<i<<" "<<j<<" "<<CSN[j]<<"\n";
                goto FindSphere;
            }
        }
        CSD[CSDN] = dr_mm;
        //std::cerr<<"W "<<i<<" "<<CSD[CSDN]<<" "<<CSDN<<"\n";
        ++CSN[CSDN];
        ++CSDN;
FindSphere:
        ;
    }
    //std::cerr<<"CoordSphereDistanceNumber "<<CSDN<<"\n";

    uint_fast32_t ntmp;
    do
    {
        find=false;
        for(uint_fast32_t j=1; j<CSDN; ++j)
        {
            if(CSD[j-1]>CSD[j])
            {
                dr_mm = CSD[j];
                CSD[j] = CSD[j-1];
                CSD[j-1] = dr_mm;

                ntmp = CSN[j];
                CSN[j] = CSN[j-1];
                CSN[j-1] = ntmp;
                find=true;
            }
        }
    }
    while(find);

    double *tempCDS = new double[MaxCSDN];
    memcpy(tempCDS, CSD, MaxCSDN*sizeof(double));
    delete[] CSD;
    CSD = tempCDS;
    uint_fast32_t *tempCSDN = new uint_fast32_t[MaxCSDN];
    memcpy(tempCSDN, CSN, MaxCSDN*sizeof(uint_fast32_t));
    delete[] CSN;
    CSN = tempCSDN;

//std::cerr<<"Q1\n";

    uint_fast32_t *SVp = new uint_fast32_t[MaxCSDN+1];
    uint_fast32_t *SVc = new uint_fast32_t[MaxCSDN];
    memset(SVc, 0, sizeof(uint_fast32_t)*MaxCSDN);
//std::cerr<<"Q2\n";
    RSN=1;
    SVp[0] = 0;
    for(uint_fast32_t i=1; i<MaxCSDN; ++i)
    {
        SVp[i] = RSN;
        RSN+=CSN[i]/2;//убираем симметричные вектора e[-k] = -e[k]
        //if(CSN[i]%2)std::cerr<<"ERRRRRRR!!!!\n\n\n";
        //std::cerr<<"Q "<<i<<" "<<RSN<<" "<<CSN[i]<<" "<<CSD[i]<<" "<<SVp[i]<<"\n";
        //std::cerr<<"Q3 "<<i<<" "<<CoordSphereParticlesNumber[i]<<" "<<SphereVectorsPosition[i]<<"\n";
    }
    //std::cin.get();
//std::cerr<<"Q3 "<<VectorNumber<<"\n";
    SVp[MaxCSDN] = RSN;
    RS = new boost::qvm::vec<double,3>[RSN];
    memset(RS,0,RSN*sizeof(boost::qvm::vec<double,3>));
    //std::cerr<<"Sphere "<<RSN<<"\n";
    //std::cerr<<boost::qvm::X(RS[22])<<" "<<boost::qvm::Y(RS[22])<<"\n";
    //std::cin.get();
    RSN = 0;
    for(uint_fast32_t i=0; i<N; ++i)
    {
        //std::cerr<<"R1\n";
        dr = R[i]-R[iCenter];
        dr_mm = boost::qvm::mag_sqr(dr);
        //std::cerr<<"SVp "<<i<<" "<<dr_mm<<" "<<boost::qvm::X(R[i]-R[iCenter])<<" "<<boost::qvm::Y(R[i]-R[iCenter])<<"\n";

        for(uint_fast32_t j=0; j<MaxCSDN; ++j)
        {
            //std::cerr<<"R1.5 "<<j<<" "<<CSD[j]<<"\n";
            //std::cin.get();
            if(fabs(dr_mm-CSD[j])<NearSphere)
            {
                //std::cerr<<"SVp "<<j<<"\n";
                //std::cin.get();
                //std::cerr<<"SVp "<<i<<" "<<j<<" "<<RSN<<" "<<dr_mm<<" "<<CSD[j];//<<" "<<boost::qvm::X(R[i]-R[iCenter])<<" "<<boost::qvm::Y(R[i]-R[iCenter])<<"\n";
                //if(j==9)
                //std::cerr<<"SVp "<<j<<" "<<SVp[j]<<" "<<SVc[j]<<" "<<SVp[j]+SVc[j]<<" "<<boost::qvm::X(RS[SVp[j]])<<" "<<boost::qvm::Y(RS[SVp[j]])<<"\n";
                for(uint_fast32_t k=SVp[j]; k<SVp[j]+SVc[j]; ++k)
                {
                    //std::cerr<<"R2\n";
                    //std::cerr<<"SVP "<<k<<" "<<boost::qvm::X(RS[k])<<" "<<boost::qvm::Y(RS[k])<<"\n";
                    //std::cerr<<"SVP "<<k<<" "<<boost::qvm::X(RS[k]+dr)<<" "<<boost::qvm::Y(RS[k]+dr)<<" "<<boost::qvm::mag_sqr(RS[k]+dr)<<"\n";
                    if(boost::qvm::mag_sqr(RS[k]+dr)<NearSphere)
                        goto SkipVector;
                }


                RS[SVp[j]+SVc[j]]=dr;
                //std::cerr<<"SVV "<<SVp[j]+SVc[j]<<" "<<boost::qvm::X(RS[SVp[j]+SVc[j]])<<" "<<boost::qvm::Y(RS[SVp[j]+SVc[j]])<<"\n";
                ++SVc[j];
                ++RSN;
                //std::cin.get();
            }/**/
SkipVector:;
            //std::cerr<<"R3\n";

        }
    }
    //std::cin.get();
    //std::cerr<<"Sphere "<<RSN<<"\n";

    /*for(uint_fast32_t i=0; i<MaxCoordSphereNumber; ++i)
    {
        std::cerr<<"Q6 "<<i<<" "<<SphereVectorsPosition[i]<<" "<<SphereVectorsPosition[i+1]<<" "<<VectorNumber<<"\n";
        for(uint_fast32_t j=SphereVectorsPosition[i]; j<SphereVectorsPosition[i+1]; ++j)
        {
            std::cerr<<i<<" "<<j<<" "<<Vectors[j].x*Vectors[j].x+Vectors[j].y*Vectors[j].y<<" "<<Vectors[j].x<<" "<<Vectors[j].y<<"\n";
        }
        std::cin.get();
    }/**/
    delete[] SVp;
    SVp = nullptr;
    delete[] SVc;
    SVc = nullptr;
    //std::cerr<<"Sphere "<<CSD<<"\n";
    //std::cin.get();
    for(uint_fast32_t i=0; i<MaxCSDN; ++i)
    {
        CSD[i] = sqrt(CSD[i]);
        std::cerr<<"Sphere "<<i<<" "<<CSD[i]<<" "<<CSD[i]*CSD[i]<<"\n";
    }
    //std::cin.get();
}

void StabilityAnalytically::decomposeVectors()
{
    for (uint_fast32_t i=0; i<RSN; ++i)
    {
        Ain[i] = boost::qvm::mag(RS[i]);
        //std::cerr<<"DV "<<RS[i].a[0]<<" "<<RS[i].a[1]<<" "<<Ain[i]<<"\n";
        _1d_A[i] = 1.0/Ain[i];
        ein[i]   = RS[i]*_1d_A[i];
    }
    memcpy(A,Ain,RSN*sizeof(double));
    memcpy(e,ein,RSN*sizeof(boost::qvm::vec<double,3>));
}

void StabilityAnalytically::createTriangularLattice()
{
//std::cerr<<"R-2 "<<n.a[0]<<"\n";
    boost::qvm::Y(n) = 2.0*uint_fast32_t(boost::qvm::X(n)*MC_1ds3);
    N = boost::qvm::X(n)*boost::qvm::Y(n);
    std::cerr<<"NumberLattice "<<N<<"\n";
///std::cerr<<"R-1\n";
    e = new boost::qvm::vec<double,3>[N];
//std::cerr<<"R0\n";
    ein = new boost::qvm::vec<double,3>[N];
//std::cerr<<"R1\n";
    ea = new boost::qvm::vec<double,3>[N];
//std::cerr<<"R2\n";
    R = new boost::qvm::vec<double,3>[N];
//std::cerr<<"R3\n";
    A = new double[N];
//std::cerr<<"R4\n";
    Ain = new double[N];
//std::cerr<<"R5\n";
    _1d_A = new double[N];
    P1 = new double[N];
    P2 = new double[N];

    //I = new uint_fast32_t[NL];
//std::cerr<<"R10\n";
    boost::qvm::vec<double,3> tmp= {0,0,0};
    uint_fast32_t i, j, iN = 0;
    for(j=0; j<boost::qvm::Y(n); ++j)
    {
        boost::qvm::X(tmp) = 0.5*boost::qvm::X(aStep)*(j%2);
        for(i=0; i<boost::qvm::X(n); ++i, ++iN)
        {
            R[iN].a[0] = i*aStep.a[0]+tmp.a[0];
            R[iN].a[1] = j*aStep.a[1]+tmp.a[1];
            R[iN].a[2] = 0;
            //std::cerr<<"PR "<<iN<<" "<<R[iN].a[0]<<" "<<R[iN].a[1]<<"\n";
        }
    }
    iCenter = boost::qvm::X(n)*(boost::qvm::Y(n)/2)+boost::qvm::X(n)/2;
    V=(R[N-1].a[0]*R[N-1].a[1])/double(N);
    _1d_V = 1.0/V;
    eXYMax = 2.0*aStep.a[1]/aStep.a[0];
    std::cerr<<"Triangular "<<N<<" "<<V<<" "<<iCenter<<"\n";
}

void StabilityAnalytically::createFCCLattice()
{
    n.a[1] = n.a[0];
    n.a[2] = n.a[0];

    n*=2;
    n+=MC_i1XYZV3;
    N = n.a[0]*n.a[1]*n.a[2]/2+1;

    e = new boost::qvm::vec<double,3>[N];
//std::cerr<<"R0\n";
    ein = new boost::qvm::vec<double,3>[N];
//std::cerr<<"R1\n";
    ea = new boost::qvm::vec<double,3>[N];
//std::cerr<<"R2\n";
    R = new boost::qvm::vec<double,3>[N];
//std::cerr<<"R3\n";
    A = new double[N];
//std::cerr<<"R4\n";
    Ain = new double[N];
//std::cerr<<"R5\n";
    _1d_A = new double[N];
    P1 = new double[N];
    P2 = new double[N];
    //I = new uint_fast32_t[NL];
    boost::qvm::vec<double,3> rCenter;
    rCenter.a[0] = n.a[0]*0.25*aStep.a[0];
    rCenter.a[1] = n.a[1]*0.25*aStep.a[1];
    rCenter.a[2] = n.a[2]*0.25*aStep.a[2];
    double drrCenterMin = 1e99, drrCenter;
    uint_fast32_t i, j, k, iN = 0;
    for(k=0; k<n.a[2]; ++k)
        for(j=0; j<n.a[1]; ++j)
            for(i=0; i<n.a[0]; ++i)
            {
                if((i%2+j%2+k%2)%2)
                    continue;

                R[iN].a[0] = i*0.5*aStep.a[0];
                R[iN].a[1] = j*0.5*aStep.a[1];
                R[iN].a[2] = k*0.5*aStep.a[2];
                //std::cerr<<"PR "<<iN<<" "<<R[iN].a[0]<<" "<<R[iN].a[1]<<" "<<R[iN].a[2]<<"\n";

                drrCenter = boost::qvm::mag_sqr(rCenter-R[iN]);
                if(drrCenter<drrCenterMin)
                {
                    drrCenterMin = drrCenter;
                    iCenter = iN;

                }
                //std::cerr<<"C "<<iN<<" "<<i<<" "<<j<<" "<<k<<" "<<P[iN].a[0]<<" "<<P[iN].a[1]<<" "<<P[iN].a[2]<<" "<<Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]<<" "<<I[Pi[iN].a[0]+Pi[iN].a[1]*n.a[0]+Pi[iN].a[2]*n.a[0]*n.a[1]]<<"\n";
                //<<" "<<ParticleCoordinate[Particles.ParticlesNumber].y
                //<<" "<<ParticleCoordinate[Particles.ParticlesNumber].z<<"\n";
                //std::cin.get();
                ++iN;
            }
    V = aStep.a[0]*aStep.a[1]*aStep.a[2]*0.25;
    //V = R1.a[0]*R2.a[1]-R1.a[1]*R2.a[0];
    //V=boost::qvm::mag(boost::qvm::cross(R[n.a[0]*(n.a[1]-1)]-R[0],));
    _1d_V = 1.0/V;
        //std::cerr<<"C "<<N<<" "<<eStep.a[0]<<" "<<eStep.a[1]<<" "<<eStep.a[2]<<"\n";
    std::cerr<<"FCC "<<N<<" "<<V<<" "<<iCenter<<"\n";
    //std::cin.get();
}

void StabilityAnalytically::calculateActualStiffness(const boost::qvm::mat<double,3,3> &dG)
{
    double EpsDeform_const = 1e-6;
    double _1d_Eps2_const = 1.0/(EpsDeform_const*EpsDeform_const);
    boost::qvm::mat<double,3,3> defGS[18], StressS[18], dStress, dGtemp = dG-MC_1T33;
    for(uint_fast32_t i=0; i<3; ++i)
    {
        for(uint_fast32_t j=0; j<3; ++j)
        {
            boost::qvm::set_identity(defGS[3*i+j]);
            defGS[3*i+j].a[i][j]   += EpsDeform_const;
            boost::qvm::set_identity(defGS[3*i+j+9]);
            defGS[3*i+j+9].a[i][j] -= EpsDeform_const;
        }
    }

    for(uint_fast32_t i=0; i<18; ++i)
    {
        defGS[i] += dGtemp;
        deformVectors(defGS[i]);
        calculateForces();
        StressS[i] = Stress;
        //std::cerr<<"Stress "<<" "<<Stress.a[0][0]<<" "<<Stress.a[0][1]<<" "<<Stress.a[1][0]<<" "<<Stress.a[1][1]<<"\n";

    }

    boost::qvm::set_zero(Stiffness);
    for(uint_fast32_t i=0; i<9; ++i)
    {
        dStress = 0.5*(StressS[i]-StressS[i+9]);
        for(uint_fast8_t j=0;j<3;++j)
            for(uint_fast8_t k=0;k<3;++k)
                if(fabs(dStress.a[j][k])<1e-8)dStress.a[j][k]=0;
        defGS[i]-=MC_1T33+dGtemp;
        //std::cerr<<dStress.a[0][0]<<" "<<dStress.a[0][1]<<" "<<dStress.a[0][2]<<"\n"
        //<<dStress.a[1][0]<<" "<<dStress.a[1][1]<<" "<<dStress.a[1][2]<<"\n"
        //<<dStress.a[2][0]<<" "<<dStress.a[2][1]<<" "<<dStress.a[2][2]<<"\n"<<"\n";
        //std::cerr<<"Stress "<<" "<<defGS[i].a[0][0]<<" "<<defGS[i].a[0][1]<<" "<<defGS[i].a[1][0]<<" "<<defGS[i].a[1][1]<<"\n";
        Stiffness += _1d_Eps2_const*tens(dStress,defGS[i]);
    }

    toVoigtNotationC3D(Stiffness,StiffnessVoigt);
    /**for(uint_fast16_t i=0; i<4; ++i)
    {
        for(uint_fast16_t j=0; j<4; ++j)
        {
            if(fabs(Stiffness.a[i][j])<1e-10)
                StiffnessVoigt.a[i][j]=0;
            std::cerr.width(7);
            std::cerr.precision(3);
            std::cerr<<Stiffness.a[i][j]*stress_const*1e-9<<" ";
        }
        std::cerr<<"\n";
    }/**/
    /*for(uint_fast16_t i=0; i<6; ++i)
    {
        for(uint_fast16_t j=0; j<6; ++j)
        {
            if(fabs(StiffnessVoigt.a[i][j])<1e-10)
                StiffnessVoigt.a[i][j]=0;
            std::cerr.width(7);
            std::cerr.precision(3);
            std::cerr<<StiffnessVoigt.a[i][j]*stress_const*1e-9<<" ";
        }
        std::cerr<<"\n";
    }

    std::cerr<<"\n";/**/
}

void StabilityAnalytically::reachSpecifiedStress(const boost::qvm::mat<double,3,3> &Svar)
{
    boost::qvm::mat<double,3,3> dS=Stress-Svar, eps_tmp, eps_tmpsum;
    boost::qvm::set_zero(eps_tmpsum);
    uint_fast32_t steps=0;
    for(;;)
    {
        calculateForces();
        dS=Stress-Svar;
        if(fabs(dS.a[0][0])*fabs(dS.a[1][1])<1e-25)
            break;

        calculateActualStiffness();
        DuctilityVoigt = boost::qvm::inverse(StiffnessVoigt);
        fromVoigtNotationS3D(DuctilityVoigt,Ductility);
        eps_tmp = -dotdot(Ductility,dS);

        eps_tmpsum += eps_tmp*(MC_1T33+eps_tmpsum);

//std::cerr<<"Stress "<<" "<<dS.a[0][0]<<" "<<dS.a[0][1]<<" "<<dS.a[1][0]<<" "<<dS.a[1][1]<<"\n";

//std::cerr<<"E "<<eps_tmp.a[0][0]<<" "<<eps_tmp.a[0][1]<<" "<<eps_tmp.a[1][0]<<" "<<eps_tmp.a[1][1]
//<<" "<<eps_tmpsum.a[0][0]<<" "<<eps_tmpsum.a[1][1]<<std::endl<<std::endl;

        eps_tmp += MC_1T33;
        deformVectors(eps_tmp);
        memcpy(ein,e,RSN*sizeof(boost::qvm::vec<double,3>));
        memcpy(Ain,A,RSN*sizeof(double));
        ++steps;


        //std::cin.get();
    }
    std::cerr<<"E "<<eps_tmpsum.a[0][0]<<" "<<eps_tmpsum.a[1][1]<<" "<<eps_tmpsum.a[2][2]<<" "<<steps<<"\n";
    //std::cin.get();
}
