#include <iostream>
#include <boost/qvm/vec.hpp>
#include "CalcBorder.h"
#include "StabilityAnalytically.h"
#include "MD.h"
#include <mpi.h>
using namespace std;
int Nproc, iproc;
int main(int argc, char *argv[])
{
    char taskdir[256]=".", filename[256]="";
    //boost::qvm::vec<float,3> v={0,0,7};
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

    if(iproc==0)
    {
        CalcBorder CB;
        // CB.createPoints3D();
        CB.createPoints3DLineBCC();
        //std::cerr<<"Q-1\n";
        CB.createIMatrix();
        //std::cin.get();
        //std::cerr<<"Q0\n";
        ///CB.checkStability();
        CB.checkStability_MPI();
    }
    else
    {
        MD MDC;
        taskSRMPIi(&MDC);
    }

    MPI_Finalize();

    //std::cerr<<"Q1\n";
    ///CB.findBorder();
    //std::cerr<<"Q2\n";
    ///CB.preciseBorderCreate();
    //std::cerr<<"Q3\n";
    ///CB.preciseBorderCalculate();
    //std::cerr<<"Q4\n";
    ///boost::qvm::mat<double,3,3> D;
    ///ElasticModules EM;
    /**strcpy(filename, taskdir);
    strcat(filename, "/Result.txt");
    std::ofstream Result_file(filename);
    for(uint_fast32_t i=0;i<CB.N;++i)
    {
        ///toElastic2D(CB.S[i],EM);
        //if(CB.Pdata[i].PointType==1)
        ///if(EM.Ex<0)
        Result_file<<CB.P[i].a[0]<<" "<<CB.P[i].a[1]<<" "<<CB.P[i].a[2]<<" "<<CB.Pdata[i].StabilitySteps<<" "<<EM.Ex<<" "<<EM.Ey<<" "<<EM.nuxy<<" "<<EM.nuyx<<" "<<EM.G<<" "<<EM.K<<"\n";
    }
    Result_file.close();
    //std::cerr<<"Q4\n";
    //std::cin.get();
    strcpy(filename, taskdir);
    strcat(filename, "/ResultBc.txt");
    std::ofstream ResultBc_file(filename);
    for(uint_fast32_t i=0;i<CB.NbLAll;++i)
    {

        ///toElastic2D(CB.Sbc[i],EM);
        //if(EM.Ex<0)
        ResultBc_file<<CB.Pbc[i].a[0]<<" "<<CB.Pbc[i].a[1]<<" "<<CB.Pbc[i].a[2]<<" "<<CB.Pbcdata[i].StabilitySteps<<" "<<i<<"\n";
    }
    ResultBc_file.close();/**/
    /*strcpy(filename, taskdir);
    strcat(filename, "/ResultBcV.txt");
    std::ofstream ResultBcL_file(filename);
    for(uint_fast32_t i=0;i<CB.NbcL;i+=2)
    {
        ResultBcL_file<<CB.Pbc[CB.PbcL[i].a[0]].a[0]<<" "<<CB.Pbc[CB.PbcL[i].a[0]].a[1]<<" "<<CB.Pbc[CB.PbcL[i].a[1]].a[0]<<" "<<CB.Pbc[CB.PbcL[i].a[1]].a[1]<<"\n"
        <<CB.Pbc[CB.PbcL[i+1].a[0]].a[0]<<" "<<CB.Pbc[CB.PbcL[i+1].a[0]].a[1]<<" "<<CB.Pbc[CB.PbcL[i+1].a[1]].a[0]<<" "<<CB.Pbc[CB.PbcL[i+1].a[1]].a[1]<<"\n";
    }
    ResultBcL_file.close();/**/
    //cout << "Hello world! " <<sizeof(int*)<<" "<<v.a[2]<< endl;
    return 0;
}
