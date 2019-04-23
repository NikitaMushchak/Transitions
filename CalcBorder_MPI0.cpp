#include "CalcBorder.h"
#include <mpi.h>
#include <time.h>
extern int_fast32_t Nproc, iproc;

void CalcBorder::checkStability_MPI()
{
    //int t1=clock(),t2;
	char filename[256]="";
    checkStabilityMPIMD(P,Pdata,S,N);

    strcpy(filename, "./ResultnStable1.txt");
    std::ofstream ResultStable1_file(filename);

	strcpy(filename, "./Resultn1.txt");
    std::ofstream Result1_file(filename);
	std::cout<<" NNNNNNNNNNNNNN = "<<N<<std::endl;
    for(uint_fast32_t i=0; i<N; ++i){

		Result1_file<<P[i].a[0] <<" "<<P[i].a[1]<<" "<<P[i].a[2]<<
			" "<<Pdata[i].StabilitySteps<<" "<<Pdata[i].StabilityTime<<" "
			<<Pdata[i].PotEnergy<<" "<<Pdata[i].NumberBonds<<" "
			<<int(Pdata[i].Stability)<<" 111111"<<"\n";


		if(int(Pdata[i].Stability)==2)
            ResultStable1_file<<P[i].a[0]-1.<<" "<<P[i].a[1]-1.<<" "<<P[i].a[2]-1<<
			" "<<Pdata[i].StabilitySteps<<" "<<Pdata[i].StabilityTime<<"\n";
    }
    ResultStable1_file.close();
	strcpy(filename, "./ResultC1.dat");
    std::ofstream ResultC1_file(filename, std::ios::binary);
    ResultC1_file.write((char*)&N,sizeof(uint_fast32_t));
    ResultC1_file.write((char*)P,N*sizeof(boost::qvm::vec<double,3>));
    ResultC1_file.write((char*)Pdata,N*sizeof(StabilityPointType));
    ResultC1_file.close();/**/
    //t2=clock();
    //std::cerr<<"TIME "<<double(t2-t1)/(2*CLOCKS_PER_SEC*112)*5000<<"\n";
	/*
	strcpy(filename, "./ResultC_done.dat");
	std::ifstream ResultCbD_file(filename, std::ios::binary);
    ResultCbD_file.read((char*)&N,sizeof(uint_fast32_t));
    ResultCbD_file.read((char*)P,N*sizeof(boost::qvm::vec<double,3>));
    ResultCbD_file.read((char*)Pdata,N*sizeof(StabilityPointType));
    ResultCbD_file.close();
	for(uint_fast32_t i=0; i<N; ++i)
    {
        Pdata[i].PointType = 0;
    }/**/


    CheckBorderMD();
    std::cerr<<"Good\n";
    //std::cin.get();

    strcpy(filename, "./ResultC.dat");
    std::ofstream ResultCb_file(filename, std::ios::binary);
    ResultCb_file.write((char*)&N,sizeof(uint_fast32_t));
    ResultCb_file.write((char*)P,N*sizeof(boost::qvm::vec<double,3>));
    ResultCb_file.write((char*)Pdata,N*sizeof(StabilityPointType));
    ResultCb_file.close();

    strcpy(filename, "./ResultC.txt");
    std::ofstream ResultC_file(filename);
    strcpy(filename, "./ResultUnStable.txt");
    std::ofstream ResultUStable_file(filename);
    strcpy(filename, "./ResultnStable.txt");
    std::ofstream ResultStable_file(filename);
    for(uint_fast32_t i = 0; i < N; ++i)
    {
        ResultC_file<<P[i].a[0]<<" "<<P[i].a[1]<<" "<<P[i].a[2]<<" "<<
		int(Pdata[i].Stability)<<" "<<Pdata[i].StabilitySteps<<" "<<
		Pdata[i].StabilityTime<<" "<<Pdata[i].PotEnergy<<" "
		<<"  a"<<"\n";

		if(int(Pdata[i].Stability)==1)
            ResultUStable_file<<P[i].a[0]<<" "<<P[i].a[1]<<" "<<P[i].a[2]<<" "<<
			Pdata[i].StabilitySteps<<" "<<Pdata[i].StabilityTime<<" "<<
			Pdata[i].PotEnergy<<"\n";

        else if(int(Pdata[i].Stability)==2)
            ResultStable_file<<P[i].a[0]<<" "<<P[i].a[1]<<" "<<P[i].a[2]<<" "<<
			Pdata[i].StabilitySteps<<" "<<Pdata[i].StabilityTime<<"\n";
    }
    ResultC_file.close();
    findBorder();
    //std::cerr<<"Q1 "<<NbL<<"\n";
    //std::cin.get();
    preciseBorderCreate();
    //std::cerr<<"Q2 "<<NbL<<"\n";
    //std::cin.get();
    strcpy(filename, "./ResultBc.txt");
    std::ofstream ResultBc_file(filename);
    strcpy(filename, "./Result2Bc.txt");
    std::ofstream Result2Bc_file(filename);
    for(uint_fast32_t i=0;i<NbL;++i)
    {

        // /toElastic2D(CB.Sbc[i],EM);
        // if(EM.Ex<0)
        ResultBc_file<<Pbc[i].a[0]-1.<<" "<<Pbc[i].a[1]-1.<<" "<<Pbc[i].a[2]-1.
			<<" "<<int(Pbcdata[i].Stability)<<" "<<Pbcdata[i].StabilitySteps
			<<" "<<Pbcdata[i].StabilityTime<<" "<<Pbcdata[i].PotEnergy<<"\n";
    }
    ResultBc_file.close();
    for(uint_fast32_t i=0;i<NbL;++i)
    {

        ///toElastic2D(CB.Sbc[i],EM);
        //if(EM.Ex<0)
        Result2Bc_file<<Pb[i].a[0]-1.<<" "<<Pb[i].a[1]-1.<<" "<<Pb[i].a[2]-1.
			<<" "<<int(Pbdata[i].Stability)<<" "<<Pbdata[i].StabilitySteps
			<<" "<<Pbdata[i].StabilityTime<<" "<<Pbcdata[i].PotEnergy<<"\n";
    }
    Result2Bc_file.close();
    //std::cerr<<"Q3\n";
    //std::cin.get();
    NbLAll = NbL;
    for(;NbL>0;)
    {
        checkStabilityMPIMD(Pbc,Pbcdata,Sbc,NbL);
        for(uint_fast32_t i=0; i<NbL; ++i)
        {
            if(Pbcdata[i].Stability==1)
                ResultUStable_file<<Pbc[i].a[0]<<" "<<Pbc[i].a[1]<<" "
				<<Pbcdata[i].StabilitySteps<<" "<<Pbcdata[i].StabilityTime<<"\n";

            else if(Pbcdata[i].Stability==2)
                ResultStable_file<<Pbc[i].a[0]<<" "<<Pbc[i].a[1]<<" "<<
				Pbcdata[i].StabilitySteps<<" "<<Pbcdata[i].StabilityTime<<"\n";
        }
        preciseBorderReCreate();
        std::cerr<<NbL<<" "<<NbLAll<<"\n";
        //std::cin.get();
    }
    ResultUStable_file.close();
    ResultStable_file.close();
    NbL = NbLAll;

    strcpy(filename, "./ResultBcS.txt");
    std::ofstream ResultBcL_file(filename);
    strcpy(filename, "./ResultBcVS.txt");
    std::ofstream ResultBcLS_file(filename);
    strcpy(filename, "./ResultBcVU.txt");
    std::ofstream ResultBcLU_file(filename);
    for(uint_fast32_t i=0;i<NbL;++i)
    {
        ResultBcL_file<<Pbc[i].a[0]<<" "<<Pbc[i].a[1]<<" "<<Pbc[i].a[2]
		<<" "<<Pbcdata[i].StabilitySteps<<"\n";

        ResultBcLS_file<<Pb[2*i].a[0]<<" "<<Pb[2*i].a[1]<<" "<<Pb[2*i].a[2]
		<<" "<<Pbdata[2*i].StabilitySteps<<"\n";

        ResultBcLU_file<<Pb[2*i+1].a[0]<<" "<<Pb[2*i+1].a[1]<<" "<<
		Pb[2*i+1].a[2]<<" "<<Pbdata[2*i+1].StabilitySteps<<"\n";
    }
    ResultBcL_file.close();
    ResultBcLS_file.close();
    ResultBcLU_file.close();/**/
    stopMPI();
}


void CalcBorder::taskSRMPI0()
{
    int flag;
    MPI_Status statusP;
    MPI_Iprobe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &flag, &statusP);
    if(flag==1)
    {
        TaskType3D Task = getNextTask();
        StabilityPointType Result;
        MPI_Status status;
        MPI_Recv(&Result,sizeof(StabilityPointType), MPI_BYTE, statusP.MPI_SOURCE, statusP.MPI_TAG, MPI_COMM_WORLD, &status);
        MPI_Send(&Task, sizeof(TaskType3D), MPI_BYTE, statusP.MPI_SOURCE, 0, MPI_COMM_WORLD);
        if(Result.i>=0)
        {
            PdataR[Result.i] = Result;
            ++NReceived;
            //std::cerr<<"Resived "<<statusP.MPI_SOURCE<<" "<<Result.i<<" "<<NReceived<<" "<<NSR<<"\n";
            //std::cerr<<"Resived "<<Result.i<<" "<<int(Result.Stability)<<" "<<Result.StabilitySteps<<"\n";
            //std::cin.get();
        }
    }
}

TaskType3D CalcBorder::getNextTask()
{
    TaskType3D Task;
	double e;
    if(NSended<NSR)
    {
        boost::qvm::set_identity(Task.D);
        Task.D.a[0][0]+=PS[NSended].a[0];
        Task.D.a[1][1]+=PS[NSended].a[1];
        Task.D.a[2][2]+=PS[NSended].a[2];

        Task.D += dG0;
        Task.Ek = 1e-6;
        //Task.Step = NStep;
		e = (Task.D.a[0][0]+Task.D.a[1][1]+Task.D.a[2][2])*MC_1d3;
		Task.Step = (NStepB-NStepA)*(e-0.4)/0.6+NStepA;
        Task.i = NSended;
        ++NSended;
    }
    else
    {
        boost::qvm::set_identity(Task.D);
        Task.Ek = 1e-6;
        Task.Step = 100;
        Task.i = -1;
    }
    return Task;
}

void taskSRMPIi(MD *MDC)
{
    char filename[256]="";
    sprintf(filename, "./result/ResultC_%08u.txt", iproc);
    std::ofstream ResultiC_file(filename);
    sprintf(filename, "./result/ResultUnStable_%08u.txt", iproc);
    std::ofstream ResultiUStable_file(filename);
    sprintf(filename, "./result/ResultnStable_%08u.txt", iproc);
    std::ofstream ResultiStable_file(filename);
    sprintf(filename, "./result/ResultC_%08u.dat", iproc);
    std::ofstream ResultiCbin_file(filename,std::ios::binary);
    TaskType3D Task;
    int flag;
    MPI_Status status;
    StabilityPointType Result;
    Result.i = -1;
    MPI_Recv(&Task, sizeof(TaskType3D), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
    //std::cerr<<"Resived "<<iproc<<" "<<Task.i<<" "<<Task.Step<<" "<<Task.D.a[0][0]<<" "<<Task.D.a[0][1]<<" "<<Task.D.a[1][0]<<" "<<Task.D.a[1][1]<<"\n";
    uint_fast32_t i;



    //MDC->reachSpecifiedStress();
    //MDC->saveLattice();
    for(;;)
    {
        if(Task.i>=0)
        {
            //MDC->eXYMax = 2.0*MDC->LV[1].a[1]*Task.D.a[1][1]/(MDC->LV[0].a[0]*Task.D.a[0][0]);
            //Task.D.a[0][1] = (MDC->eXYMax<Task.D.a[0][1])?Task.D.a[0][1]:(Task.D.a[0][1]-MDC->eXYMax*int_fast32_t(Task.D.a[0][1]/MDC->eXYMax));
            //std::cerr<<"Task.Di "<<iproc<<" "<<Task.D.a[0][0]<<" "<<Task.D.a[1][1]<<" "<<Task.D.a[2][2]<<" "<<MDC->eXYMax<<"\n";
            MDC->setMPP();
            MDC->createLattice3DT(Task.D);
            MDC->setTemperatureNormal();
            MDC->deformLattice(Task.D);
            ///MDC->saveLattice();
            MDC->createIMatrix();
            MDC->calculateStateIM();
            MDC->Ek0 = MDC->Ek;
            MDC->Ep0 = MDC->Ep;

			double PotEnergy = MDC->Ep/MDC->N;
			
            if(MDC->NM==1)goto SkipCalci;
            for(i=0; i<Task.Step; ++i)
            {

                MDC->calculateForcesIM();
                //std::cerr<<"Q10"<<"\n";
                MDC->calculateIncrements();
                //std::cerr<<"Q "<<MDC->t<<" "<<MDC->Ek*MDC->_1d_N<<" "<<MDC->Ek0*MDC->_1d_N<<" "<<MDC->Ep<<" "<<MDC->P_Cmax<<" "<<MDC->P_C<<"\n";
                //std::cin.get();
				PotEnergy = MDC->Ep/MDC->N;
                if(MDC->Ek>1.01*MDC->Ek0)
                {
					PotEnergy = MDC->Ep/MDC->N;
					//std::cerr<<"NOT Stable!!!!!\n";
                    break;
                }
            }
SkipCalci:;
            Result.i = Task.i;
            Result.Stability = (i==Task.Step)?2:1;
            Result.StabilityTime = MDC->t;
            Result.StabilitySteps = i;
			Result.NumberBonds = MDC->NM;
			Result.PotEnergy = PotEnergy;
            ResultiC_file<<Task.D.a[0][0] - 1.<<" "<<Task.D.a[1][1] - 1.<<" "<<
			Task.D.a[2][2] - 1.<<" "<<int(Result.Stability)<<" "<<
			Result.StabilitySteps<<" "<<Result.StabilityTime
			<<" "<<Result.PotEnergy<<" "<<Result.NumberBonds<<" "<<
			int(Result.Stability)<<"\n";

            if(int(Result.Stability)==1)
                ResultiUStable_file<<Task.D.a[0][0]<<" "<<Task.D.a[1][1]
				<<" "<<Task.D.a[2][2]<<" "<<Result.StabilitySteps
				<<" "<<Result.StabilityTime<<" "<<Result.PotEnergy<<"\n";

            else if(int(Result.Stability)==2)
                ResultiStable_file<<Task.D.a[0][0]<<" "<<Task.D.a[1][1]
				<<" "<<Task.D.a[2][2]<<" "<<Result.StabilitySteps
				<<" "<<Result.StabilityTime<<"\n";

            ResultiCbin_file.write((char*)&Task,sizeof(TaskType3D));
            ResultiCbin_file.write((char*)&Result,sizeof(StabilityPointType));
        }
        else
        {
            MPI_Iprobe(0, 0, MPI_COMM_WORLD, &flag, &status);
            //std::cerr<<"Try probe "<<iproc<<" "<<flag<<"\n";
            if(flag==1)
            {
                MPI_Recv(&Task, sizeof(TaskType3D), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);

                //std::cerr<<"Resv "<<iproc<<" "<<Task.i<<"\n";
                if(Task.i==-2)break;
                else if(Task.i>=0)continue;
            }
            else continue;
        }
        MPI_Send(&Result,sizeof(StabilityPointType), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
        MPI_Recv(&Task, sizeof(TaskType3D), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &status);
        ///MDC->restoreLattice();
    }
    ResultiUStable_file.close();
	ResultiStable_file.close();
	ResultiC_file.close();
	ResultiCbin_file.close();
}

void CalcBorder::stopMPI()
{
    TaskType3D Task;
    boost::qvm::set_identity(Task.D);
    Task.Ek = 1e-6;
    Task.Step = 100;
    Task.i = -2;
    MPI_Status status;
    for(int_fast32_t i=1; i<Nproc; ++i)
    {
        MPI_Send(&Task, sizeof(TaskType3D), MPI_BYTE, i, 0, MPI_COMM_WORLD);
        //std::cerr<<"Sended "<<i<<" "<<Task.i<<" "<<NReceived<<" "<<NSR<<"\n";
        //std::cin.get();
    }
}

void CalcBorder::checkStabilityMPIMD(boost::qvm::vec<double,3> *Pvar, StabilityPointType *Pdvar, boost::qvm::mat<double,3,3> *Svar, const uint_fast32_t &Nvar)
{
    //createPoints3D();
    //createIMatrix();
    char filename[256]="";
    sprintf(filename, "./result/ResultC_%08u.txt", iproc);
    std::ofstream ResultiC_file(filename);
    sprintf(filename, "./result/ResultUnStable_%08u.txt", iproc);
    std::ofstream ResultiUStable_file(filename);
    sprintf(filename, "./result/ResultnStable_%08u.txt", iproc);
    std::ofstream ResultiStable_file(filename);
    sprintf(filename, "./result/ResultC_%08u.dat", iproc);
    std::ofstream ResultiCbin_file(filename,std::ios::binary);
    MPI_Status status, statusP;
    TaskType3D Task;
    StabilityPointType Result;
    PS = Pvar;
    PdataR = Pdvar;
    NSR = Nvar;
    //int_fast32_t flag;
    //uint_fast32_t iP;
    //std::cerr<<"Q4\n";
    MDC->setMPP();
    MDC->createLattice3DT();
    //std::cerr<<"Q2\n";

    ///MDC->saveLattice();
    MDC->createIMatrix();
    //std::cerr<<"Q0\n";
    MDC->calculateActualStiffness();
    //ElasticModules3D b;toElastic3D(MDC->StiffnessVoigt, b);
    MDC->reachSpecifiedStress(); //поджатие связей
    dG0 = MDC->dG0;
    std::cerr<<"dG0 "<<dG0.a[0][0]<<" "<<dG0.a[0][1]<<" "<<dG0.a[0][2]
    <<" "<<dG0.a[1][0]<<" "<<dG0.a[1][1]<<" "<<dG0.a[1][2]<<" "<<dG0.a[2][0]<<" "<<dG0.a[2][1]<<" "<<dG0.a[2][2]<<"\n";
    //std::cerr<<"Elastic\n"
    //<<b.Ex<<" "<<b.Ey<<" "<<b.Ez<<" "<<b.K<<"\n"
    //<<b.Gxy<<" "<<b.Gyz<<" "<<b.Gzx<<"\n"
    //<<b.nuxy<<" "<<b.nuxz<<" "<<b.nuyx<<" "<<b.nuyz<<" "<<b.nuzx<<" "<<b.nuzy<<"\n";
    //std::cin.get();

    NSended = 0;
    NReceived = 0;
    for(int_fast32_t i=1; i<Nproc; ++i)
    {
        Task = getNextTask();
        MPI_Status status;
        MPI_Send(&Task, sizeof(TaskType3D), MPI_BYTE, i, 0, MPI_COMM_WORLD);
        //std::cerr<<"Sended "<<i<<" "<<Task.i<<" "<<NReceived<<" "<<NSR<<"\n";
        //std::cin.get();
    }
    uint_fast32_t i;
    for(;;)
    {
        if(NSended<NSR)
        {
            Task = getNextTask();
            //MDC->eXYMax = 2.0*MDC->LV[1].a[1]*Task.D.a[1][1]/(MDC->LV[0].a[0]*Task.D.a[0][0]);
            //Task.D.a[0][1] = (fabs(Task.D.a[0][1])<MDC->eXYMax)?Task.D.a[0][1]:(Task.D.a[0][1]-MDC->eXYMax*int_fast32_t(Task.D.a[0][1]/MDC->eXYMax));
            std::cerr<<"Task.D0 "<<NSended<<" "<<NReceived<<" "<<NSR<<" "
			<<Task.i<<" "<<Task.D.a[0][0]-1<<" "<<Task.D.a[1][1]-1
			<<" "<<Task.D.a[2][2]-1<<"\n";

            MDC->createLattice3DT(Task.D);
            MDC->setTemperatureNormal();
            MDC->deformLattice(Task.D);
            //std::cerr<<"Q8"<<"\n";
            MDC->createIMatrix();
            //std::cerr<<"Q9"<<"\n";
            MDC->calculateStateIM();
            //std::cerr<<"Stress "<<MDC->Stress.a[0][0]<<" "<<MDC->Stress.a[0][1]<<" "<<MDC->Stress.a[1][0]<<" "<<MDC->Stress.a[1][1]<<"\n";
            //std::cerr<<"E "<<MDC->Ek*MDC->_1d_N<<" "<<MDC->Ep*MDC->_1d_N<<" "<<MDC->P_Cmax<<" "<<MDC->P_C<<"\n";

            MDC->Ek0 = MDC->Ek;
            MDC->Ep0 = MDC->Ep;

			double PotEnergy = MDC->Ep / MDC->N;

			// Result.PotEnergy = Energy;
            if(MDC->NM==1)goto SkipCalc0;
            for(i=0; i<Task.Step; ++i)
            {

                MDC->calculateForcesIM();
                //std::cerr<<"Q10"<<"\n";
                MDC->calculateIncrements();
                //std::cerr<<"Q "<<MDC->t<<" "<<MDC->Ek*MDC->_1d_N<<" "<<MDC->Ek0*MDC->_1d_N<<" "<<MDC->Ep<<"\n";
                //std::cin.get();
                if(i%100==0)
                {
                    taskSRMPI0();
                }
                if(MDC->Ek>1.01*MDC->Ek0)
                {
					PotEnergy = MDC->Ep / MDC->N;
					// Result.PotEnergy = PotEnergy;
					break;
                }
            }
SkipCalc0:;
            Result.i = Task.i;
            Result.Stability = (i==Task.Step)?2:1;
            Result.StabilityTime = MDC->t;
            Result.StabilitySteps = i;
			Result.NumberBonds = MDC->NM;
			Result.PotEnergy = PotEnergy;
			// Result.PotEnergy = Energy;

            PdataR[Result.i] = Result;
            ResultiC_file<<Task.D.a[0][0]<<" "<<Task.D.a[1][1]<<" "
			<<Task.D.a[2][2]<<" "<<int(Result.Stability)<<" "<<
			Result.StabilitySteps<<" "<<Result.StabilityTime<<" "<<
			Result.PotEnergy<<" "<<Result.NumberBonds<<" num "<<"\n";

            if(int(Result.Stability)==1)
                ResultiUStable_file<<Task.D.a[0][0]<<" "<<Task.D.a[1][1]<<" "
				<<Task.D.a[2][2]<<" "<<Result.StabilitySteps<<" "
				<<Result.StabilityTime<<" "<<Result.PotEnergy<<"\n";

            else if(int(Result.Stability)==2)
                ResultiStable_file<<Task.D.a[0][0]<<" "<<Task.D.a[1][1]<<" "
				<<Task.D.a[2][2]<<" "<<Result.StabilitySteps<<" "
				<<Result.StabilityTime<<" "<<Result.PotEnergy<<"\n";

            ResultiCbin_file.write((char*)&Task,sizeof(TaskType3D));
            ResultiCbin_file.write((char*)&Result,sizeof(StabilityPointType));

            ++NReceived;

            ///MDC->restoreLattice();
            //std::cerr<<"Fin1 "<<Result.i<<" "<<int(Result.Stability)<<" "<<Result.StabilitySteps<<"\n";
            //std::cin.get();
        }
        else if(NReceived<Nvar)
        {
            //std::cerr<<"ERROR! "<<NSR<<" "<<NSended<<" "<<NReceived<<"\n";

            taskSRMPI0();
            //std::cin.get();
        }
        else
        {
            break;
        }
    }
    ResultiUStable_file.close();
	ResultiStable_file.close();
	ResultiC_file.close();
	ResultiCbin_file.close();
}

void CalcBorder::preciseBorderReCreate()
{
    boost::qvm::vec<double,3> tempP;
    StabilityPointType tempPdata;
    for(int_fast32_t i=0; i<NbL; ++i)
    {
        if(boost::qvm::mag_sqr(Pb[2*i]-Pb[2*i+1])<BPBP)
        {
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
        }
        else if(Pbcdata[i].Stability==1)
        {
            Pb[2*i+1] = Pbc[i];
            Pbc[i] = 0.5*(Pb[2*i]+Pb[2*i+1]);
            Pbdata[2*i+1] = Pbcdata[i];
        }
        else if(Pbcdata[i].Stability==2)
        {
            Pb[2*i] = Pbc[i];
            Pbc[i] = 0.5*(Pb[2*i]+Pb[2*i+1]);
            Pbdata[2*i] = Pbcdata[i];
        }
        std::cerr<<"BCL "<<boost::qvm::mag_sqr(Pb[2*i]-Pb[2*i+1])<<"\n";
    }
}
