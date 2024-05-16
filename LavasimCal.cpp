#include "LavasimCal.h"
#define DZERO 0.000001

using namespace std;

bool Lavasim::setNeibor(){
    double RScale = this->m_Para.m_RScale;
    /*
    if(this->m_Neibor)
        delete[] this->m_Neibor;
    this->m_Neibor = new int[this->m_Para.m_nRow*this->m_Para.m_nCol];
    //init randomlist
    double* randomlist = new double[this->m_Para.m_nRow*this->m_Para.m_nCol*2];
    srand(time(NULL));
    for (int i = 0; i<this->m_Para.m_nRow*this->m_Para.m_nCol*2; i++)
        randomlist[i] = (double)rand() / RAND_MAX;
    double threshold = pow(RScale, 2);
    //      map: (0,0),(0,1),(0,2)(0,3)...
    //765    (y,x) (1,0),(1,1),(1,2)(1,3)...
    //0c4          (2,0),(2,1),(2,2)(2,3)...
    //123          ......
    
    int delta_x[8] = {-1, -1, 0, +1, +1, +1, 0, -1};
    int delta_y[8] = {0, +1, +1, +1, 0, -1, -1, -1};
    int mark[8] = {1, 2, 4, 8, 16, 32, 64, 128};

    for(int i=0; i<this->m_Para.m_nRow*this->m_Para.m_nCol; ++i){
        int res = 0;
        int c_y = i/this->m_Para.m_nCol;
        int c_x = i-c_y*this->m_Para.m_nCol;
        //std::cerr<<i<<" "<<c_x<<" "<<c_y<<"INFUNC\n";
        for(int j=0; j<8; j++){
            int t_x = c_x + delta_x[j];
            int t_y = c_y + delta_y[j];
            if( !(t_x>=0 && t_x<this->m_Para.m_nCol && t_y>=0 && t_y<this->m_Para.m_nRow) )
                continue;
            //std::cerr<<t_x<<" "<<t_y<<std::endl;
            int idx = t_y*this->m_Para.m_nCol+t_x;
            double rand_cx=randomlist[i];
            double rand_cy=randomlist[i+this->m_Para.m_nRow*this->m_Para.m_nCol];
            double rand_tx=randomlist[idx];
            double rand_ty=randomlist[idx+this->m_Para.m_nRow*this->m_Para.m_nCol];
            double res_x = 1.0 + delta_x[j]*(rand_tx-rand_cx);
            double res_y = 1.0 + delta_y[j]*(rand_ty-rand_cy);
            if( (res_x*res_x+res_y*res_y)<=threshold )
                res+=mark[j];
        }
        this->m_Neibor[i]=res;
    }
    */
    return true;
}

void Lavasim::outputParam(){
    string filename = (this->m_OutputFolder+"params");
    std::ofstream ofs;
    ofs.open(filename.data(), std::ios::out);
    if(!ofs.is_open())
        return;
    time_t now=time(0);

    tm* ltime = localtime(&now);
    //asctime(ltime)
    ofs<<"Time: "<<asctime(ltime);
    ofs<<"*************************************************************"<<std::endl;
    //ctl
    ofs<<"RScale: "<<this->m_Para.m_RScale<<std::endl;
    ofs<<"CtlDelta: "<<this->m_Para.m_CtlDelta<<std::endl;
    ofs<<"nRow: "<<this->m_Para.m_nRow<<std::endl;
    ofs<<"nCol: "<<this->m_Para.m_nCol<<std::endl;
    ofs<<"CellSize: "<<this->m_Para.m_CellSize<<std::endl;
    ofs<<"TimeStepMax: "<<this->m_Para.m_TimeStepMax<<std::endl;
    ofs<<"IsAreaCorrection: "<<this->m_Para.m_AreaCorrection<<std::endl;
    ofs<<"*************************************************************"<<std::endl;
    //lava
    ofs<<"IsViscosityStatic: "<<this->m_Para.m_IsVStatic<<std::endl;
    if(this->m_Para.m_IsVStatic)
        ofs<<"Viscosity: "<<this->m_Para.m_Vis<<std::endl;
    else{
        ofs<<"ViscosityParaA: "<<this->m_Para.m_VisParaA<<std::endl;
        ofs<<"ViscosityParaB: "<<this->m_Para.m_VisParaB<<std::endl;
        ofs<<"ViscosityParaC: "<<this->m_Para.m_VisParaC<<std::endl;
    }
    ofs<<"IsDensityStatic: "<<this->m_Para.m_IsDStatic<<std::endl;
    if(this->m_Para.m_IsDStatic)
        ofs<<"Density: "<<this->m_Para.m_Den<<std::endl;
    else{
        ofs<<"DensityParaA: "<<this->m_Para.m_DenParaA<<std::endl;
        ofs<<"DensityParaB: "<<this->m_Para.m_DenParaB<<std::endl;
        ofs<<"DensityParaC: "<<this->m_Para.m_DenParaC<<std::endl;
    }
    ofs<<"IsTransform: "<<this->m_Para.m_IsTransform<<std::endl;
    if(this->m_Para.m_IsTransform)
        ofs<<"TransformTemp: "<<this->m_Para.m_TransformTemp<<std::endl;
    ofs<<"*************************************************************"<<std::endl;
    //envir
    ofs<<"Gravity: "<<this->m_Para.m_GravAcc<<std::endl;
    ofs<<"AtmTemperature: "<<this->m_Para.m_AtmTemperature<<std::endl;
    ofs<<"AtmDen: "<<this->m_Para.m_AtmDen<<std::endl;
    ofs<<"AtmSpecificHeat: "<<this->m_Para.m_AtmSpecificHeat<<std::endl;
    ofs<<"Aa: "<<this->m_Para.m_Aa<<std::endl;
    ofs<<"Ka: "<<this->m_Para.m_Ka<<std::endl;
    ofs<<"Va: "<<this->m_Para.m_Va<<std::endl;
    ofs<<"Gamma: "<<this->m_Para.m_Gamma<<std::endl;
    ofs<<"Lambda: "<<this->m_Para.m_Lambda<<std::endl;
    ofs<<"Cv: "<<this->m_Para.m_Cv<<std::endl;
    ofs<<"Emissivity: "<<this->m_Para.m_Emi<<std::endl;
    ofs<<"IsCoupled: "<<this->m_Para.m_IsCoupled<<std::endl;
    ofs<<"*************************************************************"<<std::endl;
    ofs.close();
    //output neibor_list
    filename = (this->m_OutputFolder + "neibor");
    ofs.open(filename.data(), std::ios::out);
    if(!ofs.is_open())
        return;
    ofs<<asctime(ltime);
    ofs<<this->m_Para.m_nRow<<" "<<this->m_Para.m_nCol<<std::endl;
    for(int i=0; i<this->m_Para.m_nRow; i++){
        for(int j=0; j<this->m_Para.m_nCol; j++){
            ofs<<this->m_Neibor[i*this->m_Para.m_nCol+j]<<" ";
        }
        ofs<<std::endl;
    }
    ofs.close();
    return;
}

void Lavasim::checkpoint(double timeStep, int sumTime, int sumCount) {
    const char* category[3] = {"hdata", "zdata", "tdata"};
    double* datas[3] = {this->m_HList, this->m_ZList, this->m_TList};
    for(int i=0; i<3; i++){
		string filename = (this->m_OutputFolder + category[i] + "_" + std::to_string(sumTime));// +"_" + std::to_string(sumCount));
        std::ofstream ofs;
        ofs.open(filename, std::ios::out);
        if(!ofs.is_open())
            return;
        time_t now=time(0);
        tm* ltime = localtime(&now);
        //asctime(ltime)
        ofs<<"Time: "<<asctime(ltime);
        ofs<<"Sum time: "<<sumTime<<" Sum count: "<<sumCount<<std::endl;//<<sumCount<<" currentTimeStep: "<<timeStep<<std::endl;
        ofs<<"Cell size: "<<this->m_Para.m_CellSize<<std::endl;
        ofs<<"nRow: "<<this->m_Para.m_nRow<<" nCol: "<<this->m_Para.m_nCol<<std::endl;
        for(int r=0; r<this->m_Para.m_nRow; r++){
            for(int c=0; c<this->m_Para.m_nCol; c++){
                ofs<<datas[i][r*this->m_Para.m_nCol+c]<<" ";
            }
            ofs<<std::endl;
        }
        ofs.close();
    }
}

__host__
double Lavasim::getVolcano(double nowTime){
    double* tList = this->m_VolcanoData.tList;
    double* vList = this->m_VolcanoData.vList;
    int cnt = this->m_VolcanoData.cnt;
    double effusion;
    int i = 0;
    while(i < this->m_VolcanoData.cnt && nowTime > tList[i])
        ++i;
    if(i == 0 || i == cnt)
        effusion = 0.0;
    else
        effusion = vList[i - 1] + (vList[i] - vList[i - 1]) * (nowTime - tList[i - 1]) / (tList[i] - tList[i - 1]);
    return effusion;
}

__host__
void Lavasim::initCalData(double** qflux, double** hflux, double** zData, double** hData, double** tData, double** maxQ, \
    int** vents, int** neibor, ParaData** param){
    //qflux, hflux, maxQ.
    cudaMalloc(qflux, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
    cudaMemset(*qflux, 0.0, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
    cudaMalloc(hflux, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
    cudaMemset(*hflux, 0.0, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
    cudaMalloc(maxQ, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
    cudaMemset(*maxQ, -0xaa, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
    //zData, tData, hData.
    cudaMalloc(zData, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
    cudaMemcpy(*zData, this->m_ZList, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double), cudaMemcpyHostToDevice);
    cudaMalloc(hData, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
    cudaMemcpy(*hData, this->m_HList, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double), cudaMemcpyHostToDevice);
    cudaMalloc(tData, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
    cudaMemcpy(*tData, this->m_TList, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double), cudaMemcpyHostToDevice);
    cudaMalloc(vents, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(int));
    cudaMemcpy(*vents, this->m_Vent, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(int), cudaMemcpyHostToDevice);
    cudaMalloc(neibor, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(int));
    cudaMemcpy(*neibor, this->m_Neibor, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(int), cudaMemcpyHostToDevice);
    cudaMalloc(param, sizeof(ParaData));
    cudaMemcpy(*param, &this->m_Para, sizeof(ParaData), cudaMemcpyHostToDevice);
    return;
}

__global__
void checkBoundary(int* limit_dev, int* changed, int nRow, int nCol){
    //check Boundary.
    if(changed[0])
        limit_dev[0] = (limit_dev[0]-1)>=0? (limit_dev[0]-1):0;
    if(changed[1])
        limit_dev[1] = (limit_dev[1]+1)<(nRow-1)? (limit_dev[1]+1):(nRow-1);
    if(changed[2])
        limit_dev[2] = (limit_dev[2]+1)<(nCol-1)? (limit_dev[2]+1):(nCol-1);
    if(changed[3])
        limit_dev[3] = (limit_dev[3]-1)>=0? (limit_dev[3]-1):0;
    return;
}

__host__
void Lavasim::doCal(){
    string log = this->m_OutputFolder + "lavasim.log";
    ofstream ofs;
    ofs.open(log.data(), std::ios::app);
    ofs << "Outputting params...\n";
    outputParam();
    checkpoint(0.0 , 0.0, 0);

    //initial data.
    double* qflux = NULL; // init 0.0
    double* hflux = NULL; // init 0.0
    double* zData = NULL;
    double* hData = NULL;
    double* tData = NULL;
    double* maxDt = NULL; // init 0x3f3f3f3f
    int* vents = NULL;
    int* neibor = NULL;
    ParaData* param=NULL;
    this->initCalData(&qflux, &hflux, &zData, &hData, &tData, &maxDt, &vents, &neibor, &param);
    ofs << "Finished Initialization.\n";

    int nowCount = 0;
    double nowTime = 0.0;
    double timeStepMax = this->m_Para.m_TimeStepMax;
    double timeStep = timeStepMax;
    double ctlDelta = this->m_Para.m_CtlDelta;
    double sumTime = this->m_Para.m_SumTime;
    double outputTime = this->m_Para.m_OutputTime;
    double checkTime = 0.0;

    time_t now=time(0);
    tm ltime_base = *(localtime(&now));
    ofs << "Start Computing...\n";
    ofs << "Current system time: " << asctime(&ltime_base);

    //init bounder.
    // 0-left, 1-down, 2-right, 3-up.
    int limit[4] = {this->m_Para.m_nCol, 0, 0, this->m_Para.m_nRow};
    for(int i=0; i<this->m_Para.m_nRow*this->m_Para.m_nCol; ++i){
        if(this->m_Vent[i]!=-1){
            int ir = i/this->m_Para.m_nCol;
            int ic = i - ir*this->m_Para.m_nCol;
            limit[0] = std::min<int>(limit[0], ic);
            limit[1] = std::max<int>(limit[1], ir);
            limit[2] = std::max<int>(limit[2], ic);
            limit[3] = std::min<int>(limit[3], ir);
        }
    }
    int* limit_dev=NULL;
    cudaMalloc(&limit_dev, 4*sizeof(int));
    cudaMemcpy(limit_dev, limit, 4*sizeof(int), cudaMemcpyHostToDevice);

    int* changed=NULL;
    cudaMalloc(&changed, 4*sizeof(int));
    cudaMemset(changed, 0, 4*sizeof(int));

    ofs << "Initial Boundary: left-"<<limit[0]<<" down-"<<limit[1]<<" right-"<<limit[2]<<" up-"<<limit[3]<<std::endl;
    //set up the parallel parameters.
    dim3 grid(1 +(this->m_Para.m_nRow-1)/16, 1 + (this->m_Para.m_nCol-1)/16, 1), block(16, 16, 1);

    while(nowTime<=sumTime){
        if(nowCount % 10000 == 0){
            ofs << "=======================================================\n";
            ofs <<"Now count:" << nowCount << " Now time:" << nowTime << " Timestep:" << timeStep << std::endl;
            now = time(0);
            tm* ltime_now = localtime(&now);
            ofs << "Current system time: " << asctime(ltime_now);
            int lastSec = (ltime_now->tm_mday - ltime_base.tm_mday) * 86400 + (ltime_now->tm_hour - ltime_base.tm_hour) * 3600\
                        + (ltime_now->tm_min-ltime_base.tm_min) * 60 + (ltime_now->tm_sec - ltime_base.tm_sec);
            ofs << "Lasting time: " << lastSec << "s\n";
        }
        //get effusion of vents.
        double effusion = getVolcano(nowTime);
        double lavaTemp = this->m_VolcanoData.temperature;
        calEachFlux<<<grid, block>>>(qflux, hflux, hData, zData, tData, maxDt, neibor, limit_dev, param);
        cudaDeviceSynchronize();

        getTimeStep<<<1, this->m_Para.m_nRow>>>(maxDt, this->m_Para.m_nRow, this->m_Para.m_nCol);

        cudaMemcpy(&timeStep, maxDt, sizeof(double), cudaMemcpyDeviceToHost);
        timeStep = timeStep*ctlDelta;
	if(nowCount % 1000 == 0){
		std::cerr<<nowCount<<" "<<nowTime<<" "<<timeStep<<std::endl;
	}
        timeStep = std::min<double>(timeStep, timeStepMax);

        updateData<<<grid, block>>>(qflux, hflux, zData, hData, tData, vents, limit_dev, changed, param, effusion, timeStep, lavaTemp);
        cudaDeviceSynchronize();

        checkBoundary<<<1, 1>>>(limit_dev, changed, this->m_Para.m_nRow, this->m_Para.m_nCol);

        //cudaMemcpy(limit, limit_dev, sizeof(int)*4, cudaMemcpyDeviceToHost);
        //std::cerr<<limit[0]<<" "<<limit[1]<<" "<<limit[2]<<" "<<limit[3]<<"\n";
        /*if(nowCount % 1000 == 0){
            cudaMemcpy(limit, limit_dev, sizeof(int)*4, cudaMemcpyDeviceToHost);
            std::cerr<<limit[0]<<" "<<limit[1]<<" "<<limit[2]<<" "<<limit[3]<<"\n";
        }*/

        //init the flux and maxQ array.
        cudaMemset(qflux, 0.0, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
        cudaMemset(hflux, 0.0, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
        cudaMemset(maxDt, -0xaa, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double));
        cudaMemset(changed, 0, 4*sizeof(int));

        //update time and count.
        ++nowCount;
        nowTime += timeStep;
        checkTime += timeStep;
        if(checkTime >= outputTime){
			ofs << "=======================================================\n";
            ofs << "Now count:" << nowCount << " Now time:" << nowTime << " Timestep:" << timeStep << std::endl;
            cudaMemcpy(this->m_ZList, zData, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(this->m_HList, hData, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double), cudaMemcpyDeviceToHost);
            cudaMemcpy(this->m_TList, tData, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double), cudaMemcpyDeviceToHost);
            this->checkpoint(timeStep, nowTime, nowCount);
            checkTime -= outputTime;
        }
    }

    //Finish calculation.
    cudaMemcpy(this->m_ZList, zData, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(this->m_HList, hData, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double), cudaMemcpyDeviceToHost);
    cudaMemcpy(this->m_TList, tData, this->m_Para.m_nRow*this->m_Para.m_nCol*sizeof(double), cudaMemcpyDeviceToHost);
    checkpoint(timeStep, nowTime, nowCount);
    ofs << "Finished Computing!!!\n";
    ofs.close();
    cudaFree(qflux);
    cudaFree(hflux);
    cudaFree(zData);
    cudaFree(hData);
    cudaFree(tData);
    cudaFree(maxDt);
    cudaFree(vents);
    cudaFree(neibor);
    cudaFree(limit_dev);
    cudaFree(param);
}

__device__
double calDensity(double temperature, ParaData* param){
    if(param->m_IsDStatic)
        return param->m_IsDStatic;
    else
        // A / ( B + (1672.15 - Temperature)*C )
        return param->m_DenParaA / (param->m_DenParaB + (1672.15 - temperature) * param->m_DenParaC);
}

__device__
double calCoupledHeat(double tmp){
    const double p1 = 3.48992826 * pow(10.0, -7);
    const double p2 = -0.000722366524;
    const double p3 = 0.60978315;
    const double p4 = -201.805106;
    double res = 0.0;
    if(tmp>=900)
        res += p1 * pow(tmp, 3) + p2 * pow(tmp, 2) + p3 * pow(tmp, 1) + p4;
    return res * 1000; // return with unit (W * m^-2)
}

__device__
double calViscosity(double temperature, ParaData* param){
    if(param->m_IsVStatic)
        return param->m_Vis;
    else {
        if(temperature > param->m_VisParaC)
            return pow(10.0, param->m_VisParaA + param->m_VisParaB / (temperature - param->m_VisParaC));
        else
            return pow(10.0, 50);
    }
        
}

__device__
double calSy(double temperature){
    return 0.01*(exp(0.08*(1407.15-temperature))-1);  //(Dragoni, 1989)
}

__device__
double calTempByQ(double cc, ParaData* param){
    return (param->m_DenParaB + 1672.15*param->m_DenParaC) \
            * cc / (param->m_DenParaA + param->m_DenParaC*cc);
}

__device__
void calDeltaHZ(double* hData, double* zData, double* deltaH, double* deltaZ, double* hforQtfa, bool* isOut, int ci, int ti){
    double zc = zData[ci];
    double zt = zData[ti];
    double hc = hData[ci];
    double ht = hData[ti];
    if( (zc + hc)>=(zt + ht) ){
        (*deltaH)=hc-ht;
        (*deltaZ)=zc-zt;
        (*isOut)=true;
    }
    else{
        (*deltaH)=ht-hc;
        (*deltaZ)=zt-zc;
        (*isOut)=false;
    }
    if(fabs(zc-zt)<=DZERO){
        if(hc>ht)
            (*hforQtfa) = ht;
        else
            (*hforQtfa) = hc;
    }
    else if(zc>zt){
        if(zc + hc > zt + ht)
            (*hforQtfa) = (zt+ht-zc)>0.0? (zt+ht-zc):0.0;
        else
            (*hforQtfa) = hc;
    }
    else{
        if(zc + hc > zt + ht)
            (*hforQtfa) = ht;
        else
            (*hforQtfa) = (zc+hc-zt)>0.0? (zc+hc-zt):0.0;
    }
    return;
}

__global__
void calEachFlux(double* qflux, double* hflux, double* hData, double* zData, double* tData, double* maxDt, int* neibor, int* limit, ParaData* param){
    int nRow = param->m_nRow;
    int nCol = param->m_nCol;
    double cellSize = param->m_CellSize;
    int totalNum = nRow*nCol;
    int cx = blockDim.x*blockIdx.x + threadIdx.x;
    int cy = blockDim.y*blockIdx.y + threadIdx.y;
    int ci = cy * nCol + cx;
    //check limits.
    //the cells around the upper and right bounder also requires calculation.
    if ( ci < 0 || ci > totalNum || cx<limit[0]-1 || cx>(limit[2]+1) || cy<(limit[3]-1) || cy>limit[1]+1)
        return;
    //std::cout<<"c: "<<ci<<" "<<cy<<" "<<cx<<std::endl;
    int dx[4] = {-1, -1, 0, +1};
    int dy[4] = {0, +1, +1, +1};
    //calculate q, Qm, Qtfa.
    int neibor_c = neibor[ci];
    for(int i=0; i<4; i++){
        //check neibor
        if((neibor_c & (1 << i))==0)
            continue;
        int tx = cx + dx[i];
        int ty = cy + dy[i];
        int ti = ty*nCol + tx;
        //check ranges.
        if(ti<0 || ti>=totalNum)
            continue;
        //std::cout<<"=="<<ti<<" "<<ty<<" "<<tx<<std::endl;
        bool isOut;
        double deltaH, deltaZ, HforQtfa;
        //calculate: direction_q, deltaH, deltaZ, HforQta.
        calDeltaHZ(hData, zData, &deltaH, &deltaZ, &HforQtfa, &isOut, ci, ti);
        //calculate: Hcr, q, maxDt.
        double hcr, a, q, tempDt;
        double temp, height, direct;
        //calculate: Qm, Qta.
        if(isOut){
            temp=tData[ci];
            height=hData[ci];
            direct=1.0;
        }else{
            temp=tData[ti];
            height=hData[ti];
            direct=-1.0;
        }
        //check overflow!
        //if(fabs(deltaZ-deltaH)<DZERO)
        //    continue;
        hcr = calSy(temp) * sqrt(pow(deltaZ, 2) + pow(cellSize, 2))\
                /(calDensity(temp, param) * param->m_GravAcc * (deltaZ + deltaH));
        // hcr = calSy(temp) * sqrt(pow(deltaZ, 2) + pow(cellSize, 2))\
                /(calDensity(temp, param) * param->m_GravAcc * (deltaZ - deltaH));
        a = height/hcr;
        if(a>1.0)
            q = calSy(temp) * pow(hcr, 2) * cellSize * (a*a*a-1.5*a*a+0.5) / (3.0 * calViscosity(temp, param));
        else
            q = 0.0;
        //cout<<"Temp:"<<temp<<" Sy:"<<this->calSy(temp)<<" Sqrt:"<<sqrt(pow(deltaZ, 2)+ pow(this->m_Para.m_CellSize, 2))<<\
        //        " Density:"<<this->calDensity(temp)<<" Viscosity:"<<this->calViscosity(temp)<<" Hcr:"<<hcr<<" q:"<<q<<\
        //        " deltaZ:"<<deltaZ<<" height:"<<height<<" a:"<<a<<endl;
        //update qflux array.
        atomicdoubleAdd(&qflux[ci], -direct*q);
        atomicdoubleAdd(&qflux[ti], direct*q);
        //Caculate Qta, the heat transferred between cells.
        double Qtfa = param->m_Lambda * HforQtfa * cellSize * (tData[ci]-tData[ti]);
        atomicdoubleAdd(&hflux[ci], -Qtfa);
        atomicdoubleAdd(&hflux[ti], Qtfa);

        //update Qm, the heat carried by the lava flux.
        double Qm = q * calDensity(temp, param) * param->m_Cv * temp;
        atomicdoubleAdd(&hflux[ci], -direct * Qm);
        atomicdoubleAdd(&hflux[ti], direct * Qm);
        if(q<=DZERO)
            continue;
        tempDt = height * pow(cellSize, 2) / q;
        atomicMin(&maxDt[ci], tempDt);
    }

    //calculate the gradient factor.
    double factor = 1.0;
    if(param->m_AreaCorrection) {
        double grad_xz = 0.0, grad_yz = 0.0;
        int c_x1 = cy * nCol + cx - 1, c_x2 = cy * nCol + cx + 1;
        int c_y1 = (cy - 1) * nCol, c_y2 = (cy + 1) * nCol;
        if (c_x1 >= 0 && c_x1 < totalNum && c_y1 >= 0 && c_y1 < totalNum \
 && c_x2 >= 0 && c_x2 < totalNum && c_y2 >= 0 && c_y2 < totalNum) {
            grad_xz = ((hData[c_x1] + zData[c_x1]) - (hData[c_x2] + zData[c_x2])) / (2. * cellSize);
            grad_yz = ((hData[c_y1] + zData[c_y1]) - (hData[c_y2] + zData[c_y2])) / (2. * cellSize);
        }
        factor *= (double) sqrt(pow(grad_xz, 2) + pow(grad_yz, 2) + 1);
    }
    
    if(hData[ci] >= 0.0) {
        if(param->m_IsCoupled)
            //the coupled heat computation method.
            atomicdoubleAdd(&hflux[ci], -calCoupledHeat(tData[ci]) * pow(cellSize, 2) * factor);
        else
            //qtr
            atomicdoubleAdd(&hflux[ci], -(param->m_Emi * pow(cellSize, 2) * param->m_SBConstant * (pow(tData[ci], 4)-pow(param->m_AtmTemperature, 4))) * factor);
        //Qtc with the ground
        atomicdoubleAdd(&hflux[ci], -(param->m_AtmDen * param->m_AtmSpecificHeat * param->m_Gamma *\
            cbrt((param->m_GravAcc * param->m_Aa * param->m_Ka * param->m_Ka)/param->m_Va) *\
            cbrt(pow((tData[ci] - param->m_AtmTemperature), 4))*cellSize*cellSize*factor));
    }
    return;
}

__global__
void updateData(double* qflux, double* hflux, double* zData, double* hData, double* tData, int* vents, int* limit, int* changed, ParaData* param, double effusion, double timeStep, double lavaTemp){
    int nRow = param->m_nRow;
    int nCol = param->m_nCol;
    double cellSize = param->m_CellSize;
    
    int cx = blockDim.x*blockIdx.x + threadIdx.x;
    int cy = blockDim.y*blockIdx.y + threadIdx.y;
    int ci = cy*nCol + cx;

    //skip the inactive cells
    if( !(cx>=0 && cx>=(limit[0]-1) && cx<=(limit[2]+1) && cx<nCol && cy>=0 && cy>=(limit[3]-1) && cy<=(limit[1]+1) && cy<nRow) )
        return;

    double cv = param->m_Cv;

    if(vents[ci] != -1){
        qflux[ci] += effusion;
        //Q=V*Cv*T*Dense
        double temp = lavaTemp;
        hflux[ci] += effusion * cv * temp * calDensity(temp, param);
    }
    if(qflux[ci]>=-DZERO && qflux[ci]<=DZERO)
        return;

    double dense = calDensity(tData[ci], param);
    double A = pow(cellSize, 2);
    double Qt = dense * cv * hData[ci] * tData[ci]*A;

    hData[ci] += qflux[ci] * timeStep / A;

    Qt += hflux[ci] * timeStep; //tranditional computation method
    Qt = Qt>=0.0? Qt:0.0;
    tData[ci] = calTempByQ(Qt / (cv * hData[ci] * A), param); // Original Bad Method: tData[ci] = Qt / (dense * cv * hData[ci] * A);

    if(param->m_IsTransform){
        if(tData[ci] <= param->m_TransformTemp){
            zData[ci] += hData[ci];
            hData[ci] = 0.0;
            tData[ci] = param->m_AtmTemperature>=tData[ci]? param->m_AtmTemperature:tData[ci];
        }
    }
    //check limit changing.
    if(cy==(limit[3]-1))
        changed[3]=true;
    else if(cy==(limit[1]+1))
        changed[1]=true;
    if(cx==(limit[0]-1))
        changed[0]=true;
    else if(cx==(limit[2]+1))
        changed[2]=true;
    return;
}

__global__
void getTimeStep(double* maxDt, int nRow, int nCol){
    int id = threadIdx.x;
    int idx = threadIdx.x * nCol;
    for(int i=0; i<nCol; ++i)
        maxDt[idx] = maxDt[idx]<=maxDt[idx+i]? maxDt[idx]:maxDt[idx+i];
    __syncthreads();
    for(int s=nRow/2; s>0; s>>=1){
        if(id<s){
            maxDt[idx] = maxDt[idx]<=maxDt[(id+s)*nCol]? maxDt[idx]:maxDt[(id+s)*nCol];
        }
        __syncthreads();
    }
    return;
}
