#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <cmath>
#include <iostream>
#include <algorithm> 
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

typedef std::string string;

struct VolcanoData{
    int cnt;
    double* tList = NULL;
    double* vList = NULL;
    double temperature;
    void reset(VolcanoData& cpy){
        if(this->tList)
            delete[] this->tList;
        if(this->vList)
            delete[] this->vList;
        this->cnt=cpy.cnt;
        this->temperature=cpy.temperature;
        this->tList=new double[this->cnt];
        this->vList=new double[this->cnt];
        memcpy(this->tList, cpy.tList, sizeof(double)*this->cnt);
        memcpy(this->vList, cpy.vList, sizeof(double)*this->cnt);
    };
    VolcanoData(){
        this->cnt=0;
    }
    VolcanoData(int cnt){
        this->cnt=cnt;
        this->tList=new double[this->cnt];
        this->vList=new double[this->cnt];
    }
    VolcanoData(VolcanoData& cpy){
        VolcanoData(cpy.cnt);
        this->temperature=cpy.temperature;
        memcpy(this->tList, cpy.tList, sizeof(double)*this->cnt);
        memcpy(this->vList, cpy.vList, sizeof(double)*this->cnt);
    }
    ~VolcanoData(){
        this->cnt=0;
        if(tList)
            delete[] tList;
        if(vList)
            delete[] vList;
    }
};

struct ParaData
{
    // lava
	bool m_IsVStatic;
	double m_Vis;
	double m_VisParaA;
	double m_VisParaB;
	double m_VisParaC;

	bool m_IsDStatic;
	double m_Den;
	double m_DenParaA;
	double m_DenParaB;
	double m_DenParaC;

    bool m_IsTransform;
	double m_TransformTemp;

    // envir
	const double m_SBConstant = 5.670367e-8;
    double m_Emi;
    double m_Cv;
    double m_Aa;
	double m_Ka;
	double m_Va;
	double m_Gamma;
	double m_Lambda;
    double m_AtmDen;
	double m_AtmSpecificHeat;
    double m_AtmTemperature;
    double m_GravAcc;
    bool m_IsCoupled;

    // ctl
    double m_RScale;
    double m_CtlDelta;
    double m_TimeStepMax;
    double m_SumTime;
    double m_OutputTime;
    int m_nRow;
    int m_nCol;
    double m_CellSize;
    bool m_AreaCorrection;

};

class Lavasim{
    public:
        Lavasim(){};

        string m_OutputFolder;

        ParaData m_Para;

        double* m_HList = NULL;
        double* m_ZList = NULL;
        double* m_TList = NULL;
        int* m_Vent = NULL;
        int* m_Neibor = NULL;
        int* limits = NULL;

        VolcanoData m_VolcanoData;
        //Init
        void initGrid(int row, int col, int cellSize){
            this->m_Para.m_nRow = row;
            this->m_Para.m_nCol = col;
            this->m_Para.m_CellSize = cellSize;
        };
        void initBorder(){
            this->limits = new int[4];
            // 0-left, 1-down, 2-right, 3-up.
            this->limits[0] = this->m_Para.m_nCol; this->limits[1]=0; this->limits[2]= 0; this->limits[3]=this->m_Para.m_nRow;
                //init bounder.
            for(int i = 0; i < this->m_Para.m_nRow * this->m_Para.m_nCol; ++i){
                if(this->m_Vent[i]!=-1){
                    int ir = i / this->m_Para.m_nCol;
                    int ic = i - ir * this->m_Para.m_nCol;
                    limits[0] = std::min<int>(limits[0], ic);
                    limits[1] = std::max<int>(limits[1], ir);
                    limits[2] = std::max<int>(limits[2], ic);
                    limits[3] = std::min<int>(limits[3], ir);
                }
            }
        }
        void setOutputFolder(string folder){ m_OutputFolder = folder; };
        void setRScale(double r){ m_Para.m_RScale = r; };
        void setCtlDelta(double d){ m_Para.m_CtlDelta = d; };
        void setTimeStepMax(double t){ m_Para.m_TimeStepMax = t; };
        void setSumTime(double sumTime){ m_Para.m_SumTime = sumTime; };
        void setOutputTime(double outputTime){ m_Para.m_OutputTime = outputTime; };
        void setGraveAcc(double gravAcc = 9.8) { m_Para.m_GravAcc = gravAcc; };
        void setAtmTemp(double atmTemp = 300) { m_Para.m_AtmTemperature = atmTemp; };
        void setViscosityPara(double vis) { m_Para.m_Vis = vis; m_Para.m_IsVStatic = true; };
        void setViscosityPara(double A, double B, double C) {
            m_Para.m_VisParaA = A;
            m_Para.m_VisParaB = B;
            m_Para.m_VisParaC = C;
            m_Para.m_IsVStatic = false;
        };
        void setDensityPara(double den) { m_Para.m_Den = den; m_Para.m_IsDStatic = true; };
        void setDensityPara(double A, double B, double C) {
            m_Para.m_DenParaA = A;
            m_Para.m_DenParaB = B;
            m_Para.m_DenParaC = C;
            m_Para.m_IsDStatic = false;
        };
        void setEmissivity(double emi){m_Para.m_Emi = emi;};
        void setCv(double cv){this->m_Para.m_Cv=cv; };
        void setHeatMode(bool bSet) { m_Para.m_IsCoupled = bSet; };
        void setAreaCorrection(bool b) { m_Para.m_AreaCorrection = b; };
        void setQcPara(double den, double c, double a, double k, double v, double r) {
            m_Para.m_AtmDen = den;
            m_Para.m_AtmSpecificHeat = c;
            m_Para.m_Aa = a;
            m_Para.m_Ka = k;
            m_Para.m_Va = v;
            m_Para.m_Gamma = r;
        };
        void setQtfPara(double l) { m_Para.m_Lambda = l; };
        void setTransform(double temp, bool bSet) { m_Para.m_TransformTemp = temp; m_Para.m_IsTransform = bSet; };
        bool setVolcanoData(VolcanoData& vcData, int* vcIdx){
            this->m_VolcanoData.reset(vcData);
            if(this->m_Vent)
                delete[] this->m_Vent;
            this->m_Vent = new int[this->m_Para.m_nRow*this->m_Para.m_nCol];
            memcpy(this->m_Vent, vcIdx, sizeof(int)*this->m_Para.m_nRow*this->m_Para.m_nCol);
            return true;
        };
        bool setHData(){
            setHData(0.0);
            return true;
        };
        bool setHData(double h){
            this->m_HList = new double[this->m_Para.m_nRow*this->m_Para.m_nCol];
            for(int i=0;i<this->m_Para.m_nRow*this->m_Para.m_nCol;i++)
                this->m_HList[i]=h;
            return true;
        };
        bool setHData(double* h){
            this->m_HList=new double[this->m_Para.m_nRow*this->m_Para.m_nCol];
            memcpy(this->m_HList, h, sizeof(double)*this->m_Para.m_nRow*this->m_Para.m_nCol);
            return true;
        };
        bool setZData(){
            setZData(0.0);
            return true;
        };
        bool setZData(double z){
            this->m_ZList = new double[this->m_Para.m_nRow*this->m_Para.m_nCol];
            for(int i=0;i<this->m_Para.m_nRow*this->m_Para.m_nCol;i++)
                this->m_ZList[i]=z;
            return true;
        };
        bool setZData(double* z){
            this->m_ZList = new double[this->m_Para.m_nRow*this->m_Para.m_nCol];
            memcpy(this->m_ZList, z, sizeof(double)*this->m_Para.m_nRow*this->m_Para.m_nCol);
            return true;
        };
        bool setTData(){
            setTData(m_Para.m_AtmTemperature);
            return true;
        };
        bool setTData(double t){
            this->m_TList = new double[this->m_Para.m_nRow*this->m_Para.m_nCol];
            for(int i=0;i<this->m_Para.m_nRow*this->m_Para.m_nCol;i++)
                this->m_TList[i]=t;
            return true;
        };
        bool setTData(double* t){
            this->m_TList=new double[this->m_Para.m_nRow*this->m_Para.m_nCol];
            memcpy(this->m_TList, t, sizeof(double)*this->m_Para.m_nRow*this->m_Para.m_nCol);
            return true;
        };
        bool setNeibor(int* neibor){
            this->m_Neibor=new int[this->m_Para.m_nRow*this->m_Para.m_nCol];
            memcpy(this->m_Neibor, neibor, sizeof(int)*this->m_Para.m_nRow*this->m_Para.m_nCol);
            return true;
        };
        ~Lavasim(){
            if(m_HList)
                delete[] m_HList;
            if(m_ZList)
                delete[] m_ZList;
            if(m_TList)
                delete[] m_TList;
            if(m_Vent)
                delete[] m_Vent;
            if(m_Neibor)
                delete[] m_Neibor;
        };
        bool setNeibor();
        void outputParam();
        void checkpoint(double timeStep, int sumTime, int sumCount);

        __host__ void doCal();
        __host__ void initCalData(double** qflux, double** hflux, double** zData, double** hData, double** tData, double** maxQ, int** vents, int** neibor, ParaData** param);
        __host__ double getVolcano(double nowTime);
};

__global__ void getTimeStep(double* maxDt, int nRow, int nCol);
__global__ void updateData(double* qflux, double* hflux, double* zData, double* hData, double* tData, int* vents, int* limit, int* changed, ParaData* param, double effusion, double timeStep, double lavaTemp);
__global__ void calEachFlux(double* qflux, double* hflux, double* hData, double* zData, double* tData, double* maxDt, int* neibor, int* limit, ParaData* param);
__device__ void calDeltaHZ(double* hData, double* zData, double* deltaH, double* deltaZ, double* hforQtfa, bool* isOut, int ci, int ti);
__device__ double calSy(double temperature);
__device__ double calViscosity(double temperature);
__device__ double calDensity(double temperature);
__device__ double calTempByQ(double cc, ParaData* param);
__global__ void checkBoundary(int* limit, int* changed, int nRow, int nCol);

__device__ static double atomicMin(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*) address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = ::atomicCAS(address_as_ull, assumed,
            __double_as_longlong(::fminf(val, __longlong_as_double(assumed))));
    } while (assumed != old);
    return __longlong_as_double(old);
}
__device__ inline void atomicdoubleAdd(double *address, double val)
{
       unsigned long long int ull_val = __double_as_longlong(val);
       unsigned long long int tmp0 = 0;
       unsigned long long int tmp1;
       while( (tmp1 = atomicCAS((unsigned long long int *)address, tmp0, ull_val)) != tmp0)
       {
               tmp0 = tmp1;
               ull_val = __double_as_longlong(val + __longlong_as_double(tmp1));
       }
}
