#include <iostream>
#include "Initializer.h"
#include <string>
#include <ctime>
//#include <mgl2/mgl_cf.h>

void wrap(Lavasim& ctl, Initializer& initer){
    if(!initer.loadCtlData(ctl))
        exit(0);
    ctl.setOutputFolder(initer.m_OutputFolder);
    if(!(initer.loadVolcanoData(ctl) && initer.loadEnvirData(ctl) && initer.loadLavaData(ctl)))
        exit(0);
    //init Border
    ctl.initBorder();
    //temperature
    double* data = new double[initer.m_nRow*initer.m_nCol];
    if(initer.loadData<double>(data, "tdata"))
        ctl.setTData(data);
    else
        ctl.setTData();
    //thickness
    if(initer.loadData<double>(data, "hdata"))
        ctl.setHData(data);
    else
        ctl.setHData();
    //grid
    if(initer.loadData<double>(data, "zdata"))
        ctl.setZData(data);
    else
        ctl.setZData();
    delete[] data;
    //Neibor
    int* neibor = new int[initer.m_nRow*initer.m_nCol];
    if(initer.loadData<int>(neibor, "neibor"))
        ctl.setNeibor(neibor);
    else
        ctl.setNeibor();
    delete[] neibor;
    return;
}

int main(int argc, char* argv[]) {
    Initializer initer;
    Lavasim sim;
    if(argc >= 3) {
        initer.setInputFolder(argv[1]);
        initer.setOutputFolder(argv[2]);
    }
    wrap(sim, initer);
    //std::cerr << initer.m_nRow << initer.m_nCol << initer.m_CellSize << std::endl;
    //std::cerr << sim.m_Para.m_nRow << sim.m_Para.m_nCol << sim.m_Para.m_CellSize << std::endl;

    sim.doCal();
    return 0;
}
