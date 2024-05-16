#include "Initializer.h"

Initializer::Initializer()
    : m_InputFolder("Input/")
    , m_OutputFolder("Output/") {
}

bool Initializer::setInputFolder(std::string folder) {
    this->m_InputFolder = folder;
    return true;
}

bool Initializer::setOutputFolder(std::string folder) {
    this->m_OutputFolder = folder;
    return true;
}

bool Initializer::loadCtlData(Lavasim& lsm){
    std::string log = this->m_OutputFolder + "lavasim.log";
    std::ofstream ofs;
    ofs.open(log.data(), std::ios::out);

    std::string filename = this->m_InputFolder + "ctl";
    ofs << "Loading " << filename << "..." << std::endl;
    std::ifstream ifs;
    ifs.open(filename.data(), std::ios::in);
    if(!ifs.is_open()) {
        ofs << "Error with ctl file!\n";
        return false;
    }

    char checker[30] = "ctl";
    char des[100];
    ifs >> des;
    if(strcmp(des, checker) != 0) {
        ofs << "Error with ctl checker!\n";
        return false;
    }

    int cnt;
    ifs >> cnt;
    std::map<string, double> m;
    for(int i = 0; i < cnt; ++i) {
        string key;
        double value;
        ifs >> key >> value;
        m[key] = value;
    }

    if(m.count("RScale") * m.count("TimeStepMax")\
        * m.count("SumTime") * m.count("OutputTime")\
        * m.count("CtlDelta") * m.count("CellSize")\
        * m.count("nRow") * m.count("nCol") == 0) {
        ofs << "Error with ctl data!\n";
        return false;
    }

    lsm.setRScale(m["RScale"]);
    lsm.setCtlDelta(m["CtlDelta"]);
    lsm.setTimeStepMax(m["TimeStepMax"]);

    lsm.setSumTime(m["SumTime"]);
    lsm.setOutputTime(m["OutputTime"]);
    m_nRow = m["nRow"];
    m_nCol = m["nCol"];
    lsm.initGrid(m["nRow"], m["nCol"], m["CellSize"]);
    lsm.setAreaCorrection(m["AreaCorrection"]);

    ifs.close();
    ofs.close();

    return true;
}

bool Initializer::loadVolcanoData(Lavasim& lsm) {
    std::string log = this->m_OutputFolder + "lavasim.log";
    std::ofstream ofs;
    ofs.open(log.data(), std::ios::app);

    std::string filename = this->m_InputFolder + "volcano";
    ofs << "Loading " << filename << "..." << std::endl;
    std::ifstream ifs;
    ifs.open(filename.data(), std::ios::in);
    if(!ifs.is_open()) {
        ofs << "Error with volcano file!\n";
        return false;
    }

    char checker[4][30] = {"volcano", "temperature", "rate", "location"};
    char des[100];
    ifs >> des;
    if(strcmp(des, checker[0]) != 0) {
        ofs << "Error with volcano checker!\n";
        return false;
    }
    int row, col;
    ifs >> row >> col;
    if(row != this->m_nRow || col != this->m_nCol) {
        ofs << "Error with volcano data!\n";
        return false;
    }

    ifs >> des;
    if(strcmp(des, checker[1]) != 0) {
        ofs << "Error with temperature checker!\n";
        return false;
    }
    double temperature;
    ifs >> temperature;

    ifs >> des;
    if(strcmp(des, checker[2]) != 0) {
        ofs << "Error with rate checker!\n";
        return false;
    }
    int cnt;
    ifs >> cnt;
    VolcanoData vcData(cnt);
    vcData.temperature = temperature;
    for(int i = 0; i < cnt; ++i) {
        ifs >> vcData.tList[i] >> vcData.vList[i];
        if(i > 0 && vcData.tList[i] < vcData.tList[i - 1])
            ofs << "Error with effusion data!\n";
    }
    
    ifs >> des;
    if(strcmp(des, checker[3]) != 0) {
        ofs << "Error with location checker!\n";
        return false;
    }
    int* idxMap = new int[row * col];
    int num;
    ifs >> num;
    memset(idxMap, -1, sizeof(int) * row * col);
    for(int i = 0; i < num; ++i){
        int ridx, cidx;
        ifs >> ridx >> cidx;
        idxMap[ridx * col + cidx] = i;
    }
    lsm.setVolcanoData(vcData, idxMap);
    delete[] idxMap;

    ifs.close();
    ofs.close();
    return true;
}

bool Initializer::loadEnvirData(Lavasim& lsm){
    std::string log = this->m_OutputFolder + "lavasim.log";
    std::ofstream ofs;
    ofs.open(log.data(), std::ios::app);

    std::string filename = this->m_InputFolder + "envir";
    ofs << "Loading " << filename << "..." << std::endl;
    std::ifstream ifs;
    ifs.open(filename.data(), std::ios::in);
    if(!ifs.is_open()) {
        ofs << "Error with envir file!\n";
        return false;
    }

    char checker[30] = "envir";
    char des[100];
    ifs >> des;
    if(strcmp(des, checker) != 0) {
        ofs << "Error with envir checker!\n";
        return false;
    }

    int cnt;
    ifs >> cnt;
    std::map<string, double> m;
    for(int i = 0; i < cnt; ++i) {
        string key;
        double value;
        ifs >> key >> value;
        m[key] = value;
    }

    if(m.count("AtmDen") * m.count("AtmSpecificHeat")\
        * m.count("Aa") * m.count("Ka") * m.count("Va")\
        * m.count("Gamma") * m.count("Lambda")\
        * m.count("Cv") * m.count("Emi")\
        * m.count("GravAcc") == 0) {
        ofs << "Error with envir data!\n";
        return false;
    }
    lsm.setQcPara(m["AtmDen"], m["AtmSpecificHeat"],\
        m["Aa"], m["Ka"], m["Va"], m["Gamma"]);
    lsm.setAtmTemp(m["AtmTemperature"]);
    lsm.setCv(m["Cv"]);
    lsm.setEmissivity(m["Emi"]);
    lsm.setQtfPara(m["Lambda"]);
    lsm.setGraveAcc(m["GravAcc"]);
    lsm.setHeatMode(m["IsCoupled"]);

    ifs.close();
    ofs.close();

    return true;

}

bool Initializer::loadLavaData(Lavasim& lsm){
    std::string log = this->m_OutputFolder + "lavasim.log";
    std::ofstream ofs;
    ofs.open(log.data(), std::ios::app);

    std::string filename = this->m_InputFolder + "lava";
    ofs << "Loading " << filename << "..." << std::endl;
    std::ifstream ifs;
    ifs.open(filename.data(), std::ios::in);
    if(!ifs.is_open()) {
        ofs << "Error with lava file!\n";
        return false;
    }

    char checker[30] = "lava";
    char des[100];
    ifs >> des;
    if(strcmp(des, checker) != 0) {
        ofs << "Error with lava checker!\n";
        return false;
    }

    int cnt;
    ifs >> cnt;
    std::map<string, double> m; // Create a dictionary
    for(int i = 0; i < cnt; ++i) {
        string key;
        double value;
        ifs >> key >> value;
        m[key] = value;
    }

    if((m["IsVStatic"] == 1 && m.count("StaticVis") == 0)\
        || (m["IsDStatic"] == 1 && m.count("StaticDen") == 0)\
        || m.count("VisParaA") * m.count("VisParaB") * m.count("VisParaC") == 0\
        || m.count("DenParaA") * m.count("DenParaB") * m.count("DenParaC") == 0) {
        ofs << "Error with lava data!\n";
        return false;
    }

    lsm.setTransform(m["TransTemperature"], m["IsCoolTrans"]);
    
    if(m["IsVStatic"] == 1)
        lsm.setViscosityPara(m["StaticVis"]);
    else
        lsm.setViscosityPara(m["VisParaA"], m["VisParaB"], m["VisParaC"]);
    
    if(m["IsDStatic"] == 1)
        lsm.setDensityPara(m["StaticDen"]);
    else
        lsm.setDensityPara(m["DenParaA"], m["DenParaB"], m["DenParaC"]);

    ifs.close();
    ofs.close();
    
    return true;
}

