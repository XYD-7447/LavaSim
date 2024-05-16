#include <fstream>
#include <algorithm>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include "LavasimCal.h"
//#include "visual.cpp"

typedef std::string string;

class Initializer {
    public:
        // ctl
        std::string m_InputFolder;
        std::string m_OutputFolder;
        int m_nRow;
        int m_nCol;

        Initializer();
        bool setInputFolder(std::string folder);
        bool setOutputFolder(std::string folder);
        bool loadCtlData(Lavasim& lsm);
        bool loadVolcanoData(Lavasim& lsm);
        bool loadEnvirData(Lavasim& lsm);
        bool loadLavaData(Lavasim& lsm);
        bool setPlanetEnvir(std::string envir);
        bool setLavaType(std::string type);
        
        template<typename T>
        bool loadData(T* ary, std::string name){
            string log = this->m_OutputFolder + "lavasim.log";
            std::ofstream ofs;
            ofs.open(log.data(), std::ios::app);

            std::string filename = this->m_InputFolder + name;
            ofs << "Loading " << filename << "..." << std::endl;
            std::ifstream ifs;
            ifs.open(filename.data(), std::ios::in);
            if(!ifs.is_open())
                return false;
            const char* checker = name.data();
            char des[100];
            ifs >> des;
            if(strcmp(des, checker)!=0)
                return false;
            int row, col;
            ifs >> row >> col;
            if(row != this->m_nRow || col != this->m_nCol)
                return false;
            for(int i = 0; i < row * col; ++i)
                ifs >> ary[i];
            ofs.close();
            return true;
        }
};
