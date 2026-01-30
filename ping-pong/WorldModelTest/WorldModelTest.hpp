#pragma once

#include <sg/sg.h>
#include <NetworkConfigurator.h>

#include "../ping-pong-environment.h"

class DYNAMIC_LIBRARY_EXPORTED_CLASS WorldModelTest: public IReceptors
{
    int              period;
    int              depth;
    std::ifstream    ifsspikein;
    std::ofstream    ofsout;
    int              tactin = 0;
    int              period_phase = 0;
    VECTOR<unsigned> vfl_InputSignal;
    VECTOR<STRING>   vstr_MyMeanings;
    STRING strtrans() const
    {
        int i;
        int ReceptorSection = 0;
        auto pfl = &vfl_InputSignal.front();
        BitMaskAccess bma;
        STRING strret;
        FOR_(i, GetNReceptors()) 
            if (pfl & bma++)
                strret += vstr_MyMeanings[i] + ' ';
        return strret;
    }
public:
    WorldModelTest(const pugi::xml_node& xn, int nRec);
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override;
    virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const {vstr_Meanings = vstr_MyMeanings;}
    bool operator[](int ind) const {return &vfl_InputSignal.front() & BitMaskAccess(ind);}
    void change_to(int indrec);
};
