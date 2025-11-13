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
public:
	WorldModelTest(const pugi::xml_node &xn, int nRec)
	{
        SetNReceptors(nRec);
        period = 0;
        depth = 0;
        STRING strSpikesInFile;
        GetAllParameters(xn, "period", &period, "depth", &depth, "spikes_file", &strSpikesInFile);
        ifsspikein.open(strSpikesInFile.c_str());
        ofsout.open("WorldModelTest.csv");
        vfl_InputSignal.resize(AfferentSpikeBufferSizeDW(GetNReceptors()));
    }
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override;
    bool operator[](int ind) const {return &vfl_InputSignal.front() & BitMaskAccess(ind);}
    void change_to(int indrec);
};
