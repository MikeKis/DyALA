#pragma once

#include <sg/sg.h>
#include <NetworkConfigurator.h>

#include "../ping-pong-environment.h"

class DYNAMIC_LIBRARY_EXPORTED_CLASS WorldModelTest: public IReceptors
{
    VECTOR<PAIR<int,VECTOR<PAIR<int, int> > > > vptactvp_CompleteStates;
    bool                                        bReset = true;
    Rand                                        ran;
    VECTOR<unsigned>                            vfl_InputSignal;
    unsigned                                   *pflMy;
    VECTOR<STRING>                              vstr_MyMeanings;
public:
    WorldModelTest(const pugi::xml_node &xn, int nRec);
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override;
    virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const {vstr_Meanings = vstr_MyMeanings;}
    void change_to(int indrec);
    VECTOR<PAIR<int, int> > vp_CurrentCompleteState;
    int                     ntact = 0;
    int                     nhops = 0;
};
