#pragma once

#include <sg/sg.h>
#include <NetworkConfigurator.h>

#include "../ping-pong-environment.h"

class DYNAMIC_LIBRARY_EXPORTED_CLASS WorldModelTest: public IReceptors
{
//    VECTOR<VECTOR<float> >                      vvr_States;
    VECTOR<PAIR<int,VECTOR<PAIR<int, int> > > > vptactvp_CompleteStates;
    bool                                        bReset = true;
    Rand                                        ran;
//    int                                         tactStart;
    VECTOR<unsigned>                            vfl_InputSignal;
    unsigned                                   *pflMy;

    // int              period;
    // int              depth;
    // std::ifstream    ifsspikein;
    // std::ofstream    ofsout;
    // int              tactin = 0;
    // int              period_phase = 0;
    VECTOR<STRING>   vstr_MyMeanings;
    // STRING strtrans() const
    // {
    //     int i;
    //     int ReceptorSection = 0;
    //     auto pfl = &vfl_InputSignal.front();
    //     BitMaskAccess bma;
    //     STRING strret;
    //     FOR_(i, GetNReceptors())
    //         if (pfl & bma++)
    //             strret += vstr_MyMeanings[i] + ' ';
    //     return strret;
    // }
public:
    WorldModelTest(const pugi::xml_node& xn, int nRec);
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override;
    virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const {vstr_Meanings = vstr_MyMeanings;}
    void change_to(int indrec);
    VECTOR<PAIR<int, int> > vp_CurrentCompleteState;
};
