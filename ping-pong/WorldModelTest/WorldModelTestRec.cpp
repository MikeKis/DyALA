#include <iostream>
#include "WorldModelTest.hpp"

using namespace std;

WorldModelTest *pwmtG = NULL;

const int aindReceptorSection[] = {0, 30, 60, 69, 78, 108};
const int nPeriods = 4;

WorldModelTest::WorldModelTest(const pugi::xml_node& xn, int nRec)
{
    SetNReceptors(aindReceptorSection[5] * nPeriods);
    vfl_InputSignal.resize(AfferentSpikeBufferSizeDW(GetNReceptors()));
    pflMy = &vfl_InputSignal.front();
    // ifstream ifsState("ping_pong_state.csv");
    // while (!ifsState.eof()) {
    //     string s;
    //     getline(ifsState, s);
    //     vector<float> vr_(4);
    //     sscanf(s.c_str(), "%*f,%f,%f,%f,%f", &vr_[0], &vr_[1], &vr_[2], &vr_[3]);
    //     vvr_States.push_back(vr_);
    // }
    ifstream ifsCompleteState("ping_pong_state_complete.csv");
    while (!ifsCompleteState.eof()) {
        string s;
        getline(ifsCompleteState, s);
        PAIR<int,VECTOR<PAIR<int, int> > > ptactvp_{0, VECTOR<PAIR<int, int> >(5)};
        sscanf(
               s.c_str(), 
               "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d",
               &ptactvp_.first, 
               &ptactvp_.second[0].first, 
               &ptactvp_.second[0].second, 
               &ptactvp_.second[1].first, 
               &ptactvp_.second[1].second, 
               &ptactvp_.second[2].first, 
               &ptactvp_.second[2].second, 
               &ptactvp_.second[3].first, 
               &ptactvp_.second[3].second,
               &ptactvp_.second[4].first, 
               &ptactvp_.second[4].second
              );
        vptactvp_CompleteStates.push_back(ptactvp_);
    }



    // period = 0;
    // depth = 0;
    // STRING strSpikesInFile;
    // GetAllParameters(xn, "period", &period, "depth", &depth, "spikes_file", &strSpikesInFile);
    // ifsspikein.open(strSpikesInFile.c_str());
    // ofsout.open("WorldModelTest.csv");
    const char *apchReceptorSection[] = { "x", "y", "vx", "vy", "ry" };
    int i = 0;
    int ReceptorSection = 0;
    vstr_MyMeanings.resize(GetNReceptors());
    while (ReceptorSection < sizeof(apchReceptorSection) / sizeof(apchReceptorSection[0])) {
        int section_size = aindReceptorSection[ReceptorSection + 1] - aindReceptorSection[ReceptorSection];
        for (int offset_in_section = 0; offset_in_section < section_size; ++offset_in_section)
            for (int TimePeriod = 0; TimePeriod < 4; ++TimePeriod) {
               std::stringstream ss;
               ss << apchReceptorSection[ReceptorSection];
               ss << (ReceptorSection != 2 && ReceptorSection != 3 ? offset_in_section : offset_in_section - 4) << '/' << TimePeriod;
               vstr_MyMeanings[i++] = ss.str();
           }
        ++ReceptorSection;
    }
}

RECEPTORS_SET_PARAMETERS(pchMyReceptorSectionName, nReceptors, xn)
{
	pwmtG = new WorldModelTest(xn, nReceptors);
	return pwmtG;
}

DYNAMIC_LIBRARY_ENTRY_POINT IReceptors *LoadStatus(Serializer &ser)
{
	return nullptr;
}

bool WorldModelTest::bGenerateSignals(unsigned* pfl, int bitoffset)
{
    static int ntact = 0;
    if (bReset) {
        int indstart = 0;
        cin >> indstart;
        vp_CurrentCompleteState = vptactvp_CompleteStates[indstart].second;
        bReset = false;
    }
    if (!vp_CurrentCompleteState[0].first)
        return false;
    fill(vfl_InputSignal.begin(), vfl_InputSignal.end(), 0);
    FORI(vp_CurrentCompleteState.size())
        pflMy |= BitMaskAccess((aindReceptorSection[_i] + vp_CurrentCompleteState[_i].first) * nPeriods + vp_CurrentCompleteState[_i].second);
    copy(vfl_InputSignal.begin(), vfl_InputSignal.end(), pfl);
    return ++ntact < 3000;
}

void WorldModelTest::change_to(int indrec)
{
    int ReceptorSection = (int)(upper_bound(aindReceptorSection, aindReceptorSection + sizeof(aindReceptorSection) / sizeof(aindReceptorSection[0]), indrec / nPeriods) - aindReceptorSection) - 1;
    if (ReceptorSection <= 3) {   // excluding Yracket
        int offset = indrec - aindReceptorSection[ReceptorSection] * nPeriods;
        pair<int, int> p_NewState(offset / nPeriods, offset % nPeriods);
        if (vp_CurrentCompleteState[ReceptorSection] == p_NewState)
            throw std::runtime_error("the same state predicted");
        vp_CurrentCompleteState[ReceptorSection] = p_NewState;
        vp_CurrentCompleteState[4] = vptactvp_CompleteStates[ran(vptactvp_CompleteStates.size())].second[4];   // Yracket every time is set randomly.
    }
}
