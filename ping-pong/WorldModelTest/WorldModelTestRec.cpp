#include <iostream>
#include "WorldModelTest.hpp"

using namespace std;

WorldModelTest *pwmtG = NULL;

const int aindReceptorSection[] = {0, 30, 60, 69, 78, 108};

WorldModelTest::WorldModelTest(const pugi::xml_node& xn, int nRec)
{
    SetNReceptors(nRec);
    period = 0;
    depth = 0;
    STRING strSpikesInFile;
    GetAllParameters(xn, "period", &period, "depth", &depth, "spikes_file", &strSpikesInFile);
    ifsspikein.open(strSpikesInFile.c_str());
    ofsout.open("WorldModelTest.csv");
    vfl_InputSignal.resize(AfferentSpikeBufferSizeDW(GetNReceptors()));
    const char* apchReceptorSection[] = { "x", "y", "vx", "vy", "ry" };
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
    string s;
    if (!period_phase)
        do {
            if (!getline(ifsspikein, s).good()) {
                cout << "unexpected eof in spike file\n";
                exit(-50000);
            }
        } while (GetSpikesfromText(s, vfl_InputSignal) != 5);
    else if (period_phase > depth - 30)
        fill(vfl_InputSignal.begin(), vfl_InputSignal.end(), 0);
    copy(vfl_InputSignal.begin(), vfl_InputSignal.end(), pfl);
    ofsout << strtrans() << endl;
    ofsout.flush();
    if (++period_phase == depth) {
        period_phase = 0;
        FORI(period)
            if (!getline(ifsspikein, s).good()) {
                cout << "unexpected eof in spike file\n";
                exit(-50000);
            }
    }
    return true;
}

void WorldModelTest::change_to(int indrec)
{
    int ReceptorSection = (int)(upper_bound(aindReceptorSection, aindReceptorSection + sizeof(aindReceptorSection) / sizeof(aindReceptorSection[0]), indrec / 4) - aindReceptorSection) - 1;
    BitMaskAccess bma(aindReceptorSection[ReceptorSection] * 4);
    auto pfl = &vfl_InputSignal.front();
    FORI((aindReceptorSection[ReceptorSection + 1] - aindReceptorSection[ReceptorSection]) * 4)
        pfl != bma++;
    pfl |= BitMaskAccess(indrec);
}
