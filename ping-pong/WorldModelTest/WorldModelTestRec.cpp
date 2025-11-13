#include <iostream>
#include "WorldModelTest.hpp"

using namespace std;

WorldModelTest *pwmtG = NULL;

const int aindReceptorSection[] = {0, 30, 60, 69, 78, 108};

STRING strtrans(const VECTOR<unsigned> &vfl_InputSignal)
{
    const char *apchReceptorSection[] = {"x", "y", "vx", "vy", "ry"};
    int i = 0;
    int ReceptorSection = 0;
    std::stringstream ss;
    auto pfl = &vfl_InputSignal.front();
    BitMaskAccess bma;
    while (ReceptorSection < sizeof(apchReceptorSection) / sizeof(apchReceptorSection[0])) {
        while (!(pfl & bma)) {
            ++bma;
            ++i;
            if (i == aindReceptorSection[sizeof(aindReceptorSection) / sizeof(aindReceptorSection[0]) - 1] * 4)
                return ss.str();
        }
        int indReceptor = i / 4;
        int TimePeriod = i % 4;
        ss << apchReceptorSection[ReceptorSection];
        int offset_in_section = indReceptor - aindReceptorSection[ReceptorSection];
        int section_size = aindReceptorSection[ReceptorSection + 1] - aindReceptorSection[ReceptorSection];
        ss << (offset_in_section < section_size / 2 || (section_size & 1) ? -section_size / 2 + offset_in_section : -section_size / 2 + offset_in_section + 1) << '/' << TimePeriod;
        ++i;
        ++bma;
        ++ReceptorSection;
    }
    return ss.str();
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
    ofsout << strtrans(vfl_InputSignal) << endl;
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
