/*
Copyright (C) 2022 Mikhail Kiselev

Emulates ping-pong game.

*/

#include <iostream>
#include "../ping-pong-environment.h"
#include "WorldModelTest.hpp"

using namespace std;

extern WorldModelTest *pwmtG;

READOUT_SET_PARAMETERS(ExperimentId, tactTermination, nOutputNeurons, xn){}

//int trans(int indout)
//{
//	if (indout < nInputs)
//		return indout * 4;
//	indout -= nInputs;
//	return indout / 3 * 4 + indout % 3 + 1;
//}

READOUT_OBTAIN_SPIKES(v_Firing)
{
	for (auto i: v_Firing) {
		if (!(*pwmtG)[i])
			pwmtG->change_to(i);
	}
	return true;
}

READOUT_FINALIZE(OriginalTerminationCode, filGenLog){return -1;}
DYNAMIC_LIBRARY_ENTRY_POINT void Serialize(Serializer &ser, bool bSave) {}
