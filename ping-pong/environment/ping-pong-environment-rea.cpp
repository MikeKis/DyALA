/*
Copyright (C) 2022 Mikhail Kiselev

Emulates ping-pong game.

*/

#include "ClusterBayes.hpp"
#include "ping-pong-environment.h"
#include "EnvironmentState.hpp"

using namespace std;

extern int tactStart;
extern int NeuronTimeDepth;
extern float rStateFiringFrequency;
extern int LevelDuration;
extern int LevelNeuronPeriod;
extern int minnSignificantExtraSpikes;
extern int CurrentLevel;
extern int tactLevelFixed;
extern int ntact;
extern unique_ptr<ClusterBayes> cb;
extern EnvironmentState es;

extern int nRewards;
extern int nPunishments;
extern int nRewardsTot;
extern int nPunishmentsTot;

extern bool b_forVerifier_Reward;

extern int LastRegistration;
extern int PostRewardCounter;
extern int nRecognitions;
extern int nCorr;

int nNeuronsperAction;

PING_PONG_ENVIRONMENT_EXPORT void SetParametersOut(int ExperimentId, size_t tactTermination, unsigned nOutputNeurons, const pugi::xml_node &xn) 
{
	nNeuronsperAction = (int)nOutputNeurons / 2;
	tactStart = atoi_s(xn.child("start_time").child_value());

	NeuronTimeDepth = atoi_s(xn.child("NeuronTimeDepth").child_value());
	float rNSpikesperNeuronTime = atof_s(xn.child("rNSpikesperNeuronTime").child_value());
	rStateFiringFrequency = rNSpikesperNeuronTime / NeuronTimeDepth;
	LevelDuration = atoi_s(xn.child("LevelDuration").child_value());
	LevelNeuronPeriod = max(LevelDuration / NeuronTimeDepth, 1);
	minnSignificantExtraSpikes = atoi_s(xn.child("minnSignificantExtraSpikes").child_value());

}

const float rAction = 1.F / nSpatialZones;

PING_PONG_ENVIRONMENT_EXPORT bool ObtainOutputSpikes(const vector<int> &v_Firing, int nEquilibriumPeriods)
{
	static int NoMoveTacts = 0;
	int nCommandsDown = count_if(v_Firing.begin(), v_Firing.end(), bind2nd(less<int>(), nNeuronsperAction));
	auto r = rAction * ((int)v_Firing.size() - 2 * nCommandsDown);
	if (r) {
		CurrentLevel = cb->Predict();
		tactLevelFixed = ntact;
	}
	auto rsav = *es.prRacket;
	*es.prRacket += r;
	if (*es.prRacket > 0.5F - RACKET_SIZE / 2)
		*es.prRacket = 0.5F - RACKET_SIZE / 2;
	else if (*es.prRacket < -0.5F + RACKET_SIZE / 2)
		*es.prRacket = -0.5F + RACKET_SIZE / 2;

	if (*es.prRacket == rsav) {
		if (++NoMoveTacts == 100000)
			return false;
	} else NoMoveTacts = 0;

	if (ntact && !(ntact % 200000)) {
		cout << "rew " << nRewards << " pun " << nPunishments << endl;
		nRewards = nPunishments = 0;
	}

	return true;
}

PING_PONG_ENVIRONMENT_EXPORT int Finalize(int OriginalTerminationCode) 
{
	cout << "rew " << nRewards << " pun " << nPunishments << endl;

    //double dmean, dstderr;
    //avgdis(&v_err.front(), v_err.size(), dmean, dstderr);
    //cout << "err: mean " << dmean << " stderr " << dstderr << endl;

	return nRewardsTot * 10000 / (nPunishmentsTot + nRewardsTot);
}

PING_PONG_ENVIRONMENT_EXPORT void Serialize(Serializer &ser, bool bSave) {}
