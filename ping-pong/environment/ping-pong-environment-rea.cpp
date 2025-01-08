/*
Copyright (C) 2022 Mikhail Kiselev

Emulates ping-pong game.

*/

#include <iostream>
#include "ping-pong-environment.h"
#include "EnvironmentState.hpp"

#include <ArNI/../helpers/DatTab/RF/RF.h>

using namespace std;

extern int tactStart;
extern float rStateFiringFrequency;
extern int ntact;
extern EnvironmentState es;

extern int nRewards;
extern int nPunishments;

extern int nSecondaryRewards;
extern int nSecondaryPunishments;

extern int nRewardsTot;
extern int nPunishmentsTot;

extern int LastRegistration;
extern int PostRewardCounter;
extern int nRecognitions;
extern int nCorr;
extern bool bCurrentStateOK;

extern GrowingStimulation *pgsG;

extern std::unique_ptr<RFRes> uprfr;
extern RF rf;

int nNeuronsperAction = 0;
int nNeuronsperConsequence = 0;

int nerrRF = 0;
int ncorrRF = 0;

int tactlastPunishmentPrediction = -1000000;
int tactlastRewardPrediction = -1000000;

int nPredictedPunishments = 0;
int nPredictedRewards = 0;
int nPredictedPunishmentsTot = 0;
int nPredictedRewardsTot = 0;

int PredictionPeriod = 0;

int tactPredictedPunishmentLim = -1;
int tactPredictedRewardLim = -1;

int nFalsePredictedPunishments = 0;
int nFalsePredictedRewards = 0;
int nFalsePredictedPunishmentsTot = 0;
int nFalsePredictedRewardsTot = 0;

PING_PONG_ENVIRONMENT_EXPORT void SetParametersOut(int ExperimentId, size_t tactTermination, unsigned nOutputNeurons, const pugi::xml_node &xn) 
{
	tactStart = atoi_s(xn.child("start_time").child_value());
    PredictionPeriod = atoi_s(xn.child("prediction_period").child_value());
    nNeuronsperConsequence = nOutputNeurons / 2;
}

const float rAction = 1.F / nSpatialZones;

PING_PONG_ENVIRONMENT_EXPORT bool ObtainOutputSpikes(const vector<int> &v_Firing, int nEquilibriumPeriods)
{
    static Rand rng;
    auto r = !(ntact % 100) ? (rng() < 0.5 ? -rAction : rAction) : 0.F;
	*es.prRacket += r;
	if (*es.prRacket > 0.5F - RACKET_SIZE / 2)
		*es.prRacket = 0.5F - RACKET_SIZE / 2;
	else if (*es.prRacket < -0.5F + RACKET_SIZE / 2)
		*es.prRacket = -0.5F + RACKET_SIZE / 2;

    for (auto i: v_Firing) {
        if (i < nNeuronsperConsequence) {
            tactlastPunishmentPrediction = ntact;
            if (tactPredictedPunishmentLim == -1)
                tactPredictedPunishmentLim = ntact + 3 * PredictionPeriod;
        } else {
            tactlastRewardPrediction = ntact;
            if (tactPredictedRewardLim == -1)
                tactPredictedRewardLim = ntact + 3 * PredictionPeriod;
        }
    }
    if (ntact == tactPredictedPunishmentLim) {
        tactPredictedPunishmentLim = -1;
        if (ntact >= tactStart)
            ++nFalsePredictedPunishmentsTot;
        ++nFalsePredictedPunishments;
    }
    if (ntact == tactPredictedRewardLim) {
        tactPredictedRewardLim = -1;
        if (ntact >= tactStart)
            ++nFalsePredictedRewardsTot;
        ++nFalsePredictedRewards;
    }

    if (ntact && !(ntact % 200000)) {
        float rPrecisionPunishment = nPredictedPunishments + nFalsePredictedPunishments ? (float)nPredictedPunishments / (nPredictedPunishments + nFalsePredictedPunishments) : 0.F;
        float rRecallPunishment = (float)nPredictedPunishments / nPunishments;
        float rPrecisionReward =nPredictedRewards + nFalsePredictedRewards ? (float)nPredictedRewards / (nPredictedRewards + nFalsePredictedRewards) : 0.F;
        float rRecallReward = (float)nPredictedRewards / nRewards;
        float rFPunishment = rPrecisionPunishment && rRecallPunishment ? 2 / (1 / rPrecisionPunishment + 1 / rRecallPunishment) : 0.F;
        float rFReward = rPrecisionReward && rRecallReward ? 2 / (1 / rPrecisionReward + 1 / rRecallReward) : 0.F;
        cout << "PUNISHMENTS: precision = " << rPrecisionPunishment << " recall = " << rRecallPunishment << " F = " << rFPunishment << "\nREWARD: precision = " << rPrecisionReward << " recall = " << rRecallReward << " F = " << rFReward << endl;
        nPredictedPunishments = nFalsePredictedPunishments = nPunishments = nPredictedRewards = nFalsePredictedRewards = nRewards = 0;
    }

	return true;
}

PING_PONG_ENVIRONMENT_EXPORT int Finalize(int OriginalTerminationCode) 
{
    float rPrecisionPunishment = nPredictedPunishmentsTot + nFalsePredictedPunishmentsTot ? (float)nPredictedPunishmentsTot / (nPredictedPunishmentsTot + nFalsePredictedPunishmentsTot) : 0.F;
    float rRecallPunishment = (float)nPredictedPunishmentsTot / nPunishmentsTot;
    float rPrecisionReward = nPredictedRewardsTot + nFalsePredictedRewardsTot ? (float)nPredictedRewardsTot / (nPredictedRewardsTot + nFalsePredictedRewardsTot) : 0.F;
    float rRecallReward = (float)nPredictedRewardsTot / nRewardsTot;
    float rFPunishment = rPrecisionPunishment && rRecallPunishment ? 2 / (1 / rPrecisionPunishment + 1 / rRecallPunishment) : 0.F;
    float rFReward = rPrecisionReward && rRecallReward ? 2 / (1 / rPrecisionReward + 1 / rRecallReward) : 0.F;
    return 10000 * (rFPunishment + rFReward) / 2;
}

PING_PONG_ENVIRONMENT_EXPORT void Serialize(Serializer &ser, bool bSave) {}
