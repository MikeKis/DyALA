/*
Copyright (C) 2022 Mikhail Kiselev

Emulates ping-pong game.

*/

#define _USE_MATH_DEFINES

#include <math.h>
#include <random>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>

#include "ClusterBayes.hpp"
#include "EnvironmentState.hpp"
#include "AdaptiveSpikeSource.hpp"

using namespace std;
using namespace boost::interprocess;

extern EnvironmentState es;
extern RandomNumberGenerator rng;

extern pair<float, float> prr_BallSpeed;

extern bool b_forVerifier_Reward;

int nGoalLevels;

int ntact = 0;
int tactStart;

int nRewards = 0;
int nPunishments = 0;
int nRewardsTot = 0;
int nPunishmentsTot = 0;

unique_ptr<ClusterBayes> cb;

int CurrentLevel;
bool bTargetState = false;
int tactLevelFixed = -1000; // = never

int PunishedLevel;   // used only in surrogate Bayes

int NeuronTimeDepth = 10;    // P#1
int LevelDuration;
int LevelNeuronPeriod /* = LevelDuration / NeuronTimeDepth */;

class DYNAMIC_LIBRARY_EXPORTED_CLASS Evaluator: public IReceptors
{
public:
    enum type
    {
        reward,
        punishment,

        _debug_rewnorm,
        _debug_punishment,
        _level0
    };
private:
    enum Evaluator::type typ;
    int                  TrainCounter = 0;
    int                  PeriodCounter;

protected:
    virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const override
    {
        if (typ == reward || typ == punishment) {
            vstr_Meanings.resize(1);
            vstr_Meanings.front() = typ == reward ? "REW" : "PUN";
        } else {
            vstr_Meanings.resize(nGoalLevels);
            FORI(nGoalLevels) {
                stringstream ss;
                ss << (typ == _debug_rewnorm ? "$$$rewnorm" : "$$$punishment") << _i;
                vstr_Meanings[_i] = ss.str();
            }
        }
    }
public:
    Evaluator(enum Evaluator::type t, size_t tactbeg = 0) : IReceptors(t == reward || t == punishment || t == _level0 ? 1 : nGoalLevels), typ(t) {}
    virtual bool bGenerateReceptorSignals(char *prec, size_t neuronstrsize) override
    {
        switch (typ) {
            case punishment: if (es.pprr_Ball->first < -0.5F) {
                                TrainCounter = RewardTrainLength;
                                PeriodCounter = 1;
                                if (ntact >= tactStart)
                                    ++nPunishmentsTot;
                                ++nPunishments;
                             }
                             break;
            case reward:     if (es.pprr_Ball->first == -0.5F) {
                                TrainCounter = RewardTrainLength;
                                PeriodCounter = 1;
                                if (ntact >= tactStart)
                                    ++nRewardsTot;
                                ++nRewards;
                                cb->FixReward();
                                b_forVerifier_Reward = true;
                             }
                             break;
            case _debug_rewnorm: FORI(nGoalLevels)
                                    prec[_i * neuronstrsize] = 0;
                                 if (ntact >= tactLevelFixed + NeuronTimeDepth) {
                                     if (ntact == tactLevelFixed + NeuronTimeDepth) {
                                         int NewLevel = cb->Predict();
                                         if (NewLevel < CurrentLevel) {
                                             prec[neuronstrsize * NewLevel] = 1;
                                             PunishedLevel = -1;
                                         }
                                         else if (NewLevel == CurrentLevel)
                                             PunishedLevel = -1;
                                         else PunishedLevel = CurrentLevel;
                                         CurrentLevel = NewLevel;
                                     } else if (!(ntact % NeuronTimeDepth))
                                         CurrentLevel = cb->Predict();
                                 }
                                 return true;
            case _debug_punishment: FORI(nGoalLevels)
                                        prec[_i * neuronstrsize] = 0;
                     if (ntact == tactLevelFixed + NeuronTimeDepth && PunishedLevel >= 0)
                         prec[neuronstrsize * PunishedLevel] = 1;
                                     return true;
            default: *prec = bTargetState && !(ntact % 3);
                return true;
        }
        *prec = 0;
        if (TrainCounter && !--PeriodCounter) {
            *prec = 1;
            PeriodCounter = RewardTrainPeriod;
            --TrainCounter;
        }
        return true;
    }
    virtual void Randomize(void) override {};
    virtual void SaveStatus(Serializer &ser) const override
    {
        IReceptors::SaveStatus(ser);
        ser << typ;
    }
    virtual ~Evaluator() = default;
    void LoadStatus(Serializer &ser)
    {
        IReceptors::LoadStatus(ser);
        ser >> typ;
    }
};

vector<float> vr_CurrentPhaseSpacePoint(5);

class DYNAMIC_LIBRARY_EXPORTED_CLASS rec_ping_pong: public IReceptors
{
    vector<float> vr_VelocityZoneBoundary;
    vector<AdaptiveSpikeSource> vass_;
protected:
    virtual bool bGenerateReceptorSignals(char *prec, size_t neuronstrsize) override
    {
        vector<int> vind_(6);
#define indxBall vind_[0]
#define indyBall vind_[1]
#define indvxBall vind_[2]
#define indvyBall vind_[3]
#define indRacket vind_[4]
#define indRaster vind_[5]

        UpdateWorld(vr_CurrentPhaseSpacePoint);
        indxBall = (int)((vr_CurrentPhaseSpacePoint[0] + 0.5) / (1. / nSpatialZones));
        if (indxBall == nSpatialZones)
            indxBall = nSpatialZones - 1;
        indyBall = (int)((vr_CurrentPhaseSpacePoint[1] + 0.5) / (1. / nSpatialZones));
        if (indyBall == nSpatialZones)
            indyBall = nSpatialZones - 1;
        indvxBall = (int)(lower_bound(vr_VelocityZoneBoundary.begin(), vr_VelocityZoneBoundary.end(), abs(vr_CurrentPhaseSpacePoint[2])) - vr_VelocityZoneBoundary.begin());
        indvxBall = vr_CurrentPhaseSpacePoint[2] < 0 ? nVelocityZones / 2 - indvxBall : nVelocityZones / 2 + indvxBall;
        indvyBall = (int)(lower_bound(vr_VelocityZoneBoundary.begin(), vr_VelocityZoneBoundary.end(), abs(vr_CurrentPhaseSpacePoint[3])) - vr_VelocityZoneBoundary.begin());
        indvyBall = vr_CurrentPhaseSpacePoint[3] < 0 ? nVelocityZones / 2 - indvyBall : nVelocityZones / 2 + indvyBall;
        indRacket = (int)((vr_CurrentPhaseSpacePoint[4] + 0.5) / (1. / nSpatialZones));
        if (indRacket == nSpatialZones)
            indRacket = nSpatialZones - 1;
        vector<bool> vb_Spikes(nInputs, false);
        vb_Spikes[indxBall] = vass_[indxBall].bFire();
        vb_Spikes[nSpatialZones + indyBall] = vass_[nSpatialZones + indyBall].bFire();
        vb_Spikes[nSpatialZones * 2 + indvxBall] = vass_[nSpatialZones * 2 + indvxBall].bFire();
        vb_Spikes[nSpatialZones * 2 + nVelocityZones + indvyBall] = vass_[nSpatialZones * 2 + nVelocityZones + indvyBall].bFire();
        vb_Spikes[nSpatialZones * 2 + nVelocityZones * 2 + indRacket] = vass_[nSpatialZones * 2 + nVelocityZones * 2 + indRacket].bFire();
        int indxRel = (int)((vr_CurrentPhaseSpacePoint[0] + 0.5) / rRelPosStep);
        if (indxRel < nRelPos) {
            int indyRel = (int)((vr_CurrentPhaseSpacePoint[4] - vr_CurrentPhaseSpacePoint[1] + rRelPosStep / 2) / rRelPosStep);   // Raster goes from top (higher y) to bottom - in opposite
                                                                 // direction to y axis
            if (abs(indyRel) <= (nRelPos - 1) / 2) {
                indyRel += (nRelPos - 1) / 2;
                indRaster = indyRel * nRelPos + indxRel;
                vb_Spikes[nSpatialZones * 3 + nVelocityZones * 2 + indRaster] = vass_[nSpatialZones * 3 + nVelocityZones * 2 + indRaster].bFire();;
            }
        }
        for (auto i: vb_Spikes) {
            *prec = i;
            prec += neuronstrsize;
        }

        cb->AddNewInput(vb_Spikes);

        return true;
    }
    virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const override
    {
        vstr_Meanings.resize(nInputs);
        int i = 0;
        int x, y;
        for (x = 0; x < nSpatialZones; ++x) {
            stringstream ss;
            ss << "x" << x;
            vstr_Meanings[i++] = ss.str();
        }
        for (x = 0; x < nSpatialZones; ++x) {
            stringstream ss;
            ss << "y" << x;
            vstr_Meanings[i++] = ss.str();
        }
        for (x = 0; x < nVelocityZones; ++x) {
            stringstream ss;
            ss << "vx" << x - nVelocityZones / 2;
            vstr_Meanings[i++] = ss.str();
        }
        for (x = 0; x < nVelocityZones; ++x) {
            stringstream ss;
            ss << "vy" << x - nVelocityZones / 2;
            vstr_Meanings[i++] = ss.str();
        }
        for (x = 0; x < nSpatialZones; ++x) {
            stringstream ss;
            ss << "ry" << x;
            vstr_Meanings[i++] = ss.str();
        }
        for (y = nRelPos / 2; y >= -nRelPos / 2; --y)
            for (x = 0; x < nRelPos; ++x) {
                stringstream ss;
                ss << "REL(" << x << "," << y << ")";
                vstr_Meanings[i++] = ss.str();
            }
    }
public:
    rec_ping_pong(): IReceptors(nInputs), vr_VelocityZoneBoundary((nVelocityZones - 1) / 2), vass_(nInputs)
    {
        vector<float> vr_samples(9000);
        for (auto &i: vr_samples) {
            float rBallVelocity = rMakeBallVelocity();
            float rBallMovementDirection = (float)rng(M_PI / 2);
            i = rBallVelocity * sin(rBallMovementDirection);
        }
        sort(vr_samples.begin(), vr_samples.end());
        FORI((nVelocityZones - 1) / 2)
            vr_VelocityZoneBoundary[_i] = vr_samples[vr_samples.size() / 9 + _i * 2 * vr_samples.size() / 9];
        cb.reset(new ClusterBayes);

    }
    virtual void Randomize(void) override {rng.Randomize();}
    virtual void SaveStatus(Serializer &ser) const override
    {
        IReceptors::SaveStatus(ser);
        ser << *es.pprr_Ball;
        ser << *es.prRacket;
        ser << prr_BallSpeed;
        ser << vr_VelocityZoneBoundary;
        ser << rng;
    }
    virtual ~rec_ping_pong() = default;
    void LoadStatus(Serializer &ser)
    {
        IReceptors::LoadStatus(ser);
        ser >> *es.pprr_Ball;
        es.prRacket = &((pair<pair<float, float>, float> *)es.pprr_Ball)->second;
        ser >> *es.prRacket;
        ser >> prr_BallSpeed;
        ser >> vr_VelocityZoneBoundary;
        ser >> rng;
    }
};

#ifndef BIGREALNUMBER
#define BIGREALNUMBER 1.e10F
#endif

class DYNAMIC_LIBRARY_EXPORTED_CLASS Actions: public IReceptors
{
    float rPastRY = BIGREALNUMBER;
protected:
    virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const override
    {
        vstr_Meanings.resize(2);
        vstr_Meanings[0] = "ActDown";
        vstr_Meanings[1] = "ActUp";
    }
public:
    Actions(): IReceptors(2) {}
    virtual bool bGenerateReceptorSignals(char *prec, size_t neuronstrsize) override
    {
        if (rPastRY == BIGREALNUMBER || rPastRY == vr_CurrentPhaseSpacePoint[4])
            prec[neuronstrsize] = prec[0] = 0;
        else if (rPastRY > vr_CurrentPhaseSpacePoint[4]) {
            prec[neuronstrsize] = 0;
            prec[0] = 1;
        } else {
            prec[neuronstrsize] = 1;
            prec[0] = 0;
        }
        rPastRY = vr_CurrentPhaseSpacePoint[4];
        return true;
    }
    virtual void Randomize(void) override {};
    virtual void SaveStatus(Serializer &ser) const override
    {
        IReceptors::SaveStatus(ser);
        ser << rPastRY;
    }
    virtual ~Actions() = default;
    void LoadStatus(Serializer &ser)
    {
        IReceptors::LoadStatus(ser);
        ser >> rPastRY;
    }
};

PING_PONG_ENVIRONMENT_EXPORT IReceptors *SetParametersIn(int &nReceptors, const pugi::xml_node &xn)
{
	static int CallNo = 0;
	switch (CallNo++) {
		case 0: nReceptors = 1;
			    return new Evaluator(Evaluator::punishment);
		case 1: nReceptors = 1;
			    return new Evaluator(Evaluator::reward);

        case 2: CurrentLevel = nGoalLevels = nReceptors;
				return new Evaluator(Evaluator::_debug_rewnorm);
		case 3: nReceptors = nGoalLevels;
			    return new Evaluator(Evaluator::_debug_punishment);
		case 4: nReceptors = nInputs;
			    return new rec_ping_pong;

		case 5: nReceptors = 1;
			    return new Evaluator(Evaluator::_level0);

        case 6: nReceptors = 2;
                return new Actions;
        default: cout << "Too many calls of SetParametersIn\n";
				exit(-1);
	}
}

PING_PONG_ENVIRONMENT_EXPORT IReceptors *LoadStatus(Serializer &ser)
{
	static int CallNo = 0;
	rec_ping_pong *prpp;
	Evaluator *peva;
	switch (CallNo++) {
		case 0: prpp = new rec_ping_pong;
				prpp->LoadStatus(ser);
				return prpp;
		case 1: peva = new Evaluator(Evaluator::punishment);
				peva->LoadStatus(ser);
				return peva;
		case 2: peva = new Evaluator(Evaluator::reward);
				peva->LoadStatus(ser);
				return peva;
//		case 3: papG = new AdaptivePoisson(true);
//				papG->LoadStatus(ser);
//				return papG;
		default: cout << "Too many calls of LoadStatus\n";
				exit(-1);
	}
}
