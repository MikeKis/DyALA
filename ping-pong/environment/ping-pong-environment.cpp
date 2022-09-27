/*
Copyright (C) 2021 Mikhail Kiselev
Автор Михаил Витальевич Киселев.
При модификации файла сохранение указания (со)авторства Михаила Витальевича Киселева обязательно.

Emulates signal from videocamera looking at a moving light spot.

*/

#define _USE_MATH_DEFINES

#include <utility>
#include <math.h>
#include <random>
#include <fstream>
#include <vector>
#include <iostream>
#include <memory>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include <sg/sg.h>

#include "ping-pong-environment.h"

#ifdef FOR_LINUX
#define PING_PONG_ENVIRONMENT_EXPORT extern "C"
#define DYNAMIC_LIBRARY_EXPORTED_CLASS
typedef __uint64_t UNS64;
#else
#define PING_PONG_ENVIRONMENT_EXPORT __declspec(dllexport)
#define DYNAMIC_LIBRARY_EXPORTED_CLASS __declspec(dllexport)
typedef unsigned __int64 UNS64;
#endif

using namespace std;
using namespace boost::interprocess;

const unsigned minSpotPassageTime_ms = 300;
const unsigned maxSpotPassageTime_ms = 1000;

const unsigned nSpatialZones = 30;
const unsigned nVelocityZones = 9;
const unsigned nRelPos = 5;
const unsigned nInputs = 3 * nSpatialZones + 2 * nVelocityZones + nRelPos * nRelPos;

const float rAction = 1.F / nSpatialZones;

const float rStateFiringFrequency = 0.3F;

#define SIGNAL_ON (rng() < rStateFiringFrequency ? 1 : 0)

class RandomNumberGenerator
{
	uniform_real_distribution<> urd;
	mt19937_64                  mt;
public:
	RandomNumberGenerator() : urd(0., 1.) {}
	double operator()() { return urd(mt); }
	template<class T> T operator()(T max) { return (T)((*this)() * max); }
	void Randomize(void)
	{
		mt19937_64 mtRandomized(std::random_device{}());
		mt = mtRandomized;
	}
};

RandomNumberGenerator rng;

float rMakeBallVelocity(void)
{
	int i1, i2;
	do {
		i1 = minSpotPassageTime_ms + rng(maxSpotPassageTime_ms - minSpotPassageTime_ms);
		i2 = rng(maxSpotPassageTime_ms);
	} while (i2 > i1);
	return(1.F / i1);
}

class EnvironmentState
{
	unique_ptr<shared_memory_object> shm;
	unique_ptr<mapped_region>        region;
	string                           strSharedMemoryName;
public:
	pair<float, float> *pprr_Ball;
	float              *prRacket = NULL;
	EnvironmentState()
	{
		strSharedMemoryName = ENVIRONMENT_STATE_SHARED_MEMORY_NAME;
		bool bExists = true;
		do {
			try {
				//Create a shared memory object.
				shm.reset(new shared_memory_object(open_only, strSharedMemoryName.c_str(), read_only));
				++strSharedMemoryName.front();
			}
			catch (...) {
				bExists = false;
			}
		} while (bExists);
		//Create a shared memory object.
		shm.reset(new shared_memory_object(create_only, strSharedMemoryName.c_str(), read_write));

		//Set size
		shm->truncate(sizeof(pair<pair<float, float>, float>));

		//Map the whole shared memory in this process
		region.reset(new mapped_region(*shm, read_write));
		pprr_Ball = (pair<float, float> *)region->get_address();
	}
	void ResetBall();
	~EnvironmentState() {shared_memory_object::remove(strSharedMemoryName.c_str());}
} es;

pair<float, float> prr_BallSpeed;

void EnvironmentState::ResetBall()
{
	es.pprr_Ball->first = 0.F;
	es.pprr_Ball->second = (float)(-0.5 + rng());
	float rBallVelocity = rMakeBallVelocity();
	float rBallMovementDirection = (float)rng(2 * M_PI);
	prr_BallSpeed.first = rBallVelocity * sin(rBallMovementDirection);
	prr_BallSpeed.second = rBallVelocity * cos(rBallMovementDirection);
}

void UpdateWorld(vector<float> &vr_PhaseSpacePoint)
{
	if (!es.prRacket) {
		es.prRacket = &((pair<pair<float, float>, float> *)es.pprr_Ball)->second;
		*es.prRacket = (float)(-0.5 + rng());
		es.ResetBall();
	}

	vr_PhaseSpacePoint[0] = es.pprr_Ball->first;
	vr_PhaseSpacePoint[1] = es.pprr_Ball->second;
	vr_PhaseSpacePoint[2] = prr_BallSpeed.first;
	vr_PhaseSpacePoint[3] = prr_BallSpeed.second;
	vr_PhaseSpacePoint[4] = *es.prRacket;

	if (es.pprr_Ball->first < -0.5F)   // PUNISHMENT
		es.ResetBall();
	else {
		es.pprr_Ball->first += prr_BallSpeed.first;
		if (es.pprr_Ball->first < -0.5F && es.pprr_Ball->second > *es.prRacket - RACKET_SIZE / 2 && es.pprr_Ball->second < *es.prRacket + RACKET_SIZE / 2) {
			es.pprr_Ball->first = -0.5F;   // REWARD
			prr_BallSpeed.first = -prr_BallSpeed.first;
		}
		if (es.pprr_Ball->first > 0.5F) {
			es.pprr_Ball->first = 0.5F;
			prr_BallSpeed.first = -prr_BallSpeed.first;
		}
		es.pprr_Ball->second += prr_BallSpeed.second;
		if (es.pprr_Ball->second < -0.5F) {
			es.pprr_Ball->second = -0.5F;
			prr_BallSpeed.second = -prr_BallSpeed.second;
		}
		if (es.pprr_Ball->second > 0.5F) {
			es.pprr_Ball->second = 0.5F;
			prr_BallSpeed.second = -prr_BallSpeed.second;
		}
	}
}

class DYNAMIC_LIBRARY_EXPORTED_CLASS rec_ping_pong: public IReceptors
{
	vector<float> vr_VelocityZoneBoundary;
protected:
	virtual bool bGenerateReceptorSignals(char *prec, size_t neuronstrsize) override
	{
		vector<float> vr_PhaseSpacePoint(5);
		UpdateWorld(vr_PhaseSpacePoint);
		int indxBall = (int)((vr_PhaseSpacePoint[0] + 0.5) / (1. / nSpatialZones));
		int indyBall = (int)((vr_PhaseSpacePoint[1] + 0.5) / (1. / nSpatialZones));
		int indvxBall = lower_bound(vr_VelocityZoneBoundary.begin(), vr_VelocityZoneBoundary.end(), abs(vr_PhaseSpacePoint[2])) - vr_VelocityZoneBoundary.begin();
		indvxBall = vr_PhaseSpacePoint[2] < 0 ? nVelocityZones / 2 - indvxBall : nVelocityZones / 2 + indvxBall;
		int indvyBall = lower_bound(vr_VelocityZoneBoundary.begin(), vr_VelocityZoneBoundary.end(), abs(vr_PhaseSpacePoint[3])) - vr_VelocityZoneBoundary.begin();
		indvyBall = vr_PhaseSpacePoint[3] < 0 ? nVelocityZones / 2 - indvyBall : nVelocityZones / 2 + indvyBall;
		int indRacket = (int)((vr_PhaseSpacePoint[4] + 0.5) / (1. / nSpatialZones));
		vector<char> vb_Spikes(nInputs, 0);
		vb_Spikes[indxBall] = SIGNAL_ON;
		vb_Spikes[nSpatialZones + indyBall] = SIGNAL_ON;
		vb_Spikes[nSpatialZones * 2 + indvxBall] = SIGNAL_ON;
		vb_Spikes[nSpatialZones * 2 + nVelocityZones + indvxBall] = SIGNAL_ON;
		vb_Spikes[nSpatialZones * 2 + nVelocityZones * 2 + indRacket] = SIGNAL_ON;
		if (indxBall < nRelPos) {
			int indyRel = indyBall - indRacket;
			if (abs(indyRel) <= (nRelPos - 1) / 2) {
				indyRel += (nRelPos - 1) / 2;
				int indRaster = indyRel * nRelPos + indxBall;
				vb_Spikes[nSpatialZones * 3 + nVelocityZones * 2 + indRaster] = SIGNAL_ON;
			}
		}
		for (auto i: vb_Spikes) {
			*prec = i;
			prec += neuronstrsize;
		}
		return true;
	}
public:
	rec_ping_pong(): IReceptors(nInputs), vr_VelocityZoneBoundary((nVelocityZones - 1) / 2)
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

int ntact = 0;

int RewardPunishmentBalance = 0;

class DYNAMIC_LIBRARY_EXPORTED_CLASS Evaluator: public IReceptors
{
	bool bReward;
public:
	Evaluator(bool bRew) : IReceptors(1), bReward(bRew) {}
	virtual bool bGenerateReceptorSignals(char *prec, size_t neuronstrsize) override
	{
		if (!bReward) {
			if (es.pprr_Ball->first < -0.5F) {
				*prec = 1;
				--RewardPunishmentBalance;
			} else *prec = 0;
		} else if (es.pprr_Ball->first == -0.5F) {
			*prec = 1;
			++RewardPunishmentBalance;
		} else *prec = 0;
		return true;
	}
	virtual void Randomize(void) override {};
	virtual void SaveStatus(Serializer &ser) const override
	{
		IReceptors::SaveStatus(ser);
		ser << bReward;
	}
	virtual ~Evaluator() = default;
	void LoadStatus(Serializer &ser)
	{
		IReceptors::LoadStatus(ser);
		ser >> bReward;
	}
};

PING_PONG_ENVIRONMENT_EXPORT IReceptors *SetParametersIn(int &nReceptors, const pugi::xml_node &xn)
{
	static int CallNo = 0;
	switch (CallNo++) {
		case 0: nReceptors = nInputs;
				return new rec_ping_pong;
		case 1: return new Evaluator(false);
		case 2: return new Evaluator(true);
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
		case 1: peva = new Evaluator(false);
				peva->LoadStatus(ser);
				return peva;
		case 2: peva = new Evaluator(true);
				peva->LoadStatus(ser);
				return peva;
		default: cout << "Too many calls of SetParametersIn\n";
				exit(-1);
	}
}

PING_PONG_ENVIRONMENT_EXPORT void SetParametersOut(int ExperimentId, size_t tactTermination, unsigned nOutputNeurons, const pugi::xml_node &xn) {}

PING_PONG_ENVIRONMENT_EXPORT bool ObtainOutputSpikes(const vector<int> &v_Firing, int nEquilibriumPeriods)
{
	if (v_Firing.size())
		*es.prRacket += v_Firing.front() ? rAction : -rAction;
	return true;
}

PING_PONG_ENVIRONMENT_EXPORT int Finalize(int OriginalTerminationCode) {return 5000 + RewardPunishmentBalance;}
PING_PONG_ENVIRONMENT_EXPORT void Serialize(Serializer &ser, bool bSave) {}
