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
#include <deque>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/mapped_region.hpp>

#include <sg/sg.h>
#include <NetworkConfigurator.h>

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
const int nVelocityZones = 9;
const int nRelPos = 5;
const float rRelPosStep = RACKET_SIZE / 3;   // Racket takes 3 middle positions of the nRelPos x nRelPos grid.
const unsigned nInputs = 3 * nSpatialZones + 2 * nVelocityZones + nRelPos * nRelPos;

const float rAction = 1.F / nSpatialZones;

const float rStateFiringFrequency = 0.3F;

#define SIGNAL_ON rng() < rStateFiringFrequency

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
	float rBallMovementDirection;
	do {
		rBallMovementDirection = (float)rng(2 * M_PI);
		prr_BallSpeed.first = rBallVelocity * sin(rBallMovementDirection);
	} while (prr_BallSpeed.first < 1. / (2 * maxSpotPassageTime_ms));
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

int ntact = 0;

// const double adLevelBoundaries[] = {0.27915265777848186 * 0.27915265777848186, 0.2997424323017368 * 0.2997424323017368, 0.3435725359161685 * 0.3435725359161685, 0.4109012910954389 * 0.4109012910954389};
const double adLevelBoundaries[] = {0.04, 0.16, 0.36, 0.64};
const int nGoalLevels = sizeof(adLevelBoundaries) / sizeof(adLevelBoundaries[0]);
const int LevelDuration = 100;
deque<vector<bool> > qvb_Neuron, qvb_;
const int NeuronDepth = 30;
vector<int> vn_Neuron(nInputs, 0);
vector<bool> vb_Neuron(nInputs, false);
int CurrentLevel = nGoalLevels;
const int prerewardperiod = LevelDuration * nGoalLevels;
map<vector<int>, vector<int> > mvindvn_;
int BayesModel[nGoalLevels + 1][nInputs];
int nRewardedTacts = 0;
ofstream ofsrews("rews.csv");
ofstream ofsBayes("Bayes.csv");
double d2cur = 0.;

class DYNAMIC_LIBRARY_EXPORTED_CLASS rec_ping_pong: public IReceptors
{
	vector<float> vr_VelocityZoneBoundary;
protected:
	virtual bool bGenerateReceptorSignals(char *prec, size_t neuronstrsize) override
	{
		static ofstream ofsState("ping_pong_state.csv");
		vector<float> vr_PhaseSpacePoint(5);

		vector<int> vind_(6);
#define indxBall vind_[0]
#define indyBall vind_[1]
#define indvxBall vind_[2]
#define indvyBall vind_[3]
#define indRacket vind_[4]
#define indRaster vind_[5]

		UpdateWorld(vr_PhaseSpacePoint);
		for (auto z: vr_PhaseSpacePoint)
			ofsState << z << ',';
		ofsState << endl;
		indxBall = (int)((vr_PhaseSpacePoint[0] + 0.5) / (1. / nSpatialZones));
		if (indxBall == nSpatialZones)
			indxBall = nSpatialZones - 1;
		indyBall = (int)((vr_PhaseSpacePoint[1] + 0.5) / (1. / nSpatialZones));
		if (indyBall == nSpatialZones)
			indyBall = nSpatialZones - 1;
		indvxBall = (int)(lower_bound(vr_VelocityZoneBoundary.begin(), vr_VelocityZoneBoundary.end(), abs(vr_PhaseSpacePoint[2])) - vr_VelocityZoneBoundary.begin());
		indvxBall = vr_PhaseSpacePoint[2] < 0 ? nVelocityZones / 2 - indvxBall : nVelocityZones / 2 + indvxBall;
		indvyBall = (int)(lower_bound(vr_VelocityZoneBoundary.begin(), vr_VelocityZoneBoundary.end(), abs(vr_PhaseSpacePoint[3])) - vr_VelocityZoneBoundary.begin());
		indvyBall = vr_PhaseSpacePoint[3] < 0 ? nVelocityZones / 2 - indvyBall : nVelocityZones / 2 + indvyBall;
		indRacket = (int)((vr_PhaseSpacePoint[4] + 0.5) / (1. / nSpatialZones));
		if (indRacket == nSpatialZones)
			indRacket = nSpatialZones - 1;
		vector<bool> vb_Spikes(nInputs, false);
		vb_Spikes[indxBall] = SIGNAL_ON;
		vb_Spikes[nSpatialZones + indyBall] = SIGNAL_ON;
		vb_Spikes[nSpatialZones * 2 + indvxBall] = SIGNAL_ON;
		vb_Spikes[nSpatialZones * 2 + nVelocityZones + indvyBall] = SIGNAL_ON;
		vb_Spikes[nSpatialZones * 2 + nVelocityZones * 2 + indRacket] = SIGNAL_ON;
		int indxRel = (int)((vr_PhaseSpacePoint[0] + 0.5) / rRelPosStep);
		if (indxRel < nRelPos) {
			int indyRel = (int)((vr_PhaseSpacePoint[4] - vr_PhaseSpacePoint[1] + rRelPosStep / 2) / rRelPosStep);   // Raster goes from top (higher y) to bottom - in opposite 
																 // direction to y axis
			if (abs(indyRel) <= (nRelPos - 1) / 2) {
				indyRel += (nRelPos - 1) / 2;
				indRaster = indyRel * nRelPos + indxRel;
				vb_Spikes[nSpatialZones * 3 + nVelocityZones * 2 + indRaster] = SIGNAL_ON;
			}
		}
		for (auto i: vb_Spikes) {
			*prec = i;
			prec += neuronstrsize;
		}

		double d2 = (-0.5 - vr_PhaseSpacePoint[0]) * (-0.5 - vr_PhaseSpacePoint[0]) + (vr_PhaseSpacePoint[1] - vr_PhaseSpacePoint[4]) * (vr_PhaseSpacePoint[1] - vr_PhaseSpacePoint[4]);
		CurrentLevel = lower_bound(adLevelBoundaries, adLevelBoundaries + nGoalLevels, d2) - adLevelBoundaries;

		return true;
	}
	virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const 
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

		memset(BayesModel, 0, sizeof(BayesModel));

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

int RewardPunishmentBalance = 0;
int nRewards = 0;
int nPunishments = 0;

const int RewardTrainLength = 10;
const int RewardTrainPeriod = 2;

bool b_forVerifier_Reward = false;

class DYNAMIC_LIBRARY_EXPORTED_CLASS Evaluator: public IReceptors
{
public:
	enum type
	{
		reward,
		punishment,
		_debug_rewnorm,
		_debug_punishment
	};
private:
	enum Evaluator::type typ;
	int                  TrainCounter;
	int                  PeriodCounter;

	int                  curlev = nGoalLevels;

protected:
	virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const
	{
		vstr_Meanings.resize(1);
		vstr_Meanings.front() = typ == reward ? "REW" : typ == punishment ? "PUN" : typ == _debug_rewnorm ? "$$$rewnorm" : "$$$punishment";
	}
public:
	Evaluator(enum Evaluator::type t) : IReceptors(1), typ(t) {}
	virtual bool bGenerateReceptorSignals(char *prec, size_t neuronstrsize) override
	{
		switch (typ) {
			case punishment: if (es.pprr_Ball->first < -0.5F) {
								TrainCounter = RewardTrainLength;
								PeriodCounter = 1;
								--RewardPunishmentBalance;
								++nPunishments;
							 }
							 break;
			case reward:     if (es.pprr_Ball->first == -0.5F) {
								TrainCounter = RewardTrainLength;
								PeriodCounter = 1;
								++RewardPunishmentBalance;
								++nRewards;
								b_forVerifier_Reward = true;
							 }
							 break;
			case _debug_rewnorm: *prec = CurrentLevel < curlev ? 1 : 0;
						     curlev = CurrentLevel;
							 if (*prec)
								 ofsrews << ntact << ',1,' << CurrentLevel << endl;
				             return true;
			default: *prec = CurrentLevel > curlev ? 1 : 0;
				curlev = CurrentLevel;
				if (*prec)
					ofsrews << ntact << ',0,' << CurrentLevel << endl;
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

const float rBasicPoissonFrequency = 0.0001F;   // Надо, чтобы за период допаминовой пластичности было примерно 1-2 случайных действия - не больше.
const float rMinTargetNetworkActivity = 0.01F;
const int NoisePeriod = 300000;

class DYNAMIC_LIBRARY_EXPORTED_CLASS AdaptivePoisson: public IReceptors
{
	bool bActionWasForcedbyNoise;
public:
	AdaptivePoisson() {}
	AdaptivePoisson(int nReceptors): IReceptors(nReceptors), rCurrentFrequency(rBasicPoissonFrequency / nReceptors), bActionWasForcedbyNoise(false) {}

	float rCurrentFrequency;

	virtual bool bGenerateReceptorSignals(char *prec, size_t neuronstrsize) override
	{
		rCurrentFrequency += rBasicPoissonFrequency * rMinTargetNetworkActivity;
		if (rCurrentFrequency > rBasicPoissonFrequency)
			rCurrentFrequency = rBasicPoissonFrequency;
		FORI(nReceptors) {
			if (ntact < NoisePeriod && rng() < rCurrentFrequency) {
				*prec = 1;
				bActionWasForcedbyNoise = true;
			} else *prec = 0;
			prec += neuronstrsize;
		}
		return true;
	}
	virtual void Randomize(void) override {}   // It was randomized by rec_ping_pong
	virtual void SaveStatus(Serializer &ser) const override
	{
		IReceptors::SaveStatus(ser);
		ser << rCurrentFrequency;
		ser << bActionWasForcedbyNoise;
	}
	virtual ~AdaptivePoisson() = default;

	void LoadStatus(Serializer &ser)
	{
		IReceptors::LoadStatus(ser);
		ser >> rCurrentFrequency;
		ser >> bActionWasForcedbyNoise;
	}
	void RegisterAction()
	{
		if (bActionWasForcedbyNoise)
			bActionWasForcedbyNoise = false;
		else rCurrentFrequency -= rBasicPoissonFrequency;
	}
};

AdaptivePoisson *papG;

PING_PONG_ENVIRONMENT_EXPORT IReceptors *SetParametersIn(int &nReceptors, const pugi::xml_node &xn)
{
	static int CallNo = 0;
	switch (CallNo++) {
		case 0: nReceptors = nInputs;
				return new rec_ping_pong;
		case 1: nReceptors = 1;
			    return new Evaluator(Evaluator::punishment);
		case 2: nReceptors = 1;
			    return new Evaluator(Evaluator::reward);

		case 3: nReceptors = 1;
				return new Evaluator(Evaluator::_debug_rewnorm);
		case 4: nReceptors = 1;
			    return new Evaluator(Evaluator::_debug_punishment);

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
		case 3: papG = new AdaptivePoisson(true);
				papG->LoadStatus(ser);
				return papG;
		default: cout << "Too many calls of LoadStatus\n";
				exit(-1);
	}
}

int nNeuronsperAction;

PING_PONG_ENVIRONMENT_EXPORT void SetParametersOut(int ExperimentId, size_t tactTermination, unsigned nOutputNeurons, const pugi::xml_node &xn) {nNeuronsperAction = (int)nOutputNeurons / 2;}

PING_PONG_ENVIRONMENT_EXPORT bool ObtainOutputSpikes(const vector<int> &v_Firing, int nEquilibriumPeriods)
{
	int nCommandsDown = count_if(v_Firing.begin(), v_Firing.end(), bind2nd(less<int>(), nNeuronsperAction));
	*es.prRacket += rAction * ((int)v_Firing.size() - 2 * nCommandsDown);
	if (*es.prRacket > 0.5F - RACKET_SIZE / 2)
		*es.prRacket = 0.5F - RACKET_SIZE / 2;
	else if (*es.prRacket < -0.5F + RACKET_SIZE / 2)
		*es.prRacket = -0.5F + RACKET_SIZE / 2;

	if (ntact && !(ntact % 200000)) {
		for (int z = 0; z <= nGoalLevels; ++z)
			FORI(nInputs)
				ofsBayes << BayesModel[z][_i] << (_i < nInputs - 1 ? ',' : '\n');
		ofsBayes << endl;
		cout << "rew " << nRewards << " pun " << nPunishments << endl;
		nRewards = nPunishments = 0;
	}

	return true;
}

int LastRegistration = -1000000;
int PostRewardCounter = 0;
int nRecognitions = 0;
int nCorr = 0;

PING_PONG_ENVIRONMENT_EXPORT int Finalize(int OriginalTerminationCode) 
{

	for (int z = 0; z <= nGoalLevels; ++z)
		FORI(nInputs)
			ofsBayes << BayesModel[z][_i] << (_i < nInputs - 1 ? ',' : '\n');

	cout << "rew " << nRewards << " pun " << nPunishments << endl;
	return 5000 + RewardPunishmentBalance;
}

PING_PONG_ENVIRONMENT_EXPORT void Serialize(Serializer &ser, bool bSave) {}

vector<pair<string, unsigned> > vpstrn_sec;
int nNeurons;
pair<int, int> p_LREWRange;
int indREWNORM;

PING_PONG_ENVIRONMENT_EXPORT void GetSections(const vector<pair<string, unsigned> > &vpstrn_Sections)
{
	vpstrn_sec = vpstrn_Sections;
	nNeurons = 0;
	for (const auto &i: vpstrn_Sections) {
		if (i.first == "LREW") {
			p_LREWRange.first = nNeurons;
			p_LREWRange.second = p_LREWRange.first + i.second;
		} else if (i.first == "REWNORM")
			indREWNORM = nNeurons;
		nNeurons += i.second;
	}
}

PING_PONG_ENVIRONMENT_EXPORT void ObtainSpikes(const vector<int> &v_Firing, string &strFatal)
{
	if (!PostRewardCounter && find(v_Firing.begin(), v_Firing.end(), indREWNORM) != v_Firing.end())
		++nRecognitions;
	if (!PostRewardCounter && any_of(v_Firing.begin(), v_Firing.end(), [=](int indneu) {return p_LREWRange.first <= indneu && indneu < p_LREWRange.second; }))
		LastRegistration = ntact;
	if (b_forVerifier_Reward) {
		if (ntact - LastRegistration < 100)
			++nCorr;
		PostRewardCounter = RewardTrainLength * RewardTrainPeriod;
		b_forVerifier_Reward = false;
	}
	++ntact;
	if (PostRewardCounter)
		--PostRewardCounter;
}
