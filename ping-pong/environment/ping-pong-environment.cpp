/*
Copyright (C) 2022 Mikhail Kiselev

Emulates ping-pong game.

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
#include <boost/numeric/ublas/symmetric.hpp>

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
using namespace boost::numeric;

DECLARE_UNION_OPERATORS(set<unsigned>)

//----------------------------------------------------------------------------------------------------------
// World model constants

const unsigned minSpotPassageTime_ms = 300;
const unsigned maxSpotPassageTime_ms = 1000;

const unsigned nSpatialZones = 30;
const int nVelocityZones = 9;
const int nRelPos = 5;
const float rRelPosStep = RACKET_SIZE / 3;   // Racket takes 3 middle positions of the nRelPos x nRelPos grid.
const unsigned nInputs = 3 * nSpatialZones + 2 * nVelocityZones + nRelPos * nRelPos;

const unsigned aAlternativeInputBoundaries[] = {30, 60, 69, 78, 108};

const float rAction = 1.F / nSpatialZones;
//----------------------------------------------------------------------------------------------------------

int NeuronTimeDepth = 10;    // P#1
float rStateFiringFrequency /* = 5.F / NeuronTimeDepth */;
int LevelNeuronPeriod /* = LevelDuration / NeuronTimeDepth */;
int minnSignificantExtraSpikes = 10;

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
	pair<pair<float, float>, float>  pprrr_LocalState;
public:
	pair<float, float> *pprr_Ball;
	float              *prRacket = NULL;
	EnvironmentState()
	{
		strSharedMemoryName = ENVIRONMENT_STATE_SHARED_MEMORY_NAME;
		bool bExists = true;
		try {
			shm.reset(new shared_memory_object(open_only, strSharedMemoryName.c_str(), read_write));
		} catch (...) {
			bExists = false;
		}
		if (bExists) {
			//Map the whole shared memory in this process
			region.reset(new mapped_region(*shm, read_write));
			pprr_Ball = (pair<float, float> *)region->get_address();
		} else pprr_Ball = &pprrr_LocalState.first;
	}
	void ResetBall(bool bPunishment = true);
} es;

pair<float, float> prr_BallSpeed;

void EnvironmentState::ResetBall(bool bPunishment)
{
	if (bPunishment) {
		es.pprr_Ball->first = 0.F;
		es.pprr_Ball->second = (float)(-0.5 + rng());
	}
	float rBallVelocity = rMakeBallVelocity();
	float rBallMovementDirection;
	do {
		rBallMovementDirection = (float)rng(2 * M_PI);
		prr_BallSpeed.first = rBallVelocity * sin(rBallMovementDirection);
	} while (prr_BallSpeed.first < 1. / (2 * maxSpotPassageTime_ms) || !bPunishment && prr_BallSpeed.first < 0.);
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
		es.ResetBall(true);
	else {
		bool bnegx = es.pprr_Ball->first < 0.;
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
		if (bnegx && es.pprr_Ball->first > 0.)
			es.ResetBall(false);
	}
}

int ntact = 0;
int nGoalLevels;
int CurrentLevel;
int tactLevelFixed = -1000; // = never

int PunishedLevel;   // used only in surrogate Bayes

class ClusterBayes
{
    map<vector<unsigned>, vector<double> >     mvvd_BayesModel;
    vector<double>                             vd_BayesPriors;
    deque<vector<bool> >                       qvb_RecentNeuronInput;

    deque<int>                                 q_TrueLevel;   // used only in this version of surragate Bayes neuron mechanism.

    deque<vector<bool> >                       qvb_NeuronInput;
    vector<vector<vector<bool> > >             vvvb_byLevels;
    vector<int>                                vn_Neuron;
    size_t                                     ntotMeasurements = 0;
    ublas::symmetric_matrix<vector<unsigned> > smvn_Pairs;
    map<unsigned,map<unsigned, double> >       mmd_AssociatedPairs;
    size_t nLevelMeasurements(int Level) const {return vvvb_byLevels[Level].size();}
public:
    ClusterBayes(): vvvb_byLevels(nGoalLevels + 1), vn_Neuron(nInputs, 0), smvn_Pairs(nInputs), vd_BayesPriors(nGoalLevels + 1, 0.)
    {
		FORI(nInputs)
			for (int i = _i; i < nInputs; ++i)
				smvn_Pairs(i, _i).resize(nGoalLevels + 1, 0);
    }
    void AddNewInput(const vector<bool> &vb_Spikes);
	int Predict();
	void FixReward();
    bool bClusterFeatureIsOK(const set<unsigned> &s_) const
    {
        int j;
		set<unsigned>::const_iterator i;
		const unsigned *p1 = NULL;
		FOR_ALL(i, s_) {
			auto p = upper_bound(aAlternativeInputBoundaries, aAlternativeInputBoundaries + sizeof(aAlternativeInputBoundaries) / sizeof(aAlternativeInputBoundaries[0]), *i);
			if (p == p1)
				return false;
			p1 = p;
		}
        FOR_(j, nGoalLevels) {
            auto n = count_if(vvvb_byLevels[j].begin(), vvvb_byLevels[j].end(), [&](const vector<bool> &vb_){return all_of(s_.begin(), s_.end(), [&](unsigned ind){return vb_[ind];});});
            if (n >= 4) {
                double d = 1.;
                set<unsigned>::const_iterator k;
                FORI(s_.size()) {
                    d *= smvn_Pairs(_i, _i)[j];
                    if (_i)
                        d /= nLevelMeasurements(j);
                }
                if ((n - d) * s_.size() > minnSignificantExtraSpikes)
                    break;
            }
        }
        return j < nGoalLevels;
    }
};

template<class LINKS, class CLUSTER, class MERGER> void AgglomerativeClustering(const LINKS &PairedAssociations, LIST<CLUSTER> &l_Result, const MERGER &mer)
{
    typename LIST<CLUSTER>::iterator ExistingCluster;
    for (const auto &i: PairedAssociations) {
        auto i1 = find_if(l_Result.begin(), l_Result.end(), [&](const CLUSTER &c) {return c.bIncludes(i.Node1());});
        auto i2 = find_if(l_Result.begin(), l_Result.end(), [&](const CLUSTER &c) {return c.bIncludes(i.Node2());});
        CLUSTER *pnewcluster;
        if (i1 == l_Result.end() && i2 == l_Result.end())
            l_Result.push_back(CLUSTER(i));
        else if (i1 == l_Result.end() || i2 == l_Result.end()) {
            if (i1 == l_Result.end()) {
                ExistingCluster = i2;
                auto FreeObject = i.Node1();
                auto ObjectinExistingCluster = i.Node2();
                pnewcluster = mer(*ExistingCluster, FreeObject, ObjectinExistingCluster, i);
            } else {
                ExistingCluster = i1;
                auto FreeObject = i.Node2();
                auto ObjectinExistingCluster = i.Node1();
                pnewcluster = mer(*ExistingCluster, FreeObject, ObjectinExistingCluster, i);
            }
            if (pnewcluster) {
                *ExistingCluster = *pnewcluster;
                delete pnewcluster;
            }
        } else if (i1 != i2) {
            pnewcluster = mer(*i1, *i2, i);
            if (pnewcluster) {
                *i1 = *pnewcluster;
                l_Result.erase(i2);
                delete pnewcluster;
            }
        } else i1->AddLink(i);
    }
}

struct FeaturePair: public pair<double, pair<unsigned, unsigned> >
{
    FeaturePair(double dsig, unsigned ind1, unsigned ind2)
    {
        first = dsig;
        second.first = ind1;
        second.second = ind2;
    }
    unsigned Node1() const {return second.first;}
    unsigned Node2() const {return second.second;}
};

struct ClusteredFeature: public set<unsigned>
{
    ClusteredFeature(const FeaturePair &fp)
    {
        insert(fp.Node1());
        insert(fp.Node2());
    }
    bool bIncludes(unsigned indInput) const {return find(indInput) != end();}
    void AddLink(const FeaturePair &fp){}
};

class ClusteredFeatureMerger
{
    const ClusterBayes &cb;
public:
    ClusteredFeatureMerger(const ClusterBayes &cb_): cb(cb_) {}
    ClusteredFeature *operator()(const ClusteredFeature &cf, unsigned indNewInput, unsigned indExistingInput, const FeaturePair &fp) const
    {
        ClusteredFeature *pcfnew = new ClusteredFeature(cf);
        pcfnew->insert(indNewInput);
        if (cb.bClusterFeatureIsOK(*pcfnew))
            return pcfnew;
        delete pcfnew;
        return nullptr;
    }
    ClusteredFeature *operator()(const ClusteredFeature &cf1, const ClusteredFeature &cf2, const FeaturePair &fp) const
    {
        ClusteredFeature *pcfnew = new ClusteredFeature(cf1);
        *pcfnew += cf2;
        if (cb.bClusterFeatureIsOK(*pcfnew))
            return pcfnew;
        delete pcfnew;
        return nullptr;
    }
};

ofstream ofsState("ping_pong_state.csv");
vector<float> vr_CurrentPhaseSpacePoint(5);

int TrueCurrentLevel()
{
    double d = sqrt((vr_CurrentPhaseSpacePoint[0] + 0.5) * (vr_CurrentPhaseSpacePoint[0] + 0.5) + (vr_CurrentPhaseSpacePoint[1] - vr_CurrentPhaseSpacePoint[4]) * (vr_CurrentPhaseSpacePoint[1] - vr_CurrentPhaseSpacePoint[4]));
    return min((int)(d / 0.1), 4);
}

void ClusterBayes::AddNewInput(const vector<bool> &vb_Spikes)
{
	FORI(nInputs)
		if (vb_Spikes[_i])
			++vn_Neuron[_i];
	qvb_NeuronInput.push_front(vb_Spikes);
	if (qvb_NeuronInput.size() > NeuronTimeDepth) {
		FORI(nInputs)
			if (qvb_NeuronInput.back()[_i])
				--vn_Neuron[_i];
		qvb_NeuronInput.pop_back();
	}
	if (ntact && !(ntact % NeuronTimeDepth)) {
		qvb_RecentNeuronInput.push_front(vector<bool>(nInputs, false));
		FORI(nInputs)
			if (vn_Neuron[_i])
				qvb_RecentNeuronInput.front()[_i] = true;

        q_TrueLevel.push_front(TrueCurrentLevel());

    }
}

size_t tactStart = 0;
vector<int> v_err;

int ClusterBayes::Predict()
{
//	if (CurrentLevel != 1000)
//        return TrueCurrentLevel();
    unsigned j;
    if (nLevelMeasurements(0)) {
        set<FeaturePair, greater<FeaturePair> > asspairs;
		vector<bool> vb_SingleFeatures(nInputs, false);
		for (size_t k = 0; k <= sizeof(aAlternativeInputBoundaries) / sizeof(aAlternativeInputBoundaries[0]); ++k) {
			pair<int, int> p_(!k ? 0 : aAlternativeInputBoundaries[k - 1], k < sizeof(aAlternativeInputBoundaries) / sizeof(aAlternativeInputBoundaries[0]) ? aAlternativeInputBoundaries[k] : nInputs);
			int nmax = -1;
			int indmax;
			int ind = p_.first + ntact % (p_.second - p_.first);
			FORI(p_.second - p_.first) {
				if (vn_Neuron[ind] > nmax) {
					nmax = vn_Neuron[ind];
					indmax = ind;
				}
				if (++ind == p_.second)
					ind = p_.first;
			}
			if (nmax)
				vb_SingleFeatures[indmax] = true;
		}
		FOR_(j, nInputs - 1)
            if (vb_SingleFeatures[j] && mmd_AssociatedPairs.find(j) != mmd_AssociatedPairs.end())
                for (const auto &i: mmd_AssociatedPairs[j])
					if (vb_SingleFeatures[i.first])
						asspairs.insert(FeaturePair(i.second, j, i.first));
        list<ClusteredFeature> lcf_;
        ClusteredFeatureMerger cfm(*this);
        AgglomerativeClustering(asspairs, lcf_, cfm);
        vector<map<vector<unsigned>, vector<double> >::iterator> CurrentFeatures;
        for (const auto &i: lcf_) {
            auto p_ = mvvd_BayesModel.insert(pair<vector<unsigned>, vector<double> >(vector<unsigned>(i.begin(), i.end()), vector<double>(nGoalLevels + 1, 0.)));
            CurrentFeatures.push_back(p_.first);
            if (p_.second)
                FORI(nGoalLevels + 1)
                    p_.first->second[_i] = count_if(vvvb_byLevels[_i].begin(), vvvb_byLevels[_i].end(), [=](const vector<bool> &vb_){return all_of(p_.first->first.begin(), p_.first->first.end(), [&](unsigned ind){return vb_[ind];});}) /
                                           (double)nLevelMeasurements(_i);
			for (auto k: i)
				vb_SingleFeatures[k] = false;
        }
		FOR_(j, nInputs)
			if (vb_SingleFeatures[j]) {
				auto p_ = mvvd_BayesModel.insert(pair<vector<unsigned>, vector<double> >(vector<unsigned>(1, j), vector<double>(nGoalLevels + 1, 0.)));
				CurrentFeatures.push_back(p_.first);
				if (p_.second)
					FORI(nGoalLevels + 1)
						p_.first->second[_i] = count_if(vvvb_byLevels[_i].begin(), vvvb_byLevels[_i].end(), [=](const vector<bool> &vb_) {return vb_[p_.first->first.front()];}) / (double)nLevelMeasurements(_i);
			}
        double dpmax = 0;
        int PredictedLevel = CurrentLevel;
        FORI(nGoalLevels + 1) {
            double dp = vd_BayesPriors[_i];
            for (auto i: CurrentFeatures) {
                dp *= i->second[_i];
                if (!dp)
                    break;
            }
            if (dp > dpmax) {
                dpmax = dp;
                PredictedLevel = _i;
            }
        }

		if (ofsState.is_open()) {
            ofsState << ntact << ',' << PredictedLevel << ',' << TrueCurrentLevel();
			for (auto z: vr_CurrentPhaseSpacePoint)
				ofsState << ',' << z;
			ofsState << endl;
		}
        if (ntact >= tactStart)
            v_err.push_back(PredictedLevel - CurrentLevel);

		return PredictedLevel;
    }
    return CurrentLevel;
}

void ClusterBayes::FixReward()
{
    unsigned j, lev;
    size_t l, m;
    ntotMeasurements += qvb_RecentNeuronInput.size();
//    for (lev = 0; lev < nGoalLevels && qvb_RecentNeuronInput.size(); ++lev) {
//        auto i = qvb_RecentNeuronInput.size() > LevelNeuronPeriod ? qvb_RecentNeuronInput.begin() + LevelNeuronPeriod : qvb_RecentNeuronInput.end();
//        vvvb_byLevels[lev].insert(vvvb_byLevels[lev].end(), qvb_RecentNeuronInput.begin(),  i);
//        auto k = qvb_RecentNeuronInput.begin();
//        while (k != i) {
//            FORI(nInputs)
//                if ((*k)[_i]) {
//                    ++smvn_Pairs(_i, _i)[lev];
//                    FOR_(l, _i)
//                        if ((*k)[l])
//                            ++smvn_Pairs(_i, l)[lev];
//                }
//            ++k;
//        }
//        qvb_RecentNeuronInput.erase(qvb_RecentNeuronInput.begin(), i);
//    }
//    if (lev == nGoalLevels) {
//        for (const auto &m: vvvb_byLevels[lev])
//            FORI(nInputs)
//                if (m[_i]) {
//                    ++smvn_Pairs(_i, _i)[lev];
//                    FOR_(l, _i)
//                        if (m[l])
//                            ++smvn_Pairs(_i, l)[lev];
//                }
//        vvvb_byLevels[lev].insert(vvvb_byLevels[lev].end(), qvb_RecentNeuronInput.begin(),  qvb_RecentNeuronInput.end());
//    }
//    qvb_RecentNeuronInput.clear();

    FOR_(m, q_TrueLevel.size()) {
        vvvb_byLevels[q_TrueLevel[m]].push_back(qvb_RecentNeuronInput[m]);
        FORI(nInputs)
            if (qvb_RecentNeuronInput[m][_i]) {
                ++smvn_Pairs(_i, _i)[q_TrueLevel[m]];
                FOR_(l, _i)
                    if (qvb_RecentNeuronInput[m][l])
                        ++smvn_Pairs(_i, l)[q_TrueLevel[m]];
            }
    }
    q_TrueLevel.clear();
    qvb_RecentNeuronInput.clear();

    mmd_AssociatedPairs.clear();
	FORI(nGoalLevels) {
		j = 0;
		for (auto o: aAlternativeInputBoundaries) 
			for (; j < o; ++j) {
				auto n = smvn_Pairs(j, j)[_i];
				for (l = o; l < nInputs; ++l) {
					auto n1 = smvn_Pairs(l, j)[_i];
					if (n1 > minnSignificantExtraSpikes / 2) {
						double d = n1 - n * (double)smvn_Pairs(l, l)[_i] / nLevelMeasurements(_i);
						if (d > minnSignificantExtraSpikes / 2)
							mmd_AssociatedPairs[j][(unsigned)l] = d;
					}
				}
			}
	}
    mvvd_BayesModel.clear();
    FORI(nGoalLevels + 1)
        vd_BayesPriors[_i] = nLevelMeasurements(_i) / (double)ntotMeasurements;
}

deque<vector<bool> > qvb_Neuron, qvb_;
unique_ptr<ClusterBayes> cb;
map<vector<int>, vector<int> > mvindvn_;
int nRewardedTacts = 0;
double d2cur = 0.;

class AdaptiveSpikeSource
{
    int LastTactinThisState;
    float rCurrentFrequency;
public:
    AdaptiveSpikeSource(): LastTactinThisState(-2) {}
    bool bFire()
    {
        if (LastTactinThisState < ntact - 1)
            rCurrentFrequency = rStateFiringFrequency;
        bool bret = rng() < rCurrentFrequency;
//        if (bret && rCurrentFrequency > 3.F / NeuronTimeDepth)
//            rCurrentFrequency *= 0.5F;
        LastTactinThisState = ntact;
        return bret;
    }
};

class DYNAMIC_LIBRARY_EXPORTED_CLASS rec_ping_pong: public IReceptors
{
    vector<float>               vr_VelocityZoneBoundary;
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

int nRewards = 0;
int nPunishments = 0;
int nRewardsTot = 0;
int nPunishmentsTot = 0;

const int RewardTrainLength = 1 /*10*/;
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
	Evaluator(enum Evaluator::type t, size_t tactbeg = 0) : IReceptors(t == reward || t == punishment ? 1 : nGoalLevels), typ(t) {}
	virtual bool bGenerateReceptorSignals(char *prec, size_t neuronstrsize) override
	{
		switch (typ) {
			case punishment: if (es.pprr_Ball->first < -0.5F) {
//								TrainCounter = RewardTrainLength;  no primary punishment
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
								 if (ntact == tactLevelFixed + NeuronTimeDepth) {
									int NewLevel = cb->Predict();
									if (NewLevel < CurrentLevel) {
										prec[neuronstrsize * NewLevel] = 1;
                                        PunishedLevel = -1;
									} else if (NewLevel == CurrentLevel)
                                        PunishedLevel = -1;
                                    else PunishedLevel = CurrentLevel;
                                    CurrentLevel = NewLevel;
								 }
				                 return true;
			default: FORI(nGoalLevels) 
						prec[_i * neuronstrsize] = 0;
                     if (ntact == tactLevelFixed + NeuronTimeDepth && PunishedLevel >= 0)
                         prec[neuronstrsize * PunishedLevel] = 1;
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

PING_PONG_ENVIRONMENT_EXPORT void SetParametersOut(int ExperimentId, size_t tactTermination, unsigned nOutputNeurons, const pugi::xml_node &xn) 
{
	nNeuronsperAction = (int)nOutputNeurons / 2;
	tactStart = atoi_s(xn.child("start_time").child_value());

	NeuronTimeDepth = atoi_s(xn.child("NeuronTimeDepth").child_value());
	float rNSpikesperNeuronTime = atof_s(xn.child("rNSpikesperNeuronTime").child_value());
	rStateFiringFrequency = rNSpikesperNeuronTime / NeuronTimeDepth;
	int LevelDuration = atoi_s(xn.child("LevelDuration").child_value());
	LevelNeuronPeriod = max(LevelDuration / NeuronTimeDepth, 1);
	minnSignificantExtraSpikes = atoi_s(xn.child("minnSignificantExtraSpikes").child_value());

}

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

int LastRegistration = -1000000;
int PostRewardCounter = 0;
int nRecognitions = 0;
int nCorr = 0;

PING_PONG_ENVIRONMENT_EXPORT int Finalize(int OriginalTerminationCode) 
{
	cout << "rew " << nRewards << " pun " << nPunishments << endl;

    double dmean, dstderr;
    avgdis(&v_err.front(), v_err.size(), dmean, dstderr);
    cout << "err: mean " << dmean << " stderr " << dstderr << endl;

	return nRewardsTot * 10000 / (nPunishmentsTot + nRewardsTot);
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
