/*
Copyright (C) 2022 Mikhail Kiselev

Emulates ping-pong game.

*/

#define _USE_MATH_DEFINES

#include <math.h>
#include <random>
#include <fstream>
#include <sstream>
#include <memory>

#include <ArNI/../helpers/DatTab/RF/RF.h>

#include "EnvironmentState.hpp"

//#define EXACT_STATE_EVALUATION

using namespace std;
using namespace boost::interprocess;

extern EnvironmentState es;
extern RandomNumberGenerator rng;

extern pair<float, float> prr_BallSpeed;

extern bool b_forVerifier_Reward;

extern int nerrRF;
extern int ncorrRF;

int nGoalLevels;

int ntact = -1;
int tactStart;

int nRewards = 0;
int nPunishments = 0;

int nSecondaryRewards = 0;
int nSecondaryPunishments = 0;

int nRewardsTot = 0;
int nPunishmentsTot = 0;
int InputBlockCounter = 0;

enum evaluation_type
{
    reward,
    punishment,
};

class DYNAMIC_LIBRARY_EXPORTED_CLASS Evaluator: public IReceptors
{
    enum evaluation_type typ;
    int                  TrainCounter = 0;
    int                  PeriodCounter;

protected:
    virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const override
    {
            vstr_Meanings.resize(1);
            vstr_Meanings.front() = typ == reward ? "REW" : "PUN";
    }
public:
    Evaluator(enum evaluation_type t) : typ(t) {}
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override
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
                             }
                             break;
            default:         return true;
        }
        *pfl = 0;
        if (TrainCounter && !--PeriodCounter) {
            *pfl = 1;
            PeriodCounter = RewardTrainPeriod;
            --TrainCounter;
            InputBlockCounter = afterRewardSilence;
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

vector<float> vr_CurrentPhaseSpacePoint(5);   // It does not include relative positions!

enum phase_space_coordinates
{
    BallX = 0,
    BallY = 1,
    BallVX = 2,
    BallVY = 3,
    RacketY = 4,
    CloseZone = 5,
    phase_space_dimension
};

bool bCurrentStateOK = false;
float rIntensityFactor = 0.7;

deque<pair<int, int> > aqpind_[5];   // for manual state classification
deque<pair<int, int> > qpind_all;   // for RF state classification
const int SpikeQueueDepth_tacts = 10;
const int AllSpikeQueueDepth_tacts = 30;
const int ModelRecreationPeriod_tacts = 210000;
DatTab dt(nInputs);
RF rf;
RFPar rfp;
std::unique_ptr<RFRes> uprfr;

vector<unique_ptr<DoGEncoder> > vupdog_(phase_space_dimension);

int State(double dxball, double dyball, double dvxball, double dvyball, double dyracket)
{
    if (dxball >= 0)
        return -1;
    if (dvxball >= 0)
        return -1;
//    if (dvx > 0) {  // retain it!
//        dy += (0.5 - dx) * dvy / dvx;
//        dx = 0.5;
//        dvx = -dvx;
//    }
    double ddx = dxball + 0.5;
    double dyint = dyball - ddx * dvyball / dvxball + 0.5;
    dyint -= floor(dyint / 2) * 2;
    double d = dyracket + 0.5;
    if (d - RACKET_SIZE / 2 < dyint && dyint < d + RACKET_SIZE / 2)
        return 0;
    d = 2 - d;
    return d - RACKET_SIZE / 2 < dyint && dyint < d + RACKET_SIZE / 2 ? 0 : 1;
}

class DYNAMIC_LIBRARY_EXPORTED_CLASS rec_ping_pong: public IReceptors
{
protected:
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override
    {

        ++ntact;
        UpdateWorld(vr_CurrentPhaseSpacePoint);

#ifdef LOG_STATE
        static ofstream ofsState("ping_pong_state.csv");
        if (ofsState.is_open()) {
            ofsState << ntact;
            for (auto z: vr_CurrentPhaseSpacePoint)
                ofsState << ',' << z;
            ofsState << endl;
        }
#endif

        fill(pfl, pfl + AfferentSpikeBufferSizeDW(GetNReceptors()), 0);

        for (auto &y: aqpind_)
            while (y.size() && y.front().first < ntact - SpikeQueueDepth_tacts)
                y.pop_front();
        while (qpind_all.size() && qpind_all.front().first < ntact - AllSpikeQueueDepth_tacts)
            qpind_all.pop_front();

        if (!InputBlockCounter) {
            vector<bool> vb_Spatial(nSpatialZones);
            vector<bool> vb_Velocity(nVelocityZones);
            vector<bool> vb_CloseZone(nRelPos * nRelPos);
            BitMaskAccess bma;
            (*vupdog_[BallX])(vr_CurrentPhaseSpacePoint[BallX], vb_Spatial);
            int z = 0;
            int y = 0;
            for (auto b: vb_Spatial) {
                if (b) {
                    pfl |= bma;
                    aqpind_[BallX].push_back(pair<int, int>(ntact, z));
                    qpind_all.push_back(pair<int, int>(ntact, y));
                }
                ++bma;
                ++z;
                ++y;
            }
            (*vupdog_[BallY])(vr_CurrentPhaseSpacePoint[BallY], vb_Spatial);
            z = 0;
            for (auto b: vb_Spatial) {
                if (b) {
                    pfl |= bma;
                    aqpind_[BallY].push_back(pair<int, int>(ntact, z));
                    qpind_all.push_back(pair<int, int>(ntact, y));
                }
                ++bma;
                ++z;
                ++y;
            }
            (*vupdog_[BallVX])(vr_CurrentPhaseSpacePoint[BallVX], vb_Velocity);
            z = 0;
            for (auto b: vb_Velocity) {
                if (b) {
                    pfl |= bma;
                    aqpind_[BallVX].push_back(pair<int, int>(ntact, z));
                    qpind_all.push_back(pair<int, int>(ntact, y));
                }
                ++bma;
                ++z;
                ++y;
            }
            (*vupdog_[BallVY])(vr_CurrentPhaseSpacePoint[BallVY], vb_Velocity);
            z = 0;
            for (auto b: vb_Velocity) {
                if (b) {
                    pfl |= bma;
                    aqpind_[BallVY].push_back(pair<int, int>(ntact, z));
                    qpind_all.push_back(pair<int, int>(ntact, y));
                }
                ++bma;
                ++z;
                ++y;
            }
            (*vupdog_[RacketY])(vr_CurrentPhaseSpacePoint[RacketY], vb_Spatial);
            z = 0;
            for (auto b: vb_Spatial) {
                if (b) {
                    pfl |= bma;
                    aqpind_[RacketY].push_back(pair<int, int>(ntact, z));
                    qpind_all.push_back(pair<int, int>(ntact, y));
                }
                ++bma;
                ++z;
                ++y;
            }
            (*vupdog_[CloseZone])(VECTOR<float>{vr_CurrentPhaseSpacePoint[BallX] + 0.5F, vr_CurrentPhaseSpacePoint[BallY] - vr_CurrentPhaseSpacePoint[RacketY]}, vb_CloseZone);
            for (auto b: vb_CloseZone) {
                if (b) {
                    pfl |= bma;
                    qpind_all.push_back(pair<int, int>(ntact, y));
                }
                ++bma;
                ++y;
            }
        } else --InputBlockCounter;

//        if (!(ntact % AllSpikeQueueDepth_tacts)) {
//            int CurrentStateTrue = State(vr_CurrentPhaseSpacePoint[0], vr_CurrentPhaseSpacePoint[1], vr_CurrentPhaseSpacePoint[2], vr_CurrentPhaseSpacePoint[3], vr_CurrentPhaseSpacePoint[4]);
//            if (CurrentStateTrue >= 0) {
//                vector<int> vn_spikes(nInputs, 0);
//                for (auto i: qpind_all)
//                    ++vn_spikes[i.second];
//                dt.Append(CurrentStateTrue, vn_spikes, ntact);
//            }
//            if (dt.NRecs() > 30 && !(ntact % ModelRecreationPeriod_tacts)) {
//                rfp.DominatingClass = 1;
//                uprfr.reset(dynamic_cast<RFRes *>(rf.pmlmCreateModel(&rfp, dt)));
//            }
//        }

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
    rec_ping_pong()
    {
        const int nVelocitySamples = 9000;
        vector<float> vr_vxsamples(nVelocitySamples);
        vector<float> vr_vysamples(nVelocitySamples);
        for (int i = 0; i < nVelocitySamples; ++i) {
            auto prr_v = prr_MakeBallVelocity(false);
            vr_vxsamples[i] = rng() < 0.5 ? -prr_v.first : prr_v.first;
            vr_vysamples[i] = prr_v.second;
        }
        vupdog_[BallX].reset(new DoGEncoder(-0.5, 0.5, 30, rIntensityFactor));
        vupdog_[BallY].reset(new DoGEncoder(-0.5, 0.5, 30, rIntensityFactor));
        vupdog_[BallVX].reset(new DoGEncoder(vr_vxsamples, nVelocityZones, rIntensityFactor));
        vupdog_[BallVY].reset(new DoGEncoder(vr_vysamples, nVelocityZones, rIntensityFactor));
        vupdog_[RacketY].reset(new DoGEncoder(-0.5, 0.5, 30, rIntensityFactor));
        vupdog_[CloseZone].reset(new DoGEncoder({{{0., rRelPosStep * nRelPos}, nRelPos}, {{-rRelPosStep * nRelPos / 2, rRelPosStep * nRelPos / 2}, nRelPos}}, rIntensityFactor));

    }
    virtual void Randomize(void) override {rng.Randomize();}

    // NOT WORKABLE !!

    virtual void SaveStatus(Serializer &ser) const override
    {
        IReceptors::SaveStatus(ser);
        ser << *es.pprr_Ball;
        ser << *es.prRacket;
        ser << prr_BallSpeed;
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
        ser >> rng;
    }

};

#ifndef BIGREALNUMBER
#define BIGREALNUMBER 1.e10F
#endif

class DYNAMIC_LIBRARY_EXPORTED_CLASS Actions: public IReceptors
{
    float rPastRY = BIGREALNUMBER;
    float rStep;
protected:
    virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const override
    {
        vstr_Meanings.resize(2);
        vstr_Meanings[0] = "ActDown";
        vstr_Meanings[1] = "ActUp";
    }
public:
    Actions(float rstep): rStep(rstep) {}
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override;
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

int StateChangeDelay = 0;

int StatefromSpikes()
{
    auto drestore_parameter = [&](int dim)
    {
        auto i = aqpind_[dim].rbegin();
        int n = 0;
        int lastind = -1;
        double dret = 0;
        while (i != aqpind_[dim].rend() /* && (lastind == -1 || abs(i->second - lastind) < 3)*/) {
            lastind = i->second;
            dret += vupdog_[dim]->vvg_Zone.front()[lastind].dCenter;
            ++n;
            ++i;
        }
        return dret / n;
    };
    double dx, dy, dvx, dvy, dry;

    dx = drestore_parameter(BallX);
    dy = drestore_parameter(BallY);
    dvx = drestore_parameter(BallVY);
    dvy = drestore_parameter(BallVY);
    dry = drestore_parameter(RacketY);

    return State(dx, dy, dvx, dvy, dry);
}

float rRewardSwitchThreshold = 0, rPunishmentSwitchThreshold = 0;

int StatefromRF(int CurrentState)
{
    if (!uprfr.get())
        return -1;
    vector<float> vn_spikes(nInputs, 0.F);
    for (auto i: qpind_all)
        ++vn_spikes[i.second];
    auto pr_pred = rf.pr_Apply(uprfr.get(), vn_spikes);
    return !pr_pred.first && pr_pred.second > 1 - rRewardSwitchThreshold || pr_pred.first == 1 && pr_pred.second > 1 - rPunishmentSwitchThreshold ? pr_pred.first : CurrentState;
}

bool bState(bool bRewardRequested)
{
    static int tactLastStateChange = -1000000;
    static int tactLastStateChangeTrue = -1000000;
    static bool bLastChangetoGood;
    static bool bLastChangetoGoodTrue;
    static int FormerState = -1;
    static int FormerStateTrue = -1;
    static int CurrentStateTrue, CurrentState = -1;
    if (!bRewardRequested) {   // punishment is requested first

        if (InputBlockCounter || any_of(aqpind_, aqpind_ + sizeof(aqpind_) / sizeof(aqpind_[0]), [](const deque<pair<int, int> > &qpind_){return qpind_.empty();})) {
            FormerState = FormerStateTrue = -1;
            tactLastStateChange = -1000000;
            return false;
        }

        CurrentStateTrue = State(vr_CurrentPhaseSpacePoint[0], vr_CurrentPhaseSpacePoint[1], vr_CurrentPhaseSpacePoint[2], vr_CurrentPhaseSpacePoint[3], vr_CurrentPhaseSpacePoint[4]);
        CurrentState = StatefromRF(CurrentState) /* StatefromSpikes() */;

        if (CurrentStateTrue >= 0 && CurrentState >= 0) {
            if (CurrentStateTrue == CurrentState)
                ++ncorrRF;
            else ++nerrRF;
        }

        if (CurrentStateTrue == -1) {
            FormerStateTrue = -1;
            tactLastStateChangeTrue = -1000000;
#ifdef EXACT_STATE_EVALUATION
            return false;
#endif
        }
        if (CurrentState == -1) {
            FormerState = -1;
            tactLastStateChange = -1000000;
#ifndef EXACT_STATE_EVALUATION
            return false;
#endif
        }
        if (CurrentStateTrue != -1 && FormerStateTrue != -1 && FormerStateTrue != CurrentStateTrue) {
            if (ntact - tactLastStateChangeTrue > StateChangeDelay) {
                tactLastStateChangeTrue = ntact;
                bLastChangetoGoodTrue = !CurrentStateTrue;
            } else tactLastStateChangeTrue = -1000000;
        }
        if (CurrentState != -1 && FormerState != -1 && FormerState != CurrentState) {
            if (ntact - tactLastStateChange > StateChangeDelay) {
                tactLastStateChange = ntact;
                bLastChangetoGood = !CurrentState;
            } else tactLastStateChange = -1000000;
        }
        FormerState = CurrentState;
        FormerStateTrue = CurrentStateTrue;

#ifdef EXACT_STATE_EVALUATION
        bCurrentStateOK = !CurrentStateTrue;
#else
        bCurrentStateOK = !CurrentState;
#endif

#ifdef LOG_STATE
        static ofstream ofsState("ping_pong_state_rest.csv");
        static bool bHeaderWritten = false;
        if (!bHeaderWritten) {
            ofsState << "tact,x,y,vx,vy,yr,state,statep,tactchage,tactchangep,goodchange,goodchagep,dec,decp\n";
            bHeaderWritten = true;
        }
        ofsState << ntact
                 << ','
                 << vr_CurrentPhaseSpacePoint[0]
                 << ','
                 << vr_CurrentPhaseSpacePoint[1]
                 << ','
                 << vr_CurrentPhaseSpacePoint[2]
                 << ','
                 << vr_CurrentPhaseSpacePoint[3]
                 << ','
                 << vr_CurrentPhaseSpacePoint[4]
                 << ','
                 << CurrentStateTrue
                 << ','
                 << CurrentState
                 << ','
                 << tactLastStateChangeTrue
                 << ','
                 << tactLastStateChange
                 << ','
                 << (bLastChangetoGoodTrue ? '1' : '0')
                 << ','
                 << (bLastChangetoGood ? '1' : '0')
                 << ','
                 << (ntact - tactLastStateChangeTrue == StateChangeDelay ? (bLastChangetoGoodTrue ? "1" : "-1") : "0")
                 << ','
                 << (ntact - tactLastStateChange == StateChangeDelay ? (bLastChangetoGood ? "1" : "-1") : "0")
                 << endl;
#endif

    }
#ifdef EXACT_STATE_EVALUATION
    return CurrentStateTrue != -1 && ntact - tactLastStateChangeTrue == StateChangeDelay && bLastChangetoGoodTrue == bRewardRequested;
#else
    return CurrentState != -1 && ntact - tactLastStateChange == StateChangeDelay && bLastChangetoGood == bRewardRequested;
#endif
}

class DYNAMIC_LIBRARY_EXPORTED_CLASS SecondaryEvaluator: public IReceptors
{
    bool bReward;
protected:
    virtual void GetMeanings(VECTOR<STRING>& vstr_Meanings) const override
    {
        vstr_Meanings.resize(1);
        vstr_Meanings.front() = bReward ? "SECREW" : "SECPUN";
    }
public:
    SecondaryEvaluator(bool brew) : bReward(brew)
    {
        if (!brew) {
            std::ifstream ifsRF("ping-pong.RF.bin");
            uprfr.reset(dynamic_cast<RFRes *>(rf.pmlmLoadModel(ifsRF)));
        }
    }
    virtual bool bGenerateSignals(unsigned* pfl, int bitoffset) override
    {
        *pfl = bState(bReward);

        if (*pfl) {
            if (bReward)
                ++nSecondaryRewards;
            else ++nSecondaryPunishments;

        }

        return true;
    }
    virtual void Randomize(void) override {};
    virtual void SaveStatus(Serializer& ser) const override
    {
        IReceptors::SaveStatus(ser);
        ser << bReward;
    }
    virtual ~SecondaryEvaluator() = default;
    void LoadStatus(Serializer& ser)
    {
        IReceptors::LoadStatus(ser);
        ser >> bReward;
    }
};

GrowingStimulation *pgsG = NULL;

bool Actions::bGenerateSignals(unsigned *pfl, int bitoffset)
{
    *pfl = 0;
    if (rPastRY == BIGREALNUMBER)
        rPastRY = vr_CurrentPhaseSpacePoint[4];
    else if (vr_CurrentPhaseSpacePoint[4] < rPastRY - rStep) {
        *pfl = 1;
        rPastRY = vr_CurrentPhaseSpacePoint[4];
        pgsG->reset();
    } else if (vr_CurrentPhaseSpacePoint[4] > rPastRY + rStep) {
        *pfl = 2;
        rPastRY = vr_CurrentPhaseSpacePoint[4];
        pgsG->reset();
    }
    return true;
}

RECEPTORS_SET_PARAMETERS(pchMyReceptorSectionName, nReceptors, xn)
{
    if (!strcmp(pchMyReceptorSectionName, "Punishment")) {
        if (nReceptors < 0)
            nReceptors = 1;
        else if (nReceptors != 1)
            throw std::runtime_error("Punishment - wrong input node count");
        return new Evaluator(punishment);
    } else if (!strcmp(pchMyReceptorSectionName, "Reward")) {
        if (nReceptors < 0)
            nReceptors = 1;
        else if (nReceptors != 1)
            throw std::runtime_error("Reward - wrong input node count");
        return new Evaluator(reward);
    } else if (!strcmp(pchMyReceptorSectionName, "R")) {
        if (nReceptors < 0)
            nReceptors = nInputs;
        else if (nReceptors != nInputs)
            throw std::runtime_error("R - wrong input node count");
        return new rec_ping_pong;
    } else if (!strcmp(pchMyReceptorSectionName, "Actions")) {
        if (nReceptors < 0)
            nReceptors = 2;
        else if (nReceptors != 2)
            throw std::runtime_error("Actions - wrong input node count");
        return new Actions(atof_s(xn.child_value("step")));
    } else if (!strcmp(pchMyReceptorSectionName, "SECPUN")) {
        if (nReceptors < 0)
            nReceptors = 1;
        else if (nReceptors != 1)
            throw std::runtime_error("SECPUN - wrong input node count");
        StateChangeDelay = atoi_s(xn.child_value("state_change_delay"));
        rPunishmentSwitchThreshold = atof_s(xn.child_value("punishment_switch_threshold"));
        rRewardSwitchThreshold = atof_s(xn.child_value("reward_switch_threshold"));
        return new SecondaryEvaluator(false);
    } else if (!strcmp(pchMyReceptorSectionName, "SECREW")) {
        if (nReceptors < 0)
            nReceptors = 1;
        else if (nReceptors != 1)
            throw std::runtime_error("SECREW - wrong input node count");
        return new SecondaryEvaluator(true);
    } else if (!strcmp(pchMyReceptorSectionName, "GrowingStimulation"))
        return pgsG = new GrowingStimulation(nReceptors, xn);
    else {
        cout << "ping-pong -- unknown section name: " << pchMyReceptorSectionName << endl;
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
		case 1: peva = new Evaluator(punishment);
				peva->LoadStatus(ser);
				return peva;
		case 2: peva = new Evaluator(reward);
				peva->LoadStatus(ser);
				return peva;
//		case 3: papG = new AdaptivePoisson(true);
//				papG->LoadStatus(ser);
//				return papG;
		default: cout << "Too many calls of LoadStatus\n";
				exit(-1);
	}
}

bool GrowingStimulation::bGenerateSignals(unsigned *pfl, int bitoffset)
{
    size_t *pfl64 = (size_t *)pfl;   // SAFE BECAUSE OF 64 BIT RECEPTOR ALIGNMENT!
    ++ActivityTime;
    if (ActivityTime < 0 || ntact >= tactStart)
        *pfl64 = 0;
    else {
        if (!ActivityTime)
            bh = rng() > 0.5;
        if (ActivityTime < GetNReceptors() / 2) {
            *pfl64 = (1 << (ActivityTime + 1)) - 1;
            if (bh)
                *pfl64 <<= GetNReceptors() / 2;
        } else {
            *pfl64 = 0;
            ActivityTime = -1;
        }
    }
    return true;
}
