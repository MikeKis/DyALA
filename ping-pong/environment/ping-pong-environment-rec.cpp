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

#include "EnvironmentState.hpp"

using namespace std;
using namespace boost::interprocess;

extern EnvironmentState es;
extern RandomNumberGenerator rng;

extern pair<float, float> prr_BallSpeed;

extern bool b_forVerifier_Reward;

int nGoalLevels;

int ntact = -1;
int tactStart;

int nRewards = 0;
int nPunishments = 0;
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

class DYNAMIC_LIBRARY_EXPORTED_CLASS rec_ping_pong: public IReceptors
{
    vector<unique_ptr<DoGEncoder> > vupdog_;
protected:
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override
    {

        ++ntact;
        UpdateWorld(vr_CurrentPhaseSpacePoint);

        static ofstream ofsState("ping_pong_state.csv");
        if (ofsState.is_open()) {
            ofsState << ntact;
            for (auto z: vr_CurrentPhaseSpacePoint)
                ofsState << ',' << z;
            ofsState << endl;
        }

        fill(pfl, pfl + AfferentSpikeBufferSizeDW(GetNReceptors()), 0);

        if (!InputBlockCounter) {
            vector<bool> vb_Spatial(nSpatialZones);
            vector<bool> vb_Velocity(nVelocityZones);
            vector<bool> vb_CloseZone(nRelPos * nRelPos);
            BitMaskAccess bma;
            (*vupdog_[BallX])(vr_CurrentPhaseSpacePoint[BallX], vb_Spatial);
            for (auto b: vb_Spatial) {
                if (b)
                    pfl |= bma;
                ++bma;
            }
            (*vupdog_[BallY])(vr_CurrentPhaseSpacePoint[BallY], vb_Spatial);
            for (auto b: vb_Spatial) {
                if (b)
                    pfl |= bma;
                ++bma;
            }
            (*vupdog_[BallVX])(vr_CurrentPhaseSpacePoint[BallVX], vb_Velocity);
            for (auto b: vb_Velocity) {
                if (b)
                    pfl |= bma;
                ++bma;
            }
            (*vupdog_[BallVY])(vr_CurrentPhaseSpacePoint[BallVY], vb_Velocity);
            for (auto b: vb_Velocity) {
                if (b)
                    pfl |= bma;
                ++bma;
            }
            (*vupdog_[RacketY])(vr_CurrentPhaseSpacePoint[RacketY], vb_Spatial);
            for (auto b: vb_Spatial) {
                if (b)
                    pfl |= bma;
                ++bma;
            }
            (*vupdog_[CloseZone])(VECTOR<float>{vr_CurrentPhaseSpacePoint[BallX] + 0.5F, vr_CurrentPhaseSpacePoint[BallY] - vr_CurrentPhaseSpacePoint[RacketY]}, vb_CloseZone);
            for (auto b: vb_CloseZone) {
                if (b)
                    pfl |= bma;
                ++bma;
            }
        } else --InputBlockCounter;

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
    rec_ping_pong(): vupdog_(phase_space_dimension)
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

bool bState(bool bRewardRequested)
{
    static vector<float> vr_VelocityZoneMedian;
    static int tactLastStateChange = -1000000;
    static bool bLastChangetoGood;
    static int FormerState = -1;
    if (vr_VelocityZoneMedian.empty()) {
        vector<float> vr_samples(9000, 0.F);
        for (auto &i: vr_samples) {
            float rBallVelocity = rMakeBallVelocity();
            float rBallMovementDirection = rng() * PI / 2;
            i = rBallVelocity * sin(rBallMovementDirection);
        }
        sort(vr_samples.begin(), vr_samples.end());
        vr_VelocityZoneMedian.resize((nVelocityZones - 1) / 2);
        FORI(vr_VelocityZoneMedian.size())
            vr_VelocityZoneMedian[_i] = vr_samples[(_i + 1) * 2 * vr_samples.size() / 9];
    }
    if (!bRewardRequested) {   // punishment is requested first
        double dx, dy, dvx, dvy, dry;
/*        if (curindxBall < 0 || curindyBall < 0 || curindvxBall < 0 || curindvyBall < 0 || curindRacket < 0) {
            FormerState = -1;
            tactLastStateChange = -1000000;
            return false;
        }
        dx = -0.5 + (0.5 + curindxBall) / nSpatialZones; */
                    dx = vr_CurrentPhaseSpacePoint[0];
        if (dx >= 0) {
            FormerState = -1;
            tactLastStateChange = -1000000;
            return false;
        }
//      dy = -0.5 + (0.5 + curindyBall) / nSpatialZones;
                    dy = vr_CurrentPhaseSpacePoint[1];
//      dvx = curindvxBall == nVelocityZones / 2 ? 0. : curindvxBall < nVelocityZones / 2 ? -vr_VelocityZoneMedian[nVelocityZones / 2 - curindvxBall - 1] : vr_VelocityZoneMedian[curindvxBall - nVelocityZones / 2 - 1];
                    dvx = vr_CurrentPhaseSpacePoint[2];
        if (dvx >= 0) {
            FormerState = -1;
            tactLastStateChange = -1000000;
            return false;
        }
//        dvy = curindvyBall == nVelocityZones / 2 ? 0. : curindvyBall < nVelocityZones / 2 ? -vr_VelocityZoneMedian[nVelocityZones / 2 - curindvyBall - 1] : vr_VelocityZoneMedian[curindvyBall - nVelocityZones / 2 - 1];
                    dvy = vr_CurrentPhaseSpacePoint[3];
//        dry = -0.5 + (0.5 + curindRacket) / nSpatialZones;
                    dry = vr_CurrentPhaseSpacePoint[4];
        if (dvx > 0) {  // retain it!
            dy += (0.5 - dx) * dvy / dvx;
            dx = 0.5;
            dvx = -dvx;
        }
        double ddx = dx + 0.5;
        double dyint = dy - ddx * dvy / dvx + 0.5;
        dyint -= floor(dyint / 2) * 2;
        dry += 0.5;
        int CurrentState;
        if (dry - RACKET_SIZE / 2 < dyint && dyint < dry + RACKET_SIZE / 2)
            CurrentState = 0;
        else {
            dry = 2 - dry;
            CurrentState = dry - RACKET_SIZE / 2 < dyint && dyint < dry + RACKET_SIZE / 2 ? 0 : 1;
        }
        if (FormerState != -1 && FormerState != CurrentState) {
            if (ntact - tactLastStateChange > StateChangeDelay) {
                tactLastStateChange = ntact;
                bLastChangetoGood = !CurrentState;
            } else tactLastStateChange = -1000000;
        }
        FormerState = CurrentState;
        bCurrentStateOK = !CurrentState;
    }
    return ntact - tactLastStateChange == StateChangeDelay && bLastChangetoGood == bRewardRequested;
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
    SecondaryEvaluator(bool brew) : bReward(brew) {}
    virtual bool bGenerateSignals(unsigned* pfl, int bitoffset) override
    {
        *pfl = bState(bReward);
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

class DYNAMIC_LIBRARY_EXPORTED_CLASS GrowingStimulation: public IReceptors
{
    int  ActivityTime;
    int  maxIdleTime;
    bool bh;
protected:
    virtual void GetMeanings(VECTOR<STRING> &vstr_Meanings) const override
    {
        vstr_Meanings.resize(GetNReceptors());
        FORI(GetNReceptors())
            vstr_Meanings[_i] = "sti" + str(_i);
    }
public:
    GrowingStimulation(int nrec, const pugi::xml_node &xn)
    {
        if (nrec / 2 > 32) {
            cout << "ping-pong -- Too many stimulating nodes\n";
            exit(-1);
        }
        maxIdleTime = atoi_s(xn.child_value("maxidletime"));
        reset();
    }
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override
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
    virtual void Randomize(void) override {};
    virtual void SaveStatus(Serializer &ser) const override
    {
        IReceptors::SaveStatus(ser);
    }
    virtual ~GrowingStimulation() = default;
    void LoadStatus(Serializer &ser)
    {
        IReceptors::LoadStatus(ser);
    }
    void reset(){ActivityTime = -maxIdleTime;}
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
