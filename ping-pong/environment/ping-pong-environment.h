#pragma once

#include <random>
#include <iostream>

#ifdef PINGPONGENVIRONMENT_EXPORTS
#include <sg/sg.h>
#include <NetworkConfigurator.h>
#endif

#define ENVIRONMENT_STATE_SHARED_MEMORY_NAME "ping-pong-environment.sm"

#define RACKET_SIZE 0.18F

#ifdef FOR_LINUX
#define PING_PONG_ENVIRONMENT_EXPORT extern "C"
#define DYNAMIC_LIBRARY_EXPORTED_CLASS
typedef __uint64_t UNS64;
#else
#define PING_PONG_ENVIRONMENT_EXPORT __declspec(dllexport)
#define DYNAMIC_LIBRARY_EXPORTED_CLASS __declspec(dllexport)
typedef unsigned __int64 UNS64;
#endif

const unsigned nSpatialZones = 30;
const int nVelocityZones = 30;
const float rRelPosStep = RACKET_SIZE / 3;   // Racket takes 3 middle positions of the nRelPos x nRelPos grid.
const unsigned nInputs = 3 * nSpatialZones + 2 * nVelocityZones;

void UpdateWorld(std::vector<float> &vr_PhaseSpacePoint);
float rMakeBallVelocity(void);
std::pair<float, float> prr_MakeBallVelocity(bool bPunishment);

class RandomNumberGenerator
{
    std::uniform_real_distribution<> urd;
    std::mt19937_64                  mt;
public:
    RandomNumberGenerator() : urd(0., 1.) {}
    double operator()() { return urd(mt); }
    template<class T> T operator()(T max) { return (T)((*this)() * max); }
    void Randomize(void)
    {
        std::mt19937_64 mtRandomized(std::random_device{}());
        mt = mtRandomized;
    }
};

const int RewardTrainLength = 1 /*10*/;
const int RewardTrainPeriod = 2;
const int afterRewardSilence = 30;

#ifdef PINGPONGENVIRONMENT_EXPORTS
class DYNAMIC_LIBRARY_EXPORTED_CLASS GrowingStimulation: public IReceptors
{
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
            std::cout << "ping-pong -- Too many stimulating nodes\n";
            exit(-1);
        }
        maxIdleTime = atoi_s(xn.child_value("maxidletime"));
        reset();
    }
    virtual bool bGenerateSignals(unsigned *pfl, int bitoffset) override;
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
    int  ActivityTime;
};
#endif
