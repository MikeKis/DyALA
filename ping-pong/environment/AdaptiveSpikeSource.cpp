#include "ping-pong-environment.h"
#include "AdaptiveSpikeSource.hpp"

extern int ntact;
extern RandomNumberGenerator rng;

float rStateFiringFrequency /* = 5.F / NeuronTimeDepth */;

bool AdaptiveSpikeSource::bFire()
{
    if (LastTactinThisState < ntact - 1)
        rCurrentFrequency = rStateFiringFrequency;
    bool bret = rng() < rCurrentFrequency;
//        if (bret && rCurrentFrequency > 3.F / NeuronTimeDepth)
//            rCurrentFrequency *= 0.5F;
    LastTactinThisState = ntact;
    return bret;
}
