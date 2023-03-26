#define _USE_MATH_DEFINES

#include "EnvironmentState.hpp"

using namespace std;

extern pair<float, float> prr_BallSpeed;

const unsigned minSpotPassageTime_ms = 300;
const unsigned maxSpotPassageTime_ms = 1000;

EnvironmentState es;
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
