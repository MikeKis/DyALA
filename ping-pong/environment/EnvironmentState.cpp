#define _USE_MATH_DEFINES

#include "EnvironmentState.hpp"

using namespace std;

extern int nNeuronsperAction;

extern pair<float, float> prr_BallSpeed;

const unsigned minSpotPassageTime_ms = 200;
const unsigned maxSpotPassageTime_ms = 700;

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

pair<float, float> prr_MakeBallVelocity(bool bPunishment)
{
    float rBallVelocity = rMakeBallVelocity();
    float rBallMovementDirection;
    pair<float, float> prr_ret;
    do {
        rBallMovementDirection = (float)rng(2 * M_PI);
        prr_ret.first = rBallVelocity * sin(rBallMovementDirection);
    } while (prr_ret.first < 1. / (2 * maxSpotPassageTime_ms) || !bPunishment && prr_ret.first < 0.);
    prr_ret.second = rBallVelocity * cos(rBallMovementDirection);
    return prr_ret;
}

void EnvironmentState::ResetBall(bool bPunishment)
{
    if (bPunishment) {
        es.pprr_Ball->first = 0.F;
        es.pprr_Ball->second = (float)(-0.5 + rng());
    }
    prr_BallSpeed = prr_MakeBallVelocity(bPunishment);

    // It is a special temporary regime for intermediate goal mechanism training data creation.

    if (!nNeuronsperAction)
        ResetRacket();

}

void EnvironmentState::ResetRacket() {*prRacket = (float)(-0.5 + RACKET_SIZE / 2 + rng(1. - RACKET_SIZE));}
