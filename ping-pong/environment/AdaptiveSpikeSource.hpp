#ifndef ADAPTIVESPIKESOURCE_HPP
#define ADAPTIVESPIKESOURCE_HPP


const float rStateFiringFrequency = 0.75F;

class AdaptiveSpikeSource
{
//    int LastTactinThisState;
    float rCurrentFrequency;
public:
    AdaptiveSpikeSource(): rCurrentFrequency(rStateFiringFrequency) {}
    bool bFire();
};

#endif // ADAPTIVESPIKESOURCE_HPP
