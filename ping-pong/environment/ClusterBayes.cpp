#include <fstream>
#include "ClusterBayes.hpp"
#include "ping-pong-environment.h"

using namespace std;

DECLARE_UNION_OPERATORS(set<unsigned>)

extern int nGoalLevels;
extern int LevelNeuronPeriod;
extern int CurrentLevel;
extern vector<float> vr_CurrentPhaseSpacePoint;
extern int ntact;
extern int LevelDuration;
extern bool bTargetState;
extern int tactStart;
extern int NeuronTimeDepth;

const unsigned aAlternativeInputBoundaries[] = {30, 60, 69, 78, 108};
int minnSignificantExtraSpikes = 10;

ClusterBayes::ClusterBayes(): vvvb_byLevels(nGoalLevels + 1), vn_Neuron(nInputs, 0), smvn_Pairs(nInputs), vd_BayesPriors(nGoalLevels + 1, 0.)
{
    FORI(nInputs)
        for (int i = _i; i < nInputs; ++i)
            smvn_Pairs(i, _i).resize(nGoalLevels + 1, 0);
}

bool ClusterBayes::bClusterFeatureIsOK(const set<unsigned> &s_) const
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

void ClusterBayes::FixReward()
{
    unsigned j;
    size_t l;
    ntotMeasurements += qvb_RecentNeuronInput.size();
    int lev;
    for (lev = 0; lev < nGoalLevels && qvb_RecentNeuronInput.size(); ++lev) {
        auto i = qvb_RecentNeuronInput.size() > LevelNeuronPeriod ? qvb_RecentNeuronInput.begin() + LevelNeuronPeriod : qvb_RecentNeuronInput.end();
        vvvb_byLevels[lev].insert(vvvb_byLevels[lev].end(), qvb_RecentNeuronInput.begin(),  i);
        auto k = qvb_RecentNeuronInput.begin();
        while (k != i) {
            FORI(nInputs)
                if ((*k)[_i]) {
                    ++smvn_Pairs(_i, _i)[lev];
                    FOR_(l, _i)
                        if ((*k)[l])
                            ++smvn_Pairs(_i, l)[lev];
                }
            ++k;
        }
        qvb_RecentNeuronInput.erase(qvb_RecentNeuronInput.begin(), i);
    }
    if (lev == nGoalLevels) {
        for (const auto &m: vvvb_byLevels[lev])
            FORI(nInputs)
                if (m[_i]) {
                    ++smvn_Pairs(_i, _i)[lev];
                    FOR_(l, _i)
                        if (m[l])
                            ++smvn_Pairs(_i, l)[lev];
                }
        vvvb_byLevels[lev].insert(vvvb_byLevels[lev].end(), qvb_RecentNeuronInput.begin(),  qvb_RecentNeuronInput.end());
    }
    qvb_RecentNeuronInput.clear();

    //FOR_(m, q_TrueLevel.size()) {
    //    vvvb_byLevels[q_TrueLevel[m]].push_back(qvb_RecentNeuronInput[m]);
    //        FORI(nInputs)
    //        if (qvb_RecentNeuronInput[m][_i]) {
    //            ++smvn_Pairs(_i, _i)[q_TrueLevel[m]];
    //                FOR_(l, _i)
    //                if (qvb_RecentNeuronInput[m][l])
    //                    ++smvn_Pairs(_i, l)[q_TrueLevel[m]];
    //            }
    //    }
    //q_TrueLevel.clear();
    //qvb_RecentNeuronInput.clear();

    mmd_AssociatedPairs.clear();
    FORI(nGoalLevels) {
        j = 0;
        for (auto o: aAlternativeInputBoundaries)
            for (; j < o; ++j) {
                int n = smvn_Pairs(j, j)[_i];
                for (l = o; l < nInputs; ++l) {
                    int n1 = smvn_Pairs(l, j)[_i];
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

int TrueCurrentLevel()
{
    if (vr_CurrentPhaseSpacePoint[2] >= 0)
        return 4;
    double dyintersection = vr_CurrentPhaseSpacePoint[1] - (vr_CurrentPhaseSpacePoint[0] + 0.5) * vr_CurrentPhaseSpacePoint[3] / vr_CurrentPhaseSpacePoint[2];
    if (dyintersection < vr_CurrentPhaseSpacePoint[4] - RACKET_SIZE / 2)
        dyintersection = vr_CurrentPhaseSpacePoint[4] - RACKET_SIZE / 2;
    else if (dyintersection > vr_CurrentPhaseSpacePoint[4] + RACKET_SIZE / 2)
        dyintersection = vr_CurrentPhaseSpacePoint[4] + RACKET_SIZE / 2;
    double dvs = -(vr_CurrentPhaseSpacePoint[0] + 0.5) * vr_CurrentPhaseSpacePoint[2] + (dyintersection - vr_CurrentPhaseSpacePoint[1]) * vr_CurrentPhaseSpacePoint[3];
    if (dvs <= 0.)
        return 4;
    double ds2 = (vr_CurrentPhaseSpacePoint[0] + 0.5) * (vr_CurrentPhaseSpacePoint[0] + 0.5) + (vr_CurrentPhaseSpacePoint[1] - dyintersection) * (vr_CurrentPhaseSpacePoint[1] - dyintersection);
    double dt = ds2 / dvs;
    return min((int)(dt / LevelDuration), 4);
}

ofstream ofsState("ping_pong_state.csv");

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

vector<int> v_err;

int ClusterBayes::Predict()
{
    if (CurrentLevel != 1000) {
        int PredictedLevel = TrueCurrentLevel();

        if (ofsState.is_open()) {
            ofsState << ntact << ',' << PredictedLevel;
            for (auto z : vr_CurrentPhaseSpacePoint)
                ofsState << ',' << z;
            ofsState << endl;
        }

        bTargetState = !PredictedLevel;

        return PredictedLevel;
    }

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
            v_err.push_back(PredictedLevel - TrueCurrentLevel());

        bTargetState = !PredictedLevel;
        return PredictedLevel;
    }
    return CurrentLevel;
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

//        q_TrueLevel.push_front(TrueCurrentLevel());

    }
}
