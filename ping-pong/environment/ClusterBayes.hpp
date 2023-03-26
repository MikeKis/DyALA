#ifndef CLUSTERBAYES_HPP
#define CLUSTERBAYES_HPP

#include <map>
#include <vector>
#include <deque>
#include <set>

#include <boost/numeric/ublas/symmetric.hpp>

class ClusterBayes
{
    std::map<std::vector<unsigned>, std::vector<double> >  mvvd_BayesModel;
    std::vector<double>                                    vd_BayesPriors;
    std::deque<std::vector<bool> >                         qvb_RecentNeuronInput;

//    deque<int>                                 q_TrueLevel;   // used only in this version of surragate Bayes neuron mechanism.

    std::deque<std::vector<bool> >                         qvb_NeuronInput;
    std::vector<std::vector<std::vector<bool> > >          vvvb_byLevels;
    std::vector<int>                                       vn_Neuron;
    size_t                                                 ntotMeasurements = 0;
    boost::numeric::ublas::symmetric_matrix<std::vector<unsigned> > smvn_Pairs;
    std::map<unsigned,std::map<unsigned, double> >         mmd_AssociatedPairs;
    size_t nLevelMeasurements(int Level) const {return vvvb_byLevels[Level].size();}
public:
    ClusterBayes();
    void AddNewInput(const std::vector<bool> &vb_Spikes);
    int Predict();
    void FixReward();
    bool bClusterFeatureIsOK(const std::set<unsigned> &s_) const;
};

#endif // CLUSTERBAYES_HPP
