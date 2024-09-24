//#define USE_CPU
#define EMULATOR "ArNIGPU"
#define RECOGNIZER ""

//------------------------------------------------------------------------------------------------
// $$$ COMMON DEFINES
#define MAXEMULATIONTIME 2000000
#define NINPUTS 133
#define NTESTS 3
#define SECOND_TEST_THRESHOLD 0x7fffffff   // DO NOT USE PROTOCOL POSTPROCESSOR !!
//------------------------------------------------------------------------------------------------

NetworkTestingParams ntpG{EMULATOR, RECOGNIZER, NTESTS, SECOND_TEST_THRESHOLD, NoStructuralCriterion, STARTING_EXPERIMENT_ID};   // It is because experiment ID space for GA is 1000-2000. The original experminet ID begins from 0. For local execution it is incremented
// by 1000. This incremented value is passed through -e option. But in cluster management code this task is registered by the original ID.
// Logics for CD are different!!

#ifndef ARNI_CD
//------------------------------------------------------------------------------------------------
// $$$ DEFINE CHROMOSOME

#define THRESHOLDWINC_MIN 0.0001
#define THRESHOLDWINC_MAX 1.

#define LTHRESHOLDINC_MIN 1.
#define LTHRESHOLDINC_MAX 300.

CHROMOSOME_DEFINITION
{
    INITIALIZE_GENE(THRESHOLDWINC, rng()),
    INITIALIZE_GENE(ACTIONSTEP, rng.logUniform(0.001, 0.1)),   // racket sfift causing ACTIONS firing
    INITIALIZE_GENE(MAXIDLETIME, rng.logUniform(300, 30000)),  // no-action time after which external stimalation is switched on

    // secondary evaluation neuron parameters

    INITIALIZE_GENE(SECCHARTIME, rng.logUniform(1, 30)),
    INITIALIZE_GENE(SECTHRESHOLDDECAYPERIOD, rng.logUniform(1, 3000)),
    INITIALIZE_GENE(SECTHRESHOLDINC, rng.logUniform(1, 300)),
    INITIALIZE_GENE(MINPOTSEC, rng.logUniform(0.1, 300.)),

    INITIALIZE_GENE(NNEURONSPERACTION, rng.logUniform(1, 30)),
    INITIALIZE_GENE(LCHARTIME, rng.logUniform(1, 30)),
    INITIALIZE_GENE(LTHRESHOLDDECAYPERIOD, rng.logUniform(1, 3000)),
    INITIALIZE_GENE(LTHRESHOLDINC, rng()),
    INITIALIZE_GENE(PLASTICITYTIME, rng.logUniform(3, 100)),
    INITIALIZE_GENE(MAXWEIGHT, rng.logUniform(0.3, 3.)),
    INITIALIZE_GENE(MINWEIGHT, rng.logUniform(0.003, 3.)),
    INITIALIZE_GENE(STABILITYRESOURCECHANGERATIO, rng.logUniform(0.03, 30)),
    INITIALIZE_GENE(NSILENTSYNAPSES, rng.logUniform(1, 1000)),
    INITIALIZE_GENE(MIGRATION, rng()),
    INITIALIZE_GENE(GOODSTATEINERTIA, rng.logUniform(1, 30)),
    INITIALIZE_GENE(STATECHANGEDELAY, rng.logUniform(1, 30)),
    INITIALIZE_GENE(POSWSEC, rng(10.)),
    INITIALIZE_GENE(LATERALBLOCKINGTIME, rng.logUniform(1, 100)),
    INITIALIZE_GENE(RELPUNEFFECT, rng.logUniform(0.03, 30)),
    INITIALIZE_GENE(ACTIONEVALUATIONTIME, rng.logUniform(3, 100))
};
#define NO_GENE_DEPENDENCIES
//------------------------------------------------------------------------------------------------

#ifndef FOR_LINUX
#pragma warning(disable:4244)
#pragma warning(disable:4305)
#endif

// Check that the given gene combination makes sense

bool Chromosome::bNoFatalDefect(void) const 
{
	return true;
}

int StructuralCriterion(int)
{
    throw std::logic_error("StructuralCriterion used!!!");
    return 0;
}   // not used
#else   // ARNI_CD not defined
unsigned Experiment = 3000;
//------------------------------------------------------------------------------------------------
// $$$ SET CUSTOM PARAMETER BOUNDARIES
bool bBadHyperparameters(const vector<int> &v_point, const vector<dimension> &vdim_) {return V("NLEARNINGGROUPS") < 1;}
//------------------------------------------------------------------------------------------------
#endif   // ARNI_CD defined

EVALUATOR
{
	int a, SelectedDevice;
#ifdef ARNI_CD
    if (bBadHyperparameters(v_point, vdim_))
            return 0;
#endif

// Create network configuration file and storing its contents in log file

//------------------------------------------------------------------------------------------------
// $$$ SPECIFY NNP FILE
    SETNNP(".." FOLDER_SEPARATOR "Experiments" FOLDER_SEPARATOR "ManualStateClassification.nnp");
//------------------------------------------------------------------------------------------------
// $$$ SUBSTITUTE ALL VALUES

    SUBVALUE("NNEU", nint(V(NNEURONSPERACTION)) * 2);
    SUBVALUE("INIRESOURCE0", rWeight0Resource(-V(MINWEIGHT), V(MAXWEIGHT)));
    SUBVALUE("THRESHOLDWINC", V(THRESHOLDWINC) > V(LTHRESHOLDINC) ? V(THRESHOLDWINC) : 0.);
    SUBVALUE("LTHRESHOLDINC", V(THRESHOLDWINC) > V(LTHRESHOLDINC) ? 0. : V(LTHRESHOLDINC));

    SUBV(ACTIONSTEP);   // racket sfift causing ACTIONS firing
    SUBVINT(MAXIDLETIME);  // no-action time after which external stimalation is switched on

    // secondary evaluation neuron parameters

    SUBVINT(SECCHARTIME);
    SUBVINT(SECTHRESHOLDDECAYPERIOD);
    SUBV(SECTHRESHOLDINC);
    SUBV(MINPOTSEC);

    SUBVINT(LCHARTIME);
    SUBVINT(LTHRESHOLDDECAYPERIOD);
    SUBVINT(PLASTICITYTIME);
    SUBV(MAXWEIGHT);
    SUBV(MINWEIGHT);
    SUBV(STABILITYRESOURCECHANGERATIO);
    SUBVINT(NSILENTSYNAPSES);
    SUBVINT(GOODSTATEINERTIA);
    SUBVINT(STATECHANGEDELAY);
    SUBV(POSWSEC);
    SUBVINT(LATERALBLOCKINGTIME);
    SUBV(RELPUNEFFECT);
    SUBVINT(ACTIONEVALUATIONTIME);

    if (V(MIGRATION) < 0.5)
        EXCLUDESECTION("MIGRATION");
//------------------------------------------------------------------------------------------------

#ifndef ARNI_CD
    unsigned CurrentExperiment = 1000 + ind;
#else
    unsigned CurrentExperiment = InterlockedIncrement(&Experiment);
#endif

    char buf[3000];
    sprintf(buf, ".." FOLDER_SEPARATOR "Experiments" FOLDER_SEPARATOR "%d.nnc", CurrentExperiment);
	ofstream ofs(buf);
	SAVENNC(ofs);
	ofs.close();
	mutLog.lock();
	ofsTaskLog << "-------------------------------\n";
	SAVENNC(ofsTaskLog);
	ofsTaskLog << "-------------------------------\n";
	mutLog.unlock();
	stringstream ssNNC;
	SAVENNC(ssNNC);
	NNCDONE;

	// Run simulation

    sprintf(buf, "monitoring.%d.csv", CurrentExperiment);
	string strMonitoringFile(buf);
    sprintf(buf, "-e%d", CurrentExperiment);
    int ntrials = 0;
    boost::this_thread::sleep(boost::posix_time::milliseconds((CurrentExperiment % 101) * 300));
	int sum;
	vector<string> vstr_argsEmulator;
	do {
		do {
			mutsoc.lock();
			SelectedDevice = min_element(vpn_ComputingDeviceStatus.begin(), vpn_ComputingDeviceStatus.end()) - vpn_ComputingDeviceStatus.begin();
			++vpn_ComputingDeviceStatus[SelectedDevice].first;
			mutsoc.unlock();
			mutLog.lock();
            ofsTaskLog << "Experiment " << CurrentExperiment << " - " << EMULATOR << " started at " << boost::posix_time::second_clock::local_time() << std::endl;
            mutLog.unlock();
            vstr_argsEmulator.resize(3);
			vstr_argsEmulator[0] = ".." FOLDER_SEPARATOR "Experiments";
			vstr_argsEmulator[1] = buf;
            vstr_argsEmulator[2] = "-R";
//------------------------------------------------------------------------------------------------
// $$$ SPECIFY ADDITIONAL EMULATOR COMMAND LINE OPTIONS
            vstr_argsEmulator.back() = "-RS";
            vstr_argsEmulator.push_back("-T" + str(MAXEMULATIONTIME));
//------------------------------------------------------------------------------------------------
#ifndef USE_CPU
			sprintf(buf + 100, "-C%d", SelectedDevice);
			vstr_argsEmulator.push_back(buf + 100);
#endif

			sum = 0;
            vector<int> v_ReturnedValues;
			FORI(NTESTS) {
				a = ExecuteProcess(EMULATOR, vstr_argsEmulator, ofsTaskLog);
				if (a < 0) {
					printf("\a\a\a Abnormal " EMULATOR " termination -- trying to restart few minutes later!!!\n");
					mutLog.lock();
                    ofsTaskLog << "Experiment " << CurrentExperiment << " abnormal " EMULATOR " termination -- trying to restart (" << boost::posix_time::second_clock::local_time() << ")\n";
                    mutLog.unlock();
					break;
				}
				if (a >= SHOULD_BE_LAUNCHED_ON_ANOTHER_NODE)
					break;
				sum += a;
                v_ReturnedValues.push_back(a);
			}
			mutsoc.lock();
			--vpn_ComputingDeviceStatus[SelectedDevice].first;
			mutsoc.unlock();
			mutLog.lock();
            ofsTaskLog << "Experiment " << CurrentExperiment << " - " << EMULATOR << " finished at " << boost::posix_time::second_clock::local_time() <<  " returned";
            for (auto z: v_ReturnedValues)
                ofsTaskLog << ' ' << z;
            ofsTaskLog << std::endl;
            mutLog.unlock();
			if (a >= 0)
				break;
			boost::this_thread::sleep(boost::posix_time::seconds(200));
			ntrials++;
		} while (ntrials < 30);
		if (ntrials == 30) 
			do {
				printf("CANNOT SPAWN!\n");
				boost::this_thread::sleep(boost::posix_time::seconds(10));
			} while (true);
		if (a < SHOULD_BE_LAUNCHED_ON_ANOTHER_NODE) {
			remove(strMonitoringFile.c_str());
            sprintf(buf, "%d.finsave.nns", CurrentExperiment);
			remove(buf);
			return(sum / NTESTS);
		}
		mutsoc.lock();
		++vpn_ComputingDeviceStatus[SelectedDevice].first;
		mutsoc.unlock();
		string strLog;
		time_t tim = time(NULL);
		ctime_s(buf + 200, 100, &tim);
		mutLog.lock();
        ofsTaskLog << "Experiment " << CurrentExperiment << " - task sent to another cluster node at " << boost::posix_time::second_clock::local_time() << " (device " << SelectedDevice << ")" << endl;
        mutLog.unlock();
        int ret = CoordinateCluster(a - SHOULD_BE_LAUNCHED_ON_ANOTHER_NODE, SelectedDevice, REMOTE_EXPERIMENT_ID, vstr_argsEmulator, ssNNC.str(), strLog);
        if (strLog.length()) {
			mutsoc.lock();
			--vpn_ComputingDeviceStatus[SelectedDevice].first;
			mutsoc.unlock();
			mutLog.lock();
            ofsTaskLog << "Experiment " << CurrentExperiment << " result: " << ret << endl;
            ofsTaskLog << strLog << "\n-------------------\n";
			mutLog.unlock();
			return ret;
		}
		mutLog.lock();
        ofsTaskLog << "Experiment " << CurrentExperiment << " -- too long waiting for result -- restarting on another device\n";
        mutLog.unlock();
	} while (true);
}

#ifndef ARNI_CD
double Chromosome::operator()(int ind) const {return devaluate(ind, vprvpindb_);}

void Chromosome::ApplyAdditionalCriteria(VECTOR<double> &vd_AdditionalCriterionValues) const 
{
}
#endif
