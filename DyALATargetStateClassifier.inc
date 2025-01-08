//#define USE_CPU
#define EMULATOR "ArNIGPU"
#define RECOGNIZER ""

//------------------------------------------------------------------------------------------------
// $$$ COMMON DEFINES
#define MAXEMULATIONTIME 2000000
#define LEARNINGTIME 1000000
#define NINPUTS 150
#define CAUSAL_LINK_DURATION 10
#define NCONSEQUENCES 2
#define NTESTS 3
#define SECOND_TEST_THRESHOLD 0x7fffffff   // DO NOT USE PROTOCOL POSTPROCESSOR !!
//------------------------------------------------------------------------------------------------

NetworkTestingParams ntpG{EMULATOR, RECOGNIZER, NTESTS, SECOND_TEST_THRESHOLD, NoStructuralCriterion, STARTING_EXPERIMENT_ID};   // It is because experiment ID space for GA is 1000-2000. The original experminet ID begins from 0. For local execution it is incremented
// by 1000. This incremented value is passed through -e option. But in cluster management code this task is registered by the original ID.
// Logics for CD are different!!

#ifndef ARNI_CD
//------------------------------------------------------------------------------------------------
// $$$ DEFINE CHROMOSOME

CHROMOSOME_DEFINITION
{
    INITIALIZE_GENE(NNEURONSPERCONSEQUENCE, rng.logUniform(1, 30)),
    INITIALIZE_GENE(LCHARTIME, rng.logUniform(1, 30)),
    INITIALIZE_GENE(LEARNINGRATE, rng.logUniform(0.1, 3.)),
    INITIALIZE_GENE(MAXWEIGHT, rng.logUniform(0.3, 3.)),
    INITIALIZE_GENE(MINWEIGHT, rng.logUniform(0.003, 3.)),
    INITIALIZE_GENE(STABILITYRESOURCECHANGERATIO, rng.logUniform(0.03, 30)),
    INITIALIZE_GENE(NSILENTSYNAPSES, rng.logUniform(1, 1000))
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
bool bBadHyperparameters(const vector<int> &v_point, const vector<dimension> &vdim_) {return V(NNEURONSPERACTION) < 1 || V(PUNISHMENT_SWITCH_THTRESHOLD) > 0.5 || V(REWARD_SWITCH_THTRESHOLD) > 0.5;}

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
    SETNNP(".." FOLDER_SEPARATOR "Experiments" FOLDER_SEPARATOR "TargetStateClassification.nnp");
//------------------------------------------------------------------------------------------------
// $$$ SUBSTITUTE ALL VALUES

    SUBMACRO(NINPUTS);
    SUBMACRO(LEARNINGTIME);
    SUBMACRO(CAUSAL_LINK_DURATION);
    SUBMACRO(NCONSEQUENCES);

    SUBVALUE("CAUSES", "R");
    SUBVALUE("CONSEQUENCES", "Evaluation");

    SUBVINT(NNEURONSPERCONSEQUENCE);
    SUBVALUE("NLEARNERS", nint(V(NNEURONSPERCONSEQUENCE)) * NCONSEQUENCES);
    SUBVINT(LCHARTIME);
    SUBV(LEARNINGRATE);
    SUBV(MAXWEIGHT);
    SUBV(MINWEIGHT);
    SUBV(STABILITYRESOURCECHANGERATIO);
    SUBVINT(NSILENTSYNAPSES);
    SUBVALUE("INIRESOURCE0", rWeight0Resource(-V(MINWEIGHT), V(MAXWEIGHT)));

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
            vstr_argsEmulator.back() = "-R";
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
