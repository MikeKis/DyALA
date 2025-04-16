//#define USE_CPU
#define EMULATOR "ArNIGPU"
#define RECOGNIZER ""

//------------------------------------------------------------------------------------------------
// $$$ COMMON DEFINES
#define MAXEMULATIONTIME 2000000
#define NTESTS 3
#define SECOND_TEST_THRESHOLD 0x7fffffff   // DO NOT USE PROTOCOL POSTPROCESSOR !!
#define LEARNINGGROUPSIZE 1                // because afferent connections are all-to-all
//------------------------------------------------------------------------------------------------

NetworkTestingParams ntpG{EMULATOR, RECOGNIZER, NTESTS, SECOND_TEST_THRESHOLD, NoStructuralCriterion};

//------------------------------------------------------------------------------------------------
// $$$ DEFINE OPTIMIZED VARIABLES - only variable names are used in code - to determine variable position. The array is read by an optimization control program.
OPTIMIZATION_SPACE_DEFINITION
{
    OPTIMIZED_VARIABLE_INT(MAXSTIMULATIONTIMEENTER, 3, 10),
    OPTIMIZED_VARIABLE_INT(MAXSTIMULATIONTIMECONT, 3, 10),
    OPTIMIZED_VARIABLE_INT(NLEARNINGGROUPSENTER, 1, 10),
    OPTIMIZED_VARIABLE_INT(NLEARNINGGROUPSCONT, 1, 10),
    OPTIMIZED_VARIABLE(WTHREXCESSENTER, 0.001, 0.3),
    OPTIMIZED_VARIABLE(WTHREXCESSCONT, 0.001, 0.3),
    OPTIMIZED_VARIABLE(DWPOSENTER, 0.03, 1.),
    OPTIMIZED_VARIABLE(DWPOSCONT, 0.03, 1.),
    OPTIMIZED_VARIABLE(RELDWNEGENTER, 0., 2.),
    OPTIMIZED_VARIABLE(RELDWNEGCONT, 0., 2.),
    OPTIMIZED_VARIABLE(MAXWEIGHTENTER, 0.3, 6.),
    OPTIMIZED_VARIABLE(MINWEIGHTENTER, 0.03, 6.),
    OPTIMIZED_VARIABLE(MAXWEIGHTCONT, 0.3, 6.),
    OPTIMIZED_VARIABLE(MINWEIGHTCONT, 0.03, 6.),
    OPTIMIZED_VARIABLE(STABILITYRESOURCECHANGERATIOENTER, 0.03, 30),
    OPTIMIZED_VARIABLE(STABILITYRESOURCECHANGERATIOCONT, 0.03, 30),
    OPTIMIZED_VARIABLE_INT(NSILENTSYNAPSES, 1, 3000),
    OPTIMIZED_VARIABLE_INT(CHARTIMEENTER, 2, 30),
    OPTIMIZED_VARIABLE_INT(CHARTIMECONT, 2, 30),
    OPTIMIZED_VARIABLE(HEBBIANPLASTICITYCHARTIMERATIO, 2, 10)
};

//#define NO_GENE_DEPENDENCIES - not implemented now
//------------------------------------------------------------------------------------------------

// Check that the given hyperparameter combination makes sense

bool bNoFatalDefect(const vector<float> &vr_) {return V(NLEARNINGGROUPSENTER) >= 1 && V(NLEARNINGGROUPSCONT) >= 1;}
//------------------------------------------------------------------------------------------------

int StructuralCriterion(int)
{
    throw std::logic_error("StructuralCriterion used!!!");
    return 0;
}   // not used

EVALUATOR
{
	int a, SelectedDevice;
    if (!bNoFatalDefect(vr_))
            return 0;

// Create network configuration file and storing its contents in log file

//------------------------------------------------------------------------------------------------
// $$$ SPECIFY NNP FILE
    SETNNP(".." FOLDER_SEPARATOR "Experiments" FOLDER_SEPARATOR "WorldDynamics.nnp");

//------------------------------------------------------------------------------------------------
// $$$ SUBSTITUTE ALL VALUES

    SUBVINT(MAXSTIMULATIONTIMEENTER);
    SUBVINT(MAXSTIMULATIONTIMECONT);
    SUBVINT(NLEARNINGGROUPSENTER);
    SUBVINT(NLEARNINGGROUPSCONT);
    SUBV(WTHREXCESSENTER);
    SUBV(WTHREXCESSCONT);
    SUBV(DWPOSENTER);
    SUBV(DWPOSCONT);
    SUBVALUE("RELDWNEGENTER", min(V(RELDWNEGENTER), 1.F));
    SUBVALUE("RELDWNEGCONT", min(V(RELDWNEGCONT), 1.F));
    SUBV(MAXWEIGHTENTER);
    SUBV(MINWEIGHTENTER);
    SUBV(MAXWEIGHTCONT);
    SUBV(MINWEIGHTCONT);
    SUBV(STABILITYRESOURCECHANGERATIOENTER);
    SUBV(STABILITYRESOURCECHANGERATIOCONT);
    SUBVINT(NSILENTSYNAPSES);
    SUBVINT(CHARTIMEENTER);
    SUBVINT(CHARTIMECONT);
    SUBV(HEBBIANPLASTICITYCHARTIMERATIO);

//    EXCLUDESECTION("CONVOLUTION");
//------------------------------------------------------------------------------------------------

    unsigned CurrentExperiment = STARTING_EXPERIMENT_ID + indinPack;

    char buf[3000];
    sprintf(buf, ".." FOLDER_SEPARATOR "Experiments" FOLDER_SEPARATOR "%d.nnc", CurrentExperiment);
	ofstream ofs(buf);
	SAVENNC(ofs);
	ofs.close();
	mutLog.lock();
    ofsTaskLog << "Testing network " << hex << idUnique << dec << " (local experimrnt id " << CurrentExperiment << ")" << endl;
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
			sprintf(buf + 300, "-T%d", MAXEMULATIONTIME);
            		vstr_argsEmulator.resize(4);
			vstr_argsEmulator[0] = ".." FOLDER_SEPARATOR "Experiments";
			vstr_argsEmulator[1] = buf;
			vstr_argsEmulator[2] = "-R";
			vstr_argsEmulator[3] = buf + 300;
//------------------------------------------------------------------------------------------------
// Set additional simulator command line options. 
//			vstr_argsEmulator.push_back("-f" + str(LEARNING_TIME));
//------------------------------------------------------------------------------------------------
#ifndef USE_CPU
			sprintf(buf + 100, "-C%d", SelectedDevice);
			vstr_argsEmulator.push_back(buf + 100);
#endif

			sum = 0;
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
			}
			mutsoc.lock();
			--vpn_ComputingDeviceStatus[SelectedDevice].first;
			mutsoc.unlock();
			mutLog.lock();
			ofsTaskLog << "Experiment " << CurrentExperiment << " - " << EMULATOR << " finished at " << boost::posix_time::second_clock::local_time() <<  " returned  " << a << std::endl;
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
        int ret = CoordinateCluster(a - SHOULD_BE_LAUNCHED_ON_ANOTHER_NODE, SelectedDevice, vstr_argsEmulator, ssNNC.str(), strLog);
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
