//#define USE_CPU
#define EMULATOR "ArNIGPU"
#define RECOGNIZER ""

#define LEARNING_TIME 1468200   // Because half of all time is silence.
#define STATE_PRESENTATION_TACTS 10   // Do not change: LCharTime = STATE_PRESENTATION_TACTS / 3
#define INTERSTATE_SILENCE_TACTS 10
#define NTACTPERCASE (STATE_PRESENTATION_TACTS + INTERSTATE_SILENCE_TACTS)
#define SPIKING_INPUT "inpstaticperm.txt"
#define TARGET_FILE "inpstatictargetperm.txt"
#define NINPUTS 133
#define NTESTS 3
//#define MAXNNEU 1000
#define SECOND_TEST_THRESHOLD 0x7fffffff   // DO NOT USE PROTOCOL POSTPROCESSOR !!
#define NTARGETCLASSES 1
#define NNETWORKCOPIES 1                   // It will be optimized by cooredinate descent
#define LEARNINGGROUPSIZE 1
#define TARGETBIAS 3                       // ONLY HERE!! In fact, it is ceil(THRESHOLD_BASE * 0.001) / (LCharTime * (1 - pow(1 - 1. / LCharTime, STATE_PRESENTATION_TACTS)))

NetworkTestingParams ntpG{EMULATOR, RECOGNIZER, NTESTS, SECOND_TEST_THRESHOLD, NoStructuralCriterion, 0}; 


unsigned Experiment = 3000;

int evaluate(const vector<int> &v_point, const vector<dimension> &vdim_, bool &bTerminate)
{

	if (DV("NLEARNINGGROUPS") < 1)
    		return 0;
	int a, SelectedDevice;

// Create network configuration file and storing its contents in log file

    SETNNP(".." FOLDER_SEPARATOR "Experiments" FOLDER_SEPARATOR "Classifier.nnp");

    SUBMACRO(NNETWORKCOPIES);
    SUBMACRO(LEARNING_TIME);
    SUBMACRO(SPIKING_INPUT);
    SUBMACRO(TARGET_FILE);
    SUBVALUE("NSTATES", NTARGETCLASSES);
    SUBMACRO(NINPUTS);
    SUBVALUE("NO_CLASS_LABEL", "0");
    SUBVALUE("CRITERION", "averaged_F");

    SUBVALUE("OUTPUT", "OUT");
    SUBVALUE("INPUT", "R");
    SUBVALUE("REWARD", "Target");

    SUBMACRO(NTACTPERCASE);

    SUBVALUE("NLEARNERS", NTARGETCLASSES * nint(DV("NLEARNINGGROUPS")) * LEARNINGGROUPSIZE);
    SUBVALUE("MAXOUTISI", STATE_PRESENTATION_TACTS);
    int LCharTime = STATE_PRESENTATION_TACTS / 3;
    SUBVALUE("LCHARTIME", LCharTime);
    SUBVALUE("DWNEG", DV("DWPOS") * DV("RELDWNEG"));
    
    SUBVALUE("INIRESOURCE0", rWeight0Resource(-DV("MINWEIGHT"), DV("MAXWEIGHT")));

    SUBVALUE("SPARCITY", 1.);
    SUBVALUE("DOPAMINEPLASTICITYTIME", STATE_PRESENTATION_TACTS);
    SUBVALUE("TARGETSTIMULATIONPERIOD", STATE_PRESENTATION_TACTS + INTERSTATE_SILENCE_TACTS);
    SUBMACRO(TARGETBIAS);

    SUBDIM("THRESHOLDWINC");
    SUBDIM("DWPOS");
    SUBDIM("MAXWEIGHT");
    SUBDIM("MINWEIGHT");
    SUBDIM("STABILITYRESOURCECHANGERATIO");
    SUBDIM("HEBBIANPLASTICITYCHARTIMERATIO");
    SUBDIMINT("NSILENTSYNAPSES");
    SUBDIMINT("TARGETSTIMULATIONPERIOD");

    EXCLUDESECTION("SEQUENTIAL_TRAIN");

    unsigned CurrentExperiment = InterlockedIncrement(&Experiment);

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
            vstr_argsEmulator.resize(4);
			vstr_argsEmulator[0] = ".." FOLDER_SEPARATOR "Experiments";
			vstr_argsEmulator[1] = buf;
			vstr_argsEmulator[3] = "-R";
			vstr_argsEmulator[2] = "-f" + str(LEARNING_TIME);
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
		int ret = CoordinateCluster(a - SHOULD_BE_LAUNCHED_ON_ANOTHER_NODE, SelectedDevice, CurrentExperiment, vstr_argsEmulator, ssNNC.str(), strLog);
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

