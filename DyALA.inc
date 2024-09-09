#define USE_CPU

#ifndef USE_CPU
#define EMULATOR "ArNIGPU"
#define NTESTS 3
#else
#define EMULATOR "ArNICPU"
#define NTESTS 1
#endif

#define MAXEMULATIONTIME 2000000

#define SECOND_TEST_THRESHOLD 1000000    // prohibitive value - to avoid running structural test
#define RECOGNIZER ""  // it should be retained.

NetworkTestingParams ntpG{EMULATOR, RECOGNIZER, NTESTS, SECOND_TEST_THRESHOLD, NoStructuralCriterion, 1000};   // It is because experiment ID space for GA is 1000-2000. The original experminet ID begins from 0. For local execution it is incremented
// by 1000. This incremented value is passed through -e option. But in cluster management code this task is registered by the original ID.
// Logics for CD are different!!

CHROMOSOME_DEFINITION
{
    INITIALIZE_GENE(NLEVELS, 3 + rng(5)),
    INITIALIZE_GENE(SPONTANEOUSACTIVITY, rng.logUniform(0.0001, 0.1)),
    INITIALIZE_GENE(LEVELDURATION, rng.logUniform(30, 300)),
    INITIALIZE_GENE(NGROUPSONLEVEL0, rng.logUniform(1, 10)),
    INITIALIZE_GENE(LEVELNGROUPMULTIPLIER, rng.logUniform(0.5, 2.)),
	INITIALIZE_GENE(LCHARTIME, rng.logUniform(1, 30)),
	INITIALIZE_GENE(NSILENTSYNAPSES, rng.logUniform(1, 1000)),
        INITIALIZE_GENE(DW, rng.logUniform(0.1, 10.)),
	INITIALIZE_GENE(MAXWEIGHT, rng.logUniform(1, 30)),
	INITIALIZE_GENE(MINWEIGHT, rng.logUniform(0.0001, 10)),
    INITIALIZE_GENE(STABILITYRESOURCECHANGERATIO, rng.logUniform(0.01, 10)),
	INITIALIZE_GENE(MAXISILEVELINDICATOR, rng.logUniform(3, 30)),
	INITIALIZE_GENE(NNEURONSPERACTION, rng.logUniform(1, 30)),
	INITIALIZE_GENE(LACTCHARTIME, rng.logUniform(1, 30)),
	INITIALIZE_GENE(LDOPAMINEPLASTICITYTIME, rng.logUniform(10, 30)),
    INITIALIZE_GENE(STABILITYRESOURCECHANGERATIOACT, rng.logUniform(0.03, 30)),
	INITIALIZE_GENE(NACTSILENTSYNAPSES, rng.logUniform(1, 1000)),
	INITIALIZE_GENE(MAXWEIGHTACT, rng.logUniform(1, 30)),
	INITIALIZE_GENE(MINWEIGHTACT, rng.logUniform(0.0001, 10)),
	INITIALIZE_GENE(MIGRATION, rng(2)),
	INITIALIZE_GENE(MAXACTISI, rng.logUniform(3, 30)),
        INITIALIZE_GENE(DWACT, rng.logUniform(0.003, 0.3)),
	INITIALIZE_GENE(MAXAFFERENTDELAYACT, 2 + rng(8)),
	INITIALIZE_GENE(ACTIONPOWERLEVELFACTOR, rng.logUniform(0.1, 3.)),
        INITIALIZE_GENE(BASERPLASTICITY, rng.logUniform(0.3, 10.)),
	INITIALIZE_GENE(RPLASTICITYLEVELFACTOR, rng())
};

#define NO_GENE_DEPENDENCIES

#ifndef FOR_LINUX
#pragma warning(disable:4244)
#pragma warning(disable:4305)
#endif

// Check that the given gene combination makes sense

bool Chromosome::bNoFatalDefect(void) const 
{
//	return GV(LDOPAMINEPLASTICITYTIME) > GV(NEURONTIMEDEPTH) + 2 && GV(ACTIVATIONTIME) > GV(NEURONTIMEDEPTH) + 2;
	return true;
}

int StructuralCriterion(int){return 0;}   // not used

double devaluate(int ind, const VECTOR<PAIR<float,VECTOR<PAIR<int,bool> > > > &vprvpindb_, double dConnectivityFactor = 1.)
{
	int a, SelectedDevice;

// Create network configuration file and storing its contents in log file

    SETNNP(".." FOLDER_SEPARATOR "Experiments" FOLDER_SEPARATOR "ChaoticOriginalLearning.nnp");
	SUBGENEINT(NLEVELS);
	SUBGENE(SPONTANEOUSACTIVITY);
	SUBGENEINT(LEVELDURATION);
	SUBGENEINT(NGROUPSONLEVEL0);
	SUBVALUE("LEVELNGROUPMULTIPLIER", max(GV(LEVELNGROUPMULTIPLIER), 1.F));
	SUBGENEINT(LCHARTIME);
	SUBGENEINT(NSILENTSYNAPSES);
	SUBGENE(DW);
	SUBGENE(MAXWEIGHT);
	SUBVALUE("MINWEIGHT", GV(MINWEIGHT) > 0.01 ? GV(MINWEIGHT) : 0.);
    SUBGENE(STABILITYRESOURCECHANGERATIO);
        SUBVALUE("MAXAFFERENTDELAY", 3);
	SUBGENEINT(MAXISILEVELINDICATOR);
	SUBGENEINT(NNEURONSPERACTION);
	SUBGENEINT(LACTCHARTIME);
	SUBGENEINT(LDOPAMINEPLASTICITYTIME);
    SUBGENE(STABILITYRESOURCECHANGERATIOACT);
	SUBGENEINT(NACTSILENTSYNAPSES);
	SUBGENE(MAXWEIGHTACT);
	SUBVALUE("MINWEIGHTACT", GV(MINWEIGHTACT) > 0.01 ? GV(MINWEIGHTACT) : 0.);
//	if (!GV(MIGRATION))
		EXCLUDESECTION("MIGRATION");
	SUBGENEINT(MAXACTISI);
	SUBGENE(DWACT);
	SUBGENEINT(MAXAFFERENTDELAYACT);
	SUBGENE(ACTIONPOWERLEVELFACTOR);
	SUBGENE(BASERPLASTICITY);
	SUBGENE(RPLASTICITYLEVELFACTOR);
    char buf[3000];
    sprintf(buf, ".." FOLDER_SEPARATOR "Experiments" FOLDER_SEPARATOR "%d.nnc", 1000 + ind);
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

	sprintf(buf, "monitoring.%d.csv", 1000 + ind);
	string strMonitoringFile(buf);
    sprintf(buf, "-e%d", 1000 + ind);
    int ntrials = 0;
	boost::this_thread::sleep(boost::posix_time::milliseconds((ind % 101) * 300));
	int sum;
	vector<string> vstr_argsEmulator;
	do {
		do {
			mutsoc.lock();
			SelectedDevice = min_element(vpn_ComputingDeviceStatus.begin(), vpn_ComputingDeviceStatus.end()) - vpn_ComputingDeviceStatus.begin();
			++vpn_ComputingDeviceStatus[SelectedDevice].first;
			mutsoc.unlock();
			mutLog.lock();
			ofsTaskLog << "ind" << ind << " - " << EMULATOR << " started at " << boost::posix_time::second_clock::local_time() << std::endl;
			mutLog.unlock();
			sprintf(buf + 300, "-T%d", MAXEMULATIONTIME);
			vstr_argsEmulator.resize(6);
			vstr_argsEmulator[0] = ".." FOLDER_SEPARATOR "Experiments";
			vstr_argsEmulator[1] = buf;
			vstr_argsEmulator[2] = buf + 300;
			vstr_argsEmulator[3] = "-D100000";
#ifndef USE_CPU
                        vstr_argsEmulator[4] = "-R";
#else
                        vstr_argsEmulator[4] = "-v2";
#endif
			vstr_argsEmulator[5] = "-dS30000";
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
					ofsTaskLog << "ind" << ind << " abnormal " EMULATOR " termination -- trying to restart (" << boost::posix_time::second_clock::local_time() << ")\n";
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
			ofsTaskLog << "ind" << ind << " - " << EMULATOR << " finished at " << boost::posix_time::second_clock::local_time() <<  " returned  " << a << std::endl;
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
			sprintf(buf, "spikes.%d.bin", 1000 + ind);
			string strProtocolFile(buf);
			remove(strProtocolFile.c_str());
			sprintf(buf, "%d.finsave.nns", 1000 + ind);
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
		ofsTaskLog << "ind" << ind << " - task sent to another cluster node at " << boost::posix_time::second_clock::local_time() << " (device " << SelectedDevice << ")" << endl;
		mutLog.unlock();
		int ret = CoordinateCluster(a - SHOULD_BE_LAUNCHED_ON_ANOTHER_NODE, SelectedDevice, ind, vstr_argsEmulator, ssNNC.str(), strLog);
		if (strLog.length()) {
			mutsoc.lock();
			--vpn_ComputingDeviceStatus[SelectedDevice].first;
			mutsoc.unlock();
			mutLog.lock();
			ofsTaskLog << "ind " << ind << " result: " << ret << endl;
			ofsTaskLog << strLog << "\n-------------------\n";
			mutLog.unlock();
			return ret;
		}
		mutLog.lock();
		ofsTaskLog << "ind " << ind << " -- too long waiting for result -- restarting on another device\n";
		mutLog.unlock();
	} while (true);
}

double Chromosome::operator()(int ind) const {return devaluate(ind, vprvpindb_);}

void Chromosome::ApplyAdditionalCriteria(VECTOR<double> &vd_AdditionalCriterionValues) const 
{
}
