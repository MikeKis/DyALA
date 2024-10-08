﻿//#define USE_CPU
#define EMULATOR "ArNIGPU"

#define MAXEMULATIONTIME 2000000
#define NTESTS 3

#define SECOND_TEST_THRESHOLD 1000000    // prohibitive value - to avoid running structural test
#define RECOGNIZER   // it should be retained.

CHROMOSOME_DEFINITION
{
    INITIALIZE_GENE(STABILITYRESOURCECHANGERATIO, rng.logUniform(0.03, 30)),
	INITIALIZE_GENE(NNEURONSPERACTION, rng.logUniform(10, 100)),
	INITIALIZE_GENE(LCHARTIME, rng.logUniform(10, 30)),
	INITIALIZE_GENE(LTHRESHOLDDECAYPERIOD, rng.logUniform(30, 100)),
	INITIALIZE_GENE(LTHRESHOLDINC, rng.logUniform(3, 30)),
	INITIALIZE_GENE(LPLASTICITYTIME, rng.logUniform(3, 30)),
	INITIALIZE_GENE(LDOPAMINEPLASTICITYTIME, rng.logUniform(10, 30)),
	INITIALIZE_GENE(MIGRATION, rng(2)),
	INITIALIZE_GENE(ACTIVATIONTIME, rng.logUniform(5, 30)),
	INITIALIZE_GENE(WEIGHTCHANGE, rng.logUniform(0.02, 0.5)),
	INITIALIZE_GENE(WEIGHTRELATIVEPUNISHMENT, rng.logUniform(5, 20)),
	INITIALIZE_GENE(MAXWEIGHT, rng.logUniform(2, 30)),
	INITIALIZE_GENE(DEFW, rng.logUniform(0.5, 3)),
	INITIALIZE_GENE(MAXDELAY, 3 + rng(7)),
	INITIALIZE_GENE(NSILENTINPUTS, rng.logUniform(100, 800)),
    INITIALIZE_GENE(LLATERALBLOCKINGTIME, 4 + rng(4))
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

    SETNNP(".." FOLDER_SEPARATOR "Experiments" FOLDER_SEPARATOR "ExtSNNREW.nnp");
    SUBGENE(STABILITYRESOURCECHANGERATIO);
	SUBGENEINT(NNEURONSPERACTION);
	SUBGENEINT(LCHARTIME);
	SUBGENEINT(LTHRESHOLDDECAYPERIOD);
	SUBGENEINT(LTHRESHOLDINC);
	SUBGENEINT(LPLASTICITYTIME);
	if (!GV(MIGRATION))
		EXCLUDESECTION("MIGRATION");
	SUBGENEINT(ACTIVATIONTIME);
	SUBGENEINT(LDOPAMINEPLASTICITYTIME);
	SUBGENE(WEIGHTCHANGE);
	SUBGENE(WEIGHTRELATIVEPUNISHMENT);
	SUBGENE(MAXWEIGHT);
	SUBGENE(DEFW);
	SUBGENEINT(MAXDELAY);
	SUBGENEINT(NSILENTINPUTS);
	SUBGENEINT(LLATERALBLOCKINGTIME);
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
			vstr_argsEmulator.resize(7);
			vstr_argsEmulator[0] = ".." FOLDER_SEPARATOR "Experiments";
			vstr_argsEmulator[1] = buf;
			vstr_argsEmulator[2] = buf + 300;
			vstr_argsEmulator[3] = "-S100000";
			vstr_argsEmulator[4] = "-R";
			vstr_argsEmulator[5] = "-dS30000";
			vstr_argsEmulator[6] = "-D100000";
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
			pnTasksDoneonCard[SelectedDevice]++;
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
		int ret = CoordinateCluster(a - SHOULD_BE_LAUNCHED_ON_ANOTHER_NODE, SelectedDevice, ind, vstr_argsEmulator, ssNNC.str(), strLog, pnTasksDoneonCard);
		if (strLog.length()) {
			mutsoc.lock();
			--vpn_ComputingDeviceStatus[SelectedDevice].first;
			mutsoc.unlock();
			mutLog.lock();
			ofsTaskLog << "ind " << ind << " result: " << ret << endl;
			ofsTaskLog << strLog << "\n-------------------\n";
			pnTasksDoneonCard[SelectedDevice]++;
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
