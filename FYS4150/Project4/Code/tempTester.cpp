#include "tempTester.h"

TempTester::TempTester(int L, double T0, double TN, double dT, bool random) {
    for (double T = T0; T <= TN; T += dT) {
        models.push_back( Ising(L, T, random) );
    }
}

void TempTester::calc(int cycles) {
    for (Ising &mdl : models) {
        mdl.calcEVs(cycles);
    }
}

void TempTester::write(string ofilename) {
    ofstream ofile(ofilename);
    ofile << setw(15) << "T";
    ofile << setw(15) << "E";
    ofile << setw(15) << "M";
    ofile << setw(15) << "|M|";
    ofile << setw(15) << "Cv";
    ofile << setw(15) << "X";
    ofile << setw(15) << "A" << endl;
    for (Ising &mdl : models) {
        mdl.writeEVs(ofile);
    }
    ofile.close();
}

void TempTester::calcParallell(int cycles, int processes) {
    int jobs = models.size();
    int job_size = jobs / processes;
    int extra_jobs = jobs % processes;
    int extra_jobs_start = job_size * processes;

    omp_set_num_threads(processes);
    double wtime = omp_get_wtime();
    # pragma omp parallel shared(models)
    {
        // Calculating EVs for the evenly distributed models
        for (int job = 0; job < job_size; job++) {
            models[omp_get_thread_num() * job_size + job].calcEVs(cycles);
        }
        // Calculating EVs for leftover models
        if (omp_get_thread_num() < extra_jobs) {
            models[extra_jobs_start + omp_get_thread_num()].calcEVs(cycles);
        }
    }
    wtime = omp_get_wtime() - wtime;
    cout << "Elapsed time in seconds = " << wtime << endl;
}
