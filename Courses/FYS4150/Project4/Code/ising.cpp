#include "ising.h"

Ising::Ising(int L_, double T_, bool random) {
    L = L_;
    T = T_;
    A = 0;
    // Setting up random number generator
    random_device rd;
    mt19937_64 gen(rd());
    discrete_distribution<> rspin({1, 0, 1}); // equally likely to return 0 or 2
    // Initial lattice and magnetic moment B value
    lattice = zeros<mat>(L, L);
    L2 = L * L;
    M = 0;
    for (int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++) {
            if (random) {
                lattice(x, y) = (1 - rspin(gen)); // random spin orientation, as the behaviour around chaotic spins is the most interesting
            } else {
                lattice(x, y) = 1;
            }
            M += (double) lattice(x, y);
        }
    }
    // Initial energy E. Each spin counts the energy from the interaction with the spin under and to the left
    E = 0;
    for(int x = 0; x < L; x++) {
        for (int y = 0; y < L; y++){
            E -= (double) lattice(x, y) * (lattice(periodic(x, -1), y) + lattice(x, periodic(y, -1)));
        }
    }
    // deltaE
    deltaEProb = zeros<vec>(17);
    for (int deltaE = -8; deltaE <= 8; deltaE += 4) {
        deltaEProb(deltaE + 8) = exp(-deltaE / T); // The probabilities will be at index 0, 4, 8, 12, 16, for the E arguments -8, -4, 0, 4, 8
    }
}

void Ising::calcEVs(int cycles) {
    // Setting up random number generator
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> rdist(0.0, 1.0);
    // Initializing current cycle
    int cycle = 0;
    int badCycles = min(10000, cycles / 10);
    double goodCycles = cycles - badCycles;
    // We discard the first 10% of the cycles
    for (; cycle < badCycles; cycle++) {
        flipCycle(gen, rdist);
    }
    for (; cycle < cycles; cycle++) {
        flipCycle(gen, rdist);
        // Updating expectation values each loop
        EVs(0) += E;
        EVs(1) += E * E;
        EVs(2) += M;
        EVs(3) += M * M;
        EVs(4) += fabs(M);
    }
    // Scaling down all expectation values
    EVs(0) /= (double) goodCycles * L2;
    EVs(1) /= (double) goodCycles * L2;
    EVs(2) /= (double) goodCycles * L2;
    EVs(3) /= (double) goodCycles * L2;
    EVs(4) /= (double) goodCycles * L2;
    Arate = (double) A / cycles / L2;
}

void Ising::calcWriteChange(int cycles, string filename) {
    //
    double E_EV = E;
    double M_EV = M;
    double E2_EV = E * E;
    // Setting up ofile
    ofstream ofile(filename);
    ofile << setw(15) << "E";
    ofile << setw(15) << "M";
    ofile << setw(15) << "|M|";
    ofile << setw(15) << "A";
    ofile << setw(15) << "<E>";
    ofile << setw(15) << "<|M|>";
    ofile << setw(15) << "<E2>";
    ofile << setw(15) << "T" << endl;
    ofile << setw(15) << setprecision(8) << E / L2;
    ofile << setw(15) << setprecision(8) << M / L2;
    ofile << setw(15) << setprecision(8) << fabs(M) / L2;
    ofile << setw(15) << setprecision(8) << A;
    ofile << setw(15) << setprecision(8) << E_EV / L2;
    ofile << setw(15) << setprecision(8) << M_EV / L2;
    ofile << setw(15) << setprecision(8) << E2_EV / L2;
    ofile << setw(15) << setprecision(8) << T << endl;
    // Setting up random number generator
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> rdist(0.0, 1.0);
    // Calculating the MC metropolis flips
    for (int cycle = 0; cycle < cycles; cycle++) {
        flipCycle(gen, rdist);
        E_EV += E;
        M_EV += fabs(M);
        E2_EV += E * E;
        // Writing values to file each loop
        ofile << setw(15) << setprecision(8) << (double) E / L2;
        ofile << setw(15) << setprecision(8) << (double) M / L2;
        ofile << setw(15) << setprecision(8) << (double) fabs(M) / L2;
        ofile << setw(15) << setprecision(8) << A;
        ofile << setw(15) << setprecision(8) << (double) E_EV / (cycle + 2.0) / L2;
        ofile << setw(15) << setprecision(8) << (double) M_EV / (cycle + 2.0) / L2;
        ofile << setw(15) << setprecision(8) << (double) E2_EV / (cycle + 2.0) / L2 << endl;
    }
}

void Ising::flipCycle(mt19937_64 &gen, uniform_real_distribution<double> &rdist) {
    // We suggest as many flips as there are spins in the lattice, that is L * L = L2 flips
    int x, y, deltaE;
    for (int i = 0; i < L2; i++) {
        // Suggested flip, and associated energy change
        x = (int) (rdist(gen) * (double) L);
        y = (int) (rdist(gen) * (double) L);
        deltaE = 2 * lattice(x, y) * ( lattice(x, periodic(y, -1))
                                         + lattice(periodic(x, -1), y)
                                         + lattice(x, periodic(y, 1))
                                         + lattice(periodic(x, 1), y) );
        // Test to see if we accept flip. If deltaE is negative or 0, deltaEProb is larger than 1, and the move is always accepted
        if (rdist(gen) <= deltaEProb(deltaE + 8)) {
            lattice(x, y) *= -1.0; // flip spin
            M += (double) 2 * lattice(x, y);
            E += (double) deltaE;
            A += 1;
        }
    }
}

inline int Ising::periodic(int index, int direction) {
    return (index + L + direction) % (L); //if the index goes over or under the lattice size, the index loops around
}

void Ising::writeEVs(ofstream &ofile) {
    // Writes temperature for this lattice to file
    ofile << setw(15) << setprecision(8) << T;
    //Writes all EVs to file
    ofile << setw(15) << setprecision(8) << EVs(0); // E
    ofile << setw(15) << setprecision(8) << EVs(2); // M
    ofile << setw(15) << setprecision(8) << EVs(4); // |M|
    ofile << setw(15) << setprecision(8) << (EVs(1) - EVs(0) * EVs(0) * L2) / T / T; // Cv = (<E^2> - <E>^2) T^2
    ofile << setw(15) << setprecision(8) << (EVs(3) - EVs(4) * EVs(4) * L2) / T; // X = <M^2> - <|M|>^2
    ofile << setw(15) << setprecision(8) << Arate << endl;
}
