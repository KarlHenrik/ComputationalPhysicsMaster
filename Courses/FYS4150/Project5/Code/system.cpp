#include "system.h"

// 1D system
System::System(vec &initial, double dx_, string ofilename, double t0, double D_, string sourceName_) {
    u = initial;
    dx = dx_;
    t = t0;
    t_i = 0;
    D = D_;
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    sourceName = sourceName_;
    dims = 1;
}

// 2D system
System::System(mat initial, double dx_, string ofilename, double t0, double D_) {
    U = initial;
    dx = dx_;
    t = t0;
    t_i = 0;
    D = D_;
    ofile.open(ofilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    dims = 2;
}

void System::writeState() {
    ofile << "t" << t_i << " "; // How many times the state has been written to file
    t_i++;
    if (dims == 1) {
        for (double T_i : u) {
            ofile << T_i << " ";
        }
    } else if (dims == 2) {
        int n = U.n_rows;
        for (int x = 0; x < n; x++) {
            for (int y = 0; y < n; y++) {
                ofile << U(x, y) << " ";
            }
        }
    }
    ofile << endl;
}

void System::close() {
    ofile.close();
}

double System::source(double x, double t) {
    if (sourceName == "base") {
        return heatProd(x);
    } else if (sourceName == "enrichedConst") {
        return heatProdEnrichedConst(x);
    } else if (sourceName == "enriched") {
        return heatProdEnriched(x, t);
    }
    return 0;
}

double System::heatProd(double x) {
    if (x < 0.1666666) {
        return 1.4 * Qfac;
    } else if (x < 0.3333333) {
        return 0.35 * Qfac;
    } else {
        return 0.05 * Qfac;
    }
}

double System::heatProdEnrichedConst(double x) {
    if (x < 0.1666666) {
        return 1.4 * Qfac;
    } else if (x < 0.3333333) {
        return 0.35 * Qfac;
    } else {
        return 0.55 * Qfac; // Added 0.5 extra production in the mantle
    }
}

double System::heatProdEnriched(double x, double t) {
    if (x < 0.1666666) {
        return 1.4 * Qfac;
    } else if (x < 0.333333) {
        return 0.35 * Qfac;
    } else {
        return (0.05 + 0.5 *
                            (0.4 * exp(-0.1550664833467439 * t) + //U
                             0.4 * exp(-0.0495105128971389 * t) + //Th
                             0.2 * exp(-0.5545177444479562 * t))  //K
               ) * Qfac;
    }
}
