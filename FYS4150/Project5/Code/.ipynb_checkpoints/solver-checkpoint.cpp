#include "solver.h"

void solve(System &system, string scheme, vector<int> intervals, double dt) {
    system.writeState();
    for (int t_steps : intervals) {
        schemeChooser(system, scheme, dt, t_steps);
        system.writeState();
    }
    system.close();
}

void schemeChooser(System &system, string scheme, double dt, int t_steps) {
    if (system.dims == 1) {
        if (scheme == "exp") {
            exp(system, dt, t_steps);
        } else if (scheme == "imp") {
            imp(system, dt, t_steps);
        } else if (scheme == "crank") {
            crankN(system, dt, t_steps);
        } else {
            cout << "Invalid scheme:" << scheme << endl;
        }
    } else if (system.dims == 2) {
        if (scheme == "exp") {
            exp2D(system, dt, t_steps);
        } else if (scheme == "imp") {
            imp2D(system, dt, t_steps);
        } else {
            cout << "Invalid scheme:" << scheme << endl;
        }
    } else {
        cout << "Invalid system. Dims not 1 or 2." << endl;
    }

}

void exp(System &system, double dt, int t_steps) {
    double alpha = dt / (system.dx * system.dx) * system.D;
    vec u = system.u;
    double t = system.t;
    int n = u.n_elem;
    vec new_u = zeros<vec>(n);
    new_u(0) = u(0);
    new_u(n - 1) = u(n - 1);

    for (int t_i = 0; t_i < t_steps; t_i++) {
        for (int i = 1; i < n - 1; i++) {
            new_u(i) = alpha * u(i - 1) + (1 - 2 * alpha) * u(i) + alpha * u(i + 1) + system.source(i * system.dx, t) * dt;
        }
        u.swap(new_u);
        t += dt;
    }
    system.u = u;
    system.t = t;
}

void imp(System &system, double dt, int t_steps) {
    double alpha = dt / (system.dx * system.dx) * system.D;
    vec u = system.u;
    int n = u.n_elem;
    double t = system.t;
    // Where we store the new u, "x" in the equation Ax = b
    vec x = zeros<vec>(n);
    x(0) = u(0);
    x(n - 1) = u(n - 1);
    // The right hand side "b"
    vec b = zeros(n);
    // Solving the tridiagonal matrix equation for each timestep
    // Note that the first and last elements of the vectors are never used!
    for (int t_i = 0; t_i < t_steps; t_i++) {
        // Initial value for right hand side
        for (int i = 1; i < n - 1; i++) {
            b(i) = u(i) + system.source(i * system.dx, t) * dt;
        }
        b(1) += alpha * u(0); // Accounting for boundrary conditions
        b(n - 2) += alpha * u(n - 1);
        // The tridiagonal matrix "A"
        vec diag = (1 + 2 * alpha) * ones(n);
        vec right = -alpha * ones(n);
        vec left = -alpha * ones(n);
        // Forward substitution. Updating all diagonal elements and right hand side elements, from the second counted row
        double factor;
        for (int i = 2; i < n - 1; i++) {
            factor = left(i - 1) / diag(i - 1); //the factor used for the elimination at each step
            diag(i) = diag(i) - right(i - 1) * factor;
            b(i) = b(i) - b(i - 1) * factor;
        }
        // Backward substitution. With only two diagonals, we can find the coefficients moving backwards
        x(n - 2) = b(n - 2) / diag(n - 2);
        for (int i = n - 3; i > 0; i--) {
            x(i) = (b(i) - right(i) * x(i + 1)) / diag(i);
        }
        u.swap(x);
        t += dt;
    }
    system.u = u;
    system.t = t;
}

void crankN(System &system, double dt, int t_steps) {
    double alpha = dt / (system.dx * system.dx) * system.D;
    vec u = system.u;
    int n = u.n_elem;
    double t = system.t;
    // Where we store the new u, "x" in the equation Ax = b
    vec x = zeros<vec>(n);
    x(0) = u(0);
    x(n - 1) = u(n - 1);
    // The right hand side "b"
    vec b = zeros(n);
    // Solving the tridiagonal matrix equation for each timestep
    // Note that the first and last elements of the vectors are never used!
    for (int t_i = 0; t_i < t_steps; t_i++) {
        // Initial value for right hand side
        for (int i = 1; i < n - 1; i++) {
            // The right hand side. All of the known values
            b(i) = alpha * u(i - 1) + (2 - 2 * alpha) * u(i) + alpha * u(i + 1) + dt * system.source(i * system.dx, t);
        }
        b(1) += alpha * u(0); // Accounting for boundrary conditions
        b(n - 2) += alpha * u(n - 1);
        // The tridiagonal matrix "A"
        vec diag = (2 + 2 * alpha) * ones(n);
        vec right = -alpha * ones(n);
        vec left = -alpha * ones(n);
        // Forward substitution. Updating all diagonal elements and right hand side elements, from the second counted row
        double factor;
        for (int i = 2; i < n - 1; i++) {
            factor = left(i - 1) / diag(i - 1); //the factor used for the elimination at each step
            diag(i) = diag(i) - right(i - 1) * factor;
            b(i) = b(i) - b(i - 1) * factor;
        }
        // Backward substitution. With only two diagonals, we can find the coefficients moving backwards
        x(n - 2) = b(n - 2) / diag(n - 2);
        for (int i = n - 3; i > 0; i--) {
            x(i) = (b(i) - right(i) * x(i + 1)) / diag(i);
        }
        u.swap(x);
        t += dt;
    }
    system.u = u;
    system.t = t;
}

void exp2D(System &system, double dt, int t_steps) {
    double alpha = dt / (system.dx * system.dx) * system.D;
    mat U = system.U;
    int n = U.n_rows;
    double t = system.t;
    mat new_U = U * 1; // copying U

    for (int t_i = 0; t_i < t_steps; t_i++) {
        for (int x = 1; x < n - 1; x++) {
            for (int y = 1; y < n - 1; y++) {
                new_U(x, y) = U(x, y) + alpha * ( U(x + 1, y) + U(x - 1, y) + U(x, y + 1) + U(x, y - 1) - 4 * U(x, y) );
            }
        }
        U.swap(new_U);
        t += dt;
    }
    system.U = U;
    system.t = t;
}

void imp2D(System &system, double dt, int t_steps) {
    double alpha = dt / (system.dx * system.dx) * system.D;
    mat U = system.U;
    int n = U.n_rows;
    double t = system.t;
    mat new_U = U * 1; // copying U
    mat old_U = U * 1;

    int maxIter = 10000;
    double minDiff = 1.0e-12;
    double diff;
    double tempTerm;
    int goodSteps = 0; // the number of times the algorithm converged

    for (int t_i = 0; t_i < t_steps; t_i++) {
        for (int iter = 0; iter < maxIter; iter++) {
            diff = 0;
            for (int x = 1; x < n - 1; x++) {
                for (int y = 1; y < n - 1; y++) {
                    // old_U is the previous state. it is the only known thing in the equation
                    new_U(x, y) = 1.0 / (1 + 4 * alpha) * (alpha * ( U(x + 1, y) + U(x - 1, y) + U(x, y + 1) + U(x, y - 1)) + old_U(x, y) );
                    diff += fabs(U(x, y) - new_U(x, y));
                }
            }
            U.swap(new_U); // new_U is the best guess we have to solve the equation at this time, so we put it in U to find a new best guess
            if (diff < minDiff) {
                goodSteps += 1;
                break;
            }
        }
        old_U = U * 1; // old_U is the definitive state at time t. U holds it also at this point, which will be used in output or next jacobi iteration
        t += dt;
    }
    system.U = U;
    system.t = t;
    if (goodSteps != t_steps) {
        cout << "Did not converge " << t_steps - goodSteps << " times." << endl;
    }
}
