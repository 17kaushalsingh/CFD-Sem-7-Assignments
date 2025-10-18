#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

// Initialize the grid
void initializeGrid(vector<double>& T, int nx, double T0, double T_left, double T_right) {
    for (int i = 0; i < nx; i++) {
        T[i] = T0;
    }
    T[0] = T_left;
    T[nx - 1] = T_right;
}

// FTCS step
void timeStepFTCS(const vector<double>& T, vector<double>& Tnew, int nx, double r, double T_left, double T_right) {
    Tnew = T;
    for (int i = 1; i < nx - 1; i++) {
        Tnew[i] = T[i] + r * (T[i + 1] - 2.0 * T[i] + T[i - 1]);
    }
    Tnew[0] = T_left;
    Tnew[nx - 1] = T_right;
}

// Max error between old and new
double computeError(const vector<double>& T, const vector<double>& Tnew, int nx) {
    double maxErr = 0.0;
    for (int i = 0; i < nx; i++) {
        double diff = fabs(Tnew[i] - T[i]);
        if (diff > maxErr) maxErr = diff;
    }
    return maxErr;
}

int main() {
    // Parameters
    double L = 1.0;
    double dx = 0.01;
    int nx = static_cast<int>(L / dx) + 1;
    double alpha = 2.3e-5;
    double dt = 0.005;
    double tol = 1e-4;
    int maxIt = 2000000;

    // Boundary and initial conditions
    double T_left = 0.0;
    double T_right = 100.0;
    double T0 = 30.0;

    // Stability check
    double r = alpha * dt / (dx * dx);
    if (r > 0.5) {
        cout << "Warning: FTCS may be unstable, r = " << r << endl;
    }

    // Grid
    vector<double> x(nx);
    for (int i = 0; i < nx; i++) {
        x[i] = i * dx;
    }

    // Temperature arrays
    vector<double> T(nx), Tnew(nx);
    initializeGrid(T, nx, T0, T_left, T_right);

    // Time marching
    int it = 0;
    double err = 1e9;
    while (it < maxIt) {
        it++;
        timeStepFTCS(T, Tnew, nx, r, T_left, T_right);
        err = computeError(T, Tnew, nx);
        T = Tnew;
        if (err < tol) break;
    }

    cout << "Converged in " << it << " iterations, error=" 
         << scientific << err << ", r=" << fixed << setprecision(4) << r << endl;

    // Save final x and T to CSV
    ofstream fout("final_profile.csv");
    if (!fout.is_open()) {
        cerr << "Error: Could not open output file\n";
        return 1;
    }

    // First row: x
    for (int i = 0; i < nx; i++) {
        fout << x[i];
        if (i < nx - 1) fout << ",";
    }
    fout << "\n";

    // Second row: T
    for (int i = 0; i < nx; i++) {
        fout << T[i];
        if (i < nx - 1) fout << ",";
    }
    fout << "\n";

    fout.close();
    cout << "Final profile saved to final_profile.csv\n";

    return 0;
}
