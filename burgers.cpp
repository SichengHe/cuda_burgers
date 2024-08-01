#include <iostream>
#include <vector>
#include<cmath>
#include <fstream>
using namespace std;

// // upwind flux
// std::vector<double> laxWendroffFlux(const std::vector<double>& u, double dx, double dt, double velocity) {
//     int n = u.size();
//     std::vector<double> flux(n, 0.0);

//     double alpha = velocity * dt / dx;
//     double alpha2 = alpha * alpha;

//     for (int i = 1; i < n - 1; ++i) {
//         flux[i] = u[i] - 0.5 * alpha * (u[i + 1] - u[i - 1]) 
//                          + 0.5 * alpha2 * (u[i + 1] - 2 * u[i] + u[i - 1]);
//     }

//     return flux;
// }

double compute_JST_flux(double um1, double u0, double up1, double up2){

    double flux_inv_p1, flux_inv_0, flux_inv, flux_vis;
    double delta_up15, delta_up05, delta_um05;

    // artificial viscosity
    double coeff_visc_2 = 0.5;
    double coeff_visc_4 = 0.05;

    // inviscid flux f = 0.5 u^2
    flux_inv_p1 = 0.5 * pow(up1,2);
    flux_inv_0 = 0.5 * pow(u0,2);
    flux_inv = 0.5 * (flux_inv_p1 + flux_inv_0);


    delta_up15 = up2 - up1;
    delta_up05 = up1 - u0;
    delta_um05 = u0 - um1;
    
    flux_vis = coeff_visc_2 * delta_up05 - coeff_visc_4 * (delta_up15 - 2 * delta_up05 + delta_um05);

    return  flux_inv - flux_vis;
}

double compute_dt(const std::vector<double>&u, double dx, int n, double CFL_val){

    double dt = dx / abs(2 * u[0]);
    double dt_new;

    for (int i = 0; i < n; i++){
        dt_new = dx / (abs(2 * u[i]) + 1e-8);
        if (dt > dt_new){
            dt = dt_new;
        }
    }

    return dt * CFL_val;
}

int main() {

    // constants
    int n = 1000;
    double T = 0.8;

    double xL = -1.0;
    double xR = 1.0;
    double L = xR - xL;
    double dx = L / n;

    double CFL_val = 0.8;


    double x_loc;
    double t=0.0;
    double dt;

    double um_loc_m2, um_loc_m1, um_loc_0, um_loc_p1, um_loc_p2;
    double h_p05, h_m05;

    // Initialization
    std::vector<double> u(n, 0.0);
    std::vector<double> um(n, 0.0);

    double uL = 1.0;
    double uR = 0.0;
    for (int i = 0; i < n; i++) {
        x_loc = (i + 0.5) * dx + xL;
        if (x_loc < 0.0){
            um[i] = uL;
        }
        else {
            um[i] = uR;
        }
    }

    // time stepping
    while (t < T){
        // Compute the maximum step size
        dt = compute_dt(um, dx, n, CFL_val);

        // Time integration
        for (int i = 0; i < n; i++) {
            if (i==0){
                um_loc_m2 = uL;
                um_loc_m1 = uL;
                um_loc_0 = um[i];
                um_loc_p1 = um[i+1];
                um_loc_p2 = um[i+2];
                }

            else if (i==1){
                um_loc_m2 = uL;
                um_loc_m1 = um[i-1];
                um_loc_0 = um[i];
                um_loc_p1 = um[i+1];
                um_loc_p2 = um[i+2];
            }
            else if (i==n-2){
                um_loc_m2 = um[i-2];
                um_loc_m1 = um[i-1];
                um_loc_0 = um[i];
                um_loc_p1 = um[i+1];
                um_loc_p2 = uR;
            }
            else if (i==n-1){
                um_loc_m2 = um[i-2];
                um_loc_m1 = um[i-1];
                um_loc_0 = um[i];
                um_loc_p1 = uR;
                um_loc_p2 = uR;
            }
            else{
                um_loc_m2 = um[i-2];
                um_loc_m1 = um[i-1];
                um_loc_0 = um[i];
                um_loc_p1 = um[i+1];
                um_loc_p2 = um[i+2];
            }

            h_p05 = compute_JST_flux(um_loc_m1, um_loc_0, um_loc_p1, um_loc_p2);
            h_m05 = compute_JST_flux(um_loc_m2, um_loc_m1, um_loc_0, um_loc_p1);

            // Forward Euler
            u[i] = um[i] - dt / dx * (h_p05 - h_m05);
                
        }

        // Update the history vector
        std::copy(u.begin(), u.end(), um.begin());

        // Update time
        t = t + dt;
    }

    // Create and open an output file stream
    std::ofstream outFile("./output/output_cpp.txt");

    // Check if the file is opened successfully
    if (!outFile) {
        std::cerr << "Error opening file for writing!" << std::endl;
        return 1;
    }

    // Write the array to the file
    for (int i = 0; i < n; i++) {
        outFile << (i + 0.5) * dx + xL << " " << u[i] << "\n";
    }

    // Close the file
    outFile.close();

    std::cout << "Array saved to output.txt in one column" << std::endl;

    return 0;
}