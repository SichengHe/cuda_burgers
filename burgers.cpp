#include <iostream>
#include <vector>
#include<cmath>
#include <fstream>
using namespace std;

void compute_JST_flux(const std::vector<double>& u, std::vector<double>& dudt, double uL, double uR){

    double flux_inv_p1, flux_inv_0, flux_inv, flux_vis;
    double delta_up15, delta_up05, delta_um05;

    double um_loc_m2, um_loc_m1, um_loc_0, um_loc_p1, um_loc_p2;
    double h_p05, h_m05;

    double um1, u0, up1, up2;

    // artificial viscosity
    double coeff_visc_2 = 0.5;
    double coeff_visc_4 = 0.05;

    // Get dimension
    int n = u.size();

    // Flux
    for (int i = 0; i < n; i++) {
        if (i==0){
            um_loc_m2 = uL;
            um_loc_m1 = uL;
            um_loc_0 = u[i];
            um_loc_p1 = u[i+1];
            um_loc_p2 = u[i+2];
            }

        else if (i==1){
            um_loc_m2 = uL;
            um_loc_m1 = u[i-1];
            um_loc_0 = u[i];
            um_loc_p1 = u[i+1];
            um_loc_p2 = u[i+2];
        }
        else if (i==n-2){
            um_loc_m2 = u[i-2];
            um_loc_m1 = u[i-1];
            um_loc_0 = u[i];
            um_loc_p1 = u[i+1];
            um_loc_p2 = uR;
        }
        else if (i==n-1){
            um_loc_m2 = u[i-2];
            um_loc_m1 = u[i-1];
            um_loc_0 = u[i];
            um_loc_p1 = uR;
            um_loc_p2 = uR;
        }
        else{
            um_loc_m2 = u[i-2];
            um_loc_m1 = u[i-1];
            um_loc_0 = u[i];
            um_loc_p1 = u[i+1];
            um_loc_p2 = u[i+2];
        }

        // Left neighbor flux
        um1 = um_loc_m2;
        u0 = um_loc_m1;
        up1 = um_loc_0;
        up2 = um_loc_p1;

        flux_inv_p1 = 0.5 * pow(up1,2);
        flux_inv_0 = 0.5 * pow(u0,2);
        flux_inv = 0.5 * (flux_inv_p1 + flux_inv_0);

        delta_up15 = up2 - up1;
        delta_up05 = up1 - u0;
        delta_um05 = u0 - um1;
        
        flux_vis = coeff_visc_2 * delta_up05 - coeff_visc_4 * (delta_up15 - 2 * delta_up05 + delta_um05);

        h_p05 = flux_inv - flux_vis;

        // Right neighbor flux
        um1 = um_loc_m1;
        u0 = um_loc_0;
        up1 = um_loc_p1;
        up2 = um_loc_p2;

        flux_inv_p1 = 0.5 * pow(up1,2);
        flux_inv_0 = 0.5 * pow(u0,2);
        flux_inv = 0.5 * (flux_inv_p1 + flux_inv_0);


        delta_up15 = up2 - up1;
        delta_up05 = up1 - u0;
        delta_um05 = u0 - um1;
        
        flux_vis = coeff_visc_2 * delta_up05 - coeff_visc_4 * (delta_up15 - 2 * delta_up05 + delta_um05);

        h_m05 = flux_inv - flux_vis;

        // Forward Euler
        dudt[i] = h_p05 - h_m05;
            
    }
}

void update_state(std::vector<double>&u, const std::vector<double>&dudt, double dtdx){
    // Get dimension
    int n = u.size();

    // Flux
    for (int i = 0; i < n; i++) {
        u[i] += dudt[i] * dtdx;
    }
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
    int n = pow(2, 13);
    // int n = 1000;
    double T = 0.8;

    double xL = -1.0;
    double xR = 1.0;
    double L = xR - xL;
    double dx = L / n;

    double CFL_val = 0.8;


    double x_loc;
    double t=0.0;
    double dt;
    double dtdx;

    // Initialization
    std::vector<double> u(n, 0.0);
    std::vector<double> dudt(n, 0.0);

    double uL = 1.0;
    double uR = 0.0;
    for (int i = 0; i < n; i++) {
        x_loc = (i + 0.5) * dx + xL;
        if (x_loc < 0.0){
            u[i] = uL;
        }
        else {
            u[i] = uR;
        }
    }

    // time stepping
    while (t < T){
        // Compute the maximum step size
        dt = compute_dt(u, dx, n, CFL_val);
        dtdx = dt/dx;

        compute_JST_flux(u, dudt, uL, uR);

        update_state(u, dudt, dtdx);

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