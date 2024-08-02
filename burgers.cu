#include <algorithm>
#include <cassert>
#include <iostream>
#include <vector>
#include<cmath>
#include <fstream>
using namespace std;

__global__ void compute_JST_flux(const double *__restrict u, double *__restrict dudt,
                          double uL, double uR, int n) {

    double flux_inv_p1, flux_inv_0, flux_inv, flux_vis;
    double delta_up15, delta_up05, delta_um05;

    double um_loc_m2, um_loc_m1, um_loc_0, um_loc_p1, um_loc_p2;
    double h_p05, h_m05;

    double um1, u0, up1, up2;

    // artificial viscosity
    double coeff_visc_2 = 0.5;
    double coeff_visc_4 = 0.05;

    // Calculate global thread ID
    int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

    if (tid < n){
        if (tid==0){
            um_loc_m2 = uL;
            um_loc_m1 = uL;
            um_loc_0 = u[tid];
            um_loc_p1 = u[tid+1];
            um_loc_p2 = u[tid+2];
            }

        else if (tid==1){
            um_loc_m2 = uL;
            um_loc_m1 = u[tid-1];
            um_loc_0 = u[tid];
            um_loc_p1 = u[tid+1];
            um_loc_p2 = u[tid+2];
        }
        else if (tid==n-2){
            um_loc_m2 = u[tid-2];
            um_loc_m1 = u[tid-1];
            um_loc_0 = u[tid];
            um_loc_p1 = u[tid+1];
            um_loc_p2 = uR;
        }
        else if (tid==n-1){
            um_loc_m2 = u[tid-2];
            um_loc_m1 = u[tid-1];
            um_loc_0 = u[tid];
            um_loc_p1 = uR;
            um_loc_p2 = uR;
        }
        else{
            um_loc_m2 = u[tid-2];
            um_loc_m1 = u[tid-1];
            um_loc_0 = u[tid];
            um_loc_p1 = u[tid+1];
            um_loc_p2 = u[tid+2];
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
        dudt[tid] = h_p05 - h_m05;

    }
}

__global__ void update_state(double *__restrict u, const double *__restrict dudt,double dtdx, int n){

    // Calculate global thread ID
    int tid = (blockIdx.x * blockDim.x) + threadIdx.x;

    // Flux
    // Boundary check
    if (tid < n) u[tid] = u[tid] + dudt[tid] * dtdx;

    // if (tid == 600) printf("%f",u[tid]);
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
    int n = pow(2, 16);
    // int n = pow(2, 10);
    int bytes = sizeof(double) * n;
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

    // Initialization on CPU
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

    // Allocate mem on device 
    double *d_u, *d_dudt;
    cudaMalloc(&d_u, bytes);
    cudaMalloc(&d_dudt, bytes);

    // Copy data from host to the device
    cudaMemcpy(d_u, u.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_dudt, dudt.data(), bytes, cudaMemcpyHostToDevice);

    // Threads per CTA (1024)
    int NUM_THREADS = 1 << 10;

    // CTAs per Grid
    // We need to launch at LEAST as many threads as we have elements
    // This equation pads an extra CTA to the grid if N cannot evenly be divided
    // by NUM_THREADS (e.g. N = 1025, NUM_THREADS = 1024)
    int NUM_BLOCKS = (n + NUM_THREADS - 1) / NUM_THREADS;

    // time stepping
    while (t < T){
        // Compute the maximum step size
        // dt = compute_dt(u, dx, n, CFL_val);
        // dt = 0.001 * CFL_val; // HACK
        dt = 1.5 * pow(10.0,-5) * CFL_val; // HACK

        dtdx = dt/dx;

        compute_JST_flux<<<NUM_BLOCKS, NUM_THREADS>>>(d_u, d_dudt, uL, uR, n);

        update_state<<<NUM_BLOCKS, NUM_THREADS>>>(d_u, d_dudt, dtdx, n);

        // Update time
        t = t + dt;
    }

    cudaMemcpy(u.data(), d_u, bytes, cudaMemcpyDeviceToHost);

    // Create and open an output file stream
    std::ofstream outFile("./output/output_cuda.txt");

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