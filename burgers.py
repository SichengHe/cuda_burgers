import numpy as np

# Inviscid Buergers equation solution.
# Jameson-Schmidt-Turkel (JST) scheme for inviscid flux computation 
# (with artificial viscosity), forward Euler time integration scheme. 

def compute_JST_flux(um1, u0, up1, up2):

    # inviscid flux f = 0.5 u^2
    flux_inv_p1 = 0.5 * up1 ** 2
    flux_inv_0 = 0.5 * u0 ** 2
    flux_inv = 0.5 * (flux_inv_p1 + flux_inv_0)

    # artificial viscosity
    coeff_visc_2 = 0.5
    coeff_visc_4 = 0.05

    delta_up15 = up2 - up1
    delta_up05 = up1 - u0
    delta_um05 = u0 - um1
    
    flux_vis = coeff_visc_2 * delta_up05 - coeff_visc_4 * (delta_up15 - 2 * delta_up05 + delta_um05)

    return  flux_inv - flux_vis

def compute_dt(u, dx, n, CFL_val):

    dt = dx / abs(2 * u[0])
    
    for i in range(n):
        dt_new = dx / (abs(2 * u[i]) + 1e-8)
        if dt > dt_new:
            dt = dt_new

    return dt * CFL_val

n = 2**12
# n = 65536
L = 1.0
T = 0.8

xL = -1.0
xR = 1.0
L = xR - xL
dx = L / n

CFL_val = 0.8


# Initialization
u = np.zeros(n)
um = np.zeros(n)

uL = 1.0
uR = 0.0
for i in range(n):

    x_loc = (i + 0.5) * dx + xL
    if x_loc < 0:
        um[i] = uL
    else:
        um[i] = uR


t = 0.0
while t < T:

    # Compute the maximum step size
    dt = compute_dt(um, dx, n, CFL_val)

    # Time iteration
    for i in range(n):

        if i==0:

            um_loc_m2 = uL
            um_loc_m1 = uL
            um_loc_0 = um[i]
            um_loc_p1 = um[i+1]
            um_loc_p2 = um[i+2]
        
        elif i==1: 

            um_loc_m2 = uL
            um_loc_m1 = um[i-1]
            um_loc_0 = um[i]
            um_loc_p1 = um[i+1]
            um_loc_p2 = um[i+2]

        elif (i==n-2):

            um_loc_m2 = um[i-2]
            um_loc_m1 = um[i-1]
            um_loc_0 = um[i]
            um_loc_p1 = um[i+1]
            um_loc_p2 = uR
            
        elif (i==n-1):

            um_loc_m2 = um[i-2]
            um_loc_m1 = um[i-1]
            um_loc_0 = um[i]
            um_loc_p1 = uR
            um_loc_p2 = uR

        else:
            um_loc_m2 = um[i-2]
            um_loc_m1 = um[i-1]
            um_loc_0 = um[i]
            um_loc_p1 = um[i+1]
            um_loc_p2 = um[i+2]

        h_p05 = compute_JST_flux(um_loc_m1, um_loc_0, um_loc_p1, um_loc_p2)
        h_m05 = compute_JST_flux(um_loc_m2, um_loc_m1, um_loc_0, um_loc_p1)

        # Forward Euler
        u[i] = um[i] - dt / dx * (h_p05 - h_m05)

    # Update the history vector
    um[:] = u[:]

    # Update time
    t = t + dt

x_arr = np.linspace(xL, xR, n)

data = np.zeros((n, 2))
data[:, 0] = x_arr[:]
data[:, 1] = u[:]

np.savetxt("output/output_python.txt", data)

