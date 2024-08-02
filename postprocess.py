
import matplotlib.pyplot as plt
import niceplots
import numpy as np

T = 0.8
def get_ref_sol(x, T):

    u = np.zeros_like(x)
    s = 0.5 # by Rankine-Hugnoniot condition
    for i in range(len(x)):
        if x[i]< s * T:
            u[i] = 1.0
        else: 
            u[i] = 0.0

    return u

plt.style.use(niceplots.get_style())

fig, ax = plt.subplots(2, 1)
fig.set_size_inches(12, 6)

data_cpp = np.loadtxt("./output/output_cpp.txt")
data_cuda = np.loadtxt("./output/output_cuda.txt")
data_py = np.loadtxt("./output/output_python.txt")

# Regression
ax[0].plot(data_py[:, 0], data_py[:, 1], '-o', label=r"$\text{python}$", clip_on=False)
ax[0].plot(data_cpp[:, 0], data_cpp[:, 1], '-o', label=r"$\text{C++}$", clip_on=False)
ax[0].plot(data_cuda[:, 0], data_cuda[:, 1], '-o', label=r"$\text{CUDA}$", clip_on=False)

ax[0].plot(data_cpp[:, 0], get_ref_sol(data_cpp[:, 0], T), '-', label=r"$\text{analytic}$", clip_on=False)
ax[0].set_xlabel(r"$x$")
ax[0].set_ylabel(r"$u$", rotation="horizontal", ha="right")
niceplots.adjust_spines(ax[0])
niceplots.label_line_ends(ax[0])

# Time
n_py = [2**10, 2**11, 2**12]
n_cpp = [2**13, 2**14, 2**15]
n_cuda = [2**14, 2**15, 2**16]

t_py = [2.249, 8.394, 33.682]
t_cpp = [3.372, 12.317, 47.389]
t_cuda = [4.467, 8.550, 12.729]

ax[1].plot(n_py, t_py, '-o', label=r"$\text{python}$", clip_on=False)
ax[1].plot(n_cpp, t_cpp, '-o', label=r"$\text{C++}$", clip_on=False)
ax[1].plot(n_cuda, t_cuda, '-o', label=r"$\text{CUDA}$", clip_on=False)

ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlabel(r"$x$")
ax[1].set_ylabel(r"$\text{Time, s}$", rotation="horizontal", ha="right")
niceplots.adjust_spines(ax[1])
niceplots.label_line_ends(ax[1])

plt.savefig('results.pdf', bbox_inches='tight')


plt.show()