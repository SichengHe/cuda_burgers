
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

fig, ax = plt.subplots()
fig.set_size_inches(12, 6)

data_cpp = np.loadtxt("./output/output_cpp.txt")
data_py = np.loadtxt("./output/output_python.txt")

ax.plot(data_cpp[:, 0], data_cpp[:, 1], '-o', label=r"$\text{C++ JST}$", clip_on=False)
ax.plot(data_py[:, 0], data_py[:, 1], '-o', label=r"$\text{python JST}$", clip_on=False)
ax.plot(data_cpp[:, 0], get_ref_sol(data_cpp[:, 0], T), '-', label=r"$\text{analytic}$", clip_on=False)
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$u$", rotation="horizontal", ha="right")
niceplots.adjust_spines(ax)
niceplots.label_line_ends(ax)




plt.show()