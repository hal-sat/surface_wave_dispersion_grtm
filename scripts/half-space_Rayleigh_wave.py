#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import sys
import functools


def object_func(c, args):
    vp = float(args[0])
    vs = float(args[1])

    ans = (c ** 2 / vp ** 2 - 1.0) * (c ** 2 / vs ** 2 - 1.0) - (
        1.0 - 0.5 * c ** 2 / vs ** 2
    ) ** 4
    return ans


def object_func_d(c, args):
    vp = float(args[0])
    vs = float(args[1])

    ans = (
        4 / (vp ** 2 * vs ** 2) * c ** 3
        - 2 * (1.0 / vp ** 2 + 1.0 / vs ** 2) * c
        - 4.0 * c / (vs ** 2) * (0.5 * (c / vs) ** 2 - 1) ** 3
    )
    return ans


def newton_find_root(init_x, obj_fn, obj_fn_d, args, eps=1e-10, itrmax=100):
    x_bf = init_x
    x_af = init_x
    for i in range(itrmax):
        x_af = x_bf - obj_fn(x_bf, args) / obj_fn_d(x_bf, args)
        if np.abs(obj_fn(x_af, args)) < eps:
            break
        if i == itrmax - 1:
            return np.nan
        x_bf = x_af

    return x_af


# vs = float(sys.argv[1])
# vp = float(sys.argv[2])

vs = 180.0  # m/s
vp = 1500.0  # m/s
vs = 180.0  # m/s
vp = 1500.0  # m/s
vs = 350.0  # m/s
vp = 1700.0  # m/s
args = [vp, vs]

trial_c = np.linspace(0, vs * 2, 1000)
obj_func = np.array([object_func(c, args) for c in trial_c])

opt_c = newton_find_root(vs, object_func, object_func_d, args)

z = np.zeros(trial_c.shape[0])

print("Phase velocity is {}.".format(opt_c))

fig, ax = plt.subplots()
ax.plot(trial_c, obj_func)
ax.plot(trial_c, z, "k")
ax.axvline(vs, c="orange", label = "Initial value")
ax.axvline(opt_c, c="red", label="Optimal value")
ax.set_xlabel("Phase velocity (m/s)")
ax.set_ylabel("Secular function")
ax.set_xlim([0, 2 * vs])
ax.set_xlim([vs * 0.8, vs * 1.2])

# fig.savefig("tmp.pdf")

plt.show()
