"""
Simulates a simple model of a colloidal suspension drying on a flexible cantilever.
The stress in the packed solid is assumed to be due to residual
stress arising from gel formation and shrinkage priort to solidification.
"""

import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from params import Params
import matplotlib.pyplot as plt
from math import pi

plt.rcParams.update({
    "font.size": 16,
    "figure.autolayout": True
    })


def incremental_stress(t, x, pars):
    """
    Stress in colloidal gel
    """
    # from fitting
    m = pars.m
    G_c = pars.G_c

    # compute quantities
    j_e = V_e(x)
    h_g = h_0(x, pars) - j_e * pars.t_g / pars.t_d


    tau = j_e / h_g * (t - pars.t_g / pars.t_d)
    alpha = m * h_g / j_e

    f = 1 + 1 / (1 - tau)**2
    sigma = pars.incremental * G_c / alpha * (np.exp(alpha * tau) - 1) * f

    return sigma

def V_e(x):
    """
    Evaporation rate
    """

    a = 0.42714138
    b = -0.3028127

    c = pi * a + b
    return 1 / c * (2 * a / np.sqrt(1 - x**2) + b)

def q_p(x):
    """
    Flux from the packed region
    """

    a = 0.42714138
    b = -0.3028127

    c = pi * a + b

    return -(2 * a * np.arcsin(x) - pi * a + b * x - b) / c

def h_0(x, pars):
    """
    Initial film thickness
    """
    p = 1 / (1 - 1 / np.cosh(pars.Bond**(1/2)))
    return p * (1 - np.cosh(pars.Bond**(1/2) * x) / np.cosh(pars.Bond**(1/2)))

def x_p_ode(t, x_p, pars):
    """
    ODE for x_p
    """
    J_p = pars.phi_0 / pars.phi_p

    return -q_p(x_p) / ((1 - J_p) * h_0(x_p, pars) - V_e(x_p) * t)

def solve_x_p(t_span, x_p_0, pars):
    """
    Solves for the position of the packing front
    """

    sol = solve_ivp(x_p_ode, t_span, [x_p_0], args = (pars, ), atol = 1e-12, rtol = 1e-12)

    return sol

def compute_w_pre(t, x_p, pars):
    """
    Computes the deflection before t_g_salt
    """

    Nt = len(t)
    w_all = np.zeros_like(t)

    for n in range(len(t)):

        X, H, Sigma = compute_pre_quantities(t[n], x_p[n], pars)

        X /= 2
        w_x = cumulative_trapezoid(pars.E / 2 * H * (1 + pars.alpha * H) * Sigma, X, initial=0)
        w = cumulative_trapezoid(w_x, X, initial=0)

        w_all[n] = w[-1]

    return w_all

def compute_pre_quantities(t, x_p, pars):
    """
    Computes the solid/fluid regions, solid/fluid thickness, and
    solid/fluid stress before the salt-induced gelation time
    """

    N = 1000
    xi = np.linspace(0, 1-1e-8, N)

    X_f = xi * x_p
    X_p = x_p + (1 - x_p) * xi

    # height of packed film
    h_p = pars.phi_0 / pars.phi_p * h_0(X_p, pars)
    
    # solid stress
    sigma_f = (pars.sigma_p / pars.G_0) * np.ones(N)

    X = np.r_[-X_p[::-1], -X_f[::-1], X_f, X_p]
    H = np.r_[
        h_p[::-1],
        h_0(-X_f[::-1], pars) - V_e(-X_f[::-1]) * t,
        h_0(X_f, pars) - V_e(X_f) * t,
        h_p
    ]

    Sigma = np.r_[
        sigma_f[::-1], 
        np.zeros(2*N), 
        sigma_f
    ]

    return X, H, Sigma


def compute_post_quantities(t, x_p, x_p_0, pars):
    """
    Computes the solid/gel regions, solid/gel thickness,
    solid/gel stress, and solid/gel solid fraction after the 
    salt-induced gelation time
    """


    J_p = pars.phi_0 / pars.phi_p

    N = 1000
    xi = np.linspace(0, 1-1e-8, N)

    X_f = xi * x_p
    X_p = x_p + (1 - x_p) * xi

    # fluid height and particle fraction at salt-induced gelation
    h_g = (h_0(X_f, pars) - V_e(X_f) * pars.t_g / pars.t_d)

    # film height in fluid region
    h_f = (h_0(X_f, pars) - V_e(X_f) * t)

    # elastic deformation (gelled to current)
    J_e = h_f / h_g
    
    # particle fraction in the gel
    phi_s = pars.phi_0 * h_0(X_f, pars) / h_f

    # gel stress
    sigma_f = incremental_stress(t, X_f, pars)

    # solid stress
    sigma_s = (pars.sigma_p / pars.G_0) * np.ones(N)


    # build global arrays
    X = np.r_[-X_p[::-1], -X_f[::-1], X_f, X_p]
    H = np.r_[
        J_p * h_0(-X_p[::-1], pars),
        h_f[::-1],
        h_f,
        J_p * h_0(X_p, pars)
    ]

    Sigma = np.r_[
        sigma_s[::-1], 
        sigma_f[::-1],
        sigma_f,
        sigma_s
    ]

    Phi_s = np.r_[
        pars.phi_p * np.ones(N), 
        phi_s[::-1],
        phi_s,
        pars.phi_p * np.ones(N)

    ]


    return X, H, Sigma, Phi_s


def compute_w_post(t, x_p, pars):
    """
    Computes the deflection after t_g_salt.  Here,
    G_fun is a function that computes the shear modulus
    """


    Nt = len(t)
    w_all = np.zeros_like(t)


    for n in range(Nt):

        X, H, Sigma, _ = compute_post_quantities(t[n], x_p[n], x_p[0], pars)

        X /= 2

        w_x = cumulative_trapezoid(pars.E / 2 * H * (1 + pars.alpha * H) * Sigma, X, initial=0)
        w = cumulative_trapezoid(w_x, X, initial=0)

        w_all[n] = w[-1]

    return w_all



def main():
    """
    Solve the model and plot the deflection and solidification front
    """

    fig_w, ax_w = plt.subplots()
    fig_x, ax_x = plt.subplots()

    # loop over RH values
    for rh in [20, 40, 60, 80]:

        # create params assuming c_s = 1 M
        pars = Params(rh)

        # turn off the stress in the central colloidal gel after t_g_s
        # pars.incremental = 0

        # solve pre-gelation model
        sol_pre = solve_x_p((0, pars.t_g / pars.t_d), 1-1e-6, pars)
        w_pre = compute_w_pre(sol_pre.t, sol_pre.y[0], pars)

        # solve post-gelation model
        sol_post = solve_x_p((pars.t_g / pars.t_d, pars.t_end), sol_pre.y[0, -1], pars)
        w_post = compute_w_post(sol_post.t, sol_post.y[0], pars)

        t_all = np.r_[sol_pre.t, sol_post.t]
        w_all = np.r_[w_pre, w_post]
        x_p_all = np.r_[sol_pre.y[0], sol_post.y[0]]

        """
        make some plots
        """

        # deflection
        plot = ax_w.plot
        p = plot(t_all, w_all, lw = 2, label = f'RH = {rh}')
        plot(sol_pre.t[-1], w_pre[-1], 's', ms = 6, mec = 'k', mfc = p[0].get_color(), zorder = 10)
        ax_w.set_xlabel('$t / t_d$')
        ax_w.set_ylabel('$\Delta w / h_s$')
        ax_w.legend()
        ax_w.set_ylim((0, 20))

        # packing front
        ax_x.plot(t_all, x_p_all, lw = 2, label = f'RH = {rh}', c =p[0].get_color())
        ax_x.legend()
        ax_x.set_xlabel('$t / t_d$')
        ax_x.set_ylabel('$x_p / L$')

    
    plt.show()

main()