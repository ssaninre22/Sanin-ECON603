# Patched simulation: clamp hazards, bound theta, smooth R', partial adjustment
import numpy as np
import math
import matplotlib.pyplot as plt
from dataclasses import dataclass, field
from typing import List, Tuple


# ---------- Helpers ----------
def logistic(z: float) -> float:
    if z >= 0:
        ez = math.exp(-z)
        return 1.0 / (1.0 + ez)
    else:
        ez = math.exp(z)
        return ez / (1.0 + ez)

def logistic_pdf(z: float) -> float:
    L = logistic(z)
    return L * (1.0 - L)

def clamp01(x: float) -> float:
    return max(1e-12, min(1.0, x))


# ---------- Data classes ----------
@dataclass
class MatchingParams:
    m: float
    alpha: float
    delta: float
    c: float
    xi: float

@dataclass
class ModelParams:
    x: float
    g0: float
    g1: float
    psi: float
    wbar: float
    rho_path: List[float]
    eta: float
    phi_T: float
    aN: float
    beta: float
    mu: float
    phi_bar: float
    sigma_phi: float
    VG_next: float
    b: float
    YT_bar: float
    beta_T: float

@dataclass
class ShockPaths:
    eps: List[float]
    chi: List[float]
    s_aid: List[float]

@dataclass
class Results:
    L_T: List[float] = field(default_factory=list)
    L_N: List[float] = field(default_factory=list)
    U: List[float] = field(default_factory=list)
    S: List[float] = field(default_factory=list)
    Y_T: List[float] = field(default_factory=list)
    Y_N: List[float] = field(default_factory=list)
    G: List[float] = field(default_factory=list)
    theta_T: List[float] = field(default_factory=list)
    theta_N: List[float] = field(default_factory=list)
    w_T: List[float] = field(default_factory=list)
    w_N: List[float] = field(default_factory=list)


# ---------- Matching functions (clamped) ----------
def q_of_theta(m: float, alpha: float, theta: float) -> float:
    return clamp01(m * theta**(-alpha))

def f_of_theta(m: float, alpha: float, theta: float) -> float:
    return clamp01(m * theta**(1.0 - alpha))

# ---------- Solve theta from free entry + Nash ----------
def solve_theta_given_Rprime(beta: float, b: float, sp: MatchingParams, Rprime: float,
                             theta_min: float = 1e-4, theta_max: float = 20.0) -> Tuple[float, float, float]:
    c, m, alpha, xi, delta = sp.c, sp.m, sp.alpha, sp.xi, sp.delta
    D = 1.0 - beta * (1.0 - delta) * (1.0 - xi)

    def residual(theta):
        theta = max(theta_min, min(theta_max, theta))
        # Free-entry with clamped q
        q = q_of_theta(m, alpha, theta)
        # J from free entry: c = beta q J  => J = c/(beta q)
        J = c / (beta * q)
        # Wage from Nash rule (recursive reduced form)
        w = (1.0 - xi) * b + xi * (Rprime + beta * (1.0 - delta) * J + c * theta)
        # Firm surplus identity J = R' - w + beta(1-d)J
        lhs = J
        rhs = Rprime - w + beta * (1.0 - delta) * J
        return lhs - rhs

    # Bisection within [theta_min, theta_max]
    lo, hi = theta_min, theta_max
    flo, fhi = residual(lo), residual(hi)
    # If same sign, expand hi a bit then fallback
    expanded = False
    if np.sign(flo) == np.sign(fhi):
        for _ in range(10):
            hi *= 2.0
            if hi > 1e6:
                break
            fhi = residual(hi)
            if np.sign(flo) != np.sign(fhi):
                expanded = True
                break
    if np.sign(flo) == np.sign(fhi) and not expanded:
        theta = max(theta_min, min(theta_max, 1.0))
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            fmid = residual(mid)
            if abs(fmid) < 1e-10:
                theta = mid
                break
            if np.sign(fmid) == np.sign(flo):
                lo, flo = mid, fmid
            else:
                hi, fhi = mid, fmid
        else:
            theta = 0.5 * (lo + hi)

    q = q_of_theta(m, alpha, theta)
    J = sp.c / (beta * q)
    w = (1.0 - sp.xi) * b + sp.xi * (Rprime + beta * (1.0 - sp.delta) * J + sp.c * theta)
    return theta, J, w

# ---------- Simulation ----------
def simulate(mp: ModelParams,
             T_params: MatchingParams,
             N_params: MatchingParams,
             shocks: ShockPaths,
             Lbar: float,
             L0_T: float,
             L0_N: float,
             partial_adj: float = 0.6) -> Results:
    T = len(shocks.eps)
    res = Results()

    L_T, L_N = L0_T, L0_N
    theta_T, theta_N = 1.0, 1.0
    w_T, w_N = mp.b + 0.1, mp.b + 0.1

    for t in range(T):
        eps_t = shocks.eps[t]
        chi_t = shocks.chi[t]
        s_t = shocks.s_aid[t]
        rho_t = mp.rho_path[t]

        # 1) T sales and target
        Y_T = mp.YT_bar + mp.beta_T * eps_t
        L_T_target = mp.phi_T * Y_T

        # 2) Enrollment via logit (work value proxy)
        fT = f_of_theta(T_params.m, T_params.alpha, theta_T)
        fN = f_of_theta(N_params.m, N_params.alpha, theta_N)
        V_T = fT * w_T + (1.0 - fT) * mp.b
        V_N = fN * w_N + (1.0 - fN) * mp.b
        V_W = max(V_T, V_N)

        phi_star = s_t + mp.beta * mp.VG_next - V_W
        iota = (phi_star - mp.phi_bar) / mp.sigma_phi
        S_t = mp.mu * logistic(iota)

        # 3) University demand & N fixed point target
        G_t = mp.x * mp.g0 + mp.x * mp.g1 * rho_t * S_t - chi_t * (mp.x ** mp.psi)
        kx = mp.aN / (1.0 - mp.aN * mp.eta * mp.wbar)
        L_N_target = kx * (mp.eta * mp.wbar * L_T + G_t)
        Y_N = L_N_target / mp.aN

        # 4) R' per worker (smoothed denominators)
        LT_eff = 0.5 * L_T + 0.5 * L_T_target
        LN_eff = 0.5 * L_N + 0.5 * L_N_target
        Rprime_T = Y_T / max(LT_eff, 1e-6)
        Rprime_N = Y_N / max(LN_eff, 1e-6)

        # 5) Solve θ and wages
        theta_T, J_T, w_T = solve_theta_given_Rprime(mp.beta, mp.b, T_params, Rprime_T)
        theta_N, J_N, w_N = solve_theta_given_Rprime(mp.beta, mp.b, N_params, Rprime_N)

        # 6) Vacancies to move partially toward target
        q_T = q_of_theta(T_params.m, T_params.alpha, theta_T)
        q_N = q_of_theta(N_params.m, N_params.alpha, theta_N)

        gap_T = (L_T_target - (1.0 - T_params.delta) * L_T) * partial_adj
        gap_N = (L_N_target - (1.0 - N_params.delta) * L_N) * partial_adj

        v_T = max(0.0, gap_T) / max(q_T, 1e-8)
        v_N = max(0.0, gap_N) / max(q_N, 1e-8)

        H_T = q_T * v_T
        H_N = q_N * v_N

        # 7) Update employment
        L_T_next = (1.0 - T_params.delta) * L_T + H_T
        L_N_next = (1.0 - N_params.delta) * L_N + H_N

        # 8) Aggregates
        L_E = L_T + L_N
        U = max(0.0, Lbar - L_E)

        # Store
        res.L_T.append(L_T)
        res.L_N.append(L_N)
        res.U.append(U)
        res.S.append(S_t)
        res.Y_T.append(Y_T)
        res.Y_N.append(Y_N)
        res.G.append(G_t)
        res.theta_T.append(theta_T)
        res.theta_N.append(theta_N)
        res.w_T.append(w_T)
        res.w_N.append(w_N)

        L_T, L_N = L_T_next, L_N_next

    return res


# ---------- Example run with patched code ----------
T = 1500
#eps = np.zeros(T)
eps = np.random.normal(loc=0.0, scale=1.0, size=T)
chi = np.zeros(T)
# shocks windows
#eps[8:13] = -0.2     # 2001
#eps[28:37] = -1.2    # 2008-10
#eps[48:53] = -1.5    # 2020-22
#chi[48:53] = 0.8
s_aid = 0.1 * (-eps)


x = 0.35
rho = np.full(T, 0.8 + 0.15 * x)
#rho[48:53] = 0.2

mp = ModelParams(
    x=x, g0=0.2, g1=0.8, psi=1.4, wbar=1.0,
    rho_path=rho.tolist(),
    eta=0.25 + 0.10 * x,
    phi_T=0.9,
    aN=0.5,
    beta=0.96,
    mu=1.0, phi_bar=0.0, sigma_phi=0.6, VG_next=1.2,
    b=0.4,
    YT_bar=1.0, beta_T=0.8
)

T_params = MatchingParams(m=0.5, alpha=0.5, delta=0.03, c=0.2, xi=0.5)
N_params = MatchingParams(m=0.6, alpha=0.5, delta=0.03, c=0.2, xi=0.6)

shocks = ShockPaths(eps=eps.tolist(), chi=chi.tolist(), s_aid=s_aid.tolist())

Lbar = 1.5
L0_T = 0.6
L0_N = 0.6


res = simulate(mp, T_params, N_params, shocks, Lbar, L0_T, L0_N, partial_adj=0.6)


# ---------- Plots ----------
plt.figure()
plt.plot(res.U)
plt.title("Unemployment")
plt.xlabel("t")
plt.ylabel("U")

plt.figure()
plt.plot(res.L_N, label="L_N")
plt.plot(res.L_T, label="L_T")
plt.title("Employment by Sector")
plt.xlabel("t")
plt.ylabel("Employment")
plt.legend()

plt.figure()
plt.plot(res.S)
plt.title("Enrollment share S_t (logit)")
plt.xlabel("t")
plt.ylabel("S_t")

plt.figure()
plt.plot(res.theta_N, label="theta_N")
plt.plot(res.theta_T, label="theta_T")
plt.title("Tightness by Sector")
plt.xlabel("t")
plt.ylabel("theta")
plt.legend()

plt.show()

