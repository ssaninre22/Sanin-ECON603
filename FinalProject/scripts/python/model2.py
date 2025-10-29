# ------------------------------------------------------------
# Search–Matching with EV1–Logit Enrollment
# Unemployment as stocks, Monte Carlo, and county comparisons
# ------------------------------------------------------------
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Tuple, Dict, Callable
import matplotlib.pyplot as plt


# =========================
# Helpers
# =========================
def logistic(z: float) -> float:
    if z >= 0:
        ez = np.exp(-z)
        return 1.0 / (1.0 + ez)
    else:
        ez = np.exp(z)
        return ez / (1.0 + ez)

def clamp01(x: float) -> float:
    return max(1e-12, min(1.0, x))

def q_of_theta(m: float, alpha: float, theta: float) -> float:
    """Vacancy fill prob per period in (0,1]."""
    return clamp01(m * theta**(-alpha))

def f_of_theta(m: float, alpha: float, theta: float) -> float:
    """Job-finding prob per period in (0,1]."""
    return clamp01(m * theta**(1.0 - alpha))

# =========================
# Parameters
# =========================
@dataclass
class MatchingParams:
    m: float
    alpha: float
    delta: float
    c: float
    xi: float

@dataclass
class ModelParams:
    # exposure & university demand
    x: float
    g0: float
    g1: float
    psi: float
    wbar: float
    eta: float           # MPC on N (can depend on x)
    phi_T: float         # unit labor in T
    aN: float            # unit labor in N (can depend on x)
    beta: float
    # enrollment (logit)
    mu: float
    phi_bar: float
    sigma_phi: float
    VG_next: float
    # worker outside option
    b: float
    # tradables demand
    YT_bar: float
    beta_T: float

# =========================
# Free entry + Nash (solves for θ, J, w)
# =========================
def solve_theta_given_Rprime(beta: float, b: float, sp: MatchingParams, Rprime: float,
                             theta_min: float = 1e-4, theta_max: float = 20.0,
                             max_iter: int = 200) -> Tuple[float, float, float]:
    """
    Solve the system using a residual on the firm surplus identity:
      J = R' - w + β(1-δ)J,
      with free entry (c = β q(θ) J) and Nash wage.
    Hazards are clamped to (0,1].
    """
    c, m, alpha, xi, delta = sp.c, sp.m, sp.alpha, sp.xi, sp.delta

    def residual(theta: float) -> float:
        theta = max(theta_min, min(theta_max, theta))
        q = q_of_theta(m, alpha, theta)
        J = c / (beta * q)  # from free entry
        w = (1.0 - xi) * b + xi * (Rprime + beta * (1.0 - delta) * J + c * theta)
        lhs = J
        rhs = Rprime - w + beta * (1.0 - delta) * J
        return lhs - rhs

    lo, hi = theta_min, theta_max
    flo, fhi = residual(lo), residual(hi)
    if np.sign(flo) == np.sign(fhi):
        # try expanding hi a bit
        for _ in range(10):
            hi *= 2.0
            fhi = residual(hi)
            if np.sign(flo) != np.sign(fhi):
                break
        else:
            # fallback
            theta = 1.0
            q = q_of_theta(m, alpha, theta)
            J = c / (beta * q)
            w = (1.0 - xi) * b + xi * (Rprime + beta * (1.0 - delta) * J + c * theta)
            return theta, J, w

    for _ in range(max_iter):
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
    J = c / (beta * q)
    w = (1.0 - sp.xi) * b + sp.xi * (Rprime + beta * (1.0 - sp.delta) * J + sp.c * theta)
    return theta, J, w

# =========================
# One simulation path (with unemployment stocks)
# =========================
def simulate_one_path(
    T: int,
    rng: np.random.Generator,
    mp: ModelParams,
    T_params: MatchingParams,
    N_params: MatchingParams,
    # labor forces by sector (constant)
    N_T: float,
    N_N: float,
    # initial unemployment rates by sector (fractions of sector labor force)
    u0_T: float = 0.05,
    u0_N: float = 0.05,
    # presence & closures
    rho_path: np.ndarray = None,
    chi_path: np.ndarray = None,
    # shocks
    eps_path: np.ndarray = None,
    s_aid_path: np.ndarray = None,
    # aggressiveness: how much target-gap nudges R' (tightness push)
    kappa_gap: float = 0.0,
) -> pd.DataFrame:
    """
    Simulates a single economy path with proper unemployment stocks:
      L_{s,t+1} = (1-δ_s)L_{s,t} + f_s(θ_{s,t}) U_{s,t}
      U_{s,t+1} = (1-f_s(θ_{s,t})) U_{s,t} + δ_s L_{s,t}
    Enrollment uses EV1/logit with a myopic work-value proxy.
    θ is solved from free entry + Nash each period, with an optional gap-tilt (kappa_gap).
    """
    # Defaults for exogenous paths
    if rho_path is None:
        rho_path = np.ones(T)  # full presence
    if chi_path is None:
        chi_path = np.zeros(T) # no closure
    if eps_path is None:
        eps_path = rng.normal(0.0, 1.0, size=T)  # N(0,1)
    if s_aid_path is None:
        s_aid_path = 0.0 * eps_path              # flat aid (or 0.1*(-eps_path) if desired)

    # Initialize LM stocks
    U_T = u0_T * N_T
    U_N = u0_N * N_N
    L_T = (1.0 - u0_T) * N_T
    L_N = (1.0 - u0_N) * N_N

    # Init prices/values
    theta_T, theta_N = 1.0, 1.0
    w_T, w_N = mp.b + 0.1, mp.b + 0.1

    # Results container
    out = {
        "t": [], "Y_T": [], "Y_N": [], "L_T": [], "L_N": [], "U_T": [], "U_N": [],
        "u_rate": [], "S": [], "G": [], "theta_T": [], "theta_N": [], "w_T": [], "w_N": []
    }

    # constant multiplier for N fixed point
    kx = mp.aN / (1.0 - mp.aN * mp.eta * mp.wbar)

    for t in range(T):
        eps_t = eps_path[t]
        rho_t = rho_path[t]
        chi_t = chi_path[t]
        s_t   = s_aid_path[t]

        # Tradables sales
        Y_T = mp.YT_bar + mp.beta_T * eps_t

        # Enrollment (logit) with work value based on current θ and wages
        fT = f_of_theta(T_params.m, T_params.alpha, theta_T)
        fN = f_of_theta(N_params.m, N_params.alpha, theta_N)
        V_T = fT * w_T + (1.0 - fT) * mp.b
        V_N = fN * w_N + (1.0 - fN) * mp.b
        V_W = max(V_T, V_N)

        phi_star = s_t + mp.beta * mp.VG_next - V_W
        iota = (phi_star - mp.phi_bar) / mp.sigma_phi
        S_t = mp.mu * logistic(iota)

        # University demand
        G_t = mp.x * mp.g0 + mp.x * mp.g1 * rho_t * S_t - chi_t * (mp.x ** mp.psi)

        # N fixed point target employment based on current L_T (not forced, used only to tilt R')
        L_T_target = mp.phi_T * Y_T
        L_N_target = kx * (mp.eta * mp.wbar * L_T + G_t)
        Y_N = L_N_target / mp.aN  # implied by target

        # Effective revenue-per-worker with smoothing and optional gap tilt
        LT_eff = 0.5 * L_T + 0.5 * L_T_target
        LN_eff = 0.5 * L_N + 0.5 * L_N_target
        Rprime_T = Y_T / max(LT_eff, 1e-6)
        Rprime_N = Y_N / max(LN_eff, 1e-6)

        # Optional: nudge R' upward if there is a large desired net addition (encourages tighter markets)
        gap_T = max(0.0, L_T_target - ((1.0 - T_params.delta) * L_T + fT * U_T))
        gap_N = max(0.0, L_N_target - ((1.0 - N_params.delta) * L_N + fN * U_N))
        if kappa_gap > 0.0:
            denom_T = max(L_T + U_T, 1e-6)
            denom_N = max(L_N + U_N, 1e-6)
            Rprime_T += kappa_gap * gap_T / denom_T
            Rprime_N += kappa_gap * gap_N / denom_N

        # Solve θ and wages
        theta_T, J_T, w_T = solve_theta_given_Rprime(mp.beta, mp.b, T_params, Rprime_T)
        theta_N, J_N, w_N = solve_theta_given_Rprime(mp.beta, mp.b, N_params, Rprime_N)

        # Realized hires via matching from unemployment stocks
        fT = f_of_theta(T_params.m, T_params.alpha, theta_T)
        fN = f_of_theta(N_params.m, N_params.alpha, theta_N)
        H_T = fT * U_T
        H_N = fN * U_N

        # Update stocks
        L_T_next = (1.0 - T_params.delta) * L_T + H_T
        U_T_next = (1.0 - fT) * U_T + T_params.delta * L_T

        L_N_next = (1.0 - N_params.delta) * L_N + H_N
        U_N_next = (1.0 - fN) * U_N + N_params.delta * L_N

        u_rate = (U_T + U_N) / (N_T + N_N)

        # Save
        out["t"].append(t)
        out["Y_T"].append(Y_T)
        out["Y_N"].append(Y_N)
        out["L_T"].append(L_T)
        out["L_N"].append(L_N)
        out["U_T"].append(U_T)
        out["U_N"].append(U_N)
        out["u_rate"].append(u_rate)
        out["S"].append(S_t)
        out["G"].append(G_t)
        out["theta_T"].append(theta_T)
        out["theta_N"].append(theta_N)
        out["w_T"].append(w_T)
        out["w_N"].append(w_N)

        # advance
        L_T, U_T = L_T_next, U_T_next
        L_N, U_N = L_N_next, U_N_next

    return pd.DataFrame(out)

# =========================
# Monte Carlo wrapper
# =========================
def run_monte_carlo(
    num_runs: int,
    T_total: int,
    T_keep: int,
    mp_builder: Callable[[float, int], ModelParams],
    exposure_list: List[float],
    T_params: MatchingParams,
    N_params: MatchingParams,
    # sector labor forces (can differ by county if you want)
    N_T: float = 0.8,
    N_N: float = 0.8,
    # presence/closures
    rho_rule: Callable[[float, int], np.ndarray] = lambda x, T: np.ones(T),
    chi_rule: Callable[[float, int], np.ndarray] = lambda x, T: np.zeros(T),
    # stipends rule (default: mildly countercyclical)
    s_rule: Callable[[np.ndarray], np.ndarray] = lambda eps: 0.1 * (-eps),
    # gap tilt
    kappa_gap: float = 0.0,
    seed: int = 1234
) -> Dict[float, List[pd.DataFrame]]:
    """
    For each exposure in exposure_list, runs num_runs simulations of length T_total,
    keeps the last T_keep periods, and returns a dict: exposure -> list of DataFrames.
    """
    assert T_keep <= T_total, "T_keep must be <= T_total."
    rng_master = np.random.default_rng(seed)
    out: Dict[float, List[pd.DataFrame]] = {}

    for x in exposure_list:
        out[x] = []
        for _ in range(num_runs):
            # independent seeds per run
            rng = np.random.default_rng(rng_master.integers(1, 2**31-1))
            # build model params for this exposure
            mp = mp_builder(x, T_total)
            # shocks
            eps_path = rng.normal(0.0, 1.0, size=T_total)
            chi_path = chi_rule(x, T_total)
            rho_path = rho_rule(x, T_total)
            s_aid = s_rule(eps_path)

            df = simulate_one_path(
                T=T_total, rng=rng, mp=mp,
                T_params=T_params, N_params=N_params,
                N_T=N_T, N_N=N_N,
                u0_T=0.05, u0_N=0.05,
                rho_path=rho_path, chi_path=chi_path,
                eps_path=eps_path, s_aid_path=s_aid,
                kappa_gap=kappa_gap
            )
            # keep last T_keep periods
            out[x].append(df.iloc[-T_keep:].reset_index(drop=True))
    return out

# =========================
# Comparison helper (two counties)
# =========================
def compare_two_counties(
    num_runs: int,
    T_total: int,
    T_keep: int,
    x_low: float,
    x_high: float,
    mp_builder: Callable[[float, int], ModelParams],
    T_params: MatchingParams,
    N_params: MatchingParams,
    N_T: float = 0.8,
    N_N: float = 0.8,
    rho_rule: Callable[[float, int], np.ndarray] = lambda x, T: np.ones(T),
    chi_rule: Callable[[float, int], np.ndarray] = lambda x, T: np.zeros(T),
    s_rule: Callable[[np.ndarray], np.ndarray] = lambda eps: 0.1 * (-eps),
    kappa_gap: float = 0.0,
    seed: int = 2024
) -> Tuple[List[pd.DataFrame], List[pd.DataFrame]]:
    """
    Runs Monte Carlo for two exposures (low vs high) and returns lists of final slices (DataFrames).
    """
    exposure_list = [x_low, x_high]
    res = run_monte_carlo(
        num_runs=num_runs, T_total=T_total, T_keep=T_keep,
        mp_builder=mp_builder, exposure_list=exposure_list,
        T_params=T_params, N_params=N_params,
        N_T=N_T, N_N=N_N,
        rho_rule=rho_rule, chi_rule=chi_rule, s_rule=s_rule,
        kappa_gap=kappa_gap, seed=seed
    )
    return res[x_low], res[x_high]


def plot_single_run(df, title_suffix=""):
    fig, ax = plt.subplots()
    ax.plot(df["u_rate"])
    ax.set_title(f"Unemployment rate{title_suffix}")
    ax.set_xlabel("t"); ax.set_ylabel("u_rate")

    fig, ax = plt.subplots()
    ax.plot(df["L_N"], label="L_N")
    ax.plot(df["L_T"], label="L_T")
    ax.set_title(f"Employment by sector{title_suffix}")
    ax.set_xlabel("t"); ax.set_ylabel("Employment"); ax.legend()

    fig, ax = plt.subplots()
    ax.plot(df["S"])
    ax.set_title(f"Enrollment share S_t (logit){title_suffix}")
    ax.set_xlabel("t"); ax.set_ylabel("S_t")

    fig, ax = plt.subplots()
    ax.plot(df["theta_N"], label="theta_N")
    ax.plot(df["theta_T"], label="theta_T")
    ax.set_title(f"Tightness by sector{title_suffix}")
    ax.set_xlabel("t"); ax.set_ylabel("theta"); ax.legend()

    plt.show()
    
    
# =========================
# Example: set builders & defaults (edit as needed)
# =========================
def make_model_params(x: float, T: int) -> ModelParams:
    """
    Example builder: exposure enters η(x) and aN(x).
    You can tailor g0,g1,psi,VG_next, etc. here.
    """
    eta_x = 0.25 + 0.10 * x
    aN_x  = 0.35 + 0.25 * x  # ↑ exposure ⇒ more services intensity (more labor per output)
    return ModelParams(
        x=x, g0=0.2, g1=0.8, psi=1.4, wbar=1.0,
        eta=eta_x, phi_T=0.9, aN=aN_x,
        beta=0.96,
        mu=1.0, phi_bar=0.0, sigma_phi=0.6, VG_next=1.2,
        b=0.4, YT_bar=1.0, beta_T=0.8
    )

def rho_full_presence(x: float, T: int) -> np.ndarray:
    """Full local presence throughout."""
    return np.ones(T)

def chi_none(x: float, T: int) -> np.ndarray:
    """No closure shock."""
    return np.zeros(T)

def s_countercyclical(eps: np.ndarray) -> np.ndarray:
    """Aid mildly countercyclical."""
    return 0.1 * (-eps)

# =========================
# RUN (examples)
# =========================

if __name__ == "__main__":
    # Matching parameters
    T_params = MatchingParams(m=0.5, alpha=0.5, delta=0.03, c=0.2, xi=0.5)
    N_params = MatchingParams(m=0.6, alpha=0.5, delta=0.03, c=0.2, xi=0.6)

    # Monte Carlo settings
    NUM_RUNS = 1000       # <-- your 2,000 simulations
    T_TOTAL  = 200        # total periods per run (e.g. 100 years)
    T_KEEP   = 120        # keep last 120 (e.g. 30 years quarterly)
    
    # Simple path
    rng = np.random.default_rng(0)
    df_one = simulate_one_path(
        T=200, rng=rng, mp=make_model_params(0.35, 200),
        T_params=T_params, N_params=N_params,
        N_T=0.8, N_N=0.8,
        rho_path=rho_full_presence(0.35, 200),
        chi_path=chi_none(0.35, 200),
        eps_path=rng.normal(0,1,200),
        s_aid_path=None,   # uses default (flat) if None
        kappa_gap=0.0
    )
    plot_single_run(df_one, title_suffix=" (x=0.35)")



    # Compare low vs high exposure, full presence, no closures
    low_runs, high_runs = compare_two_counties(
        num_runs=NUM_RUNS,
        T_total=T_TOTAL,
        T_keep=T_KEEP,
        x_low=0.05,
        x_high=0.20,
        mp_builder=make_model_params,
        T_params=T_params,
        N_params=N_params,
        N_T=0.8, N_N=0.8,
        rho_rule=rho_full_presence,
        chi_rule=chi_none,
        s_rule=s_countercyclical,
        kappa_gap=0.0,
        seed=42
    )

    # Example: stack Monte Carlo outputs for quick analysis
    df_low  = pd.concat(low_runs,  keys=range(NUM_RUNS), names=["run","t"]).reset_index(level=0)
    df_high = pd.concat(high_runs, keys=range(NUM_RUNS), names=["run","t"]).reset_index(level=0)


    def summarize_mc(df, col):
        # returns mean, p10, p90 by time index
        grouped = df.groupby(df.index)[col]
        mean = grouped.mean()
        p10  = grouped.quantile(0.10)
        p90  = grouped.quantile(0.90)
        return mean, p10, p90
    # The two DataFrames contain the last 120 periods of each run:
    # columns: ['run','t','Y_T','Y_N','L_T','L_N','U_T','U_N','u_rate','S','G','theta_T','theta_N','w_T','w_N']
    # You can now plot, compare means, quantiles, SDiD-style effects, etc.
    # Example summary (commented out to avoid printing huge output):
    # print(df_low.groupby("t")["u_rate"].mean())
    # print(df_high.groupby("t")["u_rate"].mean())

    # 2a) Unemployment rate comparison (mean and 10–90% bands)
    mean_low, p10_low, p90_low   = summarize_mc(df_low,  "u_rate")
    mean_high, p10_high, p90_high = summarize_mc(df_high, "u_rate")

    fig, ax = plt.subplots()
    ax.plot(mean_low.values,  label="Low exposure (mean)")
    ax.fill_between(range(len(p10_low)),  p10_low.values,  p90_low.values,  alpha=0.2)

    ax.plot(mean_high.values, label="High exposure (mean)")
    ax.fill_between(range(len(p10_high)), p10_high.values, p90_high.values, alpha=0.2)

    ax.set_title("Unemployment rate: Low vs High Exposure (Monte Carlo mean & 10–90% band)")
    ax.set_xlabel("t (last 120 periods kept)")
    ax.set_ylabel("u_rate")
    ax.legend()
    plt.show()

    # 2b Enrollment and tightness comparisons (same pattern)
    for col, label in [("S","Enrollment S_t"), ("theta_N","Non-tradables tightness")]:
        mean_low, p10_low, p90_low   = summarize_mc(df_low,  col)
        mean_high, p10_high, p90_high = summarize_mc(df_high, col)

    fig, ax = plt.subplots()
    ax.plot(mean_low.values,  label=f"Low exposure (mean)")
    ax.fill_between(range(len(p10_low)),  p10_low.values,  p90_low.values,  alpha=0.2)

    ax.plot(mean_high.values, label=f"High exposure (mean)")
    ax.fill_between(range(len(p10_high)), p10_high.values, p90_high.values, alpha=0.2)

    ax.set_title(f"{label}: Low vs High Exposure (Monte Carlo mean & 10–90% band)")
    ax.set_xlabel("t (last 120 periods kept)")
    ax.set_ylabel(label)
    ax.legend()
    plt.show()
    
    # Pick and compare specific runs (overlay a low-x run vs a high-x run)
    run_idx_low  = 100
    run_idx_high = 100
    one_low  = low_runs[run_idx_low]
    one_high = high_runs[run_idx_high]

    fig, ax = plt.subplots()
    ax.plot(one_low["u_rate"],  label=f"Low x run {run_idx_low}")
    ax.plot(one_high["u_rate"], label=f"High x run {run_idx_high}")
    ax.set_title("Unemployment rate: single-run comparison")
    ax.set_xlabel("t (last 120 periods kept)")
    ax.set_ylabel("u_rate")
    ax.legend()
    plt.show()
    
    




