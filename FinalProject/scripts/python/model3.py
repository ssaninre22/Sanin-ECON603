# ------------------------------------------------------------
# Search–Matching with EV1–Logit Enrollment
# Unemployment as stocks, Monte Carlo, synchronized shocks
# ------------------------------------------------------------
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Tuple, Dict, Callable, Optional
import matplotlib.pyplot as plt
from typing import Sequence


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
    Solve via the firm surplus identity with free entry and Nash wages:
      J = R' - w + β(1-δ)J, and c = β q(θ) J.
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
        for _ in range(10):
            hi *= 2.0
            fhi = residual(hi)
            if np.sign(flo) != np.sign(fhi):
                break
        else:
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
    J = sp.c / (beta * q)
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
    # exogenous paths (COMMON across counties in comparisons)
    eps_path: Optional[np.ndarray] = None,
    rho_path: Optional[np.ndarray] = None,
    chi_path: Optional[np.ndarray] = None,
    s_aid_path: Optional[np.ndarray] = None,
    # aggressiveness: how much target-gap nudges R' (tightness push)
    kappa_gap: float = 0.0,
) -> pd.DataFrame:
    """
    Simulates a single economy path with proper unemployment stocks:
      L_{s,t+1} = (1-δ_s)L_{s,t} + f_s(θ_{s,t}) U_{s,t}
      U_{s,t+1} = (1-f_s(θ_{s,t})) U_{s,t} + δ_s L_{s,t}
    Enrollment uses EV1/logit; θ solved from free entry + Nash each period.
    All exogenous sequences passed here can be SHARED across counties.
    """
    # Defaults for exogenous paths (if not provided)
    if eps_path is None:
        eps_path = rng.normal(0.0, 1.0, size=T)
    if rho_path is None:
        rho_path = np.ones(T)           # full presence
    if chi_path is None:
        chi_path = np.zeros(T)          # no closures
    if s_aid_path is None:
        s_aid_path = 0.0 * eps_path     # flat aid; or 0.1*(-eps_path) if desired

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

        # Tradables sales (demand)
        Y_T = mp.YT_bar + mp.beta_T * eps_t

        # Enrollment (logit) using current θ and wages as work-value proxy
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

        # N fixed point target employment (used only to tilt R')
        L_T_target = mp.phi_T * Y_T
        L_N_target = kx * (mp.eta * mp.wbar * L_T + G_t)
        Y_N = L_N_target / mp.aN

        # Effective revenue-per-worker with smoothing and optional gap tilt
        LT_eff = 0.5 * L_T + 0.5 * L_T_target
        LN_eff = 0.5 * L_N + 0.5 * L_N_target
        Rprime_T = Y_T / max(LT_eff, 1e-6)
        Rprime_N = Y_N / max(LN_eff, 1e-6)

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
# Monte Carlo with SHARED shocks across counties
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
    # common exogenous paths (rules) — SAME for all counties in a run
    build_eps: Callable[[np.random.Generator, int], np.ndarray] = lambda rng, T: rng.normal(0.0, 1.0, size=T),
    build_rho: Callable[[int], np.ndarray] = lambda T: np.ones(T),
    build_chi: Callable[[int], np.ndarray] = lambda T: np.zeros(T),
    build_s: Callable[[np.ndarray], np.ndarray] = lambda eps: 0.1 * (-eps),  # mildly countercyclical aid
    # gap tilt
    kappa_gap: float = 0.0,
    seed: int = 1234
) -> Dict[float, List[pd.DataFrame]]:
    """
    For each exposure in exposure_list and each Monte Carlo run:
      1) Draw ONE set of exogenous sequences (eps, rho, chi, s_aid).
      2) Simulate EACH county with these SAME sequences.
      3) Keep the last T_keep periods for each county.
    Returns dict: exposure -> list of DataFrames (length num_runs).
    """
    assert T_keep <= T_total, "T_keep must be <= T_total."
    rng_master = np.random.default_rng(seed)
    out: Dict[float, List[pd.DataFrame]] = {x: [] for x in exposure_list}

    for _ in range(num_runs):
        rng = np.random.default_rng(rng_master.integers(1, 2**31-1))
        # --- SHARED exogenous sequences for this run ---
        eps_path  = build_eps(rng, T_total)
        rho_path  = build_rho(T_total)                # common presence
        chi_path  = build_chi(T_total)                # common closures
        s_aid     = build_s(eps_path)                 # same mapping from eps

        for x in exposure_list:
            mp = mp_builder(x, T_total)
            df = simulate_one_path(
                T=T_total, rng=rng, mp=mp,
                T_params=T_params, N_params=N_params,
                N_T=N_T, N_N=N_N,
                u0_T=0.05, u0_N=0.05,
                eps_path=eps_path, rho_path=rho_path,
                chi_path=chi_path, s_aid_path=s_aid,
                kappa_gap=kappa_gap
            )
            out[x].append(df.iloc[-T_keep:].reset_index(drop=True))
    return out

# =========================
# Comparison helper (two counties, SHARED shocks)
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
    build_eps: Callable[[np.random.Generator, int], np.ndarray] = lambda rng, T: rng.normal(0.0, 1.0, size=T),
    build_rho: Callable[[int], np.ndarray] = lambda T: np.ones(T),
    build_chi: Callable[[int], np.ndarray] = lambda T: np.zeros(T),
    build_s: Callable[[np.ndarray], np.ndarray] = lambda eps: 0.1 * (-eps),
    kappa_gap: float = 0.0,
    seed: int = 2024
) -> Tuple[List[pd.DataFrame], List[pd.DataFrame]]:
    """
    Runs Monte Carlo for two exposures using the SAME shocks in each run.
    Returns lists of final slices (DataFrames) for the low and high exposure counties.
    """
    res = run_monte_carlo(
        num_runs=num_runs, T_total=T_total, T_keep=T_keep,
        mp_builder=mp_builder, exposure_list=[x_low, x_high],
        T_params=T_params, N_params=N_params,
        N_T=N_T, N_N=N_N,
        build_eps=build_eps, build_rho=build_rho, build_chi=build_chi, build_s=build_s,
        kappa_gap=kappa_gap, seed=seed
    )
    return res[x_low], res[x_high]

# =========================
# Single-path helper (simulate 1000, keep last 120)
# =========================
def simulate_single_county_last120(
    x: float,
    mp_builder: Callable[[float, int], ModelParams],
    T_params: MatchingParams,
    N_params: MatchingParams,
    N_T: float = 0.8,
    N_N: float = 0.8,
    T_total: int = 1000,
    T_keep: int = 120,
    seed: int = 0,
    # common path builders (so it's aligned with MC style)
    build_eps: Callable[[np.random.Generator, int], np.ndarray] = lambda rng, T: rng.normal(0.0, 1.0, size=T),
    build_rho: Callable[[int], np.ndarray] = lambda T: np.ones(T),
    build_chi: Callable[[int], np.ndarray] = lambda T: np.zeros(T),
    build_s: Callable[[np.ndarray], np.ndarray] = lambda eps: 0.1 * (-eps),
    kappa_gap: float = 0.0
) -> pd.DataFrame:
    """
    Simulate a single county for 1000 periods and return only the last 120 periods.
    """
    assert T_keep <= T_total, "T_keep must be <= T_total."
    rng = np.random.default_rng(seed)
    eps_path = build_eps(rng, T_total)
    rho_path = build_rho(T_total)
    chi_path = build_chi(T_total)
    s_aid    = build_s(eps_path)

    mp = mp_builder(x, T_total)
    df = simulate_one_path(
        T=T_total, rng=rng, mp=mp,
        T_params=T_params, N_params=N_params,
        N_T=N_T, N_N=N_N,
        u0_T=0.05, u0_N=0.05,
        eps_path=eps_path, rho_path=rho_path, chi_path=chi_path, s_aid_path=s_aid,
        kappa_gap=kappa_gap
    )
    return df.iloc[-T_keep:].reset_index(drop=True)

# =========================
# Example builders & defaults
# =========================
def make_model_params(x: float, T: int) -> ModelParams:
    """Exposure enters η(x) and aN(x)."""
    eta_x = 0.25 + 0.0 * x
    aN_x  = 0.35 + 0.0 * x  # higher exposure ⇒ more services intensity
    return ModelParams(
        x=x, g0=0.2, g1=0.8, psi=1.4, wbar=1.0,
        eta=eta_x, phi_T=0.9, aN=aN_x,
        beta=0.96,
        mu=1.0, phi_bar=0.0, sigma_phi=0.6, VG_next=1.2,
        b=0.4, YT_bar=1.0, beta_T=0.8
    )

def build_eps_normal(rng: np.random.Generator, T: int) -> np.ndarray:
    """N(0,1) shock for tradables demand."""
    return rng.normal(0.0, 1.0, size=T)

def build_rho_full(T: int) -> np.ndarray:
    """Full local presence throughout."""
    return np.ones(T)

def build_chi_none(T: int) -> np.ndarray:
    """No closure shock."""
    return np.zeros(T)

def build_s_countercyclical(eps: np.ndarray) -> np.ndarray:
    """Aid mildly countercyclical."""
    return 0.1 * (-eps)

# IRF helpers 

def make_irf_paths(
    T_total: int,
    t_shock: int,
    # demand shock in ε (negative for recession):
    eps_size: float = -2.0,
    eps_dur: int = 4,
    # presence shock (Covid-like): set rho to a lower level for a spell
    rho_level: float = None,  # e.g., 0.2 to mimic remote instruction
    rho_dur: int = 0,
    # closure shock χ (adds a negative term -χ x^ψ to G_t):
    chi_level: float = None,  # e.g., 0.8
    chi_dur: int = 0,
    base_rho: float = 1.0,
    s_rule: Callable[[np.ndarray], np.ndarray] = lambda eps: 0.1 * (-eps),
) -> Dict[str, np.ndarray]:
    """
    Construct baseline vs. shocked exogenous paths for a deterministic IRF.

    Returns a dict with: eps_base, eps_shock, rho_base, rho_shock, chi_base, chi_shock, s_base, s_shock
    """
    assert 0 <= t_shock < T_total
    eps_base = np.zeros(T_total)
    eps_shock = eps_base.copy()
    if eps_dur > 0:
        eps_shock[t_shock : min(T_total, t_shock + eps_dur)] = eps_size

    rho_base = np.full(T_total, base_rho)
    rho_shock = rho_base.copy()
    if (rho_level is not None) and (rho_dur > 0):
        rho_shock[t_shock : min(T_total, t_shock + rho_dur)] = rho_level

    chi_base = np.zeros(T_total)
    chi_shock = chi_base.copy()
    if (chi_level is not None) and (chi_dur > 0):
        chi_shock[t_shock : min(T_total, t_shock + chi_dur)] = chi_level

    s_base = s_rule(eps_base)
    s_shock = s_rule(eps_shock)

    return {
        "eps_base": eps_base, "eps_shock": eps_shock,
        "rho_base": rho_base, "rho_shock": rho_shock,
        "chi_base": chi_base, "chi_shock": chi_shock,
        "s_base": s_base, "s_shock": s_shock
    }

def compute_irf_for_exposures(
    exposures: Sequence[float],
    mp_builder: Callable[[float, int], ModelParams],
    T_params: MatchingParams,
    N_params: MatchingParams,
    N_T: float = 0.8,
    N_N: float = 0.8,
    # IRF design
    pre_burn: int = 200,   # initial burn-in before shock window
    H: int = 20,           # IRF horizon (quarters)
    t_shock: int = 200,    # time of shock (>= pre_burn)
    # shock magnitudes/durations
    eps_size: float = -2.0,
    eps_dur: int = 4,
    rho_level: float = None,
    rho_dur: int = 0,
    chi_level: float = None,
    chi_dur: int = 0,
    base_rho: float = 1.0,
    kappa_gap: float = 0.0,
    # IRF scaling
    mode: str = "diff",    # "diff" = level difference, "pct" = percent deviation
    anchors: Sequence[str] = ("u_rate","L_T","L_N","S","G","theta_N"),
    seed: int = 0
) -> Dict[float, pd.DataFrame]:
    """
    For each exposure in `exposures`, simulate:
      - baseline (no shock) and
      - shocked (ε shock +/- presence/closure shocks),
    then compute IRFs over horizon H starting at t_shock.

    Returns dict: exposure -> DataFrame with columns ['t','var','irf'] (long format).
    """
    assert t_shock >= pre_burn, "t_shock should be >= pre_burn."
    T_total = t_shock + H + 1  # run long enough to see the full horizon
    rng = np.random.default_rng(seed)

    # Build deterministic exogenous paths (shared across counties)
    paths = make_irf_paths(
        T_total=T_total, t_shock=t_shock,
        eps_size=eps_size, eps_dur=eps_dur,
        rho_level=rho_level, rho_dur=rho_dur,
        chi_level=chi_level, chi_dur=chi_dur,
        base_rho=base_rho,
        s_rule=lambda eps: 0.1 * (-eps)  # can pass in if you prefer
    )

    out: Dict[float, pd.DataFrame] = {}
    for x in exposures:
        mp = mp_builder(x, T_total)

        # Baseline
        df_b = simulate_one_path(
            T=T_total, rng=rng, mp=mp,
            T_params=T_params, N_params=N_params,
            N_T=N_T, N_N=N_N,
            u0_T=0.05, u0_N=0.05,
            eps_path=paths["eps_base"], rho_path=paths["rho_base"],
            chi_path=paths["chi_base"], s_aid_path=paths["s_base"],
            kappa_gap=kappa_gap
        )

        # Shocked
        df_s = simulate_one_path(
            T=T_total, rng=rng, mp=mp,
            T_params=T_params, N_params=N_params,
            N_T=N_T, N_N=N_N,
            u0_T=0.05, u0_N=0.05,
            eps_path=paths["eps_shock"], rho_path=paths["rho_shock"],
            chi_path=paths["chi_shock"], s_aid_path=paths["s_shock"],
            kappa_gap=kappa_gap
        )

        # Compute pre-shock anchors for percent deviations
        pre_slice = slice(t_shock - pre_burn, t_shock)  # long pre window
        base_means = df_b.loc[pre_slice, anchors].mean()

        # IRF window
        idx = range(t_shock, t_shock + H)
        records = []
        for var in anchors:
            base_series = df_b.loc[idx, var].to_numpy()
            shock_series = df_s.loc[idx, var].to_numpy()
            if mode == "pct":
                denom = np.maximum(1e-8, base_means[var])
                irf_vals = 100.0 * (shock_series - base_series) / denom
            else:
                irf_vals = shock_series - base_series
            for h, val in enumerate(irf_vals):
                records.append({"t": h, "var": var, "irf": val})

        out[x] = pd.DataFrame.from_records(records)

    return out

def plot_irf_panel(
    irfs: Dict[float, pd.DataFrame],
    vars_to_plot: Sequence[str] = ("u_rate","L_T","L_N","S","G","theta_N"),
    exposure_labels: Optional[Dict[float,str]] = None,
    title_prefix: str = "IRF",
    mode: str = "diff"
):
    """
    Simple multi-panel plot of IRFs across exposures for selected variables.
    `irfs` is the output of compute_irf_for_exposures: exposure -> long DF ['t','var','irf'].
    """
    import matplotlib.pyplot as plt
    K = len(vars_to_plot)
    ncols = 3
    nrows = int(np.ceil(K / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3*nrows), squeeze=False)

    for j, var in enumerate(vars_to_plot):
        ax = axes[j // ncols, j % ncols]
        for x, df in irfs.items():
            label = exposure_labels.get(x, f"x={x:.2f}") if exposure_labels else f"x={x:.2f}"
            sub = df[df["var"] == var].sort_values("t")
            ax.plot(sub["t"], sub["irf"], label=label)
        ax.set_title(f"{title_prefix}: {var} ({'Δ' if mode=='diff' else '%Δ'})")
        ax.set_xlabel("quarters after shock")
        ax.axhline(0.0, color="k", linewidth=0.8, alpha=0.5)
        ax.grid(True, alpha=0.2)
        if j == 0:
            ax.legend()

    # Hide any empty subplots
    for k in range(K, nrows*ncols):
        axes[k // ncols, k % ncols].axis("off")

    plt.tight_layout()
    plt.show()
    
    
# =========================
# RUN (examples)
# =========================
if __name__ == "__main__":
    # Matching parameters
    T_params = MatchingParams(m=0.5, alpha=0.5, delta=0.03, c=0.2, xi=0.5)
    N_params = MatchingParams(m=0.6, alpha=0.5, delta=0.03, c=0.2, xi=0.6)

    # --- Single county: simulate 1000, keep last 120 ---
    df_one = simulate_single_county_last120(
        x=0.35,
        mp_builder=make_model_params,
        T_params=T_params, N_params=N_params,
        N_T=0.8, N_N=0.8,
        T_total=1000, T_keep=120, seed=0,
        build_eps=build_eps_normal,
        build_rho=build_rho_full,
        build_chi=build_chi_none,
        build_s=build_s_countercyclical,
        kappa_gap=0.0
    )

    # Example single-run plots (optional)
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

    plot_single_run(df_one, title_suffix=" (x=0.35)")

    # --- Two-county comparison with SHARED shocks per run ---
    NUM_RUNS = 1000
    T_TOTAL  = 400    # simulate long horizon per run
    T_KEEP   = 120    # keep last 120 periods

    low_runs, high_runs = compare_two_counties(
        num_runs=NUM_RUNS, T_total=T_TOTAL, T_keep=T_KEEP,
        x_low=0.025, x_high=0.30,
        mp_builder=make_model_params,
        T_params=T_params, N_params=N_params,
        N_T=0.8, N_N=0.8,
        build_eps=build_eps_normal,
        build_rho=build_rho_full,
        build_chi=build_chi_none,
        build_s=build_s_countercyclical,
        kappa_gap=0.0,
        seed=42
    )

    # Stack Monte Carlo outputs (last 120 periods each)
    df_low  = pd.concat(low_runs,  keys=range(NUM_RUNS), names=["run","t"]).reset_index(level=0)
    df_high = pd.concat(high_runs, keys=range(NUM_RUNS), names=["run","t"]).reset_index(level=0)

    # Example comparison plot (mean & 10–90% bands)
    def summarize_mc(df, col):
        grouped = df.groupby(df.index)[col]
        return grouped.mean(), grouped.quantile(0.10), grouped.quantile(0.90)

    mean_low, p10_low, p90_low   = summarize_mc(df_low,  "u_rate")
    mean_high, p10_high, p90_high = summarize_mc(df_high, "u_rate")

    fig, ax = plt.subplots()
    ax.plot(mean_low.values,  label="Low exposure (mean)")
    ax.fill_between(range(len(p10_low)),  p10_low.values,  p90_low.values,  alpha=0.2)
    ax.plot(mean_high.values, label="High exposure (mean)")
    ax.fill_between(range(len(p10_high)), p10_high.values, p90_high.values, alpha=0.2)
    ax.set_title("Unemployment rate: Low vs High Exposure (MC mean & 10–90% band)")
    ax.set_xlabel("t (last 120 periods kept)")
    ax.set_ylabel("u_rate")
    ax.legend()
    plt.show()

# Choose exposures and model pieces already defined in your file:
exposures = [0.025, 0.25]  # low vs high university exposure

irfs = compute_irf_for_exposures(
    exposures=exposures,
    mp_builder=make_model_params,
    T_params=T_params, N_params=N_params,
    N_T=0.8, N_N=0.8,
    pre_burn=200, H=20, t_shock=200,
    eps_size=-2.0, eps_dur=4,         # big negative ε for one year
    rho_level=None, rho_dur=0,        # no presence shock
    chi_level=None, chi_dur=0,        # no closure shock
    base_rho=1.0,
    kappa_gap=0.0,
    mode="pct",  # or "pct"
    anchors=("u_rate","L_T","L_N","S","G","theta_N"),
    seed=0
)

plot_irf_panel(
    irfs,
    vars_to_plot=("u_rate","L_T","L_N","S","G","theta_N"),
    exposure_labels={0.05:"Low exposure", 0.25:"High exposure"},
    title_prefix="IRF to ε shock",
    mode="diff"
)

# Covid type 

irfs_covid = compute_irf_for_exposures(
    exposures=[0.05, 0.60],
    mp_builder=make_model_params,
    T_params=T_params, N_params=N_params,
    N_T=0.8, N_N=0.8,
    pre_burn=200, H=20, t_shock=200,
    eps_size=-2.0, eps_dur=4,      # demand recession
    rho_level=0.2, rho_dur=8,      # in-person presence collapse for 2 years
    chi_level=0.8, chi_dur=8,      # direct closure shock term
    base_rho=1.0,
    kappa_gap=0.0,
    mode="pct",  # show % deviations vs. pre-shock baseline
    anchors=("u_rate","L_T","L_N","S","G","theta_N"),
    seed=1
)

plot_irf_panel(
    irfs_covid,
    vars_to_plot=("u_rate","L_T","L_N","S","G","theta_N"),
    exposure_labels={0.05:"Low exposure", 0.60:"High exposure"},
    title_prefix="Covid-like IRF (ε + ρ drop + χ)",
    mode="pct"
)
