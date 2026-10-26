#!/usr/bin/env python3
"""Quantitative appendix: geomagnetic variance (Kp/Ap) x search interest ('UFO'), annual phase alignment
and the 31 October 2021 event window.

Inputs
  kp_ap_daily.csv            definitive 3-hourly Kp/ap and daily Ap, 1957-10-01 .. 2025-07-20 (GFZ via CelesTrak CSSI, parse_sw.py)
  digitised_components.csv   Prophet yearly components recovered from the author's Sept-2021 figures (digitise.py)
  digitised_trends7d.csv     Google Trends 'UFO' worldwide, hourly, 24 Oct 15:00 .. 31 Oct 15:00 (UTC+3), recovered from the screenshot
Optional
  multiTimeline.csv          Google Trends export of 'UFO', worldwide, 2004-01-01 .. 2021-09-30 (replaces the digitised UFO component)
Outputs
  results.json, fig_*.png
"""
import json, warnings
import numpy as np, pandas as pd
from scipy import signal
warnings.filterwarnings("ignore")

VINTAGE = ("2004-01-01", "2021-09-30")      # the author's estimation window (Halloween 2021 is out of sample)
N_FOURIER, P_GRID = 10, 365.0               # Prophet yearly order; period on the 365-day plotting grid
DOY_HALLOWEEN = 303                         # 31 October, 0-based day of a non-leap year
R = {}                                      # results

# ------------------------------------------------------------------ Fourier machinery
def fourier_design(t, n=N_FOURIER, period=P_GRID):
    k = np.arange(1, n + 1)
    arg = 2 * np.pi * np.outer(t, k) / period
    return np.hstack([np.cos(arg), np.sin(arg)])

def fit_fourier(y, n=N_FOURIER):
    """OLS projection of a 365-day curve on the order-n Fourier basis; returns amplitudes, phases, R^2."""
    t = np.arange(len(y)); X = np.hstack([np.ones((len(y), 1)), fourier_design(t, n)])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    a, b = beta[1:n + 1], beta[n + 1:]
    amp, ph = np.hypot(a, b), np.arctan2(b, a)
    fit = X @ beta; r2 = 1 - ((y - fit) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return dict(amp=amp, phase=ph, r2=float(r2), fit=fit)

def dot_product_coherence(sU, sA, taus):
    """C(tau) = <sU(.+tau), sA(.)> / (|sU| |sA|) on the circle; tau>0 means the UFO curve lags the Ap curve by tau days."""
    u = sU - sU.mean(); a = sA - sA.mean(); den = np.linalg.norm(u) * np.linalg.norm(a)
    return np.array([np.dot(np.roll(u, -tau), a) / den for tau in taus])

def inst_phase(s):
    """instantaneous phase of a periodic curve via the FFT-based analytic signal (exact on a full period)."""
    return np.angle(signal.hilbert(s - s.mean()))

def plv(sU, sA, tau=0):
    d = inst_phase(np.roll(sU, -tau)) - inst_phase(sA)
    return float(abs(np.exp(1j * d).mean()))

def peaks_circular(s, prominence_frac=0.10):
    ext = np.concatenate([s, s, s]); pk, _ = signal.find_peaks(ext, prominence=prominence_frac * np.ptp(s))
    return np.sort(np.unique(pk[(pk >= len(s)) & (pk < 2 * len(s))] - len(s)))

def pair_peaks(pA, pU, max_lag=45, period=365):
    """each Ap maximum is paired with the first subsequent UFO maximum within max_lag days (circular)."""
    pairs = []
    for p in pA:
        lags = (pU - p) % period; lags = lags[(lags > 0) & (lags <= max_lag)]
        if len(lags): pairs.append((int(p), int((p + lags.min()) % period), int(lags.min())))
    return pairs

def rayleigh(angles):
    n = len(angles); Rbar = abs(np.exp(1j * np.asarray(angles)).mean()); Z = n * Rbar ** 2
    p = np.exp(-Z) * (1 + (2 * Z - Z ** 2) / (4 * n) - (24 * Z - 132 * Z ** 2 + 76 * Z ** 3 - 9 * Z ** 4) / (288 * n ** 2))
    return float(Rbar), float(Z), float(min(max(p, 0), 1))

# ------------------------------------------------------------------ data
kp = pd.read_csv("kp_ap_daily.csv", index_col=0, parse_dates=True)
comp = pd.read_csv("digitised_components.csv", index_col=0)
tr7 = pd.read_csv("digitised_trends7d.csv", index_col=0, parse_dates=True)
apcols = [f"ap{i}" for i in range(8)]; kpcols = [f"Kp{i}" for i in range(8)]

# ------------------------------------------------------------------ 1. Ap yearly component, recomputed with Prophet on the definitive series
from prophet import Prophet
def prophet_yearly(series, mode="multiplicative"):
    df = pd.DataFrame({"ds": series.index, "y": series.values})
    m = Prophet(yearly_seasonality=N_FOURIER, weekly_seasonality=False, daily_seasonality=False, seasonality_mode=mode)
    m.fit(df)
    grid = pd.DataFrame({"ds": pd.date_range("2017-01-01", periods=365, freq="D")})   # Prophet's own plotting year
    y = m.predict(grid)["yearly"].values
    return 100 * y if mode == "multiplicative" else y                                   # percent of trend

ap_vint = kp.loc[VINTAGE[0]:VINTAGE[1], "Ap"].astype(float)
sA_def = prophet_yearly(ap_vint)
for fig in ("fig4", "fig5"):
    r = np.corrcoef(sA_def, comp[f"{fig}_Ap"])[0, 1]
    R[f"validation_corr_definitive_vs_{fig}_Ap"] = float(r)
    R[f"validation_rmse_definitive_vs_{fig}_Ap_pct"] = float(np.sqrt(np.mean((sA_def - comp[f"{fig}_Ap"]) ** 2)))
R["definitive_Ap_component_range_pct"] = [float(sA_def.min()), float(sA_def.max())]
R["definitive_Ap_component_at_Oct31_pct"] = float(sA_def[DOY_HALLOWEEN])
R["definitive_Ap_component_rank_of_Oct31"] = int((sA_def > sA_def[DOY_HALLOWEEN]).sum() + 1)

# optional: replace the digitised UFO component with a Prophet fit of the real Trends export
try:
    mt = pd.read_csv("multiTimeline.csv", skiprows=2); mt.columns = ["ds", "y"]
    mt["ds"] = pd.to_datetime(mt["ds"]); mt = mt.set_index("ds")["y"].astype(float)
    sU_src = prophet_yearly(mt.loc[VINTAGE[0]:VINTAGE[1]]); R["ufo_component_source"] = "multiTimeline.csv (Prophet)"
except FileNotFoundError:
    sU_src = None; R["ufo_component_source"] = "digitised from the author's figures"

# ------------------------------------------------------------------ 2. phase alignment of the annual components
taus = np.arange(-90, 91)
for fig in ("fig4", "fig5"):
    sU = comp[f"{fig}_UFO"].values if sU_src is None else sU_src
    sA = comp[f"{fig}_Ap"].values
    fU, fA = fit_fourier(sU), fit_fourier(sA)
    C = dot_product_coherence(sU, sA, taus); i = int(np.argmax(C)); tau_star = int(taus[i])
    # broad maximum: the half-height band around tau*
    band = taus[C >= C[i] - 0.5 * (C[i] - C.min())]
    pA, pU = peaks_circular(sA), peaks_circular(sU); pairs = pair_peaks(pA, pU)
    lags = np.array([l for _, _, l in pairs], float)
    Rbar, Z, p = rayleigh(2 * np.pi * lags / (P_GRID / N_FOURIER))        # on the circle of the n=10 harmonic (36.5 d)
    power = fA["amp"] ** 2
    # Halloween: nearest Ap maximum and its paired UFO maximum
    dH = (pA - DOY_HALLOWEEN + 182) % 365 - 182; pH = int(pA[np.argmin(np.abs(dH))])
    pairH = [pr for pr in pairs if pr[0] == pH]
    R[fig] = dict(
        fourier_r2_UFO=fU["r2"], fourier_r2_Ap=fA["r2"],
        amp_UFO=fU["amp"].round(2).tolist(), amp_Ap=fA["amp"].round(2).tolist(),
        phase_diff_deg=np.degrees(np.angle(np.exp(1j * (fU["phase"] - fA["phase"])))).round(1).tolist(),
        Ap_power_share_n_ge_8=float(power[7:].sum() / power.sum()),
        UFO_power_share_n_ge_8=float((fU["amp"] ** 2)[7:].sum() / (fU["amp"] ** 2).sum()),
        C0=float(C[taus == 0][0]), tau_star_days=tau_star, C_tau_star=float(C[i]),
        C_halfheight_band_days=[int(band.min()), int(band.max())],
        C_min=float(C.min()), tau_min=int(taus[np.argmin(C)]),
        PLV0=plv(sU, sA, 0), PLV_tau_star=plv(sU, sA, tau_star),
        n_Ap_peaks=len(pA), n_UFO_peaks=len(pU), n_pairs=len(pairs), pairs_doy=pairs,
        lag_mean=float(lags.mean()) if len(lags) else None, lag_sd=float(lags.std(ddof=1)) if len(lags) > 1 else None,
        rayleigh_Rbar=Rbar, rayleigh_Z=Z, rayleigh_p=p,
        Ap_peak_nearest_Oct31_doy=pH, Ap_peak_nearest_Oct31_offset_days=int(dH[np.argmin(np.abs(dH))]),
        halloween_pair=pairH[0] if pairH else None,
        sA_at_Oct31=float(sA[DOY_HALLOWEEN]), sU_at_Oct31=float(sU[DOY_HALLOWEEN]),
    )
    R[fig]["C_curve"] = dict(tau=taus.tolist(), C=C.round(4).tolist())

# ------------------------------------------------------------------ 3. climatology of Ap on the definitive record
def doy_index(idx):  # 0-based day of year with 29 Feb folded onto 28 Feb, so every year has 365 slots
    d = idx.dayofyear.values - 1; leap = np.asarray(idx.is_leap_year) & (d >= 59); return np.where(leap, d - 1, d)
def climatology(df):
    g = df["Ap"].groupby(doy_index(df.index)).mean(); return g.reindex(range(365)).values
def perm_window_p(clim, lo, hi, n_perm=20000, seed=0):
    rng = np.random.default_rng(seed); w = hi - lo + 1
    obs = np.mean([clim[(lo + k) % 365] for k in range(w)])
    null = np.array([np.mean([clim[(s + k) % 365] for k in range(w)]) for s in rng.integers(0, 365, n_perm)])
    return float(obs), float((null >= obs).mean())
for label, df in (("vintage_2004_2021", kp.loc[VINTAGE[0]:VINTAGE[1]]), ("full_1957_2025", kp)):
    clim = climatology(df); order = np.argsort(-clim)
    obsW, pW = perm_window_p(clim, DOY_HALLOWEEN - 2, DOY_HALLOWEEN + 2)
    # semiannual (Russell-McPherron) and annual harmonics of the climatology
    f2 = fit_fourier(clim, 2); peaks2 = [(k, (np.degrees(f2["phase"][k - 1]) / 360 * 365 / k) % (365 / k)) for k in (1, 2)]
    R[f"clim_{label}"] = dict(
        mean_Ap=float(df["Ap"].mean()), Ap_Oct31_mean=float(clim[DOY_HALLOWEEN]),
        Oct31_rank_of_365=int(np.where(order == DOY_HALLOWEEN)[0][0] + 1),
        Oct31_percentile=float(100 * (clim < clim[DOY_HALLOWEEN]).mean()),
        window_Oct29_Nov2_mean=obsW, window_perm_p=pW,
        harmonic_amp_pct_of_mean={"n1": float(100 * f2["amp"][0] / clim.mean()), "n2": float(100 * f2["amp"][1] / clim.mean())},
        semiannual_peak_doys=[float((365 / 2) * (-f2["phase"][1]) / (2 * np.pi) % (365 / 2) + k * 365 / 2) for k in (0, 1)],
        annual_peak_doy=float(365 * (-f2["phase"][0]) / (2 * np.pi) % 365),
    )
# share of years in which Oct 31 exceeds the year's median Ap
yrs = kp.loc["1958":"2024"]; med = yrs["Ap"].groupby(yrs.index.year).median()
o31 = yrs[(yrs.index.month == 10) & (yrs.index.day == 31)]["Ap"]; o31.index = o31.index.year
R["Oct31_above_annual_median_share"] = float((o31 > med.loc[o31.index]).mean())
R["Oct31_2003_Ap"] = int(kp.loc["2003-10-31", "Ap"]); R["Oct29_2003_Ap"] = int(kp.loc["2003-10-29", "Ap"])
# year-to-year phase stability of the n=1 and n=2 harmonics of the annual Ap cycle (vintage years)
ph = {1: [], 2: []}
for y in range(2004, 2021):
    s = kp.loc[f"{y}-01-01":f"{y}-12-31", "Ap"].values[:365].astype(float); s = s / s.mean() - 1
    f = fit_fourier(s, 2)
    for n in (1, 2): ph[n].append(f["phase"][n - 1])
R["Ap_annual_phase_stability"] = {f"n{n}": dict(zip(("Rbar", "Z", "p"), rayleigh(ph[n]))) for n in (1, 2)}

# ------------------------------------------------------------------ 4. event window 25 Oct - 1 Nov 2021
ev = kp.loc["2021-10-25":"2021-11-01"]
tab = []
for d, row in ev.iterrows():
    for i in range(8):
        tab.append(dict(t_ut=d + pd.Timedelta(hours=3 * i), Kp10=int(row[kpcols[i]]), ap=int(row[apcols[i]])))
tab = pd.DataFrame(tab).set_index("t_ut")
R["event_daily_Ap"] = {str(d.date()): int(v) for d, v in ev["Ap"].items()}
R["event_Oct31_rank_in_week_Oct25_31"] = int((ev.loc[:"2021-10-31", "Ap"] > ev.loc["2021-10-31", "Ap"]).sum() + 1)
R["event_ap_max_bins"] = [str(t) for t in tab.index[tab["ap"] == tab["ap"].max()]]; R["event_ap_max"] = int(tab["ap"].max())
# Trends hourly (UTC+3 -> UT), 3-hourly bin means aligned with the Kp bins
u = tr7["ufo_hourly"].copy(); u.index = u.index - pd.Timedelta(hours=3); u.index.name = "t_ut"
u3 = u.resample("3h", label="left", closed="left").mean()
both = pd.concat([tab["ap"], u3.rename("ufo")], axis=1).dropna()
R["event_overlap_bins"] = int(len(both)); R["event_overlap_span_ut"] = [str(both.index[0]), str(both.index[-1])]
# search maximum and its position relative to the geomagnetic maxima
t_peak_ut = u.idxmax(); R["event_search_max_ut"] = str(t_peak_ut); R["event_search_max_local"] = str(t_peak_ut + pd.Timedelta(hours=3))
R["event_ap_in_bin_of_search_max"] = int(tab.loc[tab.index[tab.index <= t_peak_ut][-1], "ap"])
R["event_Kp10_in_bin_of_search_max"] = int(tab.loc[tab.index[tab.index <= t_peak_ut][-1], "Kp10"])
first_max = tab.index[tab["ap"] == tab["ap"].max()][0]; prev_high = tab.loc[:t_peak_ut, "ap"].idxmax()
R["event_lag_h_from_preceding_ap_max_centre_to_search_max"] = float((t_peak_ut - (prev_high + pd.Timedelta(hours=1.5))) / pd.Timedelta(hours=1))
R["event_preceding_ap_max"] = dict(t_ut=str(prev_high), ap=int(tab.loc[prev_high, "ap"]), Kp10=int(tab.loc[prev_high, "Kp10"]))
R["event_lead_h_from_search_max_to_first_week_ap_max_centre"] = float(((first_max + pd.Timedelta(hours=1.5)) - t_peak_ut) / pd.Timedelta(hours=1))
# cross-correlation, 3-hourly bins (log1p ap), lags in bins; positive lag = ufo lags ap
x = np.log1p(both["ap"].values); y = both["ufo"].values; x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std()
def ccf(x, y, L):
    out = {}
    for l in range(-L, L + 1):
        out[l] = float(np.mean(x[:len(x) - l] * y[l:])) if l >= 0 else float(np.mean(x[-l:] * y[:len(y) + l]))
    return out
cc = ccf(x, y, 8); R["event_ccf_3h"] = {str(k): round(v, 3) for k, v in cc.items()}
kmax = max(cc, key=cc.get); R["event_ccf_max"] = dict(lag_bins=kmax, lag_hours=3 * kmax, r=round(cc[kmax], 3))
R["event_ccf_95pct_noise_level"] = float(1.96 / np.sqrt(len(x)))
# hourly cross-correlation (ap held constant within its bin)
ap_h = tab["ap"].reindex(pd.date_range(tab.index[0], tab.index[-1] + pd.Timedelta(hours=2), freq="h")).ffill()
bh = pd.concat([np.log1p(ap_h).rename("ap"), u.rename("ufo")], axis=1).dropna()
xh = (bh["ap"] - bh["ap"].mean()) / bh["ap"].std(); yh = (bh["ufo"] - bh["ufo"].mean()) / bh["ufo"].std()
cch = ccf(xh.values, yh.values, 24); kh = max(cch, key=cch.get)
R["event_ccf_hourly_max"] = dict(lag_hours=kh, r=round(cch[kh], 3), n=int(len(bh)), noise_95=float(1.96 / np.sqrt(len(bh))))
# bivariate VAR on the 3-hourly window; lag order by BIC; Granger tests both ways
from statsmodels.tsa.api import VAR
data = pd.DataFrame({"ap": np.log1p(both["ap"].values), "ufo": both["ufo"].values}, index=both.index)
sel = VAR(data).select_order(maxlags=4); p_bic = int(sel.bic) if int(sel.bic) > 0 else 1
res = VAR(data).fit(p_bic)
g1 = res.test_causality("ufo", ["ap"], kind="f"); g2 = res.test_causality("ap", ["ufo"], kind="f")
irf = res.irf(8).irfs   # shape (h+1, k, k): response of var j to shock in var i -> irfs[h, j, i]
R["event_VAR"] = dict(n=int(len(data)), p_bic=p_bic, bic_by_lag={str(k): float(v) for k, v in zip(range(5), sel.ics["bic"])},
                      granger_ap_to_ufo=dict(F=float(g1.test_statistic), p=float(g1.pvalue), df=[int(v) for v in np.atleast_1d(g1.df)]),
                      granger_ufo_to_ap=dict(F=float(g2.test_statistic), p=float(g2.pvalue), df=[int(v) for v in np.atleast_1d(g2.df)]),
                      coef_ufo_on_ap_lag1=float(res.params.loc["L1.ap", "ufo"]), coef_ufo_on_ufo_lag1=float(res.params.loc["L1.ufo", "ufo"]),
                      irf_ufo_to_ap_shock=[round(float(irf[h, 1, 0]), 3) for h in range(9)],
                      irf_peak_h=int(np.argmax(irf[:, 1, 0])), resid_corr=float(res.resid_corr[0, 1]))
R["trends7d_occluded_hours"] = int(tr7["occluded"].sum())
R["trends7d_dropout"] = dict(t_local=str(tr7["ufo_hourly"].idxmin()), value=float(tr7["ufo_hourly"].min()))

# ------------------------------------------------------------------ 5. crowd-tuning: Kuramoto ensemble driven by the event-window field
def kuramoto(ap_series, n=400, sigma=0.5, K0=0.5, kappa=1.0, ap_ref=27.0, dt=0.05, seed=1):
    """dtheta_i/dt = omega_i + K(t)/n * sum_j sin(theta_j - theta_i), K(t) = K0 + kappa * ap(t)/ap_ref.
    K0 is below the Gaussian critical coupling Kc = sigma*sqrt(8/pi); the field pushes the ensemble across Kc."""
    rng = np.random.default_rng(seed); om = rng.normal(0, sigma, n); th = rng.uniform(-np.pi, np.pi, n)
    hours = np.arange(0, 3 * len(ap_series), dt); apv = np.repeat(ap_series, int(round(3 / dt)))
    r = np.empty(len(hours))
    for k, h in enumerate(hours):
        K = K0 + kappa * apv[k] / ap_ref; z = np.exp(1j * th).mean(); r[k] = abs(z)
        th = th + dt * (om + K * abs(z) * np.sin(np.angle(z) - th))
    return hours, r, sigma * np.sqrt(8 / np.pi)
hrs, rr, Kc = kuramoto(tab["ap"].values)
t0 = tab.index[0]; imax = int(np.argmax(rr))
R["kuramoto"] = dict(Kc=float(Kc), K0=0.5, ap_threshold_for_Kc=float((Kc - 0.5) * 27.0),
                     r_initial=float(rr[:int(24 / 0.05)].mean()), r_max=float(rr.max()), t_r_max_ut=str(t0 + pd.Timedelta(hours=float(hrs[imax]))),
                     r_mean_Oct25_29=float(rr[:int(5 * 24 / 0.05)].mean()), r_mean_Oct31=float(rr[int(6 * 24 / 0.05):int(7 * 24 / 0.05)].mean()))

R["kuramoto"]["t_r_first_above_0p5_ut"] = str(t0 + pd.Timedelta(hours=float(hrs[int(np.argmax(rr > 0.5))])))

# ------------------------------------------------------------------ 6. supplementary quantities
# definitive Ap component: maxima and the one nearest 31 October
pk_def = peaks_circular(sA_def); dH = (pk_def - DOY_HALLOWEEN + 182) % 365 - 182
R["definitive_Ap_component_peaks_doy"] = [int(p) for p in pk_def]
R["definitive_Ap_component_peak_nearest_Oct31"] = dict(doy=int(pk_def[np.argmin(np.abs(dH))]), offset_days=int(dH[np.argmin(np.abs(dH))]),
                                                       value_pct=float(sA_def[pk_def[np.argmin(np.abs(dH))]]))
# circular mean of the instantaneous phase difference (degrees) at tau = 0 and at the anti-coherent lag
for fig in ("fig4", "fig5"):
    sU = comp[f"{fig}_UFO"].values if sU_src is None else sU_src; sA = comp[f"{fig}_Ap"].values
    for tau in (0, R[fig]["tau_min"]):
        d = inst_phase(np.roll(sU, -tau)) - inst_phase(sA)
        R[fig][f"mean_phase_diff_deg_tau{tau}"] = float(np.degrees(np.angle(np.exp(1j * d).mean())))
    # coherence against the definitive Ap component (scale-free check of the digitised pairing)
    Cd = dot_product_coherence(sU, sA_def, taus)
    R[fig]["C_vs_definitive_Ap"] = dict(C0=float(Cd[taus == 0][0]), tau_min=int(taus[np.argmin(Cd)]), C_min=float(Cd.min()),
                                        tau_max_within_45=int(taus[(taus > 0) & (taus <= 45)][np.argmax(Cd[(taus > 0) & (taus <= 45)])]),
                                        C_max_within_45=float(Cd[(taus > 0) & (taus <= 45)].max()))
# full record without 2003 (the 'Halloween storms')
k03 = kp[kp.index.year != 2003]; clim03 = climatology(k03); obsW, pW = perm_window_p(clim03, DOY_HALLOWEEN - 2, DOY_HALLOWEEN + 2)
R["clim_full_excl_2003"] = dict(Ap_Oct31_mean=float(clim03[DOY_HALLOWEEN]), Oct31_percentile=float(100 * (clim03 < clim03[DOY_HALLOWEEN]).mean()),
                                window_Oct29_Nov2_mean=obsW, window_perm_p=pW, mean_Ap=float(k03["Ap"].mean()))
# where Ap = 14 and Ap = 1 sit in the vintage distribution of daily Ap
v = kp.loc[VINTAGE[0]:VINTAGE[1], "Ap"]
R["vintage_share_days_Ap_ge_14"] = float((v >= 14).mean()); R["vintage_share_days_Ap_le_1"] = float((v <= 1).mean())
hal = kp[(kp.index.month == 10) & (kp.index.day == 31) & (kp.index.year >= 2004) & (kp.index.year <= 2024)]["Ap"]
R["halloweens_2004_2024_Ap"] = {int(y): int(a) for y, a in zip(hal.index.year, hal.values)}
R["halloween_2021_rank_among_2004_2024"] = int((hal > hal.loc["2021-10-31"]).sum() + 1)
# the 06:00 (03:00 UT) maximum against the window's diurnal profile (same hour of day, other days)
uh = tr7["ufo_hourly"].copy(); prof = uh.groupby(uh.index.hour)
h6 = uh.loc["2021-10-31 06:00"]; others = uh[(uh.index.hour == 6) & (uh.index.date != pd.Timestamp("2021-10-31").date())]
R["event_06_local_vs_same_hour_other_days"] = dict(value=float(h6), other_days_mean=float(others.mean()), other_days_sd=float(others.std(ddof=1)),
                                                   z=float((h6 - others.mean()) / others.std(ddof=1)), n_other_days=int(len(others)))
# daily means of the hourly index (local days) for the window
R["trends7d_daily_mean_local"] = {str(d): round(float(m), 1) for d, m in uh.groupby(uh.index.date).mean().items()}

json.dump(R, open("results.json", "w"), indent=1, default=float)
print(json.dumps({k: v for k, v in R.items() if k not in ("fig4", "fig5")}, indent=1, default=float))
for fig in ("fig4", "fig5"):
    print(fig, json.dumps({k: v for k, v in R[fig].items() if k != "C_curve"}, default=float))
