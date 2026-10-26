#!/usr/bin/env python3
# Listing 1. Reproducible core of Section 3 (condensed from analysis.py in the archive).
# Inputs: kp_ap_daily.csv (parse_sw.py); digitised_components.csv and
#         digitised_trends7d.csv (digitise.py).
import numpy as np, pandas as pd
from scipy import signal
from prophet import Prophet
from statsmodels.tsa.api import VAR
N, P, H = 10, 365.0, 303          # Fourier order, grid period, 31 Oct (0-based day)

def fourier_fit(y, n=N):
    t = np.arange(len(y)); k = np.arange(1, n + 1)
    arg = 2 * np.pi * np.outer(t, k) / P
    X = np.hstack([np.ones((len(y), 1)), np.cos(arg), np.sin(arg)])
    b = np.linalg.lstsq(X, y, rcond=None)[0]; fit = X @ b
    r2 = 1 - ((y - fit) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return np.hypot(b[1:n + 1], b[n + 1:]), np.arctan2(b[n + 1:], b[1:n + 1]), r2

def coherence(sU, sA, taus):      # C(tau) = <sU(.+tau), sA(.)> / (|sU| |sA|), circular
    u, a = sU - sU.mean(), sA - sA.mean()
    d = np.linalg.norm(u) * np.linalg.norm(a)
    return np.array([np.dot(np.roll(u, -t), a) / d for t in taus])

def phase(s):                     # analytic-signal phase; exact on a full period
    return np.angle(signal.hilbert(s - s.mean()))

def rayleigh(th):
    n = len(th); Rb = abs(np.exp(1j * np.asarray(th)).mean()); Z = n * Rb ** 2
    p = np.exp(-Z) * (1 + (2 * Z - Z ** 2) / (4 * n)
                      - (24 * Z - 132 * Z ** 2 + 76 * Z ** 3 - 9 * Z ** 4) / (288 * n ** 2))
    return Rb, Z, p

def peaks(s, frac=0.10):
    e = np.concatenate([s, s, s])
    p, _ = signal.find_peaks(e, prominence=frac * np.ptp(s))
    return np.unique(p[(p >= len(s)) & (p < 2 * len(s))] - len(s))

def prophet_yearly(series):       # Prophet, multiplicative, yearly order 10, % of trend
    df = pd.DataFrame({"ds": series.index, "y": series.values})
    m = Prophet(yearly_seasonality=N, weekly_seasonality=False,
                daily_seasonality=False, seasonality_mode="multiplicative").fit(df)
    grid = pd.DataFrame({"ds": pd.date_range("2017-01-01", periods=365)})
    return 100 * m.predict(grid)["yearly"].values

kp = pd.read_csv("kp_ap_daily.csv", index_col=0, parse_dates=True)
comp = pd.read_csv("digitised_components.csv", index_col=0)
tr7 = pd.read_csv("digitised_trends7d.csv", index_col=0, parse_dates=True)

# --- 3.1  recompute the Ap component on the definitive series; validate the digitisation
sA_def = prophet_yearly(kp.loc["2004-01-01":"2021-09-30", "Ap"].astype(float))
print("r(definitive, fig4 / fig5 Ap) =",
      np.corrcoef(sA_def, comp["fig4_Ap"])[0, 1].round(3),
      np.corrcoef(sA_def, comp["fig5_Ap"])[0, 1].round(3))
print("definitive Ap component, 31 Oct: %.1f%%, rank %d/365; maxima at DOY %s"
      % (sA_def[H], (sA_def > sA_def[H]).sum() + 1, peaks(sA_def)))

# --- 3.3 / 3.4  coherence, harmonic phases, instantaneous phase, peak pairing
taus = np.arange(-90, 91)
for fig in ("fig4", "fig5"):
    sU, sA = comp[f"{fig}_UFO"].values, comp[f"{fig}_Ap"].values
    C = coherence(sU, sA, taus); i = C.argmin()
    AU, phU, r2U = fourier_fit(sU); AA, phA, r2A = fourier_fit(sA)
    dphi = np.degrees(np.angle(np.exp(1j * (phU - phA))))
    pA, pU = peaks(sA), peaks(sU); lags = []
    for p in pA:                  # first UFO maximum within 45 d after each Ap maximum
        l = (pU - p) % 365; l = l[(l > 0) & (l <= 45)]
        if len(l): lags.append(l.min())
    Rb, Z, pr = rayleigh(2 * np.pi * np.array(lags) / (P / N))
    d0 = phase(sU) - phase(sA)
    print(f"{fig}: R2 {r2U:.4f}/{r2A:.4f} | C(0) {C[taus == 0][0]:+.2f} | "
          f"min C {C[i]:+.2f} at tau {taus[i]:+d} d | dphi(n=2) {dphi[1]:+.0f} deg | "
          f"PLV(0) {abs(np.exp(1j * d0).mean()):.2f}, "
          f"mean dphi {np.degrees(np.angle(np.exp(1j * d0).mean())):+.0f} deg | "
          f"pairs {len(lags)}/{len(pA)}, lag {np.mean(lags):.1f} +- "
          f"{np.std(lags, ddof=1):.1f} d, Rayleigh p {pr:.2f} | "
          f"Ap max nearest 31 Oct: DOY {pA[np.argmin(np.abs((pA - H + 182) % 365 - 182))]}")

# --- 3.2  climatology of Ap on the festival
def clim(df):
    d = df.index.dayofyear.values - 1
    d = np.where(np.asarray(df.index.is_leap_year) & (d >= 59), d - 1, d)
    return df["Ap"].groupby(d).mean().reindex(range(365)).values
windows = (("2004-2021", kp.loc["2004-01-01":"2021-09-30"]), ("1957-2025", kp),
           ("excl. 2003", kp[kp.index.year != 2003]))
for lab, df in windows:
    c = clim(df); rng = np.random.default_rng(0)
    w = lambda s: np.mean([c[(s + k) % 365] for k in range(-2, 3)])
    pW = np.mean([w(s) >= w(H) for s in rng.integers(0, 365, 20000)])
    print(f"climatology {lab}: Ap(31 Oct) {c[H]:.1f} vs mean {df['Ap'].mean():.1f}, "
          f"percentile {100 * (c < c[H]).mean():.0f}, window 29 Oct-2 Nov perm. p {pW:.3f}")

# --- 3.5  event window: 3-hourly ap (UT) against the hourly search index (UTC+3 -> UT)
ev = kp.loc["2021-10-25":"2021-11-01"]; rows = []
for d, r in ev.iterrows():
    rows += [dict(t=d + pd.Timedelta(hours=3 * i), ap=r[f"ap{i}"], Kp10=r[f"Kp{i}"])
             for i in range(8)]
tab = pd.DataFrame(rows).set_index("t")
u = tr7["ufo_hourly"].copy(); u.index = u.index - pd.Timedelta(hours=3)
both = pd.concat([tab["ap"], u.resample("3h").mean().rename("ufo")], axis=1,
                 sort=True).dropna()
x = np.log1p(both["ap"].values); y = both["ufo"].values
x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std()
cc = {l: (np.mean(x[:len(x) - l] * y[l:]) if l >= 0 else np.mean(x[-l:] * y[:len(y) + l]))
      for l in range(-8, 9)}
print("search max (UT):", u.idxmax(), "| CCF max at lag", 3 * max(cc, key=cc.get),
      "h, r =", round(max(cc.values()), 3), "| noise", round(1.96 / np.sqrt(len(x)), 3))
data = pd.DataFrame({"ap": np.log1p(both["ap"].values), "ufo": both["ufo"].values},
                    index=both.index)
p = max(int(VAR(data).select_order(4).bic), 1); res = VAR(data).fit(p)
print(f"VAR({p}): Granger ap->ufo p = "
      f"{res.test_causality('ufo', ['ap'], kind='f').pvalue:.3f}; ufo->ap p = "
      f"{res.test_causality('ap', ['ufo'], kind='f').pvalue:.3f}; "
      f"IRF peak at {3 * int(np.argmax(res.irf(8).irfs[:, 1, 0]))} h")

# --- Section 2 calibration: Kuramoto ensemble whose coupling is driven by the field
rng = np.random.default_rng(1); n, sig, K0, dt = 400, 0.5, 0.5, 0.05
om, th = rng.normal(0, sig, n), rng.uniform(-np.pi, np.pi, n)
apv = np.repeat(tab["ap"].values, int(3 / dt)); r = []
for a in apv:                     # K(t) = K0 + ap(t)/27 ; mean-field Kuramoto step
    z = np.exp(1j * th).mean(); r.append(abs(z))
    th = th + dt * (om + (K0 + a / 27) * abs(z) * np.sin(np.angle(z) - th))
r = np.array(r); Kc = sig * np.sqrt(8 / np.pi)
print(f"Kuramoto: Kc = {Kc:.2f} (crossed at ap = {(Kc - K0) * 27:.1f} nT); "
      f"r 25-29 Oct {r[:int(120 / dt)].mean():.2f}, r 31 Oct "
      f"{r[int(144 / dt):int(168 / dt)].mean():.2f}, first r > 0.5 at "
      f"{tab.index[0] + pd.Timedelta(hours=float(np.argmax(r > 0.5) * dt))} UT")
