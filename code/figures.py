"""Supplementary figures for the treatise (static PNG). Run after analysis.py."""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import dates as mdates

C1, C2, C3, INK, INK2, GRID, SURF = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": GRID, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                     "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
                     "axes.titlesize": 10, "axes.titleweight": "bold", "legend.frameon": False, "savefig.dpi": 180})
R = json.load(open("results.json"))
comp = pd.read_csv("digitised_components.csv", index_col=0)
kp = pd.read_csv("kp_ap_daily.csv", index_col=0, parse_dates=True)
tr7 = pd.read_csv("digitised_trends7d.csv", index_col=0, parse_dates=True)
DOY_H = 303
month_ticks = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
month_lab = ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"]

# ---- Figure S1: annual components
from analysis import sA_def   # noqa: E402  (re-uses the Prophet fit computed in analysis.py)
fig, ax = plt.subplots(2, 1, figsize=(8, 6.4), sharex=False)
ax[0].plot(comp.index, comp["fig5_UFO"], color=C1, lw=2, label="'UFO' search interest, yearly component (author's figure, digitised)")
ax[0].axvline(DOY_H, color=INK2, lw=1, ls=":"); ax[0].text(DOY_H + 2, ax[0].get_ylim()[1] * 0.85, "31 Oct", color=INK2, fontsize=8)
ax[0].set_ylabel("% of trend"); ax[0].set_title("Yearly components, estimation window 2004-01-01 to 2021-09-30"); ax[0].set_xticks(month_ticks); ax[0].set_xticklabels(month_lab); ax[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
ax[1].plot(comp.index, comp["fig5_Ap"], color=C2, lw=2, label="Ap, yearly component (author's figure, digitised)")
ax[1].plot(np.arange(365), sA_def, color=C3, lw=2, label="Ap, yearly component recomputed on the definitive GFZ series (Prophet, N = 10)")
ax[1].axvline(DOY_H, color=INK2, lw=1, ls=":")
ax[1].set_ylabel("% of trend"); ax[1].set_xticks(month_ticks); ax[1].set_xticklabels(month_lab); ax[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
fig.tight_layout(); fig.savefig("fig_S1_components.png"); plt.close(fig)

# ---- Figure S2: dot-product coherence
fig, ax = plt.subplots(figsize=(8, 3.6))
for key, col, lab in (("fig4", C1, "UFO (Sept-2021 figure) x Ap (same figure)"), ("fig5", C2, "UFO (pairing figure) x Ap (same figure)")):
    tau = np.array(R[key]["C_curve"]["tau"]); C = np.array(R[key]["C_curve"]["C"])
    ax.plot(tau, C, color=col, lw=2, label=lab)
    i = int(np.argmin(C)); ax.plot(tau[i], C[i], "o", color=col, ms=6); ax.annotate(f"tau = {tau[i]} d, C = {C[i]:.2f}", (tau[i], C[i]), (tau[i] + 8, C[i] - 0.06), fontsize=8, color=INK2)
ax.axhline(0, color=INK2, lw=0.8); ax.axvline(0, color=INK2, lw=0.8, ls=":")
ax.set_xlabel("lag tau (days); tau > 0: search curve lags the geomagnetic curve"); ax.set_ylabel("C(tau)")
ax.set_title("Dot-product coherence of the yearly components on the circle"); ax.legend(loc="upper right")
fig.tight_layout(); fig.savefig("fig_S2_coherence.png"); plt.close(fig)

# ---- Figure S3: event window (UT)
ev = kp.loc["2021-10-25":"2021-11-01"]
t = []; ap = []; kpv = []
for d, row in ev.iterrows():
    for i in range(8):
        t.append(d + pd.Timedelta(hours=3 * i)); ap.append(row[f"ap{i}"]); kpv.append(row[f"Kp{i}"])
t = pd.DatetimeIndex(t); ap = np.array(ap); kpv = np.array(kpv)
u = tr7["ufo_hourly"].copy(); u.index = u.index - pd.Timedelta(hours=3); occ = tr7["occluded"].values
fig, ax = plt.subplots(3, 1, figsize=(8, 7.2), gridspec_kw=dict(height_ratios=[1, 1, 0.8]))
ax[0].bar(t, ap, width=pd.Timedelta(hours=2.8), align="edge", color=C2, lw=0)
def kp_label(v):  # Kp x10 -> Bartels notation: 37 -> 4-, 40 -> 4o, 43 -> 4+
    return f"{v // 10 + 1}-" if v % 10 == 7 else (f"{v // 10}o" if v % 10 == 0 else f"{v // 10}+")
for tt, a, k in zip(t, ap, kpv):
    if a >= 15 and tt < pd.Timestamp("2021-11-01"): ax[0].text(tt + pd.Timedelta(hours=1.4), a + 0.6, "Kp " + kp_label(k), ha="center", fontsize=7, color=INK2)
ax[0].set_ylabel("ap (nT), 3-hourly"); ax[0].set_title("25 Oct - 1 Nov 2021 (UT): definitive 3-hourly ap (GFZ) and hourly 'UFO' search interest")
ax[1].plot(u.index, u.values, color=C1, lw=2, label="'UFO' hourly index (screenshot-recovered; 0-100, window-normalised)")
ax[1].plot(u.index[occ], u.values[occ], "o", ms=3, color=SURF, mec=C1, label="hours occluded by the tooltip in the source (interpolated)")
tp = u.idxmax(); ax[1].plot(tp, u.max(), "o", color=C1, ms=7); ax[1].annotate("31 Oct 03:00 UT = 06:00 UTC+3, index 100", (tp, u.max()), (tp - pd.Timedelta(hours=70), 104), fontsize=8, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
ax[1].set_ylim(0, 112); ax[1].set_ylabel("search index"); ax[1].legend(loc="lower left", fontsize=7)
for a in ax[:2]:
    a.set_xlim(pd.Timestamp("2021-10-25"), pd.Timestamp("2021-11-01 03:00")); a.xaxis.set_major_locator(mdates.DayLocator()); a.xaxis.set_major_formatter(mdates.DateFormatter("%d %b"))
cc = {int(k): v for k, v in R["event_ccf_3h"].items()}; lags = sorted(cc); vals = [cc[l] for l in lags]
ax[2].bar([3 * l for l in lags], vals, width=2.4, color=[C1 if v > 0 else INK2 for v in vals], lw=0)
nl = R["event_ccf_95pct_noise_level"]; ax[2].axhspan(-nl, nl, color=GRID, alpha=0.6, lw=0)
ax[2].axhline(0, color=INK2, lw=0.8); ax[2].set_xlabel("lag (hours); > 0: search index lags ap"); ax[2].set_ylabel("cross-correlation")
ax[2].set_title(f"Cross-correlation, 3-hourly bins (n = {R['event_overlap_bins']}); shaded: +-1.96/sqrt(n)")
fig.tight_layout(); fig.savefig("fig_S3_event.png"); plt.close(fig)

# ---- Figure S4: climatology
from analysis import climatology  # noqa: E402
def circ_smooth(x, w=15):
    k = np.ones(w) / w; return np.convolve(np.concatenate([x[-w:], x, x[:w]]), k, mode="same")[w:-w]
cv = climatology(kp.loc["2004-01-01":"2021-09-30"]); cf = climatology(kp); c03 = climatology(kp[kp.index.year != 2003])
fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(np.arange(365), circ_smooth(cf), color=C2, lw=2, label="1957-2025 (15-day circular mean)")
ax.plot(np.arange(365), circ_smooth(c03), color=C3, lw=2, label="1957-2025 without 2003")
ax.plot(np.arange(365), circ_smooth(cv), color=C1, lw=2, label="2004 - Sept 2021 (the author's window)")
ax.axvline(DOY_H, color=INK2, lw=1, ls=":"); ax.text(DOY_H + 2, ax.get_ylim()[1] * 0.92, "31 Oct", color=INK2, fontsize=8)
ax.set_xticks(month_ticks); ax.set_xticklabels(month_lab); ax.set_ylabel("mean Ap (nT)"); ax.set_title("Day-of-year climatology of Ap"); ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3)
fig.tight_layout(); fig.savefig("fig_S4_climatology.png"); plt.close(fig)
print("figures written")
