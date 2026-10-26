"""Pixel-to-data recovery of the three source figures.

fig4_components.png : Prophet yearly components, Sept-2021 vintage (UFO top, Ap bottom), composite image
fig5_pairings.png   : Prophet yearly components with the author's lead-lag pairing lines (UFO top, Ap bottom)
fig7_trends7d.png   : Google Trends 'UFO' worldwide, past 7 days, hourly, screenshot taken 31 Oct 2021 (observer UTC+3)

Calibration uses the gridlines detected in the images (tick values read from the figure labels).
Output: digitised_components.csv (day-of-year grid 0..364) and digitised_trends7d.csv (hourly).
"""
import numpy as np, pandas as pd
from PIL import Image

DOY_TICKS = np.array([0, 59, 120, 181, 243, 304, 365])  # Jan1 Mar1 May1 Jul1 Sep1 Nov1 Jan1 (non-leap year, as Prophet plots)

def load(name):
    return np.array(Image.open(name).convert("RGB")).astype(int)

def grid_rows(im, x0, x1, y0, y1, lo=170, hi=245, frac=0.5):
    reg = im[y0:y1, x0:x1]; r, g, b = reg[..., 0], reg[..., 1], reg[..., 2]
    grey = (abs(r - g) < 8) & (abs(g - b) < 8) & (r >= lo) & (r <= hi)
    rows = np.where(grey.mean(axis=1) > frac)[0] + y0
    return cluster(rows)

def grid_cols(im, x0, x1, y0, y1, lo=170, hi=245, frac=0.5):
    reg = im[y0:y1, x0:x1]; r, g, b = reg[..., 0], reg[..., 1], reg[..., 2]
    grey = (abs(r - g) < 8) & (abs(g - b) < 8) & (r >= lo) & (r <= hi)
    cols = np.where(grey.mean(axis=0) > frac)[0] + x0
    return cluster(cols)

def cluster(idx, gap=3):
    """average runs of adjacent indices (a 2-px-wide gridline -> its centre)"""
    if len(idx) == 0: return np.array([])
    out, cur = [], [idx[0]]
    for i in idx[1:]:
        if i - cur[-1] <= gap: cur.append(i)
        else: out.append(np.mean(cur)); cur = [i]
    out.append(np.mean(cur)); return np.array(out)

def trace(im, mask, x0, x1, y0, y1):
    """per-column median row of masked pixels inside the panel; NaN where none"""
    ys = np.full(x1 - x0, np.nan)
    for j, x in enumerate(range(x0, x1)):
        rows = np.where(mask[y0:y1, x])[0]
        if len(rows): ys[j] = np.median(rows) + y0
    return np.arange(x0, x1), ys

def lin_from_ticks(px, vals):
    A = np.vstack([px, np.ones_like(px)]).T
    (a, b), res, *_ = np.linalg.lstsq(A, vals, rcond=None)
    fit = a * px + b
    return a, b, np.abs(fit - vals).max()

def to_doy_grid(xpx, ypx, xa, xb, ya, yb):
    doy = xa * xpx + xb; val = ya * ypx + yb
    ok = ~np.isnan(val)
    grid = np.arange(365)
    # circular interpolation across gaps
    d, v = doy[ok], val[ok]
    order = np.argsort(d); d, v = d[order], v[order]
    dd = np.concatenate([d - 365, d, d + 365]); vv = np.concatenate([v, v, v])
    return grid, np.interp(grid, dd, vv), ok.mean()

# ---------------------------------------------------------------- Fig 5 (clean Prophet export)
def digitise_fig5():
    im = load("fig5_pairings.png"); r, g, b = im[..., 0], im[..., 1], im[..., 2]
    blue = (b > 140) & (r < 60) & (g > 80) & (g < 150)
    # blue blended with the yellow highlighter (alpha ~0.5): (128,185,89)
    blend = (r > 90) & (r < 175) & (g > 150) & (g < 220) & (b > 50) & (b < 135)
    panels = {"UFO": (50, 222, [200, 100, 0, -100, -200]), "Ap": (284, 454, [150, 100, 50, 0, -50, -100])}
    xg = grid_cols(im, 80, 647, 60, 215); assert len(xg) == 7, xg
    xa, xb, xres = lin_from_ticks(xg, DOY_TICKS)
    out = {}
    for name, (top, bot, ticks) in panels.items():
        yg = grid_rows(im, 82, 645, top + 2, bot - 2); assert len(yg) == len(ticks), (name, yg)
        ya, yb, yres = lin_from_ticks(yg, np.array(ticks, float))
        xpx, ypx = trace(im, blue | blend, 81, 647, top + 1, bot)
        grid, vals, cov = to_doy_grid(xpx, ypx, xa, xb, ya, yb)
        out[name] = vals
        print(f"fig5 {name}: x-res {xres:.2f} d, y-res {yres:.2f} %, coverage {cov:.3f}, range [{vals.min():.1f}, {vals.max():.1f}]")
    return pd.DataFrame(out, index=pd.Index(np.arange(365), name="doy"))

# ---------------------------------------------------------------- Fig 4 (composite; faint line, green overlays)
def digitise_fig4():
    im = load("fig4_components.png"); r, g, b = im[..., 0], im[..., 1], im[..., 2]
    # the line in this composite is a desaturated steel blue, e.g. (93,121,132)
    blue = (b - r > 18) & (b - g > 5) & (r < 200) & (b > 90)
    green = (g > r + 40) & (g > b + 40)
    panels = {"UFO": (21, 164, [20, 10, 0, -10, -20]), "Ap": (216, 361, [150, 100, 50, 0, -50, -100])}
    xg = grid_cols(im, 70, 560, 30, 150, frac=0.3); assert len(xg) == 7, xg
    xa, xb, xres = lin_from_ticks(xg, DOY_TICKS)
    out = {}
    for name, (top, bot, ticks) in panels.items():
        yg = grid_rows(im, 95, 540, top + 1, bot - 1, lo=150, frac=0.3)
        assert len(yg) == len(ticks), (name, yg)
        ya, yb, yres = lin_from_ticks(yg, np.array(ticks, float))
        mask = blue & ~green
        xpx, ypx = trace(im, mask, 70, 562, top + 1, bot)
        # columns occupied by the green marker lines are treated as gaps
        gcols = np.where(green[top:bot, :].any(axis=0))[0]
        ypx[np.isin(xpx, gcols)] = np.nan
        grid, vals, cov = to_doy_grid(xpx, ypx, xa, xb, ya, yb)
        out[name] = vals
        print(f"fig4 {name}: x-res {xres:.2f} d, y-res {yres:.2f} %, coverage {cov:.3f}, range [{vals.min():.1f}, {vals.max():.1f}]")
    return pd.DataFrame(out, index=pd.Index(np.arange(365), name="doy"))

# ---------------------------------------------------------------- Fig 7 (Google Trends, hourly, past 7 days)
def digitise_fig7(n_points=169):
    """169 hourly points (7 x 24 + 1). The x axis is anchored on the circled maximum, which the tooltip
    in the screenshot dates to 31 Oct 2021, 06:00 observer time (UTC+3); px/h follows from the line extent."""
    im = load("fig7_trends7d.png"); r, g, b = im[..., 0], im[..., 1], im[..., 2]
    blue = (b > 200) & (b - r > 100) & (g > 90) & (g < 190) & (r < 130)
    y0, y1 = 375, 575
    yg = grid_rows(im, 70, 1100, 380, 560, lo=200, hi=245, frac=0.5)
    dark = (r < 200) & (abs(r - g) < 10) & (abs(g - b) < 10)
    base = np.mean(np.where(dark[560:580, 100:1000].mean(axis=1) > 0.5)[0] + 560)
    top100 = yg.min(); nominal = np.array([100, 75, 50, 25, 0], float)
    ticks_px = np.concatenate([yg, [base]])
    approx = 100 * (base - ticks_px) / (base - top100)
    ticks_val = nominal[np.abs(approx[:, None] - nominal[None, :]).argmin(axis=1)]
    ya, yb, yres = lin_from_ticks(ticks_px, ticks_val)
    xpx, ypx = trace(im, blue, 40, 1110, y0, y1)
    ok = ~np.isnan(ypx); vals = ya * ypx + yb
    xs = xpx[ok]; x_start, x_end = xs.min() + 2, xs.max() - 2       # ~2px round caps
    pxph = (x_end - x_start) / (n_points - 1)
    x_peak = np.mean(xpx[ok & (vals >= 98)])                          # the circled maximum (index 100)
    hx = x_peak + (np.arange(n_points) - round((x_peak - x_start) / pxph)) * pxph
    hv = np.interp(hx, xpx[ok], vals[ok])
    # occlusion by the tooltip box: hourly points whose +-3 px neighbourhood holds no line pixel
    occ = np.array([not ok[(xpx >= x - 3) & (xpx <= x + 3)].any() for x in hx])
    i_peak = int(np.argmin(np.abs(hx - x_peak)))
    print(f"fig7: y-res {yres:.2f}, line x-range {x_start:.0f}-{x_end:.0f}, px/h {pxph:.3f}, x_peak {x_peak:.1f} -> point {i_peak}, occluded hours {occ.sum()}")
    return hx, hv, occ, i_peak

if __name__ == "__main__":
    f5 = digitise_fig5(); f4 = digitise_fig4()
    comp = pd.concat({"fig5": f5, "fig4": f4}, axis=1); comp.columns = [f"{a}_{b}" for a, b in comp.columns]
    comp.to_csv("digitised_components.csv")
    hx, hv, gap, i_peak = digitise_fig7()
    # time axis anchored on the circled maximum: 31 Oct 2021 06:00 observer time (UTC+3)
    t_peak = pd.Timestamp("2021-10-31 06:00")
    t = t_peak + pd.to_timedelta(np.arange(len(hv)) - i_peak, unit="h")
    tr = pd.DataFrame({"ufo_hourly": hv, "occluded": gap}, index=t); tr.index.name = "local_time_utc+3"
    tr.to_csv("digitised_trends7d.csv")
    print(tr.loc["2021-10-31 00:00":"2021-10-31 12:00"])
    print("argmax:", tr["ufo_hourly"].idxmax(), tr["ufo_hourly"].max())
