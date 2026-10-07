"""
Validate tinyspice two ways:

  1. against closed-form math (RC step, series RLC step), where the exact
     answer is known, including a time-step study that shows the trapezoidal
     rule is 2nd-order accurate (halve h -> error / 4) and backward Euler is
     1st-order (halve h -> error / 2);
  2. against ngspice, an established open-source SPICE, on nonlinear
     circuits (diode rectifier, switching buck converter).

Run from the repo root after building:   python tests/validate.py
Writes figures/*.png and tests/validation_results.json.
"""
import json
import os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
EXE = ROOT / ("tinyspice.exe" if os.name == "nt" else "tinyspice")
FIG = ROOT / "figures"
FIG.mkdir(exist_ok=True)

INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "axes.edgecolor": GRID, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.titlelocation": "left",
                     "axes.titleweight": "bold", "axes.titlesize": 11, "axes.labelcolor": INK2,
                     "xtick.color": INK2, "ytick.color": INK2, "legend.frameon": False, "lines.linewidth": 2})


def find_ngspice():
    for cand in (os.environ.get("NGSPICE"), shutil.which("ngspice_con"), shutil.which("ngspice"),
                 str(Path.home() / "tools/Spice64/bin/ngspice_con.exe")):
        if cand and Path(cand).exists():
            return cand
    return None


def read_csv(path):
    with open(path) as f:
        head = f.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return {h.lower(): data[:, k] for k, h in enumerate(head)}


def tiny(netlist_text, extra=()):
    with tempfile.TemporaryDirectory() as d:
        cir, out = Path(d) / "c.cir", Path(d) / "o.csv"
        cir.write_text(netlist_text)
        subprocess.run([str(EXE), str(cir), "-o", str(out), "--quiet", *extra], check=True)
        return read_csv(out)


def ngspice(netlist_text, vectors):
    exe = find_ngspice()
    if not exe:
        return None
    # same circuit; force ngspice's maximum step to the same print step
    text = re.sub(r"(?im)^\.tran\s+(\S+)\s+(\S+)\s*$", r".tran \1 \2 0 \1", netlist_text)
    text = re.sub(r"(?im)^\.end\s*$", "", text)
    with tempfile.TemporaryDirectory() as d:
        out = (Path(d) / "o.txt").as_posix()
        text += ("\n.options reltol=1e-4\n.control\nset wr_singlescale\nset wr_vecnames\nrun\n"
                 f"wrdata {out} {' '.join(vectors)}\nquit\n.endc\n.end\n")
        cir = Path(d) / "c.cir"
        cir.write_text(text)
        subprocess.run([exe, "-b", str(cir)], capture_output=True, check=True)
        data = np.loadtxt(out, skiprows=1)
    res = {"time": data[:, 0]}
    for k, v in enumerate(vectors, 1):
        res[v.lower()] = data[:, k]
    return res


results = {}

# ---------------------------------------------------------------- 1. RC step
RC = (ROOT / "examples/rc_step.cir").read_text()
r = tiny(RC)
t, v = r["time"], r["v(out)"]
exact = 1 - np.exp(-t / 1e-3)
results["rc_step_max_error_V"] = float(np.max(np.abs(v - exact)))

# order of accuracy: error vs time step, trapezoidal vs backward Euler
steps = [200e-6, 100e-6, 50e-6, 20e-6, 10e-6, 5e-6]
order = {"trap": [], "be": []}
for h in steps:
    txt = re.sub(r"(?im)^\.tran.*$", f".tran {h} 5m", RC)
    for name, extra in (("trap", ()), ("be", ("--be",))):
        rr = tiny(txt, extra)
        order[name].append(float(np.max(np.abs(rr["v(out)"] - (1 - np.exp(-rr["time"] / 1e-3))))))
slope = {k: float(np.polyfit(np.log(steps), np.log(errs), 1)[0]) for k, errs in order.items()}
results["order_of_accuracy"] = {"steps_s": steps, "max_error_V": order, "fitted_slope": slope}

# ------------------------------------------------------------ 2. series RLC
RLC = (ROOT / "examples/rlc.cir").read_text()
r2 = tiny(RLC)
R_, L_, C_ = 20.0, 10e-3, 100e-9
w0 = 1 / np.sqrt(L_ * C_)
z = R_ / 2 * np.sqrt(C_ / L_)
wd = w0 * np.sqrt(1 - z * z)
t2 = r2["time"]
exact2 = 1 - np.exp(-z * w0 * t2) * (np.cos(wd * t2) + z / np.sqrt(1 - z * z) * np.sin(wd * t2))
results["rlc_max_error_V"] = float(np.max(np.abs(r2["v(out)"] - exact2)))
results["rlc_peak_overshoot_sim_pct"] = float(100 * (r2["v(out)"].max() - 1))
results["rlc_peak_overshoot_theory_pct"] = float(100 * np.exp(-np.pi * z / np.sqrt(1 - z * z)))

# --------------------------------------------- 3./4. nonlinear vs ngspice
cmp = {}
for name, vec in (("rectifier", ["v(out)"]), ("buck", ["v(out)", "v(sw)"])):
    txt = (ROOT / f"examples/{name}.cir").read_text()
    a = tiny(txt)
    b = ngspice(txt, vec)
    if b is None:
        print("ngspice not found - skipping", name)
        continue
    tt = np.linspace(0, min(a["time"][-1], b["time"][-1]), 20000)
    ya = np.interp(tt, a["time"], a[vec[0]])
    yb = np.interp(tt, b["time"], b[vec[0]])
    rng = yb.max() - yb.min()
    cmp[name] = {"tt": tt, "ya": ya, "yb": yb}
    results[f"{name}_vs_ngspice_nrmse_pct"] = float(100 * np.sqrt(np.mean((ya - yb) ** 2)) / rng)
    results[f"{name}_vs_ngspice_max_err_pct"] = float(100 * np.max(np.abs(ya - yb)) / rng)
    if name == "buck":
        last = tt > tt[-1] - 0.5e-3
        results["buck_vout_avg_tiny_V"] = float(ya[last].mean())
        results["buck_vout_avg_ngspice_V"] = float(yb[last].mean())

(ROOT / "tests/validation_results.json").write_text(json.dumps(results, indent=2) + "\n")
for k, val in results.items():
    if not isinstance(val, dict):
        print(f"  {k:36s} {val:.4g}")
print("  order of accuracy (log-log slope):", {k: round(v, 2) for k, v in slope.items()})

# ------------------------------------------------------------------ figures
fig, ax = plt.subplots(2, 2, figsize=(12, 7.5))
a = ax[0, 0]
a.plot(t * 1e3, exact, color=INK2, lw=4, alpha=0.35, label="exact: 1 - e^(-t/RC)")
a.plot(t * 1e3, v, color=C1, lw=1.6, label="tinyspice")
a.set(title=f"RC step: max error {results['rc_step_max_error_V']:.1e} V", xlabel="time (ms)", ylabel="V")
a.legend()
a = ax[0, 1]
a.plot(t2 * 1e3, exact2, color=INK2, lw=4, alpha=0.35, label="exact (underdamped RLC)")
a.plot(t2 * 1e3, r2["v(out)"], color=C1, lw=1.6, label="tinyspice")
a.set(title=f"Series RLC: max error {results['rlc_max_error_V']:.1e} V", xlabel="time (ms)", ylabel="V")
a.legend()
for a, name, unit, title in ((ax[1, 0], "rectifier", "ms", "Diode rectifier vs ngspice"),
                             (ax[1, 1], "buck", "ms", "48 V -> 12 V buck converter vs ngspice")):
    if name not in cmp:
        continue
    d = cmp[name]
    a.plot(d["tt"] * 1e3, d["yb"], color=INK2, lw=4, alpha=0.35, label="ngspice 47")
    a.plot(d["tt"] * 1e3, d["ya"], color=C2, lw=1.4, label="tinyspice")
    a.set(title=f"{title}: NRMSE {results[name + '_vs_ngspice_nrmse_pct']:.3f} %",
          xlabel=f"time ({unit})", ylabel="output voltage (V)")
    a.legend()
fig.tight_layout()
fig.savefig(FIG / "validation.png", dpi=150)

fig, a = plt.subplots(figsize=(6.4, 4.4))
a.loglog(np.array(steps) * 1e6, order["be"], "o-", color=C2, label=f"backward Euler (slope {slope['be']:.2f})")
a.loglog(np.array(steps) * 1e6, order["trap"], "s-", color=C1, label=f"trapezoidal (slope {slope['trap']:.2f})")
a.set(xlabel="time step h (µs)", ylabel="max error vs exact (V)",
      title="Halve the step: error /2 (Euler) vs /4 (trapezoid)")
a.legend()
fig.tight_layout()
fig.savefig(FIG / "order_of_accuracy.png", dpi=150)
print("figures written to", FIG)
