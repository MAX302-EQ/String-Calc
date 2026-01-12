import math
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt

st.set_page_config(page_title="PV String Sizing", layout="wide")
st.title("PV String Sizing")

# -----------------------------
# Utilities
# -----------------------------
def to_frac_per_c(x, mode):
    return x / 100.0 if mode == "% per °C" else x

def tcell_from_noct(t_amb, noct, g):
    # Tcell = Tamb + (NOCT - 20)/800 * G
    return t_amb + ((noct - 20.0) / 800.0) * g

def module_at_tcell(voc_stc, vmp_stc, isc_stc, imp_stc,
                    beta_voc, beta_vmp, alpha_isc, alpha_imp,
                    tcell, bif_gain_pct, t_stc=25.0):
    dT = tcell - t_stc
    voc = voc_stc * (1.0 + beta_voc * dT)
    vmp = vmp_stc * (1.0 + beta_vmp * dT)
    isc = isc_stc * (1.0 + alpha_isc * dT)
    imp = imp_stc * (1.0 + alpha_imp * dT)

    # bifacial gain: boost current only (simple + practical)
    gain = 1.0 + bif_gain_pct / 100.0
    isc *= gain
    imp *= gain
    return voc, vmp, isc, imp

def mppt_effective(vmin, vmax, bottom_margin_pct, top_margin_pct):
    # Example: 200–1000 -> 220–950
    vmin_eff = vmin * (1.0 + bottom_margin_pct / 100.0)
    vmax_eff = vmax * (1.0 - top_margin_pct / 100.0)
    return vmin_eff, vmax_eff

def string_range(vmp_at_tmax, voc_at_tmin, vmin_eff, vmax_eff, abs_max_vdc):
    vals = [vmp_at_tmax, voc_at_tmin, vmin_eff, vmax_eff, abs_max_vdc]
    if any(v is None for v in vals):
        return None, None
    if vmp_at_tmax <= 0 or voc_at_tmin <= 0:
        return None, None

    n_min = math.ceil(vmin_eff / vmp_at_tmax)
    v_limit = min(vmax_eff, abs_max_vdc)
    n_max = math.floor(v_limit / voc_at_tmin)
    return n_min, n_max

def check_mppt(n, voc_at_tmin, vmp_at_tmax,
               imp_at_tmax, isc_at_tmax,
               vmin_eff, vmax_eff, abs_max_vdc,
               mppt_i_max, mppt_isc_max,
               inputs_in_parallel):
    voc_str = n * voc_at_tmin
    vmp_str = n * vmp_at_tmax

    imp_total = inputs_in_parallel * imp_at_tmax
    isc_total = inputs_in_parallel * isc_at_tmax

    ok_vmin = vmp_str >= vmin_eff
    ok_vmax = voc_str <= vmax_eff
    ok_abs  = voc_str <= abs_max_vdc
    ok_imp  = imp_total <= mppt_i_max
    ok_isc  = True if (mppt_isc_max == 0 or mppt_isc_max is None) else (isc_total <= mppt_isc_max)

    return {
        "String Voc @ Tmin (V)": round(voc_str, 2),
        "String Vmp @ Tmax (V)": round(vmp_str, 2),
        "MPPT Imp total (A)": round(imp_total, 2),
        "MPPT Isc total (A)": round(isc_total, 2),
        "Vmp >= MPPT min?": ok_vmin,
        "Voc <= MPPT max?": ok_vmax,
        "Voc <= Abs max Vdc?": ok_abs,
        "Imp_total <= Imax?": ok_imp,
        "Isc_total <= Iscmax?": ok_isc,
        "PASS": (ok_vmin and ok_vmax and ok_abs and ok_imp and ok_isc)
    }

def bar_chart(title, labels, values):
    fig = plt.figure()
    x = range(len(labels))
    plt.bar(x, values)
    plt.xticks(x, labels, rotation=20)
    plt.ylabel("Value")
    plt.title(title)
    st.pyplot(fig)

def fuse_check(n_parallel, isc_stc, i_string_fuse_max):
    # Always return SAME keys for all cases
    backfeed = (max(n_parallel, 1) - 1) * isc_stc
    fuse_min = 1.25 * isc_stc
    fuse_max = i_string_fuse_max

    fuse_required = (n_parallel > 1) and (backfeed > fuse_max)
    window_ok = fuse_min <= fuse_max

    if n_parallel <= 1:
        fuse_label = "No fuse"
        note = "Single string"
    else:
        if fuse_required and window_ok:
            fuse_label = f"Fuse {fuse_min:.1f}–{fuse_max:.1f} A"
            note = "Fusing required"
        elif fuse_required and not window_ok:
            fuse_label = "Fuse NOT possible"
            note = "❌ 1.25×Isc_STC > max series fuse"
        else:
            fuse_label = "No fuse"
            note = "Fusing not required"

    return {
        "Fuse required?": fuse_required,
        "Backfeed (n-1)*Isc_STC (A)": round(backfeed, 2),
        "Fuse min 1.25*Isc_STC (A)": round(fuse_min, 2),
        "Fuse max (A)": round(fuse_max, 2),
        "Fuse window OK?": window_ok,
        "Fuse label": fuse_label,
        "Fuse note": note,
    }


# -----------------------------
# Sidebar: site + margins
# -----------------------------
with st.sidebar:
    st.header("Site Inputs")
    tmin = st.number_input("Minimum site temperature (°C)", value=5.0, step=1.0)
    tmax = st.number_input("Maximum site temperature (°C)", value=45.0, step=1.0)
    g = st.number_input("Irradiance G (W/m²) (single value)", value=1000.0, step=50.0)

    st.header("MPPT Safety Margins")
    bottom_margin = st.number_input("Bottom margin (%) (MPPT min +)", value=10.0, step=1.0)
    top_margin = st.number_input("Top margin (%) (MPPT max -)", value=5.0, step=1.0)

    st.header("String Fuse Check")
    i_string_fuse_max = st.number_input("Module max series fuse rating (A)", value=25.0, step=1.0)

# -----------------------------
# Module input (manual but compact)
# -----------------------------
st.subheader("1) PV Module")

c1, c2, c3, c4 = st.columns(4)
with c1:
    module_name = st.text_input("Module name", value="Module-1")
    voc_stc = st.number_input("Voc STC (V)", value=49.5, step=0.1)
    vmp_stc = st.number_input("Vmp STC (V)", value=41.5, step=0.1)

with c2:
    isc_stc = st.number_input("Isc STC (A)", value=14.0, step=0.1)
    imp_stc = st.number_input("Imp STC (A)", value=13.2, step=0.1)
    noct = st.number_input("NOCT (°C)", value=45.0, step=0.5)

with c3:
    coef_mode = st.selectbox("Coefficient units", ["% per °C", "fraction per °C"])
    beta_voc_in = st.number_input("β Voc (usually -)", value=-0.29, step=0.01)
    beta_vmp_in = st.number_input("β Vmp (usually -)", value=-0.35, step=0.01)

with c4:
    alpha_isc_in = st.number_input("α Isc (usually +)", value=0.04, step=0.01)
    alpha_imp_in = st.number_input("α Imp (usually +)", value=0.03, step=0.01)
    bif_gain = st.number_input("Bifacial gain (%) (0 if mono)", value=0.0, step=1.0)

beta_voc = to_frac_per_c(beta_voc_in, coef_mode)
beta_vmp = to_frac_per_c(beta_vmp_in, coef_mode)
alpha_isc = to_frac_per_c(alpha_isc_in, coef_mode)
alpha_imp = to_frac_per_c(alpha_imp_in, coef_mode)

if voc_stc <= 0 or vmp_stc <= 0 or isc_stc <= 0 or imp_stc <= 0:
    st.error("Please enter valid module STC values (Voc, Vmp, Isc, Imp > 0).")
    st.stop()

# -----------------------------
# Inverter + MPPT cards
# -----------------------------
st.divider()
st.subheader("2) Inverter")

inv1, inv2, inv3 = st.columns(3)
with inv1:
    inv_name = st.text_input("Inverter name", value="Inverter")
with inv2:
    abs_max_vdc = st.number_input("Absolute max DC voltage (V)", value=1000.0, step=10.0)
with inv3:
    n_mppt = st.number_input("Number of MPPT", value=2, min_value=1, step=1)

if "mppt_cfg" not in st.session_state:
    st.session_state.mppt_cfg = []

while len(st.session_state.mppt_cfg) < int(n_mppt):
    st.session_state.mppt_cfg.append({
        "name": f"MPPT-{len(st.session_state.mppt_cfg)+1}",
        "vmin": 200.0,
        "vmax": 1000.0,
        "imax": 26.0,
        "iscmax": 0.0,
        "inputs": 1,
        "inputs_max": 2
    })
while len(st.session_state.mppt_cfg) > int(n_mppt):
    st.session_state.mppt_cfg.pop()

st.markdown("### MPPT settings")

for i in range(int(n_mppt)):
    cfg = st.session_state.mppt_cfg[i]

    with st.expander(f"MPPT {i+1}: {cfg['name']}", expanded=(i == 0)):
        a, b, c, d = st.columns(4)

        cfg["name"] = a.text_input(
            "MPPT name",
            value=cfg["name"],
            key=f"mppt_name_{i}"
        )

        cfg["vmin"] = b.number_input(
            "MPPT min (V)",
            value=float(cfg["vmin"]),
            step=10.0,
            key=f"mppt_vmin_{i}"
        )

        cfg["vmax"] = c.number_input(
            "MPPT max (V)",
            value=float(cfg["vmax"]),
            step=10.0,
            key=f"mppt_vmax_{i}"
        )

        cfg["imax"] = d.number_input(
            "MPPT max current Imax (A)",
            value=float(cfg["imax"]),
            step=1.0,
            key=f"mppt_imax_{i}"
        )

        e, f, g2 = st.columns(3)

        cfg["iscmax"] = e.number_input(
            "MPPT Isc max (A) (0 = ignore)",
            value=float(cfg["iscmax"]),
            step=1.0,
            key=f"mppt_iscmax_{i}"
        )

        cfg["inputs_max"] = f.number_input(
            "Physical inputs available",
            value=int(cfg["inputs_max"]),
            min_value=1,
            step=1,
            key=f"mppt_inputsmax_{i}"
        )

        # 🔒 Enforce inputs_used ≤ physical inputs
        max_allowed = int(cfg["inputs_max"])
        current_val = min(int(cfg["inputs"]), max_allowed)

        cfg["inputs"] = g2.number_input(
            "Inputs used (parallel strings)",
            value=current_val,
            min_value=1,
            max_value=max_allowed,
            step=1,
            key=f"mppt_inputs_{i}"
        )
        st.session_state.mppt_cfg[i] = cfg

# -----------------------------
# Calculations at Tmin & Tmax (your wording) + show ALL values
# -----------------------------
st.divider()
st.subheader("3) Module values at Minimum site temperature & Maximum site temperature")

tcell_min = tmin
tcell_max = tcell_from_noct(tmax, noct, g)

voc_tmin, vmp_tmin, isc_tmin, imp_tmin = module_at_tcell(
    voc_stc, vmp_stc, isc_stc, imp_stc,
    beta_voc, beta_vmp, alpha_isc, alpha_imp,
    tcell_min, bif_gain
)
voc_tmax, vmp_tmax, isc_tmax, imp_tmax = module_at_tcell(
    voc_stc, vmp_stc, isc_stc, imp_stc,
    beta_voc, beta_vmp, alpha_isc, alpha_imp,
    tcell_max, bif_gain
)

# Show as a clear table (most user-friendly)
df_wc = pd.DataFrame({
    "Parameter": ["Voc (V)", "Isc (A)", "Vmp (V)", "Imp (A)"],
    "Value": [
        round(voc_tmin, 3),   # Voc at Tmin
        round(isc_tmin, 3),   # Isc at Tmin
        round(vmp_tmax, 3),   # Vmp at Tmax
        round(imp_tmax, 3),   # Imp at Tmax
    ],
    "Temperature used": [
        f"Tmin = {tmin:.1f}°C (Tcell={tcell_min:.1f}°C)",
        f"Tmin = {tmin:.1f}°C (Tcell={tcell_min:.1f}°C)",
        f"Tmax = {tmax:.1f}°C (Tcell={tcell_max:.1f}°C)",
        f"Tmax = {tmax:.1f}°C (Tcell={tcell_max:.1f}°C)",
    ]
})

st.dataframe(df_wc, use_container_width=True)

# Helpful quick metrics (optional, but nice)
m1, m2, m3, m4 = st.columns(4)
m1.metric("Voc @ Tmin (V)", f"{voc_tmin:.2f}")
m2.metric("Vmp @ Tmax (V)", f"{vmp_tmax:.2f}")
m3.metric("Isc @ Tmin (A)", f"{isc_tmin:.2f}")
m4.metric("Imp @ Tmax (A)", f"{imp_tmax:.2f}")



# -----------------------------
# MPPT sizing + choose n
# -----------------------------
st.divider()
st.subheader("5) MPPT-wise string sizing + checks")

ranges = []
for cfg in st.session_state.mppt_cfg:
    vmin_eff, vmax_eff = mppt_effective(cfg["vmin"], cfg["vmax"], bottom_margin, top_margin)
    nmin, nmax = string_range(vmp_tmax, voc_tmin, vmin_eff, vmax_eff, abs_max_vdc)
    ranges.append((cfg["name"], vmin_eff, vmax_eff, nmin, nmax))

range_df = pd.DataFrame(ranges, columns=["MPPT", "Eff MPPT min (V)", "Eff MPPT max (V)", "n_min", "n_max"])
range_df["Feasible"] = range_df.apply(lambda r: (r["n_min"] is not None and r["n_max"] is not None and r["n_max"] >= r["n_min"]), axis=1)
st.dataframe(range_df, use_container_width=True)

feasible = range_df[range_df["Feasible"] == True]
if feasible.empty:
    st.error("No feasible string length found. Check MPPT ranges/margins, inverter max Vdc, or module inputs.")
    st.stop()

global_n_min = int(feasible["n_min"].max())
global_n_max = int(feasible["n_max"].min())

st.markdown("### Choose modules per string")
if global_n_max < global_n_min:
    st.warning("No single n fits all MPPT rows together (ranges differ). You can still test a chosen n.")
    chosen_n = st.number_input("n (modules per string)", value=22, min_value=1, step=1)
else:
    chosen_n = st.slider("n (modules per string)", min_value=global_n_min, max_value=global_n_max,
                         value=int(round((global_n_min + global_n_max) / 2)), step=1)

# checks
check_rows = []
status_map = {}
for cfg in st.session_state.mppt_cfg:
    vmin_eff, vmax_eff = mppt_effective(cfg["vmin"], cfg["vmax"], bottom_margin, top_margin)
    inputs_used = int(cfg["inputs"])
    phys = int(cfg["inputs_max"])
    phys_warn = inputs_used > phys

    result = check_mppt(
        int(chosen_n),
        voc_tmin,
        vmp_tmax,
        imp_tmax,
        isc_tmin,
        vmin_eff,
        vmax_eff,
        abs_max_vdc,
        float(cfg["imax"]),
        float(cfg["iscmax"]),
        inputs_used
    )
 # ---- Fuse check (uses Isc at Tmin, same as your sizing current basis) ----
    fuse_res = fuse_check(inputs_used,isc_stc,i_string_fuse_max)
    result.update(fuse_res)
    result["MPPT"] = cfg["name"]
    result["Inputs used"] = inputs_used
    result["Inputs max"] = phys
    result["Inputs exceed?"] = phys_warn
    check_rows.append(result)
    status_map[cfg["name"]] = bool(result["PASS"])

check_df = pd.DataFrame(check_rows)
# ---- Fuse maps for diagram use ----
fuse_label_map = dict(zip(check_df["MPPT"], check_df["Fuse label"]))
fuse_req_map = dict(zip(check_df["MPPT"], check_df["Fuse required?"]))
cols = ["MPPT", "PASS", "Inputs used", "Inputs max", "Inputs exceed?",
        "String Voc @ Tmin (V)", "String Vmp @ Tmax (V)",
        "MPPT Imp total (A)", "MPPT Isc total (A)",
        "Fuse required?", "Backfeed (n-1)*Isc_STC (A)",
        "Fuse min 1.25*Isc_STC (A)", "Fuse max (A)", "Fuse window OK?",
        "Fuse label", "Fuse note"]
st.dataframe(check_df[cols], use_container_width=True)

passed = sum(1 for v in status_map.values() if v)
total = len(status_map)
if passed == total and not check_df["Inputs exceed?"].any():
    st.success(f"✅ PASS: n={chosen_n} OK for all MPPTs ({passed}/{total}).")
elif passed == total and check_df["Inputs exceed?"].any():
    st.warning(f"⚠️ Electrical checks PASS for n={chosen_n}, but some MPPTs exceed physical inputs per MPPT.")
else:
    st.error(f"❌ FAIL: n={chosen_n} fails on some MPPTs ({passed}/{total}).")




# -----------------------------
# Diagram (colourful, clear)
# -----------------------------
st.divider()
st.subheader("6) Diagram (MPPT → Strings)")

dot = []
dot.append("digraph PV {")
dot.append("rankdir=LR;")
dot.append('node [shape=box, style="rounded,filled", fontname="Arial"];')
dot.append('edge [color="#555555"];')

# Inverter node
dot.append(
    f'INV [label="{inv_name}\\nModule: {module_name}\\nAbs Max DC: {abs_max_vdc:.0f} V", '
    f'fillcolor="#cfe8ff", color="#2b7bbb"];'
)

for i, cfg in enumerate(st.session_state.mppt_cfg, start=1):
    mppt_id = f"MPPT{i}"
    vmin_eff, vmax_eff = mppt_effective(cfg["vmin"], cfg["vmax"], bottom_margin, top_margin)

    # MPPT node
    dot.append(
        f'{mppt_id} [label="{cfg["name"]}\\n'
        f'V: {vmin_eff:.0f}–{vmax_eff:.0f} V\\n'
        f'Imax: {cfg["imax"]:.1f} A\\n'
        f'Inputs: {cfg["inputs"]}", '
        f'fillcolor="#d9f2d9", color="#2e7d32"];'
    )

    dot.append(f"INV -> {mppt_id};")

    # String node (PASS / FAIL colouring)
    # ---- String nodes: show one box per input (parallel string) ----
    inputs_used = int(cfg["inputs"])
    is_pass = status_map.get(cfg["name"], False)

    fill = "#e6f4ea" if is_pass else "#fdecea"
    border = "#2e7d32" if is_pass else "#c62828"
    status = "PASS" if is_pass else "FAIL"

     # ---- String nodes: show one box per input (parallel string) ----
    inputs_used = int(cfg["inputs"])
    is_pass = status_map.get(cfg["name"], False)

    fill = "#e6f4ea" if is_pass else "#fdecea"
    border = "#2e7d32" if is_pass else "#c62828"
    status = "PASS" if is_pass else "FAIL"

    mppt_name = cfg["name"]
    fuse_required = bool(fuse_req_map.get(mppt_name, False))
    fuse_label = str(fuse_label_map.get(mppt_name, "No fuse"))

    # Create one string box per input
    for k in range(1, inputs_used + 1):
        s_id = f"S{i}_{k}"

        # String box
        dot.append(
            f'{s_id} [shape=note, label="String {k}/{inputs_used}\\n'
            f'n = {int(chosen_n)}\\n'
            f'Voc@Tmin = {int(chosen_n)*voc_tmin:.1f} V\\n'
            f'Vmp@Tmax = {int(chosen_n)*vmp_tmax:.1f} V\\n'
            f'{status}", '
            f'fillcolor="{fill}", color="{border}", fontcolor="black"];'
        )

        # Fuse box (only if required)
        if fuse_required:
            f_id = f"F{i}_{k}"
            dot.append(
                f'{f_id} [shape=box, style="rounded,filled", '
                f'fillcolor="#fff3cd", color="#a67c00", '
                f'label="{fuse_label}"];'
            )
            dot.append(f"{mppt_id} -> {f_id};")
            dot.append(f"{f_id} -> {s_id};")
        else:
            dot.append(f"{mppt_id} -> {s_id};")

dot.append("}")

st.graphviz_chart("\n".join(dot), use_container_width=True)