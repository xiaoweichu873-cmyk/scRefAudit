"""Check both balanced-OT figures using only the three commit-pinned public XLSX files."""
import argparse
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
KEYS = ["task", "m", "eps", "B"]
NAMES = {"hao_f0": "Hao P1", "hao_f1": "Hao P2", "kang_f0": "Kang fold 1", "kang_f1": "Kang fold 2"}
RANKS = ["confidence", "entropy", "margin", "linear_radius"]


def half_up(value):
    return str(Decimal(str(value)).quantize(Decimal("0.01"), rounding=ROUND_HALF_UP))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--books", type=Path, default=ROOT / "workbooks")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    def sheet(part, name):
        return pd.read_excel(args.books / f"Supplementary_Data_Part{part}.xlsx", sheet_name=name)
    raw = sheet(2, "P2_31_timings")
    if raw.duplicated(KEYS + ["method", "rep"]).any() or len(raw) != 180:
        raise ValueError("Unexpected number of unique trial records.")
    ok = raw[raw.status == "ok"].groupby(KEYS + ["method"]).total_s.agg(["median", "count"])
    completed = ok["median"].where(ok["count"] == 3).unstack()
    ratios = completed[["reuse_AA2","reuse_chord"]].div(completed["quadratic_chord"], axis=0).sort_index()
    published = sheet(2, "P2_25_ratios_to_quadratic").set_index(KEYS).sort_index()
    pd.testing.assert_index_equal(ratios.index,published.index)
    np.testing.assert_allclose(ratios,published[ratios.columns],rtol=1e-12,equal_nan=True)
    ratios.to_csv(args.out / "figure5_ratios.csv")

    counts = sheet(3, "P3_21_review_margin_primary").set_index(["task","references"]).sort_index()
    table = sheet(1,"Table_4")
    expected = counts[RANKS].rename(columns=dict(zip(RANKS,["Confidence","Entropy","Winning margin","Linear radius"])))
    for (task,size), row in expected.iterrows():
        printed = table[(table.Task == NAMES[task]) & (table.References == size)]
        np.testing.assert_allclose(row.to_numpy(dtype=float),printed[row.index].iloc[0].to_numpy(dtype=float))
    np.testing.assert_allclose(expected.mean(), table[table.Task=="Mean"][expected.columns].iloc[0].to_numpy(dtype=float))
    counts[RANKS].to_csv(args.out / "table1.csv")
    print("Historical status counts:",raw.status.value_counts().to_dict())
    for key, group in ratios.items():
        finite=group.dropna()
        print(key,"complete",len(finite),"wins",int((finite>1).sum()),"median",half_up(finite.median()),
              "range",half_up(finite.min()),half_up(finite.max()))
    print("Top-100 means:",{c:half_up(v) for c,v in counts[RANKS].mean().items()})

    fig,axes=plt.subplots(1,2,figsize=(10,5.7),layout="constrained")
    for ax,col,label in zip(axes,ratios.columns,["A   AA2 / quadratic shared curvature","B   Shared curvature / quadratic shared curvature"]):
        s=ratios[col].dropna().sort_values()
        ticks=[f"{NAMES[t]}  ·  {m}  /  {e:g}  /  {b}" for t,m,e,b in s.index]
        ax.barh(range(len(s)),s,color=["#b26342" if v<1 else "#187578" for v in s])
        ax.set_yticks(range(len(s)),ticks,fontsize=8)
        ax.axvline(1,color="#333",ls="--",lw=.8)
        ax.set_xlabel("Ratio of median times")
        ax.set_title(label,fontsize=10)
    fig.savefig(args.out/"figure5_from_workbooks.pdf");plt.close(fig)
    fig,ax=plt.subplots(figsize=(5.2,4),layout="constrained")
    display=["confidence","margin","linear_radius"]
    for _,row in counts.iterrows():
        ax.plot(range(3),row[display],color="#acafb2",marker="o",markersize=3,lw=1)
    means=counts[display].mean()
    ax.plot(range(3),means,marker="o",lw=2.4,color="#187578")
    for x,val in enumerate(means):
        ax.annotate(half_up(val),(x,val),xytext=(0,8),textcoords="offset points",ha="center",color="#155955")
    ax.set_xticks(range(3),["Confidence","Winning margin","Linearized radius"])
    ax.set_ylim(0,113);ax.set_ylabel("Sensitive cells per 100 reviewed")
    ax.set_title("Figure 2B   Review yield on random scenarios",loc="left",fontsize=11)
    fig.savefig(args.out/"figure2B_from_workbooks.pdf");plt.close(fig)

if __name__=="__main__":
    main()
