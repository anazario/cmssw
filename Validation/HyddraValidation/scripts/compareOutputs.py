#!/usr/bin/env python3
"""
Event-by-event comparison of two HyddraSVsEXOAnalyzer ROOT outputs.

Matches events by (run, lumi, event), then compares vertex yields, kinematics,
and selection flags to verify the two producers give equivalent results.

Usage:
    python3 compareOutputs.py old.root new.root -l Old New -o comparison.pdf
    python3 compareOutputs.py old.root new.root --tree hyddraVal/Events
"""

import argparse
import sys
import warnings
warnings.filterwarnings('ignore', category=UserWarning, module='numpy')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import uproot
from matplotlib.backends.backend_pdf import PdfPages

TREE_PATH = "hyddraVal/Events"

SCALAR_BRANCHES = ['run', 'lumi', 'event', 'nLeptonTracks', 'nHyddraSV', 'Event_MET']

VECTOR_BRANCHES = [
    'HyddraSV_dxy', 'HyddraSV_dxyErr', 'HyddraSV_dxySig',
    'HyddraSV_normChi2', 'HyddraSV_pt', 'HyddraSV_eta',
    'HyddraSV_phi', 'HyddraSV_mass', 'HyddraSV_p',
    'HyddraSV_cosTheta', 'HyddraSV_decayAngle', 'HyddraSV_dR',
    'HyddraSV_passDisambiguation', 'HyddraSV_passIsolation',
    'HyddraSV_trk1Pt', 'HyddraSV_trk1Dxy', 'HyddraSV_trk1DxySig',
    'HyddraSV_trk1NormChi2',
    'HyddraSV_trk2Pt', 'HyddraSV_trk2Dxy', 'HyddraSV_trk2DxySig',
    'HyddraSV_trk2NormChi2',
]


def load_events(path, tree_path):
    f = uproot.open(path)

    # Auto-detect tree if needed — strip ROOT cycle suffix for lookup
    def _strip(k): return k.split(';')[0]
    keys_stripped = {_strip(k): k for k in f.keys()}
    if tree_path not in keys_stripped:
        candidates = [k for k in f.keys() if 'Events' in k or 'tree' in k.lower()]
        if not candidates:
            print(f"  Available keys: {f.keys()}")
            sys.exit(f"ERROR: tree '{tree_path}' not found in {path}")
        tree_path = candidates[0]
        print(f"  Auto-detected tree: {tree_path}")

    t = f[tree_path]
    avail = set(_strip(k) for k in t.keys())

    required = ['run', 'lumi', 'event']
    missing_req = [b for b in required if b not in avail]
    if missing_req:
        sys.exit(f"ERROR: required branches missing in {path}: {missing_req}")

    optional = [b for b in SCALAR_BRANCHES + VECTOR_BRANCHES
                if b not in required and b in avail]

    want = required + optional
    arrays = t.arrays(filter_name=want, library="np")
    n = len(arrays['run'])
    events = {}
    for i in range(n):
        key = (int(arrays['run'][i]), int(arrays['lumi'][i]), int(arrays['event'][i]))
        events[key] = {b: arrays[b][i] for b in want}
    return events


def match_events(ev1, ev2):
    common = sorted(set(ev1) & set(ev2))
    only1 = len(ev1) - len(common)
    only2 = len(ev2) - len(common)
    return common, only1, only2


def match_vertices(dxy_a, dxy_b, max_dxy_diff=0.5):
    """
    Greedy nearest-neighbour match on dxy. Each old vertex is matched to the
    closest unmatched new vertex (by |dxy_a - dxy_b|) within max_dxy_diff cm.
    Vertices with no close partner are left unmatched.
    """
    if len(dxy_a) == 0 or len(dxy_b) == 0:
        return [], list(range(len(dxy_a))), list(range(len(dxy_b)))

    used_b  = set()
    pairs   = []
    unmatched_a = []

    for ia in np.argsort(-dxy_a):            # process old vertices large-dxy first
        best_ib, best_d = None, max_dxy_diff
        for ib in range(len(dxy_b)):
            if ib in used_b:
                continue
            d = abs(float(dxy_a[ia]) - float(dxy_b[ib]))
            if d < best_d:
                best_d, best_ib = d, ib
        if best_ib is not None:
            pairs.append((ia, best_ib))
            used_b.add(best_ib)
        else:
            unmatched_a.append(ia)

    unmatched_b = [i for i in range(len(dxy_b)) if i not in used_b]
    return pairs, unmatched_a, unmatched_b


def collect_vertex_pairs(common_keys, ev1, ev2, max_dxy_diff=0.5):
    """Collect paired vertex quantities across all matched events."""
    fields = [b for b in VECTOR_BRANCHES
              if b in ev1[common_keys[0]] and b in ev2[common_keys[0]]]
    paired = {f: ([], []) for f in fields}
    unmatched_counts = {'extra_old': 0, 'extra_new': 0, 'matched': 0}
    nsv_pairs = []

    for key in common_keys:
        e1, e2 = ev1[key], ev2[key]
        dxy_key = 'HyddraSV_dxy'
        n1 = int(e1['nHyddraSV']) if 'nHyddraSV' in e1 else len(np.asarray(e1.get(dxy_key, [])))
        n2 = int(e2['nHyddraSV']) if 'nHyddraSV' in e2 else len(np.asarray(e2.get(dxy_key, [])))
        nsv_pairs.append((n1, n2))

        dxy1 = np.asarray(e1['HyddraSV_dxy'], dtype=float) if 'HyddraSV_dxy' in e1 else np.array([])
        dxy2 = np.asarray(e2['HyddraSV_dxy'], dtype=float) if 'HyddraSV_dxy' in e2 else np.array([])

        pairs, ua, ub = match_vertices(dxy1, dxy2, max_dxy_diff)
        unmatched_counts['matched'] += len(pairs)
        unmatched_counts['extra_old'] += len(ua)
        unmatched_counts['extra_new'] += len(ub)

        for ia, ib in pairs:
            for f in fields:
                arr1 = np.asarray(e1[f])
                arr2 = np.asarray(e2[f])
                if ia < len(arr1) and ib < len(arr2):
                    paired[f][0].append(float(arr1[ia]))
                    paired[f][1].append(float(arr2[ib]))

    for f in fields:
        paired[f] = (np.array(paired[f][0]), np.array(paired[f][1]))

    return paired, np.array(nsv_pairs), unmatched_counts


# ── Plotting helpers ────────────────────────────────────────────────────────

def residual(a, b):
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(np.abs(a) > 1e-9, (b - a) / a, b - a)


def hist_compare(ax, v1, v2, bins, label1, label2, xlabel, logy=False):
    h1, edges = np.histogram(v1[np.isfinite(v1)], bins=bins)
    h2, _     = np.histogram(v2[np.isfinite(v2)], bins=bins)
    cx = 0.5 * (edges[:-1] + edges[1:])
    ax.step(edges[:-1], h1, where='post', label=label1, color='steelblue', linewidth=1.5)
    ax.step(edges[:-1], h2, where='post', label=label2, color='tomato',    linewidth=1.5, linestyle='--')
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Vertices')
    ax.legend(fontsize=8)
    if logy:
        ax.set_yscale('log')
    ax.set_xlim(edges[0], edges[-1])


def scatter_nsv(ax, nsv_pairs, label1, label2):
    n1 = nsv_pairs[:, 0]
    n2 = nsv_pairs[:, 1]
    vmax = max(n1.max(), n2.max()) if len(n1) else 0
    vals = np.arange(0, vmax + 2)
    counts = np.zeros((len(vals), len(vals)), dtype=int)
    for a, b in zip(n1, n2):
        ia, ib = int(a), int(b)
        if ia < len(vals) and ib < len(vals):
            counts[ia, ib] += 1
    agree = np.sum(n1 == n2)
    frac  = agree / len(n1) if len(n1) else 0.0
    im = ax.imshow(counts.T, origin='lower', aspect='auto',
                   extent=[-0.5, len(vals)-0.5, -0.5, len(vals)-0.5],
                   cmap='Blues', norm=matplotlib.colors.LogNorm(vmin=1))
    plt.colorbar(im, ax=ax, label='Events')
    diag = np.arange(len(vals))
    ax.plot(diag, diag, 'k--', linewidth=0.8, alpha=0.5)
    ax.set_xlabel(f'nHyddraSV ({label1})')
    ax.set_ylabel(f'nHyddraSV ({label2})')
    ax.set_title(f'nSV agreement: {frac:.1%} ({agree}/{len(n1)} events)')


def residual_panel(ax_res, v1, v2, label, bins=50, pct_range=(-0.5, 0.5)):
    rel = residual(v1, v2)
    finite = rel[np.isfinite(rel)]
    ax_res.hist(finite, bins=bins, range=pct_range, color='mediumpurple', histtype='stepfilled', alpha=0.7)
    ax_res.axvline(0, color='k', linewidth=0.8)
    ax_res.set_xlabel(f'(New−Old)/Old  [{label}]')
    ax_res.set_ylabel('Vertex pairs')
    mu = np.mean(finite) if len(finite) else 0
    sd = np.std(finite)  if len(finite) else 0
    ax_res.set_title(f'μ={mu:.3f}  σ={sd:.3f}', fontsize=9)


def flag_bar(ax, v1, v2, name, label1, label2):
    a1 = np.asarray(v1, dtype=bool)
    a2 = np.asarray(v2, dtype=bool)
    categories = ['Both pass', 'Only Old', 'Only New', 'Both fail']
    counts = [
        np.sum( a1 &  a2),
        np.sum( a1 & ~a2),
        np.sum(~a1 &  a2),
        np.sum(~a1 & ~a2),
    ]
    colors = ['#4caf50', '#2196f3', '#ff9800', '#9e9e9e']
    bars = ax.bar(categories, counts, color=colors)
    for bar, c in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(c), ha='center', va='bottom', fontsize=8)
    total = sum(counts)
    agree = counts[0] + counts[3]
    ax.set_title(f'{name}: {agree/total:.1%} agreement' if total else name, fontsize=9)
    ax.set_ylabel('Vertex pairs')


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('file1', help='First ROOT file (old/reference)')
    parser.add_argument('file2', help='Second ROOT file (new/test)')
    parser.add_argument('-l', '--labels', nargs=2, default=['Old', 'New'], metavar=('LABEL1', 'LABEL2'))
    parser.add_argument('-o', '--output', default='comparison.pdf')
    parser.add_argument('--tree', default=TREE_PATH, help=f'TTree path (default: {TREE_PATH})')
    parser.add_argument('--dxy-match', type=float, default=0.5,
                        metavar='CM', help='Max |dxy_old - dxy_new| for vertex matching in cm (default: 0.5)')
    args = parser.parse_args()

    label1, label2 = args.labels

    print(f"Loading {args.file1} ...")
    ev1 = load_events(args.file1, args.tree)
    print(f"  {len(ev1)} events")

    print(f"Loading {args.file2} ...")
    ev2 = load_events(args.file2, args.tree)
    print(f"  {len(ev2)} events")

    branches1 = set(next(iter(ev1.values())).keys())
    branches2 = set(next(iter(ev2.values())).keys())
    only_in_1 = branches1 - branches2 - {'run','lumi','event'}
    only_in_2 = branches2 - branches1 - {'run','lumi','event'}
    if only_in_1: print(f"  Branches only in {label1}: {sorted(only_in_1)}")
    if only_in_2: print(f"  Branches only in {label2}: {sorted(only_in_2)}")

    common, only1, only2 = match_events(ev1, ev2)
    print(f"\nEvent matching:")
    print(f"  Common events : {len(common)}")
    print(f"  Only in {label1:6s}: {only1}")
    print(f"  Only in {label2:6s}: {only2}")

    if not common:
        sys.exit("No common events — check that both files ran on the same input.")

    paired, nsv_pairs, ucounts = collect_vertex_pairs(common, ev1, ev2, args.dxy_match)

    total_vtx = ucounts['matched'] + ucounts['extra_old'] + ucounts['extra_new']
    print(f"\nVertex matching (across {len(common)} common events):")
    print(f"  Matched pairs : {ucounts['matched']}")
    print(f"  Only in {label1:6s}: {ucounts['extra_old']}")
    print(f"  Only in {label2:6s}: {ucounts['extra_new']}")
    match_frac = ucounts['matched'] / total_vtx if total_vtx else 0
    print(f"  Match fraction: {match_frac:.1%}")

    # ── nLeptonTracks comparison ──────────────────────────────────────────
    has_ntk = all('nLeptonTracks' in ev1[k] and 'nLeptonTracks' in ev2[k] for k in common)
    if has_ntk:
        ntk1 = np.array([ev1[k]['nLeptonTracks'] for k in common])
        ntk2 = np.array([ev2[k]['nLeptonTracks'] for k in common])
        ntk_agree = np.sum(ntk1 == ntk2)
        print(f"\nnLeptonTracks agreement: {ntk_agree}/{len(common)} events ({ntk_agree/len(common):.1%})")
        if ntk_agree < len(common):
            diff = ntk2 - ntk1
            print(f"  Mean difference (New-Old): {np.mean(diff):.2f}  Max: {np.max(np.abs(diff))}")
    else:
        ntk1 = ntk2 = np.array([])
        print("\nnLeptonTracks: not present in one or both files, skipping")

    # ── nHyddraSV comparison ──────────────────────────────────────────────
    nsv_agree = np.sum(nsv_pairs[:, 0] == nsv_pairs[:, 1])
    print(f"\nnHyddraSV agreement: {nsv_agree}/{len(common)} events ({nsv_agree/len(common):.1%})")

    # ── Per-variable residuals summary ────────────────────────────────────
    CONT_VARS = ['HyddraSV_dxy', 'HyddraSV_normChi2', 'HyddraSV_pt', 'HyddraSV_mass',
                 'HyddraSV_trk1Pt', 'HyddraSV_trk2Pt', 'HyddraSV_trk1Dxy', 'HyddraSV_trk2Dxy']
    print(f"\nVertex-level residuals (New−Old)/Old for matched pairs [{ucounts['matched']} pairs]:")
    print(f"  {'Variable':<30} {'Mean':>10} {'Std':>10} {'|>1%|':>8} {'|>5%|':>8}")
    print(f"  {'-'*68}")
    for var in CONT_VARS:
        if var not in paired or len(paired[var][0]) == 0:
            continue
        v1, v2 = paired[var]
        rel = residual(v1, v2)
        finite = rel[np.isfinite(rel)]
        if len(finite) == 0:
            continue
        gt1  = np.sum(np.abs(finite) > 0.01)
        gt5  = np.sum(np.abs(finite) > 0.05)
        print(f"  {var:<30} {np.mean(finite):>+10.4f} {np.std(finite):>10.4f} "
              f"{gt1:>8} {gt5:>8}")

    # ── Flag agreement ────────────────────────────────────────────────────
    for flag in ['HyddraSV_passDisambiguation', 'HyddraSV_passIsolation']:
        if flag in paired and len(paired[flag][0]):
            v1f, v2f = paired[flag]
            agree = np.sum(np.array(v1f, dtype=bool) == np.array(v2f, dtype=bool))
            print(f"\n{flag}: {agree}/{len(v1f)} pairs agree ({agree/len(v1f):.1%})")

    # ── Difference diagnostics ───────────────────────────────────────────
    print(f"\n── Where do the extra {ucounts['extra_new']} New-only vertices live? ──")
    nsv_arr = np.array(nsv_pairs)
    same_mask  = nsv_arr[:,0] == nsv_arr[:,1]
    more_mask  = nsv_arr[:,1] >  nsv_arr[:,0]
    fewer_mask = nsv_arr[:,1] <  nsv_arr[:,0]
    print(f"  Events same nSV                    : {same_mask.sum()} ({same_mask.mean():.1%})")
    avg_extra = f"{(nsv_arr[more_mask,1]-nsv_arr[more_mask,0]).mean():.1f}" if more_mask.any() else "n/a"
    print(f"  Events New has MORE   (New > Old)   : {more_mask.sum()} ({more_mask.mean():.1%}), "
          f"avg extra = {avg_extra}")
    avg_missing = f"{(nsv_arr[fewer_mask,0]-nsv_arr[fewer_mask,1]).mean():.1f}" if fewer_mask.any() else "n/a"
    print(f"  Events New has FEWER  (New < Old)   : {fewer_mask.sum()} ({fewer_mask.mean():.1%}), "
          f"avg missing = {avg_missing}")

    # dxy distribution of unmatched vertices
    unmatched_old_dxy = []
    unmatched_new_dxy = []
    for key in common:
        e1, e2 = ev1[key], ev2[key]
        dxy1 = np.asarray(e1.get('HyddraSV_dxy', []), dtype=float)
        dxy2 = np.asarray(e2.get('HyddraSV_dxy', []), dtype=float)
        _, ua, ub = match_vertices(dxy1, dxy2, args.dxy_match)
        for ia in ua:
            if ia < len(dxy1):
                unmatched_old_dxy.append(float(dxy1[ia]))
        for ib in ub:
            if ib < len(dxy2):
                unmatched_new_dxy.append(float(dxy2[ib]))
    if unmatched_old_dxy:
        u = np.array(unmatched_old_dxy)
        print(f"\n  Old-only vertex dxy: "
              f"mean={u.mean():.2f} cm  median={np.median(u):.2f} cm  "
              f"max={u.max():.2f} cm  <1cm: {(u<1).mean():.1%}")
    if unmatched_new_dxy:
        u = np.array(unmatched_new_dxy)
        print(f"\n  New-only vertex dxy: "
              f"mean={u.mean():.2f} cm  median={np.median(u):.2f} cm  "
              f"max={u.max():.2f} cm  <1cm: {(u<1).mean():.1%}")

    # Flag-flip breakdown on matched pairs
    print(f"\n── Flag disagreements on {ucounts['matched']} matched pairs ──")
    for flag in ['HyddraSV_passDisambiguation', 'HyddraSV_passIsolation']:
        if flag not in paired or len(paired[flag][0]) == 0:
            continue
        f1 = np.asarray(paired[flag][0], dtype=bool)
        f2 = np.asarray(paired[flag][1], dtype=bool)
        flip_on  = np.sum(~f1 &  f2)   # Old=fail, New=pass
        flip_off = np.sum( f1 & ~f2)   # Old=pass, New=fail
        print(f"  {flag.replace('HyddraSV_',''):30s}: "
              f"Old=pass→New=fail: {flip_off:4d}   Old=fail→New=pass: {flip_on:4d}")

    # ── Plots ─────────────────────────────────────────────────────────────
    print(f"\nWriting plots to {args.output} ...")
    with PdfPages(args.output) as pdf:

        # Page 1: Event-level summary
        fig, axes = plt.subplots(2, 2, figsize=(12, 9))
        fig.suptitle(f'Event-level comparison: {label1} vs {label2}', fontsize=13)

        # nSV scatter
        scatter_nsv(axes[0, 0], nsv_pairs, label1, label2)

        # nSV distributions
        vmax = int(nsv_pairs.max()) + 1 if len(nsv_pairs) else 5
        bins_nsv = np.arange(-0.5, vmax + 1.5)
        axes[0, 1].hist(nsv_pairs[:, 0], bins=bins_nsv, alpha=0.6, label=label1, color='steelblue')
        axes[0, 1].hist(nsv_pairs[:, 1], bins=bins_nsv, alpha=0.6, label=label2, color='tomato')
        axes[0, 1].set_xlabel('nHyddraSV'); axes[0, 1].set_ylabel('Events')
        axes[0, 1].legend(); axes[0, 1].set_yscale('log')

        # nLeptonTracks
        if len(ntk1) > 0:
            vmax_tk = max(int(ntk1.max()), int(ntk2.max())) + 1
            bins_tk = np.arange(-0.5, vmax_tk + 1.5)
            axes[1, 0].hist(ntk1, bins=bins_tk, alpha=0.6, label=label1, color='steelblue')
            axes[1, 0].hist(ntk2, bins=bins_tk, alpha=0.6, label=label2, color='tomato')
            axes[1, 0].set_xlabel('nLeptonTracks'); axes[1, 0].set_ylabel('Events')
            axes[1, 0].legend()
        else:
            axes[1, 0].set_visible(False)

        # nSV difference
        diff_nsv = nsv_pairs[:, 1] - nsv_pairs[:, 0]
        dmin, dmax = int(diff_nsv.min()), int(diff_nsv.max())
        axes[1, 1].hist(diff_nsv, bins=np.arange(dmin - 0.5, dmax + 1.5),
                        color='mediumpurple', histtype='stepfilled', alpha=0.8)
        axes[1, 1].axvline(0, color='k', linewidth=0.8)
        axes[1, 1].set_xlabel(f'nHyddraSV ({label2}) − ({label1})')
        axes[1, 1].set_ylabel('Events')
        axes[1, 1].set_title(f'{nsv_agree/len(common):.1%} events agree')

        plt.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # Page 2: Kinematic distributions
        DIST_VARS = [
            ('HyddraSV_dxy',      np.linspace(0,   30, 51),  'Vertex dxy [cm]',        True),
            ('HyddraSV_normChi2', np.linspace(0,   10, 51),  'Vertex χ²/ndof',         True),
            ('HyddraSV_pt',       np.linspace(0,  100, 51),  'Vertex pT [GeV]',        True),
            ('HyddraSV_mass',     np.linspace(0,    5, 51),  'Vertex mass [GeV]',      True),
            ('HyddraSV_dR',       np.linspace(0,    5, 51),  'Track ΔR',               False),
            ('HyddraSV_cosTheta', np.linspace(-1,   1, 51),  'cos θ (vtx→PV vs p)',    False),
            ('HyddraSV_trk1Pt',   np.linspace(0,   50, 51),  'Track1 pT [GeV]',        True),
            ('HyddraSV_trk2Pt',   np.linspace(0,   50, 51),  'Track2 pT [GeV]',        True),
        ]

        # Collect all SV quantities from both files for distribution plots
        def all_sv(ev_dict, branch):
            vals = []
            for e in ev_dict.values():
                if branch in e:
                    vals.extend(np.asarray(e[branch], dtype=float).tolist())
            return np.array(vals)

        fig, axes = plt.subplots(4, 2, figsize=(12, 14))
        fig.suptitle(f'Kinematic distributions: {label1} vs {label2}', fontsize=13)
        for ax, (var, bins, xlabel, logy) in zip(axes.flat, DIST_VARS):
            if var not in paired:
                ax.set_visible(False); continue
            v_all1 = all_sv(ev1, var)
            v_all2 = all_sv(ev2, var)
            hist_compare(ax, v_all1, v_all2, bins, label1, label2, xlabel, logy)
        plt.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # Page 3: Residuals for matched vertex pairs
        RES_VARS = [
            ('HyddraSV_dxy',     'dxy',      50, (-0.05, 0.05)),
            ('HyddraSV_normChi2','χ²/ndof',  50, (-0.05, 0.05)),
            ('HyddraSV_pt',      'pT',        50, (-0.05, 0.05)),
            ('HyddraSV_mass',    'mass',      50, (-0.05, 0.05)),
            ('HyddraSV_trk1Pt',  'trk1 pT',  50, (-0.05, 0.05)),
            ('HyddraSV_trk2Pt',  'trk2 pT',  50, (-0.05, 0.05)),
            ('HyddraSV_trk1Dxy', 'trk1 dxy', 50, (-0.1,  0.1)),
            ('HyddraSV_trk2Dxy', 'trk2 dxy', 50, (-0.1,  0.1)),
        ]
        fig, axes = plt.subplots(4, 2, figsize=(12, 14))
        fig.suptitle(f'Vertex residuals (New−Old)/Old  [{ucounts["matched"]} matched pairs]', fontsize=13)
        for ax, (var, label, nbins, rng) in zip(axes.flat, RES_VARS):
            if var not in paired or len(paired[var][0]) == 0:
                ax.set_visible(False); continue
            v1, v2 = paired[var]
            residual_panel(ax, v1, v2, label, nbins, rng)
        plt.tight_layout()
        pdf.savefig(fig); plt.close(fig)

        # Page 4: Flag agreement
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        fig.suptitle(f'Selection flag agreement: {label1} vs {label2}', fontsize=13)
        for ax, flag, name in zip(axes,
                                   ['HyddraSV_passDisambiguation', 'HyddraSV_passIsolation'],
                                   ['passDisambiguation',           'passIsolation']):
            if flag in paired and len(paired[flag][0]):
                flag_bar(ax, paired[flag][0], paired[flag][1], name, label1, label2)
            else:
                ax.set_visible(False)
        plt.tight_layout()
        pdf.savefig(fig); plt.close(fig)

    print(f"Done. Output written to {args.output}")


if __name__ == '__main__':
    main()
