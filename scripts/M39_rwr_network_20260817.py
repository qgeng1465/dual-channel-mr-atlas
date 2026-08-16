#!/usr/bin/env python3
# =============================================================================
# M39_rwr_network_20260817.py - STRING gene regulatory network RWR propagation
# =============================================================================
# Seeds = 106 known effector genes (results/grid/transcript_coloc_hits.csv).
# Random Walk with Restart (r=0.7) on STRING v12 human PPI; evaluate 15 candidates'
# network proximity to the known-gene module (is a candidate close to known genes?).
# Inputs:
#   <scratch>/string/9606.protein.links.v12.0.txt.gz
#   <scratch>/string/9606.protein.info.v12.0.txt.gz
#   results/grid/transcript_coloc_hits.csv       106 known genes
#   results/candidate15_replication_20260816.csv 15 candidates
# Outputs:
#   results/m39_rwr_20260817.csv         per candidate {score, rank, module, seed_frac, perm_z, perm_p}
#   results/m39_rwr_summary_20260817.csv network summary
#   results/figures/20260817_FigS_rwr_network.png, 20260817_FigS_rwr_rank.png
# Honesty caveat: STRING is a database-derived interactome (confidence-weighted,
#   directionless); RWR proximity != functional validation.
# =============================================================================
import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

REPO = "<repo-root>"
SCR = "<scratch>"
STRING = os.path.join(SCR, "string")
FIG = os.path.join(REPO, "results", "figures")
os.makedirs(FIG, exist_ok=True)

RESTART = 0.7
TOL = 1e-8
MAXITER = 1000
N_PERM = 100
SEED_RNG = 42


def rwr(G, seeds, restart=RESTART, tol=TOL, maxiter=MAXITER):
    nodes = list(G.nodes())
    idx = {n: i for i, n in enumerate(nodes)}
    n = len(nodes)
    W = nx.to_scipy_sparse_array(G, weight="weight", format="csr").astype(float)
    deg = np.asarray(W.sum(axis=0)).ravel()
    deg[deg == 0] = 1.0
    A = W / deg  # column stochastic transition
    r0 = np.zeros(n)
    s = [idx[x] for x in seeds if x in idx]
    if not s:
        raise ValueError("no seeds in graph")
    r0[s] = 1.0 / len(s)
    p = r0.copy()
    for _ in range(maxiter):
        pn = (1 - restart) * (A @ p) + restart * r0
        if np.linalg.norm(pn - p, ord=1) < tol:
            p = pn
            break
        p = pn
    return dict(zip(nodes, p))


def main():
    print("== M39 RWR: STRING network propagation ==")
    known = pd.read_csv(os.path.join(REPO, "results/grid/transcript_coloc_hits.csv"))
    cand = pd.read_csv(os.path.join(REPO, "results/candidate15_replication_20260816.csv"))
    known_syms = sorted(set(known["symbol"].dropna().astype(str)))
    cand_syms = sorted(set(cand["symbol"].dropna().astype(str)))
    cand_syms = [g for g in cand_syms if g not in set(known_syms)]
    print(f"known seeds={len(known_syms)} candidates={len(cand_syms)}")

    info = pd.read_csv(os.path.join(STRING, "9606.protein.info.v12.0.txt.gz"),
                       sep="\t", usecols=["#string_protein_id", "preferred_name"])
    info.columns = ["string_id", "name"]
    sym2id = dict(zip(info["name"], info["string_id"]))
    print(f"STRING proteins={len(sym2id)}")

    links = pd.read_csv(os.path.join(STRING, "9606.protein.links.v12.0.txt.gz"),
                        sep=" ", compression="gzip",
                        dtype={"protein1": str, "protein2": str, "combined_score": np.int32})
    links["weight"] = links["combined_score"] / 1000.0
    print(f"raw edges={len(links)}")

    out_rows, summary = [], []
    for thr in (400, 700):
        el = links[links["combined_score"] >= thr][["protein1", "protein2", "weight"]]
        G = nx.from_pandas_edgelist(el, "protein1", "protein2", edge_attr="weight")
        G = nx.Graph(G)
        print(f"  thr={thr}: nodes={G.number_of_nodes()} edges={G.number_of_edges()}")

        seeds_all = [sym2id[g] for g in known_syms if g in sym2id]
        seeds = [sid for sid in seeds_all if sid in G]   # 种子须在阈值子图内（thr=700 会丢孤立蛋白）
        seeds_dropped = len(seeds_all) - len(seeds)
        cand_map = {g: sym2id[g] for g in cand_syms if g in sym2id}
        unmapped = [g for g in cand_syms if g not in sym2id]
        if not seeds:
            print("  ! no seeds mapped; skip")
            continue
        print(f"  seeds mapped={len(seeds)}/{len(known_syms)} (dropped {seeds_dropped} not in graph) "
              f"candidates mapped={len(cand_map)}/{len(cand_syms)} unmapped={unmapped}")

        score = rwr(G, seeds)
        node_list = list(G.nodes())
        rank_of = {n: i + 1 for i, n in enumerate(sorted(node_list, key=lambda n: -score[n]))}

        rng = np.random.default_rng(SEED_RNG)
        cand_scores_obs = np.array([score.get(sym2id[g], 0.0) for g in cand_syms if g in sym2id])
        stat_obs = cand_scores_obs.max() if len(cand_scores_obs) else 0.0
        null_stats = []
        for _ in range(N_PERM):
            rs = rng.choice(node_list, size=len(seeds), replace=False)
            sc = rwr(G, list(rs))
            cs = np.array([sc.get(sym2id[g], 0.0) for g in cand_syms if g in sym2id])
            null_stats.append(cs.max() if len(cs) else 0.0)
        null_stats = np.array(null_stats)
        perm_p = (np.sum(null_stats >= stat_obs) + 1) / (N_PERM + 1)
        perm_z = (stat_obs - null_stats.mean()) / (null_stats.std() + 1e-12)
        print(f"  RWR max-candidate stat={stat_obs:.4g} null_mean={null_stats.mean():.4g} "
              f"perm_p={perm_p:.4f} z={perm_z:.2f}")

        comm = nx.algorithms.community.louvain_communities(G, weight="weight", seed=SEED_RNG)
        node2mod = {}
        for mi, c in enumerate(comm):
            for n in c:
                node2mod[n] = mi
        seed_set = set(seeds)
        mod_info = {}
        for mi, c in enumerate(comm):
            f = len(set(c) & seed_set) / max(len(c), 1)
            mod_info[mi] = (len(c), f)

        for g in cand_syms:
            sid = sym2id.get(g)
            if sid is None:
                out_rows.append(dict(symbol=g, threshold=thr, in_network=False,
                                     rwr_score=np.nan, rank=np.nan, module=np.nan,
                                     module_n=np.nan, module_seed_frac=np.nan,
                                     perm_p=np.nan, perm_z=np.nan))
                continue
            mi = node2mod.get(sid, -1)
            mn, mf = mod_info.get(mi, (np.nan, np.nan))
            out_rows.append(dict(symbol=g, threshold=thr, in_network=True,
                                 rwr_score=score.get(sid, np.nan), rank=rank_of.get(sid, np.nan),
                                 module=mi, module_n=mn, module_seed_frac=mf,
                                 perm_p=perm_p, perm_z=perm_z))
        summary.append(dict(threshold=thr, nodes=G.number_of_nodes(), edges=G.number_of_edges(),
                            seeds_mapped=len(seeds), candidates_mapped=len(cand_map),
                            unmapped=";".join(unmapped), max_cand_stat=stat_obs,
                            null_mean=null_stats.mean(), null_sd=null_stats.std(),
                            perm_p=perm_p, perm_z=perm_z, n_modules=len(comm)))

        if thr == 700:
            keep = set(seeds) | {sid for sid in cand_map.values() if sid in G}
            nb = set()
            for n in keep:
                for nn in G.neighbors(n):
                    if G[n][nn]["weight"] >= 0.9:
                        nb.add(nn)
            sub = G.subgraph(keep | nb)
            sub = sub.subgraph(sorted(sub.nodes(), key=lambda n: -score[n])[:400])
            id2sym = {v: k for k, v in sym2id.items()}
            H = nx.relabel_nodes(sub, {n: id2sym.get(n, n) for n in sub.nodes()})
            pos = nx.spring_layout(H, seed=SEED_RNG, weight="weight", k=0.5, iterations=80)
            mod_color = {}
            for n in sub.nodes():
                mod_color[n] = f"C{node2mod.get(n, -1) % 20}"
            node_colors = [mod_color[n] for n in sub.nodes()]
            fig, ax = plt.subplots(figsize=(13, 11))
            nx.draw_networkx_edges(H, pos, ax=ax, alpha=0.15, edge_color="0.5", width=0.5)
            nx.draw_networkx_nodes(H, pos, ax=ax, node_size=22, node_color=node_colors, alpha=0.8)
            cand_nodes = [c for c in cand_map.values() if c in sub.nodes()]
            nx.draw_networkx_nodes(H, pos, ax=ax, nodelist=[id2sym[c] for c in cand_nodes],
                                   node_size=260, node_shape="*", node_color="#c00000", alpha=1.0)
            seed_sub = [s for s in seeds if s in sub.nodes()]
            nx.draw_networkx_nodes(H, pos, ax=ax, nodelist=[id2sym[s] for s in seed_sub],
                                   node_size=90, node_shape="o", node_color="#1a3a6b",
                                   edgecolors="black", linewidths=0.6, alpha=0.95)
            labels = {n: n for n in sub.nodes() if n in [id2sym[c] for c in cand_nodes]}
            nx.draw_networkx_labels(H, pos, labels=labels, font_size=9, font_color="#c00000", font_weight="bold")
            ax.legend(handles=[Line2D([0], [0], marker='*', color='w', markerfacecolor='#c00000',
                                      markersize=16, label='15 candidates'),
                               Line2D([0], [0], marker='o', color='w', markerfacecolor='#1a3a6b',
                                      markeredgecolor='black', markersize=10, label='106 known seeds')],
                      loc='lower right', fontsize=11, frameon=True)
            ax.set_title(f"STRING v12 PPI (score >= 0.9 edges), RWR from 106 known effector genes\n"
                         f"candidate max-score = {stat_obs:.3g}  perm p = {perm_p:.3f}  z = {perm_z:.2f}",
                         fontsize=13)
            ax.axis("off")
            fig.savefig(os.path.join(FIG, "20260817_FigS_rwr_network.png"), dpi=300, bbox_inches="tight")
            plt.close(fig)
            print("  saved network figure")

    df = pd.DataFrame(out_rows)
    df.to_csv(os.path.join(REPO, "results/m39_rwr_20260817.csv"), index=False)
    pd.DataFrame(summary).to_csv(os.path.join(REPO, "results/m39_rwr_summary_20260817.csv"), index=False)

    d400 = df[df["threshold"] == 400].sort_values("rwr_score", ascending=False)
    if len(d400):
        fig, ax = plt.subplots(figsize=(9, 5.5))
        y = np.arange(len(d400))[::-1]
        sc = d400["rwr_score"].fillna(0).values
        ax.barh(y, sc, color="#1a3a6b", alpha=0.85)
        ax.set_yticks(y); ax.set_yticklabels(d400["symbol"], fontsize=10)
        ax.set_xlabel("RWR score (restart=0.7, STRING v12, score>=400)")
        ax.set_title(f"15 candidate effector genes: network proximity to 106 known genes\n"
                     f"perm p = {d400['perm_p'].iloc[0]:.4f}  z = {d400['perm_z'].iloc[0]:.2f}")
        ax.axvline(0, color="0.4", lw=0.8)
        for i, p in enumerate(d400["rwr_score"].fillna(0)):
            rk = d400['rank'].iloc[i]
            ax.text(p + 0.0005, y[i], f"  rank {int(rk) if pd.notna(rk) else '-'}", va="center", fontsize=8)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, "20260817_FigS_rwr_rank.png"), dpi=300, bbox_inches="tight")
        plt.close(fig)
        print("  saved rank figure")

    print("== DONE M39 ==")
    print(df[df["threshold"] == 400][["symbol", "in_network", "rwr_score", "rank",
                                      "module_seed_frac", "perm_p"]].to_string(index=False))


if __name__ == "__main__":
    main()
