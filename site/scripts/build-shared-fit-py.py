"""Build the shared fitted AnnData ONCE per site build (not once per reference page).

Every executed example on the Python reference pages reads site/_cache/HumanDLPFC-fit.h5ad
instead of refitting. Same data and genes as the R shared fit (build-shared-fit-r.R writes
_cache/shared-genes.txt; the raw data come from _cache/HumanDLPFC.h5ad, exported by
export-h5ad.R). Cheap settings: sample_p=0.05, df_tps=2, tol=1e-2. Git-ignored; delete _cache/
to force a refit.
"""
import time
import warnings
from pathlib import Path

import anndata as ad

import spanormpy as sn

warnings.filterwarnings("ignore")
cache = Path(__file__).resolve().parent.parent / "_cache"
t0 = time.time()
adata = ad.read_h5ad(cache / "HumanDLPFC.h5ad")
genes = (cache / "shared-genes.txt").read_text().split()
adata = adata[:, genes].copy()
sn.spanorm(adata, sample_p=0.05, df_tps=2, tol=1e-2, random_state=36, verbose=False)
sn.spanorm_svg(adata, verbose=False)
sn.spanorm_pca(adata)
# not polished here: polishing invalidates the SVG results, and the polish_spanorm and
# is_polished examples polish the (unpolished) shared fit themselves
adata.write_h5ad(cache / "HumanDLPFC-fit.h5ad")
print(f"shared Python fit: {adata.n_obs} spots x {adata.n_vars} genes, {time.time() - t0:.0f} s")
