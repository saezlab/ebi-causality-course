# Data snapshots for notebook 01

These files preserve the external data used by
[01-context-specific-networks.ipynb](../../01-context-specific-networks.ipynb).
They were downloaded on **7 October 2026**. The notebook still uses its original
online loaders; it does not read these snapshots.

| File | Origin and preparation |
| --- | --- |
| [`adata.h5ad`](adata.h5ad) | Six human hepatic stellate cell RNA-seq samples from [GEO series GSE151251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151251), loaded with `decoupler.ds.hsctgfb()`. Decoupler reads the series' [processed count table](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE151251&format=file&file=GSE151251%5FHSCs%5FCtrl%2Evs%2EHSCs%5FTGFb%2Ecounts%2Etsv%2Egz), retains one row per gene name, transposes it to samples × genes, and adds sample and condition metadata. Saved as AnnData (6 samples × 58,674 genes). |
| [`collectri.parquet`](collectri.parquet) | Human CollecTRI TF–target regulons (42,990 rows), loaded with `decoupler.op.collectri(organism="human")`. This Decoupler version downloads [`CollecTRI_regulons.csv` from Zenodo record 8192729](https://zenodo.org/records/8192729), then formats resource and reference fields, removes missing and duplicate interactions, and returns a DataFrame. The loader does **not** fetch this file from the OmniPath server. |
| [`unique_receptors.parquet`](unique_receptors.parquet) | Unique target-side gene symbols (1,201 names) from `omnipath.interactions.LigRecExtra.get(genesymbols=True)`. The notebook uses these as candidate receptor annotations in its plot. The original OmniPath interaction records were reduced to a sorted, single-column DataFrame (`target_genesymbol`); interaction-level source evidence is not retained here. |
| [`df_ridden.parquet`](df_ridden.parquet) | RIDDEN receptor–gene model matrix (229 × 978), read from [`ridden_model_matrix.csv` at the pinned RIDDEN_tool commit `aa3c9d9`](https://github.com/basvaat/RIDDEN_tool/blob/aa3c9d9880f135e95a079865ee20ebec984cf54d/ridden_model/ridden_model_matrix.csv). Saved as a DataFrame with receptor names in the index and gene symbols in the columns. |
| [`pkn.parquet`](pkn.parquet) | SIGNOR-supported interactions returned by `omnipath.interactions.OmniPath.get(databases=["SIGNOR"], genesymbols=True)`, then filtered to `consensus_direction == True`, as in the notebook (60,921 rows). This is an OmniPath query result, not a direct export of the whole SIGNOR database. |

## Attribution

- **GSE151251:** V. H. Shah and U. Yaqoob, contributors, *RNA-seq from human
  hepatic stellate cells treated with TGF-β*, NCBI Gene Expression Omnibus,
  [GSE151251](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151251)
  (2020). The GEO series lists no associated publication.
- **CollecTRI:** S. Müller-Dott *et al.*, “Expanding the coverage of regulons
  from high-confidence prior knowledge for accurate estimation of transcription
  factor activities”, *Nucleic Acids Research* 51 (2023), 10934–10949,
  [doi:10.1093/nar/gkad841](https://doi.org/10.1093/nar/gkad841). Cite the
  deposited [CollecTRI dataset, doi:10.5281/zenodo.8192729](https://doi.org/10.5281/zenodo.8192729)
  for the file itself.
- **OmniPath and LigRecExtra:** D. Türei *et al.*, “Integrated intra- and
  intercellular signaling knowledge for multicellular omics analysis”,
  *Molecular Systems Biology* 17 (2021), e9923,
  [doi:10.15252/msb.20209923](https://doi.org/10.15252/msb.20209923).
  OmniPath integrates multiple source resources; consult its
  [resource information](https://omnipathdb.org/info) for the underlying
  LigRecExtra contributors and their attribution and licence terms.
- **SIGNOR:** L. Licata *et al.*, “SIGNOR 2.0, the SIGnaling Network Open
  Resource 2.0: 2019 update”, *Nucleic Acids Research* 48 (2020), D504–D510,
  [doi:10.1093/nar/gkz949](https://doi.org/10.1093/nar/gkz949). The interaction
  snapshot above was obtained through OmniPath, so cite both SIGNOR and
  OmniPath when using it.
- **RIDDEN:** S. Barsi *et al.*, “RIDDEN: Data-driven inference of receptor
  activity from transcriptomic data”, *PLOS Computational Biology* 21 (2025),
  e1013188, [doi:10.1371/journal.pcbi.1013188](https://doi.org/10.1371/journal.pcbi.1013188).
  The matrix comes from the authors' [RIDDEN_tool repository](https://github.com/basvaat/RIDDEN_tool)
  at the commit linked above.

The four `.parquet` files are pandas DataFrames. `adata.h5ad` retains the counts
and sample metadata in AnnData format. The snapshots were produced with
Decoupler 2.2.0, OmniPath Python client 1.0.12, AnnData 0.13.4, pandas 3.0.6
and PyArrow 25.0.1.
