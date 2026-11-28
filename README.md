# scRefAudit supporting data and analysis code

Research materials for *Efficient reference composition sensitivity analysis for single-cell annotation*. The manuscript and its extended methods contain the interpretation, mathematical derivations and figure legends. A Chinese guide is available in [README_zh.md](README_zh.md).

| Content | Location |
|---|---|
| Three indexed workbooks with source values, timing records and numerical checks | [workbooks](workbooks/) |
| Tissue experiments for Figure 7 | [78 prepared matrices](matched_tv_20261118/prepared/) and [source and frozen result volumes](records/) |
| Replotting scripts and source records for main and supplementary figures | [records](records/), unpack as shown below |
| Cell ranking (Figure 2B / Table 1) and solver speed (Figure 5) | [Fixed balanced OT code, prepared inputs and replay commands](https://github.com/xiaoweichu873-cmyk/scRefAudit/archive/d728fff33f0ac58840f62e2d571764fef42477c3.zip) |
| Exploratory KL-unbalanced OT study | [Derivation and benchmarks](unbalanced/Supplementary_Addendum_Unbalanced_OT.pdf) and [versioned code, exact inputs and full results](https://github.com/xiaoweichu873-cmyk/scRefAudit/tree/7969588a5194bf6b8e02256992139ae92e4c22bc/manuscript_data/unbalanced_ot_20261123) |

Download the full snapshot and extract it. To unpack further source records into their original relative directories and verify them, run the following commands from the snapshot root. The matrices in `matched_tv_20261118/prepared/` are already in place. Extraction requires Python 3.9 or later and uses the standard library.

```sh
sha256sum -c SHA256SUMS
python3 expand_records.py
```

To **redraw figures from stored values**, install the plotting dependencies in `figure_sources/requirements.txt` and run `python3 figure_sources/regenerate_figures.py`. These historical plotting scripts name their outputs using the previous numbering. Original figures 1, 2, 3, 4, 5, 6 and 7 correspond to manuscript figures 1, 3, 4, 5, 6, 7 and 2. The workbook sheets Table_1, Table_2, Table_3 and Table_4 correspond to manuscript Tables 2, 3, 4 and 1. Supplementary numbering is unchanged.

To **rerun the matched total-variation tissue analysis from prepared inputs**, see `matched_tv_20261118/README.md` after extraction. The 78 inputs cover five donor rotations across liver myeloid, lung epithelial and lung mesenchymal compartments. For the balanced immune analysis, the separate fixed download above documents fresh computation, archived-table checks and timing trials. Regeneration of the PCA inputs from original GEO expression matrices was outside the verification. The GEO accessions cited in the paper are GSE96583, GSE164378, GSE149614 and GSE135893.

The original research files retain their attribution and bundled notices. Source code for the released R and Python interfaces is in the [software repository](https://github.com/xiaoweichu873-cmyk/scRefAudit).
