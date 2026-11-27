# Manuscript supplementary result data

For **Efficient reference composition sensitivity analysis for single-cell annotation**.

These three Excel workbooks provide the tabulated results, figure source values,
timing records and numerical diagnostics accompanying the manuscript. They are
provided here for reviewer download because converting all worksheets into the
submission PDF produced thousands of pages.

## Download

| File | Contents |
| --- | --- |
| [Supplementary Data Part 1](Supplementary_Data_Part1.xlsx?raw=true) | Main and supplementary tables, matched-total-variation comparisons, and mathematical validation |
| [Supplementary Data Part 2](Supplementary_Data_Part2.xlsx?raw=true) | Immune-cell benchmarks, numerical checks, and protein measurements |
| [Supplementary Data Part 3](Supplementary_Data_Part3.xlsx?raw=true) | Tissue analyses, review prioritization, and additional diagnostics |

Download each workbook and open it in Excel, LibreOffice, or a compatible reader.
Start with the **Index** worksheet to locate the source file and corresponding
data worksheet. GitHub may not display a preview for these binary files; use
the download links above.

The **Documents** sheets preserve explanatory records. The **Original_text**
sheets retain exact source files as ordered JSON-quoted text chunks; their
reconstruction preserves decimal text, empty fields, and original line endings.
These source-record sheets are intended for machine reading, not PDF printing.

## Balanced OT experiments

Download the [versioned code and prepared inputs for the review and timing benchmarks](balanced_ot_v5/README.md). The evidence corresponds to Figure 2B / Table 1 and Figure 5 of the Technical Note.

## Scope and integrity

The workbooks are byte-for-byte copies of the files prepared for submission.
No records, failed runs, or numerical values were removed or changed for this
deposit, and no calculations were rerun. File sizes and SHA-256 checksums are
listed in [manifest.json](manifest.json).

The accompanying manuscript and Supplementary Material document contain the
methods, mathematical analyses, supplementary figures, and interpretation.
These workbooks do not constitute the complete experimental archive: they do
not supply all processed input matrices, full scenario predictions, or comparator
build materials. The main repository supplies the scRefAudit implementation.

This GitHub deposit does not yet have a permanent archive DOI.
