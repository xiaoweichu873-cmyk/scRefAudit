# Use from an existing single-cell workflow

Start with a fixed `cost` matrix whose rows are query cells and columns are
reference cells. `reference_labels[j]` must identify cost column `j`.
Keep cell identifiers aligned with the query object.

```r
fit <- scRefAudit::audit_fit(
  cost = cost,
  labels = reference_labels,
  epsilon = 0.05,
  cell_ids = rownames(cost)
)

# Class weights, not individual cells, are changed by these scenarios.
proposal <- scRefAudit::propose_scenarios(
  fit, budget = 32L, metric = "tv", radius = 0.01,
  method = "directed", seed = 1L
)
result <- scRefAudit::audit_compositions(fit, proposal$proportions, order = 1L)
review <- scRefAudit::review_cells(fit, top_n = 100L)
```

Optional metadata attachment, if `query_object` is a Seurat object and its
cells match the matrix row identifiers:

```r
metadata <- result$cells
rownames(metadata) <- metadata$cell_id
metadata$cell_id <- NULL
stopifnot(all(rownames(metadata) %in% colnames(query_object)))
query_object <- SeuratObject::AddMetaData(query_object, metadata = metadata)
```

The cost representation, distance scale and epsilon should be chosen for the
original annotation task. The package does not estimate biological reference
proportions or distinguish absolute cell counts from relative cell abundance.
Avoid using query truth labels to select these settings during evaluation.

`saveRDS(result, ...)` and `saveRDS(review, ...)` are supported because they are
ordinary R results. A fit contains a live Python object: recreate the fit in a
new R session rather than reloading its pointer from RDS.
