# WASP documentation index

## Current package documentation

These files describe the implementation in `src/wasp/` for the v2.0.0
pre-release:

1. [`../README.md`](../README.md) — installation, capabilities, minimal usage,
   sources, matching, plotting, and citation.
2. [`../README_archi.md`](../README_archi.md) — module responsibilities,
   dependencies, call graphs, and compatibility boundaries.
3. [`processing_guide.md`](processing_guide.md) — current I/O contracts, units,
   shapes, direction handling, configuration, and plotting.
4. [`../WASP_SPR_specification.md`](../WASP_SPR_specification.md) — current SPR
   regions, compatibility, prominence, energy, and decisions.
5. [`matching.md`](matching.md) — current PCSPM descriptor requirements, active
   cost, filters, dummy assignment, and return structure.

When prose conflicts with code, `src/wasp/` remains authoritative.

## Historical methodology and research material

The following documents preserve earlier analysis or manuscript iterations.
They refer to missing external validation scripts and/or historical greedy,
weighted, or threshold-based matching implementations. They must not be used as
the package API or as the current PCSPM specification:

- `REVISION_SUMMARY.md`
- `matching_implementation_vs_documentation.md`
- `methodology_partitioning_matching.md`
- `methods.md`
- `new_methods.md`
- `methods_qualisA.md`
- manuscript introduction/results/support files
- PDF/DOCX manuscript artefacts

Historical files are retained for research provenance. Their status banners
link back to this index and the current matching documentation.

## Operational status

- The current checkout contains package examples and one downstream transition
  notebook.
- `scripts/` contains no executable validation workflow.
- Stored CSV/figure outputs are generated artefacts and may use older schemas.
- Old distributions under `dist/` are not a validation source for the current
  tree.
