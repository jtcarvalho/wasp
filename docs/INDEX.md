# WASP documentation index

The implementation in `src/wasp/` is the source of truth. Current documents
describe that implementation; historical files preserve earlier research and
must not be used to infer the package API.

## Start here

1. [README](../README.md) — concise public overview, installation, minimal use,
   supported sources, citation, and project links.
2. [Project context](../PROJECT_CONTEXT.md) — authoritative scientific and
   maintenance context, invariants, known debt, and release expectations.
3. [Architecture](../README_archi.md) — modules, dependencies, call graphs, and
   data-flow boundaries.

## Current technical reference

| Topic | Authoritative document |
|---|---|
| Partitioning, merge, properties, numbering | [Partitioning reference](partitioning.md) |
| Secondary Peak Reassessment | [SPR reference](WASP_SPR_specification.md) |
| PCSPM descriptors, cost, filtering, assignment | [PCSPM matching](matching.md) |
| Readers, units, shapes, directions, configuration | [Processing and I/O](processing_guide.md) |
| Functions, signatures, returns, edge behavior | [API inventory](api_reference.md) |
| Plotting | [API inventory: plotting](api_reference.md#plotting) and [processing guide](processing_guide.md#plotting) |
| Examples, notebooks, scripts, test status | [Examples and notebooks](examples.md) |
| Naming and units | [Naming guide](../docs_naming.md) |
| Performance and numerical equivalence | [Performance guide](../docs_performance.md) |
| AI-agent rules | [AGENTS.md](../AGENTS.md) |
| AI-assisted review workflow | [CONTRIBUTING_AI.md](../CONTRIBUTING_AI.md) |

## Configuration and metadata

- `examples/config.yaml` is a notebook example; it is not automatically applied
  by the package.
- `pyproject.toml` is authoritative for build metadata, dependencies, Python
  requirement, and distribution version.
- `src/wasp/__init__.py` supplies the runtime version and top-level facade.
- `CITATION.cff` is authoritative for citing the exact v2.0.0 release and uses
  its version-specific DOI, `10.5281/zenodo.22676907`.
- `.zenodo.json` supplies Zenodo deposit metadata for version 2.0.0. It has no
  DOI field because Zenodo assigns the archive DOI. The all-version Concept DOI
  is `10.5281/zenodo.19415744`.

## Historical and internal material

- [Methodology and matching workflow](methodology_partitioning_matching.md) and
  [revision summary](REVISION_SUMMARY.md) are historical research documents.
  Their statements describe earlier external workflows, not current PCSPM.
- `TASKS.md` is local planning/history and intentionally ignored by Git.
- Manuscript drafts and earlier duplicate methodology files removed from the
  current working tree should remain available through Git history rather than
  being cited as current documentation.

## Repository artifact classes

- **A — core documentation:** README, this index, architecture, current
  partitioning/SPR/PCSPM/I/O/API/example guides, citation metadata.
- **B — development/internal:** project context, agent/contribution rules,
  naming/performance guides, and local task history.
- **C — historical/archive:** labelled methodology/revision documents and
  manuscript history.
- **D — generated/local:** `data/`, `output/`, `dist/`, `build/`, virtual
  environments, egg-info, caches, notebook checkpoints, and `.DS_Store`.
- **E — unresolved:** duplicate example notebooks and any unlabelled analysis
  artifact whose provenance is unclear.

The current checkout has no executable validation scripts or project test
suite. Stored notebook/CSV/figure outputs are not regression evidence for the
current tree.
