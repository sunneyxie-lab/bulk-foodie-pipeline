# CHANGELOG

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project adheres to [Semantic Versioning](https://semver.org/).

## [1.2.0] - 2025-12-10

### Added
- Two adjustment modes (pc, nonpc) are supported.
- Dynamically determine the optimal threshold for the conversion ratio.

## [1.1.0] - 2025-11-11

### Fixed
- Upgrade from macs2 to macs3 to fix the [ImportError when use macs2](https://github.com/macs3-project/MACS/issues/703).

## [1.1.0] - 2025-09-25

### Added

- Add TSS (Transcription Start Site) QC plot generation to assess signal around TSS regions.
- Parallelize footprint calling (fixed to 10 cores) to reduce runtime on multi-core systems and improve throughput on large datasets.

### Fixed

- Pin `igvtools` to version `2.16.2` (downgraded from `2.17.3`) to resolve TDF-generation issues observed with `2.17.3`.
