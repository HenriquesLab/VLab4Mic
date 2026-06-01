# Changelog

All notable changes to VLab4Mic are documented here in user-facing terms.

The project started as `supra_molecular_simulator`, was renamed to
`supramolsim`, and later became `vlab4mic`. Older entries preserve that package
history so releases can be traced back to the corresponding version bumps.

## [0.1.0] - 2026-06-01

Preprint submission release.

### Added
- Added and refreshed the publication-oriented example suite used to generate article and supplementary figure panels, including primary/secondary labelling, probe examples, localisation simulations, NPC quantification, depth-of-field examples, structural-integrity examples, multiple structures, and structure-parsing workflows.
- Added smoke tests for notebooks and core workflows so Colab/Jupyter entry points are checked automatically.
- Added clearer issue templates for the project scope.

### Changed
- Marked the package as preprint-ready by bumping the package version from `0.0.21` to `0.1.0`.
- Moved the user-facing MkDocs website out of the main Python package repository so the website can be hosted and built independently.
- Updated Colab notebooks and notebook CI so the public notebooks are tested closer to how users run them.
- Improved package metadata, supported Python version claims, modality tables, and API examples.

### Fixed
- Fixed the `download-structures` command so it works from an installed package context.
- Fixed structure downloads so CIF files are cached in user storage
  (`~/.vlab4mic/structures`, or `VLAB4MIC_STRUCTURE_DIR`) instead of the
  installed package directory.
- Fixed CI and packaging issues, including the missing `IPython` runtime dependency.
- Fixed runtime correctness bugs in workflows and sweep generation.
- Fixed Colab notebook problems around Google Drive setup, broken links, and notebook framing.
- Fixed structural-integrity upper-bound behaviour in the publication examples.

## [0.0.21] - 2026-04-29

Publication examples and localisation-analysis release.

### Added
- Added scripts for the main figure panels, supplementary examples, and article-ready validation workflows.
- Added localisation-table simulation for SMLM-like outputs, including localisation parameter handling in modality configuration.
- Added structure/probe visualisation improvements used by the article examples.
- Added output date stamping for experiment outputs.

### Changed
- Updated SMLM modality configuration and example scripts to support localisation outputs.
- Removed outdated notebooks in favour of the current codeless notebook workflow.
- Standardised and renamed several article figure scripts for clearer use.

### Fixed
- Fixed edge cases when localisation sampling produces no emitters.
- Fixed modality parameter updates and missing modality keys.
- Fixed primary-probe degree-of-labelling handling.

## [0.0.20] - 2026-04-09

Randomness and reproducibility patch.

### Fixed
- Fixed random-seed reset behaviour when generating image simulations, improving reproducibility across repeated runs.

## [0.0.19] - 2026-04-08

Structural-integrity parameter handling release.

### Changed
- Centralised structural-integrity setup in the experiment object.
- Updated high-level experiment functions and parameter sweeps to use the same structural-integrity setter.

## [0.0.18] - 2026-04-08

Structural-integrity bugfix release.

### Fixed
- Released a patch for structural-integrity handling immediately after `0.0.17`.

## [0.0.17] - 2026-04-08

Capsid examples and labelled-structure preview release.

### Changed
- Updated capsid scripts and configuration values.
- Improved labelled-structure preview plot limits.

## [0.0.16] - 2026-04-02

Revised structural-integrity model.

### Changed
- Revised how structural integrity is modelled so the sampled structure better matches the requested integrity percentage.
- Updated labelled-particle visualisation to show only available epitopes after structural-integrity modelling.
- Reduced test runtime by avoiding redundant experiment fixture work.

### Fixed
- Fixed stopping conditions during structural-integrity iteration.
- Fixed fracture-point sampling so each iteration samples independently.

## [0.0.15] - 2026-03-30

Custom models, positioning, and reproducibility release.

### Added
- Added support for using custom model files and image-derived positions as virtual-sample inputs.
- Added structure-axis Euler rotation controls and initial-orientation handling.
- Added global normal-direction options for structure epitopes.
- Added optional random seeds to experiments, sweeps, and high-level functions.
- Added controls for explicit rotations, orientations, offsets, and per-particle placement from high-level functions.

### Fixed
- Fixed handling of missing assembly-operation fields in CIF files.
- Fixed structure path reset behaviour when switching away from local files.
- Fixed probe visualisation with degree of labelling.

## [0.0.14] - 2026-01-15

Sweep step-size and structure-format release.

### Added
- Added package badges and coverage upload support.
- Added structure-format handling for local PDB/CIF files.

### Changed
- Changed parameter-sweep ranges to use step size rather than number of values.
- Updated scripts and tests to match the new sweep range interpretation.

### Fixed
- Fixed modality parameter updates in parameter sweeps.
- Fixed parsing fallbacks when PDB/CIF files do not contain expected protein-name metadata.

## [0.0.13] - 2026-01-12

PyPI installation and sweep usability release.

### Changed
- Updated notebooks and installation instructions to install from PyPI.
- Switched progress imports away from notebook-specific `tqdm`.
- Simplified and reduced redundant parameter-sweep tests.

### Fixed
- Fixed parameter-sweep figure saving and feedback messages.
- Fixed warning-prone plotting paths, subplot indexing, and file naming.

## [0.0.12] - 2026-01-07

Parameter-sweep channel patch.

### Added
- Added a high-level parameter-sweep test.

### Fixed
- Fixed channel-index handling in parameter sweeps.

## [0.0.11] - 2026-01-06

Expanded high-level simulation controls.

### Added
- Added expansion-factor support to high-level functions and experiment methods.
- Added an AiryScan modality configuration.

### Changed
- Updated multicolour and Figure 3 scripts to use the new controls.

### Fixed
- Fixed modality handling in `image_vsample`.

## [0.0.10] - 2025-12-12

Multichannel fluorophore release.

### Added
- Added support for two fluorophores/channels in image simulation.
- Added article example scripts for primary/secondary probes and NPC quantification.

### Fixed
- Fixed changing fluorophores in probe configurations.

## [0.0.9] - 2025-12-10

Jupyter separation release.

### Changed
- Removed Jupyter widgets from the core `vlab4mic` package so widget code can live in the companion `vlab4micjupyter` package.
- Updated notebooks to import widgets from `vlab4micjupyter`.
- Removed JupyterLab/widget dependencies from the core package dependency list.
- Updated installation instructions for the split package workflow.

### Added
- Added an option to save simulated images locally from high-level functions.

### Fixed
- Fixed notebook imports after moving widget code out of the core package.

## [0.0.8] - 2025-11-17

Examples and beginner-facing documentation release.

### Added
- Added example outputs and references for real-vs-simulation scripts.
- Added a tutorial GIF and clearer README presentation.

### Changed
- Improved script README tables and descriptions.
- Improved widget feedback text and local-structure loading descriptions.
- Removed duplicated structure loading in the main experiment flow.

## [0.0.7] - 2025-11-11

Notebook naming and cleanup release.

### Changed
- Renamed the main notebook to `VLab4Mic_main.ipynb` and fixed links to it.
- Cleaned up completed TODO comments and structure-related imports.

## [0.0.6] - 2025-11-11

Repository documentation release.

### Added
- Added Code of Conduct and initial changelog files.
- Added Google Drive mounting guidance for Colab notebooks.

### Changed
- Reworked the manual and README around codeless notebooks, Colab troubleshooting, and Python-package usage.
- Updated repository links and notebook instructions after the VLab4Mic rename.

## [0.0.5] - 2025-10-22

Colab and CIF parsing patch.

### Added
- Added Colab-facing titles to the main notebook.

### Fixed
- Fixed CIF database reading when operation lists are missing.

## [0.0.4] - 2025-10-17

Sweep widget usability release.

### Added
- Added HTML display of selected parameters in sweep widgets.
- Added methods to update sweep sliders programmatically.

### Fixed
- Fixed sweep widget slider values and manual typos.

## [0.0.3] - 2025-10-16

Reference-image workflow release.

### Added
- Added optional reference masks for sweep analysis.
- Added support for TIFF/TIF reference images.
- Added widgets for uploading reference images and masks.
- Added notebook guidance for local structure selection and advanced reference parameters.

### Fixed
- Fixed reference-sample handling when using local structure files.
- Fixed image generation and acquisition-parameter behaviour in the experiment widget.

## [0.0.2] - 2025-10-09

First public example scripts release.

### Added
- Added installation/testing instructions for TestPyPI.
- Added example scripts for imaging a virtual sample and running a parameter sweep.
- Moved example scripts into their current folder layout.

### Fixed
- Fixed mask-index handling.

## [0.0.1] - 2025-10-09

Alpha package-testing release.

### Changed
- Reset the package version from the earlier `0.1.0` baseline to `0.0.1` for alpha-stage TestPyPI deployment.

## [vlab4mic 0.1.0 baseline] - 2025-08-14

Renaming release before the alpha version reset.

### Changed
- Renamed the package from `supramolsim` to `vlab4mic`.
- Updated package metadata, coverage settings, and console entry points to use the VLab4Mic name.
- Added particle-defect controls to notebooks and improved notebook markdown.

### Fixed
- Fixed a reported issue around the renamed package workflow.

## [supramolsim 0.1.0] - 2025-05-16

Setuptools packaging baseline.

### Changed
- Deprecated the old `setup.cfg` workflow and moved package metadata into `pyproject.toml`.
- Declared setuptools package discovery, package data, dependencies, optional dependency groups, and the structure-download command.

## [2024.10.18] - 2024-10-18

### Added
- Added test fixtures for structure and coordinate-field workflows.

## [2024.10.15] - 2024-10-15

### Changed
- Moved the project from Poetry to setuptools.
- Added GitHub Actions, Makefile commands, pdoc documentation generation, and Python packaging commands.

### Fixed
- Fixed the test setup after the packaging migration.

## [2024.10.09] - 2024-10-09

### Fixed
- Fixed label display in the Jupyter widget workflow.
- Updated generated CIF ignore rules.

## [2024.10.08] - 2024-10-08

### Added
- Added image-noise transforms, image convolution tools, text/TIFF writing utilities, and multi-modality imaging workflows.
- Added widget controls for acquisition parameters and running simulations.
- Added a demo notebook and local configuration retrieval.

### Fixed
- Fixed imaging imports and multi-modality image-generation tests.

## [2024.10.07] - 2024-10-07

### Fixed
- Fixed data-visualisation imports.

## [2024.10.04] - 2024-10-04

### Added
- Added PSF generation and initial image-generation modules.
- Added workflows for creating an imaging system from field data.

## [2024.10.03] - 2024-10-03

### Added
- Added widget-driven structural-model creation.
- Added label construction from configuration files.
- Added labelled-particle creation, particle defects, coordinate-field generation, and associated tests.
- Added `scikit-learn` and test coverage support.

## [2024.10.02] - 2024-10-02

### Changed
- Updated project files and custom widget imports during the refactor.

## [2024.10.01] - 2024-10-01

### Changed
- Continued the large refactor into `generate`, `utils`, visualisation, workflow, and IO modules.
- Updated formatting and linting configuration.

## [2024.09.27] - 2024-09-27

### Added
- Added initial Jupyter widget module and `ipykernel` dependency.

### Changed
- Renamed the package from `supra_molecular_simulator` to `supramolsim`.

## [2024.06.14] - 2024-06-14

### Added
- Added templates for structures, labels, fluorophores, detectors, modalities, and PSFs.
- Added configuration metadata for the early simulator workflow.

### Changed
- Ignored generated coverage files.

## [2024.06.13] - 2024-06-13

### Changed
- Moved the project back into a `src/` package layout.
- Updated pre-commit and test tooling.
- Declared the package path in `pyproject.toml`.

## [2024.06.12] - 2024-06-12

### Changed
- Renamed the structure-download console command to `download-structures`.
- Performed early cleanup around packaging and scripts.

## [2024.06.07] - 2024-06-07

### Added
- Added the first structure downloader command.
- Added example CIF data and additional development dependencies.

## [2024.06.06] - 2024-06-06

Initial package skeleton.

### Added
- Added the first `supra_molecular_simulator` package metadata.
- Added initial dependencies, package declaration, and development/test tooling.
