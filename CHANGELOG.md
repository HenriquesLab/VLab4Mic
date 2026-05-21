# Changelog

All notable changes to VLab4Mic will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.0.21] - 2026-04-29
### Added
- SMLM support: method to simulate localisation table of emitters from ground truth positions
- Include z-position when simulating localisations
- Option to set localisation parameters from `update_modality` method
- Emitter parameters added to modality config file and imager
- Save current date in experiment output

### Fixed
- Handle case where localisation sampling yields no emitters
- Fix order of parameter check in `update_modality`
- Bugfix: structural integrity output value not initialised
- Fix missing key in modality params

### Changed
- Removed outdated example notebooks

## [0.0.20] - 2026-04-09
### Fixed
- Bugfix: resetting random seed when generating image simulations

## [0.0.19] - 2026-04-08
### Added
- Experiment: dedicated method to set structural integrity parameters
- Sweeps now use the dedicated structural integrity setup method

## [0.0.18] - 2026-04-08
### Fixed
- Bugfix in structural integrity computation

## [0.0.17] - 2026-04-08
### Fixed
- Update plot limits for labelled structure preview

## [0.0.16] - 2026-04-02
### Fixed
- Fix structural integrity condition to stop iterations correctly
- Sample fracture point each iteration in labelled particle
- Show only available epitopes after modelling structural integrity
- Ensure structural integrity percentage matches input parameter
- Fix creation of protein names dictionary for missing fields in CIF formats

## [0.0.15] - 2026-03-30
### Added
- Custom metrics support for parameter sweep (pass list of callables)
- Use image for positioning epitopes, including elevation mask
- Virtual sample: explicit option for randomised parameters
- Random seed option in high level experiment functions and sweep generator
- Support Euler rotations and initial orientation for labelled structures in virtual sample
- Functions to assign a global normal direction for structure epitopes
- Track and return number of labelled epitopes
- Show central axis when displaying assembly atoms

### Changed
- Renamed "incomplete_labelling" / "defects" to **structural integrity** throughout codebase, configs, tests and notebooks
- Renamed conjugation efficiency to **Degree of Labelling (DoL)**; DoL now used to sample conjugation sites
- Dropped support for Python 3.9; minimum numpy version set to 2.0.0

### Fixed
- Fix conditions when looking for assembly operations in creating structure
- Handle exceptions for missing fields in CIF files
- Fix probe parameterisation when specifying PDB model
- Bugfix: set structure path to None when not using local file
- Fix parameter names in experiment for single probe

## [0.0.14] - 2026-01-15
### Added
- Generate parameter iterables with stepsize instead of total number of values
- Add `structure_format` parameter; set automatically from file extension if provided
- Handle PDB format when creating CIF dictionary

### Fixed
- Return empty list if PDB/CIF does not contain protein names
- Fixed missing modality update parameters not used in parameter sweep

## [0.0.13] - 2026-01-12
### Added
- Feedback messages when running sweep
- ParameterSweep: option to treat NaN values as zero when plotting (does not affect dataframe output)
- ParameterSweep: fix subplot creation with `squeeze=False` to avoid index mismatches
- Add increasing number suffix on output filenames

### Fixed
- Fix bug when saving figures from sweep
- Resolve warning from `np.diff` call without explicit index

## [0.0.12] - 2026-01-07
### Fixed
- Fix channel index in parameter sweep

## [0.0.11] - 2026-01-06
### Added
- Expansion factor support in high level functions and experiment methods
- Add missing AiryScan modality
- Experiment: update modality and create new from template parameter
- Structure: retrieve position by residue ID
- Support input list of probes for multiple labels
- Import fluorophore configs at initialisation; method to add fluorophores to experiment
- Example scripts for multicolour, site-specific labelling and AlphaFold models

### Fixed
- Fix sequential labelling not clearing default probe
- Bugfix: changing fluorophore in probes

## [0.0.10]
- Use 2 different fluorophores/channels for image simulation

## [0.0.9]
- Removed module for jupyter widgets