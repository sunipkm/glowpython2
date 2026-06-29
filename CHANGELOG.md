# Changelog

All notable changes to glowpython2 are documented here.

## [Unreleased]
### Added
- `CHANGELOG.md` with history migrated from README
- GitHub Actions CI workflow for Linux, macOS, and Windows
- Functionality test suite with `pytest`

### Changed
- `CONTRIBUTING.md` expanded with changelog workflow and guidelines
- Changelog section removed from `README.md`
- README install instructions updated for Windows and general clarity
- `python-dateutil` dependency removed; minimum Python version bumped
- `pyproject.toml` dependencies tidied
- Moved files from tests/ to Examples & updates README to use local versions of some images

### Fixed
- Windows install:
  - gfortran now recommended to be installed via conda instead of MSYS2
  - pip forced to use shorter paths to avoid path-length issues

## [0.0.3] - 2026-01-09
### Changed
- Dataset structure changed (**breaking**)

### Fixed
- Version bumped to 0.0.3

## [0.0.1] - 2025-11-01
### Added
- Initial release
