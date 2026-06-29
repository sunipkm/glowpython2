# Contributing Guidelines

Thank you for your interest in the project. Everything is still a work in progress and
your assistance/feedback will inform future development!

## Opening an issue

If you find a bug, want a new feature, or anything else, please create an issue.
Remember to:

- **Be descriptive.**
  If you are reporting a bug, include instructions on how to replicate it.
  If you had an installation error, include the exact steps you were trying & details on your system.
- **Search for relevant issues.**
  Look through the open (and recently closed) issues to see if others have reported something similar.

## Testing

`glowpython2` is tested on Linux, MacOS and Windows.
Currently three Python versions are used (stable, latest, and experimental).

The tests are run as a GitHub action, where the source file & instructions are located in [.github/workflows/ci.yml](workflows/ci.yml)

### Installation for running tests

Because `glowpython2` uses meson-python as its build backend, 
**editable installs (`pip install -e .`) do not work correctly for running the test suite**.
Data files required by the Fortran backend are not placed where the package expects them.
Use a regular install instead:

```bash
pip install .[test]
```

**After modifying Fortran source**, you must reinstall for changes to take effect:

```bash
pip install .
```

### Running tests locally

From the repo root:

```bash
pytest tests/
```

You may run an individual test file if the tests are taking too long.

### Test layout

- **`tests/test_examples.py`** — automatically runs every script in `Examples/` as a subprocess smoke test. No maintenance needed; new examples are picked up automatically.
- **`tests/test_no_precip.py`**, **`test_maxwellian.py`** — structure and physics assertions for the top-level GLOW functions.
- **`tests/test_iri.py`** — assertions for both `Iri90` and `Iri2020`.
- **`tests/test_msis.py`** — assertions for both `NrlMsis00` and `NrlMsis21`.
- **`tests/conftest.py`** — shared fixtures (`sample_time`, `sample_lat`, `sample_lon`).

### Adding or modifying tests

- **New example script**: drop a `.py` file in `Examples/` — it will be picked up by `test_examples.py` automatically.
- **New assertion test**: add a `test_*.py` file in `tests/`, or add a test function to an existing file. Use the fixtures from `conftest.py` for the standard time/location inputs.
- **Physics assertions** should check that outputs are physically reasonable (positive densities, sane temperature ranges, expected species present) rather than matching exact numerical values, so tests stay valid across minor model updates.

## Changelog

All notable changes should be recorded in `CHANGELOG.md` before a release is tagged.

**During development:** add entries under the `[Unreleased]` section at the top of `CHANGELOG.md` as you work. Use the appropriate subsection:

- **Added** — new features or capabilities
- **Changed** — changes to existing behavior
- **Fixed** — bug fixes
- **Removed** — removed features or APIs
- **Deprecated** — features that will be removed in a future release

**At release time:**
1. Rename `[Unreleased]` to `[x.y.z] - YYYY-MM-DD` matching the version tag you're about to create.
2. Add a new empty `[Unreleased]` section above it.
3. Tag the release in git: `git tag vx.y.z`.

**Example entry:**

```markdown
## [Unreleased]
### Added
- Support for custom altitude grids

### Fixed
- Incorrect electron density below 80 km under certain solar conditions
```

Do not list every commit — only changes a user of the library would care about. Trivial internal refactors, test updates, and CI fixes can be omitted.
