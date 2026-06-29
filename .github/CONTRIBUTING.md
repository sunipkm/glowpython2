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

Should mostly be self-explanatory but reach out if you have questions/comments/issues.

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
