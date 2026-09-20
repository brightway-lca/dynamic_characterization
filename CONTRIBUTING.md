# Contributing

We welcome contributions! If you have any questions, [open an issue](https://github.com/brightway-lca/dynamic_characterization/issues) or [get in touch directly with the developers](mailto:timo.diepers@ltt.rwth-aachen.de).

This project is open source under the [BSD 3-Clause License](https://github.com/brightway-lca/dynamic_characterization/blob/main/LICENSE), and everyone taking part in it is expected to follow the [Code of Conduct](https://github.com/brightway-lca/dynamic_characterization/blob/main/CODE_OF_CONDUCT.md).

## Report bugs or errors

Something is not working as expected? You have two options:

### 🥈 Report an error

Please open a new issue in the [repository](https://github.com/brightway-lca/dynamic_characterization/issues), describing the error and where you found it. It helps to say which Python version and which version of this package you are on, what you did, what you expected, and what you got instead - ideally with a small example that reproduces it.

A member of the developer community will then take care of the issue, but it may take some time for your issue to be resolved.

### 🥇 Fix an error yourself

If you have a solution to the error, you can [create a fork](https://github.com/brightway-lca/dynamic_characterization/forks) of the repository, make your changes and [create a pull request](https://github.com/brightway-lca/dynamic_characterization/pulls). The developers will assess the changes and be eternally grateful!

## Contributing to the code, examples or documentation

If you want to contribute a new characterization function or feature, share an example, or add to the documentation, please follow the [GitHub contribution workflow (fork, branch, PR)](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests) to share your work.

It is recommended to open an issue before starting work on anything, so you can talk the approach over with the maintainers first.

### Development environment

```bash
uv sync --extra testing --extra dev
```

Or, with pip in an activated environment:

```bash
pip install -e ".[testing,dev]"
```

The optional `fair` extra pulls in the [FaIR](https://github.com/OMS-/FAIR) climate model, which only the FAIR code path needs.

### Tests

```bash
uv run pytest
```

The suite lives in `tests/` and uses [pytest](https://docs.pytest.org). CI runs it on Linux, macOS and Windows, so avoid anything platform-specific. A pull request is expected to keep the suite green and to come with tests for whatever it adds.

Code style is enforced by [pre-commit](https://pre-commit.com) (black, isort, and the standard hygiene hooks). Install it once as a git hook and it runs on every commit:

```bash
uv run pre-commit install
```

### Docstring conventions

The API reference is generated from the docstrings with [mkdocstrings](https://mkdocstrings.github.io/), using the numpydoc style. Two sections need a specific format to render properly:

- **`Examples`** (plural - `Example` is not a section name mkdocstrings knows) with the code in a fenced block, not as `>>>` doctest lines:

  ````markdown
  Examples
  --------
  ```python
  df = characterize(dynamic_inventory_df, metric="radiative_forcing", time_horizon=100)
  ```
  ````

- **`See Also`** as a markdown list, where each entry links to its target. Objects documented in this package use the [mkdocstrings cross-reference syntax](https://mkdocstrings.github.io/usage/#cross-references) (`[text][identifier]`, with the identifier being the full dotted path), everything else uses a normal markdown link:

  ```markdown
  See Also
  --------
  - [`characterize`][dynamic_characterization.dynamic_characterization.characterize]: Applies the characterization functions to a dynamic inventory.
  - [`bw_timex`](https://docs.brightway.dev/projects/bw-timex/en/latest/): Package producing the time-explicit inventories this one characterizes.
  ```

Both sections end up as admonitions in the rendered docs. Build the docs and check that no `griffe:` warnings are reported - those point at docstrings that don't parse the way they look like they should.

## Building the documentation

The docs are built with [Zensical](https://zensical.org). The pages under `docs/content/examples/` are generated from `notebooks/`, so the conversion script runs first:

```bash
uv run --with-requirements docs/requirements.txt python docs/convert_notebooks.py
uv run --with-requirements docs/requirements.txt zensical serve
```

`zensical serve` rebuilds on save and serves at [http://localhost:8000](http://localhost:8000); `zensical build` writes the static site to `site/` instead. Both report broken internal links and anchors, and a build that reports issues will not render correctly on Read the Docs either.

Never edit a page under `docs/content/examples/` - it is generated and gitignored. Edit the notebook it comes from. To publish a new notebook, add it to `NOTEBOOK_META` in `docs/convert_notebooks.py` (with an icon and tags) and to the `Examples` section of the nav in `zensical.toml`.

## Releasing a new version

The package is published in three places, and each one is triggered differently:

| Target | Triggered by | Automated |
|---|---|---|
| [PyPI](https://pypi.org/project/dynamic-characterization/) | pushing the release tag | yes, by `python-package-deploy.yml` |
| [GitHub release](https://github.com/brightway-lca/dynamic_characterization/releases) | `gh release create` | no |
| [conda `diepers` channel](https://anaconda.org/diepers/dynamic_characterization) | building the recipe locally | no |

The version number lives in exactly one place, `dynamic_characterization/__init__.py`. `pyproject.toml` reads it through `tool.setuptools.dynamic`, so bumping that one line is enough.

### 1. Prepare the release PR

Branch off `main`, never commit the bump directly:

```bash
git checkout -b release/vX.Y.Z origin/main
```

Bump `__version__` in `dynamic_characterization/__init__.py`, following [semantic versioning](https://semver.org) - new features mean a minor bump, fixes alone a patch.

Then close out the changelog. `CHANGES.md` collects entries under `## [Unreleased]` as PRs merge; turn that heading into the new version with today's date:

```markdown
## [X.Y.Z] - (YYYY-MM-DD)
```

Before opening the PR, check the section is actually complete. Entries are easy to lose when a change lands outside the usual PR flow - anything pushed straight to `main`, or shipped as a side effect of a larger branch. Compare against the commits since the last tag:

```bash
git log --oneline vLAST..origin/main
```

Every user-visible change needs a line. Open the PR, let CI pass, merge it.

### 2. Tag and create the GitHub release

Tag the **merge commit** of that PR, not your local branch tip:

```bash
git checkout main && git pull --ff-only
git tag -a vX.Y.Z -m "vX.Y.Z"
git push origin vX.Y.Z
```

!!! warning "This is the point of no return"

    Pushing the tag is the irreversible step. `python-package-deploy.yml` publishes to PyPI on any tagged ref, and a version number can never be reused there - not even after deleting the release. Make sure the tag points where you think it does before pushing it.

Release notes are the changelog section for that version, copied verbatim - no rewriting, no summarising:

```bash
awk '/^## \[X\.Y\.Z\]/{f=1;next} /^## \[/{f=0} f' CHANGES.md | sed '/^$/d' > /tmp/notes.md
gh release create vX.Y.Z --title "vX.Y.Z" --notes-file /tmp/notes.md --verify-tag
```

Watch the run and confirm the version actually landed:

```bash
gh run list --limit 3
curl -s https://pypi.org/pypi/dynamic-characterization/json | python -c "import sys,json; print(json.load(sys.stdin)['info']['version'])"
```

### 3. Publish to conda

The conda package is built and uploaded by hand from a `conda/` recipe that is **not** in the repository (it is gitignored), against a checkout that is already at the released version - the recipe builds the working tree, so make sure `main` is checked out, clean, and pulled.

One-time setup:

```bash
conda install -n base conda-build anaconda-client
anaconda login          # needs write access to the `diepers` channel
```

Then build and upload the recipe, and confirm:

```bash
conda search -c diepers dynamic_characterization
```

### 4. Check the docs rebuilt

[Read the Docs](https://dynamic-characterization.readthedocs.io) builds from `.readthedocs.yml` on every push to `main`, so the release PR already triggered it. Confirm the new version renders and that the changelog page shows the new section.
