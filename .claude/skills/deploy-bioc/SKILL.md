---
name: deploy-bioc
description: Deploy / release / publish this Bioconductor R package to GitHub. Runs the full pre-push ritual — document, update vignette, code review, R CMD check, BiocCheck — fixing issues, then bumps the version and commits & pushes. Use when asked to deploy, release, ship, publish, or "run the checks before pushing" the package, or to prepare a Bioconductor push. Works in any Bioconductor-track package.
---

# Deploy a Bioconductor package

"Deploying" means running the full pre-push ritual so GitHub CI stays green and
the package stays clean for Bioconductor review. This skill is **package-agnostic**:
it reads the package name and layout from the repo, so the same files work in any
Bioconductor-track package. Work through the steps **in order**; do not skip to
the commit.

Read the package name once — it names the BiocCheck output folder and is handy in
messages:

```bash
Rscript -e 'cat(read.dcf("DESCRIPTION")[, "Package"], "\n")'
```

The mechanical check phases are wrapped by the helper
[deploy.R](deploy.R):

```bash
Rscript .claude/skills/deploy-bioc/deploy.R <document|check|bioccheck|all>
```

It writes logs to `$DEPLOY_OUT` (default: a tempdir it prints on start) and exits
non-zero when a phase surfaces a **blocking** problem (check errors or warnings,
real BiocCheck errors/warnings). Set `DEPLOY_OUT` to a stable path so you can read
the logs:

```bash
export DEPLOY_OUT=/tmp/deploy-out
```

All paths below are relative to the repo root, and `deploy.R` must be run from
there (it uses `devtools`, which loads the source tree in place).

## Step 0 — Preflight

Confirm the push target with the user. Default to the current branch:

```bash
git rev-parse --abbrev-ref HEAD && git status --short
```

Most runs push the branch already checked out — **confirm it each time** before
the final push. Note what's already modified so you can describe it in the commit
later.

## Step 1 — Document

```bash
Rscript .claude/skills/deploy-bioc/deploy.R document
```

Regenerates `man/*.Rd` and `NAMESPACE` from roxygen comments. Review the diff —
`document()` **rewrites `NAMESPACE`**, so a stale `@export`/`@importFrom` shows up
here.

## Step 2 — Update the vignette (judgment step)

Only if a **significant new feature** shipped. Compare exported functions in
[NAMESPACE](NAMESPACE) and recent commits (`git log --oneline -15`) against the
package vignette(s) under `vignettes/*.Rmd`. If a notable new exported function or
capability isn't demonstrated, add a short section. Skip silently if nothing
significant changed — do not pad the vignette.

If you edit it, confirm it still knits before moving on (a broken vignette fails
`check` in step 4 anyway):

```bash
Rscript -e 'devtools::build_vignettes()'
```

## Step 3 — Code review

Invoke the **`/code-review`** skill on the working diff, triage the findings, and
apply fixes. (You run the skill yourself — it is deliberately *not* called from
`deploy.R`.) Re-run `document` if a fix touched roxygen.

## Step 4 — R CMD check

```bash
Rscript .claude/skills/deploy-bioc/deploy.R check
```

Read `$DEPLOY_OUT/check.log`. Fix every ERROR and WARNING, and NOTEs where
feasible; re-run until the summary reads `0 error(s) | 0 warning(s)`. This is slow
(it builds the package and runs the tests **and the vignette**).

## Step 5 — BiocCheck

```bash
Rscript .claude/skills/deploy-bioc/deploy.R bioccheck
```

Runs `BiocCheck` with `no-check-version-num = TRUE` (the manual version bump in
step 6 covers versioning). Fix all errors and warnings; resolve notes where
reasonable, else leave them. A clean BiocCheck is required for Bioconductor even
if the repo's CI (under `.github/workflows/`) doesn't run BiocCheck itself — so
run it as part of every deploy. See Gotchas for the two "errors" the helper marks
**environmental / non-blocking**.

## Step 6 — Version + NEWS

**Always** bump the **patch (last) component** of `Version:` in
[DESCRIPTION](DESCRIPTION). Bioconductor conventions:

- Under Bioc **review**, the version is `0.99.z` — bump `z` on every push.
- An **accepted** package on the Bioc **devel** branch is `x.y.z` with an odd `y`
  — also bump `z` on each change.

Updating `NEWS.md` (if the package has one) is a **judgment step, gated on the
kind of change**:

- **Write a NEWS entry only for major feature upgrades** — a new exported function
  or a notable new user-facing capability. Add a new `# <pkg> <version>` heading at
  the top, grouped under `## New Features` / `## Improvements`, and describe the
  change from the user's perspective.
- **Skip NEWS entirely for everything else** — documentation, attribution,
  vignette wording, pkgdown/config, refactors, and bug-fix-only pushes. Bump the
  version but leave `NEWS.md` untouched. Treat NEWS as a major-feature changelog,
  not a per-push log.

When unsure whether a change is "major," it isn't — default to skipping NEWS.

## Step 7 — Commit & push

Stage, commit, and push to the branch confirmed in step 0. End the commit message
with a `Co-Authored-By` trailer identifying the Claude model running the deploy,
e.g.:

```
Co-Authored-By: Claude <model> <noreply@anthropic.com>
```

If the target is the default branch (`main`/`master`), branch first per repo
policy unless the user explicitly asked to push straight to it.

## Gotchas

- **The `<pkg>.BiocCheck` output folder.** BiocCheck writes a `<pkg>.BiocCheck/`
  folder into the repo, then flags a leftover one on the *next* run
  (`checkBiocCheckOutputFolder` ERROR). `deploy.R` deletes it (via a
  `*.BiocCheck` glob) before and after each run so this never fires — but if you
  call `BiocCheck()` by hand, remove the folder afterwards and never commit it.
- **`checkSupportReg` "ERROR" is a network flake.** BiocCheck hits the Bioc
  support site to verify the maintainer's email; it fails with `HTTP 504` /
  "Unable to find your email" when offline or the site is slow. `deploy.R`
  classifies both this and the output-folder check as environmental and does
  **not** count them toward the gate. Don't chase them.
- **A broken vignette aborts `check` before findings exist.** The vignette is
  re-built inside `pkgbuild::build()`, *before* R CMD check runs, and that throws
  regardless of `error_on`. `deploy.R` catches it and reports
  `check: build aborted` — treat it as a blocking error and fix the vignette.
- **BiocCheck warnings block, not just errors.** For Bioconductor a WARNING is a
  release blocker; `deploy.R` counts check+BiocCheck warnings toward the gate too.
- **`document()` can change `NAMESPACE` unexpectedly.** Always diff it after step 1.

## Troubleshooting

- **`check: build aborted` → `could not find function "<verb>"` in a vignette
  chunk.** The vignette calls a function from an *imported* package (e.g. a bare
  ggplot2 `labs`/`theme`/`aes`) but the setup chunk only does `library(<pkg>)`.
  Imported functions are not re-exported, so a bare call can't resolve. Fix: add
  `library(<the package>)` to the vignette's setup chunk, or qualify the call
  (`pkg::fun(...)`).
- **`there is no package called '<pkg>'` from deploy.R.** Run it from the repo
  root; the helper uses `devtools`, which loads the source tree in place.
- **BiocCheck note count looks huge.** Each note lists every offending file:line;
  the helper's summary counts *distinct checks* via `res$getNum()`, which is the
  number that matters.
