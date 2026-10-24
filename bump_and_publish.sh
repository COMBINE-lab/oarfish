#!/usr/bin/env bash
# Cut an oarfish release: bump the version, date the CHANGELOG's Unreleased section, check that
# the crate builds as it will be published, and push a `[do_tag]` release commit to main. CI
# does the rest, and this script follows it:
#
#   Sanitize Cargo (.github/workflows/sanitize-cargo.yml): on a push to main whose commit
#     message has [do_tag], replaces the git dependencies with their crates.io versions,
#     commits that (off main) and tags it v<version>;
#   Release (release.yml, cargo-dist): on the tag, builds the binaries and the GitHub release,
#     with the CHANGELOG section as its notes;
#   Publish to crates.io (publish-crates-io.yml): after a successful Release, `cargo publish`.
#
# Run from a clean, up-to-date main.
set -euo pipefail

die() {
    echo "error: $*" >&2
    exit 1
}

usage() {
    cat <<'EOF'
Usage:
  ./bump_and_publish.sh <version> [--dry-run] [--no-watch] [--allow-same-version]
  ./bump_and_publish.sh <version> --watch-only

Sets the package version to <version>, renames the CHANGELOG's `## Unreleased` heading to
`## <version> - <today>` (with a new, empty `## Unreleased` above it), runs the tests, and
packages the crate with its git dependencies replaced by their crates.io versions, as the
release will be (`cargo-sanitize`, then `cargo publish --dry-run`). It then commits
`release: v<version> [do_tag]`, pushes it to main, and follows the CI runs that tag, build and
publish the release.

Options:
  --dry-run             Run the checks and show what would be done, without changing files,
                        committing or pushing.
  --no-watch            Push the release commit and stop, without following CI.
  --watch-only          Follow the CI runs of a release commit that was already pushed.
  --allow-same-version  Release the version the package is already at (e.g. to push the
                        [do_tag] commit again after the tagging run failed before tagging).
  -h, --help            Show this help message.

Examples:
  ./bump_and_publish.sh 0.11.0 --dry-run
  ./bump_and_publish.sh 0.11.0
  ./bump_and_publish.sh 0.11.0 --watch-only
EOF
}

print_cmd() {
    printf '+'
    printf ' %q' "$@"
    printf '\n'
}

run() {
    print_cmd "$@"
    if [[ "$DRY_RUN" == true ]]; then
        return 0
    fi
    "$@"
}

VERSION=""
DRY_RUN=false
WATCH=true
WATCH_ONLY=false
ALLOW_SAME_VERSION=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run) DRY_RUN=true ;;
        --no-watch) WATCH=false ;;
        --watch-only) WATCH_ONLY=true ;;
        --allow-same-version) ALLOW_SAME_VERSION=true ;;
        -h | --help)
            usage
            exit 0
            ;;
        -*) die "unknown option: $1" ;;
        *)
            [[ -z "$VERSION" ]] || die "version specified more than once"
            VERSION="$1"
            ;;
    esac
    shift
done

[[ -n "$VERSION" ]] || {
    usage
    exit 1
}

if ! [[ "$VERSION" =~ ^[0-9]+\.[0-9]+\.[0-9]+([+-][0-9A-Za-z.-]+)*$ ]]; then
    die "version must look like X.Y.Z, optionally with prerelease/build suffixes"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

CARGO_TOML="Cargo.toml"
CHANGELOG="CHANGELOG.md"
TAG="v${VERSION}"
REPO="COMBINE-lab/oarfish"
REPO_URL="https://github.com/${REPO}"
SUBJECT="release: v${VERSION} [do_tag]"
TODAY="$(date +%Y-%m-%d)"

for f in "$CARGO_TOML" "$CHANGELOG"; do
    [[ -f "$f" ]] || die "not found: $f"
done

on_crates_io() {
    local code
    code="$(curl -s -o /dev/null -w '%{http_code}' -A "oarfish-release (${REPO_URL})" \
        "https://crates.io/api/v1/crates/oarfish/${VERSION}")"
    [[ "$code" == "200" ]]
}

# The commit the tag points to (the sanitized commit, which is not on main), or nothing.
remote_tag_commit() {
    git ls-remote origin "refs/tags/${TAG}^{}" | cut -f1
}

# The id of the newest run of workflow $1 matching the `gh run list` filters that follow,
# waiting up to $WAIT_RUN seconds for it to appear.
WAIT_RUN=300
wait_for_run() {
    local wf="$1" id waited=0
    shift
    while true; do
        id="$(gh run list -R "$REPO" --workflow "$wf" "$@" --limit 1 \
            --json databaseId --jq '.[0].databaseId // empty')"
        [[ -n "$id" ]] && break
        [[ "$waited" -lt "$WAIT_RUN" ]] || return 1
        sleep 10
        waited=$((waited + 10))
    done
    echo "$id"
}

# Follows the CI runs of the release commit $1, then prints what the bioconda recipe needs.
watch_release() {
    local sha="$1" id tagged waited
    command -v gh >/dev/null || die "following the release needs the GitHub CLI (gh)"

    echo
    echo "Waiting for the Sanitize Cargo run on ${sha:0:12} (tags ${TAG})"
    id="$(wait_for_run sanitize-cargo.yml --commit "$sha")" ||
        die "no Sanitize Cargo run appeared for ${sha:0:12}; check ${REPO_URL}/actions"
    gh run watch -R "$REPO" "$id" --exit-status --compact ||
        die "the Sanitize Cargo run failed (${REPO_URL}/actions/runs/${id}); nothing was tagged.
       Fix it, then re-run with --allow-same-version to push the [do_tag] commit again"
    tagged="$(remote_tag_commit)"
    [[ -n "$tagged" ]] || die "the Sanitize Cargo run succeeded but ${TAG} is not on origin"

    echo
    echo "Waiting for the Release run on ${TAG} (${tagged:0:12})"
    id="$(wait_for_run release.yml --commit "$tagged" --event push)" ||
        die "no Release run appeared for ${TAG}; check ${REPO_URL}/actions/workflows/release.yml"
    gh run watch -R "$REPO" "$id" --exit-status --compact ||
        die "the Release run failed (${REPO_URL}/actions/runs/${id}), so nothing was published to
       crates.io. Once it is fixed and re-run, publish with:
         gh workflow run publish-crates-io.yml -R ${REPO} -f tag=${TAG}"

    # The publish run is triggered by the Release run's completion (workflow_run), so it carries
    # main's head commit rather than the tag's; crates.io itself says when it is done.
    echo
    echo "Waiting for oarfish ${VERSION} on crates.io (Publish to crates.io)"
    waited=0
    until on_crates_io; do
        [[ "$waited" -lt 1200 ]] || die "oarfish ${VERSION} is not on crates.io after 20 minutes;
       see ${REPO_URL}/actions/workflows/publish-crates-io.yml, or publish by hand with:
         gh workflow run publish-crates-io.yml -R ${REPO} -f tag=${TAG}"
        sleep 20
        waited=$((waited + 20))
    done

    echo
    echo "Released oarfish ${VERSION}: ${REPO_URL}/releases/tag/${TAG}"
    local tarball="${REPO_URL}/archive/refs/tags/${TAG}.tar.gz" sum sha256
    if command -v sha256sum >/dev/null; then sum=(sha256sum); else sum=(shasum -a 256); fi
    if sha256="$(curl -sfL "$tarball" | "${sum[@]}" | cut -d' ' -f1)" && [[ -n "$sha256" ]]; then
        echo
        echo "For the bioconda recipe:"
        echo "  url: $tarball"
        echo "  sha256: $sha256"
    fi
}

if [[ "$WATCH_ONLY" == true ]]; then
    git fetch -q origin main
    # The release commit: `release: v<version> ... [do_tag]` (earlier releases put a summary
    # between the two).
    SHA="$(git log origin/main --format='%H%x09%s' | awk -F'\t' -v p="release: v${VERSION}" '
        !found && (index($2, p " ") == 1 || $2 == p) && index($2, "[do_tag]") { print $1; found = 1 }')"
    [[ -n "$SHA" ]] || die "no '${SUBJECT}' commit on origin/main"
    watch_release "$SHA"
    exit 0
fi

# ── Pre-flight: a clean main, in sync with origin, and the release not yet out ────────────────

[[ -z "$(git status --porcelain)" ]] || die "working tree is not clean; commit or stash existing changes first"

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
[[ "$BRANCH" == "main" ]] || die "releases are cut from main (on: $BRANCH)"
git fetch -q origin main
RESUMING=false
if [[ "$(git rev-parse HEAD)" != "$(git rev-parse origin/main)" ]]; then
    # A release whose push failed leaves its release commit unpushed: main is then exactly that
    # commit ahead of origin/main, and the release resumes from it.
    if git merge-base --is-ancestor origin/main HEAD &&
        [[ "$(git rev-list --count origin/main..HEAD)" == 1 ]] &&
        [[ "$(git log -1 --format=%s)" == "$SUBJECT" ]]; then
        RESUMING=true
        echo "Resuming the release of v${VERSION} from its unpushed release commit"
    else
        die "main is not at origin/main; pull (or push) first"
    fi
fi

if git rev-parse -q --verify "refs/tags/$TAG" >/dev/null || [[ -n "$(git ls-remote --tags origin "refs/tags/$TAG")" ]]; then
    die "tag $TAG already exists (to follow its CI runs: ./bump_and_publish.sh $VERSION --watch-only)"
fi
! on_crates_io || die "oarfish $VERSION is already on crates.io"

command -v cargo-sanitize >/dev/null ||
    die "cargo-sanitize is needed to check the release manifest: cargo install cargo-sanitize"
if [[ "$WATCH" == true ]]; then
    command -v gh >/dev/null || die "following the release needs the GitHub CLI (gh); or pass --no-watch"
fi

# ── Version ─────────────────────────────────────────────────────────────────────────────────

package_version() {
    awk -F'"' '/^\[package\]/{p=1; next} /^\[/{p=0} p && /^version[[:space:]]*=/{print $2; exit}' "$CARGO_TOML"
}

CURRENT_VERSION="$(package_version)"
[[ -n "$CURRENT_VERSION" ]] || die "could not determine the [package] version from $CARGO_TOML"

SAME_VERSION=false
if [[ "$CURRENT_VERSION" == "$VERSION" ]]; then
    if [[ "$RESUMING" == true || "$ALLOW_SAME_VERSION" == true ]]; then
        SAME_VERSION=true
    else
        die "the package version is already $VERSION (pass --allow-same-version to release it as-is)"
    fi
elif [[ "$(printf '%s\n%s\n' "$CURRENT_VERSION" "$VERSION" | sort -V | tail -1)" != "$VERSION" ]]; then
    die "$VERSION is lower than the current version $CURRENT_VERSION"
fi

# ── CHANGELOG: dist takes the release notes from the section headed `## <version> - <date>` ───

VERSION_RE="${VERSION//./\\.}"
DATED_HEADING="$(grep -m1 -E "^## ${VERSION_RE} - [0-9]{4}-[0-9]{2}-[0-9]{2}[[:space:]]*$" "$CHANGELOG" || true)"
if [[ -n "$DATED_HEADING" ]]; then
    NEW_HEADING="$DATED_HEADING"
    UNRELEASED=false
else
    grep -q -E '^## Unreleased[[:space:]]*$' "$CHANGELOG" ||
        die "$CHANGELOG has no '## Unreleased' section (or '## ${VERSION} - <date>') with the release notes"
    # The Unreleased section must say something: dist uses it as the release body.
    [[ -n "$(awk '/^## Unreleased[[:space:]]*$/{p=1; next} p && /^## /{exit} p && NF' "$CHANGELOG")" ]] ||
        die "the '## Unreleased' section of $CHANGELOG is empty"
    NEW_HEADING="## ${VERSION} - ${TODAY}"
    UNRELEASED=true
fi

echo "Current version   : $CURRENT_VERSION"
echo "New version       : $VERSION"
echo "CHANGELOG section : $NEW_HEADING"
echo "Release commit    : $SUBJECT"
echo "Tag (made by CI)  : $TAG"
echo "Follow CI         : $([[ "$WATCH" == true ]] && echo yes || echo no)"
echo "Dry-run           : $([[ "$DRY_RUN" == true ]] && echo yes || echo no)"
echo

# ── Checks ──────────────────────────────────────────────────────────────────────────────────

echo "Running the tests"
cargo test --release -q

# The release publishes Cargo.toml as cargo-sanitize leaves it: git dependencies replaced by
# their crates.io versions. Package and build exactly that, so a git dependency that has moved
# past its crates.io release fails here rather than after the tag exists.
echo "Checking the sanitized package (git dependencies -> crates.io)"
WORK="$(mktemp -d)"
trap 'rm -rf "${WORK:?}"' EXIT
git archive HEAD | tar -x -C "$WORK"
if [[ "$SAME_VERSION" == false ]]; then
    awk -v ver="$VERSION" '
        /^\[package\]/ { p = 1; print; next } /^\[/ { p = 0 }
        p && /^version[[:space:]]*=/ && !done { sub(/"[^"]*"/, "\"" ver "\""); done = 1 }
        { print }
    ' "$WORK/$CARGO_TOML" >"$WORK/${CARGO_TOML}.cs_orig"
else
    cp "$WORK/$CARGO_TOML" "$WORK/${CARGO_TOML}.cs_orig"
fi
cargo-sanitize -i "$WORK/${CARGO_TOML}.cs_orig" -o "$WORK/$CARGO_TOML"
if grep -n -E '(^|[{,[:space:]])(git|branch|rev|tag)[[:space:]]*=' "$WORK/$CARGO_TOML"; then
    die "the sanitized $CARGO_TOML still has git dependencies (above)"
fi
(cd "$WORK" && CARGO_TARGET_DIR="$SCRIPT_DIR/target" cargo publish --dry-run --allow-dirty -q) ||
    die "the sanitized package does not build; a git dependency may need a crates.io release first"
echo "The sanitized package builds"
echo

# ── Edits ───────────────────────────────────────────────────────────────────────────────────

if [[ "$SAME_VERSION" == false ]]; then
    echo "Updating $CARGO_TOML: version $CURRENT_VERSION -> $VERSION"
fi
if [[ "$UNRELEASED" == true ]]; then
    echo "Updating $CHANGELOG: '## Unreleased' -> '$NEW_HEADING', with a new empty '## Unreleased' above"
fi

if [[ "$DRY_RUN" == false ]]; then
    if [[ "$SAME_VERSION" == false ]]; then
        awk -v ver="$VERSION" '
            /^\[package\]/ { p = 1; print; next } /^\[/ { p = 0 }
            p && /^version[[:space:]]*=/ && !done { sub(/"[^"]*"/, "\"" ver "\""); done = 1 }
            { print }
        ' "$CARGO_TOML" >"${CARGO_TOML}.tmp" && mv "${CARGO_TOML}.tmp" "$CARGO_TOML"
        [[ "$(package_version)" == "$VERSION" ]] || die "version update failed"
    fi
    if [[ "$UNRELEASED" == true ]]; then
        awk -v new="$NEW_HEADING" '
            !done && /^## Unreleased[[:space:]]*$/ { print "## Unreleased"; print ""; print new; done = 1; next }
            { print }
        ' "$CHANGELOG" >"${CHANGELOG}.tmp" && mv "${CHANGELOG}.tmp" "$CHANGELOG"
    fi
fi

if [[ "$RESUMING" == true ]]; then
    echo "Release commit already made"
elif [[ "$DRY_RUN" == true || -n "$(git status --porcelain)" ]]; then
    run git add "$CARGO_TOML" "$CHANGELOG"
    run git commit -m "$SUBJECT"
else
    # Already at the version, with a dated CHANGELOG (--allow-same-version): an empty commit
    # carries the [do_tag] marker.
    run git commit --allow-empty -m "$SUBJECT"
fi

# ── Push: the [do_tag] commit on main starts the release ─────────────────────────────────────

run git push origin HEAD:main

if [[ "$DRY_RUN" == true ]]; then
    echo
    echo "Dry-run complete"
elif [[ "$WATCH" == true ]]; then
    watch_release "$(git rev-parse HEAD)"
else
    echo
    echo "Pushed ${SUBJECT}. To follow CI: ./bump_and_publish.sh ${VERSION} --watch-only"
fi
