#!/bin/bash
# Bump the DESCRIPTION version, commit it, then build and publish the image.
# Usage, from anywhere:
#
#   docker/update-version.sh
#
# The bump is committed BEFORE the build on purpose: deploy-dockerhub.sh labels
# the image with the commit it came from and stamps "-dirty" on a tree with
# uncommitted changes, so an uncommitted bump publishes an image that matches no
# commit. Nothing is pushed to git; only the image goes anywhere.

set -euo pipefail

cd "$(dirname "$0")/.."
[ -f DESCRIPTION ] || { echo "No DESCRIPTION next to docker/." >&2; exit 1; }

old="$(awk '/^Version:/{print $2}' DESCRIPTION)"
[ -n "$old" ] || { echo "No Version in DESCRIPTION." >&2; exit 1; }
new="${old%.*}.$(( ${old##*.} + 1 ))"

sed -i "s/^Version: .*/Version: ${new}/" DESCRIPTION
echo "Version ${old} -> ${new}"

git add DESCRIPTION
git commit -m "chore: bump version to ${new}"

# deploy-dockerhub.sh asks y/N before it pushes. Running the whole thing IS the
# confirmation here, and there is no terminal to answer from.
printf 'y\n' | docker/deploy-dockerhub.sh
