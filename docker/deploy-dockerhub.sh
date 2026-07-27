#!/bin/bash
# Build the RiboPilot Docker image and push it to Docker Hub.
# Usage: docker/deploy-dockerhub.sh [tag]   (run from the repository root)
#
# The tag defaults to the DESCRIPTION Version, which is what new_project()
# derives its default `container` from. Keep them derived from the same source
# or the published tag and the tag every project asks for will drift apart.

repo='ribopilot'
tag="${1:-$(awk '/^Version:/{print $2}' DESCRIPTION)}"

if [ -z "$tag" ]; then
  echo "Could not determine a tag (no argument, and no Version in DESCRIPTION)."
  echo "Run from the repository root, or pass a tag explicitly."
  exit 1
fi
echo "Image tag: ${repo}:${tag}"

read -p "Enter your Docker Hub username: " username
read -s -p "Enter your Docker Hub password: " password
echo

echo "Logging into Docker Hub..."
echo "$password" | docker login -u "$username" --password-stdin
if [ $? -ne 0 ]; then
  echo "Failed to log in to Docker Hub"
  exit 1
fi

# Build the R package tarball into docker/ for the Dockerfile to COPY.
# NOTE: no devtools::document() here -- RiboPilot's NAMESPACE is hand-maintained.
Rscript -e 'devtools::build(path = "docker", vignettes = FALSE)'
if [ $? -ne 0 ]; then
  echo "Failed to build the R package tarball"
  exit 1
fi

docker build -f docker/Dockerfile --progress=plain -t ${repo}:${tag} .
if [ $? -ne 0 ]; then
  echo "Failed to build the Docker image"
  exit 1
fi

docker tag ${repo}:${tag} ${username}/${repo}:${tag}
if [ $? -ne 0 ]; then
  echo "Failed to tag the Docker image"
  exit 1
fi

docker push ${username}/${repo}:${tag}
if [ $? -ne 0 ]; then
  echo "Failed to push the Docker image to Docker Hub"
  exit 1
fi

echo "Image pushed successfully: ${username}/${repo}:${tag}"
