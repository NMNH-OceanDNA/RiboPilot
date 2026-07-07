#!/bin/bash
# Build the RiboPilot Docker image and push it to Docker Hub.
# Usage: docker/deploy-dockerhub.sh <tag>   (run from the repository root)

repo='ribopilot'
tag="${1}"

if [ -z "$tag" ]; then
  echo "You forgot to specify an image tag."
  exit 1
fi

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
