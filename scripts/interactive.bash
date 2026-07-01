IMAGE=${1:-htvs-diffsbdd:latest}
docker run -it --rm \
  -v $(pwd):/workspace \
  "${IMAGE}"
