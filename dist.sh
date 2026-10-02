#!/bin/bash
# Build the Oracode image from this repository.
#
# The Dockerfile overlays the pipeline sources onto the base image, so no
# separate docker-cp / docker-commit step is needed.
set -euo pipefail

proxy_args=()
for proxy_var in HTTP_PROXY HTTPS_PROXY NO_PROXY http_proxy https_proxy no_proxy; do
  proxy_value="${!proxy_var:-}"
  if [[ -n "$proxy_value" ]]; then
    proxy_args+=(--build-arg "$proxy_var")
  fi
done

DOCKER_BUILDKIT="${DOCKER_BUILDKIT:-1}" docker build \
  "${proxy_args[@]}" \
  --build-arg "TS_IMAGE=${TS_IMAGE:-jereiard/tsuite:5.18.1}" \
  -t "${ORACODE_IMAGE:-oracode}" \
  .

echo "Built ${ORACODE_IMAGE:-oracode}. Point ./oracode at this tag to run."
