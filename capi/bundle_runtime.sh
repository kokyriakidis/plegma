#!/usr/bin/env bash
# Copy the Mojo runtime libraries build/libplegma.so needs from plegma's pixi environment into build/, beside it:
# the library finds them by its $ORIGIN run path, so build/ alone is what a caller links against and ships (no pixi
# environment at run time). libstdc++ and libgcc_s are left to the system (the runtime needs GLIBCXX 3.4.30).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
env_lib="${CONDA_PREFIX:?run through pixi (pixi run build-capi)}/lib"
case "$(uname -s)" in
  Linux)
    deps=$(ldd build/libplegma.so | awk '$3 ~ /^\// { print $3 }') ;;
  Darwin)
    deps=$(otool -L build/libplegma.so | awk 'NR > 1 { print $1 }' | sed "s|@rpath|$env_lib|") ;;
  *) echo "bundle_runtime: unsupported platform $(uname -s)" >&2; exit 1 ;;
esac
copied=0
for lib in $deps; do
  case "$lib" in
    "$env_lib"/*) ;;
    *) continue ;;
  esac
  case "$(basename "$lib")" in
    libstdc++*|libgcc_s*) continue ;;
  esac
  cp -L "$lib" build/
  copied=$((copied + 1))
done
echo "bundle_runtime: $copied runtime libraries beside build/libplegma.so"
