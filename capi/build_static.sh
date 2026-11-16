#!/usr/bin/env bash
# build/libplegma.a: plegma's C API with the Mojo runtime functions it calls (capi/plegma_static_runtime.c), so a
# tool links plegma statically, into a static executable or (position independent) a shared library, with no Mojo
# runtime library at run time. Run through pixi: pixi run build-static.
#
# Mojo compiles for the build machine's CPU unless told otherwise. On x86-64 the library is built for x86-64-v3
# (AVX2, Intel Haswell and AMD Excavator on), so it runs on other machines; PLEGMA_TARGET_CPU chooses another
# (native: this machine's). x86-64-v2 does not build with the pinned nightly.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
target=${PLEGMA_TARGET_CPU:-}
if [ -z "$target" ] && [ "$(uname -m)" = x86_64 ]; then
  target=x86-64-v3
fi
cpu=()
if [ -n "$target" ] && [ "$target" != native ]; then
  cpu=(--target-cpu "$target")
fi
mkdir -p build
mojo build --emit object "${cpu[@]}" -D ASSERT=none -I . -I deps/sporos -I deps/karan capi/plegma_capi.mojo \
  -o build/plegma_capi.o
"${CC:-cc}" -O2 -fPIC -Wall -Wextra -c capi/plegma_static_runtime.c -o build/plegma_static_runtime.o
rm -f build/libplegma.a
ar rcs build/libplegma.a build/plegma_capi.o build/plegma_static_runtime.o
echo "build_static: build/libplegma.a for CPU ${target:-native}"
