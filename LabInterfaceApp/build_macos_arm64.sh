#!/usr/bin/env bash

set -euo pipefail

app_root="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$app_root"

if [[ "$(uname -s)" != "Darwin" ]]; then
    echo "This build must run on macOS." >&2
    exit 1
fi

if [[ "$(uname -m)" != "arm64" ]]; then
    echo "This build must run on an Apple Silicon runner, not $(uname -m)." >&2
    exit 1
fi

python_arch="$(python -c 'import platform; print(platform.machine())')"
if [[ "$python_arch" != "arm64" ]]; then
    echo "Python must be ARM64, not $python_arch." >&2
    exit 1
fi

python -m PyInstaller --clean --noconfirm LabInterface-macos-arm64.spec

app_bundle="$app_root/dist/LabInterface.app"
app_binary="$app_bundle/Contents/MacOS/LabInterface"
archive="$app_root/dist/LabInterface-macOS-arm64.zip"

if [[ ! -x "$app_binary" ]]; then
    echo "PyInstaller did not create the expected application binary." >&2
    exit 1
fi

binary_archs="$(lipo -archs "$app_binary")"
if [[ "$binary_archs" != "arm64" ]]; then
    echo "Expected an ARM64 binary, but found: $binary_archs" >&2
    exit 1
fi

codesign --verify --deep --strict --verbose=2 "$app_bundle"
rm -f "$archive"
ditto -c -k --sequesterRsrc --keepParent "$app_bundle" "$archive"

echo "Created $archive"
