#!/usr/bin/env bash

set -euo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <archive-url> <archive-md5> <archive-path> <fasta-path>" >&2
    exit 2
fi

archive_url=$1
archive_md5=$2
archive_path=$3
fasta_path=$4
curl_bin=${CURL:-curl}
md5_bin=${MD5SUM:-md5sum}

for dependency in "$curl_bin" "$md5_bin" gzip; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "download-human-reference: required command not found: $dependency" >&2
        exit 1
    fi
done

if [[ ! $archive_md5 =~ ^[[:xdigit:]]{32}$ ]]; then
    echo "download-human-reference: invalid MD5 checksum: $archive_md5" >&2
    exit 2
fi

if [[ -s $fasta_path ]]; then
    echo "Human reference already present: $fasta_path"
    exit 0
fi

mkdir -p "$(dirname "$archive_path")" "$(dirname "$fasta_path")"

verify_archive() {
    printf '%s  %s\n' "$archive_md5" "$1" | "$md5_bin" --check --status -
}

if [[ -f $archive_path ]]; then
    if ! verify_archive "$archive_path"; then
        echo "download-human-reference: checksum mismatch: $archive_path" >&2
        echo "Remove that file and rerun make download." >&2
        exit 1
    fi
else
    partial_path=${archive_path}.part
    echo "Downloading GRCh38.p14 from NCBI (approximately 1 GB compressed)..."
    "$curl_bin" \
        --fail \
        --location \
        --retry 5 \
        --continue-at - \
        --output "$partial_path" \
        "$archive_url"
    if ! verify_archive "$partial_path"; then
        echo "download-human-reference: checksum mismatch: $partial_path" >&2
        echo "Remove that partial file and rerun make download." >&2
        exit 1
    fi
    mv "$partial_path" "$archive_path"
fi

echo "Verifying and unpacking the reference (approximately 3.3 GB)..."
gzip --test "$archive_path"
gzip --decompress --stdout "$archive_path" > "${fasta_path}.part"
mv "${fasta_path}.part" "$fasta_path"
rm -f "$archive_path"

echo "Human reference ready: $fasta_path"
