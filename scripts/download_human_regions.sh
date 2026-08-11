#!/usr/bin/env bash

set -euo pipefail

if (( $# == 0 || $# % 3 != 0 )); then
    echo "Usage: $0 <url> <md5> <path> [<url> <md5> <path>]..." >&2
    exit 2
fi

curl_bin=${CURL:-curl}
md5_bin=${MD5SUM:-md5sum}

for dependency in "$curl_bin" "$md5_bin" gzip; do
    if ! command -v "$dependency" >/dev/null 2>&1; then
        echo "download-human-regions: required command not found: $dependency" >&2
        exit 1
    fi
done

verify_file() {
    local expected_md5=$1
    local path=$2
    printf '%s  %s\n' "$expected_md5" "$path" |
        "$md5_bin" --check --status -
}

download_one() {
    local url=$1
    local expected_md5=$2
    local path=$3
    local partial_path=${path}.part

    if [[ ! $expected_md5 =~ ^[[:xdigit:]]{32}$ ]]; then
        echo "download-human-regions: invalid MD5 checksum: $expected_md5" >&2
        exit 2
    fi

    mkdir -p "$(dirname "$path")"
    if [[ -f $path ]]; then
        if ! verify_file "$expected_md5" "$path"; then
            echo "download-human-regions: checksum mismatch: $path" >&2
            exit 1
        fi
    else
        echo "Downloading $url"
        "$curl_bin" \
            --fail \
            --location \
            --retry 5 \
            --output "$partial_path" \
            "$url"
        if ! verify_file "$expected_md5" "$partial_path"; then
            echo "download-human-regions: checksum mismatch: $partial_path" >&2
            exit 1
        fi
        mv "$partial_path" "$path"
    fi

    if [[ $path == *.gz ]]; then
        gzip --test "$path"
    fi
    echo "Verified human-region source: $path"
}

while [[ $# -gt 0 ]]; do
    download_one "$1" "$2" "$3"
    shift 3
done
