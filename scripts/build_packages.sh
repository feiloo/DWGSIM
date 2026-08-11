#!/usr/bin/env bash

set -euo pipefail

package_name=dwgsim
action=${1:-}
package_version=${2:-}
requested_package_root=${3:-build/packages}

die() {
    printf 'error: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "required command not found: $1"
}

case "$action" in
    dist|rpm|deb) ;;
    *) die 'usage: build_packages.sh {dist|rpm|deb} VERSION PACKAGE_ROOT' ;;
esac

[[ "$package_version" =~ ^[0-9A-Za-z][0-9A-Za-z._+~-]*$ ]] ||
    die "unsafe package version: $package_version"

require_command realpath
source_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)

is_git_worktree_root() {
    local directory=$1
    local top_level

    command -v git >/dev/null 2>&1 || return 1
    top_level=$(git -C "$directory" rev-parse --show-toplevel 2>/dev/null) || return 1
    top_level=$(realpath -e -- "$top_level") || return 1
    [[ "$top_level" == "$directory" ]]
}

case "$requested_package_root" in
    /*) package_root_candidate=$requested_package_root ;;
    *) package_root_candidate=$source_root/$requested_package_root ;;
esac
package_root_candidate=$(realpath -m -- "$package_root_candidate")

if [[ "$package_root_candidate" == "$source_root" ||
      "$package_root_candidate" != "$source_root/"* ||
      "$package_root_candidate" == "$source_root/.git" ||
      "$package_root_candidate" == "$source_root/.git/"* ]]; then
    die "PACKAGE_ROOT must be a non-.git subdirectory of $source_root"
fi

mkdir -p -- "$package_root_candidate"
package_root=$(cd -- "$package_root_candidate" && pwd -P)
[[ "$package_root" == "$source_root/"* ]] ||
    die "PACKAGE_ROOT resolves outside the source tree: $package_root"

artifact_dir=$package_root/artifacts
work_root=$package_root/work
tmp_root=$package_root/tmp
local_home=$package_root/home
cache_root=$package_root/cache

ensure_local_directory() {
    local directory=$1
    local resolved_directory

    [[ ! -L "$directory" ]] || die "managed package directory cannot be a symlink: $directory"
    [[ ! -e "$directory" || -d "$directory" ]] ||
        die "managed package path is not a directory: $directory"
    mkdir -p -- "$directory"
    resolved_directory=$(realpath -e -- "$directory")
    [[ "$resolved_directory" == "$package_root/"* ]] ||
        die "managed package directory resolves outside PACKAGE_ROOT: $directory"
}

ensure_local_directory "$artifact_dir"
ensure_local_directory "$work_root"
ensure_local_directory "$tmp_root"
ensure_local_directory "$local_home"
ensure_local_directory "$cache_root"

export HOME=$local_home
export TMPDIR=$tmp_root
export TMP=$tmp_root
export TEMP=$tmp_root
export XDG_CACHE_HOME=$cache_root
export LC_ALL=C

if [[ -z "${SOURCE_DATE_EPOCH:-}" ]]; then
    if is_git_worktree_root "$source_root"; then
        SOURCE_DATE_EPOCH=$(git -C "$source_root" log -1 --format=%ct)
    else
        SOURCE_DATE_EPOCH=0
    fi
fi
[[ "$SOURCE_DATE_EPOCH" =~ ^[0-9]+$ ]] || die 'SOURCE_DATE_EPOCH must be an integer'
export SOURCE_DATE_EPOCH

package_jobs=${PACKAGE_JOBS:-}
if [[ -z "$package_jobs" ]]; then
    package_jobs=$(getconf _NPROCESSORS_ONLN 2>/dev/null || printf '1')
fi
[[ "$package_jobs" =~ ^[1-9][0-9]*$ ]] || die 'PACKAGE_JOBS must be a positive integer'
package_cflags=${PACKAGE_CFLAGS:--g -Wall -O3}

run_dir=$(mktemp -d "$work_root/$action.XXXXXXXX")
cleanup() {
    if [[ -n "${run_dir:-}" && "$run_dir" == "$work_root/"* ]]; then
        rm -rf -- "$run_dir"
    fi
}
trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM

source_dir_name=$package_name-$package_version
source_archive_name=$source_dir_name.tar.gz
source_archive=$artifact_dir/$source_archive_name

declare -A manifest_seen=()
add_source_file() {
    local relative_path=${1#./}

    [[ -n "$relative_path" && "$relative_path" != samtools ]] || return 0
    if [[ ( -f "$source_root/$relative_path" || -L "$source_root/$relative_path" ) &&
          -z "${manifest_seen[$relative_path]+present}" ]]; then
        printf '%s\0' "$relative_path" >> "$source_manifest"
        manifest_seen[$relative_path]=1
    fi
}

build_source_manifest() {
    local file relative_path
    source_manifest=$run_dir/source-files.list
    : > "$source_manifest"

    if is_git_worktree_root "$source_root"; then
        while IFS= read -r -d '' file; do
            add_source_file "$file"
        done < <(git -C "$source_root" ls-files -z)

        if [[ ! -d "$source_root/samtools" ]] ||
           ! is_git_worktree_root "$source_root/samtools"; then
            die 'the samtools submodule is not initialized'
        fi
        while IFS= read -r -d '' file; do
            add_source_file "samtools/$file"
        done < <(git -C "$source_root/samtools" ls-files -z)
    else
        while IFS= read -r -d '' file; do
            relative_path=${file#"$source_root/"}
            case "$relative_path" in
                samtools/win32/libcurses.a|samtools/win32/libz.a) ;;
                *.o|*.a|*.so|*.so.*|*.dylib|*.exe|gmon.out|a.out|*~|\
                dwgsim|dwgsim_eval|dwgsim_mut_to_vcf|dwgsim_pileup_eval|tests/run_tests|\
                samtools/samtools|samtools/bgzip|samtools/razip) continue ;;
            esac
            add_source_file "$relative_path"
        done < <(
            find "$source_root" \
                \( -path "$source_root/.git" -o -path "$source_root/build" \
                   -o -path "$source_root/reference" -o -path "$package_root" \) -prune -o \
                \( -type f -o -type l \) -print0
        )
    fi

    # These files must be available before they have first been added to Git.
    add_source_file packaging/dwgsim.spec.in
    add_source_file scripts/build_packages.sh
}

build_dist() {
    local dist_parent=$run_dir/dist
    local dist_tree=$dist_parent/$source_dir_name
    local uncompressed_archive=$run_dir/$source_dir_name.tar
    local completed_archive=$run_dir/$source_archive_name

    require_command tar
    require_command gzip
    build_source_manifest
    mkdir -p -- "$dist_tree"
    tar -C "$source_root" --null --files-from="$source_manifest" \
        -cf "$run_dir/source-files.tar"
    tar -C "$dist_tree" -xf "$run_dir/source-files.tar"

    [[ -f "$dist_tree/Makefile" && -f "$dist_tree/LICENSE" &&
       -f "$dist_tree/samtools/Makefile" ]] ||
        die 'source staging is incomplete'

    tar --sort=name --format=gnu --mtime="@$SOURCE_DATE_EPOCH" \
        --owner=0 --group=0 --numeric-owner \
        -C "$dist_parent" -cf "$uncompressed_archive" "$source_dir_name"
    gzip -n -9 -c "$uncompressed_archive" > "$completed_archive"
    mv -f -- "$completed_archive" "$source_archive"
    printf 'Created %s\n' "$source_archive"
}

write_debian_source_control() {
    local control_file=$1

    {
        printf 'Source: %s\n' "$package_name"
        printf 'Section: science\n'
        printf 'Priority: optional\n'
        printf 'Maintainer: DWGSIM contributors <dwgsim@users.noreply.github.com>\n'
        printf 'Standards-Version: 4.6.2\n'
        printf 'Rules-Requires-Root: no\n\n'
        printf 'Package: %s\n' "$package_name"
        printf 'Architecture: any\n'
        printf 'Depends: ${shlibs:Depends}\n'
        printf 'Description: whole-genome simulator for next-generation sequencing\n'
        printf ' DWGSIM generates synthetic sequencing reads and mutation truth data.\n'
    } > "$control_file"
}

debian_dependencies() {
    local source_tree=$1
    local binary_root=debian/$package_name/usr/bin
    local dependency_output dependencies
    local dpkg_shlibdeps=${DPKG_SHLIBDEPS:-dpkg-shlibdeps}

    if command -v "$dpkg_shlibdeps" >/dev/null 2>&1; then
        mkdir -p -- "$source_tree/debian"
        write_debian_source_control "$source_tree/debian/control"
        if dependency_output=$(
            cd -- "$source_tree"
            "$dpkg_shlibdeps" -O \
                "$binary_root/dwgsim" \
                "$binary_root/dwgsim_eval" \
                "$binary_root/dwgsim_mut_to_vcf" \
                "$binary_root/dwgsim_pileup_eval"
        ); then
            dependencies=$(printf '%s\n' "$dependency_output" |
                sed -n 's/^shlibs:Depends=//p')
            if [[ -n "$dependencies" ]]; then
                printf '%s\n' "$dependencies"
                return 0
            fi
        fi
    fi

    printf 'libc6, zlib1g\n'
}

build_deb() {
    local dpkg_deb=${DPKG_DEB:-dpkg-deb}
    local dpkg_command=${DPKG:-dpkg}
    local deb_work=$run_dir/deb
    local source_tree=$deb_work/$source_dir_name
    local package_tree=$source_tree/debian/$package_name
    local control_dir=$package_tree/DEBIAN
    local doc_dir=$package_tree/usr/share/doc/$package_name
    local architecture dependencies installed_size deb_version_filename
    local reproducible_cflags
    local completed_package

    require_command "$dpkg_deb"
    require_command "$dpkg_command"
    "$dpkg_command" --validate-version "$package_version"
    architecture=${DEB_ARCH:-$("$dpkg_command" --print-architecture)}
    [[ "$architecture" =~ ^[0-9A-Za-z][0-9A-Za-z-]*$ ]] ||
        die "unsafe Debian architecture: $architecture"

    build_dist
    mkdir -p -- "$deb_work"
    tar -xzf "$source_archive" -C "$deb_work"
    reproducible_cflags="$package_cflags -ffile-prefix-map=$source_tree=. -fdebug-prefix-map=$source_tree=."
    make -C "$source_tree" -j "$package_jobs" \
        PACKAGE_VERSION="$package_version" CFLAGS="$reproducible_cflags" lib-recur
    make -C "$source_tree" -j "$package_jobs" \
        PACKAGE_VERSION="$package_version" CFLAGS="$reproducible_cflags" all

    install -d -m 0755 "$control_dir" "$package_tree/usr/bin" "$doc_dir/docs"
    for program in dwgsim dwgsim_eval dwgsim_mut_to_vcf dwgsim_pileup_eval; do
        install -m 0755 "$source_tree/$program" "$package_tree/usr/bin/$program"
    done
    install -m 0644 "$source_tree/README.md" "$source_tree/INSTALL" "$doc_dir/"
    install -m 0644 "$source_tree/LICENSE" "$doc_dir/copyright"
    install -m 0644 "$source_tree"/docs/*.md "$doc_dir/docs/"
    dependencies=$(debian_dependencies "$source_tree")

    installed_size=$(du -sk "$package_tree/usr" | awk '{print $1}')
    {
        printf 'Package: %s\n' "$package_name"
        printf 'Version: %s\n' "$package_version"
        printf 'Architecture: %s\n' "$architecture"
        printf 'Maintainer: DWGSIM contributors <dwgsim@users.noreply.github.com>\n'
        printf 'Installed-Size: %s\n' "$installed_size"
        printf 'Depends: %s\n' "$dependencies"
        printf 'Section: science\n'
        printf 'Priority: optional\n'
        printf 'Homepage: https://github.com/nh13/DWGSIM\n'
        printf 'Description: whole-genome simulator for next-generation sequencing\n'
        printf ' DWGSIM generates synthetic sequencing reads and mutation truth data.\n'
    } > "$control_dir/control"
    chmod 0644 "$control_dir/control"

    find "$package_tree" -print0 |
        xargs -0r touch -h -d "@$SOURCE_DATE_EPOCH"
    deb_version_filename=${package_version//:/_}
    completed_package=$run_dir/${package_name}_${deb_version_filename}_${architecture}.deb
    "$dpkg_deb" --root-owner-group --build "$package_tree" "$completed_package"
    mv -f -- "$completed_package" "$artifact_dir/$(basename -- "$completed_package")"
    printf 'Created %s\n' "$artifact_dir/$(basename -- "$completed_package")"
}

build_rpm() {
    local rpmbuild=${RPMBUILD:-rpmbuild}
    local rpm_work=$run_dir/rpm
    local rpm_version=${package_version//-/.}
    local spec_template=$source_root/packaging/dwgsim.spec.in
    local rendered_spec=$rpm_work/SPECS/dwgsim.spec
    local rpm_file
    local -a rpm_files=()

    require_command "$rpmbuild"
    require_command sed
    [[ -f "$spec_template" ]] || die 'missing packaging/dwgsim.spec.in'
    [[ "$rpm_version" =~ ^[0-9A-Za-z][0-9A-Za-z._+~]*$ ]] ||
        die "package version cannot be converted to an RPM version: $package_version"

    build_dist
    mkdir -p -- "$rpm_work"/{BUILD,BUILDROOT,RPMS,SOURCES,SPECS,SRPMS}
    install -m 0644 "$source_archive" "$rpm_work/SOURCES/$source_archive_name"
    sed -e "s|@UPSTREAM_VERSION@|$package_version|g" \
        -e "s|@RPM_VERSION@|$rpm_version|g" \
        "$spec_template" > "$rendered_spec"

    "$rpmbuild" -ba \
        --define "_topdir $rpm_work" \
        --define "_builddir $rpm_work/BUILD" \
        --define "_buildrootdir $rpm_work/BUILDROOT" \
        --define "_rpmdir $rpm_work/RPMS" \
        --define "_sourcedir $rpm_work/SOURCES" \
        --define "_specdir $rpm_work/SPECS" \
        --define "_srcrpmdir $rpm_work/SRPMS" \
        --define "_tmppath $tmp_root" \
        --define "_smp_mflags -j$package_jobs" \
        "$rendered_spec"

    while IFS= read -r -d '' rpm_file; do
        rpm_files+=("$rpm_file")
    done < <(find "$rpm_work/RPMS" "$rpm_work/SRPMS" -type f -name '*.rpm' -print0)
    (( ${#rpm_files[@]} > 0 )) || die 'rpmbuild did not create an RPM'
    for rpm_file in "${rpm_files[@]}"; do
        install -m 0644 "$rpm_file" "$artifact_dir/$(basename -- "$rpm_file")"
        printf 'Created %s\n' "$artifact_dir/$(basename -- "$rpm_file")"
    done
}

case "$action" in
    dist) build_dist ;;
    deb) build_deb ;;
    rpm) build_rpm ;;
esac
