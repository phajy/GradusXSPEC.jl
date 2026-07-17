#!/usr/bin/env bash
# Diagnose the GradusXSPEC build environment without changing anything.
#
# Reports PASS / WARN / FAIL for each prerequisite so a broken build produces
# one clear summary instead of a mid-build error. Read-only: it never mutates
# HEASOFT, Homebrew, or the working tree.
#
# Usage (ideally with HEADAS sourced):
#   ./check-env.sh
#
# Exit status: 0 if there are no hard failures, 1 otherwise. WARNings do not
# affect the exit status.

set -uo pipefail

cd "$(dirname "$0")" || exit 1

fails=0
warns=0

pass() { printf '  [PASS] %s\n' "$*"; }
warn() { printf '  [WARN] %s\n' "$*"; warns=$((warns + 1)); }
fail() { printf '  [FAIL] %s\n' "$*"; fails=$((fails + 1)); }
section() { printf '\n== %s ==\n' "$*"; }

section "Platform"
OS="$(uname -s)"
pass "OS: ${OS} ($(uname -m))"

section "Julia"
if command -v julia >/dev/null 2>&1; then
  jver="$(julia --version 2>/dev/null || echo '?')"
  pass "julia found: ${jver} ($(command -v julia))"
else
  fail "julia not found on PATH (needed by build-julia.sh)"
fi

if [[ -f Project.toml ]]; then
  pass "Project.toml present"
else
  fail "Project.toml missing (run from the repository root)"
fi

section "HEASOFT / XSPEC"
if [[ -n "${HEADAS:-}" ]]; then
  if [[ -d "$HEADAS" ]]; then
    pass "HEADAS set: $HEADAS"
  else
    fail "HEADAS set but directory does not exist: $HEADAS"
  fi
else
  fail "HEADAS not set (source \$HEADAS/headas-init.sh)"
fi

for cmd in initpackage hmake xspec; do
  if command -v "$cmd" >/dev/null 2>&1; then
    pass "$cmd on PATH ($(command -v "$cmd"))"
  else
    fail "$cmd not found on PATH (is HEADAS initialized?)"
  fi
done

section "Build inputs"
if [[ -f model.dat ]]; then
  pass "model.dat present (build-julia.sh has run)"
else
  warn "model.dat not found (run ./build-julia.sh before ./build-xspec.sh)"
fi

if [[ -f build/lib/libGradusXSPEC.so || -f build/lib/libGradusXSPEC.dylib ]]; then
  pass "compiled Julia library present under build/lib"
else
  warn "build/lib/libGradusXSPEC.{so,dylib} not found (run ./build-julia.sh)"
fi

if [[ -f xillverD-5.fits ]]; then
  pass "reflection table xillverD-5.fits present"
else
  warn "xillverD-5.fits not found in repo root (models need a reflection table)"
fi

# macOS-specific: the gcc@14 Fortran-library workaround for HEASOFT linking.
if [[ "$OS" == "Darwin" ]]; then
  section "macOS Fortran linking (gcc@14)"
  if command -v brew >/dev/null 2>&1; then
    pass "Homebrew found ($(command -v brew))"
    gcc_prefix="$(brew --prefix gcc@14 2>/dev/null || true)"
    if [[ -n "$gcc_prefix" && -d "$gcc_prefix" ]]; then
      pass "gcc@14 installed: $gcc_prefix"
      gfortran="${gcc_prefix}/bin/gfortran-14"
      [[ -x "$gfortran" ]] || gfortran="$(command -v gfortran-14 || true)"
      if [[ -n "$gfortran" && -x "$gfortran" ]]; then
        pass "gfortran-14 found: $gfortran"
        ok=true
        for lib in libemutls_w.a libgfortran.dylib libgcc.a; do
          path="$("$gfortran" -print-file-name="$lib" 2>/dev/null || true)"
          if [[ -z "$path" || ! -f "$path" ]]; then
            warn "gfortran cannot locate $lib (build-xspec.sh will compute the override)"
            ok=false
          fi
        done
        [[ "$ok" == true ]] && pass "gfortran reports all required Fortran libraries"
      else
        fail "gfortran-14 not found (brew install gcc@14)"
      fi
    else
      fail "gcc@14 not installed (brew install gcc@14)"
    fi
  else
    warn "Homebrew not found; XSPEC local-model linking usually needs gcc@14 on macOS"
  fi
fi

section "Summary"
if [[ $fails -gt 0 ]]; then
  printf 'FAIL: %d failure(s), %d warning(s).\n' "$fails" "$warns"
  exit 1
elif [[ $warns -gt 0 ]]; then
  printf 'OK with %d warning(s). Address warnings if a build step fails.\n' "$warns"
  exit 0
else
  printf 'All checks passed.\n'
  exit 0
fi
