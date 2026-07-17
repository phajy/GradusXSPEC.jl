#!/usr/bin/env bash
# Source the HEASoft environment, then exec the given command.
#
# The official HEASoft container configures HEADAS for interactive shells, but
# Docker `RUN`/`CMD` steps do not always inherit it. This wrapper makes the
# environment explicit and robust: it uses $HEADAS if already exported, and
# otherwise locates headas-init.sh under the install prefix.
#
# Usage:
#   docker/with-headas.sh <command> [args...]

set -euo pipefail

if [[ -z "${HEADAS:-}" ]]; then
  init="$(find /opt /usr/local /heasoft -maxdepth 4 -name headas-init.sh -print -quit 2>/dev/null || true)"
  if [[ -z "$init" ]]; then
    echo "error: HEADAS is unset and headas-init.sh could not be located" >&2
    exit 1
  fi
  HEADAS="$(cd "$(dirname "$init")" && pwd)"
  export HEADAS
fi

# shellcheck disable=SC1091
source "${HEADAS}/headas-init.sh"

exec "$@"
