#!/usr/bin/env bash

###############################################################################
# Simple shell script to generate metrics documentation with scaladoc
###############################################################################

set -euo pipefail

# Uncomment the following line to turn on debugging when running scaladoc
# JAVA_OPTS="-agentlib:jdwp=transport=dt_socket,server=y,suspend=y,address=7777"
curdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

pushd ${curdir}/../../ >/dev/null 2>&1

sources=$(find ~+/src/main/scala -name \*.scala)
cp=$(./mill show tools.runClasspath 2> /dev/null | jq --raw-output .[] | cut -f 3- -d ':' | paste -sd ":" -)

scaladoc \
    -toolcp $cp \
    -d docs \
    -doc-generator com.fulcrumgenomics.sv.internal.FgSvMetricsDoclet \
    $sources

popd 2>&1
