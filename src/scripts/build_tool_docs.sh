#!/bin/bash

###############################################################################
# Simple shell script to generate tool documentation
###############################################################################

set -euo pipefail

# Uncomment the following line to turn on debugging when running scaladoc
# JAVA_OPTS="-agentlib:jdwp=transport=dt_socket,server=y,suspend=y,address=7777"
curdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

pushd ${curdir}/../../ >/dev/null 2>&1

./mill tools.deployLocal

java -cp jars/fgsv.jar com.fulcrumgenomics.sv.internal.FgSvInternalMain BuildToolDocs -o docs/tools

popd 2>&1
