#!/bin/bash
#
# Print all of the embedded version numbers in XEphem's code and
# documentation, so the reader can see whether they match.  Particularly
# useful when doing a release.

set -e
cd "$(readlink -f $(dirname "${BASH_SOURCE[0]}"))"
cd ..
for s in archive/refs/tags 'tar xfz' 'cd xephem'
do
    grep -H "$s" INSTALL
done
grep -H PATCHLEVEL GUI/xephem/patchlevel.c
grep -H '^Version' xephem.spec
