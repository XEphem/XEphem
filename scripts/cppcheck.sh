#!/bin/sh
cppcheck --enable=all --inline-suppr --suppress=missingIncludeSystem \
--cppcheck-build-dir=../build --inconclusive --library=posix $*
