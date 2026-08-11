#!/bin/bash

set -e

TESTDATADIR=$(dirname ${0});

# Decompress the test data
pushd $TESTDATADIR;
for FILE in $(ls -1 *gz)
do 
	gunzip -c ${FILE} > $(basename ${FILE} .gz);
done
popd

mkdir tmp

# Generate the new test data
./dwgsim -z 13 -N 10000 samtools/examples/ex1.fa tmp/ex1.test

# Test the differences
for GZFILE in $(ls -1 tmp/ex1.test*gz)
do 
	gunzip $GZFILE;
    FILE=$(basename $GZFILE .gz);
	diff -q tmp/${FILE} ${TESTDATADIR}/${FILE}
done

# Remove only the top-level files decompressed by this script. Do not recurse:
# testdata also contains committed region metadata and source assets.
for GZFILE in "${TESTDATADIR}"/*.gz
do
    rm -f "${GZFILE%.gz}"
done

rm -r tmp
