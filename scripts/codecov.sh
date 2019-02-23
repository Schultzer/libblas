#!/bin/bash

# Replace 'myproject' by your project name
REPORT=$(find target/debug -maxdepth 1 -name '<blas-rs-*' -a ! -name '*.d')

for file in $REPORT; do
    mkdir -p "target/cov/$(basename $file)"
    kcov --exclude-pattern=/.cargo,/usr/lib --verify "target/cov/$(basename $file)" "$file"
done

wget -O - -q "https://codecov.io/bash" > .codecov
chmod +x .codecov
./.codecov -t $CODECOV_TOKEN
echo "Uploaded code coverage"