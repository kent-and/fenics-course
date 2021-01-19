#!/bin/sh

echo "Don't run this script as is, as it will likely upload a bunch of output files etc"
exit 1

for f in `ls | grep -v README | grep -v upload.sh`; do
    echo "Uploading $f..."
    scp -r $f fenics@fenicsproject.org:fenicsproject.org/pub/course/src/
done
