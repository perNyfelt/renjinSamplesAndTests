#!/usr/bin/env bash

basedir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

srcDir="${basedir}/src"
outDir="${basedir}/out"
libDir="${basedir}/lib"

if [[ ! -d "${libDir}" ]]; then
  echo "Fetching renjin dependencies into ${libDir}"
  ant resolve
fi

java -version
javac -version

echo "Compiling..."
if [[ ! -d "${outDir}" ]]; then
  mkdir "${outDir}"
fi
javac -g -cp "${libDir}/*" -d "${outDir}" -Xlint:unchecked $(find "${srcDir}"/* | grep .java)

echo "Running SimpleTest"
java -cp "${libDir}/*:${outDir}" se.alipsa.SimpleTest