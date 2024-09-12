#!/bin/sh
set -euo pipefail

rm -rf docs
mkdir docs

zip -r docs/jsfoot_files.zip jsfoot_files
cp -r jsfoot_files docs/

quarto render
