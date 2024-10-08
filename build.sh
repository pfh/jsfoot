#!/bin/sh
set -euo pipefail

rm -rf docs
mkdir docs

zip -r docs/jsfoot.zip jsfoot
cp -r jsfoot docs/

quarto render
