#!/bin/sh
# FORD HTML in build/doc; the user manual is llms.txt without its title line.
# Needs pip install ford. Stands in for the reference's make doc.
set -e
cd "$(dirname "$0")/.."
mkdir -p build/pages
{ printf -- '---\ntitle: User manual\n---\n\n'; sed 1d llms.txt; } > build/pages/index.md
ford ford.md
cp llms.txt build/doc/
