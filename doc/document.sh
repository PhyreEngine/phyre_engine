#!/bin/bash
top_dir=$(git rev-parse --show-toplevel)

# Exclude untracked files
git ls-files --exclude-standard --others "$top_dir/*.py" | readarray untracked

sphinx-apidoc \
    -f -e \
    -o "$top_dir/doc/source/apidoc/" \
    "$top_dir/phyre_engine" \
    "$top_dir/phyre_engine/test" \
    "${untracked[@]}"
