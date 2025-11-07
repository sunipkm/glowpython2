#!/usr/bin/env bash
rm -rf build
pip install . --config-settings=builddir=build &> build.log