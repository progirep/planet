#!/usr/bin/env sh
set -e

~/LocalLibs/caffe/build/tools/caffe train --solver=lenet_solver.prototxt $@
