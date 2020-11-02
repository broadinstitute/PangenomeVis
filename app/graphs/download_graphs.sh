#!/bin/bash

set -euxo pipefail

gsutil -m cp -r gs://broad-dsp-lrma-pgv/graphs/* .
