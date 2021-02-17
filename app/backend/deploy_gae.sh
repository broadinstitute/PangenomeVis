#!/bin/bash

cd client && yarn build && cd ..

export GOOGLE_CLOUD_PROJECT=broad-dsp-lrma
gcloud app deploy app.yaml
