#!/bin/bash

export GOOGLE_CLOUD_PROJECT=broad-dsp-lrma
gcloud app deploy app.yaml index.yaml
