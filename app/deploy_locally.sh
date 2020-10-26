#!/bin/bash

cd client && yarn build && cd ..
uvicorn pangenome_backend:app --reload
