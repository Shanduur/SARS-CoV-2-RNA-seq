#!/usr/bin/env bash

URL_REF=https://zenodo.org/record/6342228/files/ref.Rds?download=1
URL_IDX=https://zenodo.org/record/6342228/files/idx.annoy?download=1

echo "downloading references: ${URL_REF}"
curl -sSL ${URL_REF} -o ref.Rds

echo "downloading identities: ${URL_IDX}"
curl -sSL ${URL_IDX} -o idx.annoy

echo "done!"
