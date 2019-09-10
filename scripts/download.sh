#!/bin/sh

for file in $(ls gdc_manifest_tcga_*); do ./gdc-client download -t ../gdc-user-token.2019-08-13T20_11_22.889Z.txt -m $file -n 8 --log-file gdc.log -d "tcga-data/${file/gdc_manifest_tcga_/}"; done

