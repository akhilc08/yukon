#!/bin/bash
# verify.sh CAND.ec : novelty-on single-file verify, then board check against a temp copy of the board
cd /work/ec-rank
PY=${PY:-python3}
$PY -m ec_verify "$1" --no-conductor
TMP=$(mktemp -d); cp submission/board/*.ec $TMP/; cp "$1" $TMP/zz-cand.ec
$PY -m ec_verify $TMP --no-conductor | tail -8
rm -rf $TMP
