#!/bin/bash
# mkcand.sh JOBDIR P Q OUTNAME "FAMILY DESCRIPTION" [ISO]: certify the .pts record (p,q) of JOBDIR into cands/OUTNAME.ec.
# ISO=4 or 22: write the 2-isogenous curve with torsion [4] / [2,2] instead (points mapped through the isogeny, re-saturated).
set -e
cd /work/y56/tors; mkdir -p cands
D=$1; P=$2; Q=$3; OUT=$4; FAMD=$5; ISO=${6:-0}
GP="/work/ec-rank/tools/pari/bin/gp -q -s 1G -D parisizemax=4G"
cat $D/rin_0?.out.pts | grep -E "^\[$P, $Q," | head -1 > cands/$OUT.rec
[ -s cands/$OUT.rec ] || { echo "no record"; exit 1; }
DATE=$(date -u +%Y-%m-%dT%H:%MZ)
taskset -c 31 $GP <<GPEOF
read("certify.gp");
rec=eval(readstr("cands/$OUT.rec")[1]);
E=ellinit(rec[4]); P=rec[5];
read("isocert.gp");
if($ISO==4, rec=isomap(rec,[4])); if($ISO==22, rec=isomap(rec,[2,2]));
hdr=["Own find (y56/tors): family $FAMD, parameter t = $P/$Q.",
     Str("Search: y56/tors psieve/gsieve height-first Mestre-Nagao sieve + seeded ellrank; job dir /work/y56/tors/$D (s3 results/tors/$D). ISO=$ISO."),
     Str("PARI/GP ", version(), "; minimal model; ellrank + ellsaturation(E,L,100). Written $DATE.")];
n=cert(rec, "cands/$OUT.ec", hdr); print("points ",n);
GPEOF
cat cands/$OUT.ec | head -5
