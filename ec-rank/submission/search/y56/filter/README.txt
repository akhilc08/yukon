y56/filter -- can a learned filter beat the Mestre-Nagao sum S(X) at picking high-rank specialisations?
VERDICT: no, not by enough to matter. Keep using S over as many primes as you can afford. (Details below.)

FILES
  gen_feats.gp     sample family specialisations; features (a_p p<1000, S(100..1e4), disc, rootno) + rank label
  run_gen.sh       box driver for gen_feats.gp (6 workers, cores 22-27)
  data.py          loader (dedup by family+a-invariants; recomputes bad primes from the discriminant)
  train.py         hot-loop models on a_p of the first K primes: lin / quad (ridge), gbm (LightGBM), mlp
                   exp A = held-out parameter band, exp B = held-out family, *_f = same-family model,
                   meta = GBM on (S100,S316,S1000,log|D|,nbad,H[,rootno])
  add_s.gp         S(X) at X=100..1e5 (10 cutoffs) for curves already in a data file
  train_refine.py  per-family "refine S" models vs S(1e5) at matched a_p cost
  run_train.sh / run_s5.sh   box drivers; report.py aggregates; known_check.py = ledger sanity
  export.py        trained models -> model/y56_filter.h (plain C) + model/models.json
  scorer.py        python scorer (same models); bench.c = ns/candidate benchmark
  model/           y56_filter.h, models.json   (run r2; trained on all 326k labelled curves)

C API: #include "model/y56_filter.h"; int ap[168]: ap[i] = a_p of the i-th prime (2,3,5,...), 0 at bad primes.
  float y56_lin{50,100,168}(ap), y56_quad*, y56_gbm*, y56_mlp*, y56_S(ap, K)   (larger = better; compare within a family)
  Verified: C scores match python to 5e-6 (gbm thresholds are exact integer cuts on a_p).
Python: from scorer import Scorer; Scorer().score_ap(ap, "lin50") or .score_curves(ainvs, "lin50").

DATA (box, s3://ecrank-y56-115135126148/results/filter/)
  326k curves, 9 families x 2 parameter bands: Tate normal forms Z5 Z6 Z7 Z8 Z2Z4 Z2Z6 (t=u/v, max(|u|,v) <= 600 / 600..6000),
  Z2Z2 y^2=x(x+A)(x+B) (|A|,|B| <= 500 / 500..5000), random small Weierstrass "rand", Mestre ten-root
  (roots in [-30,30] / [-50,50]).  ~20k per (family,band); Mestre 9.5k.
  Label: ellrank(E,0) lower bound r (upper bound R recorded; r<R in 0.2%-20% of rows, worst for Z8/Z7/Z2Z6 band 2).
  Mestre label: #independent points from hyperellratpoints(quartic, 5000) (mestre-own/mestre.gp) -- this is
  "rank the cheap point search certifies", base 9, not the true rank (ellrank takes minutes per Mestre curve).
  S(1e5) partial sums for 120k curves (Z6, Z2Z2, rand, Z2Z4) + 11.5k Mestre.
  Targets: rank >= median + dt (median 1 for torsion/rand families, 9 for Mestre).

RESULTS (hits in top-k, summed over pools; k = k/N of each test pool)
 1. Hot loop, same K primes (exp A, held-out band, 15-16 torsion/rand pools), dt=2 (rank >= 3):
      k/N=1%:  K=50: S_K 526, lin 493, quad 484, gbm 407, mlp 333, lin_f 514 | S(1e3) 743, S(1e4) 1082
               K=168: S_K 743, lin 643, quad 540, gbm 446, lin_f 683
    held-out family (exp B), K=168: lin 606, quad 521, gbm 399 vs S_K 743.   -> learned <= S at equal cost.
    Mestre K=50, same-family lin_f: 32 vs S_K 23 (dt=2, 1%), 54 vs 37 (dt=1) -- but S(1e3) 42/57 is still better,
    and at K=168 lin_f = S_K (41 vs 42).
 2. Refine S(1e5), per family, matched a_p cost (features: S at 10 cutoffs, slope/curvature in log X, a_p p<=31,
    nbad, log|D|), band-heldout and random splits, sum over Z6/Z2Z2/rand/Z2Z4:
      dt=2 k/N=0.5%: S1e5 938, ridge 997 (+6%), gbm 980 (+4%), gbm+rootno 1011 (+8%)
      dt=2 k/N=1%:   S1e5 1554, ridge 1645 (+6%), gbm 1596 (+3%)
      dt=3: saturated (S1e5 already finds all 141 rank>=4 curves in the top 0.5%).
    Per family best case Z2Z4 +7..15%, rand +11% at one budget, Z2Z2 ~+1%.  Below the 10% bar.
    Mestre: gbm / S+corr +10..24% over S1e5 (20/20 random splits win), BUT the gain is entirely log|D| and nbad
    (a_p-only correction: +3.6%).  A 3-coefficient rule captures most of it:
        mestre:  score = S(1e5) - 0.103*log|disc| + 1.49*nbad(p<1000)      (+10..16% vs S1e5 on findable rank)
    Smaller-discriminant Mestre curves have smaller points, so this partly measures how findable the points are
    by the height-5000 quartic search, not the true rank. It still lines up with what the challenge needs
    (certifiable points, small H). Not confirmed with real ellrank labels.
    Torsion families: a consistent +1.0..1.5*nbad term (about +1..7%).
 3. meta (S100,S316,S1e3,log|D|,nbad,H): in-family beats S1e3 but not S1e4; held-out family it fails (Mestre 15 vs 57).
 4. Known-curve ledger (5586, selection-biased): Spearman(score, rank) per torsion group; quad168 >= S168 in most
    groups (e.g. trivial torsion 0.85 vs 0.74), but that did not carry over to the unbiased in-family tests above.
 Speed (c7g, loaded box, 1 core, ns/candidate): S_50 27, lin50 27, quad50 32, mlp50 227, gbm50 1417;
    K=168: S 82, lin 81, quad 83, mlp 748, gbm 1639.  A linear model costs the same as S.
 CPU sanity: in these families ellrank(E,0) costs 3-50 ms (Z7 band 2: 224 ms) and S(1e5) in gp costs ~120 ms, so at
    these heights a sieve beyond S(1e3) does not pay. Ranking quality only matters where ellrank is expensive.

REPRODUCE (box, as ubuntu; code at /work/ec-rank/submission/search/y56/filter)
  NCH=6 NPER=10000 NMES=2500 bash run_gen.sh /work/y56/filter/data 22-27   # stopped after ~2.3 rounds (326k rows)
  bash run_train.sh r2 50 100 168                     # exp A/B/meta, export header, ledger check, bench
  FAMS="Z6 Z2Z2 rand Z2Z4" NMAX=7500 bash run_s5.sh   # S(1e5) + per-family refine -> runs/refine/report.txt
  Mestre refine: same as run_s5.sh with mestre.*.txt, NMAX=3000 -> runs/refine_mestre (script runs/refine_mestre/run_s5m.sh)
  python report.py runs/r2 168 ; python known_check.py known_ap.txt model/models.json
