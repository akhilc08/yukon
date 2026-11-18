# mkec.py HITFILE OUT.ec : convert a rankrun hit file (3 lines: meta, ainvs, points) to a board-format .ec
import sys, re, datetime
meta, ai, pts = [l.strip() for l in open(sys.argv[1]) if l.strip()][:3]
m = [s.strip() for s in meta.strip('[]').split(',')]
K, p, q, H, r1, r2, j = m
a = [s.strip() for s in ai.strip('[]').split(',')]
P = re.findall(r'\[([^\[\],]+),\s*([^\[\],]+)\]', pts)
fam = {'808': 'X808 rank-9 elliptic K3 fibration: y^2 + 84(35397T^2-221T+20)xy = x(x+(1139T+115)(1269T-59)(3776T+34)(539T-5))(x+8(638T-23)(250T-1)(7811T+499)) (equation from public ICARM commentary by wgxli under curve #4078; fibre, sections, points our own)',
       '552': 'X552 rank-9 elliptic K3 fibration: y^2 + 238(19T-6)(5T-3)xy = x(x+(12T-29)(51T+22)(79T+22)(482T+29))(x+256(T-22)(83T-7)(21T-2)(7T+17)) (equation from public ICARM commentary by wgxli under curve #4061; fibre, sections, points our own)'}[K]
model = 'model of Y^2=x^3+A x^2+B x' if j == '1' else '2-isogenous model y^2=x^3-2A x^2+(A^2-4B)x'
with open(sys.argv[2], 'w') as f:
    f.write(f"# Own find: Z/2 fibre T = {p}/{q} of the {fam}.\n")
    f.write(f"# Search: y56/z2r20 height-first sieve (C ksieve: exact odd-prime + 2-adic minimal-model height bound, S(200)->S(1000)->S(30000)), PARI S(1e5), ellrank seeded with our own Q(T) sections (x = k*prod of discriminant factors), ellsaturation(E,L,100).\n")
    f.write(f"# Minimal {model} at the fibre; H = {H}, ellrank bounds [{r1},{r2}]. PARI/GP 2.17.4. Hit log: {sys.argv[1]}. Written {datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%MZ')}.\n")
    f.write(f"#CLAIM rank={r1}\n")
    f.write("E " + " ".join(a) + "\n")
    for x, y in P:
        f.write(f"P {x.strip()} {y.strip()}\n")
print(sys.argv[2], 'points', len(P))
