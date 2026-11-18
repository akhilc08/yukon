# mkec.py HITFILE OUT.ec : convert a rankek.gp hit file (meta, ainvs, points) to a board-format .ec
import sys, re, datetime
meta, ai, pts = [l.strip() for l in open(sys.argv[1]) if l.strip()][:3]
U, M, p, q, H, r1, r2, j = [s.strip().strip('"') for s in meta.strip('[]').split(',')]
a = [s.strip() for s in ai.strip('[]').split(',')]
P = re.findall(r'\[([^\[\],]+),\s*([^\[\],]+)\]', pts)
nsec = 8 if M == '0' else 9
model = 'model of y^2=x^3+2A x^2+B x' if j == '1' else "2-isogenous model y^2=x^3-4A x^2+(4A^2-4B)x"
with open(sys.argv[2], 'w') as f:
    f.write("# z2ek: Elkies-Klagsbrun 2020 Z/2 family E_u: y^2 = x^3 + 2A(t,u) x^2 + B(t,u) x (arXiv:2003.00077 sec 9 eqs 3-4), rank >= 8 over Q(t) for any u, 9 on u = 2(m^2-m-1)/(m^2+1) (App. A eqs 7-8 sections)\n")
    f.write(f"# own specialisation: u = {U}" + (f" (m = {M})" if M != '0' else "") + f", t = {p}/{q}; height-first Mestre-Nagao sieve S(200)->S(1000)->S(30000) (y56/z2ek ksieve), PARI 2.17.4 S(1e5), ellrank seeded with the {nsec} sections, ellsaturation(E,L,100)\n")
    f.write(f"# minimal {model} at the fibre, H = {H}, ellrank bounds [{r1},{r2}]. Written {datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%MZ')}\n")
    f.write("E " + " ".join(a) + "\n")
    for x, y in P:
        f.write(f"P {x.strip()} {y.strip()}\n")
    f.write(f"#CLAIM rank={r1}\n")
print(sys.argv[2], 'points', len(P))
