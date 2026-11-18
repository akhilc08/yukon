# bcgen.py FAMFILE N HB n0 d0 n1 d1 n2 d2 :
# fibres of the base change T = (n0 s^2 + n1 s + n2)/(d0 s^2 + d1 s + d2), s = p/q, |p|,q <= N.
# Prints "P Q" (T = P/Q reduced) whose cheap lower bound on the minimal-model height (both 2-isogenous
# models, H_int - log(gcd(c4^3, c6^2))/12) is below HB.
import sys, math
fam = {}
for line in open(sys.argv[1]):
    w = line.split()
    if w and w[0] in ('A', 'B'): fam[w[0]] = [int(x) for x in w[1:]]
A, B = fam['A'], fam['B']
N = int(sys.argv[2]); HB = float(sys.argv[3])
n0, d0, n1, d1, n2, d2 = map(int, sys.argv[4:10])
def hlow(c4, c6):
    if c4 == 0 and c6 == 0: return 1e9
    m = max(abs(c4) ** 3, c6 * c6); g = math.gcd(c4 ** 3, c6 * c6)
    return (math.log(m) - math.log(g)) / 12 if g else math.log(m) / 12
seen = set()
for q in range(1, N + 1):
    for p in range(-N, N + 1):
        if math.gcd(p, q) != 1: continue
        P = n0 * p * p + n1 * p * q + n2 * q * q; Q = d0 * p * p + d1 * p * q + d2 * q * q
        if Q == 0: continue
        g = math.gcd(P, Q); P //= g; Q //= g
        if Q < 0: P, Q = -P, -Q
        if (P, Q) in seen: continue
        seen.add((P, Q))
        a = sum(A[i] * P ** i * Q ** (4 - i) for i in range(5)); b = sum(B[i] * P ** i * Q ** (8 - i) for i in range(9))
        if b == 0 or a * a - 4 * b == 0: continue
        h0 = hlow(16 * (a * a - 3 * b), -64 * a ** 3 + 288 * a * b)
        h1 = hlow(16 * (a * a + 12 * b), -64 * a ** 3 + 2304 * a * b)
        if min(h0, h1) < HB: print(P, Q)
