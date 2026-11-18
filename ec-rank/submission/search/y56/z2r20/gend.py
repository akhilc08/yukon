# generate primitive quadratic forms D = a p^2 + b p q + c q^2 with |a|,|b|,|c| <= N, scaled: T = U/S
# (D(T) in U-coordinates has coefficients a, b*S, c*S^2 relative to (p, q)) ; output a b c in T-coordinates
import sys, math
N = int(sys.argv[1]); S = int(sys.argv[2])
out = sys.stdout
for a in range(0, N + 1):
    for b in range(-N, N + 1):
        for c in range(-N, N + 1):
            if a == 0 and (b < 0 or (b == 0)): continue
            if math.gcd(math.gcd(a, b), c) != 1: continue
            d = b * b - 4 * a * c
            if d >= 0 and math.isqrt(d) ** 2 == d and a != 0:
                pass  # split forms are fine (two rational ramification points)
            # D(T) with T = U/S: a U^2/S^2 + b U/S + c  ~  a U^2 + b S U + c S^2  (U = S T)
            # in T: a S^2 T^2 + b S T + c  -> coefficients (a*S*S, b*S, c)
            out.write(f"{a*S*S} {b*S} {c}\n")
