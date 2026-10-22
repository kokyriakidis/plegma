import random, sys
L, N, rate, seed = int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4])
rng = random.Random(seed)
ref = "".join(rng.choice("ACGT") for _ in range(L))
open(f"ref_{L}.fa", "w").write(f">ref\n{ref}\n")
with open(f"reads_{L}_{N}_{rate}.fa", "w") as f:
    for i in range(N):
        out = []
        for c in ref:
            r = rng.random()
            if r < rate / 3:
                out.append(rng.choice("ACGT".replace(c, "")))
            elif r < 2 * rate / 3:
                continue
            elif r < rate:
                out.append(c); out.append(rng.choice("ACGT"))
            else:
                out.append(c)
        f.write(">r%d\n%s\n" % (i, "".join(out)))
