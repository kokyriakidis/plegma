import random, subprocess, sys
binary = sys.argv[1]
def ed(a, b):
    prev = list(range(len(b) + 1))
    for i in range(1, len(a) + 1):
        cur = [i] + [0] * len(b)
        for j in range(1, len(b) + 1):
            cur[j] = min(prev[j - 1] + (a[i-1] != b[j-1]), prev[j] + 1, cur[j-1] + 1)
        prev = cur
    return prev[-1]
def mums(a, b, L):
    out = []
    for i in range(len(a)):
        for j in range(len(b)):
            if i and j and a[i-1] == b[j-1]: continue
            l = 0
            while i+l < len(a) and j+l < len(b) and a[i+l] == b[j+l]: l += 1
            if l < L: continue
            s = a[i:i+l]
            if a.count(s) == 1 and b.count(s) == 1 and (a.find(s) == i):
                # count overlapping occurrences properly
                ca = sum(1 for k in range(len(a)-l+1) if a[k:k+l] == s); cb = sum(1 for k in range(len(b)-l+1) if b[k:k+l] == s)
                if ca == 1 and cb == 1: out.append((i, j, l))
    return out
rng = random.Random(3)
bad = 0
for t in range(40):
    u = ''.join(rng.choice('ACGT') for _ in range(rng.randint(5, 12)))
    def arr(n):
        s = ''
        for _ in range(n):
            x = list(u)
            if rng.random() < .4: x[rng.randrange(len(x))] = rng.choice('ACGT')
            s += ''.join(x)
        return s
    L0 = ''.join(rng.choice('ACGT') for _ in range(30)); R0 = ''.join(rng.choice('ACGT') for _ in range(30))
    a = L0 + arr(rng.randint(5, 15)) + R0; b = L0 + arr(rng.randint(5, 15)) + R0
    line = 't\t8\t30\t%d\t30\t%d\t1\t%s\t%s\n' % (len(a) - 30, len(b) - 30, a, b)
    out = subprocess.run([binary], input=line, capture_output=True, text=True).stdout.splitlines()
    h = out[0].split('\t'); v = dict(zip(h, out[1].split('\t')))
    M = mums(a, b, 8); opt = ed(a, b)
    zero = sum(1 for (i, j, l) in M if ed(a[:i], b[:j]) + ed(a[i+l:], b[j+l:]) == opt)
    got = (int(v['opt']), int(v['n_mum']), int(v['zero_core']) + int(v['zero_flank']) + int(v['zero_cross']))
    exp = (opt, len(M), zero)
    if got != exp: bad += 1; print('MISMATCH', got, exp)
print('tests', 40, 'bad', bad)
