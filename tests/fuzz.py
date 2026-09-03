#!/usr/bin/env python3
"""Randomised differential test aimed at the vectorised nucleotide diagonal test.

Varies everything the kernel's correctness depends on: sequence lengths (which set the
diagonal span and the number of 64-lane tiles), N and masked-base placement (the 0xFF /
0x00 sentinel split), homopolymer runs (which drive the complexity weight to 1), the
identity threshold and word length (which set band_b, hence the span), and -b (the band
width itself). Any divergence from the C++ is a kernel bug.
"""
import random, subprocess, sys, os, tempfile

CPP = "./cdhit/cd-hit-est"
RS  = sys.argv[1] if len(sys.argv) > 1 else "./target/release/cd-hit-est"
N   = int(sys.argv[2]) if len(sys.argv) > 2 else 60

def make(rng, path):
    nseq = rng.randint(4, 60)
    lo   = rng.choice([30, 60, 120, 400, 900, 2000])
    hi   = lo + rng.choice([0, 20, 200, 1500])
    fams = max(1, nseq // rng.randint(2, 6))
    bases = ['ACGT'] * 7 + ['ACGTN', 'ACGTNRYKM']
    alpha = rng.choice(bases)
    recs, base = [], None
    for i in range(nseq):
        if i % fams == 0 or base is None:
            L = rng.randint(lo, hi)
            base = ''.join(rng.choice(alpha) for _ in range(L))
            if rng.random() < 0.25:                       # homopolymer run -> cpx == 1
                p = rng.randrange(max(1, len(base) - 40))
                base = base[:p] + rng.choice('ACGT') * rng.randint(8, 40) + base[p + 40:]
        s = list(base)
        for _ in range(int(len(s) * rng.uniform(0, 0.18))):
            p = rng.randrange(len(s)); s[p] = rng.choice(alpha)
        s = ''.join(s)
        if rng.random() < 0.3:                            # ragged ends -> shifts the band
            s = s[rng.randint(0, 20):len(s) - rng.randint(0, 20)]
        if rng.random() < 0.2:                            # N run -> sentinel path
            p = rng.randrange(max(1, len(s) - 12))
            s = s[:p] + 'N' * rng.randint(1, 12) + s[p + 12:]
        if len(s) >= 15:
            recs.append('>s%d\n%s' % (i, s))
    if not recs:
        return False
    open(path, 'w').write('\n'.join(recs) + '\n')
    return True

def run(binary, inp, out, opts):
    return subprocess.run([binary, '-i', inp, '-o', out, '-T', '1'] + opts,
                          capture_output=True, timeout=180)

rng = random.Random(int(os.environ.get('SEED', '1')))
tmp = tempfile.mkdtemp()
bad = 0
for t in range(N):
    inp = os.path.join(tmp, 'in.fna')
    if not make(rng, inp):
        continue
    c = rng.choice(['0.80', '0.85', '0.88', '0.90', '0.93', '0.95', '0.97', '0.99', '1.00'])
    n = rng.choice(['4', '5', '6', '7', '8', '9', '10', '11'])
    opts = ['-c', c, '-n', n]
    if rng.random() < 0.3: opts += ['-b', rng.choice(['1', '5', '10', '50', '200'])]
    if rng.random() < 0.2: opts += ['-r', '0']
    if rng.random() < 0.2: opts += ['-g', '1']
    if rng.random() < 0.2: opts += ['-p', '1']
    if rng.random() < 0.2: opts += ['-mask', 'N']
    if rng.random() < 0.15: opts += ['-aL', '0.7', '-aS', '0.7']
    rc = run(CPP, inp, os.path.join(tmp, 'c'), opts)
    rr = run(RS,  inp, os.path.join(tmp, 'r'), opts)
    ok_c, ok_r = os.path.exists(tmp + '/c.clstr'), os.path.exists(tmp + '/r.clstr')
    if ok_c != ok_r:
        print('DIVERGE(run) %s : cpp=%s rust=%s' % (' '.join(opts), ok_c, ok_r)); bad += 1
    elif ok_c and open(tmp + '/c.clstr', 'rb').read() != open(tmp + '/r.clstr', 'rb').read():
        print('DIVERGE(out) %s' % ' '.join(opts))
        os.replace(inp, '/tmp/fuzz_fail.fna'); bad += 1
        break
    for f in ('c', 'r'):
        for e in ('', '.clstr'):
            try: os.remove(tmp + '/' + f + e)
            except OSError: pass
print('%d cases, %d divergences' % (N, bad))
sys.exit(1 if bad else 0)
