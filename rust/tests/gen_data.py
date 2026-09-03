#!/usr/bin/env python3
"""Generate the synthetic FASTA/FASTQ inputs used by the differential tests."""
import random, os, sys, gzip

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
os.makedirs(OUT, exist_ok=True)
DNA = 'ACGT'
AA = 'ACDEFGHIKLMNPQRSTVWY'

def mutate(s, rate, alpha, rng, indel=0.0):
    out = []
    for c in s:
        r = rng.random()
        if r < rate:
            out.append(rng.choice(alpha))
        elif r < rate + indel:
            if rng.random() < 0.5:
                continue                      # deletion
            out.append(c); out.append(rng.choice(alpha))   # insertion
        else:
            out.append(c)
    return ''.join(out)

def revcomp(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]

def wrap(s, w=60):
    return '\n'.join(s[i:i+w] for i in range(0, len(s), w))

def write(name, recs, wrap_at=60, gz=False):
    path = os.path.join(OUT, name)
    text = ''.join('>%s\n%s\n' % (h, wrap(s, wrap_at) if wrap_at else s) for h, s in recs)
    if gz:
        with gzip.open(path, 'wt') as f: f.write(text)
    else:
        with open(path, 'w') as f: f.write(text)
    return path

# ---- dna_small: a few families, point mutations only, uniform length
rng = random.Random(11)
recs = []
for fam in range(12):
    base = ''.join(rng.choice(DNA) for _ in range(500))
    for i in range(20):
        recs.append(('d%d_%d some description here' % (fam, i), mutate(base, rng.uniform(0, .12), DNA, rng)))
rng.shuffle(recs)
write('dna_small.fna', recs)

# ---- dna_mixed: varied lengths, indels, N runs, reverse-complemented members
rng = random.Random(23)
recs = []
for fam in range(20):
    L = rng.randint(120, 1500)
    base = ''.join(rng.choice(DNA) for _ in range(L))
    for i in range(rng.randint(3, 15)):
        s = mutate(base, rng.uniform(0, .10), DNA, rng, indel=0.004)
        if rng.random() < 0.35:
            s = revcomp(s)
        if rng.random() < 0.25:                       # inject an N run
            p = rng.randrange(max(1, len(s) - 20))
            s = s[:p] + 'N' * rng.randint(1, 12) + s[p + 12:]
        if rng.random() < 0.3:                        # ragged ends
            s = s[rng.randint(0, 30):len(s) - rng.randint(0, 30)]
        recs.append(('m%d_%d' % (fam, i), s))
for i in range(40):                                   # singletons
    recs.append(('s%d' % i, ''.join(rng.choice(DNA) for _ in range(rng.randint(60, 900)))))
rng.shuffle(recs)
write('dna_mixed.fna', recs, wrap_at=70)
write('dna_mixed.fna.gz', recs, wrap_at=70, gz=True)

# ---- prot_small / prot_mixed
rng = random.Random(37)
recs = []
for fam in range(15):
    base = ''.join(rng.choice(AA) for _ in range(300))
    for i in range(rng.randint(2, 14)):
        recs.append(('p%d_%d protein description' % (fam, i),
                     mutate(base, rng.uniform(0, .35), AA, rng)))
rng.shuffle(recs)
write('prot_small.faa', recs)

rng = random.Random(41)
recs = []
for fam in range(25):
    L = rng.randint(50, 900)
    base = ''.join(rng.choice(AA) for _ in range(L))
    for i in range(rng.randint(2, 12)):
        s = mutate(base, rng.uniform(0, .45), AA, rng, indel=0.006)
        if rng.random() < 0.3:
            s = s[rng.randint(0, 20):len(s) - rng.randint(0, 20)]
        if rng.random() < 0.15:
            s = s + '*'                                # trailing stop codon
        if rng.random() < 0.1:
            p = rng.randrange(max(1, len(s) - 5))
            s = s[:p] + 'X' * rng.randint(1, 5) + s[p + 5:]
        recs.append(('q%d_%d' % (fam, i), s))
for i in range(60):
    recs.append(('t%d' % i, ''.join(rng.choice(AA) for _ in range(rng.randint(20, 600)))))
rng.shuffle(recs)
write('prot_mixed.faa', recs, wrap_at=80)

# ---- dna_dup: heavy exact duplication, short reads (the cd-hit-est common case)
rng = random.Random(53)
recs = []
pool = [''.join(rng.choice(DNA) for _ in range(150)) for _ in range(30)]
for i in range(600):
    s = rng.choice(pool)
    if rng.random() < 0.5:
        s = mutate(s, rng.uniform(0, .06), DNA, rng)
    recs.append(('r%d' % i, s))
write('dna_dup.fna', recs, wrap_at=0)

# ---- FASTQ + paired end
rng = random.Random(67)
r1, r2 = [], []
pool = [''.join(rng.choice(DNA) for _ in range(200)) for _ in range(25)]
for i in range(300):
    s = mutate(rng.choice(pool), rng.uniform(0, .05), DNA, rng)
    r1.append(('e%d' % i, s[:100]))
    r2.append(('e%d' % i, revcomp(s[-100:])))
def write_fq(name, recs):
    with open(os.path.join(OUT, name), 'w') as f:
        for h, s in recs:
            f.write('@%s\n%s\n+\n%s\n' % (h, s, 'I' * len(s)))
write_fq('pe_R1.fq', r1)
write_fq('pe_R2.fq', r2)

print('wrote', sorted(os.listdir(OUT)))
