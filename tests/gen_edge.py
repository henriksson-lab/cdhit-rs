#!/usr/bin/env python3
"""Awkward inputs: the FASTA reader is where a translation is most likely to drift."""
import os, random
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
DNA, AA = 'ACGT', 'ACDEFGHIKLMNPQRSTVWY'
rng = random.Random(97)
def w(name, text):
    open(os.path.join(OUT, name), 'w', newline='').write(text)

seq = lambda n, a=DNA: ''.join(rng.choice(a) for _ in range(n))

w('edge_empty.fna', '')
w('edge_one.fna', '>only\n%s\n' % seq(300))
w('edge_no_trailing_nl.fna', '>a\n%s\n>b\n%s' % (seq(300), seq(300)))
w('edge_crlf.fna', '>a\r\n%s\r\n>b\r\n%s\r\n' % (seq(300), seq(300)))
w('edge_blank_lines.fna', '>a\n%s\n\n\n>b\n%s\n\n' % (seq(300), seq(300)))
w('edge_lowercase.fna', '>a\n%s\n>b\n%s\n' % (seq(300).lower(), seq(300).lower()))
w('edge_dup_names.fna', ''.join('>same description\n%s\n' % seq(200) for _ in range(20)))
w('edge_long_header.fna',
  ''.join('>%s|%s\n%s\n' % ('h' * 400, i, seq(200)) for i in range(10)))
w('edge_no_space_header.fna', ''.join('>tok%d\n%s\n' % (i, seq(200)) for i in range(10)))
w('edge_tab_header.fna', ''.join('>id%d\tdesc here\n%s\n' % (i, seq(200)) for i in range(10)))
w('edge_space_after_gt.fna', ''.join('> id%d desc\n%s\n' % (i, seq(200)) for i in range(10)))
w('edge_short.fna', ''.join('>s%d\n%s\n' % (i, seq(rng.randint(1, 15))) for i in range(40)))
w('edge_one_per_line.fna', ''.join('>l%d\n%s\n' % (i, '\n'.join(seq(120))) for i in range(6)))
w('edge_iupac.fna',
  ''.join('>u%d\n%s\n' % (i, ''.join(rng.choice('ACGTNRYKMSWBDHV')
                                     for _ in range(300))) for i in range(30)))
w('edge_bad_chars.fna', '>ok\n%s\n>bad\n%s\n>ok2\n%s\n' % (seq(200), 'ACGT123!@#' * 20, seq(200)))
w('edge_empty_seq.fna', '>a\n\n>b\n%s\n>c\n\n' % seq(300))
w('edge_stops.faa', ''.join('>p%d\n%s*\n' % (i, seq(150, AA)) for i in range(20)))
w('edge_identical.fna', ''.join('>i%d\n%s\n' % (i, 'ACGTACGTAC' * 30) for i in range(50)))
w('edge_homopolymer.fna', ''.join('>h%d\n%s\n' % (i, 'A' * (200 + i)) for i in range(20)))
w('edge_all_n.fna', ''.join('>n%d\n%s\n' % (i, 'N' * 200) for i in range(10)))
print('ok')
