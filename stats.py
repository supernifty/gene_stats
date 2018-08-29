#!/usr/bin/env python
'''
  count coding sequence
'''

import argparse
import collections
import logging
import sys

def main(allowed, padding, use_max):
  logging.info('starting...')

  #bin    name    chrom   strand  txStart txEnd   cdsStart        cdsEnd  exonCount       exonStarts      exonEnds        score   name2   cdsStartStat    cdsEndStat      exonFrames
    #0       NM_001308203.1  chr1    +       66999251        67216822        67000041        67208778        22      66999251,66999928,67091529,67098752,67105459,67108492,67109226,67136677,67137626,67138963,67142686,67145360,67154830,67155872,67160121,67184976,67194946,67199430,67205017,67206340,67206954,67208755,  66999355,67000051,67091593,67098777,67105516,67108547,67109402,67136702,67137678,67139049,67142779,67145435,67154958,67155999,67160187,67185088,67195102,67199563,67205220,67206405,67207119,67216822,  0       SGIP1   cmpl    cmpl    -1,0,1,2,0,0,1,0,1,2,1,1,1,0,1,1,2,2,0,2,1,1,
    #1       NM_001323574.1  chr1    -       48998526        50489626        48999844        50489468        14      48998526,49000561,49005313,49052675,49056504,49100164,49119008,49128823,49332862,49511255,49711441,50162948,50317067,50489434,  48999965,49000588,49005410,49052838,49056657,49100276,49119123,49128913,49332902,49511472,49711536,50163109,50317190,50489626,  0       AGBL4   cmpl    cmpl    2,2,1,0,0,2,1,1,0,2,0,1,1,0,

  genes = collections.defaultdict(int)
  if allowed is not None:
    allowed_set = set(allowed)

  for line in sys.stdin:
    if line.startswith('#'):
      continue

    fields = line.strip('\n').split('\t')
    gene = fields[12]
    codingStart = int(fields[6])
    codingEnd = int(fields[7])

    if allowed is not None and gene not in allowed_set:
      continue

    if not use_max and gene in genes:
      continue

    transcript = 0
    for start, finish in zip(fields[9].split(','), fields[10].split(',')):
      if start == '' or finish == '':
        continue
      start = max(codingStart, int(start))
      finish = min(codingEnd, int(finish))
      if finish - start > 0:
        transcript += finish - start + padding * 2
    genes[gene] = max(genes[gene], transcript)

  if allowed is None:
    for gene in sorted(genes):
      sys.stdout.write('{}\t{}\n'.format(gene, genes[gene]))
  else:
    for gene in allowed:
      if gene in genes:
        sys.stdout.write('{}\t{}\n'.format(gene, genes[gene]))
      else:
        sys.stdout.write('{}\t{}\n'.format(gene, 'MISSING'))

  logging.info('done')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Count coding sequence')
  parser.add_argument('--genes', nargs='+', help='gene filter')
  parser.add_argument('--padding', default=0, type=int, help='exon padding')
  parser.add_argument('--max', action='store_true', help='max transcript size')
  parser.add_argument('--verbose', action='store_true', help='more logging')
  args = parser.parse_args()
  if args.verbose:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG)
  else:
    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.INFO)

  main(args.genes, args.padding, args.max)
