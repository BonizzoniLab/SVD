import argparse

parser = argparse.ArgumentParser('Parse VCF output from Platypus to output valid VCF.  Output is to stdout.')
parser.add_argument('-f', '--file', help="platypus vcf output file name")
opts = parser.parse_args()

read=open(opts.file)

for line in read:
	if (line.startswith('##FORMAT=<ID=NV') or line.startswith('##FORMAT=<ID=NR') or line.startswith('##INFO=<ID=NF,') or line.startswith('##INFO=<ID=NR,')):
		riga=(line.split(','))
		riga[riga.index('Number=.')]='Number=A'
		line=','.join(riga)
	print line.rstrip()