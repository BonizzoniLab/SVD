import argparse



parser = argparse.ArgumentParser('Parse VCF output from Gatk to output valid VCF.  Output is to stdout.')
parser.add_argument('-f', '--file', help="gatk vcf output file name")
opts = parser.parse_args()

read=open(opts.file)

for line in read:
    if line.startswith('##FORMAT=<ID=AD'):
        riga=(line.split(','))
        riga[riga.index('Number=.')]='Number=G'
        line=','.join(riga)
    print line.rstrip()
    