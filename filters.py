import argparse

def txt_reader(txt):
    ''' reads the dataset.txt and filter out variants that do not satisfy thresholds'''
    af=0
    dp=0
    
    for line in txt:
        line=line.rstrip()
        row=line.split('\t')
        if line.startswith('SAMPLE_ID'):
            af=row.index("AF_median")
            dp=row.index("DP_mean")
            print line
        else:
            if float(row[af])> float(1):
                row[af]=1
                
            if float(row[af])> float(opts.vaf) and float(row[af])<= float(1) and float(row[dp])>= float(opts.dp):
                print line  


def main():
    parser = argparse.ArgumentParser('Parse VCF output to output valid VCF. Output is to stdout.')

    parser.add_argument('-i', '--file', help="dataset tab delimited")
    parser.add_argument('-f', '--vaf', help="allele frequency threshold, decimal")
    parser.add_argument('-d', '--dp', help="depth of coverage threshold")
    
    global opts
    opts = parser.parse_args()

    in_file = open(opts.file)
    txt_reader(in_file)
    
main()













