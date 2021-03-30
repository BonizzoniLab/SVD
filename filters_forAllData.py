import argparse

def txt_reader(txt):
    ''' read the data set and variants that satisfy the following condition are maintained: 
    [ abs(length(REF) - length(ALT)) >= 6 & [AF_median > 0.3 OR Num_caller > 1 ] & STRBIAS_median < MaxStrandBias ] for INDELs and 
    [AF_median > 0.3 OR Num_caller > 1 ] & STRBIAS_median < MaxStrandBias] for SNPs'''
    af=0
    strb=0
    ref=''
    alt=''
    fb=0
    vd=0
    pl=0
    gk=0
    
    for line in txt:
        line=line.rstrip()
        row=line.split('\t')
        if line.startswith('SAMPLE_ID'):
            af=row.index("AF_median")
            strb=row.index("STRBIAS_median")
            ref=row.index("REF")
            alt=row.index("ALT")
            fb=row.index("CallFreebayes")
            vd=row.index("CallVardict")
            pl=row.index("CallPlatypus")
            gk=row.index("CallGatk")
            print line
            
        else:
            if row[strb]=='.':
                row[strb]=0
                
            if (opts.variantType=='INDEL' and float(row[strb]) <= float(opts.MaxStrBias) and abs(len(row[ref])-len(row[alt])) >= int(opts.Nt) and ((float(row[af]) > float(opts.MinAf)) or (int(row[fb])+int(row[vd])+int(row[pl])+int(row[gk]) >= int(opts.minCaller)))) :
                print line  
            if (opts.variantType=='SNP' and float(row[strb]) <= float(opts.MaxStrBias) and ((float(row[af]) > float(opts.MinAf)) or (int(row[fb])+int(row[vd])+int(row[pl])+int(row[gk]) >= int(opts.minCaller)))) :
                print line

def main():
    parser = argparse.ArgumentParser('Parse VCF output to output valid VCF.  Output is to stdout.')

    parser.add_argument('-i', '--file', help="dataset tab delimited")
    parser.add_argument('-f', '--MinAf', help="min allele frequency required to call a snp or an indel in AllData")
    parser.add_argument('-sb', '--MaxStrBias', help="max strand bias accepted to call a snp or an indel in AllData")
    parser.add_argument('-nt', '--Nt', help="min number of nucleotides required to call a snp or an indel in AllData", default=6)
    parser.add_argument('-calls', '--minCaller', help="min number of callers required to define a snp or an indel in AllData")
    parser.add_argument('-vt', '--variantType', help="Type of the variant analyzed")
    
    global opts
    opts = parser.parse_args()

    in_file = open(opts.file)
    txt_reader(in_file)
    
main()













