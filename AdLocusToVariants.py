import argparse

class LOCUS():
    NIRV_locus_name=''
    VirusAndGroup=''
    Viral_family=''
    Supercontig=''
    start=''
    end=''
    pass

def read_locus_info(fileinfo):
    '''method to read the RegionsInfo file containing:Supercontig     start    end    locus_name'''
    dictionary={}
    
    for line in fileinfo:
        line.rstrip()
        if line.startswith('locus_name'):
            continue
        else:
            parts = line.split("\t")
            singleLocus=LOCUS()
            singleLocus.Supercontig=parts[0].strip()
            singleLocus.start=parts[1].strip()
            singleLocus.end=parts[2].strip()
            singleLocus.locus_name=parts[3].split('\n')[0].strip()
                    
            dictionary[singleLocus.locus_name]=singleLocus        
    
    return dictionary;


def read_write_variants(vartype, locus_info_dict, in_file):
    '''method to read the indel or snp file, save variants in dictionary'''
    # read sample files of INDEL variants and select only deletions and insertions related to the specific locus selected by infovalue
    
    out_file=open(opts.outputPath+'/'+vartype+'AssignedToLocusStatistic.txt','w')
       
    for line in in_file:
        line=line.rstrip()
        row=line.split('\t')
        if line.startswith('CONTIG'):
            contig=row.index("CONTIG")
            pos=row.index("POS")
            out_file.write('LOCUS_NAME\t'+line+'\n')
        else:
            for locus, infovalue in locus_info_dict.iteritems():
                if ((str(row[contig])==infovalue.Supercontig) and (int(row[pos]) in range(int(infovalue.start),int(infovalue.end)+1))):
                    out_file.write(locus+'\t'+line+'\n')
                    
    in_file.close()
    out_file.close()
                
def main():
    #read parameteres of the script
    parser = argparse.ArgumentParser('Add Locus name to each variant found. Output is to stdout.')
    parser.add_argument('-RegionsInfo', '--RegionsInfo', help="input file containing regions locus name, viral family, Supercontig, Start and End positions tab delimited")
    parser.add_argument('-IndelVariants', '--INDELlist', help="list of the INDEL variant file")
    parser.add_argument('-SnpVariants', '--SNPlist', help="list of the SNP variant file")
    parser.add_argument('-outputPath', '--outputPath', help="Path of the output file")
    
    global opts 
    opts = parser.parse_args()
    
    #paths list of the variant file
    INDELFileList=open(opts.INDELlist)
    SNPFileList=open(opts.SNPlist)
    
    #read locus info file and gatk dp file. Save both data in dictionary structure
    in_info=open(opts.RegionsInfo)
    locus_info_dict=read_locus_info(in_info)
    read_write_variants('indel', locus_info_dict, INDELFileList)
    read_write_variants('snp', locus_info_dict, SNPFileList)
    
    INDELFileList.close()
    SNPFileList.close()
    
main()