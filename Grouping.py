import argparse
import collections    
import numpy as np

class LOCUS():
    locus_name=''
    Supercontig=''
    start=''
    end=''
    groupDict={}
    alleleDict={}
    pass
   
def GroupsEvaluation(allelesfile, group_locus_dict):
    ''' Attribution of a variation group to an allele'''
    
    #initializations
    delta='\x94'
    iota='\x99'

    
    #read allele file and assign a code: 0 for absense, A for alleles equal to the annotation, B:Z for allele with only a portion of the annotated sequence
    for line in allelesfile:
        line.rstrip()
        if line.startswith('Locus_name'):
            continue
        else:
            parts = line.split("\t")
            locus_name=parts[0]
            locus_start=parts[1]
            locus_stop=parts[2]
            d={}
            a={}
            
            for part in parts:
                if ':' in list(str(part)):
                    startallele=part.split(':')[0]
                    stopallele=part.split(':')[1].split(delta)[0].split(iota)[0]
                    
                    #case absence
                    if int(startallele) == 0 and int(stopallele) == 0:
                        samples=parts[parts.index(part)+1].split('|')
                        for s in samples:
                            if s is not '':
                                d[s.split('\n')[0]]='0'
                                
                    #case no variability, that is equal to annotation
                    elif startallele == locus_start and stopallele == locus_stop and delta not in list(str(part)) and iota not in list(str(part)):
                        samples=parts[parts.index(part)+1].split('|')
                        for s in samples:
                            if s is not '':
                                d[s.split('\n')[0]]='A'
                                a['A']=str(part)
                    elif d.values()==[] or str(chr(ord(max(d.values())))) == '0':
                        new_letter='B'
                        samples=parts[parts.index(part)+1].split('|')
                        for s in samples:
                            if s is not '':
                                d[s.split('\n')[0]]=new_letter
                                a[str(new_letter)]=str(part)
                    else:
                        new_letter=chr(ord(max(d.values()))+1)
                        samples=parts[parts.index(part)+1].split('|')
                        for s in samples:
                            if s is not '':
                                d[s.split('\n')[0]]=new_letter
                                a[str(new_letter)]=str(part)

            group_locus_dict[locus_name].groupDict=d
            group_locus_dict[locus_name].alleleDict=a

    return group_locus_dict


def WriteGroups(samplesName, group_locus_dict, snp_indel_pol_dict):
    ''' Write the output file '''
    
    outFile=open(opts.outputPath+'/SummaryOf_StructuralVariability_InSamples.txt','w')
    
    header='Locus_name\t'
    for s in sorted(samplesName):
        header=header+str(s)+'\t'
    header=header+'\t\tAllele_Legend'
    outFile.write(header+'\n')
    
    # for each locus
    for k,v in group_locus_dict.iteritems():
        # find several allele
        unique_str_allele = np.unique(v.groupDict.values())
        
        for sa in unique_str_allele:
            if str(sa) is not '0': # if locus is not always absent
                sample_with_same_str_allele=[]
                count=1
                
                # find samples supporting the same str_allele
                for sample, letter in v.groupDict.iteritems():
                    if letter == sa:
                        sample_with_same_str_allele.append(sample.split('\n')[0])
                        
                for pol, samples in snp_indel_pol_dict[k].iteritems():

                    # in loci with polymorphisms but in which samples have polymorphisms sequence like '000..00' the pedix is 0
                    sample_same_str_pol_allele=[]
                    for s_str_all in sample_with_same_str_allele:
                        if s_str_all in samples:
                            sample_same_str_pol_allele.append(s_str_all)
                    
                    if sample_same_str_pol_allele:    
                        if str(1) in list(pol):        
                            for s in sample_same_str_pol_allele:
                                old_letter=v.groupDict[s]
                                new_letter=old_letter+str(count)
                                v.groupDict[s]=new_letter
                            count=count+1
                            
                        else:
                            for s in sample_same_str_pol_allele:
                                old_letter=v.groupDict[s]
                                new_letter=old_letter+str(0)
                                v.groupDict[s]=new_letter
                                               
        string_to_print=''
        od_sample = collections.OrderedDict(sorted(v.groupDict.items()))
        for s in od_sample:
            string_to_print=string_to_print+str(od_sample[s])+'\t'
        
        string_to_print=string_to_print+'\t\t'
        
        od_allele = collections.OrderedDict(sorted(v.alleleDict.items()))
        for a in od_allele:
            string_to_print=string_to_print+str(a)+':'+str(od_allele[a])+'\t'
            
        outFile.write(k+'\t'+string_to_print+'\n')
        
    outFile.close()   


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


def read_variants(in_file):
    '''method to read the indel or snp file, save variants in dictionary'''
    
    variants_dictionary={}
    
    with in_file as f:
        lis = [x.split() for x in f]
    first_line=lis[0]
    sample_list = first_line[5:-2]
    
    j=0
    for x in zip(*lis):
        i=0
        row=[]
        
        for y in x:
            if (i==0):
                i=i+1
                continue
            
            else:
                if not y in first_line:
                    row.append(y)
                    i=i+1   
            
        variants_dictionary[first_line[j]]=row
        j=j+1
    
    in_file.close()    
    return sample_list, variants_dictionary;


def load_variants(samples_list, group_locus_dict, variants_dictionary):
    '''method to save variants in dictionary'''
    snp_indel_polimorphisms={}
    
    all_loci_var = np.array(variants_dictionary['LOCUS_NAME'])
    for locus in group_locus_dict.iterkeys():

        indicies_locus = np.where(locus == all_loci_var)
        d={}
        
        for sample in samples_list:
            all_var_sample = np.array(variants_dictionary[sample])
            locus_var_sample = all_var_sample[indicies_locus]
            
            if str(group_locus_dict[locus].groupDict[sample]) == str(0) and str(1) in locus_var_sample:
                print('Variants in 0 locus: '+str(locus)+' '+str(sample))
                print(locus_var_sample)
                locus_var_sample[np.where(locus_var_sample == '1')]='0'
                print(locus_var_sample)
                
            single_pol = ''.join(locus_var_sample)
            
            if not single_pol == '': # so there are not polimorphisms for this locus called by variant calling
                if single_pol in d.keys():
                    temp=d[single_pol]
                    temp.append(sample)
                    d[single_pol]=temp
                else:
                    d[single_pol]=[sample]
            
        snp_indel_polimorphisms[locus]=d

    return snp_indel_polimorphisms
    
    
def main():
    #read parameters of the script
    parser = argparse.ArgumentParser('Locus alleles grouped by the variations founded (INDEL, cropped left, cropped rigth and cropped central). Output is to stdout.')
    parser.add_argument('-RegionsInfo', '--RegionsInfo', help="input file containing regions locus name, viral family, Supercontig, Start and End positions tab delimited")
    parser.add_argument('-Variants', '--VariantsList', help="file of the variants")
    parser.add_argument('-AlleleFile', '--AllelesFile', help="input file containing the alleles found for each locus of interest")
    parser.add_argument('-outputPath', '--outputPath', help="Path of the output file")
    
    global opts 
    opts = parser.parse_args()
    
    # open file in input
    VarFileList=open(opts.VariantsList)
    alleles=open(opts.AllelesFile)
    in_info=open(opts.RegionsInfo)
    
    # read locus info and variants and load them in structure
    group_locus_dict=read_locus_info(in_info)
    
    # read structural different alleles and load them in structure
    group_locus_dict=GroupsEvaluation(alleles, group_locus_dict)
    
    # read variants and load them in structure
    [samplesName, variants_dictionary] = read_variants(VarFileList)
    snp_indel_pol_dict=load_variants(samplesName, group_locus_dict, variants_dictionary)
    
    # write the output file with the structural variability group
    WriteGroups(samplesName, group_locus_dict, snp_indel_pol_dict)
    
    # close file in input
    alleles.close()
    VarFileList.close()
    in_info.close()
    
main()