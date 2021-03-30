import argparse
import numpy as np

def Read_Evaluation():
    '''Read all files in list and evaluate the count of the variants'''
    variantDict={}
    sample_names=[]
    varType=''
    
    plist=open(opts.pathList)
    for line in plist:
        varType=line.split('/')[-1].split('_')[0]
        sample_names.append(line.split('/')[-3])
        
    plist.close()
    
    for varfile in open(opts.pathList):
        filename=varfile.rstrip()
        in_file=open(filename)
        
        for line in in_file:
            line=line.rstrip()
            row=line.split('\t')
            dictValue=list(np.zeros((len(sample_names),), dtype=np.int))
            
            if line.startswith('SAMPLE_ID'):
                samplename=row.index("SAMPLE_ID")
                contig=row.index("CONTIG")
                pos=row.index("POS")
                ref=row.index("REF")
                alt=row.index("ALT")
                
            else:
                ID=row[contig]+':'+row[pos]+':'+row[ref]+':'+row[alt]
                name=row[samplename]
                dictValue[sample_names.index(name)]=1
                
                if ID not in variantDict.keys():
                    variantDict[ID]=dictValue
                else:
                    val=variantDict[ID]
                    val[sample_names.index(name)]=1
                    variantDict[ID]=val
                        
        in_file.close()
    return sample_names,variantDict,varType

def WriteOutput(sample_names,varType,varDict):
    '''Write the output file including the frequency of the presence of a variant between samples analyzed'''
    names=''
    for s in sample_names:
        names=names+str(s)+'\t'
        
    outFile=open(opts.outputPath+'/'+str(varType)+'Statistics.txt','w')
    outFile.write('CONTIG\tPOS\tREF\tALT\t'+names+'COUNT'+'\t'+'FREQUENCY'+'\n')
    for k,v in varDict.iteritems():
        var=''
        idv=''
        for val in v:
            var=var+str(val)+'\t'
        for el in str(k).split(':'):
            idv=idv+str(el)+'\t'
        frequencyV=float(float(sum(v))/float(len(v)))
        sumV=sum(v)
        outFile.write(idv+var+str(sumV)+'\t'+str(frequencyV)+'\n')
        
    outFile.close()   
    
    
def main():
    #read parameteres of the script
    parser = argparse.ArgumentParser('Count the frequency of the variants of the run samples. Output is to stdout.')
    parser.add_argument('-pathList', '--pathList', help="paths of the variant files in analysis")
    parser.add_argument('-outputPath', '--outputPath', help="Path of the output file")
    
    global opts 
    opts = parser.parse_args()
    [sample_names,variantDict,varType]=Read_Evaluation()
    WriteOutput(sample_names,varType,variantDict)
    
main()