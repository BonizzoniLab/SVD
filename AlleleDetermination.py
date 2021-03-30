import argparse
import numpy as np
import collections

class Variant():
    ID='' #ID=CONTIG:POS:REF:ALT
    SampleID=''
    locus_name=''
    CONTIG=''
    POS=''
    REF=''
    ALT=''
    Call=[] #Call=[CallFB, CallVD, CallPL, CallGK]
    GT=[] #GT=[GTFB, GTVD, GTPL, GTGK]
    AF_median=0
    STRBIAS_median=0
    varType=''
    pass

class LOCUS():
    locus_name=''
    VirusAndGroup=''
    Viral_family=''
    Supercontig=''
    start=''
    end=''
    indices=[]
    DpLocusDict={}
    AllDataDict={}
    pass

class AllDataInfo():
    locus_name=''
    SampleName=''
    StartAllele=0
    StopAllele=0
    delta=u"\u0394"
    iota=u"\u0399"
    StartDeletions=[]
    StopDeletions=[]
    StartInsertions=[]
    StopInsertions=[]
    Genotypes=[]
    pass

def WriteAlleles(samplelist,locus_info_dict):
    '''write the Alleles file'''
    out_Alleles = open( opts.outputPath+'/Alleles.txt' ,"w")
    out_Alleles.write('Locus_name\tLocus_Start\tLocus_Stop\tAlleles_Name_Count\n')
    
    for locus,infovalue in locus_info_dict.iteritems():
        starts=[]
        stops=[]
        samplename=[]
        num00=''
        
        lenAnnotation=int(infovalue.end)-int(infovalue.start)
        thrSim=int(float(opts.thresholdSimilarity)*float(lenAnnotation))
        writeline=str(locus)+'\t'+str(int(infovalue.start)+1)+'\t'+str(infovalue.end)+'\t'
        od = collections.OrderedDict(sorted(infovalue.AllDataDict.items()))
        
        for sample, alldatainfo in od.iteritems():
            if (int(alldatainfo.StopAllele) - int(alldatainfo.StartAllele)) >= int(opts.minLengthAllele) and int(alldatainfo.StartAllele) != 0 and int(alldatainfo.StopAllele !=0):
                starts.append(alldatainfo.StartAllele)
                stops.append(alldatainfo.StopAllele)
                samplename.append(alldatainfo.SampleName)
            
            else:
                num00=num00+str(sample)+'|'
                
        if num00 !='':
            writeline=writeline+'0:0\t'+str(num00)+'\t'
            
        Sortedstart=sorted(starts)
        indexstart=sorted(range(len(starts)), key=lambda k: starts[k])
        Sortedstop=[]
        Sortedsamplename=[]
        
        for i in indexstart:
            Sortedstop.append(stops[i])
            Sortedsamplename.append(samplename[i])
            
        tempStart=Sortedstart
        tempStop=Sortedstop
        tempName=Sortedsamplename
        selStop=[]
        selName=[]
        minStartTemp=[]
        maxStartTemp=[]
        minStopTemp=[]
        maxStopTemp=[]
        
        while tempStart:
            uniqTempStart=np.unique(tempStart)
            minStartTemp=min(tempStart)
            maxStartTemp=min(tempStart)+thrSim
            selectedStop=[]
            selectedName=[]
            OtherStarts=[]
            OtherStop=[]
            OtherName=[]
            
            # look for the starts similarity
            for s in uniqTempStart:
                if s in range(minStartTemp, maxStartTemp+1):
                    indices = [i for i, x in enumerate(tempStart) if x == s]
                    for i in indices:
                        selectedStop.append(tempStop[i])
                        selectedName.append(tempName[i])
                        
                else:
                    indices = [i for i, x in enumerate(tempStart) if x == s]
                    for i in indices:
                        OtherStarts.append(s)
                        OtherStop.append(tempStop[i])
                        OtherName.append(tempName[i])
                        
            selStop=selectedStop
            selName=selectedName
            
            while selStop:
                
                uniqTempStop=np.unique(selStop)
                minStopTemp=min(selStop)
                maxStopTemp=min(selStop)+thrSim
                equalStop=[]
                equalName=[]
                NotEqualStop=[]
                NotEqualName=[]
                
                if maxStopTemp > int(infovalue.end):
                    maxStopTemp=int(infovalue.end)
                    
                # look for the stops similarity
                for st in uniqTempStop:
                    if st in range(minStopTemp, maxStopTemp+1):
                        indices = [i for i, x in enumerate(selStop) if x == st]
                        for i in indices:
                            equalStop.append(st)
                            equalName.append(selName[i])
                    else:
                        indices = [i for i, x in enumerate(selStop) if x == st]
                        for i in indices:
                            NotEqualStop.append(st)
                            NotEqualName.append(selName[i])
                
                tempdict={}
                for name in equalName:
                    deletions=''
                    insertions=''
                    
                    if infovalue.AllDataDict[name].StartDeletions:
                        for d in infovalue.AllDataDict[name].StartDeletions:
                            deletions=deletions+str(infovalue.AllDataDict[name].delta.encode('utf-8'))+str(d)+':'+str(infovalue.AllDataDict[name].StopDeletions[infovalue.AllDataDict[name].StartDeletions.index(d)])
                    
                    if infovalue.AllDataDict[name].StartInsertions:
                        for i in infovalue.AllDataDict[name].StartInsertions:
                            insertions=insertions+str(infovalue.AllDataDict[name].iota.encode('utf-8'))+str(i)+':'+str(infovalue.AllDataDict[name].StopInsertions[infovalue.AllDataDict[name].StartInsertions.index(i)])
                   
                    key=deletions+insertions
                    if key not in tempdict.keys():
                        tempdict[key]=name
                    else:
                        val=tempdict[key]+'|'+str(name)
                        tempdict[key]=val
                
                for k,v in tempdict.iteritems():
                    writeline=writeline+str(minStartTemp)+':'+str(max(equalStop))+k+'\t'+str(v)+'\t'
                        
                selStop=NotEqualStop        
                selName=NotEqualName
                
            tempStart=OtherStarts
            tempStop=OtherStop
            tempName=OtherName
        
        out_Alleles.write(writeline+'\n')
        
    out_Alleles.close()    
                
                   
def WriteAllDataFile_DetUndetFile(sampleList, locus_info_dict):
    '''write the AllDataTable file and the DetectedUndetected table file'''
    
    #open files in output
    out_AllData = open( opts.outputPath+'/AllData.txt' ,"w")
    out_DetUndet = open( opts.outputPath+'/Detected_Undetected.txt' ,"w")
    out_alwaysDet = open( opts.outputPath+'/LocusAlwaysDetected.txt' ,"w")
    out_alwaysUndet = open( opts.outputPath+'/LocusAlwaysUNdetected.txt' ,"w")
    
    #write the first line
    firstLineOut='Locus_name'
    out_alwaysDet.write(firstLineOut+"\n")
    out_alwaysUndet.write(firstLineOut+"\n")
    firstLineAllData=firstLineOut+'\tContig\tStartLocus\tStopLocus\t'
    
    for s in sorted(sampleList):
        firstLineOut=firstLineOut+'\t'+str(s)
        firstLineAllData=firstLineAllData+str(s)+'_totVar'+'\t'+str(s)+'_numHet'+'\t'+str(s)+'_zyg'+'\t'+str(s)+'_all'+'\t'+str(s)+'_len'+'\t'+str(s)+'_LoP'+'\t'
    out_AllData.write(firstLineAllData+"\n")
    out_DetUndet.write(firstLineOut+"\tTotDetected\tTotUndetected\n")
    
    #write the other lines
    for infokey,infovalue in locus_info_dict.iteritems():
        if infovalue.DpLocusDict:
            lineAllData=str(infokey)+'\t'+str(infovalue.Supercontig)+'\t'+str(infovalue.start)+'\t'+str(infovalue.end)+'\t'
            lineDet=str(infokey)+'\t'
            countDet=0
            countUndet=0
            od = collections.OrderedDict(sorted(infovalue.AllDataDict.items()))
            
            for AllDataValue in od.itervalues():
                
                totVarCell=0
                numHetCell=0
                zygCell='no\t'
                alleleCell=''
                lenCell=0
                deletions=''
                insertions=''
                detection=0
                gentemp=[]
                LoP='no'
                               
                if ((int(AllDataValue.StopAllele) - int(AllDataValue.StartAllele) >= int(opts.minLengthAllele)) and (int(AllDataValue.StopAllele) != 0) and (int(AllDataValue.StartAllele) != 0)):
                    lenCell=(int(AllDataValue.StopAllele) - int(AllDataValue.StartAllele))
                    #look for the existence of some deletions or insertions in this locus in this sample
                    if AllDataValue.StartDeletions:
                        for d in AllDataValue.StartDeletions:
                            deletions=deletions+str(AllDataValue.delta.encode('utf-8'))+str(d)+':'+str(AllDataValue.StopDeletions[AllDataValue.StartDeletions.index(d)])
                            lenCell=lenCell-(int(AllDataValue.StopDeletions[AllDataValue.StartDeletions.index(d)])-int(d))
                            
                    if AllDataValue.StartInsertions:
                        for i in AllDataValue.StartInsertions:
                            insertions=insertions+str(AllDataValue.iota.encode('utf-8'))+str(i)+':'+str(AllDataValue.StopInsertions[AllDataValue.StartInsertions.index(i)])
                            lenCell=lenCell+(int(AllDataValue.StopInsertions[AllDataValue.StartInsertions.index(i)])-int(i))
                            
                    alleleCell=str(AllDataValue.StartAllele)+':'+str(AllDataValue.StopAllele)+deletions+insertions+'\t'
                    
                    detection=1
                    countDet=countDet+1
                    zygCell='Homozygous\t'
                    LoP=str(0)

                    if AllDataValue.Genotypes:
                        gentemp=np.asarray(AllDataValue.Genotypes)
                        totVarCell=len(gentemp)
                        LoP=str(float(totVarCell)/float(lenCell))
                        
                        if len(gentemp[gentemp==1]):
                            numHetCell=len(gentemp[gentemp==1])
                            zygCell='Heterozygous\t'
                    
                    
                else:
                    alleleCell=str(0)+':'+str(0)+'\t'
                    countUndet=countUndet+1
                    
                lineAllData=lineAllData+str(totVarCell)+'\t'+str(numHetCell)+'\t'+zygCell+alleleCell+str(lenCell)+'\t'+LoP+'\t'
                lineDet=lineDet+str(detection)+'\t'
                
            out_AllData.write(lineAllData+"\n")
            out_DetUndet.write(lineDet+str(countDet)+'\t'+str(countUndet)+"\n")
            
            #check the NIRVs always Detected or always Undetected
            if (countDet == len(sampleList)):
                out_alwaysDet.write(str(infokey)+'\n') 
            if (countUndet == len(sampleList)):
                out_alwaysUndet.write(str(infokey)+'\n') 
               
    #close files in output        
    out_AllData.close()
    out_DetUndet.close()
    out_alwaysDet.close()
    out_alwaysUndet.close()
    
    
def AllDataDefinition(locus_info_dict, variants_dict):
      
    #for each locus
    for infokey,infovalue in locus_info_dict.iteritems():
        tempAllDataDict={}
        
        if infovalue.DpLocusDict:
            sampleList=[]
            #for each sample in that locus
            for sample,vettCoverage in infovalue.DpLocusDict.iteritems():
                
                #for each vector of coverage is defined an object ALLDataInfo
                tempAllDataOut=AllDataInfo();
                genVect=[]
                startDel=[]
                stopDel=[]
                startIn=[]
                stopIn=[]
                
                tempAllDataOut.locus_name=infovalue.locus_name;
                vect=np.asarray(map(int,vettCoverage))
                ind_extremes=np.where(vect >= int(opts.minReadsAllDef))[0]
                if ind_extremes.any():
                    tempAllDataOut.StartAllele=int(infovalue.start)+1+ind_extremes[0]
                    tempAllDataOut.StopAllele=int(infovalue.start)+1+ind_extremes[-1]
                    
                    #Control manually for loci that are, for some reason, longer than the annotation. 
                    if int(tempAllDataOut.StopAllele) > int(infovalue.end):
                        print('Control manually the following locus:')
                        print('Locus name: \t'+infokey)
                        tempAllDataOut.StopAllele=int(infovalue.end)
                            
                #info deriving from the sample file of the INDEL variants, es: pathTo/indel_allData_SampleID.txt
                strToFind=''
                strToFind=sample.split('\n')[0]
                sampleList.append(strToFind)
                tempAllDataOut.SampleName=strToFind;
                
                if sample in variants_dict.keys():
                    # selection of only deletions and insertions related to the specific locus
                    SampleVariantsDict=variants_dict[sample]
                else:
                    variants_dict[sample]={}
                    SampleVariantsDict=variants_dict[sample]
                        
                for varKey in SampleVariantsDict.iterkeys():
                    varObj=SampleVariantsDict[varKey]
                    
                    if ((varObj.CONTIG==infovalue.Supercontig) and (int(varObj.POS) in range(int(infovalue.start),int(infovalue.end)+1))):
                        num1=len(varObj.GT[varObj.GT == 1])
                        num2=len(varObj.GT[varObj.GT == 2])
                        if num1>=num2:
                            genVect.append(1)
                        else:
                            genVect.append(2)
                        
                        # only in case of INDEL  
                        if (int(len(varObj.REF))> 1 or int(len(varObj.ALT))>1):
                            if (int(len(varObj.REF))> int(len(varObj.ALT))):
                                startDel.append(int(varObj.POS))
                                stopDel.append((int(varObj.POS)+int(len(varObj.REF))))
                            else:
                                startIn.append(int(varObj.POS))
                                stopIn.append((int(varObj.POS)+int(len(varObj.ALT))))
                    
                tempAllDataOut.Genotypes=genVect
                tempAllDataOut.StartDeletions=startDel
                tempAllDataOut.StopDeletions=stopDel
                tempAllDataOut.StartInsertions=startIn
                tempAllDataOut.StopInsertions=stopIn
            
                tempAllDataDict[strToFind]=tempAllDataOut
                
        infovalue.AllDataDict=tempAllDataDict
        
    print(sampleList)
    return sampleList;


def read_variants(LocusInfoDict, VariantFileList):
    '''method to read the sample files of the SNPs and INDELs variants for each sample'''
    Variant_dictionary={}
    
    for ifile in VariantFileList:
        in_file=open(ifile.rstrip())
        
        if 'indel' in str(ifile):
            var_type_temp='INDEL'
        elif 'snp' in str(ifile):
            var_type_temp='SNP'
                
        for line in in_file:
            line=line.rstrip()
            row=line.split('\t')
            
            if line.startswith('SAMPLE_ID'):
                sample_ind=row.index("SAMPLE_ID")
                contig_ind=row.index("CONTIG")
                pos_ind=row.index("POS")
                ref_ind=row.index("REF")
                alt_ind=row.index("ALT")
                callfb_ind=row.index("CallFreebayes")
                callvd_ind=row.index("CallVardict")
                callpl_ind=row.index("CallPlatypus")
                callgk_ind=row.index("CallGatk")
                gtfb_ind=row.index("GT_Freebayes")
                gtvd_ind=row.index("GT_Vardict")
                gtpl_ind=row.index("GT_Platypus")
                gtgk_ind=row.index("GT_Gatk")
                af_med_ind=row.index("AF_median")
                sb_med_ind=row.index("STRBIAS_median")
            else:
                V=Variant()
                V.SampleID=row[sample_ind]
                V.CONTIG=row[contig_ind]
                V.POS=row[pos_ind]
                V.REF=row[ref_ind]
                V.ALT=row[alt_ind]
                V.GT=np.asarray([int(row[gtfb_ind]), int(row[gtvd_ind]), int(row[gtpl_ind]), int(row[gtgk_ind])])
                V.Call=np.asarray([int(row[callfb_ind]), int(row[callvd_ind]), int(row[callpl_ind]), int(row[callgk_ind])])
                V.ID=str(V.CONTIG)+':'+str(V.POS)+':'+str(V.REF)+':'+str(V.ALT)
                V.AF_median=row[af_med_ind]
                V.STRBIAS_median=row[sb_med_ind]
                V.varType=var_type_temp
                
                for infokey,infovalue in LocusInfoDict.iteritems():
                    if ((str(V.CONTIG)==str(infovalue.Supercontig)) and ((int(V.POS) in range(int(str(infovalue.start)),int(str(infovalue.end))+1)))):
                        V.locus_name=infokey
                if V.SampleID in Variant_dictionary:
                    Variant_dictionary[str(V.SampleID)][V.ID]=V
                else:
                    Variant_dictionary[str(V.SampleID)]={}
                    Variant_dictionary[str(V.SampleID)][V.ID]=V  
         
        in_file.close()
    
    
    return Variant_dictionary;


def  LocusSelection(matchDict,locus_info_dict,dp_dict):
    '''Selection of the locus coverage in each sample. Each locus is a dictionary containing the coverage of each sample. Each sample is a key of that dictionary'''
    #for each locus
    
    for infovalue in locus_info_dict.itervalues():
        contig_selected=[]
        index_contig_selected=[]
        final_index=[]
        cov={}
        
        for index, item in enumerate(dp_dict['Contig']):
            if (item in infovalue.Supercontig):
                contig_selected.append(dp_dict['Position'][index])
                index_contig_selected.append(index)
                 
        for pos_index, pos_item in enumerate(contig_selected):
            if (int(pos_item) in range(int(infovalue.start),int(infovalue.end)+1)):
                final_index.append(index_contig_selected[pos_index])
        
        infovalue.indices=final_index         
        
        #for each sample the coverage of all the loci is selected
        for dpkey,dpvalue in dp_dict.iteritems():
            if ((dpkey not in ['Position', 'Contig','Locus', 'Total_Depth','Average_Depth_sample']) and (final_index) and dpkey.split('Depth_for_')[1] in matchDict.keys()):
                cov[matchDict[dpkey.split('Depth_for_')[1]]]=dpvalue[int(infovalue.indices[0]):int(infovalue.indices[-1]+1)]
        	
        infovalue.DpLocusDict=cov
        
    return locus_info_dict;                 
    
    
def read_dp(filedp):
    '''method to read the depth of coverage file generated from GATKDepthOfCoverage'''
    d={}
    
    with filedp as f:
        lis = [x.split() for x in f]
    first_line=lis[0]
    
    j=0
    for x in zip(*lis):
        i=0
        row=[]
        
        for y in x:
            if (i==0):
                i=i+1
                continue
            
            else:
                row.append(y)
                i=i+1   
            
        d[first_line[j]]=row
        j=j+1
        
    contig=[]
    pos=[]       
        
    for el in d.get('Locus'):
        contig.append(el.split(':')[0])
        pos.append(el.split(':')[1])
        
    d['Contig']=contig
    d['Position']=pos
    return d;


def read_matchID(match_file):
    '''read the file of the match between the name of the sequencing data and names for the output files'''
    dictionary={}
    
    for line in match_file:
        line.rstrip()
        if ':' in list(line):
            parts = line.split(":")
            old=parts[0]
            new=parts[1].split('\n')[0]
            
            dictionary[old]=new 
            
    return dictionary;


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


def main():
    #read parameters of the script
    parser = argparse.ArgumentParser('Alldata allele extremes for each NIRV for each Sample. Output is to stdout.')
    parser.add_argument('-RegionsInfo', '--RegionsInfo', help="input file containing regions Supercontig, Start, End and locus_name tab delimited")
    parser.add_argument('-MatchID', '--MatchID', help="file including the reference to rename the sequencing data in output; oldName:newName\n")
    parser.add_argument('-GATKdp', '--GATKdp', help="GATK Depth of Coverage output file name")
    parser.add_argument('-VariantPathList', '--Variantlist', help="Paths list of the INDEL and SNP variant file")
    parser.add_argument('-minReadsAllDef', '--minReadsAllDef', help="Min number of reads required to define the extremes of the allele in AllData, default 5", default=5)
    parser.add_argument('-minLengthAllele', '--minLengthAllele', help="Min number of nucleotides required to define the allele in AllData, default 0, mining that all leghts are accepted", default=0)
    parser.add_argument('-thresholdSimilarity', '--thresholdSimilarity', help="Alleles that difference less than this threshold are considered equal in definition of the loci alleles", default=0.05)
    parser.add_argument('-outputPath', '--outputPath', help="Path of the output file")
    
    global opts 
    opts = parser.parse_args()
    
    #paths list of the variant files
    VariantList=open(opts.Variantlist)
                               
    #read locus info file and gatk dp file. Save both data in dictionary structure
    in_info=open(opts.RegionsInfo)
    LocusInfoDict=read_locus_info(in_info)
    
    match_file=open(opts.MatchID)
    MatchDict=read_matchID(match_file)
    
    in_dp=open(opts.GATKdp)
    DpFileDict=read_dp(in_dp)

    #selection of the locus indices in dp file and the sample coverages
    LocusInfoDict=LocusSelection(MatchDict,LocusInfoDict,DpFileDict)
    
    # save samples variants in dictionary structure
    Variants_dict=read_variants(LocusInfoDict, VariantList)
    
    #definition of the extremes of the alleles for each sample, for each locus
    sampleList=AllDataDefinition(LocusInfoDict,Variants_dict)
    
    #write output files
    WriteAllDataFile_DetUndetFile(sampleList, LocusInfoDict)
    
    WriteAlleles(sampleList,LocusInfoDict)
    
    VariantList.close()
    
main()
