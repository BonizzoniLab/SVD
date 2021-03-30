import argparse
import statistics
import scipy.stats as stats

''' //////////// CLASSI ////////////'''

class Caller():
	GT=''
	AO=''
	RO=''		
	AO_f=''
	AO_r=''
	DP_f=''
	DP_r=''
	DP=''
	QB=''
	Call=''
	AF=''
	StrandBias=''
	
class Freebayes(Caller):
	RO_f=''
	RO_r=''
	
class Vardict(Caller):
	RO_f=''
	RO_r=''
	ODDRATIO=''
	SSF=''
	MQ=''

class Platypus(Caller):
	RO_f=''
	RO_r=''
	pass

class Gatk(Caller):
	MQ0=''
	MQRankSum=''
	BQRankSum=''
	PhredFS=''
	StrBiasFS=''
	pass

class Features():
	GT_Freebayes='0'
	GT_Vardict='0'
	GT_Platypus='0'
	GT_Gatk='0'
	QB_Freebayes='.'
	QB_Vardict='.'
	QB_Platypus='.'
	QB_Gatk='.'
	AF_Freebayes='.'
	AF_Vardict='.'
	AF_Platypus='.'
	AF_Gatk='.'
	DP_Freebayes='.'
	DP_Vardict='.'
	DP_Platypus='.'
	DP_Gatk='.'
	CallFreebayes='0'
	CallVardict='0'
	CallPlatypus='0'
	CallGatk='0'
	STRBIAS_Freebayes='.'
	STRBIAS_Vardict='.'
	STRBIAS_Platypus='.'
	StrBiasFS_Gatk='.'
	BIAS_Vardict='.'
	ODDRATIO_Vardict='.'
	SBF_Vardict='.'
	MQ0_Gatk='.'
	MQ0F_Gatk='.'
	MQRankSum_Gatk='.'
	BQRankSum_Gatk='.'
	MQ0F_median='.'
	MQ0_median='.'
	MQ0_norm_median='.'
	MQ_Vardict='.'
	DP=float(0)
	DP_median='.'
	DP_norm_median='.'
	BQ_media='.'
	BQ_median='.'
	AF_media='.'
	AF_median='.'
	STRBIAS_media='.'
	STRBIAS_median= '.'
	
	
''' //////////// FUNZIONI ////////////'''

def get_info_freebayes(chrom,pos,ref,alt,info,format,sample,freebayes):
	'''estrae le informazioni dal vcf di freebayes'''
	
	if sample[format.index('GT')]=='1/0' or sample[format.index('GT')]=='0/1':
		freebayes.GT=1
	elif sample[format.index('GT')]=='1/1':
		freebayes.GT=2
	elif '.' in sample[format.index('GT')]:
		freebayes.GT=0
		
	if sample is not 'null' and sample[format.index('DP')] is not '0':
		try:
			freebayes.AO=float(sample[format.index('AO')])
		except:
			freebayes.AO=float(0)
		try:
			freebayes.DP=float(sample[format.index('DP')])
		except:
			freebayes.DP=0
		try:
			freebayes.RO=float(sample[format.index('RO')])
		except:
			freebayes.RO=float(0)
		try:
			freebayes.AF=freebayes.AO/freebayes.DP
		except:
			freebayes.AF='.'
	else:
		freebayes.DP=0
		freebayes.AF='.'
		freebayes.AO='.'
		freebayes.RO='.'
	
	for ind in info:
		if ind.startswith("SAF="):
			freebayes.AO_f=float(ind.split('=')[1])
		if ind.startswith("SAR="):
			freebayes.AO_r=float(ind.split('=')[1])
		if ind.startswith("SRF="):
			freebayes.RO_f=float(ind.split('=')[1])
		if ind.startswith("SRR="):
			freebayes.RO_r=float(ind.split('=')[1])

	freebayes.DP_f=float(freebayes.AO_f)+float(freebayes.RO_f)
	freebayes.DP_r=float(freebayes.AO_r)+float(freebayes.RO_r)

	try:
		freebayes.QB=float(sample[format.index('QA')])/freebayes.AO
	except:
		freebayes.QB=float(0)
	freebayes.Call=1

	
	if min((freebayes.DP_r),(freebayes.DP_f))/((freebayes.DP_r)+(freebayes.DP_f)) > 0:
		try:
			freebayes.StrandBias=1-stats.fisher_exact([[freebayes.RO_f, freebayes.RO_r], [freebayes.AO_f, freebayes.AO_r]])[1]
		except:
			freebayes.StrandBias='.'
	else:
		freebayes.StrandBias='.'


def get_info_vardict(chrom,pos,ref,alt,info,format,sample,vardict):
	'''estrae le informazioni dal vcf di vardict'''
	
	if sample[format.index('GT')]=='1/0' or sample[format.index('GT')]=='0/1':
		vardict.GT=1
	elif sample[format.index('GT')]=='1/1':
		vardict.GT=2
	elif '.' in sample[format.index('GT')]:
		vardict.GT=0
		
	if sample is not 'null' and sample[format.index('DP')] is not '0' :
		
		vardict.AO=float(sample[format.index('AD')].split(',')[1])
		vardict.RO=float(sample[format.index('AD')].split(',')[0])
		vardict.DP=float(sample[format.index('DP')])
		
		try:
			vardict.AF=float(vardict.AO/(vardict.DP))
		except:
			vardict.AF=float(0)
	else:
		vardict.AO='.'
		vardict.RO='.'
		vardict.DP=0
		vardict.AF='.'
	
	vardict.AO_f=float(sample[format.index('ALD')].split(',')[0])
	vardict.AO_r=float(sample[format.index('ALD')].split(',')[1])
	vardict.RO_f=float(sample[format.index('RD')].split(',')[0])
	vardict.RO_r=float(sample[format.index('RD')].split(',')[1])
	vardict.DP_f=vardict.AO_f+vardict.RO_f
	vardict.DP_r=vardict.AO_r+vardict.RO_r
	
	vardict.Call=1

	
	if min((vardict.DP_r),(vardict.DP_f))/((vardict.DP_r)+(vardict.DP_f)) > 0:
		try:
			vardict.StrandBias=1-stats.fisher_exact([[vardict.RO_f, vardict.RO_r], [vardict.AO_f, vardict.AO_r]])[1]
		except:
			vardict.StrandBias='.'
	else:
		vardict.StrandBias='.'
	
	for ind in info:
		if ind.startswith("QUAL="):					
			vardict.QB=float(ind.split('=')[1])
		
		if ind.startswith("MQ="):					
			vardict.MQ=float(ind.split('=')[1])		


def get_info_platypus(chrom,pos,ref,alt,info,format,sample,platypus):
	'''estrae le informazioni dal vcf di platypus'''
	
	if sample[format.index('GT')]=='1/0' or sample[format.index('GT')]=='0/1':
		platypus.GT=1
	elif sample[format.index('GT')]=='1/1':
		platypus.GT=2
	elif '.' in sample[format.index('GT')]:
		platypus.GT=0
	
	if sample is not 'null':
		
		for ind in info:
			if ind.startswith("TC="):					
				platypus.DP=float(ind.split('=')[1])
					
			if ind.startswith("TCF="):					
				platypus.DP_f=float(ind.split('=')[1])
				
			if ind.startswith("TCR="):					
				platypus.DP_r=float(ind.split('=')[1])
			
		platypus.AO=float(sample[format.index('NV')])
		platypus.RO=platypus.DP-platypus.AO
		
	if platypus.DP is not '0':
		try:
			platypus.AF=float(platypus.AO/(platypus.DP))
		except:
			platypus.AF=float(0)
	else:
		platypus.AO='.'
		platypus.RO='.'
		platypus.DP=0
		platypus.AF='.'
	
	for ind in info:
		if ind.startswith("NF="):					
			platypus.AO_f=float(ind.split('=')[1])
				
		if ind.startswith("NR="):					
			platypus.AO_r=float(ind.split('=')[1])
					
	platypus.RO_f=platypus.DP_f-platypus.AO_f
	platypus.RO_r=platypus.DP_r-platypus.AO_r
	
	platypus.Call=1
	platypus.QB='.'
	
	if int(platypus.DP) is not 0:
		if min((platypus.DP_r),(platypus.DP_f))/((platypus.DP_r)+(platypus.DP_f)) > 0:
			try:
				platypus.StrandBias=1-stats.fisher_exact([[platypus.RO_f, platypus.RO_r], [platypus.AO_f, platypus.AO_r]])[1]
			except:
				platypus.StrandBias='.'
		else:
			platypus.StrandBias='.'

			
def get_info_gatk(chrom,pos,ref,alt,info,format,sample,gatk):
	'''estrae le informazioni dal vcf di gatk'''
	
	if sample[format.index('GT')]=='1/0' or sample[format.index('GT')]=='0/1':
		gatk.GT=1
	elif sample[format.index('GT')]=='1/1':
		gatk.GT=2
	elif '.' in sample[format.index('GT')]:
		gatk.GT=0
		
	if sample is not 'null' and sample[format.index('DP')] is not '0' :
		
		gatk.DP=float(sample[format.index('DP')])
		gatk.AO=float(sample[format.index('AD')].split(',')[1])
		gatk.RO=float(sample[format.index('AD')].split(',')[0])
		
		try:
			gatk.AF=float(gatk.AO/(gatk.DP))
		except:
			gatk.AF=float(0)
	else:
		gatk.AO='.'
		gatk.RO='.'
		gatk.DP=0
		gatk.AF='.'
	
	gatk.Call=1
	gatk.QB='.'
	gatk.BQRankSum=0
	gatk.MQRankSum=0
	
	for ind in info:
		if ind.startswith("BaseQRankSum="):					
			gatk.BQRankSum=float(ind.split('=')[1])
				
		if ind.startswith("MQRankSum="):					
			gatk.MQRankSum=float(ind.split('=')[1])
		
		if ind.startswith("MQ0="):					
			gatk.MQ0=float(ind.split('=')[1])
		
		if ind.startswith("FS="):					
			gatk.PhredFS=float(ind.split('=')[1])
			
	try:		
		gatk.StrBiasFS=1-pow(10,-gatk.PhredFS/10)
	except:
		gatk.StrBiasFS='.'
		
	try:		
		gatk.MQ0F=gatk.MQ0/gatk.DP
	except:
		gatk.MQ0F='0'		


def set_features_snp(dictionary):
	'''setta i valori delle features in base alle info estratte dai vcf'''
	for variante in dictionary.keys():
		features=Features()
		varc_array=dictionary.get(variante)
		
		vett_MBQ=[]
		vett_DP=[]
		vett_MQ0=[]
		index=0

		for varcall in varc_array:
			if varcall is not "":
				vett_MBQ=vett_MBQ+[varcall.QB]
				vett_DP=vett_DP+[varcall.DP]
				
				if index == 0:

					features.GT_Freebayes=varc_array[0].GT
					features.QB_Freebayes=varc_array[0].QB
					features.AF_Freebayes=varc_array[0].AF	
					features.CallFreebayes=varc_array[0].Call
					features.STRBIAS_Freebayes=varc_array[0].StrandBias

				elif index == 1:
						
					features.GT_Vardict=varc_array[1].GT
					features.QB_Vardict=varc_array[1].QB
					features.MQ_Vardict=varc_array[1].MQ
					features.AF_Vardict=varc_array[1].AF
					features.CallVardict=varc_array[1].Call
					features.STRBIAS_Vardict=varc_array[1].StrandBias
						
				elif index == 2:
					
					features.GT_Platypus=varc_array[2].GT
					features.AF_Platypus=varc_array[2].AF
					features.CallPlatypus=varc_array[2].Call
					features.STRBIAS_Platypus=varc_array[2].StrandBias
				
				elif index == 3:
					
					features.GT_Gatk=varc_array[3].GT
					features.AF_Gatk=varc_array[3].AF
					features.CallGatk=varc_array[3].Call
					features.MQ0_Gatk=varc_array[3].MQ0
					features.MQ0F_Gatk=varc_array[3].MQ0F
					features.MQRankSum_Gatk=varc_array[3].MQRankSum
					features.BQRankSum_Gatk=varc_array[3].BQRankSum
					features.StrBiasFS_Gatk=varc_array[3].StrBiasFS
			
			index = index + 1	
		
		vett_MQ0=[features.MQ0_Gatk]
		vett_MQ0F=[features.MQ0F_Gatk]
		vett_AF_media=[features.AF_Freebayes,features.AF_Vardict,features.AF_Platypus,features.AF_Gatk]
		vett_STRB_media=[features.STRBIAS_Freebayes,features.STRBIAS_Vardict,features.STRBIAS_Platypus,features.StrBiasFS_Gatk]
		AF_med=0
		MQ0_med=0
		MQ0F_med=0
		SB_media=0
		nDP=0
		nMBQ=0
		
		i=0
		v=[]
		for dp in vett_DP:
			if dp and dp is not '':
				nDP= float(nDP)+float(dp)
				v=v+[float(dp)]
				i=i+1
		try:
			features.DP=nDP/i
			features.DP_median= statistics.median(v)
		except:
			features.DP='.'
		
		try:
			features.DP_norm_median= features.DP_median/float(opts.expectedMeanDP)
		except:
			features.DP_norm_median='.'
			
		i=0
		v=[]
		for bq in vett_MBQ:
			if bq is not '.':
				nMBQ= float(nMBQ)+float(bq)
				v=v+[float(bq)]
				i=i+1
		try:
			features.BQ_media=nMBQ/i
			features.BQ_median= statistics.median(v)
		except:
			features.BQ_media='.'

		i=0
		v=[]
		for strb in vett_STRB_media:
			if strb is not '.':
				SB_media=float(SB_media) + float(strb)
				v=v+[float(strb)]
				i=i+1
		try:
			features.STRBIAS_media= SB_media/i
			features.STRBIAS_median= statistics.median(v)
		except:
			features.STRBIAS_media='.'
			
		i=0
		v=[]
		for af in vett_AF_media:
			if af is not '.' and af is not '0':
				AF_med=float(AF_med) + float(af)
				v=v+[float(af)]
				i=i+1
		try:
			features.AF_media= AF_med/i
			features.AF_median= statistics.median(v)
		except:
			features.AF_media='.'
			
		i=0
		v=[]
		for mq in vett_MQ0:
			if mq is not '.':
				MQ0_med=float(MQ0_med) + float(mq)
				v=v+[float(mq)]
				i=i+1
		try:
			features.MQ0_media= MQ0_med/i
			features.MQ0_median= statistics.median(v)
		except:
			features.MQ0_media='.'
		
		try:
			features.MQ0_norm_median= features.MQ0_median/float(opts.expectedMeanDP)
		except:
			features.MQ0_norm_median='.'
			
		i=0
		v=[]
		for mqf in vett_MQ0F:
			if mqf is not '.':
				MQ0F_med=float(MQ0F_med) + float(mqf)
				v=v+[float(mqf)]
				i=i+1
		try:
			features.MQ0F_media= MQ0_med/i
			features.MQ0F_median= statistics.median(v)
		except:
			features.MQ0F_media='.'

		dictionary[variante]= varc_array + [features]

def switch_indel(dictionary,ID,index,chrom,pos,ref,alt,info,format,sample):
	'''tramite index richiama la funzione di estrazione delle informazioni del variant caller associato all'indice'''
	if dictionary.has_key(ID):
		vettore=dictionary[ID]
	else:
		vettore=['','','','']

	if index==0:
		freebayes=Freebayes()
		get_info_freebayes(chrom,pos,ref,alt,info,format,sample,freebayes)
		if freebayes.AF != 0.0: 
			vettore[0]=freebayes
	elif index==1:
		vardict=Vardict()
		get_info_vardict(chrom,pos,ref,alt,info,format,sample,vardict)
		if vardict.AF != 0.0 and float(vardict.MQ) >= float(opts.mq): 
			vettore[1]=vardict
	elif index==2:
		platypus=Platypus()
		get_info_platypus(chrom,pos,ref,alt,info,format,sample,platypus)
		if platypus.AF != 0.0: 
			vettore[2]=platypus
	elif index==3:
		gatk=Gatk()
		get_info_gatk(chrom,pos,ref,alt,info,format,sample,gatk)
		if gatk.AF != 0.0: 
			vettore[3]=gatk
				
	dictionary[ID]=vettore

def read(iterable,index,dictionary):
	'''legge il vcf e splitta le varie sezioni'''
	chrom=''
	alt=''
	pos=''
	ref=''
	format=''
	info=''
	sample=''
	for line in iterable:
		line.rstrip()
		if line.startswith('#'):
			continue
		else:
			ind=0
			parts = line.split("\t")
			chrom=parts[0]
			pos=parts[1]
			ref=parts[3]
			alt=parts[4]
			ID='\t'.join([chrom,pos,ref,alt])
			INFO=parts[7]
			info=INFO.split(";") 
			FORMAT=parts[8]
			format=FORMAT.split(":")
			SAMPLE=parts[9]
			sample=SAMPLE.split(":")
			
			if len(ref)==1 and len(alt)==1 or len(sample)==1 or SAMPLE.startswith('0/0:') :
				continue
			else:
				if "VD" in format:
					if sample[format.index('VD')] is '0':
						continue
					else:
						switch_indel(dictionary,ID,index,chrom,pos,ref,alt,info,format,sample)
				else:	
					switch_indel(dictionary,ID,index,chrom,pos,ref,alt,info,format,sample)

def control(dictionary):
	''' esegue un controllo sulle varianti, se non hanno variant caller che le chiama vengono eliminate'''
	for variante in dictionary.keys():
		if dictionary[variante][:4] == ['','','','']:
			del dictionary[variante]
				
def print_var_indel(dictionary):
	print "SAMPLE_ID\tCONTIG\tPOS\tREF\tALT\tCallFreebayes\tCallVardict\tCallPlatypus\tCallGatk\tGT_Freebayes\tGT_Vardict\tGT_Platypus\tGT_Gatk\tBQ_Freebayes\tBQ_Vardict\tBQ_mean\tBQ_median\tDP_mean\tDP_median\tDP_norm_median\tAF_Freebayes\tAF_Vardict\tAF_Platypus\tAF_Gatk\tAF_mean\tAF_median\tSTRBIAS_Freebayes\tSTRBIAS_Vardict\tSTRBIAS_Platypus\tStrBiasFS_Gatk\tSTRBIAS_mean\tSTRBIAS_median\tMQRankSum_Gatk\tBQRankSum_Gatk\tMQ0F_Gatk\tMQ0_Gatk\tMQ0_GATK_norm\tMQ_Vardict"
	for variante in dictionary.keys():
		features = dictionary.get(variante)[4]

		print '\t'.join([opts.sample,variante,str(features.CallFreebayes),str(features.CallVardict),str(features.CallPlatypus),str(features.CallGatk),
			str(features.GT_Freebayes),str(features.GT_Vardict),str(features.GT_Platypus),str(features.GT_Gatk),
			str(features.QB_Freebayes),str(features.QB_Vardict),str(features.BQ_media),str(features.BQ_median),str(features.DP),str(features.DP_median),str(features.DP_norm_median),
			str(features.AF_Freebayes),str(features.AF_Vardict),str(features.AF_Platypus),str(features.AF_Gatk),str(features.AF_media),str(features.AF_median),
			str(features.STRBIAS_Freebayes),str(features.STRBIAS_Vardict),str(features.STRBIAS_Platypus),str(features.StrBiasFS_Gatk),str(features.STRBIAS_media),str(features.STRBIAS_median),
			str(features.MQRankSum_Gatk),str(features.BQRankSum_Gatk),str(features.MQ0F_Gatk),str(features.MQ0_Gatk),str(features.MQ0_norm_median),str(features.MQ_Vardict)]).rstrip()
		
def main():

	
	parser = argparse.ArgumentParser('Parse VCF output from Variant callers to output a variant_dataset.txt.  Output is to stdout.')
	parser.add_argument('-fb', '--freebayes', help="Freebayes vcf output file name")
	parser.add_argument('-vd', '--vardict', help="Vardict vcf output file name")
	parser.add_argument('-pl', '--platypus', help="Platypus vcf output file name")
	parser.add_argument('-gk', '--gatk', help="Gatk vcf output file name")
	parser.add_argument('-s','--sample',help="Sample name")
	parser.add_argument('-mDP','--expectedMeanDP',help="Expected Mean Coverage")
	parser.add_argument('-mqVDThreshold','--mq',help="Threshold used in the variant calling phase to filter variants")
	
	global opts 
	opts = parser.parse_args()
	
	callers = [opts.freebayes,opts.vardict,opts.platypus,opts.gatk]
	varianti = dict() 
	index=0;
	for vcf in callers:
		in_file = open(vcf)
		read(in_file,index,varianti)
		index = index + 1
		
	set_features_snp(varianti)
	control(varianti)
	print_var_indel(varianti)
	
main()
