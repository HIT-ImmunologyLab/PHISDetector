import os,re,argparse
import numpy.random.common
import numpy.random.bounded_integers
import numpy.random.entropy
import numpy as np
import threading
import json
from Bio import Entrez, SeqIO
import multiprocessing

def get_root_path():
	pathfolders = os.environ['PATH'].split(os.pathsep)
	pathfolders.reverse()
	pathfolders.append(os.getcwd())
	pathfolders.reverse()
	CURRENTDIR = os.getcwd()
	root_path = ""
	for folder in pathfolders:
		try:
			if "PHIS" in os.listdir(folder) and "PHIS_overall_phage_to_host_two_criteria" in os.listdir(folder) and "PHIS_crispr" in os.listdir(folder) and "PHIS_prophage" in os.listdir(folder) and "PHIS_blast" in os.listdir(folder) and "PHIS_protein_protein_int" in os.listdir(folder):
				root_path = os.path.dirname(folder)
				break
		except:
			pass
	try:
		if root_path == "" and os.sep in sys.argv[0] and "PHIS" in os.listdir(sys.argv[0].rpartition(os.sep)[0]) and "PHIS_overall_phage_to_host_two_criteria" in os.listdir(sys.argv[0].rpartition(os.sep)[0]):
			root_path = os.path.dirname(sys.argv[0].rpartition(os.sep)[0])
			#os.chdir(root_path)
	except:
		pass
	if root_path == "":
		print("Error: Please add the PHIS installation directory to your $PATH environment variable before running the executable from another folder.")
		sys.exit(1)
	return root_path

def pred_orf(input_phage,outfile):
	script = os.path.join(root_path,'software','FragGeneScan','run_FragGeneScan.pl')
	command = script+" -genome="+input_phage+" -out="+outfile+" -complete=1 -train=complete -thread=20"
	os.system(command)

def getFaaFromGB(input_file,outdir):   #parse protein from genbank files in phaster
	special_pros = ['capsid','head','plate','tail','coat','portal','holin','integrase','transposase','terminase','protease','lysis','bacteriocin','tRNA']
	records = SeqIO.parse(input_file, "gb")
	file_name = os.path.basename(input_file).split()[0].split('.')[0].replace(' ','-')
	counter = 0
	outFileName = os.path.join(outdir,file_name+'.faa')
	outfiledef = os.path.join(outdir,file_name+'_protein_def')
	savefile = open(outFileName, 'w')
	savefile_protein = open(outfiledef,'w')
	for record in records:
		for feature in record.features:
			if feature.type == 'CDS':                
				location = feature.location
				if str(location).find('+') != -1:
					direction = '+'
				elif str(location).find('-') != -1:
					direction = '-'
				if '<' in str(location):
					location = str(location).replace('<','')
				if '>' in str(location):
					location = str(location).replace('>','')
				locations = re.findall("\d+\.?\d*",str(location))
				min_start = locations[0]
				max_end = locations[1]
				for loc in locations:
					if int(loc)<int(min_start):
						min_start = loc
					if int(loc)>int(max_end):
						max_end = loc
				location = min_start+'_'+max_end+'_'+direction
				counter = counter+1
				if 'product' in feature.qualifiers:
					product = feature.qualifiers['product'][0]	  
					if 'protein_id' in feature.qualifiers:
						proteinId = feature.qualifiers['protein_id'][0]
					else:
						if 'inference' in feature.qualifiers:
							strInference = str(feature.qualifiers['inference'])
							if 'RefSeq' in strInference:
								proteinId = strInference.split('RefSeq:')[1].rstrip(']').rstrip('\'')
							elif 'SwissProt' in strInference:
								proteinId = strInference.split('SwissProt:')[1].rstrip(']').rstrip('\'')
							else:
								proteinId = 'unknown'
						else:
							proteinId = 'unknown'
					if 'translation' in feature.qualifiers:
						translation = feature.qualifiers['translation'][0]
						savefile_protein.write(proteinId+'\t'+product+'\n')
						strain_id = list(strain_inf_dict.keys())[0]
						savefile.write('>' +strain_id+'|' + str(location)+ '|' + str(proteinId)+'\n')					
					# savefile.write(">"+fileID+ '_'+str(counter)+'\n')                   
						if translation[-1] == '\n':
							savefile.write(translation)
						else:
							savefile.write(translation + '\n')
	savefile.close()
	savefile_protein.close()

def diamond_blastp_nomax(file,outfile,database):
	num_threads = 20
	format = 6
	evalue = 1
	identity = 0.4
	coverage = 0.7
	diamond_path = os.path.join(root_path,'software','diamond','diamond')
	script = diamond_path+" blastp -d "+database+" -q "+file+" -f "+str(format)+" -e "+str(evalue)+" -o "+outfile+" -p "+str(num_threads)+" --id "+str(identity)+" -k 1000000"
	os.system(script)

def blastn_prophage(file,outfile):
	prophage_database = os.path.join(root_path,'database/db/database/bacteria_prophage_nucl','bacteria_prophage_nucl_db')
	format = 6
	evalue = 0.01
	blastn_path = os.path.join(root_path,'software','blast+','blastn')
	script = blastn_path+" -db "+prophage_database+" -query "+file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile+" -num_threads 30 -word_size 11 -perc_identity 70 -reward 1 -penalty -2 -gapopen 0 -gapextend 0 -max_target_seqs 100000"
	os.system(script)

def blastp_prophage(file,outfile):
	num_threads = 20
	format = 6
	database = os.path.join(root_path,'database/db/database/bacteria_prophage_protein','bacteria_prophage_protein_db')
	diamond_blastp_nomax(file,outfile,database)

def combine(lst):
	lst.sort(key=lambda obj:int(obj[0]))
	templist = lst[0]
	combinedList = []
	if len(lst)==1:
		return lst
	for newlist in lst[1:]:
		newListmin = min(int(newlist[0]),int(newlist[1]))
		newListmax = max(int(newlist[0]), int(newlist[1]))
		tempListmax = max(int(templist[0]), int(templist[1]))
		if int(newListmin) <= int(tempListmax):
			templist[1] = str(max(int(tempListmax),int(newListmax)))
		else:
			combinedList.append(templist)
			templist = newlist
	combinedList.append(templist)
	return combinedList

def parse_phage_blastn_prophage_6(outfile,outdir):
	with open(outfile) as f:
		contents = f.readlines()
	result_dict = {}
	filter_file = os.path.join(outdir,'phage_blastn_prophage_result.txt')
	outfile_dict = os.path.join(outdir,'phage_blastn_prophage_result_dict')
	f_filter = open(filter_file,'w')
	f_filter.write('phage_id\tphage_def\thit_bac_id\thit_bac_def\tprophage_id\tprophage_length\thit_length\tidentity\tcoverage\n')
	for line in contents:
		line = line.strip().split('\t')
		phage_id = line[0].split('.')[0]
		hit_prophage = line[1]
		bac_id = hit_prophage.split('|')[0].split('.')[0]
		method = hit_prophage.split('|')[-1]
		prophage_length = abs(int(hit_prophage.split('|')[1].split(':')[1])-int(hit_prophage.split('|')[1].split(':')[0])+1)
		identity = line[2]
		alignment_length = int(line[3])
		mismatch = int(line[4])
		bit_score = float(line[-1])
		start = hit_prophage.split('|')[1].split(':')[0]
		end = hit_prophage.split('|')[1].split(':')[1]
		hit_start = line[-4]
		hit_end = line[-3]
		if (float(identity) >= float(min_id_prophage)):#mismatch<=2
			if phage_id not in result_dict.keys():
				result_dict.update({phage_id:{}})
			#coverage = alignment_length/float(prophage_length)
			if bac_id not in result_dict[phage_id].keys():
				result_dict[phage_id].update({bac_id:{}})
			if method not in result_dict[phage_id][bac_id].keys():
				result_dict[phage_id][bac_id].update({method:{}})
			if start+':'+end not in result_dict[phage_id][bac_id][method].keys():
				result_dict[phage_id][bac_id][method].update({start+':'+end:[]})
			result_dict[phage_id][bac_id][method][start+':'+end].append([hit_prophage,identity,[hit_start,hit_end],bit_score])
	result_dict1 = {}
	for phage_id in result_dict.keys():
		try:
			phage_def = strain_inf_dict[phage_id.split('.')[0]]
		except:
			phage_def = phage_id	
		for bac_id in result_dict[phage_id].keys():
			bac_def = bac_inf_dict[bac_id]
			for method in result_dict[phage_id][bac_id].keys():
				identitys = []
				for region in result_dict[phage_id][bac_id][method].keys():
					hit_infos = result_dict[phage_id][bac_id][method][region]
					hit_prophage = hit_infos[0][0]
					prophage_length = abs(int(region.split(':')[1])-int(region.split(':')[0]))+1
					hit_region = [x[2] for x in hit_infos]
					identitys = [float(x[1]) for x in hit_infos]
					identity = np.mean(identitys)
					merge_region = combine(hit_region)
					hit_length = 0
					for region1 in merge_region:
						region_length = abs(int(region1[0])-int(region1[1]))+1
						hit_length = hit_length+region_length
					coverage = float(hit_length)/prophage_length
					if coverage>=float(min_cov_prophage)/100:
						if phage_id not in result_dict1.keys():
							result_dict1.update({phage_id:{}})
						if bac_id not in result_dict1[phage_id].keys():
							result_dict1[phage_id].update({bac_id:{}})
						if method not in result_dict1[phage_id][bac_id].keys():
							result_dict1[phage_id][bac_id].update({method:{}})
						result_dict1[phage_id][bac_id][method].update({region:[]})						
						result_dict1[phage_id][bac_id][method][region] = [hit_prophage,identity,coverage,hit_infos]					
					f_filter.write('\t'.join([phage_id,phage_def,bac_id,bac_def,hit_prophage,str(prophage_length),str(hit_length),str(identity),str(coverage)])+'\n')
	np.save(outfile_dict,result_dict1)
	f_filter.close()

def parse_phage_blastp_prophage_6(outfile,outdir):
	with open(outfile) as f:
		contents = f.readlines()
	result_dict = {}
	outfile_dict = os.path.join(outdir,'phage_blastp_prophage_result_dict')
	resufile = os.path.join(outdir,'phage_blastp_prophage_result.txt')
	f_result = open(resufile,'w')
	f_result.write('phage_id\tphage_def\thit_bac_id\thit_bac_def\tprophage_id\tprophage_pro_num\thit_pro_percent\thit_pro_num\n')
	for line in contents:
		line = line.strip().split('\t')
		if type == 'fasta':
			phage_id = '_'.join(line[0].split('_')[0:-3]).split('.')[0]
		else:
			phage_id = line[0].split('|')[0].split('.')[0]
		if phage_id not in result_dict.keys():
			result_dict.update({phage_id:{}})	
		hit_prophage = line[1]
		bac_id = hit_prophage.split('|')[0].split('.')[0]
		method = hit_prophage.split('|')[-2]
		prophage_region_pro_num = float(hit_prophage.split('|')[-1])
		bac_id = hit_prophage.split('|')[0]
		start = hit_prophage.split('|')[1].split(':')[0]
		end = hit_prophage.split('|')[1].split(':')[1]
		if bac_id not in result_dict[phage_id].keys():
			result_dict[phage_id].update({bac_id:{}})
		if method not in result_dict[phage_id][bac_id].keys():
			result_dict[phage_id][bac_id].update({method:{}})
		if start+':'+end not in result_dict[phage_id][bac_id][method].keys():
			result_dict[phage_id][bac_id][method].update({start+':'+end:[]})
		result_dict[phage_id][bac_id][method][start+':'+end].append([line[0],hit_prophage,prophage_region_pro_num])
	result_dict1 = {}
	for phage_id in result_dict.keys():
		try:
			phage_def = strain_inf_dict[phage_id.split('.')[0]]
		except:
			phage_def = phage_id
		for bac_id in result_dict[phage_id].keys():
			bac_def = bac_inf_dict[bac_id]
			for method,regions_infos in result_dict[phage_id][bac_id].items():
				for region,hit_prophages in regions_infos.items():
					region_pro_num = hit_prophages[0][-1]
					prophage_hit_pro_num = len(list(set([item[1] for item in hit_prophages])))
					percent = prophage_hit_pro_num/float(region_pro_num)
					#print(len(hit_prophages),region_pro_num)
					if percent>=float(min_per_prophage)/100:
						if phage_id not in result_dict1.keys():
							result_dict1.update({phage_id:{}})
						if bac_id not in result_dict1[phage_id].keys():
							result_dict1[phage_id].update({bac_id:{}})
						if method not in result_dict1[phage_id][bac_id].keys():
							result_dict1[phage_id][bac_id].update({method:{}})
						result_dict1[phage_id][bac_id][method].update({region:[]})
						result_dict1[phage_id][bac_id][method][region] = [percent,len(hit_prophages),region_pro_num,hit_prophages]
					hit_prophage_id = '|'.join(hit_prophages[0][1].split('|')[0:2]+[hit_prophages[0][1].split('|')[3]])
					f_result.write('\t'.join([phage_id,phage_def,bac_id,bac_def,hit_prophage_id,str(region_pro_num),str(percent),str(len(hit_prophages))])+'\n')
					f_result.flush()
	np.save(outfile_dict,result_dict1)
	f_result.close()

def parse_prophage_blastp_phagedb_6(outfile,outdir):
	result_dict = {}
	outfile_dict = os.path.join(outdir,'prophage_blastp_phagedb_result_dict')
	resufile = os.path.join(outdir,'prophage_blastp_phagedb_result.txt')
	f_result = open(resufile,'w')
	f_result.write('bac_id\tbac_def\tprophage_id\thit_phage_id\thit_phage_def\tprophage_pro_num\thit_pro_percent\thit_pro_num\n')
	if os.path.exists(outfile):
		with open(outfile) as f:
			contents = f.readlines()
		for line in contents:
			line = line.strip().split('\t')
			prophage_id = line[0]
			phage_pro = line[1]
			if type=='fasta':
				phage_id = '_'.join(phage_pro.split('_')[0:-3]).split('.')[0]
			else:
				phage_id = phage_pro.split('|')[0].split('.')[0]
			bac_id = prophage_id.split('|')[0].split('.')[0]
			method = prophage_id.split('|')[-2]
			prophage_region_pro_num = float(prophage_id.split('|')[-1])
			bac_id = prophage_id.split('|')[0]
			start = prophage_id.split('|')[1].split(':')[0]
			end = prophage_id.split('|')[1].split(':')[1]
			if bac_id not in result_dict.keys():
				result_dict.update({bac_id:{}})
			if phage_id not in result_dict[bac_id].keys():
				result_dict[bac_id].update({phage_id:{}})
			if method not in result_dict[bac_id][phage_id].keys():
				result_dict[bac_id][phage_id].update({method:{}})
			if start+':'+end not in result_dict[bac_id][phage_id][method].keys():
				result_dict[bac_id][phage_id][method].update({start+':'+end:[]})
			result_dict[bac_id][phage_id][method][start+':'+end].append([prophage_id,phage_pro,prophage_region_pro_num])
	result_dict1 = {}
	for bac_id in result_dict.keys():
		try:
			bac_def = strain_inf_dict[bac_id.split('.')[0]]
		except:
			bac_def = bac_id
		for phage_id in result_dict[bac_id].keys():
			phage_def = phage_inf_dict[phage_id]
			for method,regions_infos in result_dict[bac_id][phage_id].items():
				for region,hit_infos in regions_infos.items():
					region_pro_num = hit_infos[0][-1]
					prophage_hit_pro_num = len(list(set([item[0] for item in hit_infos])))
					percent = prophage_hit_pro_num/float(region_pro_num)
					if percent>=float(min_per_prophage)/100:
						if bac_id not in result_dict1.keys():
							result_dict1.update({bac_id:{}})
						if phage_id not in result_dict1[bac_id].keys():
							result_dict1[bac_id].update({phage_id:{}})
						if method not in result_dict1[bac_id][phage_id].keys():
							result_dict1[bac_id][phage_id].update({method:{}})
						result_dict1[bac_id][phage_id][method].update({region:[]})
						result_dict1[bac_id][phage_id][method][region] = [percent,len(hit_infos),region_pro_num,hit_infos]
					prophage_id = '|'.join(hit_infos[0][0].split('|')[0:2]+[hit_infos[0][0].split('|')[3]])
					f_result.write('\t'.join([bac_id,bac_def,prophage_id,phage_id,phage_def,str(region_pro_num),str(percent),str(len(hit_infos))])+'\n')
					f_result.flush()
	np.save(outfile_dict,result_dict1)
	f_result.close()

def parse_prophage_blastn_phagedb_6(outfile,outdir):
	with open(outfile) as f:
		contents = f.readlines()
	result_dict = {}
	outfile_dict = os.path.join(outdir,'prophage_blastn_phagedb_result_dict')
	resufile = os.path.join(outdir,'prophage_blastn_phagedb_result.txt')
	f_result = open(resufile,'w')
	f_result.write('bac_id\tbac_def\tprophage_id\thit_phage_id\thit_phage_def\tprophage_length\thit_length\tidentity\tcoverage\n')
	for line in contents:
		line = line.strip().split('\t')
		prophage_id = line[0]
		phage_id = line[1].split('.')[0]
		bac_id = prophage_id.split('|')[0].split('.')[0]
		method = prophage_id.split('|')[-1]
		prophage_length = abs(int(prophage_id.split('|')[1].split(':')[1])-int(prophage_id.split('|')[1].split(':')[0])+1)
		identity = line[2]
		alignment_length = int(line[3])
		mismatch = int(line[4])
		bit_score = float(line[-1])
		start = prophage_id.split('|')[1].split(':')[0]
		end = prophage_id.split('|')[1].split(':')[1]
		hit_start = line[-4]
		hit_end = line[-3]
		if float(identity) >= float(min_id_prophage):#mismatch<=2
			if bac_id not in result_dict.keys():
				result_dict.update({bac_id:{}})
			#coverage = alignment_length/float(prophage_length)
			if phage_id not in result_dict[bac_id].keys():
				result_dict[bac_id].update({phage_id:{}})
			if method not in result_dict[bac_id][phage_id].keys():
				result_dict[bac_id][phage_id].update({method:{}})
			if start+':'+end not in result_dict[bac_id][phage_id][method].keys():
				result_dict[bac_id][phage_id][method].update({start+':'+end:[]})
			result_dict[bac_id][phage_id][method][start+':'+end].append([prophage_id,identity,[hit_start,hit_end],bit_score])
	result_dict1 = {}
	for bac_id in result_dict.keys():
		try:
			bac_def = strain_inf_dict[bac_id.split('.')[0]]
		except:
			bac_def = bac_id
		for phage_id in result_dict[bac_id].keys():
			phage_def = phage_inf_dict[phage_id]
			for method in result_dict[bac_id][phage_id].keys():
				identitys = []
				for region in result_dict[bac_id][phage_id][method].keys():
					hit_infos = result_dict[bac_id][phage_id][method][region]
					hit_prophage = hit_infos[0][0]
					prophage_length = abs(int(region.split(':')[1])-int(region.split(':')[0]))+1
					hit_region = [x[2] for x in hit_infos]
					identitys = [float(x[1]) for x in hit_infos]
					identity = np.mean(identitys)
					merge_region = combine(hit_region)
					hit_length = 0
					for region1 in merge_region:
						region_length = abs(int(region1[0])-int(region1[1]))+1
						hit_length = hit_length+region_length
					coverage = float(hit_length)/prophage_length
					if coverage>=float(min_cov_prophage)/100:
						if bac_id not in result_dict1.keys():
							result_dict1.update({bac_id:{}})
						if phage_id not in result_dict1[bac_id].keys():
							result_dict1[bac_id].update({phage_id:{}})
						if method not in result_dict1[bac_id][phage_id].keys():
							result_dict1[bac_id][phage_id].update({method:{}})
						result_dict1[bac_id][phage_id][method].update({region:[]})
						result_dict1[bac_id][phage_id][method][region] = [prophage_id,identity,coverage,hit_infos]
					f_result.write('\t'.join([bac_id,bac_def,prophage_id,phage_id,phage_def,str(prophage_length),str(hit_length),str(identity),str(coverage)])+'\n')
					f_result.flush()
	np.save(outfile_dict,result_dict1)
	f_result.close()

def filter_identity_coverage(blastp_file,filter_file,target):
	f_result = open(filter_file,'w')
	with open(blastp_file) as f:
		contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		query_pro = line[0]
		hit_pro = line[1]
		identity = line[2]
		if target == 'query':
			query_start = line[-6]
			query_end = line[-5]
			hit_length = abs(int(query_start)-int(query_end))+1
			pro_length = int((abs(int(query_pro.split('|')[2].split('_')[-3])-int(query_pro.split('|')[2].split('_')[-2]))+1)/3)
		else:
			hit_start = line[-4]
			hit_end = line[-3]
			hit_length = abs(int(hit_start)-int(hit_end))+1
			pro_length = int((abs(int(hit_pro.split('|')[2].split('_')[-3])-int(hit_pro.split('|')[2].split('_')[-2]))+1)/3)
		coverage = float(hit_length)/pro_length
		if (coverage>=0.7) and ((float(identity)/100)>=0.4):
			f_result.write('\t'.join(line)+'\t'+str(coverage)+'\n')
			f_result.flush()
	f_result.close()

def phage_to_hosts(phage_file,phage_pro,outdir):
	outfile_blastn = os.path.join(outdir,'phage_blastn_prophage')
	outfile_blastp = os.path.join(outdir,'phage_blastp_prophage')
	m1 = threading.Thread(target=blastn_prophage,args=(phage_file,outfile_blastn,))
	m2 = threading.Thread(target=blastp_prophage,args=(phage_pro,outfile_blastp,))
	m1.start()
	m2.start()
	m1.join()
	m2.join()
	filter_file = os.path.join(outdir,'phage_blastp_prophage_identity_0.4_coverage_0.7')
	filter_identity_coverage(outfile_blastp,filter_file,'subject')
	outfile_blastn_dict_file = os.path.join(outdir,'phage_blastn_prophage_result_dict.npy')
	outfile_blastp_dict_file = os.path.join(outdir,'phage_blastp_prophage_result_dict.npy')
	m1 = threading.Thread(target=parse_phage_blastn_prophage_6,args=(outfile_blastn,outdir,))
	m2 = threading.Thread(target=parse_phage_blastp_prophage_6,args=(filter_file,outdir,))
	m1.start()
	m2.start()
	m1.join()
	m2.join()
	try:
		get_prophage_result(outfile_blastp_dict_file,outfile_blastn_dict_file,outdir,'p')
	except:
		print('parse prophage feature value error!')

def get_inf(file):
	with open(file) as f:
		contents = f.read().strip()
	if contents[0]=='>':
		type ='fasta'
		try:
			record = SeqIO.read(file,"fasta")
			strain_id = record.id
			strain_def = record.description
		except:
			strain_id = 'unknown'
			strain_def = 'unknown'
	else:
		type="genbank"
		record = SeqIO.read(file,"genbank")
		strain_id = record.id
		strain_def = record.description
	return strain_id,strain_def,type

def mkdir(dirname):
	command = "mkdir -p "+dirname
	os.system(command)

def getFnaFromGB(fileName,outDir,fileID):
	outFileName = outDir + '/' + fileID + '.fna'
	handle = open(fileName)
	if os.path.exists(outFileName):
		os.remove(outFileName)
	SeqIO.convert(handle, 'genbank', outFileName, 'fasta')

def get_inf(file,outdir):
	with open(file) as f:
		contents = f.read().strip()
	if contents[0]=='>':
		type ='fasta'
		try:
			record = SeqIO.read(file,"fasta")
			strain_id = record.id
			strain_def = record.description
			strain_info_dict = {strain_id.split()[0].split('.')[0]:strain_def}
		except:
			strain_info_dict = get_strain_info(file,outdir)
	else:
		type="genbank"
		record = SeqIO.read(file,"genbank")
		strain_id = record.id
		strain_def = record.description
		strain_info_dict = {strain_id.split()[0].split('.')[0]:strain_def}
	return strain_info_dict,type

def get_strain_info(file,outdir):
	#fasta or multi-fasta
	strain_file = os.path.join(outdir,'strain_inf.txt')
	f_result = open(strain_file,'w')
	strain_inf_dict = {}
	with open(file) as f:
		contents = f.read().strip()
	if '\n>' in contents:		
		for strain in contents.split('\n>'):
			strain_title = strain.split('\n')[0].strip()
			strain_id = strain_title.split()[0].strip('>')
			strain_def = strain_title.strip('>')
			strain_inf_dict.update({strain_id.split('.')[0]:strain_def})
			f_result.write(strain_id+'\t'+strain_def+'\n')
			f_result.flush()
	else:
		strain_title = contents.split('\n')[0].strip()
		strain_id = strain_title.split()[0].strip('>')
		strain_def = strain_title.strip('>')
		strain_inf_dict.update({strain_id.split('.')[0]:strain_def})
		f_result.write(strain_id+'\t'+strain_def+'\n')
		f_result.flush()
	f_result.close()
	return strain_inf_dict

def get_prophage_result(prophage_blastp_dict_file,prophage_blastn_dict_file,outdir,kind):
	prophage_blastp_dict = {}
	prophage_blastn_dict = {}
	try:
		prophage_blastp_dict = np.load(prophage_blastp_dict_file).item()
	except:
		pass
	try:
		prophage_blastn_dict = np.load(prophage_blastn_dict_file).item()
	except:
		pass
	strong_result_file = os.path.join(outdir,'prophage_true_result.txt')
	weak_result_file = os.path.join(outdir,'prophage_model_result.txt')
	all_result_file = os.path.join(outdir,'prophage_result.txt')
	f_prophage_result = open(strong_result_file,'w')
	f_prophage_model = open(weak_result_file,'w')
	f_prophage_all = open(all_result_file,'w')
	if kind=='p':
		f_prophage_result.write('phage_id\tphage_def\thost_id\thost_def\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\n')
		f_prophage_model.write('phage_id\tphage_def\thost_id\thost_def\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\n')
		f_prophage_all.write('phage_id\tphage_def\thost_id\thost_def\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\tmode\n')		
	else:
		f_prophage_result.write('bac_id\tbac_def\thit_phage_id\thit_phage_def\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\n')
		f_prophage_model.write('bac_id\tbac_def\thit_phage_id\thit_phage_def\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\n')
		f_prophage_all.write('bac_id\tbac_def\thit_phage_id\thit_phage_def\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\tmode\n')	
	all_query_list = list(set(list(prophage_blastp_dict.keys())+list(prophage_blastn_dict.keys())))
	for query_id in all_query_list:
		if query_id.split('.')[0] in strain_inf_dict.keys():
			query_def = strain_inf_dict[query_id.split('.')[0]]
		else:
			query_def = query_id
		if query_id in prophage_blastp_dict.keys():
			blastp_hit_list = list(prophage_blastp_dict[query_id].keys())
		else:
			blastp_hit_list = []
		if query_id in prophage_blastn_dict.keys():
			blastn_hit_list = list(prophage_blastn_dict[query_id].keys())
		else:
			blastn_hit_list = []
		c_hit_list = list(set(blastp_hit_list+blastn_hit_list))
		for hit_id in c_hit_list:
			query_id = query_id.split('.')[0]
			hit_id = hit_id.split('.')[0]
			try:
				hit_def = phage_inf_dict[hit_id]
			except:
				hit_def = bac_inf_dict[hit_id]
			prophage_pro_percent = 0
			prophage_identity = 0
			prophage_coverage = 0
			host_have_prophage_num = 0
			flag = 0 #this phage need to be sent to model to predict
			if query_id in prophage_blastp_dict.keys():
				if hit_id in prophage_blastp_dict[query_id].keys():
					prophage_info = prophage_blastp_dict[query_id][hit_id]
					hit_pro_percent = 0					
					for method in prophage_info.keys():
						for prophage_region,hit_info in prophage_info[method].items():
							#print(hit_info)
							if float(hit_info[0])>hit_pro_percent:
								hit_pro_percent = float(hit_info[0])
							if float(hit_info[0])>=0.7:
								flag=1
					#print(hit_pro_percent)
					prophage_pro_percent = hit_pro_percent
			if query_id in prophage_blastn_dict.keys():
				if hit_id in prophage_blastn_dict[query_id].keys():
					prophage_info = prophage_blastn_dict[query_id][hit_id]
					hit_identity = 0
					hit_coverage = 0
					for method in prophage_info.keys():
						for region,hit_infos in prophage_info[method].items():
							# print('test')
							# print(hit_infos)
							region_identity = hit_infos[1]
							region_coverage = hit_infos[2]
							if float(region_coverage)>=0.75 and float(region_identity)>=80:
								flag=1
							if hit_identity < region_identity:
								hit_identity = region_identity
							if hit_coverage < region_coverage:
								hit_coverage = region_coverage
					prophage_identity = hit_identity
					prophage_coverage = hit_coverage		
			#print(prophage_pro_percent)
			if flag == 0:
				f_prophage_model.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(prophage_pro_percent)+'\t'+str(prophage_identity)+'\t'+str(prophage_coverage)+'\n')
				f_prophage_model.flush()
				mode = "weak"
			else:
				f_prophage_result.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(prophage_pro_percent)+'\t'+str(prophage_identity)+'\t'+str(prophage_coverage)+'\n')
				f_prophage_result.flush()
				mode = "strong"
			f_prophage_all.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(prophage_pro_percent)+'\t'+str(prophage_identity)+'\t'+str(prophage_coverage)+'\t'+mode+'\n')
			f_prophage_all.flush()
	f_prophage_model.close()
	f_prophage_result.close()
	f_prophage_all.close()

if __name__=='__main__':
	multiprocessing.freeze_support()
	parser = argparse.ArgumentParser()
	parser.add_argument('--input_file',help='Path of the input phage file.\n')
	parser.add_argument('--input_protein',help='optional,Path of the input phage protein file.\n')
	parser.add_argument('--output', help='Path of the output folder.\n')
	parser.add_argument('--type', help='optional,fasta or genbank.\n')
	parser.add_argument('--min_per_prophage', help='Minimal % percentage of hit proteins on hit prophage region(default:30)\n')
	parser.add_argument('--min_id_prophage', help='Minimal % identity of hit region on hit prophage region by making blastn(default:70)\n')
	parser.add_argument('--min_cov_prophage', help='Minimal % coverage of hit region on hit prophage region by making blastn(default:30)\n')
	args = parser.parse_args()
	global type,strain_inf_dict,bac_inf_dict,root_path,min_per_prophage,min_id_prophage,min_cov_prophage
	if args.input_file:
		input_file = args.input_file
	flag = 'no'
	if args.input_protein:
		pro_file = args.input_protein
		if not os.path.exists(pro_file):
			flag = 'yes'
	else:
		flag = 'yes'
	if args.output:
		outdir = args.output
	if args.min_per_prophage:
		min_per_prophage = args.min_per_prophage
	if args.min_id_prophage:
		min_id_prophage = args.min_id_prophage
	if args.min_cov_prophage:
		min_cov_prophage = args.min_cov_prophage
	root_path = get_root_path()
	file_name = os.path.basename(input_file).split()[0].split('.')[0]
	strain_inf_dict,type = get_inf(input_file,outdir)
	if args.type:
		type = args.type
	if type=='genbank':
		global input_gb_file
		input_gb_file = input_file
	if flag == 'yes':
		faa_prefix = os.path.join(outdir,file_name)
		pro_file = faa_prefix+'.faa'
		if type == 'fasta':
			pred_orf(input_file,faa_prefix)
		else:
			getFaaFromGB(input_file,outdir)
			getFnaFromGB(input_file,outdir,file_name)
			input_file = os.path.join(outdir,file_name+'.fna')
	if args.min_per_prophage:
		min_per_prophage = args.min_per_prophage
	else:
		min_per_prophage = 30
	if args.min_id_prophage:
		min_id_prophage = args.min_id_prophage
	else:
		min_id_prophage = 70
	if args.min_cov_prophage:
		min_cov_prophage = args.min_cov_prophage
	else:
		min_cov_prophage = 30
	bac_inf_dict_file = os.path.join(root_path,'database/db/profile','bac_inf_dict.npy')
	bac_inf_dict = np.load(bac_inf_dict_file).item()
	phage_to_hosts(input_file,pro_file,outdir)
	print('finished predicting hosts in %s!'%outdir)