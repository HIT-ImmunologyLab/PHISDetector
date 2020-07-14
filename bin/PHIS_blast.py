import os
import os,threading
import numpy.random.common
import numpy.random.bounded_integers
import numpy.random.entropy
import numpy as np
import os,re,argparse
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

def pred_orf(fasta_file,faa_prefix):
	script = os.path.join(root_path,'software','FragGeneScan','run_FragGeneScan.pl')
	cmd_fragGeneScan = script+' -genome %s -out %s -complete=1 -train=complete -thread=20' % (fasta_file, faa_prefix)
	# print(cmd_fragGeneScan)
	os.system(cmd_fragGeneScan)

def getFaaFromGB(input_file,outdir,strain_id = 'no'):   #parse protein from genbank files in phaster
	special_pros = ['capsid','head','plate','tail','coat','portal','holin','integrase','transposase','terminase','protease','lysis','bacteriocin','tRNA']
	records = SeqIO.parse(input_file, "gb")
	file_name = os.path.basename(input_file).split()[0].split('.')[0]
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
						strain_id = list(strain_info_dict.keys())[0]
						savefile.write('>' +strain_id+'|' + str(location)+ '|' + str(proteinId)+'\n')
						if translation[-1] == '\n':
							savefile.write(translation)
						else:
							savefile.write(translation + '\n')
	savefile.close()
	savefile_protein.close()

def blastn_long(file,outfile,database,format,evalue):
	blastn_path = os.path.join(root_path,'software/blast+','blastn')
	num_threads = 30
	format = 6
	script = blastn_path+" -db "+database+" -query "+file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile+" -num_threads "+str(num_threads)+" -word_size 11 -reward 1 -penalty -2 -gapopen 0 -gapextend 0 -perc_identity 70 -max_target_seqs 100000"
	os.system(script)

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

def diamond_blastp_pro_int(file,outfile,database,format,evalue,identity,coverage):
	num_threads = 20
	diamond_path = os.path.join(root_path,'software','diamond','diamond')
	script = diamond_path+" blastp -d "+database+" -q "+file+" -f "+str(format)+" -e "+str(evalue)+" -o "+outfile+" -p "+str(num_threads)+" --id "+str(identity)+" -k 10000000"
	os.system(script)

def blastp(query_file,subject_file,outfile):
	format = 6
	evalue = 10
	blastp_path = os.path.join(root_path,'software/blast+','blastp')
	script = blastp_path+" -db "+subject_file+" -query "+query_file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile
	os.system(script)

def blastn_file(query_file,subject_file,outfile):
	format = 6
	evalue = 1
	blastn_path = os.path.join(root_path,'software/blast+','blastn')
	script = blastn_path+" -subject "+subject_file+" -query "+query_file+" -outfmt "+str(format)+" -evalue "+str(evalue)+" -out "+outfile
	os.system(script)

def getFnaFromGB(fileName,outDir,fileID):
	outFileName = os.path.join(outDir,fileID + '.fna')
	handle = open(fileName)
	if os.path.exists(outFileName):
		os.remove(outFileName)
	SeqIO.convert(handle, 'genbank', outFileName, 'fasta')

def parse_phage_blastp_hosts(pro_file,blastp_file,outdir):
	pro_num_dict = get_pro_num(pro_file,outdir)
	phage_blastp_bacs_result_dict = {}
	resufile = os.path.join(outdir,'phage_blastp_bac_db_result.txt')
	f_result = open(resufile,'w')
	f_result.write('phage_id\tphage_def\tbac_id\tbac_def\tphage_protein_number\tbac_protein_number\thomo_protein_pair_number\n')
	with open(blastp_file) as f:
		contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		if type == 'fasta':
			phage_id = '_'.join(line[0].split('_')[0:-3]).split('.')[0]
		else:
			phage_id = line[0].split('|')[0].split('.')[0]
		if phage_id not in phage_blastp_bacs_result_dict.keys():
			phage_blastp_bacs_result_dict.update({phage_id:{}})
		phage_pro = line[0]
		phage_pro_num = pro_num_dict[phage_id]
		hit_pros = line[1].split('|')
		for hit_bac_pro in hit_pros[0::2]:			
			hit_bac_id = '_'.join(hit_bac_pro.split('_')[0:-3]).split('.')[0]
			hit_bac_pro_num = hit_pros[hit_pros.index(hit_bac_pro)+1]
			if hit_bac_id not in phage_blastp_bacs_result_dict[phage_id].keys():
				phage_blastp_bacs_result_dict[phage_id].update({hit_bac_id:[]})
			phage_blastp_bacs_result_dict[phage_id][hit_bac_id].append([phage_pro,hit_bac_pro,phage_pro_num,hit_bac_pro_num])	
	
	phage_blastp_bacs_result_dict1 = {}	
	for phage_id in phage_blastp_bacs_result_dict.keys():
		try:
			phage_def = strain_info_dict[phage_id]
		except:
			phage_def = phage_id
		phage_pro_num = pro_num_dict[phage_id]
		for hit_bac_id in phage_blastp_bacs_result_dict[phage_id].keys():		
			try:
				hit_bac_def = bac_inf_dict[hit_bac_id]
			except:
				hit_bac_def = hit_bac_id
			homo_pro_num = len(list(set([item[0] for item in phage_blastp_bacs_result_dict[phage_id][hit_bac_id]])))
			if homo_pro_num/float(phage_pro_num)>=float(min_per_blast)/100:
				if phage_id not in phage_blastp_bacs_result_dict1.keys():
					phage_blastp_bacs_result_dict1.update({phage_id:{}})
				if hit_bac_id not in phage_blastp_bacs_result_dict1[phage_id].keys():
					phage_blastp_bacs_result_dict1[phage_id].update({hit_bac_id:phage_blastp_bacs_result_dict[phage_id][hit_bac_id]})
				f_result.write(phage_id+'\t'+phage_def+'\t'+hit_bac_id+'\t'+hit_bac_def+'\t'+str(phage_pro_num)+'\t'+str(hit_bac_pro_num)+'\t'+str(homo_pro_num)+'\n')
				f_result.flush()
	f_result.close()
	outfile = os.path.join(outdir,'phage_blastp_bac_db_result_dict')
	np.save(outfile,phage_blastp_bacs_result_dict1)

def parse_phage_blastn_hosts(file,outdir):
	phages_length_dict_file = os.path.join(outdir,'phages_length_dict.npy')
	phages_length_dict = np.load(phages_length_dict_file).item()
	phage_blastn_bacs_result_dict = {}
	filter_file = os.path.join(outdir,'phage_blastn_bac_db_result.txt')
	f_filter = open(filter_file,'w')
	f_filter.write('phage_id\tphage_def\thit_bac_id\thit_bac_def\tphage_length\thit_length\tidentity\tcoverage\n')
	with open(file) as f:
		contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		phage_id = line[0].split('.')[0]
		if phage_id in phages_length_dict.keys():			
			phage_length = phages_length_dict[phage_id]
		else:
			phage_length = 0
		hit_bac_id = line[1].split('.')[0]
		hit_phage_start = line[6]
		hit_phage_end = line[7]
		mismatch = int(line[4])
		identity = float(line[2])
		hit_length = float(line[3])
		bit_score = line[-1]
		if identity>=float(min_id_blast):
			if phage_id not in phage_blastn_bacs_result_dict.keys():
				phage_blastn_bacs_result_dict.update({phage_id:{}})
			if hit_bac_id not in phage_blastn_bacs_result_dict[phage_id].keys():
				phage_blastn_bacs_result_dict[phage_id].update({hit_bac_id:[]})
			phage_blastn_bacs_result_dict[phage_id][hit_bac_id].append([identity,[hit_phage_start,hit_phage_end],bit_score,phage_length])
	phage_blastn_bacs_result_dict1 = {}
	for phage_id in phage_blastn_bacs_result_dict.keys():
		try:
			phage_def = strain_info_dict[phage_id.split('.')[0]]
		except:
			phage_def = phage_id
		for hit_bac_id in phage_blastn_bacs_result_dict[phage_id].keys():
			try:
				hit_bac_def = bac_inf_dict[hit_bac_id.split('.')[0]]
			except:
				hit_bac_def = hit_bac_idd
			hit_region = []
			identitys = []
			for hit_record in phage_blastn_bacs_result_dict[phage_id][hit_bac_id]:
				identity = hit_record[0]
				identitys.append(identity)
				phage_hit_region = hit_record[1]
				hit_region.append(phage_hit_region)
				phage_length = hit_record[-1]
			hit_length = 0
			if len(hit_region)>0:
				merge_region = combine(hit_region)
				for region in merge_region:
					region_length = abs(int(region[0])-int(region[1]))+1
					hit_length = hit_length+region_length
				if phage_length>0:
					coverage = float(hit_length)/phage_length
				else:
					coverage = 0
				identity = np.mean(list(map(float,identitys)))
			else:
				coverage = 0
				identity = 0
				if phage_id.split('.')[0] in phages_length_dict.keys():			
					phage_length = phages_length_dict[phage_id]
				else:
					phage_length = 0
			#cutoff
			if coverage>=float(min_cov_blast)/100:
				if phage_id not in phage_blastn_bacs_result_dict1.keys():
					phage_blastn_bacs_result_dict1.update({phage_id:{}})
				if hit_bac_id not in phage_blastn_bacs_result_dict1[phage_id].keys():
					phage_blastn_bacs_result_dict1[phage_id].update({hit_bac_id:[]})
				f_filter.write('\t'.join([phage_id,phage_def,hit_bac_id,hit_bac_def,str(phage_length),str(hit_length),str(identity),str(coverage)])+'\n')
				phage_blastn_bacs_result_dict1[phage_id][hit_bac_id] = [identity,coverage,phage_length]
	f_filter.close()
	outfile = os.path.join(outdir,'phage_blastn_bac_db_result_dict')
	np.save(outfile,phage_blastn_bacs_result_dict1)

def get_pro_num(pro_file,outdir):
	pro_num_dict = {}
	with open(pro_file) as f:
		pros = f.read().strip().split('>')
	for pro in pros[1:]:
		pro_title = pro.split('\n')[0]
		if type == 'fasta':
			strain_id = '_'.join(pro_title.split('_')[0:-3]).split('.')[0]
		else:
			strain_id = pro_title.split('|')[0].split('.')[0]
		if strain_id not in pro_num_dict.keys():			
			pro_num_dict.update({strain_id:0})
		pro_num_dict[strain_id] = pro_num_dict[strain_id]+1	
	save_file = os.path.join(outdir,'protein_num_dict')
	np.save(save_file,pro_num_dict)
	return pro_num_dict

def phage_to_hosts_blastp(input_file,pro_file,outdir):
	database_blastp = os.path.join(root_path,'database/db/database/bacteria_protein_diamond_nonredundant','bacteria_protein_db')
	outfile = os.path.join(outdir,'phage_blastp_bac_db')
	diamond_blastp_pro_int(pro_file,outfile,database_blastp,6,1,0.4,0.7)
	filter_file = os.path.join(outdir,'phage_blastp_bac_db_identity_0.4_coverage_0.7')
	filter_identity_coverage(outfile,filter_file)
	parse_phage_blastp_hosts(pro_file,filter_file,outdir)

def phage_to_hosts_blastn(input_file,outdir):
	database_blastn = os.path.join(root_path,'database/db/database/bacteria_nucl','bacteria_nucl_db')
	outfile = os.path.join(outdir,'phage_blastn_bac_db')
	blastn_long(input_file,outfile,database_blastn,6,0.01)
	parse_phage_blastn_hosts(outfile,outdir)

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

def filter_identity_coverage(blastp_file,filter_file):
	f_result = open(filter_file,'w')
	with open(blastp_file) as f:
		contents = f.readlines()
	for line in contents:
		line = line.strip().split('\t')
		query_pro = line[0]
		hit_pro = line[1]
		identity = line[2]
		query_start = line[-6]
		query_end = line[-5]
		hit_length = abs(int(query_start)-int(query_end))+1
		if type=='fasta':
			pro_length = int((abs(int(query_pro.split('_')[-3])-int(query_pro.split('_')[-2]))+1)/3)
		else:
			pro_length = int((abs(int(query_pro.split('|')[1].split('_')[-3])-int(query_pro.split('|')[1].split('_')[-2]))+1)/3)
		coverage = float(hit_length)/pro_length
		if (coverage>=0.7) and ((float(identity)/100)>=0.4):
			f_result.write('\t'.join(line)+'\t'+str(coverage)+'\n')
			f_result.flush()
	f_result.close()

def get_length(file,outdir):
	#fasta or multi-fasta
	phages_length_dict = {}
	with open(file) as f:
		contents = f.read().strip()
	if '\n>' in contents:
		for phage in contents.split('\n>'):
			phage_id = phage.strip().split('\n')[0].split()[0].split('.')[0].strip('>')		
			phage_sequence = ''
			for sequence in phage.strip().split('\n')[1:]:
				phage_sequence = sequence.strip()+phage_sequence
			phage_length = len(phage_sequence)
			phages_length_dict.update({phage_id:phage_length})
	else:
		phage_id = contents.split('\n')[0].split()[0].split('.')[0].strip('>')
		phage_sequence = ''
		for sequence in contents.strip().split('\n')[1:]:
			phage_sequence = sequence.strip()+phage_sequence
		phage_length = len(phage_sequence)
		phages_length_dict.update({phage_id:phage_length})
	outfile_dict = os.path.join(outdir,'phages_length_dict')
	np.save(outfile_dict,phages_length_dict)
	return phages_length_dict

def get_blast_features(blastp_result_dict_file,blastn_result_dict_file,outdir,kind):	
	blastp_result_dict = {}
	blastn_result_dict = {}
	try:
		blastp_result_dict = np.load(blastp_result_dict_file).item()
	except:
		pass
	try:
		blastn_result_dict = np.load(blastn_result_dict_file).item()
	except:
		pass
	strong_result_file = os.path.join(outdir,'blast_true_result.txt')
	weak_result_file = os.path.join(outdir,'blast_model_result.txt')
	all_result_file = os.path.join(outdir,'blast_result.txt')
	f_blast_result = open(strong_result_file,'w')
	f_blast_model = open(weak_result_file,'w')
	f_blast_all = open(all_result_file,'w')
	if kind=='p':
		f_blast_result.write('phage_id\tphage_def\thost_id\thost_def\tblastp_score\tblastn_identity\tblastn_coverage\n')
		f_blast_model.write('phage_id\tphage_def\thost_id\thost_def\tblastp_score\tblastn_identity\tblastn_coverage\n')
		f_blast_all.write('phage_id\tphage_def\thost_id\thost_def\tblastp_score\tblastn_identity\tblastn_coverage\tmode\n')
	else:
		f_blast_result.write('bac_id\tbac_def\tphage_id\tphage_def\tblastp_score\tblastn_identity\tblastn_coverage\n')
		f_blast_model.write('bac_id\tbac_def\tphage_id\tphage_def\tblastp_score\tblastn_identity\tblastn_coverage\n')
		f_blast_all.write('bac_id\tbac_def\tphage_id\tphage_def\tblastp_score\tblastn_identity\tblastn_coverage\tmode\n')
	
	all_query_list = list(set(list(blastp_result_dict.keys())+list(blastn_result_dict.keys())))
	for query_id in all_query_list:
		if query_id.split('.')[0] in strain_info_dict.keys():
			query_def = strain_info_dict[query_id.split('.')[0]]
		else:
			query_def = query_id
		
		if query_id in blastp_result_dict.keys():
			blastp_hit_list = list(blastp_result_dict[query_id].keys())
		else:
			blastp_hit_list = []
		if query_id in blastn_result_dict.keys():
			blastn_hit_list = list(blastn_result_dict[query_id].keys())
		else:
			blastn_hit_list = []
		c_hit_list = list(set(blastp_hit_list+blastn_hit_list))
		for hit_id in c_hit_list:
			hit_id = hit_id.split('.')[0]
			try:
				hit_def = phage_inf_dict[hit_id]
			except:
				hit_def = bac_inf_dict[hit_id]
			blastp_score = 0
			blastn_identity = 0
			blastn_coverage = 0
			flag = 0#this phage need to be sent to model to predict
			if query_id in blastp_result_dict.keys():
				if hit_id in blastp_result_dict[query_id].keys():
					blastp_info = blastp_result_dict[query_id][hit_id]
					phage_pro_num = blastp_info[0][2]
					bac_pro_num = blastp_info[0][3]
					homo_pro_num = len(blastp_info)
					phage_int_pros = list(set([item[0] for item in blastp_info]))
					hit_pro_percent = float(homo_pro_num)/(float(phage_pro_num)+float(bac_pro_num)-homo_pro_num)
					blastp_score = hit_pro_percent
					if float(len(phage_int_pros))/float(phage_pro_num)>=0.7:
						flag = 1
			if query_id in blastn_result_dict.keys():
				if hit_id in blastn_result_dict[query_id].keys():
					blastn_info = blastn_result_dict[query_id][hit_id]
					blastn_identity = blastn_info[0]
					blastn_coverage = blastn_info[1]
					if (float(blastn_identity)>=80) and (float(blastn_coverage) >= 0.75):
						flag = 1
			if flag==0:
				f_blast_model.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(blastp_score)+'\t'+str(blastn_identity)+'\t'+str(blastn_coverage)+'\n')
				f_blast_model.flush()
				mode = "weak"
			else:
				f_blast_result.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(blastp_score)+'\t'+str(blastn_identity)+'\t'+str(blastn_coverage)+'\n')
				f_blast_result.flush()
				mode = "strong"
			f_blast_all.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(blastp_score)+'\t'+str(blastn_identity)+'\t'+str(blastn_coverage)+'\t'+mode+'\n')
			f_blast_all.flush()
	f_blast_model.close()
	f_blast_result.close()
	f_blast_all.close()

if __name__=='__main__':
	multiprocessing.freeze_support()
	parser = argparse.ArgumentParser()
	parser.add_argument('--input_file',help='Path of the input file.\n')
	parser.add_argument('--input_protein',help='Path of the input file.\n')
	parser.add_argument('--output', help='Path of the output file.\n')
	parser.add_argument('--type', help='Path of the output file.\n')
	parser.add_argument('--min_per_blast', help='Minimal % percentage of hit proteins on hit prophage region(default:10)\n')
	parser.add_argument('--min_id_blast', help='Minimal % identity of hit region on hit prophage region by making blastn(default:70)\n')
	parser.add_argument('--min_cov_blast', help='Minimal % coverage of hit region on hit prophage region by making blastn(default:10)\n')

	args = parser.parse_args()
	global type,strain_info_dict,bac_inf_dict,root_path,min_per_blast,min_id_blast,min_cov_blast
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
	if args.min_per_blast:
		min_per_blast = args.min_per_blast
	else:
		min_per_blast = 10
	if args.min_id_blast:
		min_id_blast = args.min_id_blast
	else:
		min_id_blast = 70
	if args.min_cov_blast:
		min_cov_blast = args.min_cov_blast
	else:
		min_cov_blast = 10
	root_path = get_root_path()
	strain_info_dict,type = get_inf(input_file,outdir)
	file_name = os.path.basename(input_file).split()[0].split('.')[0]
	if args.type:
		type = args.type
		strain_id = file_name
		strain_def = file_name
	else:
		strain_info_dict,type = get_inf(input_file,outdir)
	if flag == 'yes':
		faa_prefix = os.path.join(outdir,file_name)
		pro_file = faa_prefix+'.faa'
		if type == 'fasta':
			pred_orf(input_file,faa_prefix)
		else:
			getFaaFromGB(input_file,outdir)
			getFnaFromGB(input_file,outdir,file_name)
			input_file = os.path.join(outdir,file_name+'.fna')
	bac_inf_dict_file = os.path.join(root_path,'database/db/profile','bac_inf_dict.npy')
	bac_inf_dict = np.load(bac_inf_dict_file).item()	
	get_length(input_file,outdir)
	m1 = threading.Thread(target=phage_to_hosts_blastp,args=(input_file,pro_file,outdir,))
	m2 = threading.Thread(target=phage_to_hosts_blastn,args=(input_file,outdir,))
	m1.start()
	m2.start()
	m1.join()
	m2.join()
	blastp_result_dict_file = os.path.join(outdir,'phage_blastp_bac_db_result_dict.npy')
	blastn_result_dict_file = os.path.join(outdir,'phage_blastn_bac_db_result_dict.npy')
	get_blast_features(blastp_result_dict_file,blastn_result_dict_file,outdir,'p')
	print('finished predicting hosts in %s'%outdir)