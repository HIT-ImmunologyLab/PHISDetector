import os,threading
import numpy.random.common
import numpy.random.bounded_integers
import numpy.random.entropy
import numpy as np
import os,re,argparse
from Bio import Entrez, SeqIO
from collections import OrderedDict
import sys
import joblib
import multiprocessing
import sklearn

#get PHIS directory
def get_root_path():
	#script_path = os.path.split(os.path.realpath(__file__))[0]
	# script_path = os.path.split(os.path.realpath(sys.argv[0]))[0]
	# root_path = os.path.dirname(script_path)
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

def getFnaFromGB(fileName,outDir,fileID):
	outFileName = outDir + '/' + fileID + '.fna'
	handle = open(fileName)
	if os.path.exists(outFileName):
		os.remove(outFileName)
	SeqIO.convert(handle, 'genbank', outFileName, 'fasta')

def pred_orf(fasta_file,faa_prefix):
	script = os.path.join(root_path,'software','FragGeneScan','run_FragGeneScan.pl')
	cmd_fragGeneScan = script+' -genome %s -out %s -complete=1 -train=complete -thread=20' % (fasta_file, faa_prefix)
	print(cmd_fragGeneScan)
	os.system(cmd_fragGeneScan)

def getFaaFromGB(input_file,outdir):#parse protein from genbank files in phaster
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
						strain_id = list(strain_inf_dict.keys())[0]
						savefile.write('>' +strain_id+'|' + str(location)+ '|' + str(proteinId)+'\n')					
					# savefile.write(">"+fileID+ '_'+str(counter)+'\n')                   
						if translation[-1] == '\n':
							savefile.write(translation)
						else:
							savefile.write(translation + '\n')
	savefile.close()
	savefile_protein.close()

def predict_on_crispr(file,outdir,min_mis_crispr=2,min_cov_crispr=70):
	script = os.path.join(root_path,'bin','PHIS_crispr')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = script+" --input "+file+" --output "+outdir+" --min_mis_crispr "+str(min_mis_crispr)+" --min_cov_crispr "+str(min_cov_crispr)
	os.system(command)
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)

def predict_on_prophage(file,outdir,min_per_prophage=30,min_id_prophage=70,min_cov_prophage=30):
	script = os.path.join(root_path,'bin','PHIS_prophage')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = script+" --input_file "+file+" --output "+outdir+" --min_per_prophage "+str(min_per_prophage)+" --min_id_prophage "+str(min_id_prophage)+" --min_cov_prophage "+str(min_cov_prophage)
	os.system(command)
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)

def predict_on_blast(file,outdir,min_per_blast=10,min_id_blast=70,min_cov_blast=10):
	script = os.path.join(root_path,'bin','PHIS_blast')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = script+" --input_file "+file+" --output "+outdir+" --min_per_blast "+str(min_per_blast)+" --min_id_blast "+str(min_id_blast)+" --min_cov_blast "+str(min_cov_blast)
	os.system(command)
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)

def predict_on_pro_int(file,outdir,min_PPI=1,min_DDI=5):
	script = os.path.join(root_path,'bin','PHIS_protein_protein_int')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = script+" --input_file "+file+" --output "+outdir+" --min_PPI "+str(min_PPI)+" --min_DDI "+str(min_DDI)
	os.system(command)
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)

def run_crispr_prophage_blast(input_file,module_dir,kind):
	#crispr
	outdir_crispr = os.path.join(module_dir,'crispr')
	mkdir(outdir_crispr)
	m1 = threading.Thread(target=predict_on_crispr,args=(input_file,outdir_crispr,min_mis_crispr,min_cov_crispr))
	#prophage
	outdir_prophage = os.path.join(module_dir,'prophage')
	mkdir(outdir_prophage)
	m2 = threading.Thread(target = predict_on_prophage,args=(input_file,outdir_prophage,min_per_prophage,min_id_prophage,min_cov_prophage))
	#blast
	outdir_blast = os.path.join(module_dir,'blast')
	mkdir(outdir_blast)
	m3 = threading.Thread(target=predict_on_blast,args=(input_file,outdir_blast,min_per_blast,min_id_blast,min_cov_blast))
	m1.start()
	m2.start()
	m3.start()
	m1.join()
	m2.join()
	m3.join()

def run_genehomo_pro_int(input_file,pro_file,phage_hit_hosts_now,outdir,fasta_dir,kind):
	outdir_virhostmatcher = os.path.join(outdir,'genehomo','virhost')
	outdir_wish = os.path.join(outdir,'genehomo','wish')
	outdir_codon = os.path.join(outdir,'genehomo','codon')
	outdir_pro_int = os.path.join(outdir,'pro_pro_int')
	mkdir(outdir_virhostmatcher)
	mkdir(outdir_wish)
	mkdir(outdir_codon)
	mkdir(outdir_pro_int)
	m1 = threading.Thread(target=predict_on_virhostmatcher,args=(fasta_dir,phage_hit_hosts_now,outdir_virhostmatcher))
	m2 = threading.Thread(target=predict_on_wish,args=(fasta_dir,phage_hit_hosts_now,outdir_wish))
	m3 = threading.Thread(target=predict_on_codon,args=(input_file,phage_hit_hosts_now,outdir_codon))
	m4 = threading.Thread(target = predict_on_pro_int,args=(input_file,outdir_pro_int,min_PPI,min_DDI))
	m1.start()
	m2.start()
	m3.start()
	m4.start()
	m1.join()
	m2.join()
	m3.join()
	m4.join()

def mkdir(dirname):
	command = "mkdir -p "+dirname
	os.system(command)

def collect_result_1(input_file,pro_file,outdir,kind):
	outdir_crispr = os.path.join(outdir,'crispr')
	outdir_prophage = os.path.join(outdir,'prophage')
	outdir_pro_int = os.path.join(outdir,'pro_pro_int')
	outdir_blast = os.path.join(outdir,'blast')
	if kind == 'p':
		record_file = os.path.join(outdir,'phage_hit_hosts_crispr_prophage_pro_int_blast_dict')
		crispr_file = os.path.join(outdir_crispr,'hit_bac_spacers_dict.npy')
		prophage_file_blastp = os.path.join(outdir_prophage,'phage_blastp_prophage_result_dict.npy')
		prophage_file_blastn = os.path.join(outdir_prophage,'phage_blastn_prophage_result_dict.npy')
		pro_int_file_homolog = os.path.join(outdir_pro_int,'phage_int_bac_homo_result_dict.npy')
		pro_int_file_domain = os.path.join(outdir_pro_int,'phage_int_bac_domain_result_dict.npy')
		blastp_file = os.path.join(outdir_blast,'phage_blastp_bac_db_result_dict.npy')
		blastn_file = os.path.join(outdir_blast,'phage_blastn_bac_db_result_dict.npy')
		if os.path.exists(crispr_file):
			hit_bac_spacers_dict = np.load(crispr_file).item()
		else:
			hit_bac_spacers_dict = {}
		if os.path.exists(prophage_file_blastp):
			phage_blastp_prophage_result_dict = np.load(prophage_file_blastp).item()
		else:
			phage_blastp_prophage_result_dict = {}
		if os.path.exists(prophage_file_blastn):
			phage_blastn_prophage_result_dict = np.load(prophage_file_blastn).item()
		else:
			phage_blastn_prophage_result_dict = {}
		if os.path.exists(pro_int_file_homolog):
			phage_int_bac_homo_result_dict = np.load(pro_int_file_homolog).item()
		else:
			phage_int_bac_homo_result_dict = {}
		if os.path.exists(pro_int_file_domain):
			phage_int_bac_domain_result_dict = np.load(pro_int_file_domain).item()
		else:
			phage_int_bac_domain_result_dict = {}
		if os.path.exists(blastp_file):
			phage_blastp_bac_db_result_dict = np.load(blastp_file).item()
		else:
			phage_blastp_bac_db_result_dict = {}
		if os.path.exists(blastn_file):
			phage_blastn_bac_db_result_dict = np.load(blastn_file).item()
		else:
			phage_blastn_bac_db_result_dict = {}
		result_dict_list = [hit_bac_spacers_dict,phage_blastp_prophage_result_dict,phage_blastn_prophage_result_dict,phage_int_bac_homo_result_dict,phage_int_bac_domain_result_dict,phage_blastp_bac_db_result_dict,phage_blastn_bac_db_result_dict]
		phage_hit_bacs_dict = OrderedDict()
		for result_dict in result_dict_list:
			for phage_id in result_dict.keys():
				phage_id = phage_id.split('.')[0]
				if phage_id not in phage_hit_bacs_dict.keys():
					phage_hit_bacs_dict.update({phage_id:[]})
				phage_hit_bacs_dict[phage_id] = phage_hit_bacs_dict[phage_id]+list(result_dict[phage_id].keys())
				phage_hit_bacs_dict[phage_id] = list(set(phage_hit_bacs_dict[phage_id]))
		np.save(record_file,phage_hit_bacs_dict)
		return result_dict_list
	else:
		record_file = os.path.join(outdir,'bac_hit_phages_crispr_prophage_pro_int_blast_dict')
		crispr_file = os.path.join(outdir_crispr,'hit_phages_dict.npy')
		prophage_file_blastp = os.path.join(outdir_prophage,'prophage_blastp_phagedb_result_dict.npy')
		prophage_file_blastn = os.path.join(outdir_prophage,'prophage_blastn_phagedb_result_dict.npy')
		pro_int_file_homolog = os.path.join(outdir_pro_int,'bac_int_phage_homo_result_dict.npy')
		pro_int_file_domain = os.path.join(outdir_pro_int,'bac_int_phage_domain_result_dict.npy')
		blastp_file = os.path.join(outdir_blast,'bac_blastp_phage_db_result_dict.npy')
		blastn_file = os.path.join(outdir_blast,'bac_blastn_phage_db_result_dict.npy')
		if os.path.exists(crispr_file):
			hit_bac_spacers_dict = np.load(crispr_file).item()
		else:
			hit_bac_spacers_dict = {}
		if os.path.exists(prophage_file_blastp):
			prophage_blastp_phage_result_dict = np.load(prophage_file_blastp).item()
		else:
			prophage_blastp_phage_result_dict = {}
		if os.path.exists(prophage_file_blastn):
			prophage_blastn_phage_result_dict = np.load(prophage_file_blastn).item()
		else:
			prophage_blastn_phage_result_dict = {}
		if os.path.exists(pro_int_file_homolog):
			bac_int_phage_homo_result_dict = np.load(pro_int_file_homolog).item()
		else:
			bac_int_phage_homo_result_dict = {}
		if os.path.exists(pro_int_file_domain):
			bac_int_phage_domain_result_dict = np.load(pro_int_file_domain).item()
		else:
			bac_int_phage_domain_result_dict = {}
		if os.path.exists(blastp_file):
			bac_blastp_phage_db_result_dict = np.load(blastp_file).item()
		else:
			bac_blastp_phage_db_result_dict = {}
		if os.path.exists(blastn_file):
			bac_blastn_phage_db_result_dict = np.load(blastn_file).item()
		else:
			bac_blastn_phage_db_result_dict = {}
		result_dict_list = [hit_bac_spacers_dict,prophage_blastp_phage_result_dict,prophage_blastn_phage_result_dict,bac_int_phage_homo_result_dict,bac_int_phage_domain_result_dict,bac_blastp_phage_db_result_dict,bac_blastn_phage_db_result_dict]
		bac_hit_phages_dict = OrderedDict()
		for result_dict in result_dict_list:
			for bac_id in result_dict.keys():
				bac_id = bac_id.split('.')[0]
				if bac_id not in bac_hit_phages_dict.keys():
					bac_hit_phages_dict.update({bac_id:[]})
				bac_hit_phages_dict[bac_id] = bac_hit_phages_dict[bac_id]+list(result_dict[bac_id].keys())
				bac_hit_phages_dict[bac_id] = list(set(bac_hit_phages_dict[bac_id]))
		np.save(record_file,bac_hit_phages_dict)
		return result_dict_list

def get_pro_num(pro_file,outdir):
	pro_num_dict = {}	
	with open(pro_file) as f:
		pros = f.read().strip().split('>')
	for pro in pros:
		pro_title = pro.split('\n')[0]
		if type=='fasta':
			strain_id = '_'.join(pro_title.split('_')[0:-3]).split('.')[0]
		else:
			strain_id = pro_title.split('|')[0].split('.')[0]
		
		if strain_id not in pro_num_dict.keys():			
			pro_num_dict.update({strain_id:0})
		pro_num_dict[strain_id] = pro_num_dict[strain_id]+1	
	save_file = os.path.join(outdir,'protein_num_dict')
	np.save(save_file,pro_num_dict)
	return pro_num_dict

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
	strain_inf_dict_file = os.path.join(outdir,'strain_inf_dict')
	np.save(strain_inf_dict_file,strain_info_dict)
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

def collect_result_1(input_file,pro_file,outdir,kind):
	outdir_crispr = os.path.join(outdir,'crispr')
	outdir_prophage = os.path.join(outdir,'prophage')
	outdir_blast = os.path.join(outdir,'blast')
	
	crispr_model_result_file = os.path.join(outdir_crispr,'crispr_model_result.txt')
	prophage_model_result_file = os.path.join(outdir_prophage,'prophage_model_result.txt')
	blast_model_result_file = os.path.join(outdir_blast,'blast_model_result.txt')
	crispr_feature_dict = get_feature_value_dict(crispr_model_result_file)
	prophage_feature_dict = get_feature_value_dict(prophage_model_result_file)
	blast_feature_dict = get_feature_value_dict(blast_model_result_file)
	
	crispr_true_file = os.path.join(outdir_crispr,'crispr_true_result.txt')
	prophage_true_file = os.path.join(outdir_prophage,'prophage_true_result.txt')
	blast_true_file = os.path.join(outdir_blast,'blast_true_result.txt')
	
	crispr_true_dict = get_feature_value_dict(crispr_true_file)
	prophage_true_dict = get_feature_value_dict(prophage_true_file)
	blast_true_dict = get_feature_value_dict(blast_true_file)
	true_result_dict_list = [crispr_true_dict,prophage_true_dict,blast_true_dict]
	if kind == 'p':
		record_file = os.path.join(outdir,'phage_hit_hosts_crispr_prophage_blast_dict')
		result_dict_list = [crispr_feature_dict,prophage_feature_dict,blast_feature_dict]
		phage_hit_bacs_dict = OrderedDict()
		for result_dict in result_dict_list:
			for phage_id in result_dict.keys():
				phage_id = phage_id.split('.')[0]
				if phage_id not in phage_hit_bacs_dict.keys():
					phage_hit_bacs_dict.update({phage_id:[]})
				phage_hit_bacs_dict[phage_id] = phage_hit_bacs_dict[phage_id]+list(result_dict[phage_id].keys())
				phage_hit_bacs_dict[phage_id] = list(set(phage_hit_bacs_dict[phage_id]))
		np.save(record_file,phage_hit_bacs_dict)
		
		feature_dict_now_all = phage_hit_bacs_dict
		for result_dict in true_result_dict_list:
			for phage_id in result_dict.keys():
				phage_id = phage_id.split('.')[0]
				if phage_id not in feature_dict_now_all.keys():
					phage_hit_bacs_dict.update({phage_id:[]})
				feature_dict_now_all[phage_id] = feature_dict_now_all[phage_id]+list(result_dict[phage_id].keys())
				feature_dict_now_all[phage_id] = list(set(feature_dict_now_all[phage_id]))
				feature_dict_now_all[phage_id] = list(set(feature_dict_now_all[phage_id]).intersection(set(list(bac_inf_dict.keys()))))
		record_file_all = os.path.join(outdir,'phage_hit_hosts_crispr_prophage_blast_all_dict')		
		np.save(record_file_all,feature_dict_now_all)
		return result_dict_list,feature_dict_now_all
	else:
		record_file = os.path.join(outdir,'bac_hit_phages_crispr_prophage_blast_dict')
		result_dict_list = [crispr_feature_dict,prophage_feature_dict,blast_feature_dict]
		bac_hit_phages_dict = OrderedDict()
		for result_dict in result_dict_list:
			for bac_id in result_dict.keys():
				bac_id = bac_id.split('.')[0]
				if bac_id not in bac_hit_phages_dict.keys():
					bac_hit_phages_dict.update({bac_id:[]})
				bac_hit_phages_dict[bac_id] = bac_hit_phages_dict[bac_id]+list(result_dict[bac_id].keys())
				bac_hit_phages_dict[bac_id] = list(set(bac_hit_phages_dict[bac_id]))
		np.save(record_file,bac_hit_phages_dict)
		return result_dict_list

def predict_on_virhostmatcher(input_dir,phage_hit_hosts_now,outdir):
	script = os.path.join(root_path,'bin','PHIS_virhostmatcher')
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = script+" --input_dir "+input_dir+" --output "+outdir+" --out_dict "+phage_hit_hosts_now
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)
	os.system(command)

def predict_on_wish(input_dir,phage_hit_hosts_now,outdir):
	script = os.path.join(root_path,'bin',"PHIS_wish")
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = script+" --input_dir "+input_dir+" --output "+outdir+" --out_dict "+phage_hit_hosts_now
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)
	os.system(command) 

def predict_on_codon(file,phage_hit_hosts_now,outdir):
	script = os.path.join(root_path,'bin',"PHIS_codon")
	python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	if not os.path.exists(python_path):
		python_path = "python"
	command = script+" --input_file "+file+" --output "+outdir+" --out_dict "+phage_hit_hosts_now
	test_file = os.path.join(outdir,'script')
	with open(test_file,'w') as f:
		f.write(command)
	os.system(command)

def collect_features_virhostmatcher_wish_codon(phage_hit_hosts_now,outdir):
	virhostmatcher_dict_file = os.path.join(outdir,'virhost','virhostmatcher_result_dict.npy')
	wish_dict_file = os.path.join(outdir,'wish','wish_result_dict.npy')
	codon_dict_file = os.path.join(outdir,'codon','codon_result_dict.npy')
	#print(wish_dict_file)
	if os.path.exists(virhostmatcher_dict_file):
		virhostmatcher_dict = np.load(virhostmatcher_dict_file).item()
	if os.path.exists(wish_dict_file):
		wish_dict = np.load(wish_dict_file).item()
	if os.path.exists(codon_dict_file):
		codon_dict = np.load(codon_dict_file,encoding = 'latin1').item()
	genehomo_result_file = os.path.join(outdir,'genehomo_result.txt')
	f = open(genehomo_result_file,'w')
	f.write('phage_id\tphage_def\thost_id\thost_def\tvirhostmatcher_score\twish_score\tcodon_score\n')
	for query_id in phage_hit_hosts_now.keys():
		if query_id in strain_inf_dict.keys():
			query_def = strain_inf_dict[query_id]
		else:
			query_def = query_id
		for hit_id in phage_hit_hosts_now[query_id]:
			hit_id = hit_id.split('.')[0]
			#print(phage_id,hit_id)
			try:
				hit_def = phage_inf_dict[hit_id]
			except:
				hit_def = bac_inf_dict[hit_id]			
			try:
				virhostmatcher_score = virhostmatcher_dict[query_id][hit_id]
				# if float(virhostmatcher_score)<0.267:
				# 	virhostmatcher_score = 0
			except:
				virhostmatcher_score = 0
			try:
				wish_score = wish_dict[query_id][hit_id]
				# if float(wish_score)<-1.392:
				# 	wish_score = -2
			except:
				wish_score = -2
			try:
				codon_score = codon_dict[query_id][hit_id]
				# if float(codon_score)>0.084:
				# 	codon_score = 1
			except:
				codon_score = 1
			f.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(virhostmatcher_score)+'\t'+str(wish_score)+'\t'+str(codon_score)+'\n')
			f.flush()
	f.close()

def get_length(file,outdir,new_fna_file,max_seqs):
	#fasta or multi-fasta
	f_save = open(new_fna_file,'w')
	fasta_dir = os.path.join(outdir,'fasta_sequence')
	mkdir(fasta_dir)
	strains_length_dict = {}
	with open(file) as f:
		contents = f.read().strip()
	counter = 0
	if '\n>' in contents:
		for strain in contents.split('\n>'):
			strain_id = strain.strip().split('\n')[0].split()[0].split('.')[0].strip('>')		
			strain_dir = os.path.join(fasta_dir,strain_id)
			strain_sequence = ''
			for sequence in strain.strip().split('\n')[1:]:
				strain_sequence = sequence.strip()+strain_sequence
			strain_length = len(strain_sequence)
			strains_length_dict.update({strain_id:strain_length})
			if strain_length>=1000:
				counter = counter+1
				if counter<=int(max_seqs):
					mkdir(strain_dir)
					strain_file = os.path.join(fasta_dir,strain_id,strain_id+'.fasta')
					with open(strain_file,'w') as f:
						f.write('>'+strain.strip().strip('>')+'\n')
					f_save.write('>'+strain.strip().strip('>')+'\n')

	else:
		strain_id = contents.split('\n')[0].split()[0].split('.')[0].strip('>')
		strain_dir = os.path.join(fasta_dir,strain_id)
		mkdir(strain_dir)
		strain_file = os.path.join(fasta_dir,strain_id,strain_id+'.fasta')
		with open(strain_file,'w') as f:
			f.write('>'+contents.strip().strip('>')+'\n')
		f_save.write('>'+contents.strip().strip('>')+'\n')
		strain_sequence = ''
		for sequence in contents.strip().split('\n')[1:]:
			strain_sequence = sequence.strip()+strain_sequence
		strain_length = len(strain_sequence)
		strains_length_dict.update({strain_id:strain_length})
	outfile_dict = os.path.join(outdir,'phages_length_dict')
	np.save(outfile_dict,strains_length_dict)
	f_save.close()
	return strains_length_dict

def get_feature_value_dict(feature_file):
	feature_dict = {}
	with open(feature_file) as f:
		contents = f.readlines()
	for line in contents[1:]:
		line = line.strip().split('\t')
		query_id = line[0]
		query_def = line[1]
		hit_id = line[2]
		hit_def = line[3]
		value_list = line[4:]
		if query_id not in feature_dict.keys():
			feature_dict.update({query_id:{}})
		if hit_id not in feature_dict[query_id].keys():
			feature_dict[query_id].update({hit_id:line})
	return feature_dict

def collect_features(record_file,record_file_all,result_dict_list,module_dir,kind):
	if os.path.exists(record_file):
		hits_dict_now = np.load(record_file).item()
	else:
		hits_dict_now = {}
	if os.path.exists(record_file_all):
		hits_dict_now_all = np.load(record_file_all).item()
	else:
		hits_dict_now_all = {}
	genehomo_dir = os.path.join(module_dir,'genehomo')
	collect_features_virhostmatcher_wish_codon(hits_dict_now_all,genehomo_dir)
	
	genehomo_result_file = os.path.join(module_dir,'genehomo','genehomo_result.txt')
	pro_int_model_result_file = os.path.join(module_dir,'pro_pro_int','pro_int_result.txt')
	pro_int_feature_dict = get_feature_value_dict(pro_int_model_result_file)
	genehomo_feature_dict = get_feature_value_dict(genehomo_result_file)
	[crispr_feature_dict,prophage_feature_dict,blast_feature_dict] = result_dict_list
	
	overall_model_feature_file = os.path.join(module_dir,'overall_model_feature.txt')
	f_save = open(overall_model_feature_file,'w')
	if kind=='p':
		f_save.write('phage_id\tphage_def\thost_id\thost_def\thit_spacer_num\thit_spacer_identity\thit_spacer_coverage\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\tint_pro_num_homolog\tphage_pro_percent\tbac_pro_percent\tint_pro_num_domain\tdomain_int_phage_pro_percent\tdomain_int_bac_pro_percent\tblastp_score\tblastn_identity\tblastn_coverage\tvirhostmatcher_score\twish_score\tcondon_score\n')
		f_save.flush()
	else:
		f_save.write('bac_id\tbac_def\tphage_id\tphage_def\thit_spacer_num\thit_spacer_identity\thit_spacer_coverage\tprophage_homolog_percent\tprophage_alignment_identity\tprophage_alignment_coverage\tint_pro_num_homolog\tphage_pro_percent\tbac_pro_percent\tint_pro_num_domain\tdomain_int_phage_pro_percent\tdomain_int_bac_pro_percent\tblastp_score\tblastn_identity\tblastn_coverage\tvirhostmatcher_score\twish_score\tcondon_score\n')
		f_save.flush()	
	feature_array = []
	for query_id in hits_dict_now.keys():
		for hit_id in hits_dict_now[query_id]:
			each_array = []
			try:
				crispr_array = crispr_feature_dict[query_id][hit_id][4:]
				query_def = crispr_feature_dict[query_id][hit_id][1]
				hit_def = crispr_feature_dict[query_id][hit_id][3]
			except:
				crispr_array = [0,0,0]
			try:
				prophage_array = prophage_feature_dict[query_id][hit_id][4:]
				query_def = prophage_feature_dict[query_id][hit_id][1]
				hit_def = prophage_feature_dict[query_id][hit_id][3]
			except:
				prophage_array = [0,0,0]
			try:
				pro_int_array = pro_int_feature_dict[query_id][hit_id][4:]
				query_def = pro_int_feature_dict[query_id][hit_id][1]
				hit_def = pro_int_feature_dict[query_id][hit_id][3]
			except:
				pro_int_array = [0,0,0,0,0,0]
			try:
				blast_array = blast_feature_dict[query_id][hit_id][4:]
				query_def = blast_feature_dict[query_id][hit_id][1]
				hit_def = blast_feature_dict[query_id][hit_id][3]
			except:
				blast_array = [0,0,0]
			try:
				genehomo_array = genehomo_feature_dict[query_id][hit_id][4:]
				query_def = genehomo_feature_dict[query_id][hit_id][1]
				hit_def = genehomo_feature_dict[query_id][hit_id][3]
			except:
				genehomo_array = [0,0,0]
			each_array = [float(x) for x in crispr_array]+[float(x) for x in prophage_array]+[float(x) for x in pro_int_array]+[float(x) for x in blast_array]+[float(x) for x in genehomo_array]
			feature_array.append(each_array)
			f_save.write(query_id+'\t'+query_def+'\t'+hit_id+'\t'+hit_def+'\t'+'\t'.join(list(map(str,each_array)))+'\n')
			f_save.flush()
	f_save.close()
	return feature_array,overall_model_feature_file

def predict_model_signal(signal_result_file,special):
	signal_dict = {}
	if not os.path.exists(signal_result_file):
		return signal_dict
	with open(signal_result_file) as f:
		contents = f.readlines()
	for line in contents[1:]:
		line = line.strip().split('\t')
		query_id = line[0]
		hit_id = line[2]
		value_list = [float(item) for item in line[4:] if float(item)>0]
		if special==1:
			value_list = []
			if float(line[4])>=0.267:
				value_list.append(line[4])
			if float(line[5])>=-1.392:
				value_list.append(line[5])
			if float(line[6])<=0.084:
				value_list.append(line[6])
			#print(value_list)
		if len(value_list)>0:
			if query_id not in signal_dict.keys():
				signal_dict.update({query_id:{}})
			if hit_id not in signal_dict.keys():
				signal_dict[query_id].update({hit_id:'yes'})
	return signal_dict

def get_true_hits(module_dir,kind):
	save_file = os.path.join(module_dir,'crispr_prophage_blast_true_hit.txt')
	crispr_model_result_file = os.path.join(module_dir,'crispr','crispr_true_result.txt')
	prophage_model_result_file = os.path.join(module_dir,'prophage','prophage_true_result.txt')
	blast_model_result_file = os.path.join(module_dir,'blast','blast_true_result.txt')
	genehomo_result_file = os.path.join(module_dir,'genehomo','genehomo_result.txt')
	genehomo_signal_dict = predict_model_signal(genehomo_result_file,1)
	save_dict = {}
	if os.path.exists(crispr_model_result_file):
		with open(crispr_model_result_file) as f:
			crispr_contents = f.readlines()
		for line in crispr_contents[1:]:
			line = line.strip().split('\t')
			query_id = line[0]
			query_def = line[1]
			hit_id = line[2]
			hit_def = line[3]
			if query_id not in save_dict.keys():
				save_dict.update({query_id:{}})
			if hit_id not in save_dict[query_id].keys():
				save_dict[query_id].update({hit_id:[query_id,query_def,hit_id,hit_def,[],1]})
			save_dict[query_id][hit_id][4].append('crispr')
			try:
				genehomo_signal_dict[query_id][hit_id]
				if 'genehomo' not in save_dict[query_id][hit_id][4]:
					save_dict[query_id][hit_id][4].append('genehomo')
			except:
				pass
	
	if os.path.exists(prophage_model_result_file):
		with open(prophage_model_result_file) as f:
			prophage_contents = f.readlines()
		for line in prophage_contents[1:]:
			line = line.strip().split('\t')
			query_id = line[0]
			query_def = line[1]
			hit_id = line[2]
			hit_def = line[3]
			if query_id not in save_dict.keys():
				save_dict.update({query_id:{}})
			if hit_id not in save_dict[query_id].keys():
				save_dict[query_id].update({hit_id:[query_id,query_def,hit_id,hit_def,[],1]})
			save_dict[query_id][hit_id][4].append('prophage')
			try:
				genehomo_signal_dict[query_id][hit_id]
				if 'genehomo' not in save_dict[query_id][hit_id][4]:
					save_dict[query_id][hit_id][4].append('genehomo')
			except:
				pass

	if os.path.exists(blast_model_result_file):
		with open(blast_model_result_file) as f:
			blast_contents = f.readlines()
		for line in blast_contents[1:]:
			line = line.strip().split('\t')
			query_id = line[0]
			query_def = line[1]
			hit_id = line[2]
			hit_def = line[3]
			if query_id not in save_dict.keys():
				save_dict.update({query_id:{}})
			if hit_id not in save_dict[query_id].keys():
				save_dict[query_id].update({hit_id:[query_id,query_def,hit_id,hit_def,[],1]})
			save_dict[query_id][hit_id][4].append('blast')
			try:
				genehomo_signal_dict[query_id][hit_id]
				if 'genehomo' not in save_dict[query_id][hit_id][4]:
					save_dict[query_id][hit_id][4].append('genehomo')
			except:
				pass
	f_save = open(save_file,'w')
	if kind == 'p':
		f_save.write('phage_id\tphage_def\thost_id\thost_def\tpredict_signal\tscore\n')
	else:
		f_save.write('bac_id\tbac_def\tphage_id\tphage_def\tpredict_signal\tscore\n')	
	for query_id in save_dict.keys():
		for hit_id in save_dict[query_id].keys():
			f_save.write('\t'.join(save_dict[query_id][hit_id][0:4])+'\t'+','.join(save_dict[query_id][hit_id][4])+'\t'+str(save_dict[query_id][hit_id][-1])+'\n')
			f_save.flush()
	f_save.close()
	return save_dict

def predict_interaction_learning_model(predict_data,crispr_prophage_blast_dict,outdir,kind):
	bayes_bnb_model = joblib.load(os.path.join(root_path,'database/db/model/bayes_bnb.model'))
	bayes_gnb_model = joblib.load(os.path.join(root_path,'database/db/model/bayes_gnb.model'))
	randomforest_rfc_model = joblib.load(os.path.join(root_path,'database/db/model/randomforest_rfc.model'))
	randomforest_dtc_model = joblib.load(os.path.join(root_path,'database/db/model/randomforest_dtc.model'))
	svm_linear_model = joblib.load(os.path.join(root_path,'database/db/model/svm_linear.model'))
	svm_rbf_model = joblib.load(os.path.join(root_path,'database/db/model/svm_rbf.model'))
	logistic_model = joblib.load(os.path.join(root_path,'database/db/model/logistic.model'))
	bayes_bnb_result = bayes_bnb_model.predict(predict_data)
	bayes_bnb_result_probablity = bayes_bnb_model.predict_proba(predict_data)
	bayes_gnb_result = bayes_gnb_model.predict(predict_data)
	bayes_gnb_result_probablity = bayes_gnb_model.predict_proba(predict_data)
	randomforest_rfc_result = randomforest_rfc_model.predict(predict_data)
	randomforest_rfc_result_probablity = randomforest_rfc_model.predict_proba(predict_data)
	randomforest_dtc_result = randomforest_dtc_model.predict(predict_data)
	randomforest_dtc_result_probablity = randomforest_dtc_model.predict_proba(predict_data)
	svm_linear_result = svm_linear_model.predict(predict_data)
	svm_linear_result_probablity = svm_linear_model.predict_proba(predict_data)
	svm_rbf_result = svm_rbf_model.predict(predict_data)
	svm_rbf_result_probablity = svm_rbf_model.predict_proba(predict_data)
	logistic_result = logistic_model.predict(predict_data)
	logistic_result_probablity = logistic_model.predict_proba(predict_data)
	results = OrderedDict()
	predict_result_dict = {'bayes_bnb':[bayes_bnb_result,bayes_bnb_result_probablity],
	'bayes_gnb':[bayes_gnb_result,bayes_gnb_result_probablity],
	'randomforest_rfc':[randomforest_rfc_result,randomforest_rfc_result_probablity],
	'randomforest_dtc':[randomforest_dtc_result,randomforest_dtc_result_probablity],
	'svm_rbf':[svm_rbf_result,svm_rbf_result_probablity],
	'svm_linear':[svm_linear_result,svm_linear_result_probablity],
	'logistic':[logistic_result,logistic_result_probablity]}
	outfile = os.path.join(outdir,'predict_result_dict')
	np.save(outfile,predict_result_dict)

	outdir_results = os.path.join(outdir,'results')
	mkdir(outdir_results)
	overall_result_dict = OrderedDict()
	method_list = []
	for method,predict_values in predict_result_dict.items():
		method_list.append(method)
		out_file = os.path.join(outdir_results,method+"_result.txt")
		f_out = open(out_file,'w')
		if kind=='p':
			f_out.write('phage_id\tphage_def\thost_id\thost_def\tinteraction_type(1:interact,0:reverse)\tpredict_score(the possibilty of interaction)\n')
		else:
			f_out.write('bac_id\tbac_def\thit_phage_id\thit_phage_def\tinteraction_type(1:interact,0:reverse)\tpredict_score(the possibilty of interaction)\n')
		predict_label_list = predict_values[0]
		predict_score_list = predict_values[1]
		i = 0
		for query_id in crispr_prophage_blast_dict.keys():
			if query_id in strain_inf_dict.keys():
				query_def = strain_inf_dict[query_id]
			else:
				query_def = query_id
			hits_id = crispr_prophage_blast_dict[query_id]
			for hit_id in hits_id:
				if kind == 'p':
					hit_def = bac_inf_dict[hit_id]
				else:
					hit_def = phage_inf_dict[hit_id]
				f_out.write(query_id+'\t'+query_def+'\t'+hit_id+'\t'+hit_def+'\t'+str(predict_label_list[i])+'\t'+str(predict_score_list[i][1])+'\n')
				f_out.flush()
				if query_id not in overall_result_dict.keys():
					overall_result_dict.update({query_id:{}})
				if hit_id not in overall_result_dict[query_id].keys():
					overall_result_dict[query_id].update({hit_id:{}})
				if method not in overall_result_dict[query_id][hit_id].keys():
					overall_result_dict[query_id][hit_id].update({method:[]})
					overall_result_dict[query_id][hit_id][method]=[query_id,query_def,hit_id,hit_def,str(predict_label_list[i]),str(predict_score_list[i][1])]
				i = i+1
		f_out.close()
	
	overall_resufile = os.path.join(outdir_results,'model_result.txt')
	f_out = open(overall_resufile,'w')
	if kind=='p':
		f_out.write('phage_id\tphage_def\thost_id\thost_def\tpredict_signal\toverall_score\t'+'\t'.join(method_list)+'\n')
		f_out.flush()
	else:
		f_out.write('bac_id\tbac_def\tphage_id\tphage_def\tpredict_signal\toverall_score\t'+'\t'.join(method_list)+'\n')		
		f_out.flush()
	new_return_dict = {}
	crispr_result_file = os.path.join(module_dir,'crispr','crispr_model_result.txt')
	prophage_result_file = os.path.join(module_dir,'prophage','prophage_model_result.txt')
	pro_int_result_file = os.path.join(module_dir,'pro_int','pro_int_result.txt')
	blast_result_file = os.path.join(module_dir,'blast','blast_model_result.txt')
	genehomo_result_file = os.path.join(module_dir,'genehomo','genehomo_result.txt')
	crispr_signal_dict = predict_model_signal(crispr_result_file,0)
	prophage_signal_dict = predict_model_signal(prophage_result_file,0)

	pro_int_signal_dict = predict_model_signal(pro_int_result_file,0)
	blast_signal_dict = predict_model_signal(blast_result_file,0)
	genehomo_signal_dict = predict_model_signal(genehomo_result_file,1)
	for query_id in overall_result_dict.keys():
		for hit_id in overall_result_dict[query_id].keys():
			signal = []
			try:
				crispr_signal_dict[query_id][hit_id]
				signal.append('crispr')
			except:
				pass
			try:
				prophage_signal_dict[query_id][hit_id]
				signal.append('prophage')
			except:
				pass
			try:
				pro_int_signal_dict[query_id][hit_id]
				signal.append('protein_protein_interaction')
			except:
				pass
			try:
				blast_signal_dict[query_id][hit_id]
				signal.append('blast')
			except:
				pass
			try:
				genehomo_signal_dict[query_id][hit_id]
				signal.append('sequence_composition')
			except:
				pass
			class_list = []
			score_list = []
			if len(signal)==0:
				signal = ['none']
				#continue
			counter = 0
			for method in overall_result_dict[query_id][hit_id].keys():
				class_list.append(int(overall_result_dict[query_id][hit_id][method][-2].split('.')[0]))
				score_list.append(overall_result_dict[query_id][hit_id][method][-1])
				if float(overall_result_dict[query_id][hit_id][method][-1])>=0.8:
					counter = counter+1
			score = np.mean(list(map(float,score_list)))
			if counter>=4:#at least 4 methods support interaction
				if query_id not in new_return_dict.keys():
					new_return_dict.update({query_id:[]})
				new_return_dict[query_id].append(overall_result_dict[query_id][hit_id][method][0:-2]+[','.join(signal)]+[str(score)]+score_list)
	for query_id in new_return_dict.keys():
		hits = new_return_dict[query_id]
		hits.sort(key=lambda pair:(len(pair[4].split(',')),float(pair[5])),reverse=True)
		for hit in hits:
			f_out.write('\t'.join(list(map(str,hit)))+'\n')
			f_out.flush()
	f_out.close()
	save_model_dict_file = os.path.join(outdir_results,'model_result_dict')
	np.save(save_model_dict_file,new_return_dict)
	return new_return_dict

def get_final_result(true_result_dict,model_result_dict_file,outdir,kind):
	if os.path.exists(model_result_dict_file):
		model_result_dict = np.load(model_result_dict_file).item()
	else:
		model_result_dict = {}
	outdir = os.path.join(outdir,'results')
	mkdir(outdir)
	save_file = os.path.join(outdir,'overall_result.txt')
	f_save = open(save_file,'w')
	method_list=  ['bayes_bnb',
	'bayes_gnb',
	'randomforest_rfc',
	'randomforest_dtc',
	'svm_rbf',
	'svm_linear',
	'logistic']
	if kind=='p':
		f_save.write('phage_id\tphage_def\thost_id\thost_def\tpredict_signal\toverall_score\t'+'\t'.join(method_list)+'\n')
		f_save.flush()
	else:
		f_save.write('bac_id\tbac_def\tphage_id\tphage_def\tpredict_signal\toverall_score\t'+'\t'.join(method_list)+'\n')		
		f_save.flush()
	#print(model_result_dict)
	for query_id in true_result_dict.keys():
		for hit_id in true_result_dict[query_id].keys():
			f_save.write('\t'.join(true_result_dict[query_id][hit_id][0:4])+'\t'+','.join(true_result_dict[query_id][hit_id][4])+'\t1\t'+'\t'.join(['NA']*7)+'\n')
			f_save.flush()
	for query_id in model_result_dict.keys():
		for hit_record in model_result_dict[query_id]:
			f_save.write('\t'.join(hit_record)+'\n')
			f_save.flush()
	f_save.close()

if __name__=='__main__':
	multiprocessing.freeze_support()
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help='Path of the input file.\n')
	parser.add_argument('--output', help='Path of the output file.\n')
	parser.add_argument('--min_mis_crispr', help='Minimal mismatch of a Blastn hit on hit spacers(default: 2)\n')
	parser.add_argument('--min_cov_crispr', help='Minimal % coverage of a Blastn hit on hit spacers(default: 90)\n')
	parser.add_argument('--min_per_prophage', help='Minimal % percentage of hit proteins on hit prophage region(default:30)\n')	
	parser.add_argument('--min_id_prophage', help='Minimal % identity of hit region on hit prophage region by making blastn(default:70)\n')	
	parser.add_argument('--min_cov_prophage', help='Minimal % coverage of hit region on hit prophage region by making blastn(default:30)')	
	parser.add_argument('--min_PPI', help='Minimal PPI number of a pair of phage-host pair(default:1)\n')	
	parser.add_argument('--min_DDI', help='Minimal DDI number of a pair of phage-host pair(default:5)\n')	
	parser.add_argument('--min_per_blast', help='Minimal % percentage of hit proteins on query phage and host by making blastp(default:10)\n')	
	parser.add_argument('--min_id_blast', help='Minimal % identity of hit region on query phage and host by making blastn(default:70)\n')		
	parser.add_argument('--min_cov_blast', help='Minimal % coverage of hit region on query phage and host by making blastn(default:10)\n')		
	args = parser.parse_args()
	global type,strain_inf_dict,pro_num_dict,bac_inf_dict,phage_inf_dict,root_path,min_mis_crispr,min_cov_crispr,min_per_prophage,min_id_prophage,min_cov_prophage,min_PPI,min_DDI,min_per_blast,min_id_blast,min_cov_blast
	if args.input:
		input_file = args.input
	if args.output:
		outdir = args.output
	if args.min_mis_crispr:
		min_mis_crispr = args.min_mis_crispr
	else:
		min_mis_crispr = 2
	if args.min_cov_crispr:
		min_cov_crispr = args.min_cov_crispr
	else:
		min_cov_crispr = 70
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
	if args.min_PPI:
		min_PPI = args.min_PPI
	else:
		min_PPI = 1
	if args.min_DDI:
		min_DDI = args.min_DDI
	else:
		min_DDI = 5
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
	kind= 'p'
	root_path = get_root_path()
	bac_inf_dict_file = os.path.join(root_path,'database/db/profile','bac_inf_dict.npy')
	bac_inf_dict = np.load(bac_inf_dict_file).item()
	#create dirs
	strain_info_dir = os.path.join(outdir,'strain_info')
	module_dir = os.path.join(outdir,'module_results')
	mkdir(strain_info_dir)
	mkdir(module_dir)

	file_name = os.path.basename(input_file).split()[0].split('.')[0]
	faa_prefix = os.path.join(outdir,file_name)
	pro_file = faa_prefix+'.faa'
	strain_inf_dict,type = get_inf(input_file,strain_info_dir)
	if type == 'fasta':
		pred_orf(input_file,faa_prefix)
	else:
		getFaaFromGB(input_file,outdir)
		getFnaFromGB(input_file,outdir,file_name)
		input_file = os.path.join(outdir,file_name+'.fna')
	
	#get phage list
	pro_num_dict = get_pro_num(pro_file,strain_info_dir)
	new_fna_file = os.path.join(strain_info_dir,'filtered_phage.fna')
	phages_length_dict = get_length(input_file,strain_info_dir,new_fna_file,1000)
	if type=='fasta':
		input_file = new_fna_file
	fasta_dir = os.path.join(strain_info_dir,'fasta_sequence')
	#step1:run crispr,prophage,blast
	print('step1:start to run crispr,prophage,blast')
	run_crispr_prophage_blast(input_file,module_dir,kind)

	#step2:collect predict results through crispr,prophage,blastp
	print('step2:start to analyse crispr,prophage,blast,get the reliable hits and candidate hits!')
	result_dict_list_now,feature_dict_now_all = collect_result_1(input_file,pro_file,module_dir,kind)

	#step3: run pro_int,virhostmatcher,wish,codon_usage
	print('step3:start to run virhostmatcher,wish,codon_usage and pro_int according to step2!')
	genehomo_dir = os.path.join(module_dir,'genehomo')
	mkdir(genehomo_dir)
	if kind=='p':
		record_file = os.path.join(module_dir,'phage_hit_hosts_crispr_prophage_blast_dict.npy')
		record_file_all = os.path.join(module_dir,'phage_hit_hosts_crispr_prophage_blast_all_dict.npy')
	else:
		record_file = os.path.join(module_dir,'bac_hit_phages_crispr_prophage_blast_dict.npy')
	crispr_prophage_blast_dict = np.load(record_file).item()
	run_genehomo_pro_int(input_file,pro_file,record_file_all,module_dir,fasta_dir,kind)
	
	#step4: calculate 18 features
	print('step4:start to calculate 18 features for candidate hits')
	predict_data,overall_feature_file = collect_features(record_file,record_file_all,result_dict_list_now,module_dir,kind)

	#step5:predict model results
	print('step5:start to predict interaction for candidate hits')
	if len(predict_data)>0:
		predict_interaction_learning_model(predict_data,crispr_prophage_blast_dict,outdir,kind)

	#step6:get true hits dict
	print('step6:start to count results for reliable hits')
	true_result_dict = get_true_hits(module_dir,kind)

	#step7:merge result
	print('step7:start to count results for all hits')
	model_result_dict_file = os.path.join(outdir,'results','model_result_dict.npy')
	get_final_result(true_result_dict,model_result_dict_file,outdir,kind)

	print('PHIS finished prediction in %s'%outdir)