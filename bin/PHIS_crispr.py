import os,re,argparse
import numpy.random.common
import numpy.random.bounded_integers
import numpy.random.entropy
import numpy as np
from Bio import Entrez, SeqIO
import threading

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

def phage_blastn_spacer(phage_file,outfile,format):
	spacer_database = os.path.join(root_path,'database/db/database/merge_spacer','merge_spacer_db')
	blastn_path = os.path.join(root_path,'software/blast+','blastn')
	command = blastn_path+' -query %s -db %s -evalue 0.01 -outfmt %s -out %s -word_size 7 -dust no -soft_masking FALSE -gapopen 10 -penalty -1 -gapextend 2 -ungapped -num_threads 20 -task blastn-short' % (phage_file,spacer_database,str(format),outfile)	
	os.system(command)

def get_num(string):
	strlist = re.findall("\d+\.?\d*",str(string))
	return strlist

def parse_spacer_blastn_phage_0(file,resufile):
	with open(file) as f:
		contents = f.read().strip()
	spacers_hits = contents.split('Query=')[1:]
	hit_dict = {}
	f_result = open(resufile,'w')
	f_result.write('bac_id\tbac_def\thit_phage_id\thit_phage_def\tspacer\tspacer_length\talignment_length\tidentity\tcoverage\tmismatch\tevalue\tbit-score\n')
	print(strain_inf_dict)
	for each_spacer_hits in spacers_hits:
		if 'No hits found' in each_spacer_hits:
			continue
		spacer_name = each_spacer_hits.split('\n')[0].strip()
		bac_id = spacer_name.split('|')[-1].split('.')[0]
		try:
			bac_def = strain_inf_dict[bac_id]
		except:
			bac_def = bac_id
		spacer_length = float(spacer_name.split('|')[2])
		spacer_hits = each_spacer_hits.split('> ')[1:]
		for spacer_hit in spacer_hits:
			hit_phage = spacer_hit.split('\n')[0].strip()
			hit_phage_id = hit_phage.split()[0].split('.')[0]
			hit_phage_def = ' '.join(hit_phage.split()[1:])
			bit_score = float(get_num(spacer_hit.split('\n\n')[1].split('\n')[0])[0])
			evalue = float(spacer_hit.split('\n\n')[1].split('\n')[0].split('=')[-1].strip())
			identity = float(get_num(spacer_hit.split('\n\n')[1].split('\n')[1])[2])
			hit_length1 = int(get_num(spacer_hit.split('\n\n')[1].split('\n')[1])[0])
			hit_length = int(get_num(spacer_hit.split('\n\n')[1].split('\n')[1])[1])
			mismatch = hit_length-hit_length1
			coverage = hit_length/spacer_length
			hit_sequence_details = spacer_hit.split('\n\n')[2]
			positions = get_num(hit_sequence_details)
			hit_details = hit_sequence_details.split('\n')
			query = hit_details[0].split()[2].strip()
			subj = hit_details[-1].split()[2].strip()
			hit_info= query+'\t'+'Query'+'\t'+str(positions[0])+':'+str(positions[1])+'\n'+hit_details[1][len(hit_details[1])-len(query):]+'\n'+subj+'\t'+'Sbjct'+'\t'+str(positions[2])+':'+str(positions[3])             
			if (int(mismatch) <= 2) and (coverage>=0.7):
				if bac_id not in hit_dict.keys():
					hit_dict.update({bac_id:{}})
				if hit_phage_id not in hit_dict[bac_id].keys():
					hit_dict[bac_id].update({hit_phage_id:[]})
				hit_dict[bac_id][hit_phage_id].append([hit_phage_id,spacer_name,hit_length,bit_score,identity,coverage,mismatch,evalue,hit_info]) 
				resuline = list(map(str,[bac_id,bac_def,hit_phage_id,hit_phage_def,spacer_name,spacer_length,hit_length,identity,coverage,mismatch,evalue,bit_score]))				
				f_result.write('\t'.join(resuline)+'\n')
				f_result.flush()
	hit_dict_file = os.path.join(outdir,'hit_phages_dict')
	np.save(hit_dict_file,hit_dict)
	f_result.close()

def parse_spacer_blastn_phage_6(file,resufile):
	with open(file) as f:
		contents = f.readlines()
	hit_dict = {}
	f_result = open(resufile,'w')
	f_result.write('bac_id\tbac_def\thit_phage_id\thit_phage_def\tspacer\tspacer_length\talignment_length\tidentity\tcoverage\tmismatch\tbit-score\n')
	f_result.flush()
	for line in contents:
		line = line.strip().split('\t')
		phage_id = line[1].split('.')[0]
		phage_def = phage_inf_dict[phage_id]
		spacer_name = line[0]
		spacer_length = float(spacer_name.split('|')[2])
		identity = float(line[2])
		#hit_length = int(line[3])
		mismatch = int(line[4])
		bit_score = float(line[-1])
		query_start = int(line[-6])
		query_end = int(line[-5])
		hit_length = abs(query_start-query_end)+1
		coverage = float(hit_length)/spacer_length
		#print(line,subject_end,subject_start,spacer_length,coverage)
		bac_id = spacer_name.split('|')[-1].split('.')[0]
		try:
			bac_def = strain_inf_dict[bac_id]
		except:
			bac_def = bac_id
		evalue = float(line[-2])
		if (int(mismatch) <= 2) and (coverage>=0.7):
			if bac_id not in hit_dict.keys():
				hit_dict.update({bac_id:{}})
			if phage_id not in hit_dict[bac_id].keys():
				hit_dict[bac_id].update({hit_phage_id:[]})
			hit_dict[bac_id][phage_id].append([phage_id,spacer_name,hit_length,bit_score,identity,coverage,mismatch,evalue]) 
			resuline = list(map(str,[bac_id,bac_def,phage_id,phage_def,spacer_name,spacer_length,hit_length,identity,coverage,mismatch,evalue,bit_score]))			
			f_result.write('\t'.join(resuline)+'\n')
			f_result.flush()
	f_result.close()
	hit_dict_file = os.path.join(outdir,'hit_phages_dict')
	np.save(hit_dict_file,hit_dict)

def parse_phage_blastn_spacer_6(file,resufile):
	with open(file) as f:
		contents = f.readlines()
	hit_dict = {}
	f_result = open(resufile,'w')
	f_result.write('phage_id\tphage_def\thit_bac_id\thit_bac_def\tspacer\tspacer_length\talignment_length\tidentity\tcoverage\tmismatch\tbit-score\n')
	for line in contents:
		line = line.strip().split('\t')
		phage_id = line[0].split('.')[0]
		try:
			phage_def = strain_inf_dict[phage_id]
		except:
			phage_def = phage_id
		spacer_name = line[1]
		spacer_length = float(spacer_name.split('|')[2])
		identity = float(line[2])
		#hit_length = int(line[3])
		mismatch = int(line[4])
		bit_score = float(line[-1])
		subject_start = int(line[-4])
		subject_end = int(line[-3])
		hit_length = abs(subject_start-subject_end)+1
		coverage = float(hit_length)/spacer_length
		#print(line,subject_end,subject_start,spacer_length,coverage)
		hit_bac_id = spacer_name.split('|')[-1].split('.')[0]
		hit_bac_def = bac_inf_dict[hit_bac_id]
		evalue = float(line[-2])
		if (int(mismatch) <= 2) and (coverage>=0.7):
			if phage_id not in hit_dict.keys():
				hit_dict.update({phage_id:{}})
			if hit_bac_id not in hit_dict[phage_id].keys():
				hit_dict[phage_id].update({hit_bac_id:[]})
			hit_dict[phage_id][hit_bac_id].append([phage_id,hit_bac_id,spacer_name,bit_score,identity,coverage,mismatch,evalue]) 
			resuline = list(map(str,[phage_id,phage_def,hit_bac_id,hit_bac_def,spacer_name,spacer_length,hit_length,identity,coverage,mismatch,evalue,bit_score]))
			f_result.write('\t'.join(resuline)+'\n')
			f_result.flush()
	f_result.close()
	hit_dict_file = os.path.join(outdir,'hit_bac_spacers_dict')
	np.save(hit_dict_file,hit_dict)

def parse_phage_blastn_spacer_0(file,resufile):
	with open(file) as f:
		contents = f.read().strip()
	phage_hit_spacers = contents.split('Query=')[1:]
	hit_dict = {}
	f_result = open(resufile,'w')
	f_result.write('phage_id\tphage_def\thit_bac_id\thit_bac_def\tspacer\tspacer_length\talignment_length\tidentity\tcoverage\tmismatch\tevalue\tbit-score\n')
	for each_hits in phage_hit_spacers:
		if 'No hits found' in each_hits:
			continue
		phage_name = each_hits.split('\n')[0].strip()
		phage_id = phage_name.split()[0].split('.')[0]
		if phage_id not in hit_dict.keys():
			hit_dict.update({phage_id:{}})
		hit_spacers = each_hits.split('> ')[1:]
		for spacer_hit in hit_spacers:
			hit_spacer_name = spacer_hit.split('\n')[0].strip()
			spacer_length = float(hit_spacer_name.split('|')[2])
			hit_bac_id = hit_spacer_name.split('|')[3].split('.')[0]
			if hit_bac_id in bac_inf_dict.keys():
				hit_bac_def = bac_inf_dict[hit_bac_id]
			else:
				hit_bac_def = hit_spacer_name.split('|')[-1]
			bit_score = float(get_num(spacer_hit.split('\n\n')[1].split('\n')[0])[0])
			evalue = float(spacer_hit.split('\n\n')[1].split('\n')[0].split('=')[-1].strip())
			identity = float(get_num(spacer_hit.split('\n\n')[1].split('\n')[1])[2])
			hit_length1 = int(get_num(spacer_hit.split('\n\n')[1].split('\n')[1])[0])
			hit_length = int(get_num(spacer_hit.split('\n\n')[1].split('\n')[1])[1])
			mismatch = hit_length-hit_length1
			coverage = hit_length/spacer_length
			hit_sequence_details = spacer_hit.split('\n\n')[2]
			positions = get_num(hit_sequence_details)
			hit_details = hit_sequence_details.split('\n')
			query = hit_details[0].split()[2].strip()
			subj = hit_details[-1].split()[2].strip()
			hit_info = query+'\t'+'Query'+'\t'+str(positions[0])+':'+str(positions[1])+'\n'+hit_details[1][len(hit_details[1])-len(query):]+'\n'+subj+'\t'+'Sbjct'+'\t'+str(positions[2])+':'+str(positions[3])             
			if (mismatch <= float(min_mis_crispr)) and (coverage>=float(min_cov_crispr)/100):
				if hit_bac_id not in hit_dict[phage_id].keys():
					hit_dict[phage_id].update({hit_bac_id:[]})
				hit_dict[phage_id][hit_bac_id].append([phage_id,hit_bac_id,hit_spacer_name,bit_score,identity,coverage,mismatch,evalue,hit_info]) 
				resuline = list(map(str,[phage_id,phage_name,hit_bac_id,hit_bac_def,hit_spacer_name,spacer_length,hit_length,identity,coverage,mismatch,evalue,bit_score]))
				f_result.write('\t'.join(resuline)+'\n')
				f_result.flush()
	hit_dict_file = os.path.join(outdir,'hit_bac_spacers_dict')
	np.save(hit_dict_file,hit_dict)
	f_result.close()

def phage_to_hosts(phage_file,outdir):
	outfile = os.path.join(outdir,'phage_blastn_spacer')
	phage_blastn_spacer(phage_file,outfile,format)
	resufile = os.path.join(outdir,'hit_bac_spacers_result.txt')
	if os.path.exists(outfile):
		if os.path.getsize(outfile)>0:
			if format==6:
				parse_phage_blastn_spacer_6(outfile,resufile)
			else:
				parse_phage_blastn_spacer_0(outfile,resufile)
			crispr_result_dict_file = os.path.join(outdir,'hit_bac_spacers_dict.npy')
			get_crispr_result(crispr_result_dict_file,outdir,'p')

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

def get_crispr_result(crispr_result_file,outdir,kind):
	if os.path.exists(crispr_result_file):
		crispr_result_dict = np.load(crispr_result_file).item()
	else:
		crispr_result_dict = {}
	strong_result_file = os.path.join(outdir,'crispr_true_result.txt')
	weak_result_file = os.path.join(outdir,'crispr_model_result.txt')
	all_result_file = os.path.join(outdir,'crispr_result.txt')
	f_crispr_host = open(strong_result_file,'w')
	f_crispr_model = open(weak_result_file,'w')
	f_crispr_all = open(all_result_file,'w')
	if kind=='p':
		f_crispr_host.write('phage_id\tphage_def\thost_id\thost_def\thit_spacer_num\thit_spacer_identity\thit_spacer_coverage\n')	
		f_crispr_model.write('phage_id\tphage_def\thost_id\thost_def\thit_spacer_num\thit_spacer_identity\thit_spacer_coverage\n')
		f_crispr_all.write('phage_id\tphage_def\thost_id\thost_def\thit_spacer_num\thit_spacer_identity\thit_spacer_coverage\tmode\n')
	else:
		f_crispr_host.write('bac_id\tbac_def\tphage_id\tphage_def\thit_spacer_num\thit_spacer_identity\thit_spacer_coverage\n')	
		f_crispr_model.write('bac_id\tbac_def\tphage_id\tphage_def\thit_spacer_num\thit_spacer_identity\thit_spacer_coverage\n')
		f_crispr_all.write('bac_id\tbac_def\tphage_id\tphage_def\thit_spacer_num\thit_spacer_identity\thit_spacer_coverage\tmode\n')	
	for query_id in crispr_result_dict.keys():
		if query_id in strain_inf_dict.keys():
			query_def = strain_inf_dict[query_id]
		else:
			query_def = query_id
		for hit_id in crispr_result_dict[query_id]:
			query_id = query_id.split('.')[0]
			hit_id = hit_id.split('.')[0]
			hit_spacer_num = 0
			hit_identity = 0
			hit_coverage = 0
			spacer_names = list(set([x[2] for x in crispr_result_dict[query_id][hit_id]]))
			db_flag = 0
			try:
				hit_def = phage_inf_dict[hit_id]
			except:
				if hit_id in bac_inf_dict.keys():
					hit_def = bac_inf_dict[hit_id]
				else:
					db_flag = 1
					hit_def = spacer_names[0].split('|')[-1]
			hit_spacer_num = len(spacer_names)
			crispr_result_dict[query_id][hit_id].sort(key=lambda pair:(float(pair[3]),float(pair[4]),float(pair[5])),reverse=True)
			#get the best identity and coverage with highest bit-score
			hit_identity = crispr_result_dict[query_id][hit_id][0][4]
			hit_coverage = crispr_result_dict[query_id][hit_id][0][5]
			flag = 0
			for hit_record in crispr_result_dict[query_id][hit_id]:
				hit_spacer_name = hit_record[2]
				spacer_length = hit_spacer_name.split('|')[2]
				identity = hit_record[4]
				coverage = hit_record[5]
				mismatch = hit_record[6]
				bit_score = hit_record[3]
				evalue = hit_record[7]
				if float(coverage)>=0.95:
					flag=1
			if flag==0:
				if db_flag==0:#coverage<0.95 and is in the complete bacteria list
					f_crispr_model.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(hit_spacer_num)+'\t'+str(hit_identity)+'\t'+str(hit_coverage)+'\n')
					f_crispr_model.flush()
				mode = "weak"
			else:
				f_crispr_host.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(hit_spacer_num)+'\t'+str(hit_identity)+'\t'+str(hit_coverage)+'\n')
				f_crispr_host.flush()
				mode = "strong"
			f_crispr_all.write('\t'.join([query_id,query_def,hit_id,hit_def])+'\t'+str(hit_spacer_num)+'\t'+str(hit_identity)+'\t'+str(hit_coverage)+'\t'+mode+'\n')
			f_crispr_all.flush()
	f_crispr_model.close()
	f_crispr_host.close()

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help='''Path of the input file in genbank or fasta format for bacteria or phages
		please ensure your file name consist of word,number,_,if not,the program will probably run error!''')
	parser.add_argument('--output', help='Path of the output directory.\n')
	parser.add_argument('--type', help='optional,value:genbank or fasta.\n')
	parser.add_argument('--format', help='optional,output format of blastn result,value:0 or 6,,default:0.\n')
	parser.add_argument('--min_mis_crispr', help='Minimal mismatch of a Blastn hit on hit spacers(default: 2)\n')
	parser.add_argument('--min_cov_crispr', help='Minimal % coverage of a Blastn hit on hit spacers(default: 70)\n')
	args = parser.parse_args()
	global type,format,strain_inf_dict,min_mis_crispr,min_cov_crispr,root_path,bac_inf_dict
	format = 0
	if args.input:
		input_file = args.input
	if args.output:
		outdir = args.output
	file_name = os.path.basename(input_file).split()[0].split('.')[0]
	strain_inf_dict,type = get_inf(input_file,outdir)
	if args.format:
		format = args.format
	if type == 'genbank':
		getFnaFromGB(input_file,outdir,file_name)
		input_file = os.path.join(outdir,file_name+'.fna')
	if args.min_mis_crispr:
		min_mis_crispr = args.min_mis_crispr
	else:
		min_mis_crispr = 2
	if args.min_cov_crispr:
		min_cov_crispr = args.min_cov_crispr
	else:
		min_cov_crispr = 70
	root_path = get_root_path()
	bac_inf_dict_file = os.path.join(root_path,'database/db/profile','bac_inf_dict.npy')
	bac_inf_dict = np.load(bac_inf_dict_file).item()
	phage_to_hosts(input_file,outdir)
	print('finished predicting hosts in %s!'%outdir)