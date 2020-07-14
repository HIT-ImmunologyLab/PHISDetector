# -*- coding: utf-8 -*-
import os,argparse
import sys
from Bio import Entrez, SeqIO
import numpy.random.common
import numpy.random.bounded_integers
import numpy.random.entropy
import numpy as np
from collections import OrderedDict
def inputinstructions():
	return """PHIS 1.1.0 arguments:

Usage: PHIS [options]
Options (x is an integer number)

--input <file name>        : Query phage file path: FASTA or Multi-Fasta or GenBank file
--out <folder name>        : Output folder in which results will be stored
--signal <module name>     : Model name(default: all):
                           1.all 
                              run 5 modules:CRISPR,Prophage,Protein_protein interaction,Genetic homology and Oligonucleotide profile/sequence composition to predict the hosts of query phage
                           2.crispr
                              run CRISPR module to predict the hosts of query phage by making blastn with spacer database
                           3.prophage
                              run Prophage module to predict the hosts of query phage by making blastp and blastn with prophage database
                           4.protein_protein_interaction
                              run Protein_protein interaction module to predict the hosts of query phage by finding PPIs(protein_protein interaction) and DDIs(domain_domain interaction)
                           5.blast
                              run Genetic homology module to predict the hosts of query phage by making blastp and blastn with bacteria database
--min_mis_crispr <x>       : Minimal mismatch of a Blastn hit on hit spacers(default: 2)
--min_cov_crispr <x>       : Minimal % coverage of a Blastn hit on hit spacers(default: 70)
--min_per_prophage <x>     : Minimal % percentage of hit proteins on hit prophage region(default:30)
--min_id_prophage <x>      : Minimal % identity of hit region on hit prophage region by making blastn(default:70)
--min_cov_prophage <x>     : Minimal % coverage of hit region on hit prophage region by making blastn(default:30)
--min_PPI <x>              : Minimal PPI number of a pair of phage-host pair(default:1)
--min_DDI <x>              : Minimal DDI number of a pair of phage-host pair(default:5)
--min_per_blast <x>        : Minimal % percentage of hit proteins on query phage and host by making blastp(default:10)
--min_id_blast <x>         : Minimal % identity of hit region on query phage and host by making blastn(default:70)
--min_cov_blast <x>        : Minimal % coverage of hit region on query phage and host by making blastn(default:10)
"""

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

def predict_on_all(file,outdir,parameter_dict):
	script = os.path.join(root_path,'bin',"PHIS_overall_phage_to_host_two_criteria")
	# python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	# if not os.path.exists(python_path):
	# 	python_path = "python"
	command = script+" --input "+file+" --output "+outdir
	for par_name in parameter_dict.keys():
		if par_name!='model':
			command = command+" --"+par_name+" "+str(parameter_dict[par_name])
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)
	os.system(command)


def predict_on_crispr(file,outdir,min_mis_crispr=2,min_cov_crispr=70):
	script = os.path.join(root_path,'bin','PHIS_crispr')
	# python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	# if not os.path.exists(python_path):
	# 	python_path = "python"
	command = script+" --input "+file+" --output "+outdir+" --min_mis_crispr "+str(min_mis_crispr)+" --min_cov_crispr "+str(min_cov_crispr)
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)
	os.system(command)

def predict_on_prophage(file,outdir,min_per_prophage=30,min_id_prophage=70,min_cov_prophage=30):
	script = os.path.join(root_path,'bin','PHIS_prophage')
	#python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	# if not os.path.exists(python_path):
	# 	python_path = "python"
	command = script+" --input_file "+file+" --output "+outdir+" --min_per_prophage "+str(min_per_prophage)+" --min_id_prophage "+str(min_id_prophage)+" --min_cov_prophage "+str(min_cov_prophage)
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)
	os.system(command)

def predict_on_blast(file,outdir,min_per_blast=10,min_id_blast=70,min_cov_blast=10):
	script = os.path.join(root_path,'bin','PHIS_blast')
	# python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	# if not os.path.exists(python_path):
	# 	python_path = "python"
	command = script+" --input_file "+file+" --output "+outdir+" --min_per_blast "+str(min_per_blast)+" --min_id_blast "+str(min_id_blast)+" --min_cov_blast "+str(min_cov_blast)
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)
	os.system(command)

def predict_on_pro_int(file,outdir,min_PPI=1,min_DDI=5):
	script = os.path.join(root_path,'bin','PHIS_protein_protein_int')
	# python_path = os.path.join(root_path,'software/python/miniconda3/bin','python')
	# if not os.path.exists(python_path):
	# 	python_path = "python"
	command = script+" --input_file "+file+" --output "+outdir+" --min_PPI "+str(min_PPI)+" --min_DDI "+str(min_DDI)
	with open(os.path.join(outdir,'script'),'w') as f:
		f.write(command)
	os.system(command)

def set_default_parameter():
	min_mis_crispr = 2
	min_cov_crispr = 70
	min_per_prophage = 30
	min_id_prophage = 70
	min_cov_prophage = 30
	min_PPI = 1
	min_DDI = 5
	min_per_blast = 10
	min_id_blast = 70
	min_cov_blast = 10
	model = 'all'
	default_parameter_dict = {'min_mis_crispr':min_mis_crispr,
	'min_cov_crispr':min_cov_crispr,
	'min_per_prophage':min_per_prophage,
	'min_id_prophage':min_id_prophage,
	'min_cov_prophage':min_cov_prophage,
	'min_PPI':min_PPI,
	'min_DDI':min_DDI,
	'min_per_blast':min_per_blast,
	'min_id_blast':min_id_blast,
	'min_cov_blast':min_cov_blast,
	'model':model}
	return default_parameter_dict

def check_file(file):
	with open(file) as f:
		contents = f.read()
	if (contents[0]!='>'):
		try:
			SeqIO.parse(file,"genbank")
		except:
			print('phage file format error!please input genbank or fasta file')
			sys.exit(-1)

def check_outdir(outdir):
	if not os.path.exists(outdir):
		print('the output directory does not exists!')
		sys.exit(-1)

def check_parameter(args):
	default_parameter_dict = set_default_parameter()
	variable_dict = vars(args)
	for name in variable_dict.keys():
		if name=='input':
			continue
		if name=='output':
			continue
		if name=='model':
			if args.model:
				input_model = args.model
				if input_model not in ['all','crispr','prophage','protein_protein_interaction','blast','opc']:
					print('there no exists %s'%input_model)
					sys.exit(-1)
		else:
			if variable_dict[name]:
				try:
					cutoff = int(variable_dict[name])
				except:
					print('%s needs integer!%s is not integer!'%(name,variable_dict[name]))
					sys.exit(-1)
		if variable_dict[name]:
			default_parameter_dict[name] = variable_dict[name]
	return default_parameter_dict

def run(input_file,outdir,parameter_dict):
	if parameter_dict['model'] == 'all':
		predict_on_all(input_file,outdir,parameter_dict)	
	if parameter_dict['model'] == 'crispr':
		predict_on_crispr(input_file,outdir,parameter_dict['min_mis_crispr'],parameter_dict['min_cov_crispr'])
	if parameter_dict['model'] == 'prophage':
		predict_on_prophage(input_file,outdir,parameter_dict['min_per_prophage'],parameter_dict['min_id_prophage'],parameter_dict['min_cov_prophage'])
	if parameter_dict['model'] == 'protein_protein_interaction':
		predict_on_pro_int(input_file,outdir,parameter_dict['min_PPI'],parameter_dict['min_DDI'])
	if parameter_dict['model'] == 'blast':
		predict_on_blast(input_file,outdir,parameter_dict['min_per_blast'],parameter_dict['min_id_blast'],parameter_dict['min_cov_blast'])

def get_PHIS_overall(phage_file,outdir,kind='p'):
	strain_info_dir = os.path.join(outdir,'strain_info')
	module_dir = os.path.join(outdir,'module_results')
	overall_result_file = os.path.join(outdir,'results','overall_result.txt')
	overall_dict = OrderedDict()
	overall_file = os.path.join(outdir,'results','overall_result.txt')
	if os.path.exists(overall_file):
		if os.path.getsize(overall_file)>0:
			with open(overall_file) as f:
				contents = f.readlines()
			for line in contents[1:]:
				line = line.strip().split('\t')
				if kind =='p':
					phage_id = line[0]
					phage_def = line[1]
					host_id = line[2]
					host_def = line[3]
					predict_signal = line[4].split(',')
					score = line[5]
					score = round(float(score),4)
					if phage_id not in overall_dict.keys():
						overall_dict.update({phage_id:{}})
					c_dict = {'phage_id':phage_id,
					'phage_def':phage_def,
					'host_id':host_id,
					'host_def':host_def,
					'crispr':'None',
					'prophage':'None',
					'pro_pro_int':'None',
					'blast':'None',
					'genehomo':'None',
					'score':score}
					if 'crispr' in predict_signal:
						c_dict['crispr'] = 'CRISPR'
					if 'prophage' in predict_signal:
						c_dict['prophage'] = 'Prophage'
					if 'pro_pro_int' in predict_signal:
						c_dict['pro_pro_int'] = 'Protein_Protein_Interaction'
					if 'blast' in predict_signal:
						c_dict['blast'] = 'BLAST'
					if 'genehomo' in predict_signal:
						c_dict['genehomo'] = 'GeneHomology'
					overall_dict[phage_id].update({host_id:c_dict})

	strain_info_dict_file = os.path.join(strain_info_dir,'strain_inf_dict.npy')
	if os.path.exists(strain_info_dict_file):
		strain_info_dict = np.load(strain_info_dict_file).item()
	else:
		strain_info_dict = {}
	phage_length_dict_file = os.path.join(strain_info_dir,'phages_length_dict.npy')
	if os.path.exists(phage_length_dict_file):
		phage_length_dict = np.load(phage_length_dict_file).item()
	else:
		phage_length_dict = {}  
	phage_pro_num_dict_file = os.path.join(strain_info_dir,'protein_num_dict.npy')
	if os.path.exists(phage_pro_num_dict_file):
		phage_pro_num_dict = np.load(phage_pro_num_dict_file).item()
	else:
		phage_pro_num_dict = {}
	
	new_strain_inf_dict = {}
	for phage_id in strain_info_dict.keys():
		new_strain_inf_dict.update({phage_id:{}})
		if phage_id in overall_dict.keys():
			hosts_num = len(overall_dict[phage_id])
		else:
			hosts_num = 0
		try:
			phage_def = strain_info_dict[phage_id]
		except:
			phage_def = phage_id
		try:
			phage_length = phage_length_dict[phage_id]
		except:
			phage_length = "NA"
		try:
			phage_pro_num = phage_pro_num_dict[phage_id]
		except:
			phage_pro_num = "NA"
		new_strain_inf_dict[phage_id].update({"phage_id":phage_id,
			"phage_def":phage_def,
			"host_num":hosts_num,
			"length":phage_length,
			"pro_num":phage_pro_num})
	return new_strain_inf_dict,overall_dict

def create_results_html(phage_file,outdir):
	static_dir = os.path.join(root_path,'static')
	templeate_file = os.path.join(static_dir,'templeates','templeate.html')
	command = "cp -r "+static_dir+" "+outdir
	os.system(command)
	new_strain_inf_dict,overall_dict = get_PHIS_overall(phage_file,outdir)
	write_html = os.path.join(outdir,"PHIS_result.html")
	templeate_read_contents = open(templeate_file).read()
	result_file = os.path.join(outdir,'results','overall_result.txt')
	f_save = open(write_html,'w')
	f_save.write(templeate_read_contents.strip()+'\n')
	write_html_contents = """
	<div class="panel-group" id="accordion">
		<div class="panel panel-danger" style="text-align: left;">
			<div class="panel-heading">
			<h3 class="panel-title">      
			<a href="#">Phage Information</a>
			</h3>
			</div>
			<div class="panel-body">
        <table class="display table table-striped table-hover table-bordered" style="border:0px;">
        <thead>
        <tr>
          <th style='text-align: center;'>Phage_ID</th>
          <th style='text-align: center;'>Phage_Def</th>
          <th style='text-align: center;'>Genome_Size(bp)</th>
          <th style='text-align: center;'>Protein_Number</th>
          <th style='text-align: center;'>Host_Number</th>
        </tr>
      </thead>
      <tbody>
	"""
	write_phage_info_html = "<tr>\n"
	for phage_id in new_strain_inf_dict.keys():
		phage_info = new_strain_inf_dict[phage_id]
		write_phage_info_html = write_phage_info_html+"<td>"+phage_info['phage_id']+"</td>\n"
		write_phage_info_html = write_phage_info_html+"<td>"+phage_info['phage_def']+"</td>\n"
		write_phage_info_html = write_phage_info_html+"<td>"+str(phage_info['length'])+"</td>\n"
		write_phage_info_html = write_phage_info_html+"<td>"+str(phage_info['pro_num'])+"</td>\n"
		write_phage_info_html = write_phage_info_html+"<td>"+str(phage_info['host_num'])+"</td>\n"
	
	write_html_contents = write_html_contents+write_phage_info_html+'''</tr></tbody>
	</table>        
		<p><strong>Annotation:</strong>if phage genome size < 1000bp,the program will not consider this phage.</p>
		</div>
	</div>
	</div>

	'''
	
	write_phage_host_info = '''
	<div class="panel panel-info" style="text-align: left;">
		<div class="panel-heading">
			<h3 class="panel-title">      
			<a href="#"></a>Results Visualization</a>
			</h3>
			</div>
			<div class="panel-body">
	<div style="height: 50px;"></div>
	'''
	
	counter = 1
	for phage_id in overall_dict.keys():
		write_phage_host_info = write_phage_host_info+'''<div>
		<p style="text-align: center;border-bottom: 1px black solid">
		<strong>%s.%s</strong>
		</p>
		<table class ="display" class="table table-striped table-hover table-bordered">
		<thead>
		<tr>
			<th>Host_ID</th>
			<th>Host_Def</th>
			<th>Score</th>
			<th>CRISPR</th>
			<th>Prophage</th>
			<th>GeneHomolog</th>
			<th>Protein_Protein_Interaction</th>
			<th>Oligonucleotide profile/sequence composition</th>
		</tr>
		</thead>
		<tbody>
		'''%(str(counter),phage_id)
		for host_id in overall_dict[phage_id].keys():
			write_phage_host_info = write_phage_host_info+'''<tr style="text-align: center;">
			<td><a href="https://www.ncbi.nlm.nih.gov/nuccore/%s" target="_blank">%s</a></td>
			<td>%s</td>
			<td>%s</td>
			'''%(host_id,host_id,overall_dict[phage_id][host_id]['host_def'],overall_dict[phage_id][host_id]['score'])
			signal_list = ['crispr','prophage','blast','pro_pro_int','genehomo']
			for signal in signal_list:
				if overall_dict[phage_id][host_id][signal] != 'None':
					color = '#F05555'
					label_class = 'glyphicon glyphicon-ok-sign'
					label = 'Yes'
				else:
					color = '#9DE867'
					label_class = 'glyphicon glyphicon-remove-sign'
					label = 'No'
				write_phage_host_info = write_phage_host_info +'''
					<td style="color:%s;font-size: 17px;">
					%s<span class="%s" aria-hidden="true"></span>
					</td>
				'''%(color,label,label_class)
		write_phage_host_info = write_phage_host_info+'''
		</tbody>
		</table>
		</div>
		</div>			
      	</div>

		'''
	write_html_contents = write_html_contents+write_phage_host_info+'''
	</div>
	</div>
	</body>
	</html>
	'''	
	f_save.write(write_html_contents) 
	f_save.close()

if __name__=='__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--input',help='''Path of the input file in genbank or fasta format for bacteria or phages
		please ensure your file name consist of word,number,_,if not,the program will probably run error!''')
	parser.add_argument('--output', help='Path of the output directory.\n')
	parser.add_argument('--module', help='module name(default: all)\n')	
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
	parser.add_argument('--h',help='''print help information''')
	
	args = parser.parse_args()
	global type,root_path
	if args.h:
		print(inputinstructions())
		sys.exit(-1)
	if not args.input:
		print('Error:please input phage file!')
		print(inputinstructions())
		sys.exit(-1)
	if not args.output:
		print('Error:please input output directory!')
		print(inputinstructions())
		sys.exit(-1)
	if args.input:
		input_file = args.input
	if args.output:
		outdir = args.output
	root_path = get_root_path()
	check_file(input_file)
	check_outdir(outdir)
	parameter_dict = check_parameter(args)
	run(input_file,outdir,parameter_dict)
	if parameter_dict['module'] == 'all':
		create_results_html(input_file,outdir)