import networkx as nx
import fishertest
import itertools,sys,os,math,scipy
import numpy as np
import time
import operator,math
from operator import itemgetter
def load_go_term(graph_choice):
	graph_go_term={}
	fisher_one_occurence={}
	fisher_one_occurence2={}

	for i in ["C","P","F","R","K"]:
		graph_go_term[i]={}
		fisher_one_occurence[i]={}
		fisher_one_occurence2[i]={}
	for i in ["C","P","F","R","K"]:
		f1=open("../../web2py_test_new_version/applications/magneto/data/"+graph_choice+"/human/"+i+".txt","r")
		seq=f1.readline()
		while(seq!=""):
			seq=seq.strip().split("\t")
			graph_go_term[i][seq[0]]=seq[1:]
			seq=f1.readline()
		#f1=open("fisher_test/"+graph_choice+"/"+i+"fisher.txt","r")
		#f1=open("inf_cont_new/proteome_ic/"+i+"_count.txt","r")
		f1=open("RECTUM/"+i+"_count.txt","r")
		seq=f1.readline()
		while(seq!=""):
			seq=seq.strip().split("\t")
			fisher_one_occurence[i][seq[0]]=float(seq[1])
			seq=f1.readline()
		f1=open("fisher_test2/"+i+"fisher.txt","r")
		seq=f1.readline()
		while(seq!=""):
			seq=seq.strip().split("\t")
			fisher_one_occurence2[i][seq[0]]=float(seq[1])
			seq=f1.readline()
	return graph_go_term,fisher_one_occurence,fisher_one_occurence2

def load_input_list(file):
	f1=open(file,"r")
	seq=f1.readline()
	uniprot_list=[]
	while(seq!=""):
		seq=seq.strip()
		uniprot_list.append(seq)
		seq=f1.readline()
	return uniprot_list

def expression(file,tissue):
	f1=open(file,"r")
	seq=f1.readline()
	tissue_expr={}
	tissue_value=[]
	while(seq!=""):
		seq=seq.strip().split("\t")
		#load all in one
		tissue_expr[seq[0]]=0.00001
		for j in tissue:
			if seq[int(j)+1]!="0.0":
				tissue_value.append(float(seq[int(j)+1]))
			tissue_expr[seq[0]]=tissue_expr[seq[0]]+float(seq[int(j)+1])
		seq=f1.readline()
	for i in tissue_expr:
		tissue_expr[i]=tissue_expr[i]/len(tissue)
	return tissue_expr,tissue_value
def expression_paxdb(file):
	tissue_expr={}
	f1=open(file,"r")
	seq=f1.readline()
	seq=f1.readline()
	while(seq!=""):
		seq=seq.strip().split("\t")
		tissue_expr[seq[0]]=float(seq[2])
		seq=f1.readline()
	return tissue_expr
def mcn(nodes,graph_nodes,graph,expression,graph_choice,start_nodes):
	graph_go_term,fisher_one_occurence,fisher_one_occurence2=load_go_term(graph_choice)
	path={}
	removable=[]
	path_prob={}
	path_value={}
	path_count={}
	seek={}
	nodes_ic_value={}
	nodes_ic_path={}
	combination=list(itertools.combinations((list(set(graph_nodes).intersection(set(nodes)))),2))
	for i in nodes:
		seek[i]={}
	for i in combination:
		path_value[i]=0.0
		path_count[i]={}
		nodes_ic_value[i]=10
		nodes_ic_path[i]=[]
		f1=open("../../web2py_test_new_version/applications/magneto/data/"+graph_choice+"/path/index/"+i[0]+".txt","r")
		seq=f1.readline()
		while(seq!=""):
			seq= seq.strip().split("\t")
			if seq[0]==i[1]:
				seek[i[0]][seq[0]]=int(seq[1])
			seq=f1.readline()

	to_remove=[]

	for i in seek:
		if len(seek[i])>0:
			for j in seek[i]:
				f1=open("../../web2py_test_new_version/applications/magneto/data/"+graph_choice+"/path/"+i+".txt")
				f1.seek(seek[i][j])
				seq=f1.readline()
				key_start=seq.split("|")[0][1:].strip()
				key_end=seq.split("|")[1].strip()
				key=(key_start,key_end)
				path[key]=[]
				seq=f1.readline()
				count=0
				while(seq[0]!=">"):
					
					path_seq=seq.strip()
					node_list=seq.strip().split()[0:-2]
					#print node_list
					path[key].append(node_list)
					path_coex=[]
					path_exp=[]
					for k in range(len(node_list[0:-1])):
									
						if k==0:
							path_coex.append(graph[node_list[k]][node_list[k+1]]["capacity"])
						if k!=0:
							path_coex.append(graph[node_list[k]][node_list[k+1]]["capacity"])
							if tissue_expr.has_key(node_list[k]):
								path_exp.append(tissue_expr[node_list[k]])
							else:
								path_exp.append(0.0001)
					
					removable.append(node_list[1:-1])
					seq=f1.readline()
					if len(path_exp)==0:
						val1=1.00000
					else:
						val1=scipy.stats.mstats.gmean(path_exp)
					val2=scipy.stats.mstats.gmean(path_coex)
					total_prob=math.sqrt(val1*val2)
					path_count[key][total_prob]=node_list
					#print key,node_list,total_prob,count
					count=count+1
	#flow=[]
	#for j in path_count:
	#	flow.append(path_count[j][sorted(path_count[j])[::-1][0]])
	nodes_ic=[]
	

	#graph_go_term,fisher_one_occurence
	for i in path_count:
		if len(path_count[i])==1:
			nodes_ic_path[i]=path_count[i].values()[0]
			#flow.append(path_count[i].values()[0])
		else:
			count=0	
			#max_score=sorted(path_count[i].keys())[::-1][0]
			mean=np.mean(path_count[i].keys())
			for k in sorted(path_count[i].items(), key=itemgetter(0))[::-1]:
				
				if k[0]>=mean:
					temp_C=[]
					temp_P=[]
					temp_F=[]
					temp_R=[]
					temp_K=[]

					
					for j in k[1]:
						if graph_go_term["P"].has_key(j):
							temp_P.extend(graph_go_term["P"][j])
						if graph_go_term["K"].has_key(j):
							temp_K.extend(graph_go_term["K"][j])
						if graph_go_term["R"].has_key(j):
							temp_R.extend(graph_go_term["R"][j])
						if graph_go_term["C"].has_key(j):
							temp_C.extend(graph_go_term["C"][j])
						if graph_go_term["F"].has_key(j):
							temp_F.extend(graph_go_term["F"][j])
																											

					value_p=0.0
					value_c=0.0
					value_f=0.0
					value_r=0.0
					value_k=0.0
					c=0
					for j in set(temp_P):
						if fisher_one_occurence["P"].has_key(j):
							value_p=value_p+ fisher_one_occurence2["P"][j]
							
						else:
							c=c+1
							#value_p=value_p+ fisher_one_occurence2["P"][j]
					for j in set(temp_C):
						if fisher_one_occurence["C"].has_key(j):
							value_c=value_c+ fisher_one_occurence2["C"][j]
						else:
							c=c+1
							#value_c=value_c+ fisher_one_occurence2["C"][j]
					for j in set(temp_R):
						if fisher_one_occurence["R"].has_key(j):
							value_r=value_r+ fisher_one_occurence2["R"][j]
						else:
							c=c+1
							#value_r=value_r+ fisher_one_occurence2["R"][j]
					for j in set(temp_F):
						if fisher_one_occurence["F"].has_key(j):
							value_f=value_f+ fisher_one_occurence2["F"][j]
							
						else:
							c=c+1
							#value_r=value_r+ fisher_one_occurence2["F"][j]
					for j in set(temp_K):
						if fisher_one_occurence["K"].has_key(j):
							value_k=value_k+ fisher_one_occurence2["K"][j]
						else:
							c=c+1
							#value_r=value_r+ fisher_one_occurence2["K"][j]
					#print k[1],value_p+value_c+value_r+value_f+value_k,float(len(temp_C)+len(temp_P)+len(temp_R)+len(temp_F)+len(temp_K))
																																														
					go_value=(value_p+value_c+value_r+value_f+value_k)/float(len(temp_C)+len(temp_P)+len(temp_R)+len(temp_F)+len(temp_K))
					#go_value=(value_p+value_c+value_f)/float(len(temp_C)+len(temp_P)+len(temp_F))
					#print k[1],value_p+value_c+value_r+value_f+value_k,float(len(temp_C)+len(temp_P)+len(temp_R)+len(temp_F)+len(temp_K)),go_value
					if go_value<nodes_ic_value[i]:
						nodes_ic_value[i]=go_value
						nodes_ic_path[i]=k[1]
						#print i,k[1],go_value
					#print i,k[1],k[0],go_value,float(len(temp_C)+len(temp_P)+len(temp_R)+len(temp_F))
					count=count+1
	#for i in nodes_ic_path:
	#	print i,nodes_ic_path[i],nodes_ic_value[i]
	nodes_ic=list(set(sum(nodes_ic_path.values(),[])+start_nodes))
	
	#flow=list(set(sum(flow,[])+start_nodes))
	print "nodes_ic",len(nodes_ic)
	#print "flow",len(flow)
	fishertest.load(nodes_ic,0.05,["C","P","F","R","K","O","KDr","KDi","DB","Or","VH"],graph_choice,path_def="ic2/",single=sys.argv[3])
	#fishertest.load(flow,0.05,["C","P","F","R","K","O","KDr","KDi","DB","Or","VH"],graph_choice,path_def="flow/",single=sys.argv[3])

pa_sub_loc={0: 'Cytoskeleton (Microtubule plus end)', 1: 'Aggresome', 2: 'Golgi apparatus', 3: 'Nucleoli', 4: 'Focal Adhesions', 5: 'Cytoplasm',
                6: 'Cytoskeleton (Intermediate filaments)', 7: 'Mitochondria', 8: 'Cell Junctions', 9: 'Cytoskeleton (Microtubules)', 10: 'Microtubule organizing center',
                11: 'Cytoskeleton (Actin filaments)', 12: 'Nucleus but not nucleoli', 13: 'Nuclear membrane', 14: 'Endoplasmic reticulum', 15: 'Cytoskeleton (Cytokinetic bridge)',
                16: 'Nucleus', 17: 'Vesicles', 18: 'Centrosome', 19: 'Plasma membrane'}
pa_basal={0: 'stomach 2', 1: 'stomach 1', 2: 'seminal vesicle', 3: 'colon', 4: 'gallbladder', 5: 'esophagus', 6: 'kidney', 7: 'testis', 8: 'tonsil', 9: 'thyroid gland',
        10: 'epididymis', 11: 'adrenal gland', 12: 'endometrium 1', 13: 'endometrium 2', 14: 'salivary gland', 15: 'lateral ventricle', 16: 'smooth muscle', 17: 'heart muscle',
        18: 'pancreas', 19: 'cerebral cortex', 20: 'vagina', 21: 'lymph node', 22: 'placenta', 23: 'hippocampus', 24: 'bronchus', 25: 'parathyroid gland', 26: 'skeletal muscle',
        27: 'urinary bladder', 28: 'appendix', 29: 'liver', 30: 'rectum', 31: 'prostate', 32: 'cervix, uterine', 33: 'duodenum', 34: 'cerebellum', 35: 'fallopian tube',
        36: 'spleen', 37: 'bone marrow', 38: 'soft tissue 2', 39: 'skin 2', 40: 'small intestine', 41: 'soft tissue 1', 42: 'nasopharynx', 43: 'lung', 44: 'ovary',
        45: 'oral mucosa', 46: 'skin 1', 47: 'breast'}
pa_cancer={0: 'lymphoma', 1: 'melanoma', 2: 'glioma', 3: 'head and neck cancer', 4: 'ovarian cancer', 5: 'liver cancer', 6: 'thyroid cancer', 7: 'renal cancer', 
               8: 'skin cancer', 9: 'testis cancer', 10: 'prostate cancer', 11: 'stomach cancer', 12: 'urothelial cancer', 13: 'pancreatic cancer', 14: 'endometrial cancer',
               15: 'lung cancer', 16: 'breast cancer', 17: 'carcinoid', 18: 'cervical cancer', 19: 'colorectal cancer'}
e_mitab_2836={0: 'heart', 1: 'endometrium', 2: 'colon', 3: 'skin', 4: 'esophagus', 5: 'kidney', 6: 'thyroid', 7: 'testis', 8: 'tonsil', 9: 'adrenal gland', 10: 'salivary gland',
                  11: 'smooth muscle', 12: 'duodenum', 13: 'cerebral cortex', 14: 'lymph node', 15: 'stomach', 16: 'bladder', 17: 'adipose tissue', 18: 'skeletal muscle',
                  19: 'appendix', 20: 'liver', 21: 'rectum', 22: 'prostate', 23: 'pancreas', 24: 'placenta', 25: 'animal ovary', 26: 'fallopian tube', 27: 'spleen',
                  28: 'bone marrow', 29: 'small intestine', 30: 'lung', 31: 'gall bladder'}
e_mtab_513={0: 'lymph node', 1: 'heart', 2: 'testis', 3: 'animal ovary', 4: 'leukocyte', 5: 'skeletal muscle', 6: 'adipose tissue', 7: 'adrenal gland', 8: 'brain', 9: 'colon', 
            10: 'breast', 11: 'thyroid', 12: 'lung', 13: 'prostate', 14: 'kidney', 15: 'liver'}
e_geod_30352={0:"cerebellum",1:"frontal lobe",2:"heart",3:"kidney",4:"liver",5:"prefrontal cortex",6:"temporal lobe",7:"testis"}
e_mtab_3358_adult={0: 'heart', 1: 'skin', 2: 'locus coeruleus', 3: 'brain', 4: 'seminal vesicle', 5: 'colon', 6: 'spinal cord', 7: 'mitral valve', 8: 'substantia nigra',
                       9: 'cerebral meninges', 10: 'kidney', 11: 'putamen', 12: 'penis', 13: 'cervix', 14: 'testis', 15: 'pituitary gland', 16: 'triscuspid valve', 
                       17: 'diencephalon', 18: 'parietal lobe', 19: 'epididymis', 20: 'vagina', 21: 'pulmonary valve', 22: 'pineal gland', 23: 'globus pallidus',
                       24: 'lymph node', 25: 'placenta', 26: 'middle frontal gyrus', 27: 'hippocampus', 28: 'medulla oblongata', 29: 'tongue', 30: 'parotid gland',
                       31: 'left ventricle', 32: 'appendix', 33: 'breast', 34: 'left atrium', 35: 'artery', 36: 'prostate', 37: 'caudate nucleus', 38: 'amygdala', 
                       39: 'pancreas', 40: 'cerebellum', 41: 'animal ovary', 42: 'olfactory apparatus', 43: 'vas deferens', 44: 'spleen', 45: 'occipital lobe', 46: 'bone marrow',
                       47: 'submandibular gland', 48: 'thalamus', 49: 'occipital cortex', 50: 'smooth muscle', 51: 'uterus', 52: 'lung', 53: 'gall bladder', 54: 'dura mater',
                       55: 'middle temporal gyrus'}
e_mtab_3358_fetal={0: 'thyroid', 1: 'trachea', 2: 'eye', 3: 'diaphragm', 4: 'parietal lobe', 5: 'skeletal muscle', 6: 'occipital lobe', 7: 'throat', 8: 'small intestine',
                        9: 'colon', 10: 'uterus', 11: 'spinal cord', 12: 'skin', 13: 'stomach', 14: 'lung', 15: 'rectum', 16: 'temporal lobe', 17: 'tongue', 18: 'umbilical cord',
                        19: 'duodenum'}
GTEx={'42': 'Kidney - Cortex', '48': 'Colon - Sigmoid', '43': 'Brain - Cerebellum', '49': 'Nerve - Tibial', '24': 'Brain - Hypothalamus', '25': 'Artery - Aorta', 
	  '26': 'Prostate','27': 'Brain - Amygdala', '20': 'Small Intestine - Terminal Ileum', '21': 'Artery - Coronary', '22': 'Liver', '23': 'Esophagus - Gastroesophageal Junction',
	  '46': 'Lung', '47': 'Muscle - Skeletal', '44': 'Adipose - Visceral (Omentum)', '45': 'Adrenal Gland', '28': 'Pancreas', '29': 'Adipose - Subcutaneous', 
	  '40': 'Minor Salivary Gland', '41': 'Whole Blood', '1': 'Testis', '0': 'Thyroid', '3': 'Ovary', '2': 'Skin - Not Sun Exposed (Suprapubic)', '5': 'Vagina', 
	  '4': 'Esophagus - Muscularis', '7': 'Brain - Anterior cingulate cortex (BA24)', '6': 'Heart - Atrial Appendage', '9': 'Pituitary', '8': 'Brain - Putamen (basal ganglia)',
	  '39': 'Brain - Caudate (basal ganglia)', '38': 'Brain - Nucleus accumbens (basal ganglia)', '11': 'Cervix - Ectocervix', '10': 'Breast - Mammary Tissue', 
	  '13': 'Fallopian Tube', '12': 'Stomach', '15': 'Cervix - Endocervix', '14': 'Brain - Frontal Cortex (BA9)', '17': 'Esophagus - Mucosa', '16': 'Colon - Transverse',
	  '19': 'Brain - Cerebellar Hemisphere', '18': 'Bladder', '31': 'Spleen', '30': 'Skin - Sun Exposed (Lower leg)', '37': 'Uterus', 
	  '36': 'Brain - Spinal cord (cervical c-1)', '35': 'Artery - Tibial', '34': 'Brain - Cortex', '33': 'Heart - Left Ventricle', '32': 'Brain - Hippocampus',
	  '50': 'Brain - Substantia nigra'}

graph=nx.read_gpickle(sys.argv[1])
graph_nodes=graph.nodes()
start_nodes=load_input_list(sys.argv[2])
graph_choice="intact"
print "starting",len(start_nodes)
fishertest.load(list(set(start_nodes)),0.05,["C","P","F","R","K","O","KDr","KDi","DB","Or","VH"],graph_choice,path_def="results_old/",single="0")

nodes=list(set(start_nodes).intersection(set(graph_nodes)))

folder="tissue_expr_norm/"
for i in range(0,48,1):
	print str(i)+" "+pa_basal[i]
val=raw_input()
if "," not in val:
	val=map(int,str(val).split())
else:
	print val
	val=map(int, val.split(","))
tissue_expr,tissue_value=expression(folder+graph_choice+"/PA_basal.txt",val)

tissue_expr=expression_paxdb("RECTUM.txt")


mcn(nodes,graph.nodes(),graph,tissue_expr,graph_choice,start_nodes)

