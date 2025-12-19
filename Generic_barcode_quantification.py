from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import numpy as np
from Bio import pairwise2
import pandas as pd
from collections import Counter
from Bio import SeqIO
import csv
import gzip

def compile_experiment(folder_path, raw_barcodes_only = False, cheat_loc_tests = 5000, include_umi = False, write_umi = False, write_bc_umi_count_dict = False, shorten = False, extrap_method = 'exponential_decay', length_flex = 0, cheat_locs = None):
	# Requires 5 things in folder_path:
	# 1: a folder called 'Fastq_files' containing all the fastq files
	# 2: a file called 'Barcode_info.csv' containing each barcode and its corresponding ID name. Column names: Barcode_ID, Barcode_seq
	# 3: a file called 'Fastq_file_IDs.csv' containing each fastq filename and its corresponding descriptor. Column names: Filename, Sample_label
	# 4: a file called 'Read_tags.csv' containing each read tag and the read tag's ID. Column names: Read_tag, Tag_descriptor
	# 5: a file called 'Flanks.csv' containing the left and right sequence flanks. Column names: Left_flank, Right_flank
	# Outputs: a file called 'Quantified_reads'
	read_tag_df = pd.read_csv(folder_path + 'Read_tags.csv')
	flank_df = pd.read_csv(folder_path + 'Flanks.csv')
	print(flank_df)
	fastq_df = pd.read_csv(folder_path + 'Fastq_file_IDs.csv')
	left_flank = flank_df.Left_flank[0]
	right_flank = flank_df.Right_flank[0]
	if not raw_barcodes_only:
		barcode_info_df = pd.read_csv(folder_path + 'Barcode_info.csv')
		all_dict = dict(barcode_info_df)
	if include_umi:
		umi_all_dict = dict(barcode_info_df)
		fastq_df = pd.read_csv(folder_path + 'Fastq_file_IDs_UMI.csv')
		umi_flank_df = pd.read_csv(folder_path + 'UMI_flanks.csv')
		umi_left_flank = umi_flank_df.Left_flank[0]
		umi_right_flank = umi_flank_df.Right_flank[0]
		for i,fname in enumerate(fastq_df['Filename']):
			print('Starting filename: ',fname)
			umi_dict,barcodes,bc_umi_count_dict = read_file_umi(folder_path +'Fastq_files/' + fname + '.gz',barcode_info_df, fastq_df['Sample_label'][i],left_flank,right_flank,read_tag_df,umi_left_flank, umi_right_flank, cheat_loc_tests = cheat_loc_tests, shorten = shorten)
			all_dict.update(barcodes)
			umi_all_dict.update(bc_umi_count_dict_to_extrapolated(bc_umi_count_dict,barcodes,read_tag_df,barcode_info_df,extrap_method))
			if write_umi:
				f = open(folder_path+'umi_dict_'+fname+'.csv','w')
				writer = csv.writer(f)
				for key in umi_dict.keys():
					to_write = [key]
					for val in umi_dict[key]:
						to_write .append(val)
					# print(to_write)
					writer.writerow(to_write)
				f.close()
			if write_bc_umi_count_dict:
				bcumidf = pd.DataFrame.from_dict(bc_umi_count_dict)
				bcumidf.to_csv(folder_path+'UMI_counts_by_barcode_'+fname+'.csv')
		df = pd.DataFrame.from_dict(all_dict)
		df.to_csv(folder_path + 'Quantified_reads_UMI.csv',index = False)
		umi_df = pd.DataFrame.from_dict(umi_all_dict)
		umi_df.to_csv(folder_path + 'Quantified_reads_UMI_extrapolated.csv',index = False)

	else:
		if raw_barcodes_only:
			for i,fname in enumerate(fastq_df['Filename']):
				print('Starting filename: ',fname)
				raw_barcodes = read_file(folder_path +'Fastq_files/' + fname + '.gz',[], fastq_df['Sample_label'][i],left_flank,right_flank,read_tag_df, cheat_loc_tests, cheat_locs,raw_only = raw_barcodes_only)
				df = pd.DataFrame.from_dict(raw_barcodes)
				df.to_csv(folder_path + fastq_df['Sample_label'][i]+'_raw_barcodes.csv')
			return {}

		barcode_info_df = pd.read_csv(folder_path + 'Barcode_info.csv')
		all_dict = dict(barcode_info_df)

		# print(read_tag_df)
		# print(fastq_df)
		# print(barcode_info_df)
		for i,fname in enumerate(fastq_df['Filename']):
			print('Starting filename: ',fname)
			raw_barcodes = read_file(folder_path +'Fastq_files/' + fname + '.gz',barcode_info_df, fastq_df['Sample_label'][i],left_flank,right_flank,read_tag_df, cheat_loc_tests,cheat_locs, raw_only = raw_barcodes_only, length_flex = length_flex)
			# df = pd.DataFrame.from_dict(raw_barcodes)
			# df.to_csv(folder_path+fastq_df['Sample_label'][i]+'_raw_barcodes.csv')
			all_dict.update(raw_barcodes)
		df = pd.DataFrame.from_dict(all_dict)
		df.to_csv(folder_path + 'Quantified_reads.csv',index = False)
	return all_dict


def extrapolate_val(x,y,method):
	if method == 'exponential_decay':
		try:
			slope, intercept = np.polyfit(x,np.log(y),1, w = np.sqrt(y))
			extrap = -np.exp(intercept)/(np.exp(slope)-1)
		except:
			extrap = 0

	return extrap

def bc_umi_count_dict_to_extrapolated(count_dict, std_barcode_dict, tag_df,barcode_info_df,extrap_method, min_count = 10):
	# count_dict: a dataframe with columns as "tag_barcode" and rows as "for this tag and barcode combo, how many UMI's have this many reads"
	# min_count: only use things with this count or more
	results_dict = {key:value[:] for key,value in std_barcode_dict.items()}
	tag_dict = {}
	for i,tag in enumerate(tag_df.Read_tag):
		tag_dict[tag] = tag_df.Tag_descriptor[i]
	# read_tags = tag_df.Read_tag
	# descriptors = tag_df.Read_descriptor
	# for des in descriptors:
		# results_dict[des] = [0]*len(barcode_df.Barcode_seq)
	for key in count_dict.keys():
		allthings = count_dict[key]
		x_vals = []
		y_vals = []
		for i,val in enumerate(allthings):
			if val>=min_count:
				x_vals.append(i)
				y_vals.append(val)
		extrapolated_val = extrapolate_val(x_vals,y_vals,extrap_method)
		info = key.split('_')
		tag = info[0]
		tag_desc = tag_dict[tag]
		col = [c for c in std_barcode_dict.keys() if tag_desc in c][0]
		barcode = info[1]
		try:
			row = list(barcode_info_df['Barcode_seq']).index(barcode)
			results_dict[col][row] = extrapolated_val
		except:
			pass
	return results_dict
		# Tag_descriptor
		# Read_tag

def read_file_umi(record_fname, barcode_df, fastq_id, left_flank, right_flank, tag_df, left_flank_umi, right_flank_umi, cheat_loc_tests, shorten = False):
	# gzip.open("practicezip.fasta.gz", "rt")
	# records = list(SeqIO.parse(record_fname,'fastq'))
	records = list(SeqIO.parse(gzip.open(record_fname, "rt"),'fastq'))
	if shorten:
		records = records[:5*cheat_loc_tests]
	umi_dict, raw_barcode_dict, raw_barcode_count_dict = read_records_to_raw_barcode_dict_umi(records, left_flank, right_flank, tag_df, cheat_loc_tests, left_flank_umi, right_flank_umi)
	bc_with_counts, bc_count_dict = raw_barcodes_to_barcode_dict_umi(raw_barcode_dict, barcode_df, fastq_id, raw_barcode_count_dict)
	return umi_dict,bc_with_counts, bc_count_dict

def read_records_to_raw_barcode_dict_umi(records, left_flank, right_flank, tag_df, cheat_loc_tests, left_flank_umi, right_flank_umi, max_umi_count = 100, min_phred = 10):
	cheat_coords = []
	for record in records[:cheat_loc_tests]:
	    bc_loc = raw_barcode_loc(record, left_flank, right_flank)
	    if len(bc_loc) == 2:
	        cheat_coords.append((bc_loc[0],bc_loc[1]))
	cheat_coords = set(cheat_coords)
	print('Cheat coords: ',cheat_coords)

	cheat_coords_umi = []
	for record in records[:cheat_loc_tests]:
	    bc_loc = raw_barcode_loc(record, left_flank_umi, right_flank_umi)
	    if len(bc_loc) == 2:
	        cheat_coords_umi.append((bc_loc[0],bc_loc[1]))
	cheat_coords_umi = set(cheat_coords_umi)
	print('Cheat coords UMI: ',cheat_coords_umi)

	tag_dict = {}
	for i,tag in enumerate(tag_df.Read_tag):
		tag_dict[tag] = tag_df.Tag_descriptor[i]

	umi_dict = {}
	raw_barcodes = {}
	for tag in tag_dict.keys():
		raw_barcodes[tag_dict[tag]] = {}
	for record in records:
		for cheat in cheat_coords:
			tag, raw_barcode = read_to_tagged_barcode(record, left_flank, right_flank, tag_dict.keys(), cheat_left = cheat[0], cheat_right = cheat[1], min_phred = min_phred)
			if not (tag=='' or raw_barcode ==''):
				for cheat_umi in cheat_coords_umi:
					_, umi = read_to_tagged_barcode(record, left_flank_umi, right_flank_umi, tag_dict.keys(), cheat_left = cheat_umi[0], cheat_right = cheat_umi[1])
					if not umi=='':
						if not umi in umi_dict.keys():
							umi_dict[umi] = [(tag,raw_barcode)]
						else:
							umi_dict[umi].append((tag,raw_barcode))
					else:
						pass
						# if np.random.randint(1000)<3:
						# 	print('Couldnt find UMI for sequence: ',record.seq)
						# 	print('Flanks: ',left_flank_umi,' and ',right_flank_umi)
						# if raw_barcode in raw_barcodes[tag].keys():
							# raw_barcodes[tag][raw_barcode].append(umi)
						# else:
							# raw_barcodes[tag][raw_barcode] = [umi]
	good_umi_dict = {}
	barcode_count_dict = {}
	for key in umi_dict.keys():
		if len(set(umi_dict[key]))==1:
			tag, raw_barcode = umi_dict[key][0]
			if raw_barcode in raw_barcodes[tag_dict[tag]].keys():
				raw_barcodes[tag_dict[tag]][raw_barcode].append(key)
			else:
				raw_barcodes[tag_dict[tag]][raw_barcode] = [key]
			# raw_barcodes[umi_dict[key][0][0]][umi_dict[key][0][1]].append(key)
			good_umi_dict[key] = umi_dict[key]
			bcd_key = tag+'_'+raw_barcode
			if not bcd_key in barcode_count_dict.keys():
				barcode_count_dict[bcd_key] = [0]*max_umi_count
			barcode_count_dict[bcd_key][min(max_umi_count-1,len(umi_dict[key]))] += 1
		# else:
			# print('blah multiple UMI maps: ',umi_dict[key])
	return good_umi_dict, raw_barcodes, barcode_count_dict

def raw_barcodes_to_barcode_dict_umi(raw_barcode_dict, barcode_df, fastq_id, raw_barcode_count_dict):
	# currently doesn't do any error correcting!
	# Barcode_df: 2 columns: 'Barcode_ID', 'Barcode_seq'
	to_return = {}
	for tag in raw_barcode_dict.keys():
		barcode_counts = []
		for barcode in barcode_df.Barcode_seq:
			if barcode in raw_barcode_dict[tag].keys():
				barcode_counts.append(len(set(raw_barcode_dict[tag][barcode])))
			else:
				barcode_counts.append(0)
		to_return[fastq_id+'_'+tag] = barcode_counts
	to_del = []
	for key in raw_barcode_count_dict.keys():
		keep = False
		for barcode in barcode_df.Barcode_seq:
			if barcode in key:
				keep = True
		if not keep:
			to_del.append(key)
	for key in to_del:	
		raw_barcode_count_dict.pop(key)
	return to_return, raw_barcode_count_dict



def read_file(record_fname, barcode_df, fastq_id, left_flank, right_flank, tag_df, cheat_loc_tests,cheat_locs, raw_only = False, length_flex = 0):
	# records = list(SeqIO.parse(record_fname,'fastq'))
	records = list(SeqIO.parse(gzip.open(record_fname, "rt"),'fastq'))
	raw_barcode_dict = read_records_to_raw_barcode_dict(records, left_flank, right_flank, tag_df, cheat_loc_tests, length_flex,cheat_locs)
	if raw_only:
		return raw_barcode_dict
	bc_with_counts = raw_barcodes_to_barcode_dict(raw_barcode_dict, barcode_df, fastq_id)
	return bc_with_counts

def read_records_to_raw_barcode_dict(records, left_flank, right_flank, tag_df, cheat_loc_tests, length_flex,cheat_locs):
	if cheat_locs == None:
		cheat_coords = []
		for record in records[:cheat_loc_tests]:
		    bc_loc = raw_barcode_loc(record, left_flank, right_flank)
		    if len(bc_loc) == 2:
		    	if length_flex>0:
		    		for ch in range(max(bc_loc[0],bc_loc[1]-length_flex),bc_loc[1]+length_flex):
		    			cheat_coords.append((bc_loc[0],ch))
		    	else:
			        cheat_coords.append((bc_loc[0],bc_loc[1]))
		cheat_coords = set(cheat_coords)
	else:
		cheat_coords = set(cheat_locs)
	print('cheat_coords: ',cheat_coords)

	tag_dict = {}
	for i,tag in enumerate(tag_df.Read_tag):
		tag_dict[tag] = tag_df.Tag_descriptor[i]

	raw_barcodes = {}
	for tag in tag_dict.keys():
		raw_barcodes[tag_dict[tag]] = []
	for record in records:
		for cheat in cheat_coords:
			tag, raw_barcode = read_to_tagged_barcode(record, left_flank, right_flank, tag_dict.keys(), cheat_left = cheat[0], cheat_right = cheat[1])
			if not (tag=='' or raw_barcode ==''):
				# print(tag)
				# print(tag_dict)
				raw_barcodes[tag_dict[tag]].append(raw_barcode)
	# print(raw_barcodes)
	to_return = {}
	for key in tag_dict.keys():
		# print(raw_barcodes[key])
		to_return[tag_dict[key]] = dict(Counter(raw_barcodes[tag_dict[key]]))
	return(to_return)

def raw_barcodes_to_barcode_dict(raw_barcode_dict, barcode_df, fastq_id):
	# currently doesn't do any error correcting!
	# Barcode_df: 2 columns: 'Barcode_ID', 'Barcode_seq'
	to_return = {}
	for tag in raw_barcode_dict.keys():
		barcode_counts = []
		for barcode in barcode_df.Barcode_seq:
			if barcode in raw_barcode_dict[tag].keys():
				barcode_counts.append(raw_barcode_dict[tag][barcode])
			else:
				barcode_counts.append(0)
		to_return[fastq_id+'_'+tag] = barcode_counts
	return to_return

def raw_barcode_loc(seq_record, left_flank, right_flank):
	to_return = []
	raw_barcode = ''
	phred_scores = seq_record.letter_annotations['phred_quality']
	read = seq_record.seq
	strread = str(read)
	fwd_seqs = []
	rev_seqs = []
	fwd_seqs = seq_between_flanks_forward(read, left_flank, right_flank)
	rev_seqs = seq_between_flanks_forward(read, Seq(right_flank).reverse_complement(), Seq(left_flank).reverse_complement())
	# print(rev_seqs)
	to_return = fwd_seqs
	if len(rev_seqs)>1:
		# print('here')
		# print(read)
		# print(rev_seqs)
		to_return = [-1*rev_seqs[0],-1*rev_seqs[1]]
	return to_return

def read_to_tagged_barcode(seq_record, left_flank, right_flank, read_tags, max_errors = 2, min_phred = 10, cheat_left = 0, cheat_right = 0):
	# tag_len = 
	strread = str(seq_record.seq)
	readtag = strread[:len(list(read_tags)[0])]
	tag = ''
	for temp_tag in read_tags:
		if temp_tag == readtag:
			tag = temp_tag
	if tag == '':
		return '',[]
	return tag, read_to_raw_barcode(seq_record, left_flank, right_flank, max_errors = max_errors, min_phred = 10, cheat_left = cheat_left, cheat_right = cheat_right)

def read_to_raw_barcode(seq_record, left_flank, right_flank, max_errors = 2, min_phred = 10, cheat_left = 0, cheat_right = 0):
	raw_barcode = ''
	phred_scores = seq_record.letter_annotations['phred_quality']
	read = seq_record.seq
	strread = str(read)
	fwd_seqs = []
	rev_seqs = []
	if cheat_left>0:
		if len(strread) > cheat_right + len(right_flank):
			if strread[cheat_left-len(left_flank):cheat_left] == left_flank and strread[cheat_right:cheat_right + len(right_flank)] == right_flank:
				fwd_seqs = [cheat_left, cheat_right]
	elif cheat_left<0:
		cheat_left = -cheat_left
		cheat_right = -cheat_right
		# print(cheat_left)
		# print(cheat_right)
		if len(strread) > cheat_right + len(right_flank):
			if strread[cheat_left-len(right_flank):cheat_left] == str(Seq(right_flank).reverse_complement()) and strread[cheat_right:cheat_right + len(left_flank)] == str(Seq(left_flank).reverse_complement()):
				rev_seqs = [cheat_left, cheat_right]
	else:
		fwd_seqs = seq_between_flanks_forward(read, left_flank, right_flank)
		# print(degenerate_code)
		rev_seqs = seq_between_flanks_forward(read, Seq(right_flank).reverse_complement(), Seq(left_flank).reverse_complement())
	if len(fwd_seqs)+len(rev_seqs)>4:
		print('Too many sequences: fwd_seqs = ',fwd_seqs,' and rev_seqs = ',rev_seqs)
	elif len(fwd_seqs) == 2:
		# print('here!')
		barcode_phred = phred_scores[fwd_seqs[0]:fwd_seqs[1]]
		# Now just create a single bad datapoint if barcode_phred is (inexplicably) empty:
		if len(barcode_phred)<1:
			barcode_phred = [min_phred-5]
		# print(read)
		# print(phred_scores)
		# print(barcode_phred)
	elif len(rev_seqs) == 2:
		barcode_phred = phred_scores[rev_seqs[0]:rev_seqs[1]]
		# Now just create a single bad datapoint if barcode_phred is (inexplicably) empty:
		if len(barcode_phred)<1:
			barcode_phred = [min_phred-5]
	if len(fwd_seqs)+len(rev_seqs)==2:
		# print('fwd_seqs: ',fwd_seqs)
		if min(barcode_phred) >= min_phred:
			if len(fwd_seqs) == 2:
				raw_barcode = read[fwd_seqs[0]:fwd_seqs[1]]
			elif len(rev_seqs) == 2:
				raw_barcode = read[rev_seqs[0]:rev_seqs[1]].reverse_complement()
	return str(raw_barcode)

# 49, 64 for reverse
# 109, 124 for forward

def seq_between_flanks_forward(read, left_flank, right_flank, max_flank_errors = 0, print_multiple_instances = False):
	# Returns sequence and also where in the sequence it was found (to check for local phred score)
	left_loc = seq_loc(read, left_flank, max_flank_errors)
	if len(left_loc)>0:
		right_loc = seq_loc(read, right_flank, max_flank_errors)
		if len(right_loc) > 0:
			if len(right_loc)>1:
				if print_multiple_instances:
					print('Multiple instances of right flank: ',right_flank)
			elif len(left_loc)>1:
				if print_multiple_instances:
					print('Multiple instances of left flank: ', left_flank)
			else:
				left_edge = left_loc[0]+len(left_flank)
				right_edge = right_loc[0]
				return [left_edge, right_edge]
	return []


def seq_loc(read, target_seq, max_errors):
	# Returns sequence and also where in the sequence it was found (to check for local phred score)
	seqs_to_return = []
	readstr = str(read)
	codestr = str(target_seq)
	match_score = 1
	target_score = len(target_seq) - max_errors
	# mismatch_score = float(len(degenerate_code))/max_errors + 0.001
	mismatch_score = 0
	# Note: NO GAPS ALLOWED

	score_matrix = [[0]*(len(codestr)+1) for i in range(len(readstr)+1)]
	for i in range(len(readstr)):
		for j in range(len(codestr)):
			score_matrix[i+1][j+1] = score_matrix[i][j]
			# print(match_dict[codestr[j]])
			if readstr[i] == codestr[j]:
				score_matrix[i+1][j+1] += match_score
	for i in range(len(target_seq),len(readstr)):
		if score_matrix[i][len(codestr)] >= target_score:
			# print(i)
			seqs_to_return.append(i-len(target_seq))
	return seqs_to_return











