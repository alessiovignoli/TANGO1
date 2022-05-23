#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always

# take a look at phobius_short_prediction_field_retriever.py  it is called to create the most important inputs for this script
# take a look at the print error for the call of the function



import sys
import os
sys.path.append(os.path.abspath("/home/alessio/Desktop/erasmus-internship/scripts"))
from phobius_short_prediction_field_retriever import phobius_short_pred_field_selecter



def number_and_pos_computer(short_pred_file, input_field):
	safety_check = ["s", "n"]
	if input_field not in safety_check:
		print("the allowed keywords are: ['s', 'n']",  file=sys.stderr)
		raise SystemExit
	dict_of_num = {}
	dict_of_pos = {}
	with open(short_pred_file, 'r') as infile:
		for line in infile:
			seq_id, list_boundaries = phobius_short_pred_field_selecter(line, input_field)
			#print(seq_id, list_boundaries)
			#print(len(list_boundaries))
			if len(list_boundaries) in dict_of_num:
				dict_of_num[len(list_boundaries)] += 1
			else:
				dict_of_num[len(list_boundaries)] = 1
			other_feature = None
			if input_field=="s":
				other_feature = "n"
			else:
				other_feature = "s"
			bubba, list_boundaries_other = phobius_short_pred_field_selecter(line, other_feature)
			#print(list_boundaries_other, len(list_boundaries_other))
			for n, elem in enumerate(list_boundaries):
				for i, comapare_to in enumerate(list_boundaries_other):
					if elem[0] > comapare_to[1]:
						#print('boundary :', elem, '\tleft_extremity :', elem[0],  '\tright_extremity :', comapare_to[1], '\ti :', i)
						continue
					else:
						#print('####  boundary :', elem, '\tleft_extremity :', elem[0],  '\tright_extremity :', comapare_to[1], '\ti :', i, '\tn :', n, 'count :', n+i+1)
						if (n+i+1) in dict_of_pos:
							dict_of_pos[(n+i+1)] += 1
						else:
							dict_of_pos[(n+i+1)] = 1
						break

	print(dict(sorted(dict_of_num.items(), key=lambda item: item[0])))
	print(dict(sorted(dict_of_pos.items(), key=lambda item: item[0])))
	


if __name__ == "__main__":
	try:
		input_short_pred_file = sys.argv[1]
		feature = sys.argv[2]
	except Exception:
		print('Program usage: text.py <a phobius short output redirection, where on each line is present the prediction and position of the latter on the sequence> <the one letter code for the TM, possible = [s, n], take a look at phobius_short_prediction_field_retriever.py  error and help message for more details, \n the scripts computes how many times a feature present in the prediction and if it is the first second ecc.. of the TM in the protein>', file=sys.stderr)
		raise SystemExit
	else:
		number_and_pos_computer(input_short_pred_file, feature)
