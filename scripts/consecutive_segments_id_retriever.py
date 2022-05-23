#!/usr/bin/env python3


# tm = trans membrane            always

# take a look at phobius_short_prediction_field_retriever.py  it is called to create the most important inputs for this script
# take a look at the print error for the call of the function



import sys
import os
sys.path.append(os.path.abspath("/home/alessio/Desktop/erasmus-internship/scripts"))
from phobius_short_prediction_field_retriever import phobius_short_pred_field_selecter


def consecutive_fields_retriever(input_txt, output_file, field, number_consec=2, distanciator=25):
	#print(input_txt, output_file, field, number_consec, distanciator)
	with open(input_txt, 'r') as in_txt, open(output_file, 'w') as outfile:
		for line in in_txt:
			seq_id, list_boundaries = phobius_short_pred_field_selecter(line, field)
			#print(list_boundaries)
			if len(list_boundaries) == int(number_consec):	# here can set to more than x but changes ahve to be made also downstream
				#print(list_boundaries)
				distance_below_threshold = True
				for n, feature_ranges in enumerate(list_boundaries[:-1]):
					#print(n, feature_ranges)
					if (list_boundaries[(n+1)][0] - feature_ranges[1]) <= int(distanciator):	#checking distance between istances
						#print("yes, distance calculi:", list_boundaries[(n+1)][0], "-", feature_ranges[1] )
						continue
					else:
						#print("no, distance calculi:", list_boundaries[(n+1)][0], "-", feature_ranges[1] )
						distance_below_threshold = False
				#print(distance_below_threshold)
				if distance_below_threshold:
					outfile.write(seq_id + "\n")



if __name__ == "__main__":
	input_field = None
	try:
		input_txt_file = sys.argv[1]
		output_filename = sys.argv[2]
		input_field = sys.argv[3]
		input_number_consec = sys.argv[4]
		spatiator = sys.argv[5]
	except Exception:
		if input_field is not None:
			consecutive_fields_retriever(input_txt_file, output_filename, input_field)
		else:
			print('Program usage: text.py <a phobius short output redirection file with ids and prediction on same line, with n lines> <the nameor path for the output file> <the letter code for selecting with featur of the prediction to look for, the scripts check whether in the prediction there are exactly input_number_consec number of consecutive requested elements the threshold isset to 25 but it can be changed with spatiator>  <optional field, the number in digit of exact instances to obe found in the prediction, default 2, this field rwquires for spatiator to be specified to (following)> <optional field, the threshold to consider consecutive two features default 25, this field requires the input_number_consec to be specified as well>', file=sys.stderr)
			raise SystemExit
	else:
		consecutive_fields_retriever(input_txt_file, output_filename, input_field, input_number_consec, spatiator)
