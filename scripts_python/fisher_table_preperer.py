#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always




import sys


def values_for_fisher_computer(in_set, total_inset, backgoundset, total_backgr, column_interest=0):
	#print('bubba', inset, total_inset, backgoundset, total_backgr, column_interest)
	with open(in_set, 'r') as inset:
		for inline in inset:
			if inline[0] == '#':
				continue
			else:
				annotation = inline.split()
				a_val = int(annotation[int(column_interest)])
				del annotation[int(column_interest)]
				#print('a_val:', a_val, 'annot:', annotation)
				c_val = 0
				str_annotation = ' '.join(annotation)
				with open(backgoundset, 'r') as backgound:
					for backline  in backgound:
						if annotation[0] in backline and annotation[-1] in backline:
							back_annotation = backline.split()
							tmp = int(back_annotation[int(column_interest)])
							del back_annotation[int(column_interest)]
							if annotation == back_annotation:	# exactstring match without spaces or tab
								c_val = tmp
								#print('c_val_back:', c_val, backline, end='')
								break
				if c_val == 0:
					print('this annotation has not been found in the backgound set:\t', str_annotation, '\tc value was set to zero')
				if (int(total_inset) - a_val) < 0:
					print('this annotation has avalue greater than the total given:\t', str_annotation, '\tvalue:\t', a_val, '\ttotal input set given:\t', total_inset, file=sys.stderr)
					raise SystemExit
				if (int(total_backgr) - c_val) < 0:
					print('this annotation has avalue greater than the total given:\t', str_annotation, '\tvalue:\t', c_val, '\ttotal backgroun given:\t', total_backgr, file=sys.stderr)
					raise SystemExit
				#print('a_val:', a_val, 'b_val:', (int(total_inset) - a_val), 'c_val:', c_val, 'd_val:', (int(total_backgr) - c_val), 'annot:', annotation)
				output_filename = in_set.split('/')[-1].split('.')[0] + backgoundset.split('/')[-1].split('.')[0] + '.ftable'
				with open(output_filename, 'a') as outfile:
					outfile.write(str(a_val) + ' ' + str(int(total_inset) - a_val) + ' ' + str(c_val) + ' ' + str(int(total_backgr) - c_val) + ' ' + str_annotation + '\n')




if __name__ == "__main__":
	total_background = None
	try:
		input_set_file = sys.argv[1]
		total_inputset = sys.argv[2]
		background_set = sys.argv[3]
		total_background = sys.argv[4]
		column_interest = sys.argv[5]
	except Exception:
		if total_background is not None: 
			values_for_fisher_computer(input_set_file, total_inputset, background_set, total_background)
		else:
			print('Program usage: text.py <a space or tab separated file with on a given column has integers associated to a annotation (string) on that same line, aline would look like -> 24410 region disordered, the annotation that is everything but the column selected (default 0, first column) is used to search on the background file, so that the values in the input file are associated by the annotation to thoose in the background, when not found the script wil print to screen and assign a 0 value to the third output number, the search is an exact string match one, so is up to the user to make sure the annotations are written in the same way in the two files, if you need to generate the files above look at    <uniprot_annotation_domain_processer.py >   > < the total as integer necessary to compute the non-annotation value. on every line of the input file there is a value specific for a feature, but fisher needs two features and in this script the second feature is the absence of the feature itself, shortly   total - value    will give the second number/feature value  > <background file with the same architecture of the input one (same relevant column)> < total associted to the values found in the bacgrounds > <optional flag, column to take into consideration to find the integer associated to an annotation, this follows python annotation where 0 is the first column> <the output file will have the same number of lines of the input set file but in this format each line: value_feature_in_inset value_not-freature_in_inset value_same-feature_in_background value_not_same-feature_in_background feature name, like -> 62 21 34 54 disorder region, take also a look at	 annotation_enrichment_analizer.nf 	pipeline >', file=sys.stderr)
			raise SystemExit
	else:
		values_for_fisher_computer(input_set_file, total_inputset, background_set, total_background, column_interest)
