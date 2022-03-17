#!/usr/bin/env python3

# pp = posterior probability     always
# tm = trans membrane            always




import sys
from scipy.stats import fisher_exact


def fisherer(input_table, mode_keyword):
	print(input_table, mode_keyword)
	check_modes = ["two-sided", "greater", "less"]
	if mode_keyword not in check_modes:
		print('the allowed keyword for modes of fidsher exact test are : "two-sided", "greater", "less", check for script error message for more details', file=sys.stderr)
		raise SystemExit
	with open(input_table) as infile:
		for line in infile:
			table = []
			table.append(line.split()[:2])
			table.append(line.split()[2:4])
			annotation_name = ' '.join(line.split()[4:])
			print(table, annotation_name)
			oddsr, p = fisher_exact(table, alternative=mode_keyword)
			print(oddsr, p)
	


if __name__ == "__main__":
	try:
		input_file = sys.argv[1]
		analysis_mode = sys.argv[2]
	except Exception:
		print('Program usage: text.py <a one or multiline file that has on one line the first 4 fields as integers values space separated and from the fifth on the name of the feature, like: 1 2 3 4 proline rich, where hte first value is left top value (a) second is top right (b) third is bottom left (c) and right bottom the fourth (d), is user responasability to set correctly up take a look at pipeline    fisher_exact_test.nf    for more info> <second mandatory field resonsible for type of analisis to be done with fischer, as extreme as "two-sided" default of scipy module that is the probability a random table would have a probability equal to or less than the probability of the input table, "greater" is the probability that a random table has x >= a, which means with left top/ first number higher than what is given as input, "less" the one-sided p-value is the probability that a random table has x <= a >', file=sys.stderr)
		raise SystemExit
	else:
		fisherer(input_file, analysis_mode)
