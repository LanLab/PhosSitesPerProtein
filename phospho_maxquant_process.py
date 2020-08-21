import re
import sys
import argparse

############################################################## FUNCTIONS


def checkAndIncr(prob, countVar, styProbTh):

    if float(prob) >= styProbTh:
        countVar = countVar + 1

    return (countVar)

def getCountInPeptide(peptide, styProbTh):
    s = doTheCount('S', peptide, styProbTh)
    t = doTheCount('T', peptide, styProbTh)
    y = doTheCount('Y', peptide, styProbTh)

    return (s, t, y)


def doTheCount(residue, seq, styProbTh):
    count = 0
    val = re.findall(residue + '\(([0-9\.]+)\)', seq)
    if len(val)>0:
        if type(val) is list:
            for v in val:
                    # print(seq + '\t' + v)
                count = checkAndIncr(v, count, styProbTh)
        else:
            count = checkAndIncr(val.group(1), count, styProbTh)
    return (count)


def residual_count(dict_seq, aa, styProbTh):
    y_count = 0
    for seq in dict_seq:
        val = re.findall(aa+'\(([0-9\.]+)\)', seq)
        if len(val)>0:
            if type(val) is list:
                for v in val:
                    # print(seq + '\t' + v)
                    y_count = checkAndIncr(v, y_count, styProbTh)
            else:
                y_count = checkAndIncr(val.group(1), y_count, styProbTh)
    return (y_count)

############################################################## TOP LEVEL


def countAndPrint(fn_maxquant, protColumn, locProbColumn, locProbTh, styProbColumn, styProbTh, fn_out):

	dict_prot = dict()
	dict_seq = dict()

	with open(fn_maxquant, 'r+') as fh:
		isHeader = True
		for line in fh:
			if isHeader == True:
				isHeader = False
				continue

			line = re.sub('[\n\r]$', '', line)
			arr = line.split("\t")

			# print(arr[protColumn]);
			if re.match('^CON', arr[protColumn]) or re.match('^[\s]*$', arr[protColumn]):
				continue


			if float(arr[locProbColumn]) >= locProbTh:

				if (arr[protColumn], arr[styProbColumn]) not in dict_prot:
					(s, t, y) = getCountInPeptide(arr[styProbColumn], styProbTh)
					dict_prot[(arr[protColumn], arr[styProbColumn])] = (s, t, y) # 'S\t' + str(s) + '\t' + 'T' + '\t' + str(t) + '\t' + str(y)
				else:
					# print("Already added " + arr[protColumn] + " " + arr[styProbColumn])
					pass

				if arr[styProbColumn] in dict_seq:
					dict_seq[arr[styProbColumn]] = dict_seq[arr[styProbColumn]] + 1
				else:
					dict_seq[arr[styProbColumn]] = 1


	        #print(line)
	        #break


	##### Printing the outputs:
	s_count = residual_count(dict_seq, 'S', styProbTh)
	t_count = residual_count(dict_seq, 'T', styProbTh)
	y_count = residual_count(dict_seq, 'Y', styProbTh)

	print("S: " + str(s_count))
	print("T: " + str(t_count))
	print("Y: " + str(y_count))

	# print(dict_seq);

	print("======")
	print("Writing to file " + fn_out);
	fhandle = open(fn_out, 'w+')
	# print(dict_prot);
	for (prot, peptide) in dict_prot:
		(ser, thr, tyr) = dict_prot[(prot, peptide)]
		fhandle.write(prot + "\t" + peptide + "\t" +  'S\t' + str(ser) + '\t' + 'T' + '\t' + str(thr) + '\t' + 'Y' + '\t' + str(tyr) + '\n')

	fhandle.close()

############################################################## MAIN

def main():
	# usage = "python " + sys.argv[0] + " <styFile.txt> <sty_colnum>"

	# if len(sys.argv) != 2:
	# 	sys.exit("Error: incorrect number of inputs\n" + usage + '\n\n')

	# countAndPrint(sys.argv[1])


	parser = argparse.ArgumentParser(description='Extract S, T and Y site counts (which are above thresholds) from MaxQuant data.')
	parser.add_argument('--maxQuantInFile', nargs=1, required=True, help='An output file from MaxQuant, as input to this script.')

	parser.add_argument('--protColumn', nargs=1, default=0, type=int, help="Column containing the protein identifier. Default is 0.")

	parser.add_argument('--locProbColumn', nargs=1, default=5, type=int, help="Column containing the localization probability. Default is 5.")
	parser.add_argument('--locProbTh', nargs=1, default=0.75, type=float, help="Threshold to filter localization probability (between 0 and 1). Default is 0.75.")

	parser.add_argument('--styProbColumn', nargs=1, default=29, type=int, help="Column containing the phospho-(STY)-probabilities. Default is 29.")
	parser.add_argument('--styProbTh', nargs=1, default=0.75, type=float, help="Threshold to filter phospho-(STY)-probabilities (between 0 and 1). Default is 0.75.")

	parser.add_argument('--outfile', nargs=1, default="Output.txt", help="Output file to save the results to. Default is Output.txt.")
	# parser.add_argument('--maxQuantInFile', nargs=1, required=True, help='An output file from MaxQuant, as input to this script.')

	args = parser.parse_args()

	if isinstance(args.protColumn, list):
		print("is list");
		args.protColumn = args.protColumn[0];

	if isinstance(args.locProbColumn, list):
		args.locProbColumn = args.locProbColumn[0];

	if isinstance(args.locProbTh, list):
		args.locProbTh = args.locProbTh[0];

	if isinstance(args.styProbColumn, list):
		args.styProbColumn = args.styProbColumn[0];

	if isinstance(args.styProbTh, list):
		args.styProbTh = args.styProbTh[0];

	if isinstance(args.outfile, list):
		args.outfile = args.outfile[0];

	if not (args.locProbTh >= 0 and args.locProbTh <= 1):
		#sys.stderr.write()
		raise argparse.ArgumentTypeError("Error: Argument to --locProbTh should be between 0 and 1.")

	if not (args.styProbTh >= 0 and args.styProbTh <= 1):
		#sys.stderr.write()
		raise argparse.ArgumentTypeError("Error: Argument to --styProbTh should be between 0 and 1.")

	print(isinstance(args.protColumn, list))

	countAndPrint(args.maxQuantInFile[0], args.protColumn, args.locProbColumn, args.locProbTh, args.styProbColumn, args.styProbTh, args.outfile)

if __name__ == '__main__':
	main()





# LIST_SEQ = list()


"""
if len(sys.argv) != 3:
    sys.stderr.write("Incorrect number of inputs " + '\n' + usage + '\n\n')
    sys.exit()
"""
