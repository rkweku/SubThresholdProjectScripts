# Use the summary stats and the external binary search to annotate the 
# positions of SNPs from a summary stats file

import csv
import subprocess
import sys

###############################################################################
# Summary stats filename
summaryStatsFilename = sys.argv[1]

# The SNP ID index in the summary stats (Index starting at 0)
snpIndex = int(sys.argv[2])

pvalueIndex = int(sys.argv[3])

outputFilename = sys.argv[4]

# The name of the variation file
variationFilename = "/mnt/isilon/sfgi/hammondrk/hg19Data/variation/hg19.snp151.sortedByName.txt"

binarySearchPath = "/home/hammondrk/pts-line-bisect/pts_lbsearch"

###############################################################################

def main():
    # Open the output file to write the data to
    with open(outputFilename, "w") as g_out:

        with open(summaryStatsFilename, "r") as f_in:
            reader = csv.reader(f_in, delimiter = "\t")
            next(reader)

            for row in reader:
                snpID = row[snpIndex].split(":")[0]
                pvalue = row[pvalueIndex]

                command = "%s -p %s %s |  head -1 | awk -v OFS='\t' '{ "\
                    "print $1,$2,$3; }'" % (binarySearchPath,
                    variationFilename, snpID)

                results = subprocess.run(command, shell = True,
                    stdout=subprocess.PIPE)

                toWrite = results.stdout.decode("UTF-8").strip()
                toWrite += "\t%s" % pvalue

                g_out.write("%s\n" % toWrite)

if __name__ == "__main__":
    main()
