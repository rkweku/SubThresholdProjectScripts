# Report the sensitivity and specificity of this constrained method for the
# loci identified by each cell type, as well as the pooled gene set from
# all provided cell types. Additionally, identify the average sensitivity and
# specificity using 10,000 randomly sampled populations where the population
# size is the number of loci identified at the corresponding band level of
# the pooled constrained loci. Finally, identify the sensitivity and 
# specificity of the negative constrained dataset, where we assess the loci
# that are not in regions of OCR or whose gene promoters do not physically
# contact an open SNP

import csv
import math
import numpy
import random
import sys

###############################################################################

# A constrained count file
# Example dataset: "Figure 4-Source Data 1.csv"
filteredFilename = sys.argv[1]

# An unconstrained count file
# Example dataset: "Figure 4-Source Data 1.csv"
unfilteredFilename = sys.argv[2]

# Name of the output file (in csv format)
outputFilename = sys.argv[3]

###############################################################################

def randomSampling(filteredRatesDict, samplingRatesDict):
    # Initialize the population dictionary for each band consisting of 1s
    # and 0s
    samplePopulationDict = {}
    # Initialize two variables to track the sensitivities and specificities
    # of the bootstrapped populations
    trackingSensitivities = {}
    trackingSpecificities = {}
    for bandLevel in samplingRatesDict:
        samplePopulationDict[bandLevel] = []
        trackingSensitivities[bandLevel] = []
        trackingSpecificities[bandLevel] = []
        for i in range(samplingRatesDict[bandLevel][0]):
            samplePopulationDict[bandLevel].append(1)
        for i in range(samplingRatesDict[bandLevel][1]):
            samplePopulationDict[bandLevel].append(0)

    # Create 10,000 bootstrapped populations to identify the average
    # sensitivity and specificity for each sampling rates dictionary using
    # a population size of the corresponding filteredRatesDict band level
    for i in range(10000):
        subsamplesDict = {}
        prevBand = ""
        for bandLevel in samplePopulationDict:
            # Get the average population size from this band level using the
            # prediction data
            populationSize = math.ceil(sum(
                filteredRatesDict["All Cells"][bandLevel]))
            samplingListCopy = samplePopulationDict[bandLevel].copy()

            if(prevBand):
                subsamplesDict[bandLevel] = subsamplesDict[prevBand].copy()
                for entry in subsamplesDict[prevBand]:
                    samplingListCopy.remove(entry)
                populationSize = populationSize - len(subsamplesDict[prevBand])

            else:
                subsamplesDict[bandLevel] = []

            subsamplesDict[bandLevel].extend(random.sample(samplingListCopy,
                populationSize))
            prevBand = bandLevel

            positiveCount = samplingRatesDict[bandLevel][0]
            negativeCount = samplingRatesDict[bandLevel][1]

            tp = sum(subsamplesDict[bandLevel])
            fp = len(subsamplesDict[bandLevel]) - tp

            if(positiveCount > 0):
                sensitivity = tp / positiveCount
            else:
                sensitivity = 0

            if(negativeCount > 0):
                specificity = 1 - (fp / negativeCount)
            else:
                specificity = 0

            trackingSensitivities[bandLevel].append(sensitivity)
            trackingSpecificities[bandLevel].append(specificity)

    ratesDict = {}
    for bandLevel in trackingSensitivities:
        sensitivity = numpy.mean(trackingSensitivities[bandLevel])
        specificity = numpy.mean(trackingSpecificities[bandLevel])
        ratesDict[bandLevel] = (sensitivity, specificity)

    return(ratesDict)

def getRatesFromPPVFile(filteredFilename, unfilteredFilename):
    """Get the sensitivity and specifificity from the PPV files

    """

    filteredRatesDict = {}
    unfilteredRatesDict = {}

    # Read through the PPV file for the constraned file
    with open(filteredFilename) as f:
        reader = csv.reader(f, delimiter = ",")
        header = next(reader)[1:]

        for row in reader:
            bandLevel = row[0]
            
            for i in range(len(row[1:])):
                cellType = header[i]

                if(cellType not in filteredRatesDict):
                    filteredRatesDict[cellType] = {}

                # Find the number of positive and negative predictions for
                # the SNP clusters in this band level and cell type
                ppv = row[1+i]
                numPositive = int(ppv.split("/")[0])
                numNegative = int(ppv.split("/")[1]) - numPositive

                filteredRatesDict[cellType][bandLevel] = (numPositive,
                    numNegative)

    # Read through the PPV file for the unconstrained file
    with open(unfilteredFilename) as f:
        # Skip the header here
        next(f)
        reader = csv.reader(f, delimiter = ",")

        # For each band level, find the number of true positives and true
        # negatives to be used for tp, fp, fn, and tn values later
        for row in reader:
            countsList = []
            bandLevel = row[0]
            numPositive = int(row[1].split("/")[0])
            N = int(row[1].split("/")[1])
            numNegative = N - numPositive

            unfilteredRatesDict[bandLevel] = (numPositive, numNegative)

    # Using the constrained and unconstrained values, identify the set of
    # SNP clusters that would be part of the NEGATIVE set. i.e. remove the
    # SPNs that survive the filters in either cell type to plot those values
    # Additionally, output this data to its own file for separate analysis, if
    # necessary
    negativeRatesDict = {}
    for cellType in filteredRatesDict:
        negativeRatesDict[cellType] = {}

        for bandLevel in filteredRatesDict[cellType]:
            unfilteredPos = unfilteredRatesDict[bandLevel][0]
            unfilteredNeg = unfilteredRatesDict[bandLevel][1]
            filteredPos = filteredRatesDict[cellType][bandLevel][0]
            filteredNeg = filteredRatesDict[cellType][bandLevel][1]
            negativeRatesDict[cellType][bandLevel] = (
                unfilteredPos - filteredPos, unfilteredNeg - filteredNeg) 

    # Using the true values of the unconstrained data, find the sensitivity
    # and specificity for each band and cell type
    ratesDict = {}
    for cellType in filteredRatesDict:
        if(cellType not in ratesDict):
            ratesDict[cellType] = {}
            ratesDict["Negative %s" % cellType] = {}
        for bandLevel, values in filteredRatesDict[cellType].items():
            totalPositives = unfilteredRatesDict[bandLevel][0]
            totalNegatives = unfilteredRatesDict[bandLevel][1]

            tp = values[0]
            fp = values[1]
            fn = totalPositives - tp
            tn = totalNegatives - fp

            if(tp > 0 or fn > 0):
                sensitivity = tp / (tp + fn)
            else:
                sensitivity = 0

            if(tn > 0 or fp > 0):
                specificity = tn / (tn + fp)
            else:
                specificity = 0

            ratesDict[cellType][bandLevel] = (sensitivity, specificity)

    # This code generates the sampled data to get the sensitivity and 
    # specificity of randomly sampled data in each suggestive zone
    ratesDict["Naive Sampling"] = randomSampling(filteredRatesDict,
        unfilteredRatesDict)

    # Get the sensitivity and specificity of the randomly sampled population
    # of the negative set of data for each cell type
    for cellType in negativeRatesDict:
        ratesDict["Negative %s" % cellType] = randomSampling(filteredRatesDict,
            negativeRatesDict[cellType])

    return(ratesDict)

def main():
    ratesDict = getRatesFromPPVFile(filteredFilename, unfilteredFilename)

    with open(outputFilename, "w") as g_out:
        writer = csv.writer(g_out, delimiter = ",")

        writer.writerow(["Cell Type", "Band Level", "Sensitivity",
            "Specificity"])
        for cellType in ratesDict:
            for bandLevel in ratesDict[cellType]:
                sensitivity = ratesDict[cellType][bandLevel][0]
                specificity = ratesDict[cellType][bandLevel][1]
                writer.writerow([cellType, bandLevel, sensitivity,
                    specificity])

if __name__ == "__main__":
    main()

