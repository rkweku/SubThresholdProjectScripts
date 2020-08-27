# Create empirical distribution plot depicting the distribution of PPV of
# a randomly sampled set of loci identified from unconstrained dataset and
# where the PPV of each constrained dataset would fall within that distribution

import csv
import math
import random
import sys

import matplotlib.pyplot as plt 
import matplotlib.patches as mpatches
import numpy as np
import scipy.stats as st
import seaborn as sns

from statsmodels.distributions.empirical_distribution import ECDF
from sklearn.metrics import auc

###############################################################################

# File with the counts of TP vs FP in the constrained dataset
# Example file would be "Figure 2-Source Data 1"
filteredFilename = sys.argv[1]

# File with the counts of the TP vs FP in the unconstrained datset (assuming 
# all subthreshold signals would reach GWS by the future GWAS, only those in
# near LD [r2>0.8] with a signal that is GWS in a future GWAS are TP)
# Example file would be "Figure 2-Source Data 2"
unfilteredFilename = sys.argv[2]

# Prefix for the output file. 4 files per run, hence why we use a prefix here
# This can be anything. Output file will be "outputPrefix_5x10^-N_PPV.png"
# where N is either -7, -6. -5, or -4
outputPrefix = sys.argv[3]

###############################################################################

def getRatesFromPPVFile(filteredFilename, unfilteredFilename):
    """Get the TPR and FPR from the values within the filtered file and the
    unfiltered file

    """

    filteredRatesDict = {}

    with open(filteredFilename) as f:
        reader = csv.reader(f, delimiter = ",")
        header = next(reader)[1:]

        for row in reader:
            bandLevel = row[0]
            
            
            for i in range(len(row[1:])):
                cellType = header[i]
                if(cellType not in filteredRatesDict):
                    filteredRatesDict[cellType] = {}

                ppv = row[1+i]
                numPositive = int(ppv.split("/")[0])
                numNegative = int(ppv.split("/")[1]) - numPositive

                filteredRatesDict[cellType][bandLevel] = (numPositive, numNegative)

    ratesDict = {}
    for cellType in filteredRatesDict:
        ratesDict[cellType] = ([0],[0])
        numPositive = filteredRatesDict[cellType]["5x10^-4"][0]
        numNegative = filteredRatesDict[cellType]["5x10^-4"][1]

        for bandLevel in filteredRatesDict[cellType]:
            truePositive = filteredRatesDict[cellType][bandLevel][0]
            falsePositive = filteredRatesDict[cellType][bandLevel][1]

            if(numPositive):
                tpr = float(truePositive / numPositive)
            else:
                tpr = 0
            if(numNegative):
                fpr = float(falsePositive / numNegative)
            else:
                fpr = 0
            
            ratesDict[cellType][0].append(tpr)
            ratesDict[cellType][1].append(fpr)

    countsDict = {}

    with open(unfilteredFilename) as f:
        # Skip the header here
        next(f)
        reader = csv.reader(f, delimiter = ",")

        for row in reader:
            countsList = []
            bandLevel = row[0]
            numPositive = int(row[1].split("/")[0])
            N = int(row[1].split("/")[1])
            numNegative = N - numPositive

            countsDict[bandLevel] = (numPositive, numNegative)

    samplesDict = {}
    for bandLevel in countsDict:
        samplesDict[bandLevel] = []
        for i in range(countsDict[bandLevel][0]):
            samplesDict[bandLevel].append(1)
        for i in range(countsDict[bandLevel][1]):
            samplesDict[bandLevel].append(0)

    trackingPositives = {}
    trackingNegatives = {}
    trackingPPVs = {"Naive Sampling": {}}

    for bandLevel in samplesDict:
        trackingPositives[bandLevel] = []
        trackingNegatives[bandLevel] = []
        trackingPPVs["Naive Sampling"][bandLevel] = []

    # Code to get a bootstrapped distribution
    for i in range(10000):
        for bandLevel in samplesDict:
            # Bootstrap a population from this band level
            bootstrappedPopulation = np.random.choice(
                samplesDict[bandLevel], size = 10000, replace = True)
            randomSamplingList = np.random.choice(list(bootstrappedPopulation),
                1000, replace=False)
            tp = sum(randomSamplingList)
            fp = len(randomSamplingList) - tp
            trackingPPVs["Naive Sampling"][bandLevel].append(tp / (tp + fp))

    # This code generates the sampled data to create the distribution
    for i in range(10000):
        subsamplesDict = {}
        prevBand = ""
        for bandLevel in samplesDict:
            # Get the average population size from this band level using the
            # prediction data
            populationSize = math.ceil(sum(filteredRatesDict["All Cells"][bandLevel]))
            samplingListCopy = samplesDict[bandLevel].copy()

            if(prevBand):
                subsamplesDict[bandLevel] = subsamplesDict[prevBand].copy()
                for entry in subsamplesDict[prevBand]:
                    samplingListCopy.remove(entry)
                populationSize = populationSize - len(subsamplesDict[prevBand])

            else:
                subsamplesDict[bandLevel] = []

            subsamplesDict[bandLevel].extend(random.sample(samplingListCopy, populationSize)) 
            prevBand = bandLevel

        positiveCount = sum(subsamplesDict["5x10^-4"])
        negativeCount = len(subsamplesDict["5x10^-4"]) - positiveCount

        for bandLevel, samplesList in subsamplesDict.items():
            tp = sum(samplesList)
            fp = len(samplesList) - tp

            if(positiveCount):
                tpr = tp / positiveCount
            else:
                tpr = 0

            if(negativeCount):
                fpr = fp / negativeCount
            else:
                fpr = 0

            trackingPositives[bandLevel].append(tpr)
            trackingNegatives[bandLevel].append(fpr)

    ratesDict["Naive Sampling"] = ([0], [0])
    for bandLevel in trackingPositives:
        trackingPositivesList = trackingPositives[bandLevel]
        trackingNegativesList = trackingNegatives[bandLevel]
        ratesDict["Naive Sampling"][0].append(np.mean(trackingPositivesList))
        ratesDict["Naive Sampling"][1].append(np.mean(trackingNegativesList))

    for cellType in filteredRatesDict:
        trackingPPVs[cellType] = {}
        for bandLevel in filteredRatesDict[cellType]:
            tp = filteredRatesDict[cellType][bandLevel][0]
            fp = filteredRatesDict[cellType][bandLevel][1]
            if(tp+fp > 0):
                trackingPPVs[cellType][bandLevel] = (tp / (tp + fp))
            else:
                trackingPPVs[cellType][bandLevel] = 0

    return(ratesDict, trackingPPVs)

def main():
    filteredRatesDict, ppvs = getRatesFromPPVFile(filteredFilename, unfilteredFilename)

    # Code to plot probability distribution
    for bandLevel in ppvs["Naive Sampling"]:
        data = ppvs["Naive Sampling"][bandLevel]
        print (np.percentile(data, [.01, 50, 99.99]))
#        print(np.mean(data))
#        CIs = st.t.interval(0.95, len(data) - 1, loc=np.mean(data),
#            scale=st.sem(data))
#        print(bandLevel, CIs)        
        
        plt.figure(figsize = (11,5))
        p = sns.distplot(ppvs["Naive Sampling"][bandLevel])
        xMin = min(ppvs["Naive Sampling"][bandLevel])
        xMax = 0
        for cellType in ppvs:
            if(cellType == "Naive Sampling"):
                continue
            if(ppvs[cellType][bandLevel] > xMax):
                xMax = ppvs[cellType][bandLevel]

        #fig, ax = plt.subplots()
        p.set_title("Empirical Distribution of PPVs from Naively Sampled "\
            "Data %s" % bandLevel, fontsize = 16)
        p.set_xlabel("Positive Predictive Value", fontsize = 14)
        p.set_xlim(xMin-.1, xMax+.1)
        p.tick_params(labelsize = 12)
#        p.set_xticklabels(p.set_xticks(), size=24)
#        p.set_yticklabels(p.get_yticks(), size=24)
        p.set_ylabel("Frequency", fontsize = 14)
#        p.set(title = "Empirical Distribution of PPVs from Naively Sampled "\
#            "Data %s" % bandLevel, xlabel = "Positive Predictive Value",
#            ylabel = "Frequency", xlim = (xMin - .1, xMax + .1))
        for cellType in ppvs:
            if(cellType == "Naive Sampling" or cellType == "All Cells"):
                continue
            if(cellType == "Adipose"):
                plt.axvline(x=ppvs[cellType][bandLevel], ymin = 0,
                    ymax = .75, color = "red", label = "Adipose", 
                    linewidth = 3, linestyle = "dashed")
            elif(cellType == "HypothalamicNeurons"):
                plt.axvline(x=ppvs[cellType][bandLevel], ymin = 0,
                    ymax = .75, color = "orange", label = "Hypothalamic Neurons", 
                    linewidth = 3, linestyle = "dotted")
        plt.rc("legend", fontsize = 12)
        plt.legend()
        fig = p.get_figure()

        fig.savefig("%s_%s-PPV.png" % (outputPrefix, bandLevel))
        fig.clf()

if __name__ == "__main__":
    main()

