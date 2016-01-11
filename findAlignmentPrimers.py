from __future__ import division
import argparse
import numpy
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pygal 
from pygal.style import Style

def main():
    options = parseArguments()
    aln = AlignIO.read(open(options.alignment),"fasta")
    primers = findPrimers(aln, options)
    # primers is a list of [(start,stop),consensus,% conserved over primer] #conservation can be used for ranking?
    primerRegionSequences = getPrimerRegions(primers, options, aln)
    # Idea would now be that each primer region could be used with primer3 to design a primer to that region.
    
    
    
    
def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-alignment', help='Multiple sequence alignment in FASTA format', required=True)
    parser.add_argument('-maxprimer', help='Maximum primer length in nucleotides. Default=30', default=30, type=int)
    parser.add_argument('-similarity', help='Minimum percent similarity between sequences in the primer region. Default=95', default=95, type=int)
    parser.add_argument('-degeneracy', help='Percentage of bases in an alignment column which can be different before a degeneracy is called. Default=5', default=5, type=int)
    parser.add_argument('-plot', help='Produces an SVG file of a plot of the conservation across the alignment. Filename as argument, otherwise default is plot.svg',default='plot.svg')
    parser.add_argument('-maxgaps', help='Maximum number of alignment gaps that should be included in the output primer region. Alternatively, edit your alignment to remove single insertions. Default=3',default=4,type=int)
    args = parser.parse_args()
    return args

"""
getPrimerRegions(primers)
Loops through the primers returned by findPrimers() and calculates the start and stop
positions of each primer region as many of the primers will overlap by 1bp.
The regions are returned as consensus (read: degenerate) seqRecord objects where the seq.id shows
the coordinates in the main alignment.
Arguments:
primers     List of primers [(start,stop), consensus, %conserved for whole primer]
regionPositions List of tuples (start,stop) for each primer region in the main alignment object
alignment   AlignIO object of the main alignment being used.
Returns:
consensusSequences  List of SeqRecord objects of the consensus sequences at each primer region.
"""
def getPrimerRegions(primers, options, alignment):
    consensusSequences = []
    
    # Pull out all the (start,end) tuples from primers
    starts = []
    for primer in primers:
        starts.append(primer[0][0])
    
    coordinates = []
    # Go through the starts. Check if the next number is i+1. If not, then output slice
    sliceStart = 0
    primerLength = options.maxprimer
    for i in range(0,len(starts)):
        start = starts[i]
        try:
            nextStart = starts[i+1]
        except:
            break
        if (nextStart != start+1):
            coordinates.append((sliceStart,start+primerLength))
            sliceStart = nextStart
    
    # Create the slices
    for pair in coordinates:
        slicedAlignment = alignment[:,pair[0]:pair[1]]
        consensus = Seq(calculateConsensusSequence(slicedAlignment,options.degeneracy))
        seq = SeqRecord(consensus)
        seq.id = str(pair[0]) + "-" + str(pair[1])
        consensusSequences.append(seq)
    
    return consensusSequences

"""
findPrimers()
Arguments:
alignment   Multiple sequence alignment object (BioPython)
options     The output of the parseArguments() function

Returns:
primers List of degenerate primer sequences and their locations (5'->3') on the given alignment
"""
def findPrimers(alignment, options):
    slidingWindow = options.maxprimer #Maximum primer size
    similarityThreshold = options.similarity
    degeneracyThreshold = options.degeneracy
    primers = []
    
    maxBaseFrequencyPercentages = [] #Holder for column percentages (TEST PLOT).
    highSimilaritySlices = [] #Holder for potential primers
    # Go through the alignment L-R, slicing out the sliding window and assessing
    # it for the base frequencies in each column in the slice.
    alignmentLength = alignment.get_alignment_length()
    for i in range(0,alignmentLength-slidingWindow):
        # Slice the alignment
        sliced = alignment[:, i:i+slidingWindow]
        aCounts = [] 
        gCounts = []
        cCounts = []
        tCounts = []
        nCounts = []
        dashCounts = []
        sliceMaxBaseFrequencyPercentages = [] # %age of the most freq base in a column
        #Loop through the slice and calculate the base frequencies for each column
        for j in range(0,sliced.get_alignment_length()):
            column = sliced[:,j].upper() #upper !imp
            a = column.count('A')/len(column)*100
            g = column.count('G')/len(column)*100
            c = column.count('C')/len(column)*100
            t = column.count('T')/len(column)*100
            n = column.count('N')/len(column)*100
            dash = column.count('-')/len(column)*100
            aCounts.append(a)
            gCounts.append(g)
            cCounts.append(c)
            tCounts.append(t)
            nCounts.append(n)
            dashCounts.append(dash)
            # Grab the percentage of the most abundant base in the column and append
            # to the slice holding list.
            maxBaseFrequencyPercentage = max(a,g,c,t,n,dash)
            sliceMaxBaseFrequencyPercentages.append(maxBaseFrequencyPercentage)
            #print i,a,c,g,t,n,dash,maxBaseFrequencyPercentage
        # Calculate the mean maxBaseFrequencyPercentage for the slice
        sliceMeanMaxBaseFrequencyPercentage = numpy.mean(sliceMaxBaseFrequencyPercentages)
        # add to the holding list for plotting. Tuple for plotting
        maxBaseFrequencyPercentages.append((i,sliceMeanMaxBaseFrequencyPercentage)) 
        # If this slice mean percentage is above the similarity threshold put it in the holder
        if (sliceMeanMaxBaseFrequencyPercentage >= similarityThreshold):
            highSimilaritySlices.append(((i,i+slidingWindow),sliceMeanMaxBaseFrequencyPercentage))
    
    # Go through the highSimilaritySlices and calculate the consensus sequences for each
    for i in highSimilaritySlices:
        sliceMeanMaxBaseFrequencyPercentage = i[1]
        alignmentSlice = alignment[:, i[0][0]:i[0][1]]
        alignmentSummary = AlignInfo.SummaryInfo(alignmentSlice)
        consensus = calculateConsensusSequence(alignmentSlice, degeneracyThreshold)
        # Calculate the number of gaps ('-') in the consensus sequence. If more than the
        # requested number don't report this sequence
        gaps = consensus.count('-')
        if (gaps <= options.maxgaps):
            #primerData = [(alignment start, alignment end),consensus,percentageConserved]
            primerData = [(i[0][0],i[0][1]),consensus,"{0:.2f}".format(i[1])]
            primers.append(primerData)

    # Print a test plot of the MaxBaseFrequencyPercentage at each position in the alignment
    # if requested
    if (options.plot):
        producePercentagesPlot(maxBaseFrequencyPercentages, options.plot, slidingWindow, similarityThreshold)
        
    #return the primers list
    return primers

"""
calculateConsensusSequence()
Arguments:
aln   Multiple sequence alignment object (BioPython)
threshold  options.degeneracy. Int. The percentage of two or more bases that is required in an alignment column to trigger the calling of a degeneracy. 

Returns:
consensus   String. The consensus sequence of the alignment provided.
"""
def calculateConsensusSequence(aln,threshold):
    # Degeneracy lookup. The single bases are required because the later collapse of the
    # bases required in the consensus calculation using the threshold might result in a
    # clear base call.
    degeneracies = {
        'AG':'R', 'CT':'Y', 'AC':'M', 'CG':'S', 'AT':'W', 'GT':'K', 'ACG':'V',
        'AGT':'D', 'ACT':'H', 'CGT':'B', 'ACGT':'N', 'A':'A', 'C':'C', 'T':'T','G':'G'
    } 
    consensus = ""
    for i in range(0,aln.get_alignment_length()):
        column = aln[:,i].upper() #upper !imp
        a = column.count('A')/len(column)*100
        g = column.count('G')/len(column)*100
        c = column.count('C')/len(column)*100
        t = column.count('T')/len(column)*100
        n = column.count('N')/len(column)*100
        dash = column.count('-')/len(column)*100
        # If the column is clearly above the % threshold, call the base
        if (a >= (100-threshold)): consensus += 'A'
        elif (g >= (100-threshold)): consensus += 'G'
        elif (c >= (100-threshold)): consensus += 'C'
        elif (t >= (100-threshold)): consensus += 'T'
        elif (n >= (100-threshold)): consensus += 'N'
        elif (dash >= (100-threshold)): consensus += '-'
        else: #Otherwise, this column requires a degeneracy
            column = column.translate(None,'-') #remove dashes
            # Recalculate the percentages without the dashes
            a = column.count('A')/len(column)*100
            g = column.count('G')/len(column)*100
            c = column.count('C')/len(column)*100
            t = column.count('T')/len(column)*100
            # If any of the %bases are above the threshold, add to consensusBases
            # for calculation of the degenerate base call
            degenerateBases = [] 
            if (a >= threshold): degenerateBases.append('A')
            if (g >= threshold): degenerateBases.append('G')
            if (c >= threshold): degenerateBases.append('C')
            if (t >= threshold): degenerateBases.append('T')
            degenerateBases = ''.join(sorted(set(degenerateBases))) #collapse string
            if degenerateBases in degeneracies:
                degeneracy = degeneracies[degenerateBases]
            else:
                degeneracy = 'X' #undetermined
            consensus += degeneracy
    return consensus

"""
producePercentagesPlot()
Arguments:
percentages List of tuples. The maximumBaseFrequencyPercentages for each slice of the alignment using a sliding window of the maximum primer length.
outfile String. Name of the output svg file.
slidingWindow   options.maxprimer. Int. The maximum size of the primer region required.
similarity  options.similarity. Int. The similarity threshold chosen by the user that is used to select conserved regions. 
Output:
file    SVG image (open in browser) of the conservation along the length of the alignment. 
"""
def producePercentagesPlot(percentages,outfile,slidingWindow,similarity):
    similarityLine = [] #To put a line for the similarity threshold on the plot
    for i in range(0,len(percentages)):
        similarityLine.append((i,similarity)) #i for XY
    custom_style = getPygalCustomStyle()
    chart = pygal.XY(
        style=custom_style,
        show_dots=False,
        fill=False,
        legend_at_bottom=True,
        legend_font_size=12)
    chart.title = "Conservation across the alignment with a sliding window of " + str(slidingWindow) + "nt"
    chart.x_title = "Alignment position"
    chart.y_title = "Percent conservation"
    chart.add("Conservation", percentages)
    chart.add("Similarity threshold for primer selection",similarityLine)
    chart.render_to_file(outfile)

"""
getPygalCustomStyle()
Returns a custom colour style for the PyGal plotting library. 
"""    
def getPygalCustomStyle():
    style = Style(
        background = '#fff',
        plot_background = '#fff',
        foreground='#000',
        foreground_light = '#545454',
        opacity = '0.9',
        opacity_hover = '1',
        transition = '400ms ease-in',
        colors = ( '#88bbdd','#ff1100','#66bb44','#ffee00','#0000ff', '#999', '#000')
    )
    return style

if __name__ == '__main__':
    main()