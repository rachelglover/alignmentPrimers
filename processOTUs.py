# For the diatom alignment taxonomy assessments
# To be run on the QIIME server so no BioPython!

import argparse
import sys
from collections import defaultdict

def main():
    options = parseArguments()
    # The OTUs and taxonomy files are indexed by the OTU number. I'd rather have
    # an output of each actual sequence in the alignment slice and it's taxonomic
    # assignment
    sequences = defaultdict(lambda: defaultdict(dict))

    # Grab the names of the sequences in the alignment file
    alignment_sequences = []
    for line in open(options.alignment, "rU"):
        if line.startswith(">"):
            line = line.rstrip()
            sample = line.strip(">")
            alignment_sequences.append(sample)

    # Load in all the correct taxonomies to the sequences dict
    for line in open(options.alltax, "rU"):
        line = line.rstrip()
        if (line.startswith("Strain")):
            pass
        else:
            linelist = line.split('\t')
            seqname = linelist[0]
            if seqname in alignment_sequences:
                taxonomy = linelist[1]
                taxonomylist = taxonomy.split(';')
                sequences[seqname]["correct_taxonomy"]["full"] = taxonomy
                sequences[seqname]["correct_taxonomy"]["class"] = taxonomylist[0]
                sequences[seqname]["correct_taxonomy"]["family"] = taxonomylist[1]
                sequences[seqname]["correct_taxonomy"]["genus"] = taxonomylist[2]
                sequences[seqname]["correct_taxonomy"]["species"] = taxonomylist[3]
                sequences[seqname]["correct_taxonomy"]["strain"] = taxonomylist[4]
                # Note: "strain" for the DTM_composite is the seqname

    # Now go through the OTU taxonomy and create a lookup
    otu_taxonomies = {}
    for line in open(options.tax, "rU"):
        line = line.rstrip()
        linelist = line.split('\t')
        otu = int(linelist[0])
        taxonomy = linelist[1]
        otu_taxonomies[otu] = taxonomy

    # Now go through the otus and go through each of the samples and assign the actual taxonomy
    for line in open(options.otus, "rU"):
        line = line.rstrip()
        linelist = line.split('\t')
        otu = int(linelist[0])
        linelist.pop(0) #linelist now only contains sequence ids
        # Grab the taxonomy for this OTU
        otu_tax = otu_taxonomies[otu]
        if (otu_tax.startswith("No blast hit")):
            otu_tax = "NULL;NULL;NULL;NULL;NULL;"
        otu_taxlist = otu_tax.split(';')
        # Assign this OTU taxonomy to all sequences associated with this OTU
        for seqname in linelist:
            sequences[seqname]["actual_taxonomy"]["full"] = otu_tax
            sequences[seqname]["actual_taxonomy"]["class"] = otu_taxlist[0]
            sequences[seqname]["actual_taxonomy"]["family"] = otu_taxlist[1]
            sequences[seqname]["actual_taxonomy"]["genus"] = otu_taxlist[2]
            sequences[seqname]["actual_taxonomy"]["species"] = otu_taxlist[3]
            sequences[seqname]["actual_taxonomy"]["strain"] = otu_taxlist[4]

    # Not all DTM taxonomy sequences will have been in the original alignment
    for seqname in sequences:
        try:
            actual = sequences[seqname]["actual_taxonomy"]["full"]
        except:
            sequences[seqname]["actual_taxonomy"]["full"] = "not-in-slice"
            sequences[seqname]["actual_taxonomy"]["class"] = "not-in-slice"
            sequences[seqname]["actual_taxonomy"]["family"] = "not-in-slice"
            sequences[seqname]["actual_taxonomy"]["genus"] = "not-in-slice"
            sequences[seqname]["actual_taxonomy"]["species"] = "not-in-slice"
            sequences[seqname]["actual_taxonomy"]["strain"] = "not-in-slice"

    for seqname in sequences:
        match = "no_match"
        #sequentially get more specific on the match between correct/actual
        if (sequences[seqname]["correct_taxonomy"]["class"] == sequences[seqname]["actual_taxonomy"]["class"]):
            match = "class"
        if (sequences[seqname]["correct_taxonomy"]["family"] == sequences[seqname]["actual_taxonomy"]["family"]):
            match = "family"
        if (sequences[seqname]["correct_taxonomy"]["genus"] == sequences[seqname]["actual_taxonomy"]["genus"]):
            match = "genus"
        if (sequences[seqname]["correct_taxonomy"]["species"] == sequences[seqname]["actual_taxonomy"]["species"]):
            match = "species"
        if (sequences[seqname]["correct_taxonomy"]["strain"] == sequences[seqname]["actual_taxonomy"][
            "strain"]):
            match = "strain"
        print seqname, sequences[seqname]["correct_taxonomy"]["full"], sequences[seqname]["actual_taxonomy"]["full"], match

def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-otus', help='The picked otus TEXT file from QIIME. This is the output of pick_otus.py', required=True)
    parser.add_argument('-tax', help='The OTU taxonomy assignments from QIIME. This is the output of assign_taxonomy.py', required=True)
    parser.add_argument('-alltax', help='All the sequence taxonomy assignments from the main taxonomy file input used in assign_taxonomy.py')
    parser.add_argument('-alignment', help='FASTA alignment slice', required=True)
    return parser.parse_args()

if __name__ == '__main__':
    main()

