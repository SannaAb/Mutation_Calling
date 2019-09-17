#!/usr/bin/py 

import sys 
import pandas as pd
import os 

# This script creates a matrix per each position that the user defines. Is the position mutated or not? Cords inputed as chrom:start-stop, vcf are inputed as an array 

def readdata(): 
    """
    Reading the data
    """
    #coords = "?" 
    coords = sys.argv[1]
    vcfs = sys.argv[2:]
    return (coords,vcfs)



def createthematrix(coords,vcfs):
    """
    creates the matrix, goes through each position in the the vcf files 
    """    
    startcoord = int(coords.split(":")[1].split("-")[0])
    endcoord = int(coords.split(":")[1].split("-")[-1])
    dictmut = {}
    for i in vcfs:
        with open(i, "r") as invcf: 
            for line in invcf: 
                line = line.strip()
                if "#" in line: 
                    continue 
                elif line.split("\t")[0] == coords.split(":")[0] and int(line.split("\t")[1]) >= startcoord and int(line.split("\t")[1]) <= endcoord:
                    if i in dictmut:
                        dictmut[i].append(int(line.split("\t")[1]))
                    else: 
                        dictmut[i] = [int(line.split("\t")[1])]
    endcoord+=1 
    dictForallpositions = {}
    rangecoord = range(startcoord,endcoord)
    for key, value in dictmut.iteritems():
        print "Investigating... " + key
        for i in rangecoord:
            if i in value: # We are at the position containing a mutation 
                 if key in dictForallpositions:
                     dictForallpositions[key].append('1')
                 else:
                     dictForallpositions[key] = ['1']
            else: 
                if key in dictForallpositions:
                    dictForallpositions[key].append('0')
                else: 
                    dictForallpositions[key] = ['0']
    dfObj = pd.DataFrame(dictForallpositions, index=rangecoord) # Create a matrix from the dictionary 
    outfile = coords.split(":")[0] + "_" + str(startcoord) + "_" + str(endcoord-1) + ".MutationTable.txt"  
    dfObj.to_csv(outfile, sep="\t")

    return outfile

def plotting(outfile): 
    """
    This part makes the plotting, creates an Rscript that we run 
    """

    RoutScript = outfile.split(".MutationTable.txt")[0] + ".R"
    outpdf = outfile.split(".MutationTable.txt")[0] + ".pdf"

    with open(RoutScript, "w") as Rcode:
        print >> Rcode, """
#!/usr/bin/R
library(dendextend)
library(ape)
data<-read.table("%s", sep = \"\t\", header = T, row.names=1,check.names=FALSE)
data<-t(data)
rownames(data)<-gsub(x=row.names(data),pattern=".Sorted.Passed.vcf", replacement="")
d <- dist(data, method = "binary")
hc <- hclust(d, method = "ward.D2")
dend <- as.dendrogram(hc)
pdf("%s",width=16, height=16) 
plot(as.phylo(dend), type = "fan", cex=0.8)
dev.off()
        """ %(outfile,outpdf)

    command = "Rscript %s" %RoutScript

    os.system(command)

def main(): 
    (coords,vcfs)=readdata()
    outfile=createthematrix(coords,vcfs)
    plotting(outfile)

main()

