from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast import NCBIWWW
from Bio.Align.Applications import MafftCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
from os import path
from Bio.SeqUtils import IUPACData
import numpy as math
print(math.log2(15/30))

for record in SeqIO.parse('alb.fasta', "fasta"):
    originalGene = record

searchResults = []
similarGenes = []
if not path.exists("similarData.xml"):
    result_handler = NCBIWWW.qblast(
        program="blastp",
        database="swissprot",
        sequence=originalGene.seq,
        perc_ident=80,
        entrez_query='mammals[Organism]')

    with open("similarData.xml", "w") as out_handle:
        out_handle.write(result_handler.read())
    result_handler.close()
searchResults = SearchIO.read("similarData.xml", "blast-xml")

with open("similarData.fasta", "w") as output:
    for hit in searchResults.hits:
        similarGenes.append(hit)
        SeqIO.write(hit.hsps[0].hit, output, "fasta")

scriptLocation = "mafft-win\mafft.bat"
cliTool = MafftCommandline(scriptLocation, input="similarData.fasta")
stdout, _ = cliTool()
with open("alignedSequences.fasta", "w") as output:
    output.write(stdout)

alignment = AlignIO.read("alignedSequences.fasta", "fasta")
alignment._alphabet = IUPACData.protein_letters
constructor = DistanceTreeConstructor(DistanceCalculator('blosum62'))
tree = constructor.build_tree(alignment)
Phylo.draw(tree)

length = len(alignment[0])
maxIndex = -1
minIndex = -1
maxIC = -10e4
minIC = 10e4
ic = []
indexes = []

for i in range(length - 14):
    column = alignment[:, i:i+15]
    currentIC = 0

    for j in range(15):
        row = column[:, j]
        aminoAcidCount = {}
        for letter in IUPACData.protein_letters:
            aminoAcidCount[letter] = 0
        for letter in IUPACData.protein_letters:
            aminoAcidCount[letter] = aminoAcidCount[letter] + row.count(letter)

        for count in aminoAcidCount.values():
            if count == 0:
                continue
            currentIC += math.log2(count/len(column))
    
    ic.append(currentIC)
    indexes.append(i)
    
    if currentIC > maxIC:
        maxIC = currentIC
        maxIndex = i
    if currentIC < minIC:
        minIC = currentIC
        minIndex = i

plt.plot(indexes, ic)
plt.ylabel('information content')
plt.xlabel('index')
plt.show()

print("Similar index: " + str(maxIndex))
print("Similar: " + alignment[0].seq[maxIndex:maxIndex+15])
print("Differing index: " + str(minIndex))
print("Differing sequence: " + alignment[0].seq[minIndex:minIndex+15])