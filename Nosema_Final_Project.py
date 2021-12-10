#!/usr/bin/env python3

################## Understanding gene lengths and identifying proteins with sequences >4000 for Nosema ceranae

from Bio import SeqIO

sizes = [len(rec) for rec in SeqIO.parse("MicrosporidiaDB-55_NceranaePA08_1199_AnnotatedCDSs.fasta", "fasta")]

len(sizes), min(sizes), max(sizes)

print("The number of gene sequences for Nosema ceranae is : " + str(len(sizes)))
print("The smallest gene sequence length for Nosema ceranae is : " + str(min(sizes)))
print("The largest gene sequence length for Nosema ceranae is : " + str(max(sizes)))

# We see that Nosema ceranae has a high max(sizes)

# Figuring out how many sequences are greater than 4,000 bp
  
# initializing k
k = 4000

count = 0
for i in sizes :
    if i > k :
        count = count + 1
  
# printing the intersection 
print ("The sequence lengths greater than 4000 for Nosema ceranae : " + str(count))

# initializing j
j = 400

count = 0
for i in sizes :
    if i < j :
        count = count + 1
  
# printing the intersection 
print ("The sequence lengths less than 400 for Nosema ceranae : " + str(count))

from Bio import SeqIO
for seq_record in SeqIO.parse("MicrosporidiaDB-55_NceranaePA08_1199_AnnotatedCDSs.fasta", format = "fasta"):
    if len(seq_record) > 4000:
        print(seq_record.description)
                     
import pylab
import matplotlib.pyplot as plt

plt.hist(sizes, bins=20)

plt.title(
    "%i Nosema ceranae sequences\nLengths %i to %i" % (len(sizes), min(sizes), 4000)
)
# Set cut off to 4000 to clean up the plots, only a few sequences fell outside of this

plt.xlabel("Sequence length (bp)")
plt.ylabel("Count")
plt.xlim([100, 4000])
plt.show()
plt.savefig("NCHistogram.png")
plt.close() #closes figure

############################################################################# GC content for Nosema ceranae

from Bio import SeqIO
from Bio.SeqUtils import GC

gc_values = sorted(GC(rec.seq) for rec in SeqIO.parse("MicrosporidiaDB-55_NceranaePA08_1199_AnnotatedCDSs.fasta", "fasta"))

import numpy as np

print("The average GC content for Nosema ceranae is:\n", np.mean(gc_values))

print("The standard deviation for GC content for Nosema ceranae is:\n", np.std(gc_values))

print("The lowest GC content for Nosema ceranae is:\n", np.min(gc_values))

print("The highest GC content for Nosema ceranae is:\n", np.max(gc_values))

import pylab
import matplotlib.pyplot as plt

plt.plot(gc_values)
plt.title(
    "%i Nosema ceranae sequences\nGC%% %0.1f to %0.1f"
    % (len(gc_values), min(gc_values), max(gc_values))
)
plt.xlabel("Genes")
plt.ylabel("GC%")
plt.show()
plt.savefig("GCNCeranae.png") 
plt.close()

################ Understanding gene lengths and identifying proteins with sequences >4000 for Nosema bombycis

from Bio import SeqIO

sizes = [len(rec) for rec in SeqIO.parse("MicrosporidiaDB-55_NbombycisCQ1_AnnotatedCDSs.fasta", "fasta")]

len(sizes), min(sizes), max(sizes)

print("The number of gene sequences for Nosema bombycis is : " + str(len(sizes)))
print("The smallest gene sequence length for Nosema bombycis is : " + str(min(sizes)))
print("The largest gene sequence length for Nosema bombycis is : " + str(max(sizes)))

# Figuring out how many sequences are greater than 4,000 bp
  
# initializing k
k = 4000

count = 0
for i in sizes :
    if i > k :
        count = count + 1

# printing the intersection 
print ("The sequence lengths greater than 4000 for Nosema bombycis : " + str(count))

# initializing j
j = 400

count = 0
for i in sizes :
    if i < j :
        count = count + 1
  
# printing the intersection 
print ("The sequence lengths less than 400 for Nosema bombycis : " + str(count))

from Bio import SeqIO
for seq_record in SeqIO.parse("MicrosporidiaDB-55_NbombycisCQ1_AnnotatedCDSs.fasta", format = "fasta"):
    if len(seq_record) > 4000:
        print(seq_record.description)


import pylab
import matplotlib.pyplot as plt

plt.hist(sizes, bins=20)

plt.title(
    "%i Nosema bombycis sequences\nLengths %i to %i" % (len(sizes), min(sizes), 4000)
)
    
plt.xlabel("Sequence length (bp)")
plt.ylabel("Count")
plt.xlim([100, 4000])
plt.show()
plt.savefig("NBHistogram.png") 
plt.close()

############################################################################# GC content for Nosema bombycis


from Bio import SeqIO
from Bio.SeqUtils import GC

gc_values = sorted(GC(rec.seq) for rec in SeqIO.parse("MicrosporidiaDB-55_NbombycisCQ1_AnnotatedCDSs.fasta", "fasta"))

import numpy as np

print("The average GC content for Nosema bombysics is:\n", np.mean(gc_values))

print("The standard deviation for GC content for Nosema bombysics is:\n", np.std(gc_values))

print("The lowest GC content for Nosema bombysics is:\n", np.min(gc_values))

print("The highest GC content for Nosema bombysics is:\n", np.max(gc_values))

import pylab
import matplotlib.pyplot as plt

plt.plot(gc_values)
plt.title(
    "%i Nosema bombycis sequences\nGC%% %0.1f to %0.1f"
    % (len(gc_values), min(gc_values), max(gc_values))
)
plt.xlabel("Genes")
plt.ylabel("GC%")
plt.show()
plt.savefig("GCNBombycis.png")
plt.close()

########################################################################## verify that all plots were created

#verify that NCHistogram.png was created

from pathlib import Path

path_to_file = 'NCHistogram.png'
path = Path(path_to_file)

if path.is_file():
    print(f'The file {path_to_file} exists')
else:
    print(f'The file {path_to_file} does not exist')

#verify that NBHistogram.png was created
    
from pathlib import Path

path_to_file1 = 'NBHistogram.png'
path = Path(path_to_file1)

if path.is_file():
    print(f'The file {path_to_file1} exists')
else:
    print(f'The file {path_to_file1} does not exist')

#verify that GCNCeranae.png was created

from pathlib import Path

path_to_file2 = 'GCNCeranae.png'
path = Path(path_to_file2)

if path.is_file():
    print(f'The file {path_to_file2} exists')
else:
    print(f'The file {path_to_file2} does not exist')

#verify that GCNBombycis.png was created

from pathlib import Path

path_to_file3 = 'GCNBombycis.png'
path = Path(path_to_file3)

if path.is_file():
    print(f'The file {path_to_file3} exists')
else:
    print(f'The file {path_to_file3} does not exist')
    
print('Printing in a Nutshell', end='\n * ')
print('Calling Print', end='\n * ')
print('Separating Multiple Arguments', end='\n * ')
print('Preventing Line Breaks')
