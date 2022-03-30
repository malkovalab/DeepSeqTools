import sys
sys.path.append('./.local/lib/python3.6/site-packages')

import pysam
import regex as re
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from joblib import Parallel, delayed

def dataFrameImport(directory,datatype = 'sam'):
  frames_list_samples = []
  sample_names = []

  counter=1
  for filename in os.listdir(directory):
      if filename.endswith(f".{datatype}"):
          sample_names.append(filename[:8].strip())
          path = os.path.join(directory, filename)
          globals()[f"sample{counter}"] = pd.DataFrame(({'name': x.qname, 'seq': x.seq, 'qual': x.query_qualities} for x in pysam.Samfile(path).fetch()))
          print(f'dataframe {sample_names[-1]} has {len(globals()[f"sample{counter}"])} total rows')
          frames_list_samples.append(globals()[f"sample{counter}"])
          counter+=1
  print(f"found {len(frames_list_samples)} samples in {directory}")
  print(len(sample_names), sample_names)

  return(frames_list_samples, sample_names)

def matchSequence(row, sequence):
  if (sequence in row):
    return(True)
  else:
    return(False)

def findReadLength(row):
  return(len(row))

def trimseq(row, f_seq, r_seq, consensus):

  if row.forward_primer == True:
    f_seq_index = row.seq.find(f_seq)
    if f_seq_index == -1:
      print("no forward primer found")
      row.forward_primer = False
    row.seq = row.seq[f_seq_index:]

  if row.reverse_primer == True:
    r_seq_index = row.seq.find(r_seq)
    if r_seq_index == -1:
      print("no reverse primer found")
      row.reverse_primer = False
    row.seq = row.seq[:r_seq_index+len(r_seq)]
  if (row.forward_primer == False) & (row.reverse_primer == False):
    row.seq = "no_primers"

  if re.search('('+row.seq+')'+'{e<=2}', consensus):
    row.seq = "no_excision"
  return(row)

def classify(row, consensus, classification):

  if row.classification != "other":
    return row

  if re.search('('+row.seq+')'+'{e<=2}', consensus):
  #if row.seq in consensus:
    row.classification = classification
  return row

def run_sample(frames, i):

  print("starting analysis")
  #find the length of each read
  frames[i]['read_length'] = frames[i].apply(lambda row: findReadLength(row.seq), axis=1)
  #remove reads that are 220bp or shorter
  frames[i] = frames[i][frames[i]["read_length"] > 220 ]

  #find if there is a forward or reverse primer in each read
  frames[i]['forward_primer'] = frames[i].apply(lambda row: matchSequence(row.seq, forward_primer), axis=1)
  frames[i]['reverse_primer'] = frames[i].apply(lambda row: matchSequence(row.seq, reverse_primer), axis=1)
  #take only the reads that either contain the fw or rv primer
  frames[i] = frames[i][( frames[i]['forward_primer'] == True) | (frames[i]['reverse_primer'] == True) ]

  frames[i] = frames[i].apply(lambda row: trimseq(row, forward_primer, reverse_primer, lys2_consensus), axis=1)

  print(frames[i][frames[i]['seq'] == "no_excision"])

  #find the length of each read again after trimming and remove shorter reads
  frames[i]['read_length'] = frames[i].apply(lambda row: findReadLength(row.seq), axis=1)
  frames[i] = frames[i][frames[i]["read_length"] > 220]

  #count up reads that have the same sequences under a new column "read_count"
  frames[i] = frames[i].groupby(["seq", 'forward_primer', 'reverse_primer','read_length']).count().reset_index()
  frames[i] = frames[i].rename(columns={"name": "read_count"})
  
  #create a column for classification of the new read
  frames[i]["classification"] = "other"

  #classify reads by finding if they match a known excision event 
  frames[i] = frames[i].apply(lambda row: classify(row, ex_full, "full_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_type1, "type1_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_type2, "type2_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_J1, "J1_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_type2_v2, "type2_v2_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_type5, "type5_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_type2_v1, "type2_v1_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, flipped_spacer, "flipped_spacer"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_internal1, "internal_excision1"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_J2, "J2_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_J3, "J3_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_J4, "J4_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_J5, "J5_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_J6, "J6_excision"), axis=1)
  frames[i] = frames[i].apply(lambda row: classify(row, ex_J7, "J7_excision"), axis=1)

  print(frames[i])

  #show a graph of reads by

  graph_order=["full_excision","type1_excision","type5_excision", "type2_excision", "type2_v1_excision", "type2_v2_excision", "J1_excision", "flipped_spacer", "internal_excision1", "other" ]

  fig, ax = plt.subplots(figsize=(16,8))
  ax.set(yscale="log")

  frames[i] = frames[i].astype({'read_count': 'int'})

  ax = sns.stripplot(x="classification", y="read_count", data=frames[i], jitter=0.25, order=graph_order, alpha=.75) #, size="read_count", sizes=(10,500)

  plt.savefig(f'figures/{names[i]}.png', transparent=False)
  frames[i].to_csv(f'{names[i]}.csv', index=True)

################################################################################

forward_primer="GTTCGTACCCCTCTCGAGAATA"
reverse_primer="GGTCTTTCAGATGAGAAGTGGATGG"
lys2_consensus="GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGGCCTGACTCTTATACACAAGTAGCGTCCTGAACGGAACCTTTCCCGTTTTCCAGGATCTGATCTTCCATGTTAGGAGGTCACATGGAAGATCAGATCCTGGAAAACGGGAAAGGTTCCGTTCAGGACGCTACTTGTGTATAAGAGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_type2 =     "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGCTACTTGTGTATAAGAGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_type1 =     "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGGCCTGACTCTTATACACAAGTAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_full =      "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_J1 =        "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGACGCTACTTGTGTATAAGAGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_type2_v2 =  "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_type5 =     "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGGCCTGACTCTTATACACAAGTAGCGTCCTGAACGGAACCTTTCCCGTTTTCCAGGATAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_type2_v1 =  "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGATCCTGGAAAACGGGAAAGGTTCCGTTCAGGACGCTACTTGTGTATAAGAGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
flipped_spacer="GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGGCCTGACTCTTATACACAAGTAGCGTCCTGAACGGAACCTTTCCCGTTTTCCAGGATCTGATCTTCCATGTGACCTCCTAACATGGAAGATCAGATCCTGGAAAACGGGAAAGGTTCCGTTCAGGACGCTACTTGTGTATAAGAGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_internal1 = "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGGCCTGACTCTTATACACAAGTAGCGTCCTGGAAAACGGGAAAGGTTCCGTTCAGGACGCTACTTGTGTATAAGAGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_J2 =        "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTATAAGAGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_J3 =        "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGGCCTGACTCTTATACACAAGTAGCGTCCTGAACGGAACCTTTCCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_J4 =        "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGGCCTGACTCTTATACACAAGTAGCGTCCTGAACGGAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_J5 =        "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCACGTCAGGGCCTGACTCTTATACACAAGTAGCGTCCTGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_J6 =        "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGTGTTTGCCCGTTCAGGACGCTACTTGTGTATAAGAGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"
ex_J7 =        "GTTCGTACCCCTCTCGAGAATATTTTGTTGAACCTAATAGTGCCGAAGGAAAAACAACAATTAATGTGTTTGTTACCGGTGTCACAGGATTTCTGGGCTCCTACATCCTTGCAGATTTGTTAGGACGTTCTCCAAAGAACTACAGTTTCAAAGGTTCCGTTCAGGACGCTACTTGTGTATAAGAGTCAGCGTCAGGGCCAAGGATGAAGAAGCTGCATTTGCAAGATTACAAAAGGCAGGTATCACCTATGGTACTTGGAACGAAAAATTTGCCTCAAATATTAAAGTTGTATTAGGCGATTTATCTAAAAGCCAATTTGGTCTTTCAGATGAGAAGTGGATGG"

################################################################################

frames, names =  dataFrameImport('30-614151413/00_fastq/filtered_bams', 'bam')

jobs = len(frames)

Parallel(n_jobs=18)(
    delayed(run_sample)(frames, i) for i in range(jobs)
    )