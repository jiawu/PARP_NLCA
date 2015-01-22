import os
import glob
import pdb
import re
import pickle

messed_up_jobs = []
list_of_dirs = os.listdir("/projects/p20519/jia_output/parp_permutation_test/")
entries = []
list_of_dicts = []
T0T3list_of_dicts = []
T3T6list_of_dicts = []
T6T9list_of_dicts = []
T9T24list_of_dicts = []
average_list = []

for dir in list_of_dirs:
  try:
    fullpath = "/projects/p20519/jia_output/parp_permutation_test/" + dir + "/CNO/Scaffold.complete.sif"
    #print fullpath
    with open(fullpath, "r") as input_file:
      for line in input_file:
        entry = line.rstrip('\n')
        entries.append(entry)
    T0T3_path = "/projects/p20519/jia_output/parp_permutation_test/" + dir + "/CNO/ReportFromT0ToT3/EvolutionFit.txt"
    T3T6_path = "/projects/p20519/jia_output/parp_permutation_test/" + dir + "/CNO/ReportFromT3ToT6/EvolutionFit.txt"
    T6T9_path = "/projects/p20519/jia_output/parp_permutation_test/" + dir + "/CNO/ReportFromT6ToT9/EvolutionFit.txt"
    T9T24_path = "/projects/p20519/jia_output/parp_permutation_test/" + dir + "/CNO/ReportFromT9ToT24/EvolutionFit.txt"

    with open(T0T3_path,'rb') as fh:
      offs = -100
      while True:
        fh.seek(offs, 2)
        lines = fh.readlines()
        if len(lines)>1:
          last = lines[-1]
          break
        offs *= 2
    #print last
    entry = last.rstrip('\n').split("\t")
    score = entry[1]
    T0T3_score = score
    temp_dict = {'dir': dir, 'score': score, 'timepoint':'T0T3'}
    T0T3list_of_dicts.append(temp_dict)

    with open(T3T6_path,'rb') as fh:
      offs = -100
      while True:
        fh.seek(offs, 2)
        lines = fh.readlines()
        if len(lines)>1:
          last = lines[-1]
          break
        offs *= 2
    #print last
    entry = last.rstrip('\n').split("\t")
    score = entry[1]
    T3T6_score = score
    temp_dict = {'dir': dir, 'score': score, 'timepoint':'T3T6'}
    T3T6list_of_dicts.append(temp_dict)

    with open(T6T9_path,'rb') as fh:
      offs = -100
      while True:
        fh.seek(offs, 2)
        lines = fh.readlines()
        if len(lines)>1:
          last = lines[-1]
          break
        offs *= 2
    #print last
    entry = last.rstrip('\n').split("\t")
    T6T9_score = score
    score = entry[1]
    temp_dict = {'dir': dir, 'score': score, 'timepoint':'T6T9'}
    T6T9list_of_dicts.append(temp_dict)

    with open(T9T24_path,'rb') as fh:
      offs = -100
      while True:
        fh.seek(offs, 2)
        lines = fh.readlines()
        if len(lines)>1:
          last = lines[-1]
          break
        offs *= 2
    #print last
    entry = last.rstrip('\n').split("\t")
    T9T24_score = score
    score = entry[1]
    temp_dict = {'dir': dir, 'score': score, 'timepoint':'T9T24'}
    T9T24list_of_dicts.append(temp_dict)

    average_score = (float(T0T3_score) + float(T3T6_score) + float(T6T9_score) + float(T9T24_score))/4
    temp_dict = {'dir': dir, 'average_score' : average_score}
    average_list.append(temp_dict)
    
  except IOError:
    print "Error: does not appear to exist."
    messed_up_jobs.append(dir)

T0T3newlist = sorted(T0T3list_of_dicts, key=lambda k: k['score'])
T3T6newlist = sorted(T3T6list_of_dicts, key=lambda k: k['score'])
T6T9newlist = sorted(T6T9list_of_dicts, key=lambda k: k['score'])
T9T24newlist = sorted(T9T24list_of_dicts, key=lambda k: k['score'])
average_newlist = sorted(average_list, key=lambda k: k['average_score'])


list_of_top_networks = []
list_of_top_networks.append(T0T3newlist[0:5])
list_of_top_networks.append(T3T6newlist[0:5])
list_of_top_networks.append(T6T9newlist[0:5])
list_of_top_networks.append(T9T24newlist[0:5])

list_of_top_average_networks = []
list_of_top_average_networks.append(average_newlist)
import itertools
merged = list(itertools.chain(*list_of_top_networks))
best_networks = []
for network in merged:
  best_networks.append(network['dir'])

list_of_all_edges = []
list_of_all_edges.append(average_newlist)
all_merged = list(itertools.chain(*list_of_all_edges))
all_networks = []
for network in all_merged:
  all_networks.append(network['dir'])
all_entries = []
for dir in all_networks:
  try:
    fullpath = "/projects/p20519/jia_output/parp_permutation_test/" + dir + "/CNO/Scaffold.complete.sif"
    #print fullpath
    with open(fullpath, "r") as input_file:
      for line in input_file:
        entry = line.rstrip('\n')
        all_entries.append(entry)
  except IOError:
    print "Error: does not appear to exist."
    messed_up_jobs.append(dir)

entries2=[]
for dir in best_networks:
  try:
    fullpath = "/projects/p20519/jia_output/parp_permutation_test/" + dir + "/CNO/Scaffold.complete.sif"
    #print fullpath
    with open(fullpath, "r") as input_file:
      for line in input_file:
        entry = line.rstrip('\n')
        entries2.append(entry)
  except IOError:
    print "Error: does not appear to exist."
    messed_up_jobs.append(dir)
