# Firstname Lastname
# NetID
# COMP 182 Spring 2021 - Homework 5, Problem 4

# You can import any standard library, as well as Numpy and Matplotlib.
# You can use helper functions from provided.py, and autograder.py,
# but they have to be copied over here.

# Your code here...
from typing import Tuple
from collections import *
from copy import *

g1 = {0: {1: 20, 2: 4, 3: 20}, 1: {2: 2, 5: 16}, 2: {3: 8, 4: 20}, 3: {4: 4, 5: 8}, 4: {1: 4}, 5: {}}

x = {}
x[2] = {}
x[2][2] = 5
#print(x)

### Functions for use in Problem 4

def infer_transmap(gen_data, epi_data, patient_id):
    """
        Infers a transmission map based on genetic
        and epidemiological data rooted at patient_id
        
        Arguments:
        gen_data -- filename with genetic data for each patient
        epi_data -- filename with epidemiological data for each patient
        patient_id -- the id of the 'patient 0'
        
        Returns:
        The most likely transmission map for the given scenario as the RDMST 
        of a weighted, directed, complete digraph
        """
    
    complete_digraph = construct_complete_weighted_digraph(gen_data, epi_data)
    return compute_rdmst(complete_digraph, patient_id)


def read_patient_sequences(filename):
    """
        Turns the bacterial DNA sequences (obtained from patients) into a list containing tuples of
        (patient ID, sequence).
        
        Arguments:
        filename -- the input file containing the sequences
        
        Returns:
        A list of (patient ID, sequence) tuples.
        """
    sequences = []
    with open(filename) as f:
        line_num = 0
        for line in f:
            if len(line) > 5:
                patient_num, sequence = line.split("\t")
                sequences.append( (int(patient_num), ''.join(e for e in sequence if e.isalnum())) )
    return sequences

def read_patient_traces(filename):
    """
        Reads the epidemiological data file and computes the pairwise epidemiological distances between patients
        
        Arguments:
        filename -- the input file containing the sequences
        
        Returns:
        A dictionary of dictionaries where dict[i][j] is the
        epidemiological distance between i and j.
    """
    trace_data = []
    patient_ids = []
    first_line = True
    with open(filename) as f:
        for line in f:
            if first_line:
                patient_ids = line.split()
                patient_ids = list(map(int, patient_ids))
                first_line = False
            elif len(line) > 5:
                trace_data.append(line.rstrip('\n'))
    return compute_pairwise_epi_distances(trace_data, patient_ids)

def compute_pairwise_gen_distances(sequences, distance_function):
    """
        Computes the pairwise genetic distances between patients (patients' isolate genomes)
        
        Arguments:
        sequences -- a list of sequences that correspond with patient id's
        distance_function -- the distance function to apply to compute the weight of the 
        edges in the returned graph
        
        Returns:
        A dictionary of dictionaries where gdist[i][j] is the
        genetic distance between i and j.
        """
    gdist = {}
    cultures = {}
    
    # Count the number of differences of each sequence
    for i in range(len(sequences)):
        patient_id = sequences[i][0]
        seq = sequences[i][1]
        if patient_id in cultures:
            cultures[patient_id].append(seq)
        else:
            cultures[patient_id] = [seq]
            gdist[patient_id] = {}
    # Add the minimum sequence score to the graph
    for pat1 in range(1, max(cultures.keys()) + 1):
        for pat2 in range(pat1 + 1, max(cultures.keys()) + 1):
            min_score = float("inf")
            for seq1 in cultures[pat1]:
                for seq2 in cultures[pat2]:
                    score = distance_function(seq1, seq2)
                    if score < min_score:
                        min_score = score
            gdist[pat1][pat2] = min_score
            gdist[pat2][pat1] = min_score
    return gdist



### HELPER FUNCTIONS. ###

def find_first_positives(trace_data):
    """
        Finds the first positive test date of each patient
        in the trace data.
        Arguments:
        trace_data -- a list of data pertaining to location
        and first positive test date
        Returns:
        A dictionary with patient id's as keys and first positive
        test date as values. The date numbering starts from 0 and
        the patient numbering starts from 1.
        """
    first_pos = {}
    for pat in range(len(trace_data[0])):
        first_pos[pat + 1] = None
        for date in range(len(trace_data)):
            if trace_data[date][pat].endswith(".5"):
                first_pos[pat + 1] = date
                break
    return first_pos



def compute_epi_distance(pid1, pid2, trace_data, first_pos1, first_pos2, patient_ids):
    """
        Computes the epidemiological distance between two patients.
        
        Arguments:
        pid1 -- the assumed donor's index in trace data
        pid2 -- the assumed recipient's index in trace data
        trace_data -- data for days of overlap and first positive cultures
        first_pos1 -- the first positive test day for pid1
        first_pos2 -- the first positive test day for pid2
        patient_ids -- an ordered list of the patient IDs given in the text file
        
        Returns:
        Finds the epidemiological distance from patient 1 to
        patient 2.
        """
    first_overlap = -1
    assumed_trans_date = -1
    pid1 = patient_ids.index(pid1)
    pid2 = patient_ids.index(pid2)
    # Find the first overlap of the two patients
    for day in range(len(trace_data)):
        if (trace_data[day][pid1] == trace_data[day][pid2]) & \
            (trace_data[day][pid1] != "0"):
            first_overlap = day
            break
    if (first_pos2 < first_overlap) | (first_overlap < 0):
        return len(trace_data) * 2 + 1
    # Find the assumed transmission date from patient 1 to patient 2
    for day in range(first_pos2, -1, -1):
        if (trace_data[day][pid1] == trace_data[day][pid2]) & \
            (trace_data[day][pid1] != "0"):
            assumed_trans_date = day
            break
    sc_recip = first_pos2 - assumed_trans_date

    if first_pos1 < assumed_trans_date:
        sc_donor = 0
    else:
        sc_donor = first_pos1 - assumed_trans_date
    return sc_donor + sc_recip



def compute_pairwise_epi_distances(trace_data, patient_ids):
    """
        Turns the patient trace data into a dictionary of pairwise 
        epidemiological distances.
        
        Arguments:
        trace_data -- a list of strings with patient trace data
        patient_ids -- ordered list of patient IDs to expect
        
        Returns:
        A dictionary of dictionaries where edist[i][j] is the
        epidemiological distance between i and j.
        """
    edist = {}
    proc_data = []
    # Reformat the trace data
    for i in range(len(trace_data)):
        temp = trace_data[i].split()[::-1]
        proc_data.append(temp)
    # Find first positive test days and remove the indication from the data
    first_pos = find_first_positives(proc_data)
    for pid in first_pos:
        day = first_pos[pid]
        proc_data[day][pid - 1] = proc_data[day][pid - 1].replace(".5", "")
    # Find the epidemiological distance between the two patients and add it
    # to the graph
    for pid1 in patient_ids:
        edist[pid1] = {}
        for pid2 in patient_ids:
            if pid1 != pid2:
                epi_dist = compute_epi_distance(pid1, pid2, proc_data,
                                                first_pos[pid1], first_pos[pid2], patient_ids)
                edist[pid1][pid2] = epi_dist
    return edist


def compute_genetic_distance(seq1, seq2):
    """
    Input: two sequences (from the list built by read_patient_sequences)
            (patient ID, sequence)
    Output: hamming distance
    """
    hamming_distance = 0
    sequence1 = seq1[1]
    sequence2 = seq2[1]
    #print(sequence1)
    #print(sequence2)
    count = 0
    # print(seq)
    while count < len(sequence1):
        if sequence1[count] != sequence2[count]:
    #         print(seq1[count])
    #         print(seq2[count])
             hamming_distance += 1
        count += 1
    return hamming_distance

seq1 = [0, 0, 1, 0, 1]
seq2 = [1, 0, 1, 0, 0]
#print(compute_genetic_distance(seq1, seq2))

ham_dist_gen_data = read_patient_sequences('patient_sequences.txt')

#print(ham_dist_gen_data[0])
#print(ham_dist_gen_data[1])

print(compute_genetic_distance(ham_dist_gen_data[0], ham_dist_gen_data[1]))
# for patient_data1 in ham_dist_gen_data:
#     for patient_data2 in ham_dist_gen_data:
#         if patient_data1 != patient_data2:
#             print(compute_genetic_distance(patient_data1, patient_data1))

def construct_complete_weighted_digraph(genetic_data, epidemiological_data):
    """
    Write a function construct_complete_weighted_digraph that takes as arguments
    the filenames of the genetic data and epidemiological data (in this order)
    and returns a complete, weighted, digraph whose nodes are the patients (use the
    patient id's for node labels) and whose edge weights are based on Equation (1)
    in the Homework 5 description. Read the epidemiological data
    from the dictionary built by read_patient_traces.

    """
    di_graph = {}

    
    #A list of (patient ID, sequence) tuples.
    #[(ID, [0101])]
    # Returns:
    #     A list of (patient ID, sequence) tuples.
    gen_seq = read_patient_sequences(genetic_data)
    print("gen_seq: " + str(gen_seq))
    # Returns:
    #     A dictionary of dictionaries where gdist[a][b] is the
    #     genetic distance between a and b.
    gen_dist = compute_pairwise_gen_distances(gen_seq, compute_genetic_distance)
    print("gen_dist: " + str(gen_dist))

    # Returns:
    #     A dictionary of dictionaries where dict[a][b] is the
    #     epidemiological distance between a and b.
    epi_dist = read_patient_traces(epidemiological_data)
    print("epi_dist: " + str(epi_dist))

    di_graph = deepcopy(epi_dist)

    for patient_a, contacts in epi_dist.items():
        for patient_b, epi_dist_cur in contacts.items():
            #print("patient a: " + str(patient_a))
            #print("patient b: " + str(patient_b))
            #print(gen_dist[patient_a][patient_b])
            gen_dist_cur = gen_dist.get(patient_a).get(patient_b)
            final_cur_dist = gen_dist_cur + ((999 * (epi_dist_cur / max(epi_dist))) / (10**5))
            print("values: " + str(list(epi_dist.values().values())))
            print("MAX E: " + str(max(epi_dist.values())))
            #print(final_cur_dist)
            #print(di_graph)
            di_graph[5][2] = 3
            #print(di_graph)
            di_graph[patient_a][patient_b] = final_cur_dist


    return di_graph


print("CONSTRUCT COMPLETE WEIGHTED DIGRAPH: " + str(construct_complete_weighted_digraph('patient_sequences.txt', 'patient_traces.txt')))

#print(infer_transmap('patient_sequences.txt', 'patient_traces.txt', 1))