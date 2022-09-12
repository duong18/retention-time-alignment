#!/usr/bin/env python
# coding: utf-8

# In[3]:
# FASTER IMPLEMENTATION AS USING NEIGHBOR CANDIDATES partners_before, partners_after
# FASTER THAN graph_alignment_sigmul_fast_time_interval_distr_no_loop.py


import statistics
import networkx as nx
from pqdict import pqdict
from collections import defaultdict
from sortedcontainers import SortedList
import time
start_time = time.time()


delta_RT = 5.0  # minutes
#delta_med = 0.1 # minutes
delta_med = 0.5     # minutes
delta_common_run = 3
N_touch = 75        #50        #25       #100       # 75  # 50
TIME_INTERVAL = 10      # to cover each 10 minutes


# In[4]:


print "step -2"
selected_features = []
#with open("/Users/dtn074/quantification/selected_features.txt") as fi:
#with open("/Users/dtn074/quantification/29tissues_selected_features_sigmul.txt") as fi:
with open("/Users/dtn074/quantification/29tissues_selected_features_msgf_mapping_ver_04.txt") as fi:        # from the msgf+ search, using the mapping script of version 0.4
    #header = fi.readline()
    while (True):
        line = fi.readline()
        if (line == ""):
            break

        S = line.rstrip()
        selected_features.append(S)


# In[5]:


print "step -1"
import csv
import sys
import random
import scipy.stats
maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10 
    # as long as the OverflowError occurs.

    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)


#ids must be a defaultdict(list)
#filenames = {}
rows = []
#task_ids = ["2d4b38a97a6040beaec15396063a35be", "a7b993bcdcfa476db46f4a64fda05423", "b1fec8f49f834c1d8ace9425f53c06ae", "18ab0e92291a4cb5b87a57e45b7d5e0a", "d1e74871bdf643d28705a0a6b867b2c1", "52e8078eac1c4a1fa2585e971c914dc9", "be4ff0a4253c419b9940e52756cd409f", "756be9fc4c3245dab703ab78352f2604", "cf6b74ba327d4064b7ca72b1a804f158", "f6938a28c15d4dabb977c24ea134d40f", "a3810c14cab64dbfb378cb770ba42a71", "9433b4a7a2434831a1f801b320a850c2", "3b68327d25024b3787c5ec375864fbcc", "cbac6695d40b4db08fe54c025151a4a2", "6e0607be949342079b568459de5adb6f", "850409d10e7f44c49adba7d04090f447", "a79238227ad04abeb82e028432aa59e0", "9442043f110d4a2d803f0a4f18527f00", "2a73c4b9604344059b026b0b212b27cf", "e9d58777d2ce48f8a84a30b78e3faf30", "bc36aaeb17c44a27a879a3be52854402", "b90c3fae3d8245c38e133fe22d29153f", "4b5bb7012e3b4a5eb1144622a220b458", "c753f461c3574c18a5f58ea5c7b48c3e", "611f5368627146c3962a74aa7a43c440", "2c6b3bee214f4e46af08ad4900388c79", "e9250556f0e1457eac22b5a07939cd09", "c2cbb70f685d43e3944c6bee948a1add", "f0dcb283613940f19bf07f0f7dcb2762", "e0abc2d33884455caacd43ba0e118390"]
sample_folders = ["adrenal_gland_matched", "appendix_vermiform_appendix_matched", "bone_marrow_matched", "brain_matched", "colon_matched", "duodenum_matched", "endometrium_uterine_endometrium_matched", "esophagus_matched", "fallopian_tube_oviduct_matched", "fat_adipose_tissue_matched", "gallbladder_matched", "heart_matched", "kidney_matched", "liver_matched", "lung_matched", "lymph_node_matched", "ovary_matched", "pancreas_matched", "pituitary_hypophysis_matched", "placenta_matched", "prostate_matched", "rectum_matched", "salivary_gland_matched", "small_intestine_matched", "smooth_muscle_matched", "spleen_matched", "stomach_matched", "testis_matched", "thyroid_matched", "urinary_bladder_matched"]

#for task in task_ids:
#for sf in sample_folders:
#sf = args.sample_folder
if (True):
    #mouse
    #evidence_file = "/data/beta-proteomics2/tasks/84d33e2940f144e09eca07979def50bc/msstats_input_evidence_label_free/688caf299fb542e3af57e1269aff2471"
    #29 tissues
    #evidence_file = "/data/beta-proteomics2/tasks/fc67dbb61e8b46bb86aa5b92a3c14ed7/msstats_input_evidence_label_free/c3562edd59a745fcbadc7c131a896d0d"
    #evidence_file = "/data/ccms-data/tasks/jswertz/fc67dbb61e8b46bb86aa5b92a3c14ed7/msstats_input_evidence_label_free/c3562edd59a745fcbadc7c131a896d0d"
    evidence_file = "/Users/dtn074/quantification/feature_mapping/evidence_files_concated.tsv"
    #swath
    #evidence_file = "/data/beta-proteomics2/tasks/5ecb619bb6a347b7bab2e1549ae53cc9/msstats_input_evidence_label_free/4607066630c344579fd6d1d6af51221c"
    #plasma
    #evidence_file = "/data/beta-proteomics2/tasks/e8fa802dca7d4ea69ff5cc18b69cc77b/msstats_input_evidence_label_free/d08f570803a042e0a5d4b629e91a5e59"
    with open(evidence_file) as f:
        nextline = f.readline()
        headers = nextline.rstrip().split('\t')
        mztab_dict = csv.DictReader(f, fieldnames = headers, delimiter = '\t')
        for r in list(mztab_dict):
            rows.append(r)

#rows = rows[:100000]
print("step 0: --- %s seconds ---" % (time.time() - start_time))
feat_dict = {}

#inf = []
pID = {}
lID = {}

selection_set = set(selected_features)
sigmul_feat_inf = {}
print len(selection_set)
all_runs = set([])
all_colon_runs = set(["01308_H04_P013387_B00_N32_R1", "01308_H03_P013387_B00_N24_R1", "01308_H02_P013387_B00_N16_R1", "01308_H01_P013387_B00_N08_R1", "01308_G04_P013387_B00_N31_R1", "01308_G03_P013387_B00_N23_R1", "01308_G02_P013387_B00_N15_R1", "01308_G01_P013387_B00_N07_R1", "01308_F04_P013387_B00_N30_R1", "01308_F03_P013387_B00_N22_R1", "01308_F02_P013387_B00_N14_R1", "01308_F01_P013387_B00_N06_R1", "01308_E04_P013387_B00_N29_R1", "01308_E03_P013387_B00_N21_R1", "01308_E02_P013387_B00_N13_R1", "01308_E01_P013387_B00_N05_R1", "01308_D05_P013387_B00_N36_R1", "01308_D04_P013387_B00_N28_R1", "01308_D03_P013387_B00_N20_R1", "01308_D02_P013387_B00_N12_R1", "01308_D01_P013387_B00_N04_R1", "01308_C05_P013387_B00_N35_R1", "01308_C04_P013387_B00_N27_R1", "01308_C03_P013387_B00_N19_R1", "01308_C02_P013387_B00_N11_R1", "01308_C01_P013387_B00_N03_R1", "01308_B05_P013387_B00_N34_R1", "01308_B04_P013387_B00_N26_R1", "01308_B03_P013387_B00_N18_R1", "01308_B02_P013387_B00_N10_R1", "01308_B01_P013387_B00_N02_R1", "01308_A05_P013387_B00_N33_R1", "01308_A04_P013387_B00_N25_R1", "01308_A03_P013387_B00_N17_R1", "01308_A02_P013387_B00_N09_R1", "01308_A01_P013387_B00_N01_R1"])
N_colon_rows = 0

for i in range(len(rows)):
    row = rows[i]
    protein = row['Proteins']
    peptide = row['Modified sequence']
    
    charge = int(row["Charge"])
    precursor = peptide + ',' + str(charge)
    RT_start = float(row["Retention time start"])
    RT_end = float(row["Retention time end"])
    RT = float(row["Retention time"])
    N_iso_peaks = int(row["Number of isotopic peaks"])
    file_name = row["Raw file"]
    mz_value = float(row["m/z"])
    intensity = float(row["Intensity"])
    f_idx = int(row["Idx"])
    pID[f_idx] = precursor
    lID[f_idx] = i
    #inf.append((precursor, RT, intensity, f_idx))
    #inf.append((precursor, RT, intensity, f_idx, RT_start, RT_end, mz_value))
    
    if (file_name not in all_colon_runs):    
        all_runs.add(file_name)
        
        if (precursor in selection_set):
            if ((precursor, file_name) in sigmul_feat_inf):
                sigmul_feat_inf[(precursor, file_name)].append((intensity, RT))
            else:
                sigmul_feat_inf[(precursor, file_name)] = [(intensity, RT)]
    else:
        N_colon_rows = N_colon_rows + 1
        
print "N_colon_rows =", N_colon_rows

pre_inf = {}
for (precursor, file_name) in sigmul_feat_inf:
    (max_intensity, best_RT) = max(sigmul_feat_inf[(precursor, file_name)])
    if (precursor in pre_inf):
        pre_inf[precursor].append((file_name, best_RT))
    else:
        pre_inf[precursor] = [(file_name, best_RT)]
for precursor in pre_inf:
    pre_inf[precursor] = sorted(pre_inf[precursor])

touch = {}  # for each key of precursor, store the list of runs that the precursor "touches"
for (p, L) in pre_inf.items():
    touch[p] = set([])
    for (file_name, RT) in L:
        touch[p].add(file_name)


# In[6]:


print("step 1: --- %s seconds ---" % (time.time() - start_time))
node_labels = list(pre_inf.keys())
count_and_label = []  # store pairs of #runs and labels
for p in node_labels:
    n_p = len(pre_inf[p])
    count_and_label.append((n_p, p))
count_and_label = sorted(count_and_label, reverse = True)  # sorted node labels by #runs as we expect labels of high values of runs would have more chance to be selected as the nodes of heaviest weight
node_labels = []
for (n_p, p) in count_and_label:
    node_labels.append(p)
    
n = len(node_labels)
print "n =", n

avg_RT = []
med_RT = []
belong_nodes = {}
for i in range(n):
    precursor = node_labels[i]
    temp = []
    summ = 0.0
    for (f, t) in pre_inf[precursor]:
        summ = summ + t
        temp.append(t)
        if f in belong_nodes:       # store all nodes that appear on run f and their RTs
            belong_nodes[f].append((t, i))
        else:
            belong_nodes[f] = [(t, i)]
    avg_RT.append(summ / len(pre_inf[precursor]))
    med_RT.append(statistics.median(temp))
    
    
partners_before = defaultdict(set)
partners_after = defaultdict(set)
for file in belong_nodes:
    belong_nodes[file] = sorted(belong_nodes[file])
    temp = belong_nodes[file]
    N_nodes = len(temp)
    
    for i in range(N_nodes-1):
        (RT_u, u) = temp[i]
        pu = node_labels[u]
        for j in range(i+1, N_nodes):
            (RT_v, v) = temp[j]
            pv = node_labels[v]
            if RT_v >= RT_u + delta_RT:
                break
            
            if v not in partners_after[u]:
                if (RT_u < RT_v) and (avg_RT[u] < avg_RT[v]) and (len(touch[pu].intersection(touch[pv])) >= delta_common_run):
                    partners_before[v].add(u)       # 'before' as avg_RT[u] < avg_RT[v]
                    partners_after[u].add(v)
                    
print("step 2: --- %s seconds ---" % (time.time() - start_time))
#print "len(all_pairs) =", len(all_pairs)


def get_weight(pi, pj, pre_inf):
    ni = len(pre_inf[pi])
    nj = len(pre_inf[pj])
    ii = 0
    jj = 0
    temp = []
    
    while ((ii < ni) and (jj < nj)):  # consider all runs, which both pi and pj appear
        (file_name_i, RT_i) = pre_inf[pi][ii]
        (file_name_j, RT_j) = pre_inf[pj][jj]
        if (file_name_i == file_name_j):
            if ((RT_i >= RT_j) or (RT_j >= RT_i + delta_RT)):
                break
            temp.append(RT_j - RT_i)
            ii = ii + 1
            jj = jj + 1
        elif (file_name_i < file_name_j):
            ii = ii + 1
        else:
            jj = jj + 1
    
    if (len(temp) >= delta_common_run):
        if ((ii >= ni) or (jj >= nj)):
            med = statistics.median(temp)
            count = 0

            for v in temp:
                if (abs(med - v) <= delta_med):
                    count = count + 1
                    
            if (count >= max(delta_common_run, len(temp)/2)):
                return (med, count)

    return (-1.0, -1)


G = nx.DiGraph()
max_weight = -1  # the weight of the heaviest edge
i = 0

while (i < n):
    pi = node_labels[i]
    ni = len(pre_inf[pi])
    
    if (max_weight < ni):
        j = 0
        while (j < n):
            pj = node_labels[j]
            nj = len(pre_inf[pj])
            if (i != j) and (max_weight < nj) and (avg_RT[i] < avg_RT[j]):
                (label, weight) = get_weight(pi, pj, pre_inf)
                if (weight > max_weight):
                    max_weight = weight
                    first_edge = [i, j]                    
            j = j + 1
    
    i = i + 1

print("step 3: --- %s seconds ---" % (time.time() - start_time))
print "first_edge =", first_edge, node_labels[first_edge[0]], node_labels[first_edge[1]]
print "max_weight =", max_weight

ref_set = set(first_edge)
coverage = defaultdict(int)
for v in first_edge:
    precursor = node_labels[v]
    #ti = int(med_RT[v] / TIME_INTERVAL)     #TODO fix the med RT
    for (run, rt) in pre_inf[precursor]:
        ti = int(rt / TIME_INTERVAL)
    #for r in touch[p]:      # list all runs that p touches
        coverage[(run, ti)] = coverage[(run, ti)] + 1
        
        
def get_predecessors(v, node_labels, pre_inf, partners_before, G):
    #n = len(node_labels)
    pv = node_labels[v]
    neighbors = []
    
    #for u in range(n):
    #    if (u != v) and (avg_RT[u] < avg_RT[v]):
    for u in partners_before[v]:
            pu = node_labels[u]
            (label, weight) = get_weight(pu, pv, pre_inf)
            if (weight > 0):
                G.add_edge(u, v, label=label, weight=weight)
                neighbors.append(u)
    
    return set(neighbors)


def get_successors(v, node_labels, pre_inf, partners_after, G):
    #n = len(node_labels)
    pv = node_labels[v]
    neighbors = []
    
    #for w in range(n):
    #    if (v != w) and (avg_RT[v] < avg_RT[w]):
    for w in partners_after[v]:
            pw = node_labels[w]
            (label, weight) = get_weight(pv, pw, pre_inf)
            if (weight > 0):
                G.add_edge(v, w, label=label, weight=weight)
                neighbors.append(w)
    
    return set(neighbors)


'''
def get_consensus_weight(v, preds, succs, ref_set, cRT, G):
    infer_RTs = []
    for u in preds[v].intersection(ref_set):
        infer_RTs.append((cRT[u]+G[u][v]['label'], u))
    for w in succs[v].intersection(ref_set):
        infer_RTs.append((cRT[w]-G[v][w]['label'], w))
    
    infer_RTs = sorted(infer_RTs)
    N = len(infer_RTs)
    consensus_weight = -1
    for i in range(N):
        (cRT_u, u) = infer_RTs[i]
        ll = i - 1
        while (ll >= 0) and (infer_RTs[i][0]-infer_RTs[ll][0] <= delta_med):
            ll = ll - 1
            
        rr = i + 1
        while (rr < N) and (infer_RTs[rr][0]-infer_RTs[i][0] <= delta_med):
            rr = rr + 1
        
        weight_u = 0
        for j in range(ll+1, rr):       # iterate all neighbor reference nodes of v within the tolerance
            w = infer_RTs[j][1]
            if w in preds[v]:
                weight_u = weight_u + G[w][v]['weight']
            if w in succs[v]:
                weight_u = weight_u + G[v][w]['weight']
        
        if weight_u > consensus_weight:
            consensus_weight = weight_u
            if u in preds[v]:
                cRT[v] = cRT_u + G[u][v]['label']
            if u in succs[v]:
                cRT[v] = cRT_u - G[v][u]['label']
    
    return consensus_weight
'''


def get_consensus_weight(v, storage_w, preds, succs, ref_set, cRT, G):
    #if v in consens_w:      # if the consensus weight of v has been calculated
    #   return consens_w[v]
    
    # calculate the consensus weight from scratch
    infer_RTs = []
    for u in preds[v].intersection(ref_set):
        infer_RTs.append((cRT[u]+G[u][v]['label'], u))
    for w in succs[v].intersection(ref_set):
        infer_RTs.append((cRT[w]-G[v][w]['label'], w))
    
    infer_RTs = sorted(infer_RTs)
    N = len(infer_RTs)
    consensus_weight = -1
    temp_w = []
    for i in range(N):
        (cRT_uv, u) = infer_RTs[i]
        ll = i - 1
        while (ll >= 0) and (infer_RTs[i][0]-infer_RTs[ll][0] <= delta_med):
            ll = ll - 1
            
        rr = i + 1
        while (rr < N) and (infer_RTs[rr][0]-infer_RTs[i][0] <= delta_med):
            rr = rr + 1
        
        weight_uv = 0
        for j in range(ll+1, rr):       # iterate all neighbor reference nodes of v within the tolerance
            w = infer_RTs[j][1]
            if w in preds[v]:
                weight_uv = weight_uv + G[w][v]['weight']
            if w in succs[v]:
                weight_uv = weight_uv + G[v][w]['weight']
        
        if u in preds[v]:
            temp_w.append([cRT_uv, weight_uv, G[u][v]['weight']])
        elif u in succs[v]:
            temp_w.append([cRT_uv, weight_uv, G[v][u]['weight']])
            
        if weight_uv > consensus_weight:
            consensus_weight = weight_uv
            cRT[v] = cRT_uv
            '''
            if u in preds[v]:
                cRT[v] = cRT_u + G[u][v]['label']
            if u in succs[v]:
                cRT[v] = cRT_u - G[v][u]['label']
            '''
    
    #consens_w[v] = consensus_weight
    storage_w[v] = SortedList(temp_w)
    
    return consensus_weight


def update_consensus_weight(u, v, heap_weight, storage_w, preds, succs, cRT, G):       # update consensus weight for u after v is added to the reference set
    if u in preds[v]:
        cRT_vu = cRT[v] - G[u][v]['label']
        w_vu = G[u][v]['weight']
    elif u in succs[v]:
        cRT_vu = cRT[v] + G[v][u]['label']
        w_vu = G[v][u]['weight']
        
    it = storage_w[u].irange([cRT_vu-delta_med, -1, -1], [cRT_vu+delta_med, maxInt, maxInt])        # the weight sum is at the second place in the triple
    consensus_weight = heap_weight      # the weight assigned before v is added to the reference set
    new_w = 0
    for i in it:
        #print("i =", i)
        #print("type(i) =", type(i))
        new_w = new_w + i[2]
        i[1] = i[1] + w_vu
        if consensus_weight < i[1]:
            consensus_weight = i[1]
            cRT[u] = i[0]
    new_w = new_w + w_vu
    
    if consensus_weight < new_w:
        consensus_weight = new_w
        cRT[u] = cRT_vu
    storage_w[u].add([cRT_vu, new_w, w_vu])
    
    return consensus_weight
    
    
preds = {}
succs = {}
first_candidates = set([])  # 1-hop neighbor nodes of the nodes of heaviest weight

for v in ref_set:
    preds[v] = get_predecessors(v, node_labels, pre_inf, partners_before, G)
    succs[v] = get_successors(v, node_labels, pre_inf, partners_after, G)
    first_candidates = first_candidates.union(preds[v]).union(succs[v])

# assign the consensus RTs for the first two references
cRT = {}        # the the final consensus RTs
[i, j] = first_edge
pi = node_labels[i]
ni = len(pre_inf[pi])
temp_RT = []
for ii in range(ni):
    (file_name_i, RT_i) = pre_inf[pi][ii]
    temp_RT.append(RT_i)
cRT[i] = statistics.median(temp_RT)     # the consensus RT of the first reference is set as the median of RTs on all of its runs
cRT[j] = cRT[i] + G[i][j]['label']

candidates = set(range(n)).difference(ref_set)
pq = pqdict({}, reverse=True)       # using "reverse" for max heaps
#consens_w = {}
storage_w = {}

for v in first_candidates.difference(ref_set):  # u -> v -> w
    if v not in preds:
        preds[v] = get_predecessors(v, node_labels, pre_inf, partners_before, G)
    if v not in succs:
        succs[v] = get_successors(v, node_labels, pre_inf, partners_after, G)
    
    consensus_weight = get_consensus_weight(v, storage_w, preds, succs, ref_set, cRT, G)
    pq.additem(v, consensus_weight)
    
    '''
    refined_weight = 0
    N_in = 0
    N_out = 0
    
    for u in preds[v].intersection(ref_set):
        #if ((u in ref_set) and (check_consistency(u, v))):
        #if (u in ref_set):
        refined_weight = refined_weight + G[u][v]['weight']
        N_in = N_in + 1
            
    for w in succs[v].intersection(ref_set):
        #if ((w in ref_set) and (check_consistency(v, w))):
        #if (w in ref_set):
        refined_weight = refined_weight + G[v][w]['weight']
        N_out = N_out + 1
    
    #if ((N_in > 0) and (N_out > 0)):
    if ((N_in > 0) or (N_out > 0)):
        pq.additem(v, refined_weight)
    '''

print "len(pq) =", len(pq)

while (len(pq) > 0):
    #(v, refined_weight) = pq.popitem()
    (v, consensus_weight) = pq.popitem()
    candidates.remove(v)
    N_unsatisfied_run = 0
    p = node_labels[v]
    #ti = int(med_RT[v] / TIME_INTERVAL)
    
    #print v, node_labels[v], preds[v].intersection(ref_set)
    #print v, node_labels[v], succs[v].intersection(ref_set)

    #for r in touch[p]:
    #    coverage[(r, ti)] = coverage[(r, ti)] + 1
    
    for (run, rt) in pre_inf[p]:
        ti = int(rt / TIME_INTERVAL)
        if coverage[(run, ti)] < N_touch:
            N_unsatisfied_run = N_unsatisfied_run + 1
        coverage[(run, ti)] = coverage[(run, ti)] + 1
    
    '''
    for r in touch[p]:
        if (coverage[(r, ti)] < N_touch):
            N_unsatisfied_run = N_unsatisfied_run + 1
    '''
    
    #if ((N_unsatisfied_run > 0) and (len(preds.intersection(ref_set)) > 0) and (len(succs.intersection(ref_set)) > 0)):
    if N_unsatisfied_run > 0:       # add v if it improves the coverity
        ref_set.add(v)  # v is selected to the reference node set
        if len(ref_set) % 100 == 0:
            print("len(ref_set) =", len(ref_set))
            print v, node_labels[v]
        
        # update weights for in-comming and out-going nodes to/from v
        if v not in preds:
            preds[v] = get_predecessors(v, node_labels, pre_inf, partners_before, G)
        if v not in succs:
            succs[v] = get_successors(v, node_labels, pre_inf, partners_after, G)
        
        for u in preds[v].union(succs[v]).intersection(candidates):
            if u not in preds:
                preds[u] = get_predecessors(u, node_labels, pre_inf, partners_before, G)
            if u not in succs:
                succs[u] = get_successors(u, node_labels, pre_inf, partners_after, G)
                
            #consensus_weight = get_consensus_weight(u, preds, succs, ref_set, cRT, G)
            #consensus_weight = update_consensus_weight(u, v, consens_w, storage_w, preds, succs, cRT, G)
            if u in pq:
                consensus_weight = update_consensus_weight(u, v, pq[u], storage_w, preds, succs, cRT, G)
                pq.updateitem(u, consensus_weight)
            else:
                consensus_weight = get_consensus_weight(u, storage_w, preds, succs, ref_set, cRT, G)
                pq[u] = consensus_weight
        
        '''
        for u in preds[v].intersection(candidates):  # and (check_consistency(u, v)):
            if (u in pq):
                pq.updateitem(u, pq[u] + G[u][v]['weight'])
            else:
                pq[u] = G[u][v]['weight']
                
        for w in succs[v].intersection(candidates):
            if (w in pq):
                pq.updateitem(w, pq[w] + G[v][w]['weight'])
            else:
                pq[w] = G[v][w]['weight']
                
        ref_set.add(v)  # v is selected to the reference node set
        print v, node_labels[v]
        #print v, node_labels[v], preds[v].intersection(ref_set)
        #print v, node_labels[v], succs[v].intersection(ref_set)
        
        for r in touch[p]:
            coverage[(r, ti)] = coverage[(r, ti)] + 1
        '''
        
        
for (k, c) in coverage.items():
    print "k, c =", k, c

print("step 4: --- %s seconds ---" % (time.time() - start_time))
print "G.number_of_nodes() =", G.number_of_nodes()
print "G.number_of_edges() =", G.number_of_edges()
#consens_RT = {}
N_ref_edges = 0

for (i, j) in G.edges():
    if (i in ref_set) and (j in ref_set):
        print i, j, node_labels[i], node_labels[j]
        N_ref_edges = N_ref_edges + 1
        pi = node_labels[i]
        ni = len(pre_inf[pi])
        pj = node_labels[j]
        nj = len(pre_inf[pj])

print("step 5: --- %s seconds ---" % (time.time() - start_time))
print "len(ref_set) =", len(ref_set)
print "N_ref_edges =", N_ref_edges


with open('/Users/dtn074/quantification/ref_set_consensus_weight_N_touch_75_ver2.txt', 'w') as fo:
    for r in ref_set:
        #fo.write(str(r) + '\t' + node_labels[r] + '\t' + str(statistics.median(consens_RT[r])) + '\n')
        fo.write(str(r) + '\t' + node_labels[r] + '\t' + str(cRT[r]) + '\n')
        
print("step 7: --- %s seconds ---" % (time.time() - start_time))

with open('/Users/dtn074/quantification/graph_info_consensus_weight_N_touch_75_ver2.txt', 'w') as fo:
    for u in G.nodes():
	fo.write(str(u) + '\t' + node_labels[u] + '\n')
        #fo.write(str(r) + '\t' + node_labels[r] + '\t' + str(statistics.median(consens_RT[r])) + '\n')
    fo.write("-----------\n")
    
    for (u, v) in G.edges():
        fo.write(str(u) + '\t' + str(v) + '\n')
