#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import argparse
import csv
import sys
import random
random.seed(30)
import scipy.stats
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, mean_absolute_error
from statistics import mean
maxInt = sys.maxsize

while True:
    # decrease the maxInt value by factor 10 
    # as long as the OverflowError occurs.

    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)


# In[ ]:


K = 11  # the half windows of neighbors to predict


# In[ ]:


def arguments():
    parser = argparse.ArgumentParser(description='Do alignment and calculate the RT error for 29 tissues, each by a job')
    parser.add_argument('-e','--evidence', type = str, help='The path to the evidence file for doing the alignment')
    parser.add_argument('-r','--reference', type = str, help='The path to the reference file')
    parser.add_argument('-o','--output', type = str, help='The path to the output file which will contain the alignment result')
    #parser.set_defaults(variant_comparisons=False)
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()


args = arguments()
args_evidence = args.evidence
args_reference = args.reference
args_output = args.output

mul_precursor = set([])

run_A = "01075_A03_P010693_S00_N17_R1"
run_B = "01075_B03_P010693_S00_N18_R1"
run_B = "01075_C03_P010693_S00_N19_R1"
run_B = "01075_D03_P010693_S00_N20_R1"
run_B = "01075_E03_P010693_S00_N21_R1"
run_B = "01075_F03_P010693_S00_N22_R1"
run_B = "01075_G03_P010693_S00_N23_R1"


selected_features = []
with open("/Users/dtn074/quantification/plasma_selected_features_Es.txt") as fi:
    while (True):
        line = fi.readline()
        if (line == ""):
            break

        S = line.rstrip()
        selected_features.append(S)
selection_set = set(selected_features)


'''
with open("/Users/db/Downloads/complete_29tissues_triple_information_150.out") as fi:
#with open("/Users/db/Downloads/triple_information.out") as fi:
    header = fi.readline()
    while (True):
        line = fi.readline()
        if (line == ""):
            break

        S = line.rstrip().split()
        file_name = S[1]
        
        num = int(S[2])
        if (num > 1):
            mul_precursor.add(S[0])
            
print "len(mul_precursor) =", len(mul_precursor)
'''

# In[ ]:


#ids must be a defaultdict(list)
#filenames = {}
rows = []
#task_ids = ["2d4b38a97a6040beaec15396063a35be", "a7b993bcdcfa476db46f4a64fda05423", "b1fec8f49f834c1d8ace9425f53c06ae", "18ab0e92291a4cb5b87a57e45b7d5e0a", "d1e74871bdf643d28705a0a6b867b2c1", "52e8078eac1c4a1fa2585e971c914dc9", "be4ff0a4253c419b9940e52756cd409f", "756be9fc4c3245dab703ab78352f2604", "cf6b74ba327d4064b7ca72b1a804f158", "f6938a28c15d4dabb977c24ea134d40f", "a3810c14cab64dbfb378cb770ba42a71", "9433b4a7a2434831a1f801b320a850c2", "3b68327d25024b3787c5ec375864fbcc", "cbac6695d40b4db08fe54c025151a4a2", "6e0607be949342079b568459de5adb6f", "850409d10e7f44c49adba7d04090f447", "a79238227ad04abeb82e028432aa59e0", "9442043f110d4a2d803f0a4f18527f00", "2a73c4b9604344059b026b0b212b27cf", "e9d58777d2ce48f8a84a30b78e3faf30", "bc36aaeb17c44a27a879a3be52854402", "b90c3fae3d8245c38e133fe22d29153f", "4b5bb7012e3b4a5eb1144622a220b458", "c753f461c3574c18a5f58ea5c7b48c3e", "611f5368627146c3962a74aa7a43c440", "2c6b3bee214f4e46af08ad4900388c79", "e9250556f0e1457eac22b5a07939cd09", "c2cbb70f685d43e3944c6bee948a1add", "f0dcb283613940f19bf07f0f7dcb2762", "e0abc2d33884455caacd43ba0e118390"]
sample_folders = ["adrenal_gland_matched", "appendix_vermiform_appendix_matched", "bone_marrow_matched", "brain_matched", "colon_matched", "duodenum_matched", "endometrium_uterine_endometrium_matched", "esophagus_matched", "fallopian_tube_oviduct_matched", "fat_adipose_tissue_matched", "gallbladder_matched", "heart_matched", "kidney_matched", "liver_matched", "lung_matched", "lymph_node_matched", "ovary_matched", "pancreas_matched", "pituitary_hypophysis_matched", "placenta_matched", "prostate_matched", "rectum_matched", "salivary_gland_matched", "small_intestine_matched", "smooth_muscle_matched", "spleen_matched", "stomach_matched", "testis_matched", "thyroid_matched", "urinary_bladder_matched"]
all_colon_runs = set(["01308_H04_P013387_B00_N32_R1", "01308_H03_P013387_B00_N24_R1", "01308_H02_P013387_B00_N16_R1", "01308_H01_P013387_B00_N08_R1", "01308_G04_P013387_B00_N31_R1", "01308_G03_P013387_B00_N23_R1", "01308_G02_P013387_B00_N15_R1", "01308_G01_P013387_B00_N07_R1", "01308_F04_P013387_B00_N30_R1", "01308_F03_P013387_B00_N22_R1", "01308_F02_P013387_B00_N14_R1", "01308_F01_P013387_B00_N06_R1", "01308_E04_P013387_B00_N29_R1", "01308_E03_P013387_B00_N21_R1", "01308_E02_P013387_B00_N13_R1", "01308_E01_P013387_B00_N05_R1", "01308_D05_P013387_B00_N36_R1", "01308_D04_P013387_B00_N28_R1", "01308_D03_P013387_B00_N20_R1", "01308_D02_P013387_B00_N12_R1", "01308_D01_P013387_B00_N04_R1", "01308_C05_P013387_B00_N35_R1", "01308_C04_P013387_B00_N27_R1", "01308_C03_P013387_B00_N19_R1", "01308_C02_P013387_B00_N11_R1", "01308_C01_P013387_B00_N03_R1", "01308_B05_P013387_B00_N34_R1", "01308_B04_P013387_B00_N26_R1", "01308_B03_P013387_B00_N18_R1", "01308_B02_P013387_B00_N10_R1", "01308_B01_P013387_B00_N02_R1", "01308_A05_P013387_B00_N33_R1", "01308_A04_P013387_B00_N25_R1", "01308_A03_P013387_B00_N17_R1", "01308_A02_P013387_B00_N09_R1", "01308_A01_P013387_B00_N01_R1"])
#all_colon_runs = set(["01308_C01_P013387_B00_N03_R1", "01308_E01_P013387_B00_N05_R1"])
#all_colon_runs = set(["01308_F04_P013387_B00_N30_R1", "01308_D03_P013387_B00_N20_R1"])
#all_colon_runs = set(["01308_E03_P013387_B00_N21_R1", "01308_D03_P013387_B00_N20_R1"])
#all_colon_runs = set(["01308_B03_P013387_B00_N18_R1", "01308_C03_P013387_B00_N19_R1"])
#all_colon_runs = set(["01308_D01_P013387_B00_N04_R1", "01308_F01_P013387_B00_N06_R1"])

#for task in task_ids:
#for sf in sample_folders:
#sf = args.sample_folder
if (True):
    #evidence_file = '/data/beta-proteomics2/tasks/' + task + '/mapped_features/evidence.tsv'
    #evidence_file = "/Users/dtn074/quantification/feature_mapping/" + sf + "/evidence.tsv"
    #evidence_file = "/Users/db/Downloads/colon_matched/evidence.tsv"
    #evidence_file = "/Users/db/Downloads/unaligned_29tissues_evidence.csv"  # Evidence file from Merge Maestro IDs + unaligned Maxquant
    #evidence_file = "/Users/db/Downloads/evidence_files_concated_100ppm.tsv"  # Evidence file from MSGF+ IDs + unaligned Maxquant
    #evidence_file = "/Users/db/Downloads/pan_human_evidence_files_concated.tsv"  # Evidence file from MSGF+ IDs + unaligned Maxquant
    evidence_file = args_evidence
    with open(evidence_file) as f:
        nextline = f.readline()
        headers = nextline.rstrip().split('\t')
        mztab_dict = csv.DictReader(f, fieldnames = headers, delimiter = '\t')
        for r in list(mztab_dict):
            rows.append(r)

c_rows = list(rows)
max_intensity = {}
split_count = {}

for i in range(len(c_rows)):
    row = c_rows[i]
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
    
    if (precursor, file_name) in max_intensity:
        max_intensity[(precursor, file_name)] = max(max_intensity[(precursor, file_name)], intensity)
        split_count[(precursor, file_name)] = split_count[(precursor, file_name)] + 1
    else:
        max_intensity[(precursor, file_name)] = intensity
        split_count[(precursor, file_name)] = 1
        
rows = []
max_pos = {}
for i in range(len(c_rows)):
    row = c_rows[i]
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
    
    if (precursor, file_name) not in max_pos:
        if max_intensity[(precursor, file_name)] == intensity:
            rows.append(row)
            max_pos[(precursor, file_name)] = i


# In[ ]:


row_A = []
dict_B = {}
inf = []
pID = {}
lID = {}

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
    inf.append((precursor, RT, intensity, f_idx))
    #inf.append((precursor, RT, intensity, f_idx, RT_start, RT_end, mz_value))
    
    #if ((file_name in all_colon_runs) and ('+' not in peptide) and (';' not in peptide)):
    if ('+' not in peptide) and (';' not in peptide):
        if (file_name == run_A):
            #if precursor in selection_set:
            #if (precursor not in mul_precursor):
            if True:
                row_A.append((RT, i))
        #if (file_name == run_B):
        if (file_name in dict_B):
            #if precursor in selection_set:
            #if (precursor not in mul_precursor):
            if True:
                dict_B[file_name].append((RT, i))
        else:
            #if precursor in selection_set:
            #if (precursor not in mul_precursor):
            if True:
                dict_B[file_name] = [(RT, i)]


print "evidence_file =", evidence_file
print dict_B.keys()
print len(dict_B.keys())


# In[ ]:


iRT = {}
#with open("/Users/db/Downloads/ref_set.txt") as fi:
#with open("/Users/db/Downloads/[single_feature]ref_set.txt") as fi:
#with open("/Users/db/Downloads/ref_set_faster.txt") as fi:
#with open("/Users/db/Downloads/ref_set_faster_N_touch_100.txt") as fi:
#with open("/Users/db/Downloads/ref_set_msgfIDs_N_touch_100.txt") as fi:
#with open("/Users/db/Downloads/ref_set_msgfIDs_time_interval_distr_N_touch_75_N_loop_10000.txt") as fi:
with open(args_reference) as fi:
#with open("/Users/db/Downloads/ref_set_iRT.txt") as fi:
#with open("/Users/db/Downloads/ref_set_msgfIDs_time_interval_distr_N_touch_50_no_loop.txt") as fi:
    
#with open("/Users/db/Downloads/ref_set_29tissues_80.txt") as fi:
#with open("/Users/db/Downloads/ref_set_29tissues_100.txt") as fi:
#with open("/Users/db/Downloads/[heaviest_path]ref_set.txt") as fi:
    #header = fi.readline()
    while (True):
        line = fi.readline()
        if (line == ""):
            break

        S = line.rstrip().split('\t')
        iRT[S[1]] = float(S[2])
print iRT


# In[ ]:


def find_time_correspondence(inf, iRT, row_A):  # iRT time is in the former
    n_A = len(row_A)
    #print "n_A =", n_A
    #for t in row_A:
    #    print t

    A_time = []
    p_matched = set([])
    idx_matched = []
    fim_A = set([])

    for i in range(n_A):
        iA = row_A[i][1]
        p_A = inf[iA][0]
        RT_A = inf[iA][1]
        intensity_A = inf[iA][2]
        f_idx_A = inf[iA][3]

        if (p_A in iRT):
            #print p_A, RT_A, iRT[p_A]
            #A_time.append((RT_A, iRT[p_A], 1.0, 1.0))
            A_time.append((iRT[p_A], RT_A, 1.0, 1.0))
            p_matched.add(p_A)
            fim_A.add(f_idx_A)

    A_time = sorted(A_time, reverse=True)
    #print "p_matched =", p_matched
    #print "len(A_time) =", len(A_time)
    #print A_time[0]

    #for (tA, tB, intensity_A, intensity_B) in A_time:
    #    print tA, tB

    return A_time


def find_alignment_time(RT_A, A_time):  # A_time should be sorted Increasingly before; Predict from left to right (RT_A in on the left)
    #mini = 1000000.0
    L = 0
    R = len(A_time) - 1
    pos = len(A_time) - 1
    #print L, R

    while (L <= R):
        M = (L + R) / 2
        if (A_time[M][0] >= RT_A):
            pos = M
            R = M - 1
        else:
            L = M + 1
            
    distances = []
    for j in range(max(0, pos-K), min(pos+K-1, len(A_time))):
        d = abs(A_time[j][0] - RT_A)
        if (d < 0.000001):
            d = 0.000001
        distances.append(d)
    mean_distance = mean(distances)
    
    weight = 0.0
    summ = 0.0
    for j in range(max(0, pos-K), min(pos+K-1, len(A_time))):
        d = abs(A_time[j][0] - RT_A)
        if (d < 0.000001):
            d = 0.000001
        d = d + mean_distance
        
        #weight = weight + 1.0 / abs(A_time[j][0] - RT_A)
        #summ = summ + 1.0/ abs(A_time[j][0] - RT_A) * (A_time[j][1] - A_time[j][0])
        #weight = weight + A_time[j][2] / abs(A_time[j][0] - RT_A)
        #summ = summ + A_time[j][2] / abs(A_time[j][0] - RT_A) * (A_time[j][1] - A_time[j][0])
        weight = weight + A_time[j][2] / d  # intensity of the left
        #summ = summ + (A_time[j][2] / d) * (A_time[j][1] - A_time[j][0])
        summ = summ + (A_time[j][2] / d) * A_time[j][1]

    #RT_B = RT_A + summ/weight
    RT_B = summ/weight
    
    """
    print "RT_A, A_time[pos], pos, len(A_time) =", RT_A, A_time[pos], pos, len(A_time)
    print "summ, weight =", summ, weight
    print "A_time[(pos-5):pos] =", A_time[(pos-5):pos]
    print "A_time[pos:(pos+5)] =", A_time[pos:(pos+5)]
    """

    return RT_B


def longest_increasing_sequence(A):
    A = sorted(A)
    n = len(A)
    if (n <= 1):
        return A

    F = []
    prev = []

    for i in range(n):
        F.append(0)
        prev.append(-1)

    for i in range(n):
        (tAi, tBi, intensity_Ai, intensity_Bi) = A[i]
        maxi = 0
        pos = -1

        for j in range(i):
            (tAj, tBj, intensity_Aj, intensity_Bj) = A[j]
            if ((tAi >= tAj) and (tBi >= tBj)):
            #if ((tAi > tAj) and (tBi > tBj)):
                if (maxi < F[j]):
                    maxi = F[j]
                    pos = j

        F[i] = maxi + 1
        prev[i] = pos

    maxi = 0
    for i in range(n):
        if (F[i] > maxi):
            maxi = F[i]
            pos = i

    r_seq = []
    while (pos >= 0):
        (tA, tB, intensity_A, intensity_B) = A[pos]
        r_seq.append((tA, tB, intensity_A, intensity_B))
        pos = prev[pos]

    r_seq = sorted(r_seq)

    return r_seq


def pairs_reverse(A_time):
    result = []
    for (tA, tB, intensity_A, intensity_B) in A_time:
        result.append((tB, tA, intensity_B, intensity_A))
    
    return result


# #with open('/Users/db/Downloads/MSE.out', 'w') as fo:
# #with open('/Users/db/Downloads/[heaviest_path]MSE.out', 'w') as fo:
# def MSE_calculation(fo, inf, iRT, run_A, row_A, dict_B):
#     print "Calculating MSEs for", run_A
#     
#     ori_MSE = []
#     from_A_MSE = []
#     from_ref_MSE = []
#     
#     for (run_B, row_B) in dict_B.items():
#         print len(row_A)
#         print len(row_B)
# 
#         s_A = set([])
#         for (RT, i) in row_A:
#             p_A = inf[i][0]
#             s_A.add(p_A)
# 
#         s_B = set([])
#         for (RT, i) in row_B:
#             p_B = inf[i][0]
#             s_B.add(p_B)
# 
#         #print len(s_A)
#         #print len(s_B)                
# 
# 
#         I = list(s_A.intersection(s_B))
#         print "len(I) =", len(I)
#         #print I
# 
#         random.seed(30)
#         random.shuffle(I)
#         n_m_I = int(0.7 * len(I))  # change to different percentages
#         print "n_m_I =", n_m_I
#         test_precursors = set(I[n_m_I:])
#         #print "len(test_precursors) =", len(test_precursors)
# 
#         I_B = []
#         for (RT, iB) in row_B:
#             p_B = inf[iB][0]
#             if (p_B in iRT):
#                 I_B.append(p_B)
#         print "len(I_B) =", len(I_B)
#         #print I_B
# 
#         random.shuffle(I_B)
#         n_m_I_B = int(0.7 * len(I_B))  # change to different percentages
#         print "n_m_I_B =", n_m_I_B
#         test_precursors_B = set(I_B[n_m_I_B:])
#         print "len(I_B[n_m_I_B:]) =", len(I_B[n_m_I_B:])
#         print "len(test_precursors_B) =", len(test_precursors_B)
# 
#         t_row_A = []
#         x_A = {}
#         for (RT, iA) in row_A:
#             p_A = inf[iA][0]
#             if (p_A not in test_precursors):
#                 t_row_A.append((RT, iA))
#             x_A[p_A] = RT
# 
#         t_row_B = []
#         x_B = {}
#         for (RT, iB) in row_B:
#             p_B = inf[iB][0]
#             #if (p_B not in test_precursors):
#             if (p_B not in test_precursors_B):
#                 t_row_B.append((RT, iB))
#             x_B[p_B] = RT
# 
#         row_A = sorted(t_row_A)
#         row_B = sorted(t_row_B)
#         n_A = len(row_A)
#         n_B = len(row_B)
#         print "n_A =", n_A
#         print "n_B =", n_B
#         #for t in row_A:
#         #    print t
#         #print row_B
# 
# 
#         A_time = find_time_correspondence(inf, iRT, row_A)
#         #A_time = sorted(A_time, reverse=True)
# 
#         B_time = find_time_correspondence(inf, iRT, row_B)
#         #B_time = sorted(B_time, reverse=True)
# 
#         #print "p_matched =", p_matched
#         #print A_time[0]
# 
#         #for (tA, tB, intensity_A, intensity_B) in A_time:
#         #    print tA, tB
# 
#         A_time = sorted(longest_increasing_sequence(A_time))
#         B_time = sorted(longest_increasing_sequence(B_time))
#         print "len(A_time) =", len(A_time)
#         print "len(B_time) =", len(B_time)
# 
# 
# 
#         A_to_B = []
#         for p in set(I).difference(test_precursors):
#             A_to_B.append((x_A[p], x_B[p], 1.0, 1.0))
#         #A_to_B = sorted(A_to_B, reverse=True)
# 
#         A_to_B = sorted(longest_increasing_sequence(A_to_B))
#         print "len(A_to_B) =", len(A_to_B)
#         #for (tA, tB, intensity_A, intensity_B) in A_to_B:
#         #    print tA, tB
# 
# 
# 
#         print "len(test_precursors) =", len(test_precursors)
#         print "len(test_precursors_B) =", len(test_precursors_B)
# 
#         ori_A = []
#         ori_B = []
# 
#         if (len(test_precursors) >= 30):
#             N_test = len(test_precursors)
#             for p in test_precursors:
#                 ori_A.append(x_A[p])
#                 ori_B.append(x_B[p])
#                 print "ori_A, ori_B =", p, ori_A[-1], ori_B[-1]
#             MSE = mean_absolute_error(ori_A, ori_B)
#         else:           
#             MSE = -1.0
#         print "run_A, run_B =", run_A, run_B        
#         print "mean_absolute_error(ori_A, ori_B) =", MSE        
#         ori_MSE.append(MSE)
#         fo.write("run_A, run_B =" + '\t' + run_A + '\t' + run_B + '\t' + str(MSE))
#         
# 
#         ori_B = []
#         predict_from_A = []
# 
#         if (len(test_precursors) >= 30):
#             N_test = len(test_precursors)
#             for p in test_precursors:
#                 ori_B.append(x_B[p])
#                 #predict_from_A.append(reg_A_to_B.predict(np.array([[x_A[p], 0.0]])))
#                 predict_from_A.append(find_alignment_time(x_A[p], A_to_B))
#                 print "ori_B, predict_from_A =", p, ori_B[-1], predict_from_A[-1]
#             MSE = mean_absolute_error(predict_from_A, ori_B)
#         else:            
#             MSE = -1.0
#         print "run_A, run_B =", run_A, run_B        
#         print "mean_absolute_error(predict_from_A, ori_B) =", MSE
#         from_A_MSE.append(MSE)
#         fo.write('\t' + str(MSE))
#         
# 
#         ori_B = []
#         ori_iRT = []
#         predict_from_ref = []        
#         
#         if (len(test_precursors_B) >= 30):
#             N_test = len(test_precursors_B)
#             print "run_A, run_B =", run_A, run_B, "precursor, ori RT from A, ori RT from B, RT from reference, predicted RT"
#             #for p in test_precursors_B:
#             for p in I_B:
#                 ori_B.append(x_B[p])
#                 ori_iRT.append(iRT[p])
#                 #predict_from_ref.append(reg_B.predict(np.array([[iRT[p], 0.0]])))
# 
#                 predict_from_ref.append(find_alignment_time(iRT[p], B_time))
#                 if p in x_A:
#                     print "ori_B, ori_iRT, predict_from_ref =", p, x_A[p], ori_B[-1], ori_iRT[-1], predict_from_ref[-1]
#                 else:
#                     print "ori_B, ori_iRT, predict_from_ref =", p, "-1.0", ori_B[-1], ori_iRT[-1], predict_from_ref[-1]
#             MSE = mean_absolute_error(predict_from_ref, ori_B)
#             iRT_MSE = mean_absolute_error(ori_iRT, ori_B)
#         else:            
#             MSE = -1.0
#             iRT_MSE = -1.0
#         print "run_A, run_B =", run_A, run_B
#         #print "mean_squared_error(predict_from_ref, ori_B) =", mean_squared_error(predict_from_ref, ori_B)
#         print "mean_absolute_error(predict_from_ref, ori_B) =", MSE
#         #fo.write("run_A, run_B =" + '\t' + run_A + '\t' + run_B + '\t' + "mean_squared_error(predict_from_ref, ori_B) =" + '\t' + str(MSE0) + '\t' + str(MSE1) + '\n')
#         #fo.write("run_A, run_B =" + '\t' + run_A + '\t' + run_B + '\t' + "mean_squared_error(predict_from_ref, ori_B) =" + '\t' + str(MSE) + '\n')
#         #fo.write("run_A, run_B =" + '\t' + run_A + '\t' + run_B + '\t' + "mean_absolute_error(predict_from_ref, ori_B) =" + '\t' + str(MSE) + '\n')
#         from_ref_MSE.append(MSE)
#         fo.write('\t' + str(MSE) + '\t' + str(iRT_MSE) + '\n')
#         
#     return ori_MSE, from_A_MSE, from_ref_MSE

# In[ ]:


#with open('/Users/db/Downloads/MSE.out', 'w') as fo:
#with open('/Users/db/Downloads/[heaviest_path]MSE.out', 'w') as fo:
def MSE_calculation(fo, fh, inf, iRT, run_A, row_A, dict_B):
    print "Calculating MSEs for", run_A
    
    ori_row_A = list(row_A)
    ori_MSE = []
    from_A_MSE = []
    from_ref_MSE = []
    
    for (run_B, row_B) in dict_B.items():
        print len(row_A)
        print len(row_B)
        
        row_A = list(ori_row_A)
        s_A = set([])
        for (RT, i) in row_A:
            p_A = inf[i][0]
            s_A.add(p_A)

        s_B = set([])
        for (RT, i) in row_B:
            p_B = inf[i][0]
            s_B.add(p_B)

        #print len(s_A)
        #print len(s_B)                


        I = list(s_A.intersection(s_B))
        print "len(I) =", len(I)
        #print I
        
        random.shuffle(I)
        n_m_I = int(0.7 * len(I))  # change to different percentages
        print "n_m_I =", n_m_I
        test_precursors = set(I[n_m_I:])
        #print "len(test_precursors) =", len(test_precursors)

        I_B = []
        for (RT, iB) in row_B:
            p_B = inf[iB][0]
            if (p_B in iRT):
                I_B.append(p_B)
        print "len(I_B) =", len(I_B)
        #print I_B

        random.shuffle(I_B)
        n_m_I_B = int(0.7 * len(I_B))  # change to different percentages
        print "n_m_I_B =", n_m_I_B
        test_precursors_B = set(I_B[n_m_I_B:])
        print "len(I_B[n_m_I_B:]) =", len(I_B[n_m_I_B:])
        print "len(test_precursors_B) =", len(test_precursors_B)

        t_row_A = []
        x_A = {}
        for (RT, iA) in row_A:
            p_A = inf[iA][0]
            if (p_A not in test_precursors):
                t_row_A.append((RT, iA))
            x_A[p_A] = RT

        t_row_B = []
        x_B = {}
        for (RT, iB) in row_B:
            p_B = inf[iB][0]
            #if (p_B not in test_precursors):
            if (p_B not in test_precursors_B):
                t_row_B.append((RT, iB))
            x_B[p_B] = RT

        row_A = sorted(t_row_A)
        row_B = sorted(t_row_B)
        n_A = len(row_A)
        n_B = len(row_B)
        print "n_A =", n_A
        print "n_B =", n_B
        #for t in row_A:
        #    print t
        #print row_B


        A_time = find_time_correspondence(inf, iRT, row_A)
        #A_time = sorted(A_time, reverse=True)

        B_time = find_time_correspondence(inf, iRT, row_B)
        #B_time = sorted(B_time, reverse=True)

        #print "p_matched =", p_matched
        #print A_time[0]

        #for (tA, tB, intensity_A, intensity_B) in A_time:
        #    print tA, tB

        A_time = sorted(longest_increasing_sequence(A_time))
        B_time = sorted(longest_increasing_sequence(B_time))
        print "len(A_time) =", len(A_time)
        print "len(B_time) =", len(B_time)



        A_to_B = []
        for p in set(I).difference(test_precursors):
            A_to_B.append((x_A[p], x_B[p], 1.0, 1.0))
        #A_to_B = sorted(A_to_B, reverse=True)

        A_to_B = sorted(longest_increasing_sequence(A_to_B))
        print "len(A_to_B) =", len(A_to_B)
        #for (tA, tB, intensity_A, intensity_B) in A_to_B:
        #    print tA, tB



        print "len(test_precursors) =", len(test_precursors)
        print "len(test_precursors_B) =", len(test_precursors_B)

        ori_A = []
        ori_B = []

        if (len(test_precursors) >= 30):
            N_test = len(test_precursors)
            for p in test_precursors:
                ori_A.append(x_A[p])
                ori_B.append(x_B[p])
                #print "ori_A, ori_B =", p, ori_A[-1], ori_B[-1]
            MSE = mean_absolute_error(ori_A, ori_B)
        else:           
            MSE = -1.0
        print "run_A, run_B =", run_A, run_B        
        print "mean_absolute_error(ori_A, ori_B) =", MSE        
        ori_MSE.append(MSE)
        fo.write("run_A, run_B =" + '\t' + run_A + '\t' + run_B + '\t' + str(MSE))
        

        ori_B = []
        predict_from_A = []
        align_for_A = []
        align_for_B = []
        '''__________'''
        print "pairs_reverse(A_time)"
        #print pairs_reverse(A_time)
        print "pairs_reverse(B_time)"
        #print pairs_reverse(B_time)

        if (len(test_precursors) >= 30):
            N_test = len(test_precursors)
            for p in test_precursors:                
                ori_B.append(x_B[p])
                #predict_from_A.append(reg_A_to_B.predict(np.array([[x_A[p], 0.0]])))
                predict_from_A.append(find_alignment_time(x_A[p], A_to_B))
                #print "ori_B, predict_from_A =", p, ori_B[-1], predict_from_A[-1]
                
                
                align_for_A.append(find_alignment_time(x_A[p], pairs_reverse(A_time)))
                align_for_B.append(find_alignment_time(x_B[p], pairs_reverse(B_time)))
                
                if p in iRT:
                    fh.write(p + '\t' + str(split_count[(p, run_A)]) + '\t' + str(split_count[(p, run_B)]) + '\t' + run_A + '\t' + run_B + '\t' + str(x_A[p]) + '\t' + str(x_B[p]) + '\t' + str(align_for_A[-1]) + '\t' + str(align_for_B[-1]) + '\t' + str(iRT[p]) + '\n')
                else:
                    fh.write(p + '\t' + str(split_count[(p, run_A)]) + '\t' + str(split_count[(p, run_B)]) + '\t' + run_A + '\t' + run_B + '\t' + str(x_A[p]) + '\t' + str(x_B[p]) + '\t' + str(align_for_A[-1]) + '\t' + str(align_for_B[-1]) + '\t' + "-1.0" + '\n')
                """                
                if p in iRT:
                    align_for_A.append(find_alignment_time(x_A[p], pairs_reverse(A_time)))
                    align_for_B.append(iRT[p])
                    ori_B.append(find_alignment_time(x_B[p], pairs_reverse(B_time)))
                    predict_from_A.append(iRT[p])
                    
                    print "p =", p
                    print "align_for_A[-1] =", align_for_A[-1]
                    print "ori_B[-1] =", ori_B[-1]
                    print "iRT[p] =", iRT[p]
                    print "find_alignment_time(iRT[p], A_time) =", find_alignment_time(iRT[p], A_time)
                    print "find_alignment_time(iRT[p], B_time) =", find_alignment_time(iRT[p], B_time)
                    print "x_A[p] =", x_A[p]
                    print "x_B[p] =", x_B[p]
                """    
                
            MSE = mean_absolute_error(predict_from_A, ori_B)
            MAE_after_align = mean_absolute_error(align_for_A, align_for_B)
        else:            
            MSE = -1.0
            MAE_after_align = -1.0
        print "run_A, run_B =", run_A, run_B        
        print "mean_absolute_error(predict_from_A, ori_B) =", MSE
        print "mean_absolute_error(align_for_A, align_for_B) =", MAE_after_align
        from_A_MSE.append(MSE)
        fo.write('\t' + str(MAE_after_align) + '\t' + str(MSE))
        

        ori_B = []
        ori_iRT = []
        predict_from_ref = []        
        
        if (len(test_precursors_B) >= 30):
            N_test = len(test_precursors_B)
            print "run_A, run_B =", run_A, run_B, "precursor, ori RT from A, ori RT from B, RT from reference, predicted RT"
            for p in test_precursors_B:
            #for p in I_B:
                ori_B.append(x_B[p])
                ori_iRT.append(iRT[p])
                #predict_from_ref.append(reg_B.predict(np.array([[iRT[p], 0.0]])))

                predict_from_ref.append(find_alignment_time(iRT[p], B_time))
                """
                if p in x_A:
                    print "ori_B, ori_iRT, predict_from_ref =", p, x_A[p], ori_B[-1], ori_iRT[-1], predict_from_ref[-1]
                else:
                    print "ori_B, ori_iRT, predict_from_ref =", p, "-1.0", ori_B[-1], ori_iRT[-1], predict_from_ref[-1]
                """
            MSE = mean_absolute_error(predict_from_ref, ori_B)
            iRT_MSE = mean_absolute_error(ori_iRT, ori_B)
        else:            
            MSE = -1.0
            iRT_MSE = -1.0
        print "run_A, run_B =", run_A, run_B
        #print "mean_squared_error(predict_from_ref, ori_B) =", mean_squared_error(predict_from_ref, ori_B)
        print "mean_absolute_error(predict_from_ref, ori_B) =", MSE
        #fo.write("run_A, run_B =" + '\t' + run_A + '\t' + run_B + '\t' + "mean_squared_error(predict_from_ref, ori_B) =" + '\t' + str(MSE0) + '\t' + str(MSE1) + '\n')
        #fo.write("run_A, run_B =" + '\t' + run_A + '\t' + run_B + '\t' + "mean_squared_error(predict_from_ref, ori_B) =" + '\t' + str(MSE) + '\n')
        #fo.write("run_A, run_B =" + '\t' + run_A + '\t' + run_B + '\t' + "mean_absolute_error(predict_from_ref, ori_B) =" + '\t' + str(MSE) + '\n')
        from_ref_MSE.append(MSE)
        fo.write('\t' + str(MSE) + '\t' + str(iRT_MSE) + '\n')
        
    return ori_MSE, from_A_MSE, from_ref_MSE


# In[ ]:


#name_set_A = np.random.choice(dict_B.keys(), size=20, replace=False)
name_set_A = list(dict_B.keys())
random.shuffle(name_set_A)
#fo = open('/Users/db/Downloads/[heaviest_path]MSE.out', 'w')
#fo = open('/Users/db/Downloads/[singmul_feat_ref]MSE.out', 'w')
#fo = open('/Users/db/Downloads/[singmul_feat_ref]MSE_for_sing_only.out', 'w')
#fo = open('/Users/db/Downloads/[sing_feat_ref]MSE_for_sing_only.out', 'w')
#fo = open('/Users/db/Downloads/[singmul_feat_ref_faster]MAE.out', 'w')
#fo = open('/Users/db/Downloads/[singmul_feat_ref_faster_N_touch_100]mean_absolute_error.out', 'w')
#fo = open('/Users/db/Downloads/[singmul_feat_ref_set_msgfIDs_N_touch_100]mean_absolute_error.out', 'w')
#fo = open('/Users/db/Downloads/[singmul_feat_ref_set_msgfIDs_time_interval_distr_N_touch_75_N_loop_10000]mean_absolute_error.out', 'w')
#fo = open('/Users/db/Downloads/[pan_human_singmul_feat_ref_set_msgfIDs_time_interval_distr_N_touch_75_N_loop_10000]mean_absolute_error.out', 'w')
fo = open(args_output, 'w')
#fo = open('/Users/db/Downloads/[pan_human_singmul_feat_ref_set_iRT]mean_absolute_error.out', 'w')
#fo = open('/Users/db/Downloads/[singmul_feat_ref_set_msgfIDs_time_interval_distr_N_touch_50_no_loop]mean_absolute_error.out', 'w')
fo.write("Run A\tRun B\tA and B originally\tA and B after aligning\tFrom A\tFrom Refs with alignment\tFrom Refs without alignment\n")
fh = open(args_output.replace(".txt", ".his"), 'w')
fh.write("precursor" + '\t' + "split count on run A" +  '\t' + "split count on run B" + '\t' + "run A" + '\t' + "run B" + '\t' + "RT on run A" + '\t' + "RT on run B" + '\t' + "aligned RT for A" + '\t' + "aligned RT for B" + '\t' + "ref RT" + '\n')

for run_A in name_set_A:
    row_A = []
    
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
        
        #if ((file_name in all_colon_runs) and ('+' not in peptide) and (';' not in peptide)):
        if ('+' not in peptide) and (';' not in peptide):
            if (file_name == run_A):
                #if precursor in selection_set:
                #if (precursor not in mul_precursor):
                if True:
                    row_A.append((RT, i))
    
    ori_MSE, from_A_MSE, from_ref_MSE = MSE_calculation(fo, fh, inf, iRT, run_A, row_A, dict_B)
    print "run_A =", run_A
    
    temp = [x for x in ori_MSE if (x > 0)]
    if len(temp) >= 5:
        print "ori_MSE =", mean(temp)
    else:
        print "ori_MSE =", -1.0
        
    temp = [x for x in from_A_MSE if (x > 0)]
    if len(temp) >= 5:
        print "from_A_MSE =", mean(temp)
    else:
        print "from_A_MSE =", -1.0
        
    temp = [x for x in from_ref_MSE if (x > 0)]
    if len(temp) >= 5:
        print "from_ref_MSE =", mean(temp)
    else:
        print "from_ref_MSE =", -1.0
        
fo.close()
fh.close()


# In[ ]:





# # end here
