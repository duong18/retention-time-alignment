# NOTE: This code can be used for alignment results for BOTH modification and unmodification features. The defaut is for unmodification only

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


K = 11  # the half windows (the half number) of neighbors to predict
TRAIN_PERCENT = 0.7     # the precentage to split the set of shared precursors when aligning a pair of runs
MIN_TEST = 30   # the minimum required number of shared precursors of a run pair to calculate the RT alignment errors
eps = 0.000001  # small epsilon value to prevent dividing by zero
MIN_REPORT = 5  # for printout purpose only
                # the minimum required number of RT differences to calculate the mean absolute error of RTs, if a run pair does not have enough return -1 instead


def arguments():
    parser = argparse.ArgumentParser(description='Do alignment and calculate the RT errors for run pairs in the evidence file')
    parser.add_argument('-e','--evidence', type = str, help='The path to the evidence file for doing this alignment')
    parser.add_argument('-r','--reference', type = str, help='The path to the reference file')
    parser.add_argument('-o','--output', type = str, help='The path to the output file which will contain the alignment result')
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()


def read_evidence_file(args_evidence):
    rows = []
    evidence_file = args_evidence
    with open(evidence_file) as f:
        nextline = f.readline()
        headers = nextline.rstrip().split('\t')
        mztab_dict = csv.DictReader(f, fieldnames = headers, delimiter = '\t')
        for r in list(mztab_dict):
            rows.append(r)
            
    return rows


def get_highest_intensity_features(rows):
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
    return rows, split_count


def index_runs(rows):
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

        # in case we care unmod features only, the feature will be ignored if one or more of the ID sequences contain mods
        c_peptide = peptide.replace("C+57.", "C")
        if ('+' in c_peptide) or ('-' in c_peptide):        
            continue
        
        if ';' not in peptide:      # should not be an ambiguous ID
            if (file_name in dict_B):
                dict_B[file_name].append((RT, i))
            else:
                dict_B[file_name] = [(RT, i)]
                
    print(dict_B.keys())
    print(len(dict_B.keys()))
    return dict_B, inf


def read_references(args_reference):  
    iRT = {}
    with open(args_reference) as fi:
        while (True):
            line = fi.readline()
            if (line == ""):
                break

            S = line.rstrip().split('\t')
            iRT[S[1]] = float(S[2])
            
    print(iRT)
    return iRT


def find_time_correspondence(inf, iRT, row_A):  # return the correspondence between the reference RTs and the RTs on the actual run of the same precursors
                                                # iRT time is in the former values of the pairs returned
    n_A = len(row_A)
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
            A_time.append((iRT[p_A], RT_A, 1.0, 1.0))
            p_matched.add(p_A)
            fim_A.add(f_idx_A)

    A_time = sorted(A_time, reverse=True)
    return A_time


def find_alignment_time(RT_A, A_time):  # calculate the aligned RT for RT_A, given the time correspondence A_time
                                        # A_time should be sorted increasingly before; Predict from left to right (RT_A in on the left)
    L = 0
    R = len(A_time) - 1
    pos = len(A_time) - 1
    
    while (L <= R):
        M = (L + R) // 2
        if A_time[M][0] >= RT_A:
            pos = M
            R = M - 1
        else:
            L = M + 1
            
    distances = []
    for j in range(max(0, pos-K), min(pos+K-1, len(A_time))):
        d = max(abs(A_time[j][0] - RT_A), eps)
        distances.append(d)
    mean_distance = mean(distances)
    
    weight = 0.0
    summ = 0.0
    for j in range(max(0, pos-K), min(pos+K-1, len(A_time))):
        d = max(abs(A_time[j][0] - RT_A), eps)
        d = d + mean_distance
        weight = weight + A_time[j][2] / d  # intensity of the left
        summ = summ + (A_time[j][2] / d) * A_time[j][1]

    RT_B = summ/weight
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


def MAE_calculation(fo, fh, inf, iRT, run_A, row_A, dict_B, split_count):    # report the alignment errors (Mean Absolute Errors) when aligning run_A to others runs in the evidence file
    print("Calculating MAEs for", run_A)
    ori_row_A = list(row_A)
    ori_MAE = []
    from_A_MAE = []
    from_ref_MAE = []
    
    for (run_B, row_B) in dict_B.items():        
        row_A = list(ori_row_A)
        s_A = set([])
        for (RT, i) in row_A:
            p_A = inf[i][0]
            s_A.add(p_A)
        s_B = set([])
        for (RT, i) in row_B:
            p_B = inf[i][0]
            s_B.add(p_B)
            
        I = list(s_A.intersection(s_B))        
        random.shuffle(I)
        n_m_I = int(TRAIN_PERCENT * len(I))  # OK to change to different percentages
        test_precursors = set(I[n_m_I:])
        
        I_B = []
        for (RT, iB) in row_B:
            p_B = inf[iB][0]
            if (p_B in iRT):
                I_B.append(p_B)
        random.shuffle(I_B)
        n_m_I_B = int(TRAIN_PERCENT * len(I_B))  # OK to change to different percentages
        test_precursors_B = set(I_B[n_m_I_B:])
        
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
            if (p_B not in test_precursors):
            #if (p_B not in test_precursors_B):
                t_row_B.append((RT, iB))
            x_B[p_B] = RT

        row_A = sorted(t_row_A)
        row_B = sorted(t_row_B)
        n_A = len(row_A)
        n_B = len(row_B)
        A_time = find_time_correspondence(inf, iRT, row_A)
        B_time = find_time_correspondence(inf, iRT, row_B)
        A_time = sorted(longest_increasing_sequence(A_time))
        B_time = sorted(longest_increasing_sequence(B_time))        
        A_to_B = []
        for p in set(I).difference(test_precursors):
            A_to_B.append((x_A[p], x_B[p], 1.0, 1.0))
        A_to_B = sorted(longest_increasing_sequence(A_to_B))
        
        ori_A = []
        ori_B = []
        if (len(test_precursors) >= MIN_TEST):
            N_test = len(test_precursors)
            for p in test_precursors:
                ori_A.append(x_A[p])
                ori_B.append(x_B[p])
            MAE = mean_absolute_error(ori_A, ori_B)
        else:           
            MAE = -1.0        
        ori_MAE.append(MAE)
        fo.write(run_A + '\t' + run_B + '\t' + str(MAE))
        
        ori_B = []
        predict_from_A = []
        align_for_A = []
        align_for_B = []
        if (len(test_precursors) >= MIN_TEST):
            N_test = len(test_precursors)
            for p in test_precursors:                
                ori_B.append(x_B[p])
                predict_from_A.append(find_alignment_time(x_A[p], A_to_B))             
                align_for_A.append(find_alignment_time(x_A[p], pairs_reverse(A_time)))
                align_for_B.append(find_alignment_time(x_B[p], pairs_reverse(B_time)))
                
                if p in iRT:
                    fh.write(p + '\t' + str(split_count[(p, run_A)]) + '\t' + str(split_count[(p, run_B)]) + '\t' + run_A + '\t' + run_B + '\t' + str(x_A[p]) + '\t' + str(x_B[p]) + '\t' + str(align_for_A[-1]) + '\t' + str(align_for_B[-1]) + '\t' + str(iRT[p]) + '\n')
                else:
                    fh.write(p + '\t' + str(split_count[(p, run_A)]) + '\t' + str(split_count[(p, run_B)]) + '\t' + run_A + '\t' + run_B + '\t' + str(x_A[p]) + '\t' + str(x_B[p]) + '\t' + str(align_for_A[-1]) + '\t' + str(align_for_B[-1]) + '\t' + "-1.0" + '\n')
                    
            MAE = mean_absolute_error(predict_from_A, ori_B)
            MAE_after_align = mean_absolute_error(align_for_A, align_for_B)
        else:            
            MAE = -1.0
            MAE_after_align = -1.0
        from_A_MAE.append(MAE)
        fo.write('\t' + str(MAE_after_align) + '\t' + str(MAE))        

        ori_B = []
        ori_iRT = []
        predict_from_ref = []
        if (len(test_precursors_B) >= MIN_TEST):
            N_test = len(test_precursors_B)
            for p in test_precursors_B:
            #for p in I_B:
                ori_B.append(x_B[p])
                ori_iRT.append(iRT[p])
                predict_from_ref.append(find_alignment_time(iRT[p], B_time))
            MAE = mean_absolute_error(predict_from_ref, ori_B)
            iRT_MAE = mean_absolute_error(ori_iRT, ori_B)
        else:            
            MAE = -1.0
            iRT_MAE = -1.0
        from_ref_MAE.append(MAE)
        fo.write('\t' + str(MAE) + '\t' + str(iRT_MAE) + '\n')
        
    return ori_MAE, from_A_MAE, from_ref_MAE


def main():
    args = arguments()
    args_evidence = args.evidence
    args_reference = args.reference
    args_output = args.output

    rows = read_evidence_file(args_evidence)
    rows, split_count = get_highest_intensity_features(rows)
    dict_B, inf = index_runs(rows)
    iRT = read_references(args_reference)

    
    name_set_A = list(dict_B.keys())
    random.shuffle(name_set_A)
    fo = open(args_output, 'w')
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
            
            # in case we care unmod features only, the feature will be ignored if one or more of the ID sequences contain mods
            c_peptide = peptide.replace("C+57.", "C")
            if ('+' in c_peptide) or ('-' in c_peptide):        
                continue
            
            if ';' not in peptide:      # should not be an ambiguous ID
                if (file_name == run_A):
                    row_A.append((RT, i))
        
        ori_MAE, from_A_MAE, from_ref_MAE = MAE_calculation(fo, fh, inf, iRT, run_A, row_A, dict_B, split_count)
        print("run_A =", run_A)    
        temp = [x for x in ori_MAE if (x > 0)]
        if len(temp) >= MIN_REPORT:
            print("ori_MAE =", mean(temp))
        else:
            print("ori_MAE =", -1.0)
            
        temp = [x for x in from_A_MAE if (x > 0)]
        if len(temp) >= MIN_REPORT:
            print("from_A_MAE =", mean(temp))
        else:
            print("from_A_MAE =", -1.0)
            
        temp = [x for x in from_ref_MAE if (x > 0)]
        if len(temp) >= MIN_REPORT:
            print("from_ref_MAE =", mean(temp))
        else:
            print("from_ref_MAE =", -1.0)
            
    fo.close()
    fh.close()


if __name__ == "__main__":
    main()
