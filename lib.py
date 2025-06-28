import argparse
import numpy as np
import os
import time
import sys

from ham_param import *

if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *

def get_mask_arguments():
    st_ubc = 'Nawar Malhis (2024), the University of British Columbia'
    parser = argparse.ArgumentParser(description=f"Masking Homologous Annotations identified by ham.py. {st_ubc}.")
    parser.add_argument('-in', '--in_file', type=str, help='Input file in annotated fasta format',
                        required=True)
    parser.add_argument('-p', "--path", type=str,
                        help="Path for input, results, and output files", required=True)

    arguments = parser.parse_args()
    return arguments


def get_resolve_arguments():
    st_ubc = 'Nawar Malhis (2024), the University of British Columbia'
    parser = argparse.ArgumentParser(description=f"Resolving Annotation Conflict identified by hac.py. {st_ubc}.")
    parser.add_argument('-in', '--in_file', type=str, help='Input file in annotated fasta format',
                        required=True)
    parser.add_argument('-p', "--path", type=str,
                        help="Path for input, results, and output files", required=True)
    parser.add_argument('-pr', '--priority', type=str, choices=['10', '01', '-'], default='10',
                        help="Priority for resolving conflicts")
    arguments = parser.parse_args()
    return arguments



def get_hac_arguments():
    st_ubc = 'The University of British Columbia'
    xx = 'HAC (Homology Annotation Conflict) identifies conflict in annotation between homologous amino acid sequences.'
    parser = argparse.ArgumentParser(description=f"{xx} Nawar Malhis (2024), {st_ubc}.")
    parser.add_argument('-in', '--in_file', type=str, help='Input file in annotated fasta format',
                        required=True)
    parser.add_argument('-p', "--path", type=str, default='./',
                        help="Path for input & output files, default='./'")
    parser.add_argument('-ico', "--identity_cut_off", type=float, default=80,
                        help="Identity cut off, default: 80")
    parser.add_argument('-msz', "--minimum_aligned_size", type=float, default=10,
                        help="Minimum aligned region size, default: 10 amino acids")
    parser.add_argument('-num_threads', "--number_of_threads", type=float, default=8,
                        help="The number of threads used by BLASTP, default: 8")

    arguments = parser.parse_args()
    return arguments


def get_ham_arguments():
    st_name = 'By Nawar Malhis (2024), the University of British Columbia'
    seq_extra = 'identifies cross-annotation homologies between two sets of annotated amino acid sequences.'
    parser = argparse.ArgumentParser(description=f"HAM (Measuring Annotated Homology) {seq_extra}\n{st_name}")
    parser.add_argument('-in1', '--in_file1', type=str, help='Input file1 in annotated fasta format',
                        required=True)
    parser.add_argument('-in2', '--in_file2', type=str, help='Input file2 in annotated fasta format',
                        required=True)
    parser.add_argument('-p', "--path", type=str, default='./',
                        help="Path for input & output files, default='./'")
    parser.add_argument('-ico', "--identity_cut_off", type=float, default=80,
                        help="Identity cut off, default: 80")
    parser.add_argument('-msz', "--minimum_aligned_size", type=float, default=10,
                        help="Minimum aligned region size, default: 10 amino acids")
    parser.add_argument('-num_threads', "--number_of_threads", type=float, default=8,
                        help="The number of threads used by BLASTP, default: 8")
    arguments = parser.parse_args()
    return arguments


def load_data_1file(arg):
    # data = {}
    counts = [0, 0]
    pf1 = f"{arg.path}/{arg.in_file}"
    if1 = arg.in_file.split('.')[0]
    af = aff_load_simple(in_file=pf1)  # aff_load(in_file=pf1)

    tg = 'ANN'
    if len(af['metadata']['tags_dict']) > 0:
        tg = list(af['metadata']['tags_dict'].keys())[0]
    for ac in af['data']:
        for m in af['data'][ac][tg]:
            if m in ['0', '1']:
                counts[int(m)] += 1

    return af, counts, if1


def load_data_2files(arg):
    data2 = {'a_fasta': {}}
    pf1 = f"{arg.path}/{arg.in_file1}"
    pf2 = f"{arg.path}/{arg.in_file2}"
    if1 = arg.in_file1.split('.')[0]
    if2 = arg.in_file2.split('.')[0]
    data2['a_fasta'][if1] = aff_load_simple(in_file=pf1)
    data2['a_fasta'][if2] = aff_load_simple(in_file=pf2)
    for df in data2['a_fasta']:
        tg = list(data2['a_fasta'][df]['metadata']['tags_dict'].keys())[0]
        for ac in data2['a_fasta'][df]['data']:
            data2['a_fasta'][df]['data'][ac][tg] = data2['a_fasta'][df]['data'][ac][tg].replace('x', '-')
    add_data2_metadata(data2)
    return data2


def add_data2_metadata(data2):
    [if1, if2] = list(data2['a_fasta'].keys())
    if len(data2['a_fasta'][if1]['data']) > len(data2['a_fasta'][if2]['data']):
        db_id = if1
        sq_id = if2
    else:
        db_id = if2
        sq_id = if1
    data2['metadata'] = {'DB': db_id, 'SEQ': sq_id}
    return


def compute_ham(arg, data2):
    if not os.path.exists(f"{arg.path}/results/"):
        os.system(f"mkdir {arg.path}/results/")
    if not os.path.exists(f"tmp_ham"):
        os.system(f"mkdir tmp_ham")
    db_id = data2['metadata']['DB']
    sq_id = data2['metadata']['SEQ']
    db_data = data2['a_fasta'][db_id]['data']
    sq_data = data2['a_fasta'][sq_id]['data']
    make_db(db_data)
    details = open(f"{arg.path}/results/ham-details-{sq_id}-{db_id}.tsv", 'w')
    print(f"# identity cutoff:\t{arg.identity_cut_off}\n# minimum_aligned_size:\t{arg.minimum_aligned_size}",
          file=details)
    sq_tag = list(data2['a_fasta'][sq_id]['metadata']['tags_dict'].keys())[0]
    db_tag = list(data2['a_fasta'][db_id]['metadata']['tags_dict'].keys())[0]
    print(f"# AC({sq_id})\tsize\tpos\tAA\t{sq_tag}\tAC({db_id})\tsize\tpos\tAA\t{db_tag}\tidentity\tlength",
          file=details)
    for xx in ['0', '1', '-']:
        data2['a_fasta'][db_id]['metadata']['tags_dict'][f'H{xx}'] = f'{xx} annotation'
        data2['a_fasta'][sq_id]['metadata']['tags_dict'][f'H{xx}'] = f'{xx} annotation'
        for d_ac in db_data:
            db_data[d_ac][f'H{xx}'] = ['.'] * len(db_data[d_ac][db_tag])
        for q_ac in sq_data:
            sq_data[q_ac][f'H{xx}'] = ['.'] * len(sq_data[q_ac][sq_tag])

    for q_ac in sq_data:
        q_seq = sq_data[q_ac]['seq']
        q_mask = sq_data[q_ac][sq_tag]  # .replace('x', '-')
        q_size = len(q_mask)
        align_of6(q_ac, q_seq, num_threads=arg.number_of_threads)
        with open("tmp_ham/tmp_out.tsv", 'r') as fin:
            for line in fin:
                lst = line.strip().split()
                identity = float(lst[2])
                sz = int(lst[3])
                g_open = int(lst[5])
                if identity > arg.identity_cut_off and g_open == 0 and sz >= arg.minimum_aligned_size:
                    d_ac = lst[1]
                    if d_ac in db_data:
                        d_mask = db_data[d_ac][db_tag]  # .replace('x', '-')
                        d_size = len(d_mask)
                        d_seq = db_data[d_ac]['seq']
                        d_st = int(lst[8]) - 1
                        q_st = int(lst[6]) - 1
                        for i in range(sz):
                            qi = i + q_st
                            di = i + d_st
                            dm = d_mask[di]
                            qm = q_mask[qi]
                            if dm not in ['0', '1', '-'] or qm not in ['0', '1', '-']:
                                print(f"{d_ac}\t{di}\t{dm}\t{q_ac}\t{qi}\t{qm}")
                                continue
                            sq_data[q_ac][f'H{dm}'][qi] = dm
                            db_data[d_ac][f'H{qm}'][di] = qm
                            print(f"{q_ac}\t{q_size:,}\t{qi + 1:,}\t{q_seq[qi]}\t{qm}\t", end='', file=details)
                            print(f"{d_ac}\t{d_size:,}\t{di + 1:,}\t{d_seq[di]}\t{dm}\t", end='', file=details)
                            print(f"{identity:.2f}\t{sz:,}", file=details)

    for d_ac in db_data:
        for xx in ['0', '1', '-']:
            db_data[d_ac][f'H{xx}'] = ''.join(db_data[d_ac][f'H{xx}'])  # ['0'] * len(db_data[d_ac][db_tag])
    for q_ac in sq_data:
        for xx in ['0', '1', '-']:
            sq_data[q_ac][f'H{xx}'] = ''.join(sq_data[q_ac][f'H{xx}'])  # ['0'] * len(sq_data[q_ac][sq_tag])

    details.close()
    os.system("rm -r tmp_ham")
    return


def make_db(db_data):
    if not os.path.exists('tmp_ham'):
        os.system("mkdir tmp_ham")
    if not os.path.exists('tmp_ham/DB'):
        os.system("mkdir tmp_ham/DB")
    f_db = 'tmp_ham/db.fasta'
    db = 'tmp_ham/DB/db'
    with open(f_db, 'w') as fout:
        for ac in db_data:
            print(f">{ac}\n{db_data[ac]['seq']}", file=fout)
    cmd = f"makeblastdb -in {f_db} -out {db} -parse_seqids -dbtype prot > tmp_ham/junk"
    os.system(cmd)


def align_of6(ac, seq, num_threads=8):
    with open('tmp_ham/tmp_in.fasta', 'w') as f_out:
        print(f">{ac}\n{seq}", file=f_out)
    cmd = f"blastp -db tmp_ham/DB/db -query tmp_ham/tmp_in.fasta"
    cmd = f"{cmd} -evalue 1000 -outfmt 6 -word_size 3 -num_threads {int(num_threads)} -out tmp_ham/tmp_out.tsv"
    # print(cmd, flush=True)
    os.system(cmd)
    os.system("rm tmp_ham/tmp_in.fasta")


def get_homology_file_names(names):
    _hf_names = {}
    for fn in names:
        _hf_names[fn] = {}
        for fn2 in names:
            if fn2 != fn:
                _hf_names[fn]['to'] = fn2
                _hf_names[fn]['out_file'] = f"{fn}-homology-to-{fn2}.af"
                break
    return _hf_names


def hac_cross_stat(af, counts, d_name=''):
    sz = len(af['data'])
    # counts = [0, 0]
    tg = list(af['metadata']['tags_dict'])[0]
    homology = np.zeros((2, 2), dtype='int32')
    for ac in af['data']:
        mask = af['data'][ac][tg]
        # cnt[0] += mask.count('0')
        # cnt[1] += mask.count('1')
        for i in range(len(mask)):
            if mask[i] in ['0', '1']:
                qi = int(mask[i])
                for di in [0, 1]:
                    if af['data'][ac][f'H{di}'][i] == str(di):
                        homology[qi][di] += 1
    ret = '# Annotation conflict summery:\n'
    ret = f"{ret}#   {d_name} '0' count:\t{counts[0]:,}.\n"
    ret = f"{ret}#   {d_name} '1' count:\t{counts[1]:,}.\n#\n"
    ret = f"{ret}#   The number of {d_name} '0' residues homologous to at least one {d_name} '0' is"
    ret = f"{ret} {homology[0][0]:,}, which is {homology[0][0]/counts[0]:.2%} of {d_name} '0'.\n"
    ret = f"{ret}#   The number of {d_name} '0' residues homologous to at least one {d_name} '1' is"
    ret = f"{ret} {homology[0][1]:,}, which is {homology[0][1]/counts[0]:.2%} of {d_name} '0'.\n"
    ret = f"{ret}#   The number of {d_name} '1' residues homologous to at least one {d_name} '0' is"
    ret = f"{ret} {homology[1][0]:,}, which is {homology[1][0]/counts[1]:.2%} of {d_name} '1'.\n"
    ret = f"{ret}#   The number of {d_name} '1' residues homologous to at least one {d_name} '1' is"
    ret = f"{ret} {homology[1][1]:,}, which is {homology[1][1]/counts[1]:.2%} of {d_name} '1'.\n#"
    return ret


def ham_cross_stat(af, f1, f2):
    # if not af['metadata']['counts']:
    aff_gen_counts(af)
    # print(list(af['metadata']['counts']), flush=True)
    # exit(0)
    ky_list = list(af['metadata']['tags_dict'].keys())
    counts = np.zeros((3, 3), dtype='int32')
    for ac in af['data']:
        for m, h0, h1, h_ in zip(af['data'][ac][ky_list[0]], af['data'][ac]['H0'],
                                 af['data'][ac]['H1'], af['data'][ac]['H-']):
            if m == '-':
                im = 2
            else:
                im = int(m)
            if h0 == '0':
                counts[im][0] += 1
            if h1 == '1':
                counts[im][1] += 1
            if h_ == '-':
                counts[im][2] += 1

    m_totals = np.array([af['metadata']['counts']['tags_dict'][ky_list[0]]['0'],
                         af['metadata']['counts']['tags_dict'][ky_list[0]]['1'],
                         af['metadata']['counts']['tags_dict'][ky_list[0]]['-']], dtype='int32')
    tg = ky_list[0]  # af['metadata']['tags_dict'][0]
    total = m_totals[0:2].sum()
    h_total = counts[0:2, 0:2].sum()
    h0_disc = f"‘0’ for residues that are homologous to ‘0’ annotated residues in {f2}, otherwise, it is ‘.’."
    h1_disc = f"‘1’ for residues that are homologous to ‘1’ annotated residues in {f2}, otherwise, it is ‘.’."
    h__disc = f"‘-’ for residues that are homologous to ‘-’ annotated residues in {f2}, otherwise, it is ‘.’."
    af['metadata']['tags_dict']['H0'] = h0_disc
    af['metadata']['tags_dict']['H1'] = h1_disc
    af['metadata']['tags_dict']['H-'] = h__disc
    ret = f"# HAM Cross-annotations\n# ------------\t-------\t{f1}(0)\t{f1}(0)\t{f1}(1)\t{f1}(1)\n"
    ret = ret + f"# ------------\th_total\t{f2}(0)\t{f2}(1)\t{f2}(0)\t{f2}(1)\n"
    ret = ret + f"# Counts ....:\t{h_total:,}\t{counts[0][0]:,}\t{counts[0][1]:,}\t{counts[1][0]:,}\t{counts[1][1]:,}\n"
    ret = ret + f"# Percentages:\t{float(h_total)/total:.2%}\t{float(counts[0][0]) / m_totals[0]:.2%}\t"
    ret = ret + f"{float(counts[0][1]) / m_totals[0]:.2%}\t{float(counts[1][0]) / m_totals[1]:.2%}\t"
    ret = ret + f"{float(counts[1][1]) / m_totals[1]:.2%}\n#\n"
    ret = ret + f"# The cross-annotation 'Counts' provides the counts of {f1} residues homologous to {f2} for all\n"
    ret = ret + f"#    possible '0' and '1' annotations in {f1} and {f2}.\n#\n"
    ret = ret + f"# The cross-annotation 'Percentages' provides the percentages of each {f1} class residues\n"
    ret = ret + f"#    homologous to {f2} for all possible '0' and '1' annotations in {f1} and {f2}.\n#\n"
    ret = ret + "# Examples:\n"
    c0 = af['metadata']['counts']['tags_dict'][tg]['0']
    ret = ret + f"#  1) {h_total/total:.2%} ({h_total:,}) of the {f1} total annotated residues ({total:,}) are "
    ret = ret + f"homologous to {f2}.\n"
    ret = ret + f"#  2) {counts[0][0]/c0:.2%} ({counts[0][0]:,}) of the {f1} '0' annotated residues ({c0:,}) are "
    ret = ret + f"homologous to {f2} '0' annotated residues."
    return ret


def compute_hac(arg, af, d_name):
    cnt = 0
    path = arg.path
    identity_cut = arg.identity_cut_off
    msz = arg.minimum_aligned_size  # minimum_aligned_size
    num_threads = arg.number_of_threads
    if not os.path.exists(f"{path}/results/"):
        os.system(f"mkdir {path}/results/")
    make_db(af['data'])
    tg = list(af['metadata']['tags_dict'])[0]
    for ac in af['data']:
        mask_sz = af['data'][ac][tg]
        af['data'][ac]['H0'] = ['.'] * len(mask_sz)
        af['data'][ac]['H1'] = ['.'] * len(mask_sz)

    details = open(f"{path}/results/hac-details-{d_name}.tsv", 'w')
    print('# P1: Protein 1\n# P2: Protein 2\n# P1_AC\tP1_index\tP1_AA\tP1_class', end='\t', file=details)
    print('P2_AC\tP2_index\tP2_AA\tP2_class\tIdentity', file=details)
    start_time = time.time()
    # test_list = ['PDB:1a3b_I', 'PDB:1a4w_I', 'PDB:1d4p_H', 'PDB:1fpc_I', 'PDB:1ghv_I', 'PDB:1lhf_I', 'PDB:1vit_I']
    # test_list = ['PDB:1a3b_I']
    # test_list = ['PDB:1a4w_I']
    for q_ac in af['data']:
        # print(q_ac, flush=True)
        q_seq = af['data'][q_ac]['seq']
        q_mask = af['data'][q_ac][tg]
        align_of6(q_ac, q_seq, num_threads=num_threads)
        # exit(0)
        with open("tmp_ham/tmp_out.tsv", 'r') as fin:
            for line in fin:
                line = line.strip()
                lst = line.split()
                d_ac = lst[1]
                if q_ac != d_ac:
                    identity = float(lst[2])
                    sz = int(lst[3])
                    g_open = int(lst[5])
                    if identity > identity_cut and g_open == 0 and sz >= msz:
                        d_mask = af['data'][d_ac][tg]
                        d_seq = af['data'][d_ac]['seq']
                        d_st = int(lst[8]) - 1
                        q_st = int(lst[6]) - 1
                        for i in range(sz):
                            qi = i + q_st
                            di = i + d_st
                            if d_mask[di] in ['0', '1'] and q_mask[qi] in ['0', '1']:
                                dmi = int(d_mask[di])
                                qmi = int(q_mask[qi])
                                af['data'][q_ac][f'H{dmi}'][qi] = d_mask[di]
                                af['data'][d_ac][f'H{qmi}'][di] = q_mask[qi]
                                print(f"{q_ac}\t{qi+1}\t{q_seq[qi]}\t{q_mask[qi]}\t", end='', file=details)
                                print(f"{d_ac}\t{di+1}\t{d_seq[di]}\t{d_mask[di]}\t{identity}", file=details)
        cnt += 1
        if cnt % 100 == 0:
            print(f"{float(cnt) / len(af['data']):.0%} completed in {float(time.time() - start_time):5,.0f} seconds",
                  flush=True)
    details.close()
    tmp = "‘0’ for residues that are homologous to ‘0’ annotated residues, otherwise, it is ‘.’."
    af['metadata']['tags_dict']['H0'] = tmp
    tmp = "‘1’ for residues that are homologous to ‘1’ annotated residues, otherwise, it is ‘.’."
    af['metadata']['tags_dict']['H1'] = tmp
    for ac in af['data']:
        af['data'][ac]['H0'] = ''.join(af['data'][ac]['H0'])
        af['data'][ac]['H1'] = ''.join(af['data'][ac]['H1'])
    os.system("rm -r tmp_ham")
    process_time = f"completed in {float(time.time() - start_time):,.0f} seconds"
    print(f"100% {process_time}", flush=True)
    return process_time


def save_hac(arg, af, d_name, process_time, homology, counts):
    path = arg.path
    iden_cut = arg.identity_cut_off
    msz = arg.minimum_aligned_size
    h_file = f"{path}/results/{d_name}_annotation_conflict.af"
    with open(h_file, 'w') as fout:
        print_hac_header(d_name=d_name, fout=fout, process_time=process_time, seq_count=len(af['data']), counts=counts,
                         homology=homology, cut=iden_cut, msz=msz)
        for ac in af['data']:
            h0 = ''.join(af['data'][ac]['H'][0])
            h1 = ''.join(af['data'][ac]['H'][1])
            if '0' in h0 or '1' in h1:
                print(f">{ac}", file=fout)
                print(f"{af['data'][ac]['seq']}", file=fout)
                print(f"{af['data'][ac]['mask']}", file=fout)
                print(f"{h0}", file=fout)
                print(f"{h1}", file=fout)


def print_hac_header(d_name, fout, process_time, seq_count, counts, homology, cut, msz):
    ht = 'residues homologous to at least one'
    print(f"# {d_name} annotation conflict at {cut}% identity cut off", end='', file=fout)
    print(f" and minimum aligned size of {msz} amino acids\n#", file=fout)
    print('# Annotation conflict summery:\n#', file=fout)
    print(f'# {d_name} includes {seq_count:,} sequences {process_time}', file=fout)
    print(f"# {d_name}[0] count:\t{counts[0]:,}\n# {d_name}[1] count:\t{counts[1]:,}", file=fout)

    x = homology[0][0]
    xp = x / counts[0]
    print('#', file=fout)
    print(f'# the number of {d_name}[0] {ht} {d_name}[0] is {x:,}, which is {xp:.2%} of {d_name}[0]', file=fout)
    x = homology[0][1]
    xp = x / counts[0]
    print(f'# the number of {d_name}[0] {ht} {d_name}[1] is {x:,}, which is {xp:.2%} of {d_name}[0]', file=fout)
    x = homology[1][0]
    xp = x / counts[1]
    print(f'# the number of {d_name}[1] {ht} {d_name}[0] is {x:,}, which is {xp:.2%} of {d_name}[1]', file=fout)
    x = homology[1][1]
    xp = x / counts[1]
    print(f'# the number of {d_name}[1] {ht} {d_name}[1] is {x:,}, which is {xp:.2%} of {d_name}[1]', file=fout)
    print("#", file=fout)
    print(f'# Format:\n# 1. >ID\n# 2. Amino acid sequence\n# 3. Annotation', file=fout)
    ex_st = 'stands for no aligned sequences with the target annotation'
    print(f"# 4. Homologous to '0'\n# 5. Homologous to '1'\n#\n# '.' {ex_st}\n#", file=fout)


def get_hac_top(d_name, sz, min_identity, min_size):
    ret = f"# Homology Annotation Conflict (HAC) for {d_name}, {sz:,} sequences."
    ret = f"{ret}\n#\tHomology minimum identity: {min_identity}%\n#\tHomology minimum length: {min_size} residues\n#"
    return ret


def hac_resolve_conflict(afc, priority=None):
    tg = list(afc['metadata']['tags_dict'])[0]
    del afc['metadata']['tags_dict']['H0']
    del afc['metadata']['tags_dict']['H1']
    for ac in afc['data']:
        mask_list = list(afc['data'][ac][tg])
        for ii in range(len(mask_list)):
            if priority != '-':
                if mask_list[ii] == priority[0]:
                    continue
                m_set = {afc['data'][ac][tg][ii], afc['data'][ac]['H1'][ii], afc['data'][ac]['H0'][ii]}
                for pp in priority:
                    if pp in m_set:
                        mask_list[ii] = pp
                        break
            else:
                m_set = {afc['data'][ac][tg][ii], afc['data'][ac]['H1'][ii], afc['data'][ac]['H0'][ii]}
                if '0' in m_set and '1' in m_set:
                    mask_list[ii] = '-'
        afc['data'][ac][tg] = ''.join(mask_list)

