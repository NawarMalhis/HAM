from lib import get_mask_arguments
# from lib.annotated_fasta import annotated_fasta_load, annotated_fasta_save
import os
import sys
from ham_param import *
if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *


if __name__ == '__main__':
    # target_file = 'TS2022p.af1'
    # data_path = 'data'
    arg = get_mask_arguments()
    target_file = arg.in_file
    data_path = arg.path
    results_path = f'{data_path}/results/'
    target_name = target_file.split('.')[0]
    af = annotated_fasta_load(f"{data_path}/{target_file}")
    tg = af['metadata']['tags'][0]
    hm_files = [x for x in os.listdir(results_path) if x[-3:] == '.af' and x.startswith(f"{target_name}-homology-to-")]
    for h_file in hm_files:
        h_name = h_file.split('-')[-1].split('.')[0]
        afh = annotated_fasta_load(f"{results_path}{h_file}")
        print(f"Masking {target_name} homology to {h_name}")
        for ac in af['data']:
            tg_list = list(af['data'][ac][tg])
            for ii in range(len(af['data'][ac]['seq'])):
                if afh['data'][ac]['H0'][ii] == '0' or afh['data'][ac]['H1'][ii] == '1':
                    tg_list[ii] = '-'
            af['data'][ac][tg] = ''.join(tg_list)
    annotated_fasta_save(af, f"{data_path}/ham-masked-{target_name}.af")
