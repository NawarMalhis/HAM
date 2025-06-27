# from annotated_fasta import annotated_fasta_load, annotated_fasta_save
from lib import hac_resolve_conflict, get_resolve_arguments
import os
import sys
from ham_param import *

if aff_path not in sys.path:
    sys.path.append(aff_path)

from annotated_fasta import *


if __name__ == '__main__':
    arg = get_resolve_arguments()
    target_file = arg.in_file
    data_path = arg.path
    priority = arg.priority
    # target_file = arg.in_file
    # data_path = arg.path
    # priority = arg.priority
    results_path = f'{data_path}/results/'
    target_name = target_file.split('.')[0]
    # af = annotated_fasta_load(f"{data_path}/{target_file}")
    afc = aff_load2(f"{results_path}{target_name}-annotation-conflict.af")
    hac_resolve_conflict(afc, priority=priority)
    ex = priority
    if priority == '-':
        ex = 'masked'
    aff_save2(afc, f"{data_path}/{target_name}-resolved-{ex}.af")
