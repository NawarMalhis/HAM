from lib import *


if __name__ == "__main__":
    arg = get_ham_arguments()
    data2 = load_data_2files(arg)
    compute_ham(arg, data2)
    hf_names = get_homology_file_names(names=list(data2['a_fasta'].keys()))
    for af_dta in data2['a_fasta']:
        _f_name = f"{arg.path}/results/{hf_names[af_dta]['out_file']}"
        _header_top = f"# Minimum identity: {arg.identity_cut_off}%\n# Minimum match: {arg.minimum_aligned_size}"
        _header_top = f"{_header_top} residues\n#"
        _header_bottom = ham_cross_stat(data2['a_fasta'][af_dta], f1=af_dta, f2=hf_names[af_dta]['to'])
        # annotated_fasta_remove_empty(data2['a_fasta'][af_dta], 'H1')
        data2['metadata']['data_name'] = af_dta
        aff_save2(af=data2['a_fasta'][af_dta], f_name=_f_name, header_top=_header_top, header_bottom=_header_bottom)
