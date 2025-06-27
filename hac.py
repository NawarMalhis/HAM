from lib import *


# https://www.cell.com/current-biology/home
if __name__ == '__main__':
    # _p = 'data'
    arg = get_hac_arguments()
    af, counts, d_name = load_data_1file(arg)
    process_time = compute_hac(arg, af, d_name)
    sz = len(af['data'])
    _header_top = get_hac_top(d_name, sz, min_identity=arg.identity_cut_off, min_size=arg.minimum_aligned_size)
    header_bottom = hac_cross_stat(af, counts=counts, d_name=d_name)
    f_name = f"{arg.path}/results/{d_name}-annotation-conflict.af"
    aff_save2(af=af, f_name=f_name, header_top=_header_top, header_bottom=header_bottom)


    # print(homology_2x2)
    # save_hac(arg, data, d_name, process_time, homology, counts)
