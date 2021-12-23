# Function to display help text. Also acts as information how to run the script.


def f_display_help():
    print(
        '''
########################################################################################################################
#                                                                                                                      #
#  Function: filter_vcf                                                                                                #
#  This function parses a vcf file containing SNPs and filters SNPs according to parameters given. The function was    #
#  Designed to find allelic SNPs present in an in group sub population of samples, but absent in a out group           #
#  sub population. Groups where only a subset of samples have to be hit can also be defined.                           #
#                                                                                                                      #
#  Author: Pascal Schlaepfer                                                                                           #
#  Date: 2021/06/21                                                                                                    #
#  Version: 1.0.0                                                                                                      #
#                                                                                                                      #
########################################################################################################################

Usage:
python filter_vcf.py -in vcf_file_name [-argument argument_value]

Input:
 - Mandatory:
    -in [vcf_file]: Input vcf file that should be searched for a SNP variant of interest.

 - Optional:
    -p [path]: Path where input file is stored. If empty, script uses current working directory to look up input file.
    -out [file]: Output file name, where results should be stored (in tsv format). (default results_filtering_vcf.tsv)
    -rf [path]: Path where results should be stored. If not given, current working directory is used.
    -mtrc: Sets the parameter minimum absolute read coverage of a SNP in a sample of interest. (default is 20)
    -mnrc: Sets the parameter maximum absolute reads of a SNP in an out group sub population, so that if this number of 
           reads is found for the out group sample, it is assumed to be noise and the SNP is still flagged as
           interesting. (default is 2).
    -mrrc: Sets the parameter minimum relative read coverage of a variant of a SNP of interest in an in group sub
           population. (default is 0.33)
    -mq: Sets the parameter minimum quality a SNP has to reach to be flagged as interesting. (default is 100)
    -i: Defines the samples that should be counted as in group sub population. Sample ids have to be separated by ";".
    -ni: Defines the number of samples given in -i (can be lower than actual number).
    -o: Defines the samples that should be counted as out group sub population. Sample ids have to be separated by ";".
    -no: Defines the number of samples given in -o (can be lower than actual number).
    -fi: Defines the samples that should be counted as optional in group sub population. Sample ids have to be 
    separated by ";".
    -nfi: Defines the number of samples given in -fi (can be lower than actual number).
    -fo: Defines the samples that should be counted as optional out group sub population. Sample ids have to be 
         separated by ";".
    -nfo: Defines the number of samples given in -fo (can be lower than actual number).
    
 - Auxiliary:
   -h: Display this help message
   -v [True, False]: Set verbose to True or False and give screen output.
   -s TME204 or -s TME14 or -s 91-02324 Special mode to reproduce SNP findings of the associated manuscript.
        '''
    )


# Main function that performs filtering.


def filter_vcf(file_in, file_out, min_reads, in_group_allele_min_ratio, out_group_max_reads, min_snp_quality,
               in_group_samples, min_number_hits_in_group, out_group_samples, min_number_hits_out_group,
               optional_in_group_samples, min_number_hits_optional_in_group, optional_out_group_samples,
               min_number_hits_optional_out_group, print_verbose):
    # Print out arguments to screen
    if print_verbose:
        print('VCF file to read: ', file_in)
        print('Output file to be written to: ', file_out)
        print('Minimal required reads in samples of interest: ', min_reads)
        print('Minimal read ratio for a variant in ingroup(s): ', in_group_allele_min_ratio)
        print('Maximum number of reads for a variant in outgroup: ', out_group_max_reads)
        print('Minimum SNP quality: ', min_snp_quality)
        print('Sample ids for in group: ')
        print(in_group_samples)
        print('Minimum number of in group samples required: ', min_number_hits_in_group)
        print('Sample ids for out group: ')
        print(out_group_samples)
        print('Minimum number of out group samples required: ', min_number_hits_out_group)
        print('Sample ids for optional in group: ')
        print(optional_in_group_samples)
        print('Minimum number of optional in group samples required: ', min_number_hits_optional_in_group)
        print('Sample ids for optional out group: ')
        print(optional_out_group_samples)
        print('Minimum number of optional out groups samples required: ', min_number_hits_optional_out_group)

    # Headers for position, alternative alleles, format and quality
    chromosome_header = 'CHROM'
    position_header = 'POS'
    alternative_variant_header = 'ALT'
    format_header = 'FORMAT'
    dp_format_header = 'DP'
    ad_format_header = 'AD'
    quality_header = 'QUAL'

    # Test if vcf file can be opened
    if print_verbose:
        print('Test if files can be opened and written to')
    try:
        vfh = open(file_in, 'r')
        vfh.close()
    except IOError:
        print('File ' + file_in + ' can not be read.')
        sys.exit()

    # Try to write to output file
    try:
        vfh = open(file_out, 'w')
        vfh.close()
    except IOError:
        print('File ' + file_out + ' can not be written to.')
        sys.exit()

    # If tests are passed, open files
    fh_in = open(file_in, 'r')
    fh_out = open(file_out, 'w')

    # Prepare storing place for samples
    if print_verbose:
        print('Prepare storing space')
    in_group_samples_id = []
    for _ in in_group_samples:
        in_group_samples_id.append(0)

    out_group_samples_id = []
    for _ in out_group_samples:
        out_group_samples_id.append(0)

    optional_in_group_samples_id = []
    for _ in optional_in_group_samples:
        optional_in_group_samples_id.append(0)

    optional_out_group_samples_id = []
    for _ in optional_out_group_samples:
        optional_out_group_samples_id.append(0)

    # Initialize ids
    chromosome_header_id = -1
    position_header_id = -1
    alternative_variant_header_id = -1
    format_header_id = -1
    dp_format_header_id = -1
    ad_format_header_id = -1
    quality_header_id = -1
    quality = -1

    found_position = False

    # Go through file line by line.
    if print_verbose:
        print('Parse information part of vcf')
    for vcf_line_id, vcf_line in enumerate(fh_in):
        if vcf_line.startswith('##'):
            if print_verbose:
                print('Ignoring comment line starting with "##"')
        elif vcf_line.startswith('#'):  # Header line
            if print_verbose:
                print('Reading in header line')
            vcf_header = vcf_line.lstrip('#').rstrip('\n').split('\t')
            if print_verbose:
                print('Header:')
                print(vcf_header)
            fh_out.write(vcf_line.lstrip('#'))
            for vh, vcf_header_item in enumerate(vcf_header):
                if vcf_header_item == position_header:
                    position_header_id = vh
                elif vcf_header_item == chromosome_header:
                    chromosome_header_id = vh
                elif vcf_header_item == alternative_variant_header:
                    alternative_variant_header_id = vh
                elif vcf_header_item == format_header:
                    format_header_id = vh
                elif vcf_header_item == quality_header:
                    quality_header_id = vh
                else:
                    for vs, sample in enumerate(in_group_samples):
                        if vcf_header_item == sample:
                            in_group_samples_id[vs] = vh
                    for vs, sample in enumerate(out_group_samples):
                        if vcf_header_item == sample:
                            out_group_samples_id[vs] = vh
                    for vs, sample in enumerate(optional_in_group_samples):
                        if vcf_header_item == sample:
                            optional_in_group_samples_id[vs] = vh
                    for vs, sample in enumerate(optional_out_group_samples):
                        if vcf_header_item == sample:
                            optional_out_group_samples_id[vs] = vh
        else:
            vcf_content = vcf_line.lstrip('#').rstrip('\n').split('\t')
            n_all = 1
            all_complete = []
            all_abs_complete = []
            all_tot_complete = []
            position_header_id_to_show = ''
            chromosome_header_id_to_show = ''
            for vc, vcf_content_item in enumerate(vcf_content):
                if vc == format_header_id:
                    format_header_split = vcf_content_item.split(':')
                    for vf, format_item in enumerate(format_header_split):
                        if format_item == dp_format_header:
                            dp_format_header_id = vf
                        elif format_item == ad_format_header:
                            ad_format_header_id = vf
                if vc == alternative_variant_header_id:
                    format_header_split = vcf_content_item.split(',')
                    for _ in format_header_split:
                        n_all += 1
                if vc == chromosome_header_id:
                    chromosome_header_id_to_show = vcf_content_item
                if vc == position_header_id:
                    position_header_id_to_show = vcf_content_item
                if vc == quality_header_id:
                    quality = float(vcf_content_item)
            if print_verbose:
                if vcf_line_id % 1000 == 0:
                    print('At line ' + str(vcf_line_id) + ': ' + chromosome_header_id_to_show + ', ' +
                          position_header_id_to_show)

            # If quality is not met, skip
            if quality >= min_snp_quality:
                # Extract information per organism
                for id_samples in in_group_samples_id:
                    all_tot, all_ratios, all_abs = f_extract_info_variety(vcf_content, id_samples, n_all,
                                                                          dp_format_header_id, ad_format_header_id)
                    all_complete.append(all_ratios)
                    all_tot_complete.append(all_tot)
                    all_abs_complete.append(all_abs)
                for id_samples in out_group_samples_id:
                    all_tot, all_ratios, all_abs = f_extract_info_variety(vcf_content, id_samples, n_all,
                                                                          dp_format_header_id, ad_format_header_id)
                    all_complete.append(all_ratios)
                    all_tot_complete.append(all_tot)
                    all_abs_complete.append(all_abs)
                for id_samples in optional_in_group_samples_id:
                    all_tot, all_ratios, all_abs = f_extract_info_variety(vcf_content, id_samples, n_all,
                                                                          dp_format_header_id, ad_format_header_id)
                    all_complete.append(all_ratios)
                    all_tot_complete.append(all_tot)
                    all_abs_complete.append(all_abs)
                for id_samples in optional_out_group_samples_id:
                    all_tot, all_ratios, all_abs = f_extract_info_variety(vcf_content, id_samples, n_all,
                                                                          dp_format_header_id, ad_format_header_id)
                    all_complete.append(all_ratios)
                    all_tot_complete.append(all_tot)
                    all_abs_complete.append(all_abs)

                # Set up check
                is_a_true_hit = False

                # Check ingroup for minimal reads
                min_hits_reached = []
                count_in_group_min_reads = 0
                for vi_in in range(len(in_group_samples_id)):
                    if all_tot_complete[vi_in] >= min_reads:
                        min_hits_reached.append(True)
                        count_in_group_min_reads += 1
                    else:
                        min_hits_reached.append(False)

                # Check outgroup for minimal reads
                count_out_group_min_reads = 0
                for vi_out in range(len(out_group_samples_id)):
                    if all_tot_complete[vi_out] >= min_reads:
                        min_hits_reached.append(True)
                        count_out_group_min_reads += 1
                    else:
                        min_hits_reached.append(False)

                # Check optional ingroup for minimal reads
                count_opt_in_group_min_reads = 0
                for vi_opt_in in range(len(optional_in_group_samples_id)):
                    if all_tot_complete[vi_opt_in] >= min_reads:
                        min_hits_reached.append(True)
                        count_opt_in_group_min_reads += 1
                    else:
                        min_hits_reached.append(False)

                # Check outgroup for minimal reads
                count_opt_out_group_min_reads = 0
                for vi_opt_out in range(len(optional_out_group_samples_id)):
                    if all_tot_complete[vi_opt_out] >= min_reads:
                        min_hits_reached.append(True)
                        count_opt_out_group_min_reads += 1
                    else:
                        min_hits_reached.append(False)

                # If first check was successful, do second check
                if count_in_group_min_reads >= min_number_hits_in_group and \
                        count_out_group_min_reads >= min_number_hits_out_group and \
                        count_opt_in_group_min_reads >= min_number_hits_optional_in_group and \
                        count_opt_out_group_min_reads >= min_number_hits_optional_out_group:
                    all_complete_t = list(map(list, zip(*all_complete)))
                    all_abs_complete_t = list(map(list, zip(*all_abs_complete)))
                    for va, all_complete_t_line in enumerate(all_complete_t):
                        all_abs_complete_t_line = all_abs_complete_t[va]
                        count_in_group = 0
                        count_out_group = 0
                        count_opt_in_group = 0
                        count_opt_out_group = 0
                        for vai, all_complete_t_line_item in enumerate(all_complete_t_line):
                            if vai < len(in_group_samples_id) and \
                                    all_complete_t_line_item >= in_group_allele_min_ratio and \
                                    min_hits_reached[vai]:
                                count_in_group += 1
                            elif len(in_group_samples_id) <= vai < len(in_group_samples_id) + \
                                    len(out_group_samples_id) and \
                                    all_abs_complete_t_line[vai] <= out_group_max_reads and \
                                    min_hits_reached[vai]:
                                count_out_group += 1
                            elif len(in_group_samples_id) + len(out_group_samples_id) <= vai < \
                                    len(in_group_samples_id) + len(out_group_samples_id) + \
                                    len(optional_in_group_samples_id) and \
                                    all_complete_t_line_item >= in_group_allele_min_ratio and \
                                    min_hits_reached[vai]:
                                count_opt_in_group += 1
                            elif len(in_group_samples_id) + len(out_group_samples_id) + \
                                    len(optional_in_group_samples_id) <= vai < len(in_group_samples_id) + \
                                    len(out_group_samples_id) + len(optional_in_group_samples_id) + \
                                    len(optional_out_group_samples_id) and \
                                    all_abs_complete_t_line[vai] <= out_group_max_reads \
                                    and min_hits_reached[vai]:
                                count_opt_out_group += 1
                        if count_in_group >= min_number_hits_in_group and \
                                count_out_group >= min_number_hits_out_group and \
                                count_opt_in_group >= min_number_hits_optional_in_group and \
                                count_opt_out_group >= min_number_hits_optional_out_group:
                            is_a_true_hit = True
                if is_a_true_hit:
                    fh_out.write(vcf_line)
    fh_in.close()
    fh_out.close()


# Function to extract information from the SNP record of a Sample in a VCF file


def f_extract_info_variety(vcf_content, sample_column, n_all, dp_format_header_id, ad_format_header_id):
    all_ratios = []
    all_abs = []
    all_tot = 0
    for vc, vcf_content_item in enumerate(vcf_content):
        if vc == sample_column:
            for all_id in range(n_all):
                all_ratios.append(0.0)
                all_abs.append(0)
            if not vcf_content_item == '.':
                format_header_split = vcf_content_item.split(':')
                all_tot = int(format_header_split[dp_format_header_id])
                all_ratios = []
                all_abs = []
                for format_header_split_item in format_header_split[ad_format_header_id].split(','):
                    all_ratios.append(float(format_header_split_item) / float(all_tot))
                    all_abs.append(int(format_header_split_item))
    return all_tot, all_ratios, all_abs


# Entry point for script when  running from command line.
if __name__ == "__main__":
    # Test if sys is installed in python
    try:
        import sys
    except ImportError:
        print('Error, module sys is required.')
        exit()
        sys.exit()

    try:
        import os
    except ImportError:
        print('Error, module os is required.')
        sys.exit()

    # Default parameters:
    main_display_help = False
    main_arg = 1
    main_path = ''
    main_result_path = ''
    main_file = ''
    main_comb_in = [[]]
    main_n_comb_in = 0
    main_comb_out = [[]]
    main_n_comb_out = 0
    main_comb_opt_in = [[]]
    main_n_comb_opt_in = 0
    main_comb_opt_out = [[]]
    main_n_comb_opt_out = 0
    main_min_abs_reads = 20
    main_min_rel_reads_hit = 0.33
    main_max_noise_level = 2
    main_min_quality = 100
    main_output_name = ''
    display_help = False
    allow_change = True
    main_verbose = False
    while main_arg < len(sys.argv):
        if sys.argv[main_arg] == '-h':
            if allow_change:
                display_help = True
        elif sys.argv[main_arg] == '-mtrc':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_min_abs_reads = float(sys.argv[main_arg + 1])
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-mnrc':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_max_noise_level = float(sys.argv[main_arg + 1])
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-mrrc':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_min_rel_reads_hit = float(sys.argv[main_arg + 1])
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-mq':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_min_quality = float(sys.argv[main_arg + 1])
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-ni':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_n_comb_in = int(sys.argv[main_arg + 1])
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-no':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_n_comb_out = int(sys.argv[main_arg + 1])
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-nfi':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_n_comb_opt_in = int(sys.argv[main_arg + 1])
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-nfo':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_n_comb_opt_out = int(sys.argv[main_arg + 1])
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-i':
            if main_arg + 1 < len(sys.argv):
                main_comb_in_string = sys.argv[main_arg + 1]
                main_comb_in = [main_comb_in_string.split(';')]
            else:
                print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-o':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_comb_out_string = sys.argv[main_arg + 1]
                    main_comb_out = [main_comb_out_string.split(';')]
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-fi':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_comb_opt_in_string = sys.argv[main_arg + 1]
                    main_comb_opt_in = [main_comb_opt_in_string.split(';')]
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-fo':
            if allow_change:
                if main_arg + 1 < len(sys.argv):
                    main_comb_opt_out_string = sys.argv[main_arg + 1]
                    main_comb_opt_out = [main_comb_opt_out_string.split(';')]
                else:
                    print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                    sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-out':
            if main_arg + 1 < len(sys.argv):
                main_output_name = sys.argv[main_arg + 1]
            else:
                print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-rf':
            if main_arg + 1 < len(sys.argv):
                main_result_path = sys.argv[main_arg + 1]
            else:
                print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-in':
            if main_arg + 1 < len(sys.argv):
                main_file = sys.argv[main_arg + 1]
            else:
                print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-p':
            if main_arg + 1 < len(sys.argv):
                main_path = sys.argv[main_arg + 1]
            else:
                print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-v':
            if main_arg + 1 < len(sys.argv):
                if sys.argv[main_arg + 1] == 'True' or sys.argv[main_arg + 1] == 'true' or \
                        sys.argv[main_arg + 1] == '1' or sys.argv[main_arg + 1] == 't':
                    main_verbose = True
                elif sys.argv[main_arg + 1] == 'False' or sys.argv[main_arg + 1] == 'false' or \
                        sys.argv[main_arg + 1] == '0' or sys.argv[main_arg + 1] == 'f':
                    main_verbose = False
                else:
                    print('Unknown state of verbose: ' + sys.argv[main_arg + 1])
            else:
                print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                sys.exit()
            main_arg += 1
        elif sys.argv[main_arg] == '-s':
            allow_change = False
            if main_arg + 1 < len(sys.argv):
                if sys.argv[main_arg + 1] == 'TME204':
                    main_comb_in = ['TME204-F1-3594-R', 'TME204-F1-3639-R', 'TME204-F1-3641-R']
                    main_n_comb_in = 3
                    main_comb_out = ['TME204-F1-3591-S', 'TME204-F1-3596-S', 'TME204-F1-3604-S', 'TME204-F1-3610-S',
                                     'TME204-F1-3636-S']
                    main_n_comb_out = 5
                    main_comb_opt_in = []
                    main_n_comb_opt_in = 0
                    main_comb_opt_out = []
                    main_n_comb_opt_out = 0
                    main_min_abs_reads = 20
                    main_min_rel_reads_hit = 0.20
                    main_max_noise_level = 2
                    main_output_name = 'Results_for_manuscript_TME204.tsv'
                elif sys.argv[main_arg + 1] == 'TME14':
                    main_comb_in = ['TME14-3548-R', 'TME14-3555-R', 'TME14-3565-R', 'TME14-3558-R']
                    main_n_comb_in = 4
                    main_comb_out = ['TME14-3552-S', 'TME14-3557-S']
                    main_n_comb_out = 2
                    main_comb_opt_in = []
                    main_n_comb_opt_in = 0
                    main_comb_opt_out = ['60444-FEC-A-S3', '60444-FEC-B-C1', 'TME3-FEC-A-C2', 'TME3-FEC-B-C4',
                                         'TME7-FEC-A', 'TME7-FEC-B', 'TME7-FEC-C', 'TME8-OES-A', 'TME8-OES-B',
                                         'TME9-OES-A', 'TME9-OES-B', 'TME14-WT', 'TME204-FEC', 'TME204-F1-3591-S',
                                         'TME204-F1-3596-S', 'TME204-F1-3604-S', 'TME204-F1-3610-S',
                                         'TME204-F1-3636-S', 'TME419-FEC-A', 'TME419-FEC-B']
                    main_n_comb_opt_out = 11
                    main_min_abs_reads = 20
                    main_min_rel_reads_hit = 0.20
                    main_max_noise_level = 2
                    main_output_name = 'Results_for_manuscript_TME14.tsv'
                elif sys.argv[main_arg + 1] == '91-02324':
                    main_comb_in = ['91_02324-WT']
                    main_n_comb_in = 1
                    main_comb_out = ['TME14-3552-S', 'TME14-3557-S', 'TME204-F1-3591-S', 'TME204-F1-3596-S',
                                     'TME204-F1-3604-S', 'TME204-F1-3610-S', 'TME204-F1-3636-S', '60444-WT']
                    main_n_comb_out = 8
                    main_comb_opt_in = []
                    main_n_comb_opt_in = 0
                    main_comb_opt_out = ['60444-FEC-A-S3', '60444-FEC-B-C1', 'TME3-FEC-A-C2', 'TME3-FEC-B-C4',
                                         'TME7-FEC-A', 'TME7-FEC-B', 'TME7-FEC-C', 'TME8-OES-A', 'TME8-OES-B',
                                         'TME9-WT', 'TME9-OES-A', 'TME9-OES-B', 'TME14-WT', 'TME204-FEC',
                                         'TME419-FEC-A', 'TME419-FEC-B']
                    main_n_comb_opt_out = 14
                    main_min_abs_reads = 20
                    main_min_rel_reads_hit = 0.20
                    main_max_noise_level = 2
                    main_output_name = 'Results_for_manuscript_TMS9102324.tsv'
                else:
                    print('Option ' + sys.argv[main_arg] + ' ' + sys.argv[main_arg + 1] + ' unknown.')
                    sys.exit()
            else:
                print(sys.argv[main_arg] + ' needs to be followed by an argument value.')
                sys.exit()
            main_arg += 1
        else:
            print('Option ' + sys.argv[main_arg] + ' unknown.')
            sys.exit()
        main_arg += 1

    if display_help:
        f_display_help()
    else:
        missing_info = False
        if main_path == '':
            print('Path to vcf file missing (-p). Using current path ' + os.getcwd())
        if main_result_path == '':
            print('Path to store results file missing (-rf). Using current path ' + os.getcwd())
        if main_file == '':
            print('Mandatory input vcf file not given (-in).')
            missing_info = True
        if main_output_name == '':
            print('Output file name not given (-out). Using results_filtering_vcf.tsv')
            main_output_name = 'results_filtering_vcf.tsv'

        if missing_info:
            sys.exit()

        filter_vcf(main_path + main_file, main_result_path + main_output_name, main_min_abs_reads,
                   main_min_rel_reads_hit, main_max_noise_level, main_min_quality, main_comb_in, main_n_comb_in,
                   main_comb_out, main_n_comb_out, main_comb_opt_in, main_n_comb_opt_in, main_comb_opt_out,
                   main_n_comb_opt_out, main_verbose)
        print('Filtering performed successfully.\nResults can be found in ' + main_result_path + main_output_name)
