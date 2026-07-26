#!/usr/bin/env python3
import os
import sys
import argparse

def parse_info(info_str):
    info_dict = {}
    for item in info_str.split(';'):
        if '=' in item:
            k, v = item.split('=', 1)
            info_dict[k] = v
        else:
            info_dict[item] = True
    return info_dict

def get_format_dict(format_keys, sample_str):
    keys = format_keys.split(':')
    vals = sample_str.split(':')
    return dict(zip(keys, vals))

def parse_num_list(val_str, mode='max'):
    if not val_str or val_str == '.':
        return None
    try:
        vals = [float(x) for x in val_str.replace('|', ',').split(',') if x != '.']
        if not vals:
            return None
        return max(vals) if mode == 'max' else min(vals)
    except ValueError:
        return None

def classify_variant_type(ref, alt_str):
    """
    Classifies a variant record as 'snv', 'indel', or 'mnv'.
    Handles multi-allelic ALTs (comma-separated).
    """
    alts = alt_str.split(',')
    types = set()
    ref_len = len(ref)

    for alt in alts:
        alt_len = len(alt)
        if ref_len == 1 and alt_len == 1:
            types.add('snv')
        elif ref_len != alt_len:
            types.add('indel')
        else:
            types.add('mnv')

    # Priority return for filtering selection
    if 'indel' in types:
        return 'indel'
    elif 'mnv' in types:
        return 'mnv'
    return 'snv'

def filter_vcf_record(fields, args):
    ref = fields[3]
    alt = fields[4]
    var_type = classify_variant_type(ref, alt)

    # Variant Type Selection Filter
    if args.var_type and args.var_type != 'all':
        if args.var_type == 'snv' and var_type != 'snv':
            return False, var_type
        elif args.var_type == 'indel' and var_type != 'indel':
            return False, var_type
        elif args.var_type == 'mnv' and var_type != 'mnv':
            return False, var_type

    chrom = fields[0]
    qual_str = fields[5]
    filter_col = fields[6]
    info_str = fields[7]
    info = parse_info(info_str)

    # 1. Chromosome & Standard Checks
    if args.chrom and chrom != args.chrom:
        return False, var_type

    if args.pass_only and filter_col != 'PASS':
        return False, var_type

    # 2. QUAL Filter (Column 6)
    if args.min_qual is not None:
        try:
            if float(qual_str) < args.min_qual: return False, var_type
        except ValueError:
            return False, var_type

    # 3. Mapping Quality Filters (MQ / MMQ / MAQ)
    if args.min_mq is not None:
        mq_val = parse_num_list(info.get('MQ'), 'min')
        if mq_val is None or mq_val < args.min_mq: return False, var_type

    if args.min_mmq is not None:
        mmq_val = parse_num_list(info.get('MMQ'), 'min')
        if mmq_val is None or mmq_val < args.min_mmq: return False, var_type

    # 4. Total Depth Filters (INFO/DP)
    if args.min_info_dp is not None:
        if 'DP' not in info or float(info['DP']) < args.min_info_dp: return False, var_type

    # 5. Mutect2 INFO Field Filters
    if args.min_tlod is not None:
        val = parse_num_list(info.get('TLOD'), 'max')
        if val is None or val < args.min_tlod: return False, var_type

    if args.min_nlod is not None:
        val = parse_num_list(info.get('NLOD'), 'max')
        if val is None or val < args.min_nlod: return False, var_type

    if args.min_nalod is not None:
        val = parse_num_list(info.get('NALOD'), 'max')
        if val is None or val < args.min_nalod: return False, var_type

    if args.max_ecnt is not None:
        if 'ECNT' not in info or int(info['ECNT']) > args.max_ecnt: return False, var_type

    if args.min_mbq is not None:
        val = parse_num_list(info.get('MBQ'), 'min')
        if val is None or val < args.min_mbq: return False, var_type

    if args.min_mfrl is not None:
        val = parse_num_list(info.get('MFRL'), 'min')
        if val is None or val < args.min_mfrl: return False, var_type

    if args.min_mpos is not None:
        val = parse_num_list(info.get('MPOS'), 'min')
        if val is None or val < args.min_mpos: return False, var_type

    if args.exclude_str and 'STR' in info:
        return False, var_type

    # 6. FORMAT / Sample Level Filters (Tumor Sample - Col 10)
    if len(fields) >= 10:
        fmt = get_format_dict(fields[8], fields[9])

        if args.min_sample_dp is not None:
            val = parse_num_list(fmt.get('DP'), 'max')
            if val is None or val < args.min_sample_dp: return False, var_type

        if args.min_af is not None:
            val = parse_num_list(fmt.get('AF'), 'max')
            if val is None or val < args.min_af: return False, var_type

        if args.min_alt_ad is not None:
            ad_str = fmt.get('AD')
            if not ad_str: return False, var_type
            ad_vals = [int(x) for x in ad_str.split(',') if x != '.']
            alt_reads = max(ad_vals[1:]) if len(ad_vals) > 1 else 0
            if alt_reads < args.min_alt_ad: return False, var_type

        if args.min_f1r2 is not None:
            val = parse_num_list(fmt.get('F1R2'), 'max')
            if val is None or val < args.min_f1r2: return False, var_type

        if args.min_f2r1 is not None:
            val = parse_num_list(fmt.get('F2R1'), 'max')
            if val is None or val < args.min_f2r1: return False, var_type

        if args.require_phased and 'PGT' not in fmt and '|' not in fmt.get('GT', ''):
            return False, var_type

    return True, var_type

def run_pipeline(args):
    input_path = args.input
    if not os.path.isfile(input_path):
        print(f"Error: File '{input_path}' not found.", file=sys.stderr)
        sys.exit(1)

    input_dir, file_name = os.path.split(input_path)
    base_name, ext = os.path.splitext(file_name)
    if ext == ".gz":
        base_name, ext2 = os.path.splitext(base_name)
        ext = ext2 + ext

    # Setup Out Files (Split Mode vs Single File Mode)
    split_mode = args.split_types
    writers = {}
    header_lines = []

    print("\n" + "=" * 65)
    print("      COMPREHENSIVE SOMATIC VCF FILTERING PIPELINE      ")
    print("=" * 65)
    print(f"Input VCF: {input_path}")

    # Read Header
    with open(input_path, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                header_lines.append(line)
            else:
                break

    if split_mode:
        snv_path = os.path.join(input_dir, f"{base_name}_snvs{ext}")
        indel_path = os.path.join(input_dir, f"{base_name}_indels{ext}")
        mnv_path = os.path.join(input_dir, f"{base_name}_mnvs{ext}")
        
        writers['snv'] = open(snv_path, 'w')
        writers['indel'] = open(indel_path, 'w')
        writers['mnv'] = open(mnv_path, 'w')
        
        for w in writers.values():
            w.writelines(header_lines)
            
        print(f"Mode:      Splitting into separate SNV/INDEL files")
        print(f"Outputs:   \n - {snv_path}\n - {indel_path}\n - {mnv_path}\n")
    else:
        if args.output:
            output_path = args.output if os.path.isabs(args.output) else os.path.join(input_dir, args.output)
        else:
            tag = args.tag if args.tag else "filtered"
            type_suffix = f"_{args.var_type}" if args.var_type and args.var_type != 'all' else ""
            output_path = os.path.join(input_dir, f"{base_name}_{tag}{type_suffix}{ext}")

        out_fh = open(output_path, 'w')
        out_fh.writelines(header_lines)
        print(f"Output VCF: {output_path}\n")

    total, passed = 0, 0
    counts = {'snv': 0, 'indel': 0, 'mnv': 0}

    with open(input_path, 'r') as infile:
        for line in infile:
            if line.startswith('#'):
                continue
            total += 1
            fields = line.strip().split('\t')
            keep, vtype = filter_vcf_record(fields, args)

            if keep:
                passed += 1
                counts[vtype] += 1
                if split_mode:
                    writers[vtype].write(line)
                else:
                    out_fh.write(line)

    # Close handles
    if split_mode:
        for w in writers.values():
            w.close()
    else:
        out_fh.close()

    print(f"Total Variants Evaluated: {total}")
    print(f"Variants Retained:        {passed}")
    print(f"  └─ SNVs Passed:         {counts['snv']}")
    print(f"  └─ Indels Passed:       {counts['indel']}")
    print(f"  └─ MNVs Passed:         {counts['mnv']}")
    print("=" * 65 + "\n")

def interactive_mode():
    print("\n--- Interactive Filter Setup ---")
    input_path = input("Path to Input VCF: ").strip()

    print("\nVariant Type Filtering:")
    print(" 1. All Variants (Keep SNVs + Indels + MNVs)")
    print(" 2. Only SNVs")
    print(" 3. Only Indels")
    print(" 4. Split into separate files (_snvs.vcf and _indels.vcf)")
    vtype_choice = input("Select Option (1-4) [default: 1]: ").strip()

    var_type = 'all'
    split_types = False
    if vtype_choice == '2': var_type = 'snv'
    elif vtype_choice == '3': var_type = 'indel'
    elif vtype_choice == '4': split_types = True

    out_name = None
    tag = "filtered"
    if not split_types:
        out_name = input("\nCustom Output File Name [ENTER to use suffix tag]: ").strip() or None
        if not out_name:
            tag = input("Output Tag/Suffix [default: 'filtered']: ").strip() or "filtered"

    def prompt(label, val_type=float):
        res = input(f"{label}: ").strip()
        return val_type(res) if res else None

    args = argparse.Namespace(
        input=input_path,
        output=out_name,
        tag=tag,
        var_type=var_type,
        split_types=split_types,
        chrom=input("Chrom (e.g. chr2) [ENTER to skip]: ").strip() or None,
        pass_only=input("Keep ONLY 'PASS' variants? (y/N): ").lower().startswith('y'),
        min_qual=prompt("Min QUAL Score (Col 6) [ENTER to skip]"),
        min_mq=prompt("Min Standard Mapping Quality (MQ) [e.g. 50]"),
        min_mmq=prompt("Min Median Mapping Quality (MMQ) [e.g. 50]"),
        min_mbq=prompt("Min Median Base Quality (MBQ) [e.g. 30]"),
        min_info_dp=prompt("Min Total INFO Depth (DP) [e.g. 15]", int),
        min_sample_dp=prompt("Min Tumor Sample Depth (FORMAT DP) [e.g. 10]", int),
        min_tlod=prompt("Min Tumor LOD (TLOD) [e.g. 6.0]"),
        min_nlod=prompt("Min Normal LOD (NLOD) [e.g. 3.0]"),
        min_nalod=prompt("Min Normal Alt LOD (NALOD)"),
        max_ecnt=prompt("Max Event Count (ECNT) [e.g. 3]", int),
        min_mfrl=prompt("Min Fragment Length (MFRL)"),
        min_mpos=prompt("Min Distance from Read End (MPOS) [e.g. 10]"),
        exclude_str=input("Exclude Short Tandem Repeats (STR)? (y/N): ").lower().startswith('y'),
        min_af=prompt("Min Tumor Allele Frequency (AF) [e.g. 0.10]"),
        min_alt_ad=prompt("Min Alt Read Count (AD alt) [e.g. 3]", int),
        min_f1r2=prompt("Min F1R2 Read Count", int),
        min_f2r1=prompt("Min F2R1 Read Count", int),
        require_phased=input("Keep ONLY Phased Variants (PGT/PID)? (y/N): ").lower().startswith('y')
    )
    run_pipeline(args)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Full VCF parameter filtering & SNV/Indel separator tool.")
    parser.add_argument("-i", "--input", required=False, help="Input VCF file path")
    parser.add_argument("-o", "--output", help="Exact output filename")
    parser.add_argument("-t", "--tag", help="Output suffix tag (e.g., 'filtered')")
    parser.add_argument("--var-type", choices=['all', 'snv', 'indel', 'mnv'], default='all', help="Filter for specific variant type")
    parser.add_argument("--split-types", action="store_true", help="Split output into separate _snvs.vcf, _indels.vcf, and _mnvs.vcf files")
    parser.add_argument("--chrom", help="Chromosome filter")
    parser.add_argument("--pass-only", action="store_true", help="Keep PASS only")
    parser.add_argument("--min-qual", type=float, help="Min QUAL score")
    parser.add_argument("--min-mq", type=float, help="Min Mapping Quality (INFO/MQ)")
    parser.add_argument("--min-mmq", type=float, help="Min Median Mapping Quality (INFO/MMQ)")
    parser.add_argument("--min-mbq", type=float, help="Min Median Base Quality (INFO/MBQ)")
    parser.add_argument("--min-info-dp", type=float, help="Min total INFO depth")
    parser.add_argument("--min-sample-dp", type=int, help="Min sample FORMAT depth")
    parser.add_argument("--min-tlod", type=float, help="Min TLOD")
    parser.add_argument("--min-nlod", type=float, help="Min NLOD")
    parser.add_argument("--min-nalod", type=float, help="Min NALOD")
    parser.add_argument("--max-ecnt", type=int, help="Max ECNT")
    parser.add_argument("--min-mfrl", type=float, help="Min MFRL")
    parser.add_argument("--min-mpos", type=float, help="Min MPOS")
    parser.add_argument("--exclude-str", action="store_true", help="Exclude STRs")
    parser.add_argument("--min-af", type=float, help="Min Tumor AF")
    parser.add_argument("--min-alt-ad", type=int, help="Min Alt AD")
    parser.add_argument("--min-f1r2", type=int, help="Min F1R2")
    parser.add_argument("--min-f2r1", type=int, help="Min F2R1")
    parser.add_argument("--require-phased", action="store_true", help="Require phasing")

    if len(sys.argv) == 1:
        interactive_mode()
    else:
        args = parser.parse_args()
        if not args.input:
            print("Error: -i/--input parameter is required.", file=sys.stderr)
            sys.exit(1)
        run_pipeline(args)
