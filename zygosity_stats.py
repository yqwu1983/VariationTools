#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import cStringIO
from vcfstats.zygosity import SampleZygosity
from vcfstats.zygosity import VariantZygosity


def main():
    """
    The main function that is called when the script is launched.
    """
    parser = argparse.ArgumentParser(
        description='Get zygosity statistics for each variant.')

    parser.add_argument('vcf_file', help='a VCF file of genomic '
                                         'variants')
    parser.add_argument('output_file', help='the output file name')
    parser.add_argument('--samples', action='store_true',
                        help='get heterozygosity statistics for each '
                             'sample')

    args = parser.parse_args()

    output = cStringIO.StringIO()

    if args.samples:
        variant_statistics = SampleZygosity(args.vcf_file)
        template = '\t'.join(['{}'] * 4) + '\n'
        output.write('SAMPLE\t#HOM_REF\t#HET\t#HOM_ALT\n')
        for sample in variant_statistics.samples():
            output.write(template.format(*sample))
    else:
        variant_statistics = VariantZygosity(args.vcf_file)
        template = '\t'.join(['{}'] * 6) + '\n'
        output.write('SAMPLE\tCHROM\tPOS\t#HOM_REF\t#HET\t#HOM_ALT\n')
        for variant in variant_statistics.variants():
                output.write(template.format(*variant))

    with open(args.output_file, 'w') as output_file:
        output_file.write(output.getvalue())
    output.close()


if __name__ == '__main__':
    main()
