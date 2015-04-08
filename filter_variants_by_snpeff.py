#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import vcf
from snpeff import SnpEffAnn


def read_effect_list(filename):
    """
    Given a file name, read a list of effects from it and return
    their set.

    :param filename: a name of a file with a list of variant effects
    :type filename: str
    :return: a set of variant effects
    :rtype: set
    """
    result = []
    with open(filename) as input_file:
        for line in input_file:
            line = line.rstrip()
            result.append(line)
    return set(result)

def main():
    """
    This is the main function of the script.
    """
    parser = argparse.ArgumentParser(description='Filter variants by '
                                                 'their snpEff '
                                                 'annotation.')
    parser.add_argument('vcf_file', help='a VCF file of annotated '
                                         'variants')
    parser.add_argument('effect_file', help='a text file of variant '
                                            'effects')
    parser.add_argument('output_file', help='the output VCF file')
    args = parser.parse_args()

    current_effects = read_effect_list(args.effect_file)

    vcf_reader = vcf.Reader(filename=args.vcf_file)
    vcf_writer = vcf.Writer(open(args.output_file, 'w'), vcf_reader)
    for record in vcf_reader:
        annotation = SnpEffAnn(record)
        if annotation.effects[0].effect in current_effects:
            record.INFO['ANN'] = ['|'.join(annotation.effects[0])]
            vcf_writer.write_record(record)

if __name__ == '__main__':
    main()
