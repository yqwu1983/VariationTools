#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import bcftools.roh


def parse_command_line_args():
    """
    Parse command-line agruments of the script.

    :return: a tuple of command-line arguments
    """
    parser = argparse.ArgumentParser(
        description='Produce a BED track of ROHs from bcftools roh '
                    'output')
    parser.add_argument('roh_file', help='a file produced by '
                                         'bcftools roh')
    parser.add_argument('output_file', help='a name of the output '
                                            'BED file')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_command_line_args()
    template = '\t'.join(['{}']*3) + '\n'
    with open(args.output_file, 'w') as output:
        for region in bcftools.roh.Run(args.roh_file).regions():
            # we output a BED file, so start coordinates must start with
            # 0 and end coordinates must be shifted by 1; variant
            # positions are 1-based and their ends are not shifted,
            # so we substract 1 from start coordinates
            output.write(template.format(region.seq,
                                         region.start - 1,
                                         region.end))
