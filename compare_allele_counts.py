#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import logging
from vcftools.stats import AlleleComparator

logging.basicConfig()
logger = logging.getLogger(__name__)


def main():
    """
    This is the main function that is called when the script is
    launched.
    """
    parser = argparse.ArgumentParser(
        description='Compare allele counts in the specified vcftools '
                    '.frq.count files.')

    parser.add_argument('first', help='the first file of allele '
                                      'counts')
    parser.add_argument('second', help='the second file of allele '
                                       'counts')
    parser.add_argument('output', help='the output file name')

    parser.add_argument('--debug', action='store_true',
                        help='enable the debug mode (diasnostic '
                             'messages are printed)')

    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.INFO)

    logger.info('{0} - {1}: loading started'.format(args.first,
                                                    args.second))
    comparator = AlleleComparator(args.first, args.second)
    logger.info('{0} - {1}: loading completed'.format(args.first,
                                                      args.second))

    template = '\t'.join(['{}'] * 8) + '\n'

    logger.info('{0} - {1}: writing started'.format(args.first,
                                                    args.second))

    with open(args.output, 'w') as output_file:
        for i in comparator.comparisons():
            for j in i:
                output_file.write(template.format(*j))

    logger.info('{0} - {1}: writing completed'.format(args.first,
                                                      args.second))


if __name__ == '__main__':
    main()
