#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import argparse
import logging
import vcf
from collections import namedtuple

logging.basicConfig()
logger = logging.getLogger(__name__)

GenotypeCounts = namedtuple('GenotypeCounts',
                            ('ref_hom', 'ref_alt_het', 'alt_het_hom'))


class MultiVCFReader(object):
    """
    This class implements obtaining variants from a set of VCF files
    each of that stores variants located on a single chromosome. For
    example, VCF files of the 1000 Genomes project are organized in
    this way.
    """
    def __init__(self, vcf_files):
        """
        Given a list of VCF files, create a MultiVCFReader object
        corresponding to them.

        :param vcf_files: a list of VCF file paths
        :type vcf_files: list
        :return: a MultiVCFReader object providing access to variants
            in the specified VCF files
        :rtype: MultiVCFReader
        """
        self.__readers = {}
        for i in vcf_files:
            first_variant = vcf.Reader(open(i)).next()
            reader_chrom = first_variant.CHROM
            if reader_chrom not in self.__readers:
                self.__readers[reader_chrom] = vcf.Reader(open(i))
            else:
                logger.error('multiple VCF files of the chromosome {'
                             '}'.format(reader_chrom))
                raise Exception

    def get_variant(self, chrom, pos):
        """
        Get a record for a variant situated at the specified position.

        :param chrom: a variant chromosome
        :param pos: a variant position
        :type chrom: str
        :type pos: int
        :return: a variant record representing a variant at the
            specified position or None if there is no variant at the
            position
        :rtype: vcf.model._Record
        """
        result = None
        if chrom in self.__readers:
            variant = list(self.__readers[chrom].fetch(chrom, pos, pos))
            result = variant[0] if variant else None

        return result


def count_genotypes(record, samples=None):
    """
    Given a record from a VCF file, return genotype counts for it.

    :param record: a record from a VCF file
    :param samples: a list of sample names
    :type record: vcf.model._Record
    :type samples: list
    :return: a tuple of genotype counts
    :rtype: GenotypeCounts
    """
    ref_hom = ref_alt_het = alt_het_hom = 0
    for i in record.samples:
        if samples is None or i.sample in samples:
            genotype = sorted(map(int, i.gt_alleles))
            if sum(genotype) == 0:
                ref_hom += 1
            elif genotype[0] == 0:
                ref_alt_het += 1
            else:
                alt_het_hom += 1

    return GenotypeCounts(ref_hom=ref_hom, ref_alt_het=ref_alt_het,
                          alt_het_hom=alt_het_hom)


def main():
    """
    This is the main function that is called when the script is
    launched.
    """
    parser = argparse.ArgumentParser(description='Analyze genotype '
                                                 'counts.')
    parser.add_argument('vcf_file',
                        help='a VCF file of genomic variants')
    parser.add_argument('--compare', nargs='+',
                        help='compare genotype counts with variants '
                             'from 1000 Genomes')
    parser.add_argument('--samples',
                        help='a file with sample names for comparison')
    parser.add_argument('output', help='an output file')

    args = parser.parse_args()

    samples = None
    if args.samples:
        samples = []
        for line in open(args.samples):
            samples.append(line.rstrip())

    comparator = MultiVCFReader(args.compare) if args.compare else None
    template = '\t'.join(['{}'] * 8) + '\n'
    template_compared = '\t'.join(['{}'] * 12) + '\n'

    with open(args.output, 'w') as output_file:
        for variant in vcf.Reader(open(args.vcf_file)):
            genotype_counts = count_genotypes(variant)
            variant_type = 'SNV' if variant.is_snp else 'INDEL'
            variant_len = max(len(variant.REF),
                              reduce(max, map(len, variant.ALT)))
            if comparator:
                compared_variant = comparator.get_variant(
                    variant.CHROM, variant.POS)
                if compared_variant:
                    compared_genotype_counts = count_genotypes(
                        compared_variant, samples)
                    compared_aaf = \
                        (compared_genotype_counts.alt_het_hom * 2 +
                         compared_genotype_counts.ref_alt_het) / \
                        float(sum(compared_genotype_counts)*2)
                else:
                    # the variant is not present in the specified VCF
                    # files
                    compared_genotype_counts = GenotypeCounts(
                        *(['NA'] * 3))
                    compared_aaf = 'NA'
                output_file.write(template_compared.format(
                    variant.CHROM, variant.POS,
                    variant.aaf[0],
                    compared_aaf,
                    variant_type, variant_len,
                    *genotype_counts + compared_genotype_counts))
            else:
                output_file.write(template.format(
                    variant.CHROM, variant.POS,
                    variant.aaf,
                    variant_type, variant_len,
                    *genotype_counts))


if __name__ == '__main__':
    main()
