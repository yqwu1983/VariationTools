#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

import vcf
from collections import namedtuple


class VariantZygosity(object):
    """
    This class implements routines to collect zygosity statistics for
    each variant from a VCF file.
    """
    __stat_record_names = ('id', 'chr', 'pos', 'hom_ref', 'het',
                           'hom_alt')

    ZygosityStatRecord = namedtuple('ZygosityStatRecord',
                                    __stat_record_names)

    def __init__(self, filename):
        """
        Create a VariantZygosity object from the specified VCF file.

        :param filename: a name of a VCF file
        :type filename: str
        :return: a VariantZygosity object created from the specified
            VCF file
        :rtype VariantZygosity
        """
        self.__filename = filename

    def variants(self):
        """
        Return the iterator that iterates through variants and
        returns zygosity statistics for each of them.

        :return: the itetator that iterates through variants of the
            VCF file the object was created from
        """
        for variant in vcf.Reader(open(self.__filename)):
            record = VariantZygosity.ZygosityStatRecord(
                id=variant.ID if variant.ID else '.',
                chr=variant.CHROM,
                pos=variant.POS,
                hom_ref=variant.num_hom_ref,
                het=variant.num_het,
                hom_alt=variant.num_hom_alt
            )
            yield record


class SampleZygosity(object):
    """
    This class implements routines to collect zygosity statistics for
    each sample from a VCF file.
    """
    __stat_record_names = ('individual', 'hom_ref', 'het', 'hom_alt')

    ZygosityStatRecord = namedtuple('ZygosityStatRecord',
                                    __stat_record_names)

    def __init__(self, filename):
        """
        Create a SampleZygosity object from the specified VCF file.

        :param filename: a name of a VCF file
        :type filename: str
        :return: a SampleZygosity object created from the specified
            VCF file
        :rtype: SampleZygosity
        """
        self.__filename = filename

    def samples(self):
        """
        Return the iterator that iterates through samples and returns
        zygosity statistics for each of them.

        :return: the iterator that iterates through samples of the
            VCF file the object was created from
        """
        all_samples = vcf.Reader(open(self.__filename)).samples
        stats = {sample: [0, 0, 0] for sample in all_samples}

        for variant in vcf.Reader(open(self.__filename)):
            for i in variant.samples:
                stats[i.sample][i.gt_type] += 1

        for sample in sorted(all_samples):
            yield SampleZygosity.ZygosityStatRecord(sample,
                                                    *stats[sample])
