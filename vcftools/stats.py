#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from collections import namedtuple
from vcftools.frqcount import FrqCount


class AlleleComparator(object):
    """
    This class implements routines to compare allele counts between
    two populations.
    """

    __comparison_features = ('chrom', 'pos', 'ref', 'alt',
                             'first_ref_counts',
                             'first_alt_counts',
                             'second_ref_counts',
                             'second_alt_counts')

    ComparisonLine = namedtuple('AlleleComparisonLine',
                                __comparison_features)

    def __init__(self, first, second):
        """
        Create an AlleleComparator object given names of two allele
        count files.

        :param first: a name of the first file of allele counts
        :param second: a name of the second file of allele counts
        :type first: str
        :type second: str
        :return: an AlleleComparator object comparing allele counts
            from the specified pair of files
        :rtype: AlleleComparator
        """
        self.__first = FrqCount(first, first.endswith('.gz'))
        self.__second = FrqCount(second, second.endswith('.gz'))

    @staticmethod
    def __combine_allele_counts(x, y):
        """
        Given two allele frequency records, produce a ComparisonLine
        object for them.

        :param x: the first allele count object
        :param y: the second allele count object
        :type x: FrqCount.FrqCountLine
        :type y: FrqCount.FrqCountLine
        :return: the ComparisonLine object representing the allele
            count comparison result
        :rtype: AlleleComparator.ComparisonLine
        """
        x_alt = set(x.alt.keys())
        y_alt = set(y.alt.keys())

        common_alt_alleles = x_alt.intersection(y_alt)
        x_unique_alt = x_alt.difference(y_alt)
        y_unique_alt = y_alt.difference(x_alt)

        result = []
        for i in common_alt_alleles:
            result.append(AlleleComparator.ComparisonLine(
                chrom=x.chrom,
                pos=x.pos,
                ref=x.ref.keys()[0],
                alt=i,
                first_ref_counts=x.ref.values()[0],
                first_alt_counts=x.alt[i],
                second_ref_counts=y.ref.values()[0],
                second_alt_counts=y.alt[i]
            ))
        for i in x_unique_alt:
            result.append(AlleleComparator.ComparisonLine(
                chrom=x.chrom,
                pos=x.pos,
                ref=x.ref.keys()[0],
                alt=i,
                first_ref_counts=x.ref.values()[0],
                first_alt_counts=x.alt[i],
                second_ref_counts=y.ref.values()[0],
                second_alt_counts=0
            ))
        for i in y_unique_alt:
            result.append(AlleleComparator.ComparisonLine(
                chrom=x.chrom,
                pos=x.pos,
                ref=x.ref.keys()[0],
                alt=i,
                first_ref_counts=x.ref.values()[0],
                first_alt_counts=0,
                second_ref_counts=y.ref.values()[0],
                second_alt_counts=y.alt[i]
            ))

        return result

    def comparisons(self):
        """
        Produces the iterator to comparisons of allele counts.
        """
        first = self.__first.records()
        second = self.__second.records()

        common_chromosomes = set(first.keys()).intersection(
            set(second.keys()))

        for i in common_chromosomes:
            for pos in sorted(first[i].keys()):
                if pos in second[i]:
                    yield self.__combine_allele_counts(first[i][pos],
                                                       second[i][pos])
