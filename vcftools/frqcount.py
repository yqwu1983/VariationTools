#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from collections import namedtuple
import gzip
import logging

logging.basicConfig()
logger = logging.getLogger(__name__)


class VcfToolsError(Exception):
    """
    This class represent any exception raised by a class in the
    vcftools module
    """
    pass


class FrqCount(object):
    """
    This class provides routines to process allele frequency files
    (*.frq.count) produced by VCFtools using its --counts option.
    """

    __numeric_fields = (1, 2, 3)

    __frq_count_names = ('chrom', 'pos', 'n_alleles', 'n_chr',
                         'ref', 'alt')

    FrqCountLine = namedtuple('FrqCountLine', __frq_count_names)

    @staticmethod
    def __parse_allele_counts(line, lineno):
        """
        Given a line of allele counts from an allele frequency file
        produced by VCFtools, parse it to a dictionary which keys are
        alleles and values are their counts.

        :param line: a line of allele counts to be parsed
        :param lineno: the current line number; it is used when
            displaying error messages
        :type line: str
        :type lineno: int
        :return: the dictionary representing allele counts from the
            specified line; its keys are allels and values are their
            counts
        """
        result = {}

        for allele_record in line.split():
            try:
                allele, count = allele_record.rsplit(':', 1)
            except ValueError:
                logger.error('line {0}: the incorrect allele record {'
                             '1}'.format(lineno, allele_record))
                raise VcfToolsError
            # check if the count value is a valid integer
            try:
                count = int(count)
            except ValueError:
                logger.error('line {0}: the incorrect allele count '
                             'value {1}'.format(lineno, count))
                raise VcfToolsError
            # check if the count allele value is unique, otherwise
            # raise the exception
            if allele not in result:
                result[allele] = count
            else:
                logger.error('line {0}: multiple allele counts for '
                             'the allele {1}'.format(lineno, allele))
                raise VcfToolsError

        return result

    @staticmethod
    def __parse_frq_count_line(line, lineno):
        """
        Given a line from a vcftools allele frequency file, parse it
        and return a FrqCountLine object.

        :param line: a line from a vcftools frequency count file
        :param lineno: the current line number; it is used when
            displaying error messages
        :type line: str
        :type lineno: int
        :return: a named tuple of the FrqCountLine class representing
            the specified line
        """
        line_parts = line.split(None, 5)

        # check if the line contains all required fields
        if len(line_parts) < 6:
            logger.error('line {0}: the incorrect number of '
                         'fields'.format(lineno))
            raise VcfToolsError

        # convert numeric values
        for i in FrqCount.__numeric_fields:
            try:
                line_parts[i] = int(line_parts[i])
            except ValueError:
                logger.error('line {0}: the incorrect numeric value {'
                             '1}'.format(lineno, line_parts[i]))
                raise VcfToolsError

        # parse the parts containing allele counts
        for i in (4, 5):
            line_parts[i] = FrqCount.__parse_allele_counts(
                line_parts[i], lineno)

        return FrqCount.FrqCountLine(*line_parts)

    def __init__(self, filename, gzipped=False):
        """
        Create an FrqCount object from the specified file.

        :param filename: a name of a vcftools frequency count file
        :param gzipped: is the specified file gzipped or not
        :type filename: str
        :type gzipped: bool
        :return: an FrqCount object containing frequency counts from
            the specified file
        :rtype: FrqCount
        """
        self.__records = dict()

        # We store records obtained from the specified file in a
        # dictionary which keys are chromosome names to sort them by
        # two keys: first by their position on a chromosome and
        # second by the chromosome name. Finally the tuple of sorted
        # allele frequency counts in produced.
        temp_records = dict()
        file_handler = open if not gzipped else gzip.open
        with file_handler(filename) as count_file:
            lineno = 0
            for line in count_file:
                line = line.rstrip()
                lineno += 1
                if line.startswith('CHROM'):
                    # skip the comment line
                    continue
                new_record = self.__parse_frq_count_line(line, lineno)
                if new_record.chrom in temp_records:
                    temp_records[new_record.chrom][
                        new_record.pos] = new_record
                else:
                    temp_records[new_record.chrom] = dict()
                    temp_records[new_record.chrom][new_record.pos] \
                        = new_record

        self.__records = temp_records

    def records(self):
        """
        Return the iterator to iterate through the allele frequency
        records.

        :return: the allele count record iterator
        """
        return self.__records
