#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

from collections import namedtuple


class Reader(object):
    """
    This class implements a reader for output files produced by
    bcftools roh.
    """
    Variant = namedtuple('RohVariant', ('seq', 'pos', 'state',
                                        'quality'))

    def __init__(self, filename):
        """
        Creates the object to read ROH data from the specified file.

        :param filename: a name of a file produced by bcftools roh
        :type filename: str
        """
        self.__filename = filename

    def variants(self):
        """
        Iterate through variants that were processed by bcftools roh.

        :return: a variant with its ROH annotation
        :rtype: RohVariant
        """
        with open(self.__filename) as input_file:
            for line in input_file:
                if line.startswith('#'):
                    # skip comments
                    continue
                line = line.rstrip()
                line_parts = line.split('\t')
                line_parts[1] = int(line_parts[1])
                line_parts[2] = line_parts[2] == '1'
                yield Reader.Variant(*line_parts)


class Run(object):
    """
    This class implements routines to obtain runs of homozygosity
    from output files produced by bcftools roh
    """
    Region = namedtuple('Region', ('seq', 'start', 'end'))

    def __init__(self, filename):
        """
        Given a name of an output file produced by bcftools roh,
        create a Run object from it.

        :param filename: a name of a file produced by bcftools roh
        :type filename: str
        """
        self.__filename = filename

    def regions(self):
        """
        Iterate through runs of homozygosity.

        :return: a run of homozygosity, that is, a region in a genome
        :rtype: Region
        """
        # we utilize a simple finite-state machine
        curr_seq = None
        start_variant = None
        end_variant = None
        reader = Reader(self.__filename)
        for variant in reader.variants():
            if variant.seq != curr_seq:
                # we moved to a new sequence, check if there is an ROH
                # to be output
                if end_variant is not None:
                    yield Run.Region(
                        seq=curr_seq,
                        start=start_variant.pos,
                        end=end_variant.pos)
                    start_variant = None
                    end_variant = None
                # change the current sequence
                curr_seq = variant.seq
                if variant.state:
                    start_variant = variant
                else:
                    start_variant = None
            else:
                if variant.state:
                    # we read an ROH variant
                    if start_variant is None:
                        start_variant = variant
                    else:
                        end_variant = variant
                else:
                    # we read a non-ROH variant and check if a
                    # ROH is detected
                    if start_variant is not None:
                        if end_variant is None:
                            # there was a single ROH variant
                            # surrounded by non-ROH ones, it does
                            # not form an ROH region
                            start_variant = None
                        else:
                            # there were at least two ROH variants
                            # without any non-ROH variants within
                            # them
                            yield Run.Region(
                                seq=curr_seq,
                                start=start_variant.pos,
                                end=end_variant.pos)
                            start_variant = None
                            end_variant = None
