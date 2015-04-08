#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2015 by Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

"""
This file contains routines to process VCF files which variant
effects were annotated by snpEff.
"""

import logging
import vcf
from collections import namedtuple

logging.basicConfig()
logger = logging.getLogger(__name__)


class SnpEffAnnError(Exception):
    """
    This class represents exceptions that may occur when running
    routines related to snpEff annotations.
    """
    pass


class SnpEffAnn(object):
    """
    This class implements routines to process snpEff annotations of
    variants in a VCF file.
    """

    __snpeff_all_effects = ('chromosome_number_variation',
                            'exon_loss_variant',
                            'frameshift_variant',
                            'stop_gained',
                            'stop_lost',
                            'start_lost',
                            'splice_acceptor_variant',
                            'splice_donor_variant',
                            'rare_amino_acid_variant',
                            'missense_variant',
                            'inframe_insertion',
                            'disruptive_inframe_insertion',
                            'inframe_deletion',
                            'disruptive_inframe_deletion',
                            '5_prime_UTR_truncation+exon_loss_variant',
                            '3_prime_UTR_truncation+exon_loss',
                            '5_prime_UTR_truncation',
                            '3_prime_UTR_truncation',
                            'splice_branch_variant',
                            'splice_region_variant',
                            'splice_branch_variant',
                            'stop_retained_variant',
                            'initiator_codon_variant',
                            'synonymous_variant',
                            'initiator_codon_variant+non_canonical_start_codon',
                            'stop_retained_variant',
                            'coding_sequence_variant',
                            '5_prime_UTR_variant',
                            '3_prime_UTR_variant',
                            '5_prime_UTR_premature_start_codon_gain_variant',
                            'upstream_gene_variant',
                            'downstream_gene_variant',
                            'TF_binding_site_variant',
                            'regulatory_region_variant',
                            'miRNA',
                            'transcript',
                            'custom',
                            'sequence_feature',
                            'conserved_intron_variant',
                            'intron_variant',
                            'intragenic_variant',
                            'conserved_intergenic_variant',
                            'intergenic_region',
                            'coding_sequence_variant',
                            'non_coding_exon_variant',
                            'nc_transcript_variant',
                            'gene_variant',
                            'chromosome')

    __snpeff_ann_keys = dict(zip(__snpeff_all_effects,
                                 range(len(__snpeff_all_effects))))

    __snpeff_ann_fields = ('allele',
                           'effect',
                           'impact',
                           'gene_name',
                           'gene_id',
                           'feature_type',
                           'feature_id',
                           'transcript_biotype',
                           'exon_position',
                           'variant',
                           'protein_variant',
                           'cdna_position',
                           'cds_position',
                           'protein_position',
                           'feature_distance',
                           'info')

    SnpEffRecord = namedtuple('SnpEffRecord', __snpeff_ann_fields)

    def __parse_snpeff_ann(self):
        """
        Given an snpEff annotation line from the INFO field of a VCF
        file, parse it to the SnpEffRecord named tuple object.

        :return: a named tuple of the SnpEffRecord type representing
            the specified annotation
        :rtype: SnpEffAnnotation.SnpEffRecord
        """
        if 'ANN' not in self.__vcf_record.INFO:
            logger.error('no snpEff annotation in the record')
            raise SnpEffAnnError
        else:
            annotations = map(lambda x: x.split('|'),
                              self.__vcf_record.INFO['ANN'])
            # process composite effects, if any
            for i in xrange(len(annotations)):
                annotations[i][1] = annotations[i][1].split('&')[0]
            result = map(lambda x: SnpEffAnn.SnpEffRecord(*x),
                         annotations)
            # check the effects
            for i in xrange(len(result)):
                if result[i].effect not in \
                        SnpEffAnn.__snpeff_all_effects:
                    logger.error('missing effect {}'.format(
                        result[i].effect))
                    raise SnpEffAnnError
            return result

    def __sort_effects(self):
        """
        Sort annotated effects of the variants based on their impact.
        """
        self.effects.sort(key=lambda x: self.__snpeff_ann_keys[
            x.effect])

    def __init__(self, vcf_record):
        """
        Given a record from a VCF file, create an SnpEffAnnotation
        object from it.

        :param vcf_record: a VCF record
        :type vcf_record: vcf.model._Record
        """
        self.__vcf_record = vcf_record
        self.effects = self.__parse_snpeff_ann()
        self.__sort_effects()
