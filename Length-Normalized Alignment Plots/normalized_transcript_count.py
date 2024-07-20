import json

import functools

import matplotlib.pyplot as plt


print()


blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'
print(f'{blast_output_file_path=}')
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output!')
print()


reports = [item['report'] for item in blast_output['BlastOutput2']]
print(f'{len(reports)=}')
print()


all_searches = [report['results']['search'] for report in reports]
print(f'{len(all_searches)=}')
print()


def is_to_CY1(hit):
    return hit['description'][0]['title'].upper() == 'CY1'


def is_for_CY1_read(search):
    hits = search['hits']
    assert len(hits) <= 1
    return len(hits) == 1 and is_to_CY1(hits[0])


searches_for_CY1_reads = list(filter(is_for_CY1_read, all_searches))
print(f'{len(searches_for_CY1_reads)=}')
print()


def has_plus_hit_strand(hsp):
    return hsp['hit_strand'] == 'Plus'


def only_has_plus_hsps(hit):
    hsps = hit['hsps']
    assert len(hsps) > 0
    return all([has_plus_hit_strand(hsp) for hsp in hsps])


def is_for_plus_strand_CY1_read(search):
    hits = search['hits']
    assert len(hits) <= 1
    return is_for_CY1_read(search) and only_has_plus_hsps(hits[0])


searches_for_plus_strand_CY1_reads = list(filter(is_for_plus_strand_CY1_read, searches_for_CY1_reads))
print(f'{len(searches_for_plus_strand_CY1_reads)=}')
print()


def is_for_full_length_plus_strand_CY1_gRNA_read(search):
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert is_to_CY1(hit)
    hsps = hit['hsps']
    assert len(hsps) > 0
    first_hsp = hsps[0]
    return len(hsps) == 1 \
        and has_plus_hit_strand(first_hsp) \
        and first_hsp['hit_from'] <= 31 \
        and 2692 - 10 <= first_hsp['hit_to']


searches_for_full_length_plus_strand_CY1_gRNA_reads = list(filter(
    is_for_full_length_plus_strand_CY1_gRNA_read,
    searches_for_plus_strand_CY1_reads,
))
print(f'{len(searches_for_full_length_plus_strand_CY1_gRNA_reads)=}')
print()


def is_for_plus_strand_CY1_F281_read(search):
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert is_to_CY1(hit)
    hsps = hit['hsps']
    assert len(hsps) > 0
    hsp1 = hsps[0]
    return len(hsps) == 1 \
        and has_plus_hit_strand(hsp1) \
        and hsp1['hit_from'] <= 31 \
        and 281 - 10 <= hsp1['hit_to'] <= 281 + 10


searches_for_plus_strand_CY1_F281_reads = list(filter(
    is_for_plus_strand_CY1_F281_read,
    searches_for_plus_strand_CY1_reads,
))
print(f'{len(searches_for_plus_strand_CY1_F281_reads)=}')
print()


def is_for_plus_strand_CY1_F671_read(search):
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert is_to_CY1(hit)
    hsps = hit['hsps']
    assert len(hsps) > 0
    hsp1 = hsps[0]
    return len(hsps) == 1 \
        and has_plus_hit_strand(hsp1) \
        and hsp1['hit_from'] <= 31 \
        and 671 - 10 <= hsp1['hit_to'] <= 671 + 10


searches_for_plus_strand_CY1_F671_reads = list(filter(
    is_for_plus_strand_CY1_F671_read,
    searches_for_plus_strand_CY1_reads,
))
print(f'{len(searches_for_plus_strand_CY1_F671_reads)=}')
print()


def cmp_query_froms(hsp1, hsp2):
    return hsp1['query_from'] - hsp2['query_from']


def is_for_two_segment_read(search):
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert is_to_CY1(hit)
    hsps = hit['hsps']
    return len(hsps) == 2


searches_for_two_segment_reads = list(filter(is_for_two_segment_read, searches_for_CY1_reads))
print(f'{len(searches_for_two_segment_reads)=}')
print()


def is_for_plus_strand_CY1_DRNA_read(search):
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert is_to_CY1(hit)
    hsps = hit['hsps']
    assert len(hsps) == 2
    hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
    hsp1 = hsps[0]
    hsp2 = hsps[1]
    return has_plus_hit_strand(hsp1) \
        and has_plus_hit_strand(hsp2) \
        and hsp1['hit_from'] <= 31 \
        and 2692 - 10 <= hsp2['hit_to'] \
        and 944 - 10 <= hsp1['hit_to'] + (hsp2['hit_to'] - hsp2['hit_from']) <= 944 + 10


searches_for_plus_strand_CY1_DRNA_reads = list(filter(is_for_plus_strand_CY1_DRNA_read, searches_for_two_segment_reads))
print(f'{len(searches_for_plus_strand_CY1_DRNA_reads)=}')
print()


def is_for_F442_read(search):
    """Returns True if the search is for a F442 read."""
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert is_to_CY1(hit)
    hsps = hit['hsps']
    assert len(hsps) > 0
    hsp1 = hsps[0]
    return len(hsps) == 1 \
        and has_plus_hit_strand(hsp1) \
        and hsp1['hit_from'] <= 31 \
        and 442 - 10 <= hsp1['hit_to'] <= 442 + 10


searches_for_F442_reads = list(filter(is_for_F442_read, searches_for_plus_strand_CY1_reads))
print(f'{len(searches_for_F442_reads)=}')
print()


def is_for_F944_read(search):
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert is_to_CY1(hit)
    hsps = hit['hsps']
    assert len(hsps) > 0
    hsp1 = hsps[0]
    return len(hsps) == 1 \
        and has_plus_hit_strand(hsp1) \
        and hsp1['hit_from'] <= 31 \
        and 944 - 10 <= hsp1['hit_to'] <= 944 + 10


searches_for_F944_reads = list(filter(is_for_F944_read, searches_for_plus_strand_CY1_reads))
print(f'{len(searches_for_F944_reads)=}')
print()
