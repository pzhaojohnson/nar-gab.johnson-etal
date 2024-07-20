import json

import matplotlib.pyplot as plt

import random

import functools


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


def is_for_plus_strand_read(search):
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert is_to_CY1(hit)
    hsps = hit['hsps']
    assert len(hsps) > 0
    return all([hsp['hit_strand'] == 'Plus' for hsp in hsps])


searches_for_plus_strand_reads = list(filter(is_for_plus_strand_read, searches_for_CY1_reads))
print(f'{len(searches_for_plus_strand_reads)=}')
print()


length_multiplied = []
for search in searches_for_plus_strand_reads:
    read_length = search['query_len']
    for i in range(read_length):
        length_multiplied.append(search)
print(f'{len(length_multiplied)=}')
print()


def cmp_read_lengths(search1, search2):
    return search1['query_len'] - search2['query_len']


def num_aligned_positions(search):
    assert is_for_plus_strand_read(search)
    hit = search['hits'][0]
    hsps = hit['hsps']
    return sum([hsp['hit_to'] - hsp['hit_from'] + 1 for hsp in hsps])


def cmp_num_aligned_positions(search1, search2):
    return num_aligned_positions(search1) - num_aligned_positions(search2)


def is_for_full_length_gRNA_read(search):
    assert is_for_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    return len(hsps) == 1 \
        and hsps[0]['hit_from'] <= 31 \
        and 2692 - 10 <= hsps[0]['hit_to']


def is_for_F281_read(search):
    assert is_for_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    return len(hsps) == 1 \
        and hsps[0]['hit_from'] <= 31 \
        and 281 - 10 <= hsps[0]['hit_to'] <= 281 + 10


def is_for_F442_read(search):
    assert is_for_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    return len(hsps) == 1 \
        and hsps[0]['hit_from'] <= 31 \
        and 442 - 10 <= hsps[0]['hit_to'] <= 442 + 10


def is_for_F671_read(search):
    assert is_for_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    return len(hsps) == 1 \
        and hsps[0]['hit_from'] <= 31 \
        and 671 - 10 <= hsps[0]['hit_to'] <= 671 + 10


def is_for_F1070_read(search):
    assert is_for_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    return len(hsps) == 1 \
        and hsps[0]['hit_from'] <= 31 \
        and 1070 - 10 <= hsps[0]['hit_to'] <= 1070 + 10


def cmp_query_froms(hsp1, hsp2):
    return hsp1['query_from'] - hsp2['query_from']


def is_for_DRNA_read(search):
    assert is_for_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
    if len(hsps) == 2:
        assert hsps[0]['query_from'] < hsps[1]['query_from']
    return len(hsps) == 2 \
        and hsps[0]['hit_from'] <= 31 \
        and hsps[1]['hit_to'] >= 2692 - 10 \
        and 944 - 10 <= hsps[0]['hit_to'] + (hsps[1]['hit_to'] - hsps[1]['hit_from']) <= 944 + 10


print(f'{len(searches_for_plus_strand_reads)=}')
print()

print(f'{len(list(filter(is_for_full_length_gRNA_read, searches_for_plus_strand_reads)))=}')
print()

print(f'{len(list(filter(is_for_F281_read, searches_for_plus_strand_reads)))=}')
print()

print(f'{len(list(filter(is_for_F442_read, searches_for_plus_strand_reads)))=}')
print()

print(f'{len(list(filter(is_for_F671_read, searches_for_plus_strand_reads)))=}')
print()

print(f'{len(list(filter(is_for_F1070_read, searches_for_plus_strand_reads)))=}')
print()

print(f'{len(list(filter(is_for_DRNA_read, searches_for_plus_strand_reads)))=}')
print()



"""
fig, ax = plt.subplots()

random.shuffle(length_multiplied)

subset = length_multiplied[:int(1.5e4)]
subset.sort(key=functools.cmp_to_key(cmp_num_aligned_positions), reverse=True)

y = 0

for search in subset:
    if is_for_full_length_gRNA_read(search):
        color = 'red'
    elif is_for_F281_read(search):
        color = 'dodgerblue'
    elif is_for_F442_read(search):
        color = 'magenta'
    elif is_for_F671_read(search):
        color = 'orange'
    elif is_for_F1070_read(search):
        color = 'blue'
    elif is_for_DRNA_read(search):
        color = 'lime'
    else:
        color = 'black'
    alpha = 0.05
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert is_to_CY1(hit)
    hsps = hit['hsps']
    for hsp in hsps:
        plt.plot([hsp['hit_from'], hsp['hit_to']], [y, y], color=color, alpha=alpha, linewidth=1)
    y += 1

plt.show()
"""
