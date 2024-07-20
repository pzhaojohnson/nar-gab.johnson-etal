import json

import matplotlib.pyplot as plt

from random import shuffle


print()


#blast_output_file_path = 'blast_output_ivt_cy1_gRNA.json'
blast_output_file_path = 'blast_output_cy1_pemv2_pfbv_in_line.json'

#blast_output_file_path = 'blast_output_cy1_nb_2wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_2wpi_root.json'
#blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_6wpi_root.json'

print('BLAST output file path: ' + blast_output_file_path);
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output.')
print()


reports = [output['report'] for output in blast_output['BlastOutput2']]
print('Reports: ' + str(len(reports)))
print()


searches = [report['results']['search'] for report in reports]
print('Searches: ' + str(len(searches)))
print()


def has_a_hit(search):
    hits = search['hits']
    assert type(hits) == list
    return len(hits) > 0


searches_with_a_hit = list(filter(has_a_hit, searches))
print('Searches with a hit: ' + str(len(searches_with_a_hit)))
print()


def is_cy1_hit(hit):
    return hit['description'][0]['title'].lower() == 'cy1'


for search in searches_with_a_hit:
    hits = search['hits']
    assert type(hits) == list
    assert len(hits) == 1
print('All searches have at most one hit.')
print()


for search in searches_with_a_hit:
    hits = search['hits']
    assert type(hits) == list
    for hit in hits:
        assert is_cy1_hit(hit)
print('All hits are CY1.')
print()


def has_plus_hit_strand(hsp):
    return hsp['hit_strand'] == 'Plus'


def has_minus_hit_strand(hsp):
    return hsp['hit_strand'] == 'Minus'


def only_has_plus_hsps(hit):
    hsps = hit['hsps']
    assert type(hsps) == list
    plus_hsps = list(filter(has_plus_hit_strand, hsps))
    return len(plus_hsps) == len(hsps)


def only_has_minus_hsps(hit):
    hsps = hit['hsps']
    assert type(hsps) == list
    minus_hsps = list(filter(has_minus_hit_strand, hsps))
    return len(minus_hsps) == len(hsps)


def is_for_plus_strand_read(search):
    hits = search['hits']
    assert type(hits) == list
    assert len(hits) == 1
    hit = hits[0]
    return only_has_plus_hsps(hit)


def is_for_minus_strand_read(search):
    hits = search['hits']
    assert type(hits) == list
    assert len(hits) == 1
    hit = hits[0]
    return only_has_minus_hsps(hit)


def is_for_plus_minus_hybrid_read(search):
    hits = search['hits']
    assert type(hits) == list
    assert len(hits) == 1
    hit = hits[0]
    return len(hit['hsps']) >= 2 and not only_has_plus_hsps(hit) and not only_has_minus_hsps(hit)


searches_for_plus_strand_reads = list(filter(is_for_plus_strand_read, searches_with_a_hit))
print('Searches for plus-strand reads: ' + str(len(searches_for_plus_strand_reads)))
print()

searches_for_minus_strand_reads = list(filter(is_for_minus_strand_read, searches_with_a_hit))
print('Searches for minus-strand reads: ' + str(len(searches_for_minus_strand_reads)))
print()

searches_for_plus_minus_hybrid_reads = list(filter(is_for_plus_minus_hybrid_read, searches_with_a_hit))
print('Searches for plus-minus hybrid reads: ' + str(len(searches_for_plus_minus_hybrid_reads)))
print()


searches_accounted_for = (
    searches_for_plus_strand_reads
    + searches_for_minus_strand_reads
    + searches_for_plus_minus_hybrid_reads
)

# check that all searches with a hit are accounted for
assert len(searches_accounted_for) == len(searches_with_a_hit)
for search in searches_with_a_hit:
    assert search in searches_accounted_for
print('All searches with a hit are accounted for.')
print()


def covered_hit_positions(hsp):
    hit_from = hsp['hit_from']
    hit_to = hsp['hit_to']
    assert type(hit_from) == int
    assert type(hit_to) == int
    return [p for p in range(hit_from, hit_to + 1)]


def unique_covered_hit_positions(hit):
    hsps = hit['hsps']
    assert type(hsps) == list
    return set([p for hsp in hsps for p in covered_hit_positions(hsp)])


#"""
print('Creating coverage map... (This might take a little bit...)')
print()
fig, ax = plt.subplots()

plt.hist(
    [p for search in searches_for_plus_strand_reads for p in unique_covered_hit_positions(search['hits'][0])],
    bins=[i + 0.5 for i in range(0, 2692 + 1)],
    color='black',
)

print('Showing coverage map... (This might take a little bit...)')
print()
plt.show()
#"""
