import json

import matplotlib.pyplot as plt

import random

import functools


print()


blast_output_file_path = 'all_reads_cy1_nb_6wpi_leaf_28S_rRNA_blast.json'
#blast_output_file_path = 'all_reads_cy1_nb_6wpi_leaf_18S_rRNA_blast.json'

#blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'

#blast_output_file_path = 'blast_output_ivt_cy1_gRNA.json'

print('BLAST output file path: ' + blast_output_file_path)
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output.')
print()


reports = [item['report'] for item in blast_output['BlastOutput2']]
print('Reports: ' + str(len(reports)))
print()

searches = [report['results']['search'] for report in reports]
print('Searches: ' + str(len(searches)))
print()


class ReadLength:
    @staticmethod
    def for_(search):
        """Returns the length of the read for the search."""
        read_length = search['query_len']
        assert type(read_length) == int
        return read_length


def cmp_read_lengths(search1, search2):
    return ReadLength.for_(search1) - ReadLength.for_(search2)


class Hits:
    @staticmethod
    def for_(search):
        """Returns the hits for the search."""
        hits = search['hits']
        assert type(hits) == list
        return hits


class Hit:
    @staticmethod
    def for_(search):
        """Returns the single hit for the search.

        Raises if the search does not have exactly one hit.
        """
        hits = Hits.for_(search)
        assert len(hits) == 1
        return hits[0]


class NumHits:
    @staticmethod
    def for_(search):
        """Returns the number of hits for the search."""
        hits = Hits.for_(search)
        return len(hits)


def has_a_hit(search):
    return NumHits.for_(search) > 0


searches_with_a_hit = list(filter(has_a_hit, searches))
print('Searches with a hit: ' + str(len(searches_with_a_hit)))
print()


for search in searches_with_a_hit:
    assert NumHits.for_(search) == 1
print('All searches have at most one hit.')
print()


def is_for_plus_strand_read(search):
    hit = Hit.for_(search)
    hsps = Hsps.of(hit)
    plus_hsps = PlusHsps.of(hit)
    return len(plus_hsps) == len(hsps)


def is_for_minus_strand_read(search):
    hit = Hit.for_(search)
    hsps = Hsps.of(hit)
    minus_hsps = MinusHsps.of(hit)
    return len(minus_hsps) == len(hsps)


def is_for_plus_minus_hybrid_read(search):
    return not is_for_plus_strand_read(search) and not is_for_minus_strand_read(search)


class Hsps:
    @staticmethod
    def of(hit):
        """Returns the hsps of the search."""
        hsps = hit['hsps']
        assert type(hsps) == list
        return hsps


class PlusHsps:
    @staticmethod
    def of(hit):
        """Returns the hsps of the hit with a plus hit strand."""
        hsps = Hsps.of(hit)
        return list(filter(has_plus_hit_strand, hsps))


class MinusHsps:
    @staticmethod
    def of(hit):
        """Returns the hsps of the hit with a minus hit strand."""
        hsps = Hsps.of(hit)
        return list(filter(has_minus_hit_strand, hsps))


class UniqueHitPositions:
    @staticmethod
    def covered_by(hit):
        """Returns a set of all unique hit positions covered by the hsps of the hit."""
        hsps = Hsps.of(hit)
        return set([p for hsp in hsps for p in HitPositions.covered_by(hsp)])


class MinHitPosition:
    @staticmethod
    def covered_by(hit):
        """Returns the minimum hit position covered by the hsps of the hit."""
        return min(UniqueHitPositions.covered_by(hit))


class MaxHitPosition:
    @staticmethod
    def covered_by(hit):
        """Returns the maximum hit position covered by the hsps of the hit."""
        return max(UniqueHitPositions.covered_by(hit))


class QueryFrom:
    @staticmethod
    def of(hsp):
        """Returns the query-from position of the hsp."""
        query_from = hsp['query_from']
        assert type(query_from) == int
        return query_from


class QueryTo:
    @staticmethod
    def of(hsp):
        """Returns the query-to position of the hsp."""
        query_to = hsp['query_to']
        assert type(query_to) == int
        return query_to


class HitFrom:
    @staticmethod
    def of(hsp):
        """Returns the hit-from position of the hsp."""
        hit_from = hsp['hit_from']
        assert type(hit_from) == int
        return hit_from


class HitTo:
    @staticmethod
    def of(hsp):
        """Returns the hit-to position of the hsp."""
        hit_to = hsp['hit_to']
        assert type(hit_to) == int
        return hit_to


def has_plus_hit_strand(hsp):
    return hsp['hit_strand'].lower() == 'plus'


def has_minus_hit_strand(hsp):
    return hsp['hit_strand'].lower() == 'minus'


class HitPositions:
    @staticmethod
    def covered_by(hsp):
        """Returns a list of all hit positions covered by the hsp."""
        hit_from = HitFrom.of(hsp)
        hit_to = HitTo.of(hsp)
        return [p for p in range(min(hit_from, hit_to), max(hit_from, hit_to) + 1)]


def is_to_ntomentosiformis_28S_rRNA(hit):
    """Returns True if the hit is to N. tomentosiformis 28S rRNA and False otherwise."""
    hit_name = hit['description'][0]['title'].lower()
    return 'nicotiana tomentosiformis' in hit_name and '28s ribosomal rna' in hit_name


def is_to_ntomentosiformis_18S_rRNA(hit):
    """Returns True if the hit is to N. tomentosiformis 18S rRNA and False otherwise."""
    hit_name = hit['description'][0]['title'].lower()
    return 'nicotiana tomentosiformis' in hit_name and '18s ribosomal rna' in hit_name


def is_to_cy1(hit):
    """Returns True if the hit is to CY1 gRNA and False otherwise."""
    return hit['description'][0]['title'].lower() == 'cy1'


for search in searches_with_a_hit:
    hit = Hit.for_(search)
    assert is_to_ntomentosiformis_28S_rRNA(hit)
print('All hits are to N. tomentosiformis 28S rRNA.')
print()


searches_for_plus_strand_reads = list(filter(is_for_plus_strand_read, searches_with_a_hit))
print('Searches for plus-strand reads: ' + str(len(searches_for_plus_strand_reads)))
print()

searches_for_minus_strand_reads = list(filter(is_for_minus_strand_read, searches_with_a_hit))
print('Searches for minus-strand reads: ' + str(len(searches_for_minus_strand_reads)))
print()

searches_for_plus_minus_hybrid_reads = list(filter(is_for_plus_minus_hybrid_read, searches_with_a_hit))
print('Searches for plus-minus hybrid reads: ' + str(len(searches_for_plus_minus_hybrid_reads)))
print()


searches_for_one_segment_reads = list(filter(
    lambda search : len(Hsps.of(Hit.for_(search))) == 1,
    searches_with_a_hit,
))
print('Searches for one segment reads: ' + str(len(searches_for_one_segment_reads)))
print()


searches_for_two_segment_reads = list(filter(
    lambda search : len(Hsps.of(Hit.for_(search))) == 2,
    searches_with_a_hit,
))
print('Searches for two segment reads: ' + str(len(searches_for_two_segment_reads)))
print()


#"""
fig, ax = plt.subplots()

plt.hist(
    [ReadLength.for_(search) for search in searches_with_a_hit],
    bins=225,
    color='black',
)

plt.show()
#"""


"""
fig, ax = plt.subplots()

random.shuffle(searches_for_one_segment_reads)

plt.hist(
    [p for search in searches_for_one_segment_reads[:5000] for p in UniqueHitPositions.covered_by(Hit.for_(search))],
    [i + 0.5 for i in range(0, 3381 + 1)],
    color='black',
)

ax.set_xticks([1, 500, 1000, 1500, 2000, 2500, 3000, 3381])

plt.xticks(rotation=90)

plt.show()
#"""


"""
fig, ax = plt.subplots()

plt.hist(
    [MinHitPosition.covered_by(Hit.for_(search)) for search in searches_with_a_hit],
    bins=100,
    color='black',
)

plt.show()
#"""


"""
fig, ax = plt.subplots()

plt.hist(
    [MaxHitPosition.covered_by(Hit.for_(search)) for search in searches_with_a_hit],
    bins=250,
    color='black',
)

plt.show()
#"""


"""
fig, ax = plt.subplots()

random.shuffle(searches_with_a_hit)

searches_with_a_hit.sort(key=functools.cmp_to_key(cmp_read_lengths), reverse=True)

i = 1

for search in searches_with_a_hit:
    hit = Hit.for_(search)
    hsps = Hsps.of(hit)

    for hsp in hsps:
        plt.plot([HitFrom.of(hsp), HitTo.of(hsp)], [i, i], color='black', alpha=0.001)

    i += 1

plt.show()
#"""
