import json

import matplotlib.pyplot as plt

import functools


class ReadLength:
    @staticmethod
    def for_(search):
        """Returns the length of the read for the search."""
        read_length = search['query_len']
        assert type(read_length) == int
        assert read_length > 0
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
        return len(Hits.for_(search))


def has_a_hit(search):
    return NumHits.for_(search) > 0


def has_exactly_one_hit(search):
    return NumHits.for_(search) == 1


def is_for_plus_strand_read(search):
    hit = Hit.for_(search)
    assert NumHsps.for_(hit) > 0
    return NumPlusHsps.for_(hit) == NumHsps.for_(hit)


def is_for_minus_strand_read(search):
    hit = Hit.for_(search)
    assert NumHsps.for_(hit) > 0
    return NumMinusHsps.for_(hit) == NumHsps.for_(hit)


def is_for_plus_minus_hybrid_read(search):
    hit = Hit.for_(search)
    assert NumHsps.for_(hit) > 0
    return not is_for_plus_strand_read(search) and not is_for_minus_strand_read(search)


def is_for_type_i_minus_strand_read(search):
    assert is_for_minus_strand_read(search)
    read_length = ReadLength.for_(search)
    hit = Hit.for_(search)
    first_hsp = FirstHsp.in_(hit)
    return QueryFrom.for_(first_hsp) / read_length <= 0.05


def is_for_type_ii_minus_strand_read(search):
    assert is_for_minus_strand_read(search)
    read_length = ReadLength.for_(search)
    hit = Hit.for_(search)
    first_hsp = FirstHsp.in_(hit)
    return 0.4 <= QueryFrom.for_(first_hsp) / read_length <= 0.53


class Hsps:
    @staticmethod
    def for_(hit):
        """Returns the hsps for the hit."""
        hsps = hit['hsps']
        assert type(hsps) == list
        # all hits should have at least one hsp
        assert len(hsps) > 0
        return hsps


class NumHsps:
    @staticmethod
    def for_(hit):
        """Returns the number of hsps for the hit."""
        return len(Hsps.for_(hit))


class PlusHsps:
    @staticmethod
    def for_(hit):
        """Returns the hsps for the hit with a plus hit-strand."""
        hsps = Hsps.for_(hit)
        return list(filter(has_plus_hit_strand, hsps))


class NumPlusHsps:
    @staticmethod
    def for_(hit):
        """Returns the number of hsps with a plus hit-strand for the hit."""
        return len(PlusHsps.for_(hit))


class MinusHsps:
    @staticmethod
    def for_(hit):
        """Returns the hsps for the hit with a minus hit-strand."""
        hsps = Hsps.for_(hit)
        return list(filter(has_minus_hit_strand, hsps))


class NumMinusHsps:
    @staticmethod
    def for_(hit):
        """Returns the number of hsps with a minus hit-strand for the hit."""
        return len(MinusHsps.for_(hit))


class FirstHsp:
    @staticmethod
    def in_(hit):
        """Returns the first hsp in the hit (by query-from position)."""
        hsps = Hsps.for_(hit)
        # don't forget to sort!
        hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
        assert len(hsps) > 0
        return hsps[0]


def is_to_cy1(hit):
    """Returns True if the hit is to CY1 and False otherwise."""
    return hit['description'][0]['title'] == 'CY1'


class QueryFrom:
    @staticmethod
    def for_(hsp):
        """Returns the query-from position for the hsp."""
        query_from = hsp['query_from']
        assert type(query_from) == int
        return query_from


def cmp_query_froms(hsp1, hsp2):
    return QueryFrom.for_(hsp1) - QueryFrom.for_(hsp2)


class QueryTo:
    @staticmethod
    def for_(hsp):
        """Returns the query-to position for the hsp."""
        query_to = hsp['query_to']
        assert type(query_to) == int
        return query_to


def cmp_query_tos(hsp1, hsp2):
    return QueryTo.for_(hsp1) - QueryTo.for_(hsp2)


def has_plus_hit_strand(hsp):
    return hsp['hit_strand'] == 'Plus'


def has_minus_hit_strand(hsp):
    return hsp['hit_strand'] == 'Minus'


print()


blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'

print(f'{blast_output_file_path=}')
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output.')
print()


reports = [item['report'] for item in blast_output['BlastOutput2']]
print(f'{len(reports)=}')
print()


searches = [report['results']['search'] for report in reports]
print(f'{len(searches)=}')
print()


searches_with_a_hit = list(filter(has_a_hit, searches))
print(f'{len(searches_with_a_hit)=}')
print()


for search in searches_with_a_hit:
    assert has_exactly_one_hit(search)
print('All searches have at most one hit.')
print()


for search in searches_with_a_hit:
    hit = Hit.for_(search)
    assert is_to_cy1(hit)
print('All hits are to CY1.')
print()


searches_for_plus_strand_reads = list(filter(is_for_plus_strand_read, searches_with_a_hit))
print(f'{len(searches_for_plus_strand_reads)=}')
print()

searches_for_minus_strand_reads = list(filter(is_for_minus_strand_read, searches_with_a_hit))
print(f'{len(searches_for_minus_strand_reads)=}')
print()

searches_for_plus_minus_hybrid_reads = list(filter(is_for_plus_minus_hybrid_read, searches_with_a_hit))
print(f'{len(searches_for_plus_minus_hybrid_reads)=}')
print()


assert len(searches_with_a_hit) == len([
    *searches_for_plus_strand_reads,
    *searches_for_minus_strand_reads,
    *searches_for_plus_minus_hybrid_reads,
])

for search in searches_for_minus_strand_reads:
    assert search not in searches_for_plus_strand_reads
    assert search not in searches_for_plus_minus_hybrid_reads

for search in searches_for_plus_minus_hybrid_reads:
    assert search not in searches_for_plus_strand_reads
    assert search not in searches_for_minus_strand_reads


searches_for_type_i_minus_strand_reads = list(filter(
    is_for_type_i_minus_strand_read,
    searches_for_minus_strand_reads,
))
print(f'{len(searches_for_type_i_minus_strand_reads)=}')
print()

searches_for_type_ii_minus_strand_reads = list(filter(
    is_for_type_ii_minus_strand_read,
    searches_for_minus_strand_reads,
))
print(f'{len(searches_for_type_ii_minus_strand_reads)=}')
print()


for search in searches_for_type_i_minus_strand_reads:
    assert search not in searches_for_type_ii_minus_strand_reads


#"""
fig, ax = plt.subplots()

searches_for_type_ii_minus_strand_reads.sort(key=functools.cmp_to_key(cmp_read_lengths), reverse=True)

i = 1

for search in searches_for_type_ii_minus_strand_reads:
    read_length = ReadLength.for_(search)

    hit = Hit.for_(search)

    first_hsp = FirstHsp.in_(hit)
    first_hsp_query_from = QueryFrom.for_(first_hsp)

    plt.plot([
        1 - first_hsp_query_from,
        read_length - first_hsp_query_from,
    ], [i, i], color='gray')

    for hsp in Hsps.for_(hit):
        plt.plot([
            QueryFrom.for_(hsp) - first_hsp_query_from,
            QueryTo.for_(hsp) - first_hsp_query_from,
        ], [i, i], color='red')

    i += 1

ax.set_xticks([-2692, 0, 2691])

plt.show()
#"""
