import json

import matplotlib.pyplot as plt

import functools


print()


blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'

print('BLAST output file path: ' + blast_output_file_path)
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output.')
print()


class Reports:
    @staticmethod
    def in_blast_output(blast_output):
        """Returns all reports in the BLAST output."""
        assert type(blast_output['BlastOutput2']) == list
        return [item['report'] for item in blast_output['BlastOutput2']]


class Searches:
    @staticmethod
    def in_blast_output(blast_output):
        """Returns all searches in the BLAST output."""
        reports = Reports.in_blast_output(blast_output)
        return [report['results']['search'] for report in reports]


class ReadLength:
    @staticmethod
    def for_search(search):
        """Returns the read length for the search."""
        read_length = search['query_len']
        assert type(read_length) == int
        return read_length


def is_for_plus_strand_read(search):
    """Returns True if the search is for a plus-strand read and False otherwise."""
    hit = Hit.for_search(search)
    return len(PlusHsps.for_hit(hit)) == len(Hsps.for_hit(hit))


def is_for_minus_strand_read(search):
    """Returns True if the search is for a minus-strand read and False otherwise."""
    hit = Hit.for_search(search)
    return len(MinusHsps.for_hit(hit)) == len(Hsps.for_hit(hit))


def is_for_type_i_minus_strand_read(search):
    hit = Hit.for_search(search)
    first_hsp = FirstHsp.of(hit)
    query_from = QueryFrom.for_hsp(first_hsp)
    read_length = ReadLength.for_search(search)
    return is_for_minus_strand_read(search) and query_from / read_length <= 0.05


def is_for_type_ii_minus_strand_read(search):
    hit = Hit.for_search(search)
    first_hsp = FirstHsp.of(hit)
    query_from = QueryFrom.for_hsp(first_hsp)
    read_length = ReadLength.for_search(search)
    return is_for_minus_strand_read(search) and 0.4 <= query_from / read_length <= 0.53


class Hits:
    @staticmethod
    def for_search(search):
        """Returns the hits for the search."""
        hits = search['hits']
        assert type(hits) == list
        return hits


class Hit:
    @staticmethod
    def for_search(search):
        """Returns the single hit for the search.

        Raises if the search does not have exactly one hit.
        """
        hits = Hits.for_search(search)
        assert len(hits) == 1
        return hits[0]


class Hsps:
    @staticmethod
    def for_hit(hit):
        """Returns the hsps for the hit."""
        hsps = hit['hsps']
        assert type(hsps) == list
        return hsps


class PlusHsps:
    @staticmethod
    def for_hit(hit):
        """Returns hsps with a plus hit strand for the hit."""
        hsps = Hsps.for_hit(hit)
        return list(filter(has_plus_hit_strand, hsps))


class MinusHsps:
    @staticmethod
    def for_hit(hit):
        """Returns hsps with a minus hit strand for the hit."""
        hsps = Hsps.for_hit(hit)
        return list(filter(has_minus_hit_strand, hsps))


class FirstHsp:
    @staticmethod
    def of(hit):
        """Returns the first hsp of the hit (i.e., the hsp with the smallest query-from position)."""
        hsps = Hsps.for_hit(hit)
        assert len(hsps) > 0
        # don't forget to sort!
        hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
        return hsps[0]


class LastHsp:
    @staticmethod
    def of(hit):
        """Returns the last hsp of the hit (i.e., the hsp with the largest query-from position)."""
        hsps = Hsps.for_hit(hit)
        assert len(hsps) > 0
        # don't forget to sort!
        hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
        return hsps[-1]


class QueryFrom:
    @staticmethod
    def for_hsp(hsp):
        """Returns the query-from position for the hsp."""
        query_from = hsp['query_from']
        assert type(query_from) == int
        return query_from


def cmp_query_froms(hsp1, hsp2):
    return QueryFrom.for_hsp(hsp1) - QueryFrom.for_hsp(hsp2)


class HitFrom:
    @staticmethod
    def for_hsp(hsp):
        """Returns the hit-from position for the hsp."""
        hit_from = hsp['hit_from']
        assert type(hit_from) == int
        return hit_from


class HitTo:
    @staticmethod
    def for_hsp(hsp):
        """Returns the hit-to position for the hsp."""
        hit_to = hsp['hit_to']
        assert type(hit_to) == int
        return hit_to


def has_plus_hit_strand(hsp):
    return hsp['hit_strand'].lower() == 'plus'


def has_minus_hit_strand(hsp):
    return hsp['hit_strand'].lower() == 'minus'


all_searches = Searches.in_blast_output(blast_output)
print('All searches: ' + str(len(all_searches)))
print()


searches_with_a_hit = list(filter(
    lambda search : len(Hits.for_search(search)) > 0,
    all_searches,
))
print('Searches with a hit: ' + str(len(searches_with_a_hit)))
print()


for search in searches_with_a_hit:
    assert len(Hits.for_search(search)) == 1
print('All searches have at most one hit.')
print()


for search in searches_with_a_hit:
    hit = Hit.for_search(search)
    assert hit['description'][0]['title'].lower() == 'cy1'
print('All hits are CY1.')
print()


searches_for_plus_strand_reads = list(filter(is_for_plus_strand_read, searches_with_a_hit))
print('Searches for plus-strand reads: ' + str(len(searches_for_plus_strand_reads)))
print()

searches_for_type_i_minus_strand_reads = list(filter(is_for_type_i_minus_strand_read, searches_with_a_hit))
print('Searches for type I minus-strand reads: ' + str(len(searches_for_type_i_minus_strand_reads)))
print()

searches_for_type_ii_minus_strand_reads = list(filter(is_for_type_ii_minus_strand_read, searches_with_a_hit))
print('Searches for type II minus-strand reads: ' + str(len(searches_for_type_ii_minus_strand_reads)))
print()


fig, ax = plt.subplots()

plt.hist(
    [HitFrom.for_hsp(FirstHsp.of(Hit.for_search(search))) for search in searches_for_plus_strand_reads],
    #[HitTo.for_hsp(LastHsp.of(Hit.for_search(search))) for search in searches_for_plus_strand_reads],
    #[HitFrom.for_hsp(FirstHsp.of(Hit.for_search(search))) for search in searches_for_type_i_minus_strand_reads],
    #[HitTo.for_hsp(LastHsp.of(Hit.for_search(search))) for search in searches_for_type_i_minus_strand_reads],
    bins=135,
    color='black',
)

ax.set_xlim([1 - 150, 2692 + 150])

ax.set_xticks([1, 74, 281, 455, 671, 1086, 2420, 2631, 2692])

plt.xticks(rotation=90)

plt.show()
