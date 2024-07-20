import json

import matplotlib.pyplot as plt

import functools


class ReadID:
    @staticmethod
    def for_(search):
        """Returns the ID of the read for the search."""
        read_id = search['query_title']
        assert type(read_id) == str
        # all read IDs should be UUIDs
        assert len(read_id) == 36
        return read_id


class ReadLength:
    @staticmethod
    def for_(search):
        """Returns the length of the read for the search."""
        read_length = search['query_len']
        assert type(read_length) == int
        assert read_length > 0
        return read_length


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


def is_for_two_segment_plus_minus_hybrid_read(search):
    assert is_for_plus_minus_hybrid_read(search)
    hit = Hit.for_(search)
    hsps = Hsps.for_(hit)
    return len(hsps) == 2


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


class LastHsp:
    @staticmethod
    def in_(hit):
        """Returns the last hsp in the hit (by query-to position)."""
        hsps = Hsps.for_(hit)
        # don't forget to sort!
        hsps.sort(key=functools.cmp_to_key(cmp_query_tos))
        assert len(hsps) > 0
        return hsps[-1]


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


class HitFrom:
    @staticmethod
    def for_(hsp):
        """Returns the hit-from position for the hsp."""
        hit_from = hsp['hit_from']
        assert type(hit_from) == int
        return hit_from


class HitTo:
    @staticmethod
    def for_(hsp):
        """Returns the hit-to position for the hsp."""
        hit_to = hsp['hit_to']
        assert type(hit_to) == int
        return hit_to


def has_plus_hit_strand(hsp):
    return hsp['hit_strand'] == 'Plus'


def has_minus_hit_strand(hsp):
    return hsp['hit_strand'] == 'Minus'


class MountainPlotHeight:
    @staticmethod
    def of(folding):
        """Returns the mountain plot height of a folding (given in dot-bracket notation)."""
        mountain_plot_height = 0
        h = 0
        for c in folding:
            if c == '(':
                h += 1
                mountain_plot_height = max(mountain_plot_height, h)
            elif c == ')':
                h -= 1
            elif c != '.':
                raise Exception(f'Unrecognized character in dot-bracket notation: "{c}".')
        assert h == 0
        return mountain_plot_height


class NormalizedMountainPlotHeight:
    @staticmethod
    def of(folding):
        """Returns the normalized mountain plot height of the folding.

        (Normalized to the length of the folding.)
        """
        assert len(folding) > 0
        return MountainPlotHeight.of(folding) / len(folding)


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


searches_for_two_segment_plus_minus_hybrid_reads = list(filter(
    is_for_two_segment_plus_minus_hybrid_read,
    searches_for_plus_minus_hybrid_reads,
))
print(f'{len(searches_for_two_segment_plus_minus_hybrid_reads)=}')
print()


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


folded_cy1_reads_file_path = 'folded_cy1_reads_cy1_nb_6wpi_leaf.fasta'

print(f'{folded_cy1_reads_file_path=}')
print()


cy1_read_foldings = {}

with open(folded_cy1_reads_file_path, 'r') as f:
    lines = f.read().splitlines()
    for i in range(0, len(lines), 3):
        assert lines[i][0] == '>'
        read_id = lines[i][1:]
        # should be a UUID
        assert len(read_id) == 36
        read_seq = lines[i + 1]
        assert len(read_seq) > 0
        # in dot-bracket notation
        folding = lines[i + 2]
        # should have a trailing delta G value
        assert len(folding) > len(read_seq) + 3
        assert folding[len(read_seq)] == ' '
        # remove tailing delta G
        folding = folding[:len(read_seq)]
        assert len(folding) == len(read_seq)
        cy1_read_foldings[read_id] = folding

print(f'{len(cy1_read_foldings)=}')
print()


for search in searches_for_plus_strand_reads:
    assert ReadID.for_(search) in cy1_read_foldings

plus_strand_read_foldings = list(map(
    lambda search : cy1_read_foldings[ReadID.for_(search)],
    searches_for_plus_strand_reads,
))
print(f'{len(plus_strand_read_foldings)=}')
print()

plus_strand_read_norm_mtn_plot_heights = (
    [NormalizedMountainPlotHeight.of(folding) for folding in plus_strand_read_foldings]
)
print(f'{len(plus_strand_read_norm_mtn_plot_heights)=}')
print()


for search in searches_for_plus_minus_hybrid_reads:
    assert ReadID.for_(search) in cy1_read_foldings

plus_minus_hybrid_read_foldings = list(map(
    lambda search : cy1_read_foldings[ReadID.for_(search)],
    searches_for_plus_minus_hybrid_reads,
))
print(f'{len(plus_minus_hybrid_read_foldings)=}')
print()

plus_minus_hybrid_read_norm_mtn_plot_heights = (
    [NormalizedMountainPlotHeight.of(folding) for folding in plus_minus_hybrid_read_foldings]
)
print(f'{len(plus_minus_hybrid_read_norm_mtn_plot_heights)=}')
print()


for search in searches_for_two_segment_plus_minus_hybrid_reads:
    assert ReadID.for_(search) in cy1_read_foldings

two_segment_plus_minus_hybrid_read_foldings = (
    [cy1_read_foldings[ReadID.for_(search)] for search in searches_for_two_segment_plus_minus_hybrid_reads]
)
print(f'{len(two_segment_plus_minus_hybrid_read_foldings)=}')
print()

two_segment_plus_minus_hybrid_read_norm_mtn_plot_heights = (
    [NormalizedMountainPlotHeight.of(folding) for folding in two_segment_plus_minus_hybrid_read_foldings]
)
print(f'{len(two_segment_plus_minus_hybrid_read_norm_mtn_plot_heights)=}')
print()


for search in searches_for_type_i_minus_strand_reads:
    assert ReadID.for_(search) in cy1_read_foldings

type_i_minus_strand_read_foldings = (
    [cy1_read_foldings[ReadID.for_(search)] for search in searches_for_type_i_minus_strand_reads]
)
print(f'{len(type_i_minus_strand_read_foldings)=}')
print()

type_i_minus_strand_read_norm_mtn_plot_heights = (
    [NormalizedMountainPlotHeight.of(folding) for folding in type_i_minus_strand_read_foldings]
)
print(f'{len(type_i_minus_strand_read_norm_mtn_plot_heights)=}')
print()


for search in searches_for_type_ii_minus_strand_reads:
    assert ReadID.for_(search) in cy1_read_foldings

type_ii_minus_strand_read_foldings = (
    [cy1_read_foldings[ReadID.for_(search)] for search in searches_for_type_ii_minus_strand_reads]
)
print(f'{len(type_ii_minus_strand_read_foldings)=}')
print()

type_ii_minus_strand_read_norm_mtn_plot_heights = (
    [NormalizedMountainPlotHeight.of(folding) for folding in type_ii_minus_strand_read_foldings]
)
print(f'{len(type_ii_minus_strand_read_norm_mtn_plot_heights)=}')
print()


#"""
fig, ax = plt.subplots()

parts = plt.violinplot(
    #[*two_segment_plus_minus_hybrid_read_norm_mtn_plot_heights, *type_ii_minus_strand_read_norm_mtn_plot_heights],
    #plus_strand_read_norm_mtn_plot_heights,
    #type_i_minus_strand_read_norm_mtn_plot_heights,
    #two_segment_plus_minus_hybrid_read_norm_mtn_plot_heights,
    type_ii_minus_strand_read_norm_mtn_plot_heights,
    showextrema=False,
)

for pc in parts['bodies']:
    pc.set_facecolor('black')
    pc.set_edgecolor('black')
    pc.set_alpha(1)

ax.set_ylim([-0.1, 0.6])
ax.set_yticks([-0.1, -0.05, 0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6])

plt.show()
#"""
