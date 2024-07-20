import json

import matplotlib.pyplot as plt

import functools


class Reports:
    @staticmethod
    def in_(blast_output):
        """Returns all reports in the BLAST output."""
        return [item['report'] for item in blast_output['BlastOutput2']]


class Searches:
    @staticmethod
    def in_(blast_output):
        """Returns all searches in the BLAST output."""
        reports = Reports.in_(blast_output)
        return [report['results']['search'] for report in reports]


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
    """Returns True if the search is for a plus-strand read and False otherwise."""
    hit = Hit.for_(search)
    return NumPlusHsps.for_(hit) == NumHsps.for_(hit)


def is_for_minus_strand_read(search):
    """Returns True if the search is for a minus-strand read and False otherwise."""
    hit = Hit.for_(search)
    return NumMinusHsps.for_(hit) == NumHsps.for_(hit)


def is_for_plus_minus_hybrid_read(search):
    """Returns True if the search is for a plus-minus hybrid read and False otherwise."""
    hit = Hit.for_(search)
    assert NumHsps.for_(hit) > 0
    return not is_for_plus_strand_read(search) and not is_for_minus_strand_read(search)


class FirstAlignedReadPosition:
    @staticmethod
    def for_(search):
        """Returns the first read position aligned in the single hit for the search.

        Raises if the search does not have exactly one hit.
        """
        hit = Hit.for_(search)
        return min(UniqueAlignedReadPositions.for_(hit))


def cmp_first_aligned_read_positions(search1, search2):
    return FirstAlignedReadPosition.for_(search1) - FirstAlignedReadPosition.for_(search2)


class Hsps:
    @staticmethod
    def of(hit):
        """Returns the hsps of the hit."""
        hsps = hit['hsps']
        assert type(hsps) == list
        return hsps


class FirstHsp:
    @staticmethod
    def in_(hit):
        """Returns the first hsp in the hit (by query-from position)."""
        hsps = Hsps.of(hit)
        hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
        assert len(hsps) > 0
        return hsps[0]


class NumHsps:
    @staticmethod
    def for_(hit):
        """Returns the number of hsps for the hit."""
        return len(Hsps.of(hit))


class PlusHsps:
    @staticmethod
    def for_(hit):
        """Returns the hsps with a plus hit-strand for the hit."""
        return list(filter(has_plus_hit_strand, Hsps.of(hit)))


class NumPlusHsps:
    @staticmethod
    def for_(hit):
        """Returns the number of hsps with a plus hit-strand for the hit."""
        return len(PlusHsps.for_(hit))


class MinusHsps:
    @staticmethod
    def for_(hit):
        """Returns the hsps with a minus hit-strand for the hit."""
        return list(filter(has_minus_hit_strand, Hsps.of(hit)))


class NumMinusHsps:
    @staticmethod
    def for_(hit):
        """Returns the number of hsps with a minus hit-strand for the hit."""
        return len(MinusHsps.for_(hit))


def is_to_cy1(hit):
    """Returns True if the hit is to CY1 and False otherwise."""
    return hit['description'][0]['title'] == 'CY1'


class QueryFrom:
    @staticmethod
    def of(hsp):
        """Returns the query-from position of the hsp."""
        query_from = hsp['query_from']
        assert type(query_from) == int
        return query_from


def cmp_query_froms(hsp1, hsp2):
    return QueryFrom.of(hsp1) - QueryFrom.of(hsp2)


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
    return hsp['hit_strand'] == 'Plus'


def has_minus_hit_strand(hsp):
    return hsp['hit_strand'] == 'Minus'


class AlignedReadPositions:
    @staticmethod
    def in_(hsp):
        """Returns the read positions that are aligned in the hsp."""
        query_from = QueryFrom.of(hsp)
        query_to = QueryTo.of(hsp)
        assert query_from < query_to
        return [p for p in range(query_from, query_to + 1)]


class UniqueAlignedReadPositions:
    @staticmethod
    def for_(hit):
        """Returns the set of unique read positions that are aligned in the hsps of the hit."""
        hsps = Hsps.of(hit)
        return set([p for hsp in hsps for p in AlignedReadPositions.in_(hsp)])


class NumUniqueAlignedReadPositions:
    @staticmethod
    def for_(hit):
        """Returns the number of unique read positions that are aligned in the hsps of the hit."""
        return len(UniqueAlignedReadPositions.for_(hit))


print()


blast_to_cy1_output_file_path = 'blast_to_cy1_output_cy1_nb_6wpi_leaf.json'
print(f'{blast_to_cy1_output_file_path=}')
print()


with open(blast_to_cy1_output_file_path, 'r') as f:
    blast_to_cy1_output = json.loads(f.read())
print('Successfully parsed BLAST-to-CY1 output.')
print()


blast_to_nb_transcripts_output_file_path = 'blast_to_nb_transcripts_output_cy1_nb_6wpi_leaf.json'
print(f'{blast_to_nb_transcripts_output_file_path=}')
print()


with open(blast_to_nb_transcripts_output_file_path, 'r') as f:
    blast_to_nb_transcripts_output = json.loads(f.read())
print('Successfully parsed BLAST-to-NB-transcripts output.')
print()


blast_to_nb_genome_output_file_path = 'blast_to_nb_genome_output_cy1_nb_6wpi_leaf.json'
print(f'{blast_to_nb_genome_output_file_path=}')
print()


with open(blast_to_nb_genome_output_file_path, 'r') as f:
    blast_to_nb_genome_output = json.loads(f.read())
print('Successfully parsed BLAST-to-NB-genome output.')
print()


searches_with_a_cy1_hit = list(filter(has_a_hit, Searches.in_(blast_to_cy1_output)))
assert all(has_exactly_one_hit(search) for search in searches_with_a_cy1_hit)
assert all(is_to_cy1(Hit.for_(search)) for search in searches_with_a_cy1_hit)
print(f'{len(searches_with_a_cy1_hit)=}')
print()


cy1_read_ids = set([ReadID.for_(search) for search in searches_with_a_cy1_hit])
print(f'{len(cy1_read_ids)=}')
print()


searches_with_an_nb_transcripts_hit = list(filter(has_a_hit, Searches.in_(blast_to_nb_transcripts_output)))
print(f'{len(searches_with_an_nb_transcripts_hit)=}')
print()


nb_transcript_read_ids = list(map(lambda search : ReadID.for_(search), searches_with_an_nb_transcripts_hit))
cy1_nb_transcript_chimeric_read_ids = list(filter(lambda read_id : read_id in cy1_read_ids, nb_transcript_read_ids))
print(f'{len(cy1_nb_transcript_chimeric_read_ids)=}')
print()


searches_with_an_nb_genome_hit = list(filter(has_a_hit, Searches.in_(blast_to_nb_genome_output)))
print(f'{len(searches_with_an_nb_genome_hit)=}')
print()


nb_genome_read_ids = list(map(lambda search : ReadID.for_(search), searches_with_an_nb_genome_hit))
cy1_nb_genome_chimeric_read_ids = list(filter(lambda read_id : read_id in cy1_read_ids, nb_genome_read_ids))
print(f'{len(cy1_nb_genome_chimeric_read_ids)=}')
print()


unique_chimeric_read_ids = set([*cy1_nb_transcript_chimeric_read_ids, *cy1_nb_genome_chimeric_read_ids])
print(f'{len(unique_chimeric_read_ids)=}')
print()


searches_to_cy1_for_cy1_nb_transcript_chimeric_reads = list(map(
    lambda read_id : next(search for search in searches_with_a_cy1_hit if ReadID.for_(search) == read_id),
    cy1_nb_transcript_chimeric_read_ids,
))
print(f'{len(searches_to_cy1_for_cy1_nb_transcript_chimeric_reads)=}')
print()


searches_to_cy1_for_cy1_nb_genome_chimeric_reads = list(map(
    lambda read_id : next(search for search in searches_with_a_cy1_hit if ReadID.for_(search) == read_id),
    cy1_nb_genome_chimeric_read_ids,
))
print(f'{len(searches_to_cy1_for_cy1_nb_genome_chimeric_reads)=}')
print()


num_plus_strand_cy1_nb_genome_chimeric_reads = len(list(filter(is_for_plus_strand_read, searches_to_cy1_for_cy1_nb_genome_chimeric_reads)))
print(f'{num_plus_strand_cy1_nb_genome_chimeric_reads=}')
print()

num_minus_strand_cy1_nb_genome_chimeric_reads = len(list(filter(is_for_minus_strand_read, searches_to_cy1_for_cy1_nb_genome_chimeric_reads)))
print(f'{num_minus_strand_cy1_nb_genome_chimeric_reads=}')
print()

num_plus_minus_hybrid_cy1_nb_genome_chimeric_reads = len(list(filter(is_for_plus_minus_hybrid_read, searches_to_cy1_for_cy1_nb_genome_chimeric_reads)))
print(f'{num_plus_minus_hybrid_cy1_nb_genome_chimeric_reads=}')
print()

assert len(searches_to_cy1_for_cy1_nb_genome_chimeric_reads) == (
    num_plus_strand_cy1_nb_genome_chimeric_reads
    + num_minus_strand_cy1_nb_genome_chimeric_reads
    + num_plus_minus_hybrid_cy1_nb_genome_chimeric_reads
)


for search in searches_to_cy1_for_cy1_nb_genome_chimeric_reads:
    print(ReadID.for_(search))
print()


cy1_searches_to_plot = searches_to_cy1_for_cy1_nb_genome_chimeric_reads
print(f'{len(cy1_searches_to_plot)=}')
print()

cy1_searches_to_plot.sort(key=functools.cmp_to_key(cmp_read_lengths), reverse=True)
cy1_searches_to_plot.sort(key=functools.cmp_to_key(cmp_first_aligned_read_positions))

"""
fig, ax = plt.subplots()

i = 1

for cy1_search in cy1_searches_to_plot:
    cy1_hit = Hit.for_(cy1_search)

    read_length = ReadLength.for_(cy1_search)

    plt.plot([
        1 - QueryFrom.of(FirstHsp.in_(cy1_hit)),
        read_length - QueryFrom.of(FirstHsp.in_(cy1_hit)),
    ], [i, i], color='gray', linewidth=0.5)

    for hsp in Hsps.of(cy1_hit):
        plt.plot([
            QueryFrom.of(hsp) - QueryFrom.of(FirstHsp.in_(cy1_hit)),
            QueryTo.of(hsp) - QueryFrom.of(FirstHsp.in_(cy1_hit)),
        ], [i, i], color='red', linewidth=0.5)

    i += 1

plt.show()
#"""


"""
fig, ax = plt.subplots()

i = 1

for cy1_search in cy1_searches_to_plot:
    for hsp in Hsps.of(Hit.for_(cy1_search)):
        plt.plot([HitFrom.of(hsp), HitTo.of(hsp)], [i, i], color='red', linewidth=0.5)
    i += 1

ax.set_xticks([-299, 1, 281, 671, 2420, 2692, 2992])

plt.show()
#"""
