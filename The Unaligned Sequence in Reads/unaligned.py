import json

import matplotlib.pyplot as plt


print()


class Search:
    """A BLAST search."""

    def __init__(self, raw):
        self.raw = raw

    @property
    def read_length(self):
        """Returns the length of the read for the search."""
        read_length = self.raw['query_len']
        assert type(read_length) == int
        assert read_length > 0
        return read_length

    @property
    def hits(self):
        """Returns the hits for the search."""
        hits = self.raw['hits']
        assert type(hits) == list
        return list(map(lambda hit : Hit(hit), hits))

    @property
    def hit(self):
        """Returns the single hit for the search.

        Raises if the search does not have exactly one hit.
        """
        assert len(self.hits) == 1
        return self.hits[0]

    @property
    def num_hits(self):
        """The number of hits for the search."""
        return len(self.hits)

    @property
    def unique_aligned_read_positions(self):
        """The set of unique read positions that aligned in the single hit for the search.

        Raises if the search does not have exactly one hit.
        """
        assert self.num_hits == 1
        return self.hit.unique_aligned_query_positions

    @property
    def num_unique_aligned_read_positions(self):
        """The number of unique read positions that aligned in the single hit for the search.

        Raises if the search does not have exactly one hit.
        """
        assert self.num_hits == 1
        return len(self.unique_aligned_read_positions)

    @property
    def num_unique_unaligned_read_positions(self):
        """The number of unique read positions that did not align in the single hit for the search.

        Raises if the search does not have exactly one hit.
        """
        assert self.num_hits == 1
        return self.read_length - self.num_unique_aligned_read_positions


def cmp_read_lengths(search1, search2):
    return search1.read_length - search2.read_length


class Hit:
    """A hit for a search."""

    def __init__(self, raw):
        self.raw = raw

    @property
    def hsps(self):
        """The hsps for the hit."""
        hsps = self.raw['hsps']
        assert type(hsps) == list
        assert len(hsps) > 0
        return list(map(lambda hsp : Hsp(hsp), hsps))

    @property
    def unique_aligned_query_positions(self):
        """The set of unique query positions that are aligned in the hsps for the hit."""
        return set([p for hsp in self.hsps for p in hsp.unique_aligned_query_positions])


class Hsp:
    """An hsp for a hit."""

    def __init__(self, raw):
        self.raw = raw

    @property
    def query_from(self):
        """The query-from position for the hsp."""
        query_from = self.raw['query_from']
        assert type(query_from) == int
        assert query_from >= 1
        return query_from

    @property
    def query_to(self):
        """The query-to position for the hsp."""
        query_to = self.raw['query_to']
        assert type(query_to) == int
        assert query_to >= 1
        return query_to

    @property
    def min_query_position(self):
        """The smallest query position aligned in the hsp."""
        return min(self.query_from, self.query_to)

    @property
    def max_query_position(self):
        """The largest query position aligned in the hsp."""
        return max(self.query_from, self.query_to)

    @property
    def aligned_query_positions(self):
        """The query positions that are aligned in the hsp."""
        return list(range(self.min_query_position, self.max_query_position, 1))

    @property
    def unique_aligned_query_positions(self):
        """The set of unique query positions that are aligned in the hsp."""
        return set(self.aligned_query_positions)


print()


#blast_output_file_path = 'blast_output_ivt_cy1_gRNA.json'

blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'

print('BLAST output file path: ' + blast_output_file_path)
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output.')
print()


reports = [thing['report'] for thing in blast_output['BlastOutput2']]
print(f'{len(reports)=}')
print()


all_searches = [Search(report['results']['search']) for report in reports]
print(f'{len(all_searches)=}')
print()


searches_with_a_hit = list(filter(lambda search : search.num_hits > 0, all_searches))
print(f'{len(searches_with_a_hit)=}')
print()


assert all(map(lambda search : search.num_hits == 1, searches_with_a_hit))
print('All searches have at most one hit.')
print()


fig, ax = plt.subplots()

plt.hist(
    [search.num_unique_unaligned_read_positions for search in searches_with_a_hit],
    color='black',
    bins=500,
)

plt.show()
