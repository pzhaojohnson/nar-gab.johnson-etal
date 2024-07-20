import json

import matplotlib.pyplot as plt


print()


class Reports:
    @staticmethod
    def in_(blast_output):
        return [item['report'] for item in blast_output['BlastOutput2']]


class Searches:
    @staticmethod
    def in_(blast_output):
        return [report['results']['search'] for report in Reports.in_(blast_output)]


class ReadLength:
    @staticmethod
    def for_(search):
        read_length = search['query_len']
        assert type(read_length) == int
        assert read_length > 0
        return read_length


class Hits:
    @staticmethod
    def for_(search):
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


class Hsps:
    @staticmethod
    def for_(hit):
        """Returns the hsps for the hit."""
        hsps = hit['hsps']
        assert type(hsps) == list
        assert len(hsps) > 0
        return hsps


class MaxQueryFrom:
    @staticmethod
    def for_(hsps):
        """Returns the maximum query-from position for the hsps."""
        return max(list(map(lambda hsp : QueryFrom.for_(hsp), hsps)))


class MaxQueryTo:
    @staticmethod
    def for_(hsps):
        """Returns the maximum query-to position for the hsps."""
        return max(list(map(lambda hsp : QueryTo.for_(hsp), hsps)))


class MaxQueryPosition:
    @staticmethod
    def for_(hsps):
        """Returns the maximum query position for the hsps."""
        return max([MaxQueryFrom.for_(hsps), MaxQueryTo.for_(hsps)])


class QueryFrom:
    @staticmethod
    def for_(hsp):
        """Returns the query-from position for the hsp."""
        query_from = hsp['query_from']
        assert type(query_from) == int
        return query_from


class QueryTo:
    @staticmethod
    def for_(hsp):
        """Returns the query-to position for the hsp."""
        query_to = hsp['query_to']
        assert type(query_to) == int
        return query_to


blast_output_file_path = 'blast_to_rubisco_large_output_cy1_nb_2wpi_leaf.json'
print(f'{blast_output_file_path}')
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output.')
print()


all_searches = Searches.in_(blast_output)
print(f'{len(all_searches)=}')
print()


searches_with_a_hit = list(filter(lambda search : len(Hits.for_(search)) > 0, all_searches))
print(f'{len(searches_with_a_hit)=}')
print()


for search in searches_with_a_hit:
    assert len(Hits.for_(search)) == 1
print('All searches have at most one hit.')
print()


fig, ax = plt.subplots()

plt.hist(
    [ReadLength.for_(search) - MaxQueryPosition.for_(Hsps.for_(Hit.for_(search))) for search in searches_with_a_hit],
    color='black',
    bins=150,
)

plt.show()
