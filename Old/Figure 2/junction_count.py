import json

import functools


print()


blast_output_file_path = 'blast_output_cy1_nb_2wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'

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
        return [item['report'] for item in blast_output['BlastOutput2']]


class Searches:
    @staticmethod
    def in_blast_output(blast_output):
        """Returns all searches in the BLAST output."""
        reports = Reports.in_blast_output(blast_output)
        return [report['results']['search'] for report in reports]


class Hits:
    @staticmethod
    def for_search(search):
        """Returns all hits for the search."""
        hits = search['hits']
        assert type(hits) == list
        return hits


class NumHits:
    @staticmethod
    def for_search(search):
        """Returns the number of hits for the search."""
        return len(Hits.for_search(search))


def has_a_hit(search):
    """Returns True if the search has at least one and False otherwise."""
    return NumHits.for_search(search) > 0


def has_exactly_one_hit(search):
    """Returns True if the search has exactly one hit and False otherwise."""
    return NumHits.for_search(search) == 1


def has_a_cy1_hit(search):
    """Returns True if the search has at least one CY1 hit and False otherwise."""
    hits = Hits.for_search(search)
    cy1_hits = list(filter(is_cy1_hit, hits))
    return len(cy1_hits) > 0


class Hit:
    @staticmethod
    def of(search):
        """Returns the single hit of the search.

        Raises if the search does not have exactly one hit."""
        hits = Hits.for_search(search)
        assert len(hits) == 1
        return hits[0]


def is_cy1_hit(hit):
    """Returns True if the hit is a CY1 hit and False otherwise."""
    return hit['description'][0]['title'].lower() == 'cy1'


class Hsps:
    @staticmethod
    def of(hit):
        """Returns the hsps of the hit."""
        hsps = hit['hsps']
        assert type(hsps) == list
        return hsps


class QueryFrom:
    @staticmethod
    def of(hsp):
        """Returns the query-from position of the hsp."""
        query_from = hsp['query_from']
        assert type(query_from) == int
        return query_from


def cmp_query_froms(hsp1, hsp2):
    return QueryFrom.of(hsp1) - QueryFrom.of(hsp2)


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


all_searches = Searches.in_blast_output(blast_output)
print('All searches: ' + str(len(all_searches)))
print()

searches_with_a_hit = list(filter(has_a_hit, all_searches))
print('Searches with a hit: ' + str(len(searches_with_a_hit)))
print()

for search in searches_with_a_hit:
    assert has_exactly_one_hit(search)
print('All searches have at most one hit.')
print()

for search in searches_with_a_hit:
    assert has_a_cy1_hit(search)
print('All hits are CY1.')
print()


junctions = []

for search in searches_with_a_hit:
    hit = Hit.of(search)

    # don't forget to sort by read position!
    hsps = Hsps.of(hit)
    hsps.sort(key=functools.cmp_to_key(cmp_query_froms))

    for i in range(len(hsps) - 1):
        hsp1 = hsps[i]
        hsp2 = hsps[i + 1]
        junctions.append(str(HitTo.of(hsp1)) + ' / ' + str(HitFrom.of(hsp2)))

print('Junctions found: ' + str(len(junctions)))
print()


unique_junctions = set(junctions)
print('Unique junctions: ' + str(len(unique_junctions)))
print()


junction_counts = {}

for junc in unique_junctions:
    occurrences = list(filter(lambda j : j == junc, junctions))
    junction_counts[junc] = len(occurrences)

print('Junction counts:')
print()

for junc, ct in junction_counts.items():
    print(junc + '\t' + str(ct))
print()
