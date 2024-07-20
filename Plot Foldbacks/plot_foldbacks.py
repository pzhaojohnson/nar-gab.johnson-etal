import json

import matplotlib.pyplot as plt

import functools


class Search:
    def __init__(self, data):
        self.data = data

    @property
    def read_length(self):
        return self.data['query_len']

    @property
    def hits(self):
        return [Hit(hit) for hit in self.data['hits']]

    @property
    def num_hits(self):
        return len(self.hits)

    @property
    def CY1_hit(self):
        CY1_hits = list(filter(lambda hit : hit.isto_CY1(), self.hits))
        assert len(CY1_hits) == 1
        return CY1_hits[0]

    def hasa_CY1_hit(self):
        return self.num_hits == 1 and self.hits[0].isto_CY1()

    def is_type_I_minus(self):
        assert self.num_hits == 1
        assert self.CY1_hit.is_minus()
        return self.CY1_hit.min_aligned_pos / self.read_length <= 0.05

    def is_type_II_minus(self):
        assert self.num_hits == 1
        assert self.CY1_hit.is_minus()
        return 0.4 <= self.CY1_hit.min_aligned_pos / self.read_length <= 0.53


def cmp_read_lengths(search1, search2):
    return search1.read_length - search2.read_length


class Hit:
    def __init__(self, data):
        self.data = data

    def isto_CY1(self):
        return self.data['description'][0]['title'].upper() == 'CY1'

    @property
    def hsps(self):
        return [Hsp(hsp) for hsp in self.data['hsps']]

    @property
    def num_hsps(self):
        return len(self.hsps)

    @property
    def hsps_sortedby_query_from(self):
        hsps = self.hsps
        hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
        return hsps

    @property
    def minus_hsps(self):
        return list(filter(lambda hsp : hsp.is_minus(), self.hsps))

    @property
    def minus_hsps_sortedby_query_from(self):
        minus_hsps = self.minus_hsps
        minus_hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
        return minus_hsps

    @property
    def min_aligned_pos(self):
        return self.hsps_sortedby_query_from[0].query_from

    def is_plus(self):
        return all([hsp.is_plus() for hsp in self.hsps])

    def is_minus(self):
        return all([hsp.is_minus() for hsp in self.hsps])

    def is_plus_minus(self):
        return len(self.hsps) > 1 and not self.is_plus() and not self.is_minus()


class Hsp:
    def __init__(self, data):
        self.data = data

    def is_plus(self):
        return self.data['hit_strand'] == 'Plus'

    def is_minus(self):
        return self.data['hit_strand'] == 'Minus'

    @property
    def query_from(self):
        return self.data['query_from']

    @property
    def query_to(self):
        return self.data['query_to']


def cmp_query_froms(hsp1, hsp2):
    return hsp1.query_from - hsp2.query_from


print()


blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'
print(f'{blast_output_file_path=}')
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output!')
print()


reports = [item['report'] for item in blast_output['BlastOutput2']]

all_searches = [Search(report['results']['search']) for report in reports]
print(f'{len(all_searches)=}')
print()


for search in all_searches:
    assert search.num_hits <= 1


CY1_searches = list(filter(lambda search : search.hasa_CY1_hit(), all_searches))
print(f'{len(CY1_searches)=}')
print()


plus_searches = list(filter(lambda search : search.CY1_hit.is_plus(), CY1_searches))
print(f'{len(plus_searches)=}')
print()

minus_searches = list(filter(lambda search : search.CY1_hit.is_minus(), CY1_searches))
print(f'{len(minus_searches)=}')
print()

plus_minus_searches = list(filter(lambda search : search.CY1_hit.is_plus_minus(), CY1_searches))
print(f'{len(plus_minus_searches)=}')
print()


type_I_minus_searches = list(filter(lambda search : search.is_type_I_minus(), minus_searches))
print(f'{len(type_I_minus_searches)=}')
print()

type_II_minus_searches = list(filter(lambda search : search.is_type_II_minus(), minus_searches))
print(f'{len(type_II_minus_searches)=}')
print()


foldback_searches = [
    *list(filter(lambda search : search.CY1_hit.num_hsps == 2, plus_minus_searches)),
    *type_II_minus_searches,
]

print(f'{len(foldback_searches)=}')
print()


fig, ax = plt.subplots()

foldback_searches.sort(key=functools.cmp_to_key(cmp_read_lengths), reverse=True)

y = 1

for search in foldback_searches:
    hit = search.CY1_hit
    first_minus_hsp = hit.minus_hsps_sortedby_query_from[0]
    plt.plot([1 - first_minus_hsp.query_from, search.read_length - first_minus_hsp.query_from], [y, y], color='gray')
    for hsp in hit.minus_hsps:
        plt.plot([hsp.query_from - first_minus_hsp.query_from, hsp.query_to - first_minus_hsp.query_from], [y, y], color='red')
    y += 1

ax.set_xticks([-2692, 0, 2691])

plt.show()
