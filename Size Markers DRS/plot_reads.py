import json

from matplotlib import pyplot as plt

import functools

import random


class Search:
    def __init__(self, data):
        self.data = data

    @property
    def hits(self):
        return [Hit(hit) for hit in self.data['hits']]

    @property
    def num_hits(self):
        return len(self.hits) == 1

    def has_CY1_hit(self):
        return len(list(filter(lambda hit : hit.isto_CY1(), self.hits))) > 0

    @property
    def CY1_hit(self):
        assert self.num_hits == 1
        assert self.hits[0].isto_CY1()
        return self.hits[0]

    @property
    def unique_covered_positions(self):
        return set([p for hsp in self.CY1_hit.hsps for p in hsp.covered_positions])


def cmp_num_unique_covered_positions(search1, search2):
    return len(search1.unique_covered_positions) - len(search2.unique_covered_positions)


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
    def hsp(self):
        assert self.num_hsps == 1
        return self.hsps[0]

    @property
    def first_hsp(self):
        assert self.num_hsps > 0
        return self.hsps[0]

    @property
    def second_hsp(self):
        assert self.num_hsps > 1
        return self.hsps[1]

    def matches_F281(self):
        return self.num_hsps == 1 and self.hsp.matches_F281()

    def matches_F671(self):
        return self.num_hsps == 1 and self.hsp.matches_F671()

    def matches_DRNA(self):
        return self.num_hsps == 2 \
            and self.first_hsp.matches_F671() \
            and self.second_hsp.matches_3prime_251()

    def matches_F1600(self):
        return self.num_hsps == 1 and self.hsp.matches_F1600()

    def matches_CY1_gRNA(self):
        return self.num_hsps == 1 and self.hsp.matches_CY1_gRNA()

    @property
    def plotted_color(self):
        if self.matches_F281():
            return 'dodgerblue'
        elif self.matches_F671():
            return 'darkorange'
        elif self.matches_DRNA():
            return 'lime'
        elif self.matches_F1600():
            return 'turquoise'
        elif self.matches_CY1_gRNA():
            return 'red'
        else:
            return 'black'


def are_within(n1, n2, max_diff):
    return abs(n1 - n2) <= max_diff


class Hsp:
    def __init__(self, data):
        self.data = data

    @property
    def hit_from(self):
        return self.data['hit_from']

    @property
    def hit_to(self):
        return self.data['hit_to']

    @property
    def covered_positions(self):
        return [p for p in range(self.hit_from, self.hit_to + 1)]

    def is_5prime_coterminal(self):
        return are_within(self.hit_from, 1, 30)

    def is_3prime_coterminal(self):
        return are_within(self.hit_to, 2692, 10)

    def matches_F281(self):
        return self.is_5prime_coterminal() and are_within(self.hit_to, 281, 10)

    def matches_F671(self):
        return self.is_5prime_coterminal() and are_within(self.hit_to, 671, 10)

    def matches_3prime_251(self):
        return self.is_3prime_coterminal() and are_within(self.hit_from, 2442, 10)

    def matches_F1600(self):
        return self.is_5prime_coterminal() and are_within(self.hit_to, 1600, 10)

    def matches_CY1_gRNA(self):
        return self.is_5prime_coterminal() and self.is_3prime_coterminal()


print()


blast_output_file_path = 'blastn_IVT_size_markers_HAC.json'
print(f'{blast_output_file_path=}')
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Parsed BLAST output.')
print()


reports = [item['report'] for item in blast_output['BlastOutput2']]
print(f'{len(reports)=}')
print()


searches = [Search(report['results']['search']) for report in reports]
print(f'{len(searches)=}')
print()


CY1_searches = list(filter(lambda search : search.has_CY1_hit(), searches))
print(f'{len(CY1_searches)=}')
print()


F281_searches = list(filter(lambda search : search.CY1_hit.matches_F281(), CY1_searches))
print(f'{100 * len(F281_searches) / len(CY1_searches)=}')
print()

F671_searches = list(filter(lambda search : search.CY1_hit.matches_F671(), CY1_searches))
print(f'{100 * len(F671_searches) / len(CY1_searches)=}')
print()

DRNA_searches = list(filter(lambda search : search.CY1_hit.matches_DRNA(), CY1_searches))
print(f'{100 * len(DRNA_searches) / len(CY1_searches)=}')
print()

F1600_searches = list(filter(lambda search : search.CY1_hit.matches_F1600(), CY1_searches))
print(f'{100 * len(F1600_searches) / len(CY1_searches)=}')
print()

CY1_gRNA_searches = list(filter(lambda search : search.CY1_hit.matches_CY1_gRNA(), CY1_searches))
print(f'{100 * len(CY1_gRNA_searches) / len(CY1_searches)=}')
print()


fig, ax = plt.subplots()

searches_to_plot = random.sample(CY1_searches, 10000)

searches_to_plot.sort(key=functools.cmp_to_key(cmp_num_unique_covered_positions), reverse=True)

i = 1

for search in searches_to_plot:
    color = search.CY1_hit.plotted_color
    alpha = 0.05 if color == 'black' else 1
    for hsp in search.CY1_hit.hsps:
        plt.plot([hsp.hit_from, hsp.hit_to], [i, i], color=color, alpha=alpha)
    i += 1

ax.set_xticks([1, 281, 442, 671, 1600, 2442, 2692])

plt.show()
