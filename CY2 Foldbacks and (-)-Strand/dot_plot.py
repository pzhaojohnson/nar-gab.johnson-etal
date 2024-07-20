import json

import matplotlib.pyplot as plt


print()


class Search:
    def __init__(self, data):
        self.data = data

    @property
    def hits(self):
        return [Hit(hit) for hit in self.data['hits']]

    @property
    def CY2_hit(self):
        CY2_hits = list(filter(lambda hit : hit.isto_CY2(), self.hits))
        assert len(CY2_hits) == 1
        return CY2_hits[0]

    @property
    def num_hits(self):
        return len(self.hits)

    def has_exactly_two_segments(self):
        return len(self.CY2_hit.hsps) == 2


class Hit:
    def __init__(self, data):
        self.data = data

    def isto_CY2(self):
        return self.data['description'][0]['title'].upper() == 'CY2'

    @property
    def hsps(self):
        return [Hsp(hsp) for hsp in self.data['hsps']]

    @property
    def plus_hsps(self):
        return list(filter(lambda hsp : hsp.is_plus(), self.hsps))

    def is_plus_minus(self):
        return len(self.plus_hsps) < len(self.hsps) and len(self.plus_hsps) > 0


class Hsp:
    def __init__(self, data):
        self.data = data

    def is_plus(self):
        return self.data['hit_strand'] == 'Plus'

    def is_minus(self):
        return self.data['hit_strand'] == 'Minus'

    @property
    def hit_from(self):
        return self.data['hit_from']

    @property
    def hit_to(self):
        return self.data['hit_to']

    @property
    def query_from(self):
        return self.data['query_from']

    @property
    def query_to(self):
        return self.data['query_to']


blast_output_file_path = 'blast_output_cy2_nb_14wpi_leaf.json'

with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())


reports = [item['report'] for item in blast_output['BlastOutput2']]

all_searches = [Search(report['results']['search']) for report in reports]
print(f'{len(all_searches)=}')
print()


CY2_searches = list(filter(lambda search : search.num_hits == 1 and search.hits[0].isto_CY2(), all_searches))
print(f'{len(CY2_searches)=}')
print()


plus_minus_searches = list(filter(lambda search : search.CY2_hit.is_plus_minus(), CY2_searches))
print(f'{len(plus_minus_searches)=}')
print()

two_segment_plus_minus_searches = list(filter(lambda search : search.has_exactly_two_segments(), plus_minus_searches))
print(f'{len(two_segment_plus_minus_searches)=}')
print()


"""
fig, ax = plt.subplots()

for search in plus_minus_searches:
    hit = search.CY2_hit
    for hsp in hit.hsps:
        plt.plot([hsp.query_from, hsp.query_to], [hsp.hit_from, hsp.hit_to], color='black', alpha=0.2)

#ax.set_xticks([1, 281, 442, 671, 2983, 2 * 2983, (2 * 2983) + 300])
#ax.set_yticks([1, 281, 442, 671, 2068, 2710, 2983])

plt.show()
#"""


#"""
fig, ax = plt.subplots()

for search in two_segment_plus_minus_searches:
    hit = search.CY2_hit
    hsps = hit.hsps
    assert len(hsps) == 2
    # sort hsps
    hsps = hsps if hsps[0].query_from < hsps[1].query_from else [hsps[1], hsps[0]]
    plt.plot(
        [hsps[0].query_from - hsps[1].query_from, hsps[0].query_to - hsps[1].query_from],
        [hsps[0].hit_from - hsps[1].hit_from, hsps[0].hit_to - hsps[1].hit_from],
        color='black',
        alpha=0.2,
    )
    plt.plot(
        [hsps[1].query_from - hsps[1].query_from, hsps[1].query_to - hsps[1].query_from],
        [hsps[1].hit_from - hsps[1].hit_from, hsps[1].hit_to - hsps[1].hit_from],
        color='red',
        alpha=0.2,
    )

ax.set_xticks([-2983, -1, 2982, 2982 + 300])
ax.set_yticks([-2982, 0])

plt.show()
#"""

