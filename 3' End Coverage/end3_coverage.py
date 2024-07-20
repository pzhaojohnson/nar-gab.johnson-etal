import json

import matplotlib.pyplot as plt


print()


class Search:
    def __init__(self, data):
        self.data = data

    @property
    def read_length(self):
        return self.data['query_len']

    @property
    def hit(self):
        hits = self.data['hits']
        assert len(hits) == 1
        return Hit(hits[0])

    @property
    def num_hits(self):
        return len(self.data['hits'])

    def is_viral(self):
        return self.hit.data['description'][0]['title'].upper() in ['CY1', 'CY2']

    def is_plus(self):
        return all([hsp.is_plus() for hsp in self.hit.hsps])

    def is_minus(self):
        return all([hsp.is_minus() for hsp in self.hit.hsps])

    def is_plus_minus_hybrid(self):
        return not self.is_plus() and not self.is_minus()

    def is_type_I_minus(self):
        return self.is_minus() \
            and self.hit.min_queryp / self.read_length <= 0.05

    def is_foldback(self):
        if self.is_plus_minus_hybrid():
            return self.hit.num_hsps == 2
        elif self.is_minus():
            return 0.4 <= self.hit.min_queryp / self.read_length <= 0.53
        else:
            return False


class Hit:
    def __init__(self, data):
        self.data = data

    @property
    def hsps(self):
        return [Hsp(hsp) for hsp in self.data['hsps']]

    @property
    def num_hsps(self):
        return len(self.hsps)

    @property
    def plus_hsps(self):
        return list(filter(lambda hsp : hsp.is_plus(), self.hsps))

    @property
    def num_plus_hsps(self):
        return len(self.plus_hsps)

    @property
    def minus_hsps(self):
        return list(filter(lambda hsp: hsp.is_minus(), self.hsps))

    @property
    def num_minus_hsps(self):
        return len(self.minus_hsps)

    @property
    def min_queryp(self):
        query_froms = [hsp.query_from for hsp in self.hsps]
        query_tos = [hsp.query_to for hsp in self.hsps]
        return min([*query_froms, *query_tos])


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

    @property
    def hit_from(self):
        return self.data['hit_from']

    @property
    def hit_to(self):
        return self.data['hit_to']

    @property
    def covered_pos(self):
        if self.is_plus():
            assert self.hit_from < self.hit_to
            return [p for p in range(self.hit_from, self.hit_to + 1)]
        else:
            assert self.hit_to < self.hit_from
            return [p for p in range(self.hit_to, self.hit_from + 1)]


#blast_output_file_path = 'blast_output_cy1_nb_2wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_2wpi_root.json'
#blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_6wpi_root.json'

blast_output_file_path = 'blast_output_cy2_nb_14wpi_leaf.json'

#blast_output_file_path = 'blast_output_cy2_hemp_leaf.json'

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


searches_that_hit = list(filter(lambda search : search.num_hits > 0, all_searches))
print(f'{len(searches_that_hit)=}')
print()


viral_searches = list(filter(lambda search : search.is_viral(), searches_that_hit))
print(f'{len(viral_searches)=}')
print()


plus_searches = list(filter(lambda search : search.is_plus(), viral_searches))
print(f'{len(plus_searches)=}')
print()

type_I_minus_searches = list(filter(lambda search : search.is_type_I_minus(), viral_searches))
print(f'{len(type_I_minus_searches)=}')
print()

foldback_searches = list(filter(lambda search : search.is_foldback(), viral_searches))
print(f'{len(foldback_searches)=}')
print()


def unique_covered_pos(hsps):
    return set([p for hsp in hsps for p in hsp.covered_pos])


fig, ax = plt.subplots()

plt.hist(
    #[p for search in plus_searches for p in unique_covered_pos(search.hit.hsps)],
    #[p for search in type_I_minus_searches for p in unique_covered_pos(search.hit.hsps)],
    [p for search in foldback_searches for p in unique_covered_pos(search.hit.minus_hsps)],
    color='black',
    #bins=[i + 0.5 for i in range(0, 2692 + 1, 1)],
    bins=[i + 0.5 for i in range(0, 2983 + 1, 1)],
)

#ax.set_xticks([1, 3, 4, 8, 10, 14, 15, 20, 23, 25, 30, 35])

ax.set_xticks([2068, 2068 + 3, 2068 + 9, 2068 + 34])

#ax.set_xticks([2692 - 34, 2692 - 3, 2692])
#ax.set_xticks([2983 - 34, 2983 - 3, 2983])

#plt.gca().invert_xaxis()

plt.show()
