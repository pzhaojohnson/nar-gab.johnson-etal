import json

import functools

import matplotlib.pyplot as plt


print()


class Search:
    def __init__(self, data):
        self.data = data

    @property
    def read_length(self):
        return self.data['query_len']

    @property
    def num_hits(self):
        return len(self.data['hits'])

    def has_a_hit(self):
        return self.num_hits > 0

    @property
    def hit(self):
        hits = self.data['hits']
        assert len(hits) == 1
        return Hit(hits[0])

    def is_minus(self):
        return self.hit.is_minus()

    def is_type_i_minus(self):
        return self.is_minus() \
            and self.hit.hsps_sorted_by_query_from[0].query_from / self.read_length <= 0.05

    @property
    def num_unique_covered_pos(self):
        return len(self.hit.unique_covered_pos)


def cmp_num_unique_covered_pos(search1, search2):
    return search1.num_unique_covered_pos - search2.num_unique_covered_pos


class Hit:
    def __init__(self, data):
        self.data = data

    @property
    def hsps_sorted_by_query_from(self):
        hsps = [Hsp(hsp) for hsp in self.data['hsps']]
        hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
        return hsps

    def is_minus(self):
        return all([hsp.is_minus() for hsp in self.hsps_sorted_by_query_from])

    @property
    def unique_covered_pos(self):
        return set([p for hsp in self.hsps_sorted_by_query_from for p in hsp.covered_pos])


class Hsp:
    def __init__(self, data):
        self.data = data

    @property
    def query_from(self):
        return self.data['query_from']

    @property
    def hit_from(self):
        return self.data['hit_from']

    @property
    def hit_to(self):
        return self.data['hit_to']

    def is_minus(self):
        return self.data['hit_strand'] == 'Minus'

    @property
    def covered_pos(self):
        # the below range only works correctly if the hsp is minus!
        assert self.is_minus()
        return [p for p in range(self.hit_to, self.hit_from + 1, 1)]


def cmp_query_froms(hsp1, hsp2):
    return hsp1.query_from - hsp2.query_from


def are_within(num1, num2, maxDiff):
    return abs(num1 - num2) <= maxDiff


class HspRange:
    def __init__(self, hit_from, hit_to):
        self.hit_from = hit_from
        self.hit_to = hit_to

    def closely_matches(self, hsp):
        return are_within(self.hit_from, hsp.hit_from, 30) and are_within(self.hit_to, hsp.hit_to, 30)


def matches_full_length_gRNA(search):
    assert search.is_minus()
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(2983, 1).closely_matches(hsps[0])


def matches_F281(search):
    assert search.is_minus()
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(281, 1).closely_matches(hsps[0])


def matches_F442(search):
    assert search.is_minus()
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(442, 1).closely_matches(hsps[0])


def matches_F671(search):
    assert search.is_minus()
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(671, 1).closely_matches(hsps[0])


def matches_DRNA(search):
    assert search.is_minus()
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 2 \
        and are_within(2983, hsps[0].hit_from, 30) \
        and are_within(1, hsps[1].hit_to, 30) \
        and search.num_unique_covered_pos < 2000


def matches_sgRNA(search):
    assert search.is_minus()
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(2983, 2068).closely_matches(hsps[0])


def plotted_color(search):
    if matches_full_length_gRNA(search):
        return 'red'
    elif matches_F281(search):
        return 'dodgerblue'
    elif matches_F442(search):
        return 'magenta'
    elif matches_F671(search):
        return 'darkorange'
    elif matches_DRNA(search):
        return 'lime'
    elif matches_sgRNA(search):
        return 'mediumturquoise'
    else:
        return 'black'


blast_output_file_path = 'blast_output_cy2_nb_14wpi_leaf.json'

print(f'{blast_output_file_path=}')
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output!')
print()


reports = [item['report'] for item in blast_output['BlastOutput2']]
print(f'{len(reports)=}')
print()


all_searches = [Search(report['results']['search']) for report in reports]
print(f'{len(all_searches)=}')
print()


searches_with_a_hit = list(filter(lambda search : search.has_a_hit(), all_searches))
print(f'{len(searches_with_a_hit)=}')
print()


for search in searches_with_a_hit:
    assert search.num_hits == 1
print('All searches have at most one hit.')
print()


minus_searches = list(filter(lambda search : search.is_minus(), searches_with_a_hit))
print(f'{len(minus_searches)=}')
print()


type_i_minus_searches = list(filter(lambda search : search.is_type_i_minus(), minus_searches))
print(f'{len(type_i_minus_searches)=}')
print()


print(f'{100 * len(list(filter(matches_full_length_gRNA, type_i_minus_searches))) / len(type_i_minus_searches)=}')
print()

print(f'{100 * len(list(filter(matches_sgRNA, type_i_minus_searches))) / len(type_i_minus_searches)=}')
print()

print(f'{100 * len(list(filter(matches_DRNA, type_i_minus_searches))) / len(type_i_minus_searches)=}')
print()

print(f'{100 * len(list(filter(matches_F671, type_i_minus_searches))) / len(type_i_minus_searches)=}')
print()

print(f'{100 * len(list(filter(matches_F442, type_i_minus_searches))) / len(type_i_minus_searches)=}')
print()

print(f'{100 * len(list(filter(matches_F281, type_i_minus_searches))) / len(type_i_minus_searches)=}')
print()


"""
fig, ax = plt.subplots()

type_i_minus_searches.sort(key=functools.cmp_to_key(cmp_num_unique_covered_pos), reverse=True)

y = 1

for search in type_i_minus_searches:
    color = plotted_color(search)
    alpha = 1 if color == 'black' else 1
    hsps = search.hit.hsps_sorted_by_query_from
    for hsp in hsps:
        plt.plot([hsp.hit_from, hsp.hit_to], [y, y], color=color, alpha=alpha)
    y += 1

ax.set_xticks([1, 281, 442, 671, 2068, 2710, 2983])

plt.show()
#"""


#"""
fig, ax = plt.subplots()

plt.hist(
    [p for search in type_i_minus_searches for p in search.hit.unique_covered_pos],
    color='black',
    bins=[i + 0.5 for i in range(0, 2983 + 1, 1)],
)

ax.set_xticks([1, 281, 442, 671, 2068, 2710, 2983])

plt.show()
#"""
