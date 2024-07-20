import json

import functools

import matplotlib.pyplot as plt


print()


class Search:
    def __init__(self, data):
        self.data = data

    @property
    def num_hits(self):
        return len(self.data['hits'])

    def has_a_hit(self):
        return self.num_hits > 0

    @property
    def hit(self):
        """Returns the single hit for the search.

        Raises if the search does not have exactly one hit.
        """
        hits = [Hit(hit) for hit in self.data['hits']]
        assert len(hits) == 1
        return hits[0]

    def is_plus(self):
        """Returns True if and only if the single hit for the search is plus."""
        return self.hit.is_plus()

    def is_minus(self):
        """Returns True if and only if the single hit for the search is minus."""
        return self.hit.is_minus()

    @property
    def num_unique_covered_pos(self):
        """The number of unique hit positions covered by the hit for the search."""
        return len(self.hit.unique_covered_pos)


def cmp_num_unique_covered_pos(search1, search2):
    return search1.num_unique_covered_pos - search2.num_unique_covered_pos


class Hit:
    def __init__(self, data):
        self.data = data

    def is_to_CY1(self):
        return self.data['description'][0]['title'].upper() == 'CY1'

    def is_to_CY2(self):
        return self.data['description'][0]['title'].upper() == 'CY2'

    def is_to_rubisco_large(self):
        return 'JF419563.1' in self.data['description'][0]['title']

    @property
    def hsps_sorted_by_query_from(self):
        hsps = [Hsp(hsp) for hsp in self.data['hsps']]
        hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
        return hsps

    def is_plus(self):
        """Returns True if and only if all hsps of the hit are plus."""
        return all([hsp.is_plus() for hsp in self.hsps_sorted_by_query_from])

    def is_minus(self):
        """Returns True if and only if all hsps of the hit are minus."""
        return all([hsp.is_minus() for hsp in self.hsps_sorted_by_query_from])

    @property
    def unique_covered_pos(self):
        """A set of the unique hit positions covered by the hsps of the hit."""
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

    def is_plus(self):
        return self.data['hit_strand'] == 'Plus'

    def is_minus(self):
        return self.data['hit_strand'] == 'Minus'

    @property
    def covered_pos(self):
        """A list of the hit positions covered by the hsp."""
        # hit-from must be less than or equal to hit-to for the below range to work!
        assert self.is_plus()
        return [p for p in range(self.hit_from, self.hit_to + 1, 1)]


def cmp_query_froms(hsp1, hsp2):
    return hsp1.query_from - hsp2.query_from


def are_within(num1, num2, maxDiff):
    """Returns True if and only if the absolute difference between the two number is less than or equal to the max difference."""
    return abs(num1 - num2) <= maxDiff


class HspRange:
    def __init__(self, end5, end3):
        """The range of an hsp defined by 5' and 3' terminal end positions."""
        self.end5 = end5
        self.end3 = end3

    def matches(self, hsp):
        """Returns True if and only if:
            - The 5'end position of the hsp is within 30 nt of the 5'end position of this hsp range.
            - And the 3'end position of the hsp is within 10 nt of the 3'end position of this hsp range.

        (5'ends are given a bit more leeway due to the relative imprecision of 5'end sequencing with DRS.)
        """
        assert hsp.is_plus()
        return are_within(hsp.hit_from, self.end5, 30) and are_within(hsp.hit_to, self.end3, 10)


def matches_full_length_gRNA(search):
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(1, 2692).matches(hsps[0])


def matches_F281(search):
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(1, 281).matches(hsps[0])


def matches_F442(search):
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(1, 442).matches(hsps[0])


def matches_F671(search):
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(1, 671).matches(hsps[0])


def matches_DRNA(search):
    assert search.is_plus()
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 2 \
        and are_within(1, hsps[0].hit_from, 30) \
        and are_within(2692, hsps[1].hit_to, 10) \
        and search.num_unique_covered_pos < 2000


def matches_full_length_rubisco_large(search):
    hsps = search.hit.hsps_sorted_by_query_from
    return len(hsps) == 1 \
        and HspRange(1, 1434).matches(hsps[0])


def plotted_color(search):
    #return 'red' if matches_full_length_rubisco_large(search) else 'black'
    #"""
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
    else:
        return 'black'
    #"""


#blast_output_file_path = 'blast_output_cy1_nb_2wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_2wpi_root.json'
blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_6wpi_root.json'

#blast_output_file_path = 'blast_to_rubisco_large_output_cy1_nb_2wpi_leaf.json'
#blast_output_file_path = 'blast_to_rubisco_large_output_cy1_nb_2wpi_root.json'
#blast_output_file_path = 'blast_to_rubisco_large_output_cy1_nb_6wpi_leaf.json'
#blast_output_file_path = 'blast_to_rubisco_large_output_cy1_nb_6wpi_root.json'

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
    assert search.hit.is_to_CY1()
    #assert search.hit.is_to_rubisco_large()
print('All hits are to CY1.')
#print('All hits are to rubisco large.')
print()


plus_searches = list(filter(lambda search : search.is_plus(), searches_with_a_hit))
print(f'{len(plus_searches)=}')
print()

minus_searches = list(filter(lambda search : search.is_minus(), searches_with_a_hit))
print(f'{len(minus_searches)=}')
print()


#print(f'{100 * len(list(filter(matches_full_length_rubisco_large, plus_searches))) / len(plus_searches)=}')
#print()

#"""
print(f'{100 * len(list(filter(matches_full_length_gRNA, plus_searches))) / len(plus_searches)=}')
print()

print(f'{100 * len(list(filter(matches_DRNA, plus_searches))) / len(plus_searches)=}')
print()

print(f'{100 * len(list(filter(matches_F671, plus_searches))) / len(plus_searches)=}')
print()

print(f'{100 * len(list(filter(matches_F442, plus_searches))) / len(plus_searches)=}')
print()

print(f'{100 * len(list(filter(matches_F281, plus_searches))) / len(plus_searches)=}')
print()
#"""


#"""
fig, ax = plt.subplots()

plus_searches.sort(key=functools.cmp_to_key(cmp_num_unique_covered_pos), reverse=True)

y = 1

for search in plus_searches:
    color = 'black'
    #color = plotted_color(search)
    alpha = 0.05 if color == 'black' else 1
    hsps = search.hit.hsps_sorted_by_query_from
    for hsp in hsps:
        plt.plot([hsp.hit_from, hsp.hit_to], [y, y], color=color, alpha=alpha)
    y += 1

ax.set_xticks([1, 281, 442, 671, 2331, 2385, 2420, 2580, 2692])
#ax.set_xticks([1, 500, 1000, 1434])

plt.show()
#"""
