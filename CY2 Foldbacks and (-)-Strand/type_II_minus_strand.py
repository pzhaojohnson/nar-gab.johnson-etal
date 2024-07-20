import json

import matplotlib.pyplot as plt


print()


blast_output_file_path = 'blast_output_cy2_nb_14wpi_leaf.json'
print(f'{blast_output_file_path=}')
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output!')
print()


reports = [item['report'] for item in blast_output['BlastOutput2']]

all_searches = [report['results']['search'] for report in reports]
print(f'{len(all_searches)=}')
print()


def isfor_CY2(search):
    hits = search['hits']
    assert len(hits) <= 1
    return len(hits) == 1 and hits[0]['description'][0]['title'].upper() == 'CY2'


CY2_searches = list(filter(isfor_CY2, all_searches))
print(f'{len(CY2_searches)=}')
print()


def is_plus(search):
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    return all([hsp['hit_strand'] == 'Plus' for hsp in hit['hsps']])


def is_minus(search):
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    return all([hsp['hit_strand'] == 'Minus' for hsp in hit['hsps']])


plus_searches = list(filter(is_plus, CY2_searches))
print(f'{len(plus_searches)=}')
print()

minus_searches = list(filter(is_minus, CY2_searches))
print(f'{len(minus_searches)=}')
print()


def normalized_upstream_length(minus_search):
    hits = minus_search['hits']
    assert len(hits) == 1
    hit = hits[0]
    assert all([hsp['hit_strand'] == 'Minus' for hsp in hit['hsps']])
    return min([hsp['query_from'] for hsp in hit['hsps']]) / minus_search['query_len']


type_II_minus_searches = list(filter(
    lambda search : 0.4 <= normalized_upstream_length(search) <= 0.53,
    minus_searches,
))
print(f'{len(type_II_minus_searches)=}')
print()


"""
fig, ax = plt.subplots()

plt.hist(
    [normalized_upstream_length(search) for search in minus_searches],
    color='black',
    bins=100,
)

ax.set_xticks([-0.15, 0, 0.05, 0.4, 0.53, 1, 1.15])

plt.show()
#"""



#"""
fig, ax = plt.subplots()

y = 1

for search in type_II_minus_searches:
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    hsps = hit['hsps']
    assert all([hsp['hit_strand'] == 'Minus' for hsp in hsps])
    min_queryp = min([hsp['query_from'] for hsp in hsps])
    max_queryp = max([hsp['query_to'] for hsp in hsps])
    plt.plot([1 - min_queryp, max_queryp - min_queryp], [y, y], color='gray')
    for hsp in hsps:
        plt.plot([hsp['query_from'] - min_queryp, hsp['query_to'] - min_queryp], [y, y], color='red')
    y += 1

ax.set_xticks([-2983 - 300, -2983, 0, 2982, 2982 + 300])

plt.show()
#"""
