import json


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
    return len(hits) == 1 \
        and hits[0]['description'][0]['title'].upper() == 'CY2'


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


def is_plus_minus_hybrid(search):
    hits = search['hits']
    assert len(hits) == 1
    return not is_plus(search) and not is_minus(search)


plus_searches = list(filter(is_plus, CY2_searches))
print(f'{len(plus_searches)=}')
print()

minus_searches = list(filter(is_minus, CY2_searches))
print(f'{len(minus_searches)=}')
print()

plus_minus_hybrid_searches = list(filter(is_plus_minus_hybrid, CY2_searches))
print(f'{len(plus_minus_hybrid_searches)=}')
print()


segment_counts = {}

for search in plus_minus_hybrid_searches:
    hits = search['hits']
    assert len(hits) == 1
    hit = hits[0]
    num_segments = str(len(hit['hsps']))
    prev_count = 0 if num_segments not in segment_counts else segment_counts[num_segments]
    segment_counts[num_segments] = prev_count + 1

print(f'{segment_counts=}')
print()
