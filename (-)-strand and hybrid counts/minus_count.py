import json

import functools


print()


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
    def CY_hit(self):
        CY_hits = list(filter(lambda hit : hit.isto_CY1() or hit.isto_CY2(), self.hits))
        assert len(CY_hits) == 1
        return CY_hits[0]

    def hasa_CY_hit(self):
        try:
            self.CY_hit
            return True
        except:
            return False

    def is_plus(self):
        return all([hsp.is_plus() for hsp in self.CY_hit.hsps])

    def is_minus(self):
        return all([hsp.is_minus() for hsp in self.CY_hit.hsps]) \
            and self.CY_hit.min_queryp / self.read_length <= 0.05

    def is_foldback(self):
        if all([hsp.is_minus() for hsp in self.CY_hit.hsps]):
            return 0.4 <= self.CY_hit.min_queryp / self.read_length <= 0.53
        else:
            return self.CY_hit.num_hsps == 2 and self.CY_hit.num_plus_hsps == 1 and self.CY_hit.num_minus_hsps == 1

    def is_plus_minus_hybrid(self):
        if self.is_foldback():
            return True
        else:
            return self.CY_hit.num_plus_hsps > 0 and self.CY_hit.num_minus_hsps > 0


class Hit:
    def __init__(self, data):
        self.data = data

    def isto_CY1(self):
        return self.data['description'][0]['title'].upper() == 'CY1'

    def isto_CY2(self):
        return self.data['description'][0]['title'].upper() == 'CY2'

    @property
    def hsps(self):
        return [Hsp(hsp) for hsp in self.data['hsps']]

    @property
    def hsps_sortedby_query_from(self):
        hsps = self.hsps
        hsps.sort(key=functools.cmp_to_key(cmp_query_froms))
        return hsps

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
        return list(filter(lambda hsp : hsp.is_minus(), self.hsps))

    @property
    def num_minus_hsps(self):
        return len(self.minus_hsps)

    @property
    def min_queryp(self):
        assert all([hsp.query_from < hsp.query_to for hsp in self.hsps])
        return self.hsps_sortedby_query_from[0].query_from


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


#blast_output_file_path = 'blast_output_ivt_cy1_gRNA.json'

#blast_output_file_path = 'blast_output_cy1_nb_2wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_2wpi_root.json'
#blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_6wpi_root.json'

#blast_output_file_path = 'blast_output_cy2_nb_14wpi_leaf.json'
blast_output_file_path = 'blast_output_cy2_hemp_leaf.json'

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


CY_searches = list(filter(lambda search : search.hasa_CY_hit(), all_searches))
print(f'{len(CY_searches)=}')
print()


plus_searches = list(filter(lambda search : search.is_plus(), CY_searches))
print(f'{len(plus_searches)=}')
print(f'{100 * len(plus_searches) / len(CY_searches)=}')
print()

minus_searches = list(filter(lambda search : search.is_minus(), CY_searches))
print(f'{len(minus_searches)=}')
print(f'{100 * len(minus_searches) / len(CY_searches)=}')
print()

foldback_searches = list(filter(lambda search : search.is_foldback(), CY_searches))
print(f'{len(foldback_searches)=}')
print(f'{100 * len(foldback_searches) / len(CY_searches)=}')
print()

plus_minus_hybrid_searches = list(filter(lambda search : search.is_plus_minus_hybrid(), CY_searches))
print(f'{len(plus_minus_hybrid_searches)=}')
print(f'{100 * len(plus_minus_hybrid_searches) / len(CY_searches)=}')
print()
