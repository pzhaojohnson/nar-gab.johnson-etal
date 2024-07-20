import json

import matplotlib.pyplot as plt

import random


print()


class Search:
    def __init__(self, data):
        self.data = data

    @property
    def hits(self):
        return [Hit(hit) for hit in self.data['hits']]

    @property
    def CY1_hit(self):
        CY1_hits = list(filter(lambda hit : hit.isto_CY1(), self.hits))
        assert len(CY1_hits) == 1
        return CY1_hits[0]

    def hasa_CY1_hit(self):
        try:
            self.CY1_hit
            return True
        except:
            return False


class Hit:
    def __init__(self, data):
        self.data = data

    def isto_CY1(self):
        return self.data['description'][0]['title'].upper() == 'CY1'

    @property
    def hsps(self):
        return [Hsp(hsp) for hsp in self.data['hsps']]

    def is_plus(self):
        return all([hsp.is_plus() for hsp in self.hsps])

    @property
    def unique_covered_pos(self):
        return set([p for hsp in self.hsps for p in hsp.covered_pos])


class Hsp:
    def __init__(self, data):
        self.data = data

    @property
    def hit_from(self):
        return self.data['hit_from']

    @property
    def hit_to(self):
        return self.data['hit_to']

    def is_plus(self):
        return self.data['hit_strand'] == 'Plus'

    @property
    def covered_pos(self):
        assert self.is_plus()
        return [p for p in range(self.hit_from, self.hit_to + 1, 1)]


#blast_output_file_path = 'blast_output_ivt_cy1_gRNA.json'
#blast_output_file_path = 'blast_output_cy1_pemv2_pfbv_in_line.json'

#blast_output_file_path = 'blast_output_cy1_nb_2wpi_leaf.json'
#blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'
blast_output_file_path = 'blast_output_cy1_nb_6wpi_root.json'


print('BLAST output file path: ' + blast_output_file_path);
print()


with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())
print('Successfully parsed BLAST output.')
print()


reports = [output['report'] for output in blast_output['BlastOutput2']]
print('Reports: ' + str(len(reports)))
print()


all_searches = [Search(report['results']['search']) for report in reports]
print('All Searches: ' + str(len(all_searches)))
print()


CY1_searches = list(filter(lambda search : search.hasa_CY1_hit(), all_searches))
print(f'{len(CY1_searches)=}')
print()


plus_searches = list(filter(lambda search : search.CY1_hit.is_plus(), CY1_searches))
print(f'{len(plus_searches)=}')
print()


fig, ax = plt.subplots()

plt.hist(
    [p for search in plus_searches for p in search.CY1_hit.unique_covered_pos],
    color='black',
    bins=[i + 0.5 for i in range(0, 2692 + 1, 1)],
)

plt.show()
