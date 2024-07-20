import json


print()


class BLASTOutput:
    def __init__(self, data):
        self.data = data

    @property
    def reports(self):
        """All reports in the BLAST output."""
        return [item['report'] for item in self.data['BlastOutput2']]

    @property
    def searches(self):
        """All searches in the BLAST output."""
        searches = [Search(report['results']['search']) for report in self.reports]
        assert len(searches) > 0
        return searches


class Search:
    def __init__(self, data):
        self.data = data

    @property
    def hits(self):
        """All hits for the search."""
        return [Hit(hit) for hit in self.data['hits']]

    @property
    def hit(self):
        """The single hit for the search.

        Raises if the search does not have exactly one hit.
        """
        assert len(self.hits) == 1
        return self.hits[0]


class Hit:
    def __init__(self, data):
        self.data = data

    @property
    def hsps(self):
        """The hsps of the hit."""
        return [Hsp(hsp) for hsp in self.data['hsps']]

    @property
    def num_hsps(self):
        """The number of hsps for the hit."""
        return len(self.hsps)


class Hsp:
    def __init__(self, data):
        self.data = data

    @property
    def query_from(self):
        query_from = self.data['query_from']
        assert type(query_from) == int
        return query_from


blast_output_file_path = 'blast_output_cy1_nb_6wpi_leaf.json'

with open(blast_output_file_path, 'r') as f:
    blast_output = json.loads(f.read())

reports = [item['report'] for item in blast_output['BlastOutput2']]

searches = [report['results']['search'] for report in reports]

for search in searches:
    assert len(search['hits']) <= 1

def is_for_a_cy1_read(search):
    hits = search['hits']
    return len(hits) == 1 and hits[0]['description'][0]['title'].upper() == 'CY1'

searches_for_cy1_reads = list(filter(is_for_a_cy1_read, searches))
print(f'{len(searches_for_cy1_reads)=}')
print()

def is_for_a_plus_strand_read(search):
    assert is_for_a_cy1_read(search)
    hsps = search['hits'][0]['hsps']
    return all([hsp['hit_strand'] == 'Plus' for hsp in hsps])

searches_for_plus_strand_reads = list(filter(is_for_a_plus_strand_read, searches_for_cy1_reads))
print(f'{len(searches_for_plus_strand_reads)=}')
print()


def is_for_an_F281_read(search):
    assert is_for_a_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    return len(hsps) == 1 and hsps[0]['hit_from'] <= 31 and 281 - 10 <= hsps[0]['hit_to'] <= 281 + 10

searches_for_F281_reads = list(filter(is_for_an_F281_read, searches_for_plus_strand_reads))
print(f'{len(searches_for_F281_reads)=}')
print()


def is_for_an_F671_read(search):
    assert is_for_a_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    return len(hsps) == 1 and hsps[0]['hit_from'] <= 31 and 671 - 10 <= hsps[0]['hit_to'] <= 671 + 10

searches_for_F671_reads = list(filter(is_for_an_F671_read, searches_for_plus_strand_reads))
print(f'{len(searches_for_F671_reads)=}')
print()



def is_for_a_gRNA_read(search):
    assert is_for_a_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    return len(hsps) == 1 and hsps[0]['hit_from'] <= 31 and 2692 - 10 <= hsps[0]['hit_to']

searches_for_gRNA_reads = list(filter(is_for_a_gRNA_read, searches_for_plus_strand_reads))
print(f'{len(searches_for_gRNA_reads)=}')
print()



def is_for_a_DRNA_read(search):
    assert is_for_a_plus_strand_read(search)
    hsps = search['hits'][0]['hsps']
    return len(hsps) == 2 and hsps[0]['hit_from'] <= 31 and 2692 - 10 <= hsps[1]['hit_to']

searches_for_DRNA_reads = list(filter(is_for_a_DRNA_read, searches_for_plus_strand_reads))
print(f'{len(searches_for_DRNA_reads)=}')
print()
