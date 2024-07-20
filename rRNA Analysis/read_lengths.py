import json

import matplotlib.pyplot as plt


print()


reads_file_path = 'all_reads_cy1_nb_6wpi_leaf.fastq'

print('Reads file path: ' + reads_file_path)
print()


reads = {}

with open(reads_file_path, 'r') as f:
    lines = f.read().splitlines()
    for i in range(0, len(lines), 4):
        assert lines[i][0] == '@'
        read_id = lines[i][1:37]
        assert len(read_id) == 36
        read_seq = lines[i + 1]
        assert len(read_seq) > 0
        for c in read_seq:
            assert c in 'AUGC'
        reads[read_id] = read_seq

print('Reads found: ' + str(len(reads)))
print()


fig, ax = plt.subplots()

plt.hist(
    [len(read_seq) for read_id, read_seq in reads.items()],
    bins=200,
    color='black',
)

ax.set_xticks([0, 150, 225, 500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000])

plt.xticks(rotation=90)

plt.show()
