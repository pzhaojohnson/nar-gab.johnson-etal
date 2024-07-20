print()


reads_file_path = 'all_reads_ivt_cy1_gRNA.fastq'

#reads_file_path = 'all_reads_cy1_nb_2wpi_leaf.fastq'
#reads_file_path = 'all_reads_cy1_nb_2wpi_root.fastq'
#reads_file_path = 'all_reads_cy1_nb_6wpi_leaf.fastq'
#reads_file_path = 'all_reads_cy1_nb_6wpi_root.fastq'

print('Reads file path: ' + reads_file_path)
print()


with open(reads_file_path, 'r') as f:
    lines = f.read().splitlines()
print('Successfully read in reads file.')
print()


reads = {}
for i in range(0, len(lines), 4):
    assert lines[i][0] == '@'
    assert lines[i][37] == ' '
    read_id = lines[i][1:37]
    assert len(read_id) == 36
    read_seq = lines[i + 1]
    assert len(read_seq) > 0
    reads[read_id] = read_seq
print('Reads found: ' + str(len(reads)))
print()
