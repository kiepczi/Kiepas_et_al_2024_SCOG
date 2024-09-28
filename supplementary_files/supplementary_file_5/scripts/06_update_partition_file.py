
"""Our concatanated alignment was trimmed using trimAL with an -automated1 parameter, and columns that were fully conserved 
where also removed.
This script was used to update the partition file for modeltest-ng. 

The relationship between the columns in the old and new alignment can be reported with the -colnumbering paratemer
This can be used to update the partition file, which will be used as an imput for modeltest-ng.

Thing to remember: trimAl indexing starts at 0, where modeltestng at 1!"""



def update_part_file(colmatchfile, partfile):
    """Update the modeltest-ng partition file.
    """

    #Loading the trimal -colnumbering output 
    alg_col_match = open(colmatchfile, 'r')
    alg_col_match = [int(_)+1 for _ in alg_col_match.read().replace("#ColumnsMap\t", '').replace(' ', '').replace('\n', '').split(',')] #Update the column numbers by adding 1 (trimal indexing starts at 0, and modeltest-ng at 1)

    #Loading the parition file
    part_file = open(partfile, 'r')
    part_file = part_file.readlines()


    partition = []
    last_position = 0
    #Getting the length of the alignment from the partition file
    part_file_len = max([int(_.replace('\n', '').split('-')[-1]) for _ in part_file])
    #Getting the list of columns that were removed with trimAl
    removed_col = [int(_) for _ in range(1, part_file_len+1) if int(_) not in alg_col_match]
    
 
    for _ in part_file:
        start = int(_.split('= ')[-1].split('-')[0]) #start of the current partition
        end = int(_.split('= ')[-1].split('-')[1]) #end of the current partition
        OG = _.split(', ')[1].split(' =')[0] # the partition name 
        col_match = len([_ for _ in range(start, end+1) if _ not in removed_col]) #here we can get the total number of positions that remained in the new alignment
        updated_pos = last_position + col_match #updating the position
        partition.append(f"DNA, {OG} = {last_position+1}-{updated_pos}")
        last_position = updated_pos


    return partition


partition_file_no_gaps = update_part_file('../output/alignments/concatenated/no_gaps_concatenated_columns.txt', '../output/alignments/concatenated/initial_concatenated_modeltest.part')

#we can save the partition
with open("../output/alignments/concatenated/concatenated_modeltest_fixed_positions.part", "w") as text_file:
    text_file.write('\n'.join(partition_file_no_gaps))

