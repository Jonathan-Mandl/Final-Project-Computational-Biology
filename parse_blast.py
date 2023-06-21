from Bio.SearchIO._legacy import NCBIStandalone

def parse_blast_file():
    # Open the BLAST output file
    result_handle = open("nr_blast_output.txt", 'r')
    blast_parser = NCBIStandalone.BlastParser() # object for parsing blast file output
    blast_iterator = NCBIStandalone.Iterator(result_handle, blast_parser) # iterator for blast parser output

    # iterate over blast records
    for blast_record in blast_iterator:
        # extract query id
        query_id = blast_record.query.split()[0]
        # iterate over alignments in for blast record. only take first alignment
        for alignment in blast_record.alignments:
            hit_id = alignment.title.split()[0][1:]  # Extract the hit ID from the alignment title
            first_hsp = alignment.hsps[0]  # stores alignment object
            # call align() method to find different indexes using the alignment
            diff_indexes = align(first_hsp.query_start, first_hsp.query,first_hsp.sbjct)
            # display alignment statistics
            print("Query: human protein ", query_id)
            print("Hit: pan troglodytes ", hit_id)
            print("Alignment length:", first_hsp.align_length)
            print("E-value:", first_hsp.expect)
            print("Score:", first_hsp.score)
            print("Sequence identity:", first_hsp.identities)
            print("Query start:", first_hsp.query_start)
            print("Alignment:")
            print(first_hsp.query)
            print(first_hsp.sbjct)
            print("-" * 50)

            # write different indexes to a file
            with open('diff_indexes_updated.txt', 'a') as f:
                row = f"{query_id}\t{hit_id}\t{','.join(map(str, diff_indexes))}\n"
                f.write(row)
            break

def align(human_start,human_seq, chimp_seq):

    diff_indexes = [] # list of different indexes between the 2 proteins.
    human_index = human_start-1
    # iterate over the aligned residues in each row
    for j in range(len(human_seq)):
        res1 = human_seq[j]
        res2 = chimp_seq[j]
        # if there is no gap in human protein, increase human index by 1
        if res1 != "-":
            human_index += 1
        # check for a difference between the residues
        if res1 != res2 and res1 != "-" and res2 != "-":
            # if there is a difference between residues, add human index to diff_indexes array
            diff_indexes.append(human_index)
    return diff_indexes