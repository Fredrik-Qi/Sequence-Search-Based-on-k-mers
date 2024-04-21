'''
Sequence Mapping Program

This Python script performs sequence mapping using k-mer and seed-and-extend approaches. 
It reads input FASTA files containing reference and read sequences, while builds k-string libraries, 
and generates seeds from first k mers of read sequences and matches to the library of references. 
The program then extends these seeds with a threshold of success matching rate to identify mappings 
and produces three output files, including a histogram plot of mapping counts, a txt file containing 
information on reference coverage and another txt file containing mapping information.

Author: Zhouhui Qi
Date: 2024-02-14
Example:
  python script_name.py --k 10 reference.fasta reads.fasta output/


'''
import argparse
import matplotlib.pyplot as plt
def check_dna(seq):
    """
    Function: Check if a DNA sequence contains valid characters.
    
    Parameters:
    - seq (str): The DNA sequence to be checked.
    
    Error handling:
    - ValueError: If the sequence contains at least one invalid character.
    
    Notes:
    - This function is adapted from "lab material," Lab 2, "Loading a 
        sequence database from multiple files" Part in this course.
    """
    if len(set(seq) - set({'T', 'A', 'G', 'C'})) > 0:
        raise ValueError(f"DNA sequnece {seq} contains "
                         f"at least one invalid character")

def read_fasta_to_dict(file,k):
    """
    Function: Read a FASTA file and return a dictionary of 
    sequences with reads/ref IDs in keys and sequences in values.

    Parameters:
    - file (str): The filename of the FASTA file.
    - k (int): The number of k-mer, which is also the minimum 
        length required for each sequence.

    Returns:
    - dict: A dictionary containing sequence IDs as keys and 
        corresponding sequences as values.

    Error handling:
    - FileNotFoundError: If the specified file does not exist.
    - ValueError: If a sequence is shorter than the specified 
        length (k) or contains invalid characters.

    Notes:
    - The function assumes that the input file is in FASTA format.
    - Sequences are checked for valid DNA characters using the 
        check_dna function.
    - If a sequence is shorter than k, a ValueError is raised.
    """
    sequences = {}
    try:
        with open(file, 'r') as file:
            lines = file.readlines()
            sequence_id = None
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    sequence_id = line[1:]
                    sequences[sequence_id] = ''
                else:
                    check_dna(line)
                    if len(line) < k:
                        raise ValueError(f"Sequence {sequence_id} "
                                        f"is shorter than k ({k})")
                    sequences[sequence_id] += line
    except FileNotFoundError:
        print(f"The file {file} you entered does not exist." 
              f"Please check it again.")
        raise SystemExit 
    return sequences

def calculate_length(sequences_dict):
    """
    Function: Calculate the length of each sequence in a dictionary.

    Parameters:
    - sequences_dict (dict): A dictionary containing sequence IDs 
        as keys and sequences as values.

    Returns:
    - dict: A dictionary containing sequence IDs as keys and their 
        respective lengths as values.

    Notes:
    - This function can facilitate the direct indexing of the corresponding 
        read and ref lengths based on read_id and ref_id when calculating 
        coverage (calculate_reference_coverage())
    """
    sequences_length_dict = {}
    for sequence_id, sequence in sequences_dict.items():
        sequences_length_dict[sequence_id] = len(sequence)
    return sequences_length_dict

def build_library(k, sequences_dict):
    """
    Fuction: Build a k-mer library for a given set of sequences.

    Parameters:
    - k (int): The length of the k-mer.
    - sequences_dict (dict): A dictionary containing sequence IDs 
        as keys and sequences as values.

    Returns:
    - dict: A dictionary where each sequence ID is associated 
        with a dictionary of k-string positions and k-strings.

    Example:
    >>> sequences_dict = {'seq1': 'ATCG', 'seq2': 
                        'GATCGATC', 'seq3': 'CGAT'}
    >>> build_library(3, sequences_dict)
    # Returns {'seq1': {0: 'ATC', 1: 'TCG'}, 
        'seq2': {0: 'GAT', 1: 'ATC', 2: 'TCG', 
        3: 'CGA', 4: 'GAT', 5: 'ATC', 6: 'TCG'}, 
        'seq3': {0: 'CGA', 1: 'GAT'}}
    """
    library = {}
    for sequence_id, sequence in sequences_dict.items():
        kstring_position = {}  
        for i in range(len(sequence) - k + 1):
            substring = sequence[i:i+k]
            kstring_position[i] = substring
        library[sequence_id] = kstring_position
    return library

def compare_strings(read_str, ref_str):
    """
    Function:Compare two strings character by character.

    Parameters:
    - read_str (str): The first string to be compared.
    - ref_str (str): The second string to be compared.

    Returns:
    - tuple: A tuple containing the count of matches, 
        count of mismatches, and matching rate.
    """
    match_count = 0
    mismatch_count = 0

    for i in range(len(read_str)):
        if read_str[i] == ref_str[i]:
            match_count += 1
        else:
            mismatch_count += 1

    success_rate = match_count / len(read_str)
    return match_count, mismatch_count, success_rate

def generate_seed(onereadin_read_library, 
                onerefin_ref_library,extension_threshold):
    """
    Function: Generate seed alignments based on comparing the first
        k-mer of read with reference reference k-strings library.
        If the matching rate of the alignment is greater than
        the threshold, then this first k-mer of read will be 
        treated as a seed and stored in a dictionary.

    Parameters:
    - onereadin_read_library (dict): A dictionary containing k-mers 
        of a read, with k-mers position as keys and k-mers as values.
    - onerefin_ref_library (dict): A dictionary containing k-mers  
        of a reference, with k-mers position as keys and k-mers as values.
    - extension_threshold (float): The minimum matching rate required 
        for a seed alignment.

    Returns:
    - dict: A dictionary containing seed alignments with relevant information.
        The keys are tuples of (read_position, ref_position), and the values
        are dictionaries with matching informations, like mismatch count.

    Example:
    >>> onereadin_read_library = {0: 'ATCG'}
    >>> onerefin_ref_library = {1: 'ATCG', 2: 'GATC'}
    >>> extension_threshold = 0.8
    >>> generate_seed(onereadin_read_library, 
                onerefin_ref_library, extension_threshold)
    # Returns {(0, 1): {'read_value': 'ATCG', 'ref_value': 'ATCG', 
                'match_count': 4, 'mismatch_count': 0}}
    """
    results = {}
    for ref_position, ref_kstring in onerefin_ref_library.items():
        match_count, mismatch_count, success_rate = compare_strings(
                        onereadin_read_library[0],ref_kstring)
        if success_rate > extension_threshold and ref_position:
            key = (0, ref_position)
            results[key] = {
                'read_value': onereadin_read_library[0],
                'ref_value': ref_kstring,
                'match_count': match_count,
                'mismatch_count': mismatch_count,
            }
            
    return results
def end_ref_error():
    """
    Function: Print an error message when the mapping is outside
    """
    print("Error: Mappng outside the border of the reference")

def extend_seeds(result, onereadin_read_library, 
    onerefin_ref_library, k,extension_threshold):
    """
    Function: Extend seed alignments based on comparing the following
        k-mers of read and reference. The main structure of this 
        function is composed of while and else statements. The 
        else code block will only be executed when the expression 
        of while returns False. If the loop is interrupted by the 
        break statement, which means the alignment is not right, 
        the expression won't return False, so the else code block 
        won't be executed. 
        The while statement is mainly responsible for comparing 
        the complete k-mer, with k as the distance of each movement
        of the reading frame; while the else statement performs the 
        comparison of the incomplete k-mer at the end of the read,
        finally achieving alignment of complete read sequences

    Parameters:
    - result (dict): A dictionary containing informations of seed 
        alignments.
    - onereadin_read_library (dict): A dictionary containing k-mers 
        of a read, with k-mers position as keys and k-mers as values.
    - onerefin_ref_library (dict): A dictionary containing k-mers  
        of a reference, with k-mers position as keys and k-mers as 
        values.
    - extension_threshold (float): The minimum matching rate required 
        for a seed alignment.
    - k (int): The length of the k-mer.

    Returns:
    - dict: A dictionary containing extended seed alignments with 
        relevant information, with read_position, ref_position as 
        key and match_count, mismatch_count as values.

    Example:
    extended_results = {{(0, 4): {'match_count': 9, 'mismatch_count': 0}}
        {(0, 18): {'match_count': 9, 'mismatch_count': 0}}}
    """
    extended_results = {}
    # Calculate the total number of matches, which isthe length of read
    sum_count = len(onereadin_read_library)+k-1
    for key, seed_info in result.items():
        read_position, ref_position = key
        original_read_postion, original_ref_position = read_position, ref_position
        match_count, mismatch_count = seed_info[
            'match_count'], seed_info['mismatch_count']
        sum_match = 0
                
        # Extend the match based on seeds with complete k-mer
        while (
            read_position + k in onereadin_read_library
            and ref_position + k in onerefin_ref_library
            ):
            match_count_extension, mismatch_count_extension, success_rate_extension = compare_strings(
                onereadin_read_library[read_position + k],
                onerefin_ref_library[ref_position + k]
            )
            # If the matching rate is less than the threshold, exit the loop
            if success_rate_extension < extension_threshold:
                # After break, skip the following else and go directly to the next for loop
                break 
            
            # If the matching rate is greater than the threshold, accumulate results
            match_count += match_count_extension
            mismatch_count += mismatch_count_extension
            sum_match = match_count + mismatch_count
            
            read_position += k
            ref_position += k
        
        # If the previous while loop is executed successfully, it will 
        # finally return False, so that all complete k-mers are matched,
        # then the remaining sequence will be matched in the else statement
        else:
            # If read can be divided into an integer number of complete 
            #k-mers, store the final comparison result
            if sum_match == sum_count:
                extended_results[key] = {
                    'match_count': match_count,
                    'mismatch_count': mismatch_count,
                }
            
            # If (read_position + k) is not within the scope of onereadin_read_library
            # or (ref_position + k is) not within the scope of onerefin_ref_library
            elif (
                read_position + k not in onereadin_read_library
                or ref_position + k not in onerefin_ref_library
                ):
                # Report error if the seed matches the the end of the ref 
                # andthe length of the ref is not enough to extend the match.
                if original_ref_position > (
                    len(onerefin_ref_library)-len(onereadin_read_library)):  
                    end_ref_error()
                    break
                # Compare the remaining sequences
                remaining_sequence = onereadin_read_library[
                    original_read_postion+sum_count-k][-(sum_count % k):]
                remaining_ref_sequence = onerefin_ref_library[
                    original_ref_position+sum_count-k][-(sum_count % k):]
                match_count_extension, mismatch_count_extension, success_rate_extension = compare_strings(
                    remaining_sequence,remaining_ref_sequence)
                
                if success_rate_extension < extension_threshold:
                    break

                match_count += match_count_extension
                mismatch_count += mismatch_count_extension
                sum_string = match_count + mismatch_count

                if sum_string == sum_count:
                    extended_results[key] = {
                        'match_count': match_count,
                        'mismatch_count': mismatch_count,
                    }
    return extended_results

def generate_seperated_dict_and_map(read_library,ref_library,k,
            extension_threshold,read_sequnce_length,ref_sequnce_length):
    """
    Function: Traverse each read and reference and use generate_seed() 
        and extend_seed() to get the full mapping results of reads and 
        references.

    Parameters:
    - read_library (dict): dict: A dictionary of read containing k-string 
        positions and k-strings.
    - ref_library (dict): A dictionary of reference containing k-string 
        positions and k-strings.
    - k (int): The length of the k-mer.
    - extension_threshold (float): The matching rate required for extending
        a seed alignment.
    - read_sequence_length (dict): A dictionary containing read IDs as keys 
        and corresponding sequence lengths as values.
    - ref_sequence_length (dict): A dictionary containing reference IDs as 
        keys and corresponding sequence lengths as values.

    Returns:
    - dict: A dictionary containing the tuples (read_id,ref_id) as keys and 
        extended_result as values to store the full mapping results.
    """
    full_mapping_results = {}
    for key1 in read_library.keys():
        onereadin_read_library = read_library[key1]
        for key2 in ref_library.keys():
            if read_sequnce_length[key1] > ref_sequnce_length[key2]:
                raise ValueError(f'The length of read sequnce {key1} is '
                 f'longer than the length of reference sequnce {key2}.')

            onerefin_ref_library = ref_library[key2]
            result = generate_seed(onereadin_read_library, 
                onerefin_ref_library,extension_threshold)
            extended_result = extend_seeds(result, onereadin_read_library, 
                onerefin_ref_library, k,extension_threshold)
            # If extended_result is not empty, store it in full_mapping_results
            if extended_result: 
                full_mapping_results[(key1,key2)] = extended_result
    return full_mapping_results

def draw_bar_plot(dict_ref_info, myoutfile):
    """
    Function: Draw a bar plot based on the occurrences 
        of reads mappings for each reference.

    Parameters:
    - dict_ref_info (dict): A dictionary containing the
        information about the number of mappings for each read.
    - myoutfile (str): The folder for the output files.

    Example:
    - dict_ref_info = {'read1': 10, 'read2': 20, 'read3': 15}
    """
    occurrences = list(dict_ref_info.values())
    plt.hist(occurrences, bins=list(range(51)))

    x_ticks_positions = list(range(0, 51, 2))
    x_ticks_labels = list(range(0, 51, 2))
    plt.xticks(x_ticks_positions, x_ticks_labels)
    
    plt.xlabel('Countänd Nr. best locations')
    plt.ylabel('Count')
    plt.savefig(myoutfile+"reference_coverage.pdf")
    plt.show()
    
def calculate_reference_coverage(dict_read_ref_info, outfile,
                        read_sequnce_length,ref_sequnce_length):
    """
    Function:Calculate reference coverage based on read-reference 
        mapping information.

    Parameters:
    - dict_read_ref_info (dict): A dictionary containing information 
        about read-reference mappings.
    - outfile (str): The folder for the output files 
    - read_sequence_length (dict): A dictionary containing read IDs 
        as keys and corresponding sequence lengths as values.
    - ref_sequence_length (dict): A dictionary containing reference 
        IDs as keys and corresponding sequence lengths as values.

    Returns:
    - dict: A dictionary containing reference coverage values.

    Example:
    >>> dict_read_ref_info = {'read1': ['ref1', 'ref2'], 
                            'read2': ['ref2', 'ref3']}
    >>> m_i = {'read1':3, 'read2':3, 'read3':3, 
                'read4':4, 'read5':5}
    >>> m_ij = {
        'read1': {'ref1':1, 'ref2':1, 'ref3':1},
        'read2': {'ref1':1, 'ref2':1, 'ref3':1}}
    >>> c_j = {'ref1': 0.5, 'ref2': 0.5, 'ref3': 0.5}
    >>> read_sequence_length = {'read1': 4, 'read2': 6}
    >>> ref_sequence_length = {'ref1': 5, 'ref2': 4, 'ref3': 7}
    """
    c_j = {}
    readplus_m_ij_of_m_i = {}
    m_i = {}
    m_ij = {}
    
    for read_ref, ref_list in dict_read_ref_info.items():
        m_i[read_ref] = len(ref_list)
        m_ij.setdefault(read_ref, {})  
        for ref in ref_list:
            m_ij[read_ref][ref] = m_ij[read_ref].get(ref, 0) + 1

    for read, ref_counts in m_ij.items():
        for ref, count_ij in ref_counts.items():
            readplus_m_ij_of_m_i[ref] = readplus_m_ij_of_m_i.get(ref, 0) + \
                count_ij / m_i[read] * read_sequnce_length[read]
    for ref in readplus_m_ij_of_m_i:
        c_j[ref] = readplus_m_ij_of_m_i[ref] / ref_sequnce_length[ref]
    with open(outfile+'reference_coverage.txt', 'w') as output_file:
        for ref, c_j_value in c_j.items():
            output_file.write(f"{ref} , {c_j_value}\n")

    return c_j

def format_covert_full_mapping_result(full_mapping_results):
    """
    Function: convert full mapping results to a 
        simplified format.

    Parameters:
    - full_mapping_results (dict): A dictionary 
        containing full mapping results.

    Returns:
    - dict: A dictionary containing full mapping 
        results to a simplified format.
    """
    full_mapping_results_converted = {}
    for key, value in full_mapping_results.items():
            read_info, ref_info = key
            if read_info not in full_mapping_results_converted:
                full_mapping_results_converted[read_info] = {}

            for pos, pos_info in value.items():
                ref_start_position, mapping_start_position = pos
                mismatch_count = pos_info['mismatch_count']
                full_mapping_results_converted[read_info][ref_info] = [
                    mapping_start_position, mismatch_count]
    return full_mapping_results_converted

def filter_best_match(best_match_dict):
    """
    Function: Filter the full match dictionary to keep only entries 
        with the minimum or one of the minimum mismatch counts.

    Parameters:
    - best_match_dict (dict): A dictionary containing full match results 
        for each read-reference pair.

    Returns:
    - dict: A filtered dictionary containing entries with the minimum 
        or one of the minimum mismatch counts.
    """
    # all_min_match_counts stores the minimum mismatch_count of each read
    all_min_match_counts = {}
    for read_info, values in best_match_dict.items():
        for value in values.values():
            if read_info not in all_min_match_counts:
                all_min_match_counts[read_info] = value[1]
            else:
                if value[1] < all_min_match_counts[read_info]:
                    all_min_match_counts[read_info] = value[1]
    #Create a list to be deleted
    to_delete = []
    # Keep only the entries with the minimum or one of the minimum match_count
    for read_info, values in best_match_dict.items():
        for ref_info, value in values.items():
            if value[1] != all_min_match_counts[read_info]:
                to_delete.append((read_info, ref_info))

    # Delete non-minimum read-ref pairs
    for read_info, ref_info in to_delete:
        del best_match_dict[read_info][ref_info]

    return best_match_dict

def generate_outfile1(full_mapping_results,outfile):
    '''
    Function: Generate the first output file containing 
        formatted mapping information.

    Parameters:
    - full_mapping_results (dict): A dictionary containing 
        full mapping results.
    - outfile (str): The folder for the output files

    Returns:
    - tuple: A tuple containing dictionaries with information 
        about mappings.
    '''
    dict_ref_info = {}
    dict_read_ref_info = {}
    with open(outfile+"best_mapping.txt", "w") as output_file_1:
        
        full_mapping_results_converted = format_covert_full_mapping_result(
            full_mapping_results)
        full_mapping_results_converted = filter_best_match(
            full_mapping_results_converted)

        for read_info, value in full_mapping_results_converted.items():
            for ref_info, value in value.items():
                mapping_start_position = value[0]
                mismatch_count = value[1]
                if read_info not in dict_ref_info:
                    dict_ref_info[read_info] = 0
                    dict_read_ref_info[read_info] = []
                dict_ref_info[read_info] += 1
                dict_read_ref_info[read_info].append(ref_info)
              
                output_string = (f"{read_info},{ref_info},"
                    f"{mapping_start_position},{mismatch_count}\n")

                output_file_1.write(output_string)
    return dict_ref_info, dict_read_ref_info

def main():
    '''
    This is the main body of the program.
    '''
    
    parser = argparse.ArgumentParser(description='Read map')
    # Add parameters
    parser.add_argument('--k', type=int, 
        help='Size of k-mer for seeds generation and library comparison')
    parser.add_argument('infile_ref', type=str, 
        help='Input FASTA file containing reference sequences.')
    parser.add_argument('infile_read', type=str, 
        help='Input FASTA file containing read sequences.')
    parser.add_argument('outfile', type=str, 
        help=('folder for output files'))
    parser.add_argument('--match_threshold', type=float, 
        help=('Threshold of success matching rate for each '
            'k-mer (e.g., 0.8)'), default=0.7)
    args = parser.parse_args()
    
    # Data processing
    fasta_read = read_fasta_to_dict(args.infile_read,args.k)
    fasta_ref = read_fasta_to_dict(args.infile_ref,args.k)
    read_sequnce_length = calculate_length(fasta_read)
    ref_sequnce_length = calculate_length(fasta_ref)
    read_library = build_library(args.k, fasta_read)
    ref_library = build_library(args.k, fasta_ref)
    full_mapping_results = generate_seperated_dict_and_map(
        read_library,ref_library,args.k,args.match_threshold,
        read_sequnce_length,ref_sequnce_length)
    dict_ref_info, dict_read_ref_info = generate_outfile1(
        full_mapping_results, args.outfile)
    
    # Plot the histogram
    draw_bar_plot(dict_ref_info, args.outfile)
    calculate_reference_coverage(dict_read_ref_info,
        args.outfile,read_sequnce_length,ref_sequnce_length)

if __name__ == '__main__':
    main()
















