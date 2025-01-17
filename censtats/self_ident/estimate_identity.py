# ModDotPlot
# https://github.com/marbl/ModDotPlot/commit/50ecda4eff91acd00584090afd380d4a355be7aa
import math
import numpy as np
import mmh3


def removeAmbiguousBases(mod_list: set[int], k: int) -> set[int]:
    # Ambiguous IUPAC codes
    bases_to_remove = ["R", "Y", "M", "K", "S", "W", "H", "B", "V", "D", "N"]
    kmers_to_remove = set()
    for i in range(len(bases_to_remove)):
        result_string = str(bases_to_remove[i]) * k
        kmers_to_remove.add(mmh3.hash(result_string))
    mod_set = set(mod_list)
    # Remove homopolymers of ambiguous nucleotides
    mod_set.difference_update(kmers_to_remove)
    return mod_set


def createSelfMatrix(
    sequence: list[int],
    window_size: int,
    delta: float,
    k: int,
    identity: float,
    ambiguous: bool,
    modimizer: int,
) -> np.ndarray:
    """
    Create self-identity matrix.

    Args:
    * sequence
            * Sequence as list of mmh3 kmers.
    * window_size
            * Window size.
    * delta
            * Fraction of neighboring partition to include in identity estimation. Must be between 0 and 1, use > 0.5 is not recommended.
    * k
            * kmer length
    * identity
            * Identity cutoff threshold.
    * ambiguous
            * Preserve diagonal when handling strings of ambiguous homopolymers (eg. long runs of N's).
    * modimizer
            * Modimizer sketch size.
            * A lower value will reduce the number of modimizers, but will increase performance.
            * Must be less than window size.

    Returns
    * 2D ndarray of identity values.
    """
    sequence_length = len(sequence)
    seq_sparsity = round(window_size / modimizer)
    if seq_sparsity <= modimizer:
        seq_sparsity = 2 ** int(math.log2(seq_sparsity))
    else:
        seq_sparsity = 2 ** (int(math.log2(seq_sparsity - 1)) + 1)
    sketch_size = round(window_size / seq_sparsity)

    no_neighbors = partitionOverlaps(sequence, window_size, 0, sequence_length, k)
    if delta > 0:
        neighbors = partitionOverlaps(sequence, window_size, delta, sequence_length, k)
    else:
        neighbors = no_neighbors

    neighbors_mods = convertToModimizers(
        neighbors, seq_sparsity, ambiguous, k, sketch_size
    )
    no_neighbors_mods = convertToModimizers(
        no_neighbors, seq_sparsity, ambiguous, k, sketch_size
    )
    matrix = selfContainmentMatrix(
        no_neighbors_mods, neighbors_mods, k, identity, ambiguous
    )
    return matrix


def partitionOverlaps(
    lst: list[int], win: int, delta: float, seq_len: int, k: int
) -> list[list[int]]:
    kmer_list = []
    kmer_to_genomic_coordinate_offset = win - k + 1
    delta_offset = win * delta

    # set the first window to contain win - k + 1 kmers.
    starting_end_index = int(round(kmer_to_genomic_coordinate_offset + delta_offset))
    kmer_list.append(lst[0:starting_end_index])
    counter = win - k + 1

    # set normal windows
    while counter <= (seq_len - win):
        start_index = counter + 1
        end_index = win + counter + 1
        delta_start_index = int(round(start_index - delta_offset))
        delta_end_index = int(round(end_index + delta_offset))
        if delta_end_index > seq_len:
            delta_end_index = seq_len
        try:
            kmer_list.append(lst[delta_start_index:delta_end_index])
        except Exception as e:
            print(e)
            kmer_list.append(lst[delta_start_index:seq_len])
        counter += win

    # set the last window to get the remainder
    if counter <= seq_len - 2:
        final_start_index = int(round(counter + 1 - delta_offset))
        kmer_list.append(lst[final_start_index:seq_len])

    # Test that last value was added on correctly

    assert kmer_list[-1][-1] == lst[-1]
    return kmer_list


def populateModimizers(
    partition: list[int], sparsity: int, ambiguous: bool, expectation: int, k: int
) -> set[int]:
    mod_set = set()
    for kmer in partition:
        if kmer % sparsity == 0:
            mod_set.add(kmer)
    if not ambiguous:
        mod_set = removeAmbiguousBases(mod_set, k)
    if (len(mod_set) < round(expectation / 2)) and (sparsity > 1):
        populateModimizers(partition, sparsity // 2, ambiguous, expectation, k)
    return mod_set


def convertToModimizers(
    kmer_list: list[list[int]], sparsity: int, ambiguous: bool, k: int, expectation: int
) -> list[set[int]]:
    mod_total = []
    for partition in kmer_list:
        mod_set = populateModimizers(partition, sparsity, ambiguous, expectation, k)
        mod_total.append(mod_set)
    return mod_total


def convertMatrixToBed(
    matrix: np.ndarray,
    window_size: int,
    id_threshold: float,
    x_name: str,
    y_name: str,
    self_identity: bool,
) -> list[tuple[str, int, int, str, int, int, float]]:
    """
    # Returns
    * "#query_name"
    * "query_start"
    * "query_end"
    * "reference_name"
    * "reference_start"
    * "reference_end"
    * "perID_by_events"
    """
    bed = []
    rows, cols = matrix.shape
    for x in range(rows):
        for y in range(cols):
            value = matrix[x, y]
            if (not self_identity) or (self_identity and x <= y):
                if value >= id_threshold / 100:
                    start_x = x * window_size + 1
                    end_x = (x + 1) * window_size
                    start_y = y * window_size + 1
                    end_y = (y + 1) * window_size

                    bed.append(
                        (
                            x_name,
                            int(start_x),
                            int(end_x),
                            y_name,
                            int(start_y),
                            int(end_y),
                            float(value),
                        )
                    )
    return bed


def binomial_distance(containment_value: float, kmer_value: int) -> float:
    """
    Calculate the binomial distance based on containment and kmer values.

    Args:
        containment_value (float): The containment value.
        kmer_value (int): The k-mer value.

    Returns:
        float: The binomial distance.
    """
    return math.pow(containment_value, 1.0 / kmer_value)


def containment_neighbors(
    set1: set[int],
    set2: set[int],
    set3: set[int],
    set4: set[int],
    identity: float,
    k: int,
) -> float:
    """
    Calculate the containment neighbors based on four sets and an identity threshold.

    Args:
    * set1:
            * The first set.
    * set2:
            * The second set.
    * set3:
            * The third set.
    * set4:
            * The fourth set.
    * identity:
            * The identity threshold.
    * k:
            * Kmer value.

    Returns:
    * The containment neighbors value.
    """
    len_a = len(set1)
    len_b = len(set2)

    intersection_a_b_prime = len(set1 & set4)
    if len_a != 0:
        containment_a_b_prime = intersection_a_b_prime / len_a
    else:
        # If len_a is zero, handle it by setting containment_a_b_prime to a default value
        containment_a_b_prime = 0

    if binomial_distance(containment_a_b_prime, k) < identity / 100:
        return 0.0

    else:
        intersection_a_prime_b = len(set2 & set3)
        if len_b != 0:
            containment_a_prime_b = intersection_a_prime_b / len_b
        else:
            # If len_a is zero, handle it by setting containment_a_b_prime to a default value
            containment_a_prime_b = 0

        return max(containment_a_b_prime, containment_a_prime_b)


def selfContainmentMatrix(
    mod_set: list[set],
    mod_set_neighbors: list[set],
    k: int,
    identity: float,
    ambiguous: bool,
) -> np.ndarray:
    """
    Create a self-containment matrix based on containment similarity calculations.

    Args:
    * mod_set
            * A list of sets representing elements.
    * mod_set_neighbors
            * A list of sets representing neighbors for each element.
    * k
            * A parameter for containment similarity calculation.

    Returns:
        np.ndarray: A NumPy array representing the self-containment matrix.
    """
    n = len(mod_set)
    containment_matrix = np.empty((n, n))

    for w in range(n):
        containment_matrix[w, w] = 100.0
        if len(mod_set[w]) == 0 and not ambiguous:
            containment_matrix[w, w] = 0

        for r in range(w + 1, n):
            c_hat = binomial_distance(
                containment_neighbors(
                    mod_set[w],
                    mod_set[r],
                    mod_set_neighbors[w],
                    mod_set_neighbors[r],
                    identity,
                    k,
                ),
                k,
            )
            containment_matrix[r, w] = c_hat * 100.0
            containment_matrix[w, r] = c_hat * 100.0

    return containment_matrix
