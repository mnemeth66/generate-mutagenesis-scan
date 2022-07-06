from pathlib import Path
from enum import Enum
from typing import Union, Optional
import csv

from latch import small_task, workflow
from latch.types import LatchFile, LatchDir
from Bio import SeqIO

class Methods(Enum):
    alanine = 'alanine'
    deep = 'deep'

alphabet = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

@small_task
def create_vars(
    sequence: Union[str, LatchFile],
    output_dir: Optional[LatchDir],
    range_start: int,
    range_end: int,
    method: Methods, 
) -> LatchFile:
    """
    Create variants from the input file.
    """
    output_filename = 'variants.csv'
    if isinstance(sequence, LatchFile):
        local_sequence = Path(sequence).resolve()
        output_filename = local_sequence.stem + '_variants.csv'
        for record in SeqIO.parse(local_sequence, 'fasta'):
            sequence = str(record.seq)
            break
    vars = []

    # Generate the scan
    # Currently only doing single-site mutagenesis: alanine scanning and deep scanning
    if range_start < 0:
        range_start += len(sequence)
    range_start = max(0, range_start)
    if range_end < 0:
        range_end += len(sequence)
    range_end = min(len(sequence) - 1, range_end)
    assert range_start <= range_end, 'range_start must be less than range_end'
    for i in range(range_start, range_end + 1):
        aa = sequence[i]
        if method == Methods.alanine:
            vars.append(aa + str(i) + 'A')
        elif method == Methods.deep:
            for a in alphabet:
                if a != aa:
                    vars.append(aa + str(i) + a)
                    
    # Create a csv with vars in the 'mutant' column)
    with open(f'/root/{output_filename}', 'w') as f:
        writer = csv.writer(f)
        writer.writerow(['mutant'])
        for var in vars:
            writer.writerow([var])
    output_dir = output_dir.remote_path if output_dir else 'latch:///'
    return LatchFile(f'/root/{output_filename}', f'{output_dir}/{output_filename}')

@workflow
def generate_scan(
    sequence: Union[str, LatchFile],
    range_start: int = 0,
    range_end: int = -1,
    method: Methods = Methods.deep,
    output_dir: Optional[LatchDir] = None,
) -> LatchFile:
    """
    Generate a Mutagenesis Scan 
    ----
    # Mutational Scanning
    Mutational scanning is a set of methods developed to analyzed the effects
    of single of combinatorial amino acid mutations on protein structure and 
    function. 
    Historically, this began with alanine scanning methods, in which experiments
    would replace amino acids one at a time to alanine, a relatively inert amino 
    acid. However, as sequencing methods became more powerful, so called deep
    mutational scanning methods were developed, in which biologists could 
    feasibly test the effects of much larger libraries of single-site mutations,
    not just restricted to alanine.

    # This Workflow
    This workflow is a the start to a computational mutational scanning workflow.
    It takes a fasta file or string representing a protein sequence as an input,
    and allows the user to initiate a mutational scanning workflow, by selecting
    their mutagenesis method of choice.
    Currently, the options are only alanine scanning and full deep mutational scanning
    over a specified range, and the workflow simply outputs a csv of variants to use
    in downstream workflows.

    ## Inputs
    - sequence: A fasta file or string representing a protein sequence.
    - range: The range over the protein to mutate.
        - This is 0 indexed, so the first amino acid is 0, the second is 1, etc.
    - range_start: The start of the range to mutate. Default is 0.
    - range_end: The end of the range to mutate. Default is the length of the sequence. Negative indexing supported.
    - method: The method of mutagenesis to use.
        - Alanine scanning 
        - Deep mutational scanning. 
    - output_dir: The directory to output the csv file to.

    ## Outputs
    - {sequence}_scan.csv: A csv file containing the variants created

    ## Other Tools
    - [dms-view.github.io](DMS-View is a site used for visualizing the results of a DMS Scan)

    __metadata__:
        display_name: Generate Mutagenesis Scan
        author:
            name: Matthew Nemeth
            email: mnemeth6@berkeley.edu
            github:
        repository: 
        license:
            id: MIT

    Args:
        sequence: 
            A fasta file or string representing a protein sequence.
            __metadata__:
                display_name: Sequence
        range_start:
            The range over the protein to mutate.
            __metadata__:
                display_name: Range Start
                description: 0-indexed.
        range_end:
            The end of the range to mutate.
            __metadata__:
                display_name: Range End
                description: Can use negative index. Must be larger than range_start (after accounting for negative indexing).
        method:
            The method of mutagenesis to use.
            __metadata__:
                display_name: Method
        output_dir:
            The directory to output the csv file to.
            __metadata__:
                display_name: Output Directory
    """
    return create_vars(
        sequence = sequence,
        output_dir = output_dir,
        range_start = range_start,
        range_end = range_end,
        method = method
    )
