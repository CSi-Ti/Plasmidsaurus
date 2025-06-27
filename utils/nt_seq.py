import subprocess
from collections import Counter
import numpy as np
from pathlib import Path
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO

# Define motifs
kozak = "GCCACCATG"
tcra_motif = "ATCCAGAACCCTGACCCTGCCGT"
tcrb_motif = "TTTGAGCCATCAGAAGCAGA"
plx_r1 = "GTGGTTTAGTAA"
plx_f = "GCTAACTGTCGGGATCAACAAGTTTG"


def trim_plx(seq):
    start = seq.find(plx_f) + len(plx_f)
    end = seq.find(plx_r1)
    if start != -1 and end != -1 and end > start:
        return seq[start:end]
    return seq

def trim_kozak(seq):
    start = seq.find(kozak)
    if start != -1:
        return seq[start + 6:]
    else:
        return None

def trim_translate(seq):
    start = seq.find("ATG")
    if start == -1:
        return None
    stop_codons = ["TAA", "TAG", "TGA"]
    seq_length = len(seq)

    for i in range(start, seq_length - 2, 3):
        codon = seq[i:i + 3]
        if codon in stop_codons:
            return seq[start:i + 3]
    return None


def is_tcra(seq):
    return tcra_motif in seq

def is_tcrb(seq):
    return tcrb_motif in seq

def get_alignment_matrix_mafft(seqs):
    input_file = Path(__file__).parent.parent / 'data' / 'tmp.fasta'
    output_file = Path(__file__).parent.parent / 'data' / 'tmp.aln'
    records = [SeqRecord(seq, id=f"seq{i}", description="") for i, seq in enumerate(seqs)]
    SeqIO.write(records, input_file, "fasta")

    # subprocess.run(f'mafft --auto --thread 20 {input_file} > {output_file}', shell=True)
    subprocess.run(f'mafft --6merpair --retree 1 --maxiterate 0 --thread 20 {input_file} > {output_file}', shell=True)

    alignment = AlignIO.read(output_file, "fasta")
    aligned_strs = [str(rec.seq).upper() for rec in alignment]
    alignment_matrix = np.array([list(seq) for seq in aligned_strs])
    return alignment_matrix


def get_consensus_sequence(seqs, threshold=0.6):
    if len(seqs) == 0:
        return [], []
    # alignment_matrix = get_alignment_matrix(seqs)
    alignment_matrix = get_alignment_matrix_mafft(seqs)
    unmatched_mask = np.ones(len(seqs), dtype=bool)
    consensus_seqs = []
    consensus_matches = []
    match_percent = 1
    while match_percent > 0.01:
        min_count = threshold * np.sum(unmatched_mask)
        consensus = []
        for col in range(alignment_matrix.shape[1]):
            column = alignment_matrix[unmatched_mask, col]
            for k, v in Counter(column).most_common():
                if k == '-':
                    continue
                if v >= min_count:
                    consensus.append(k)
                    break
        consensus_seq = ''.join(consensus)
        consensus_seq = trim_translate(consensus_seq)
        if consensus_seq is None:
            break
        consensus_match_count = 0
        for i in range(len(alignment_matrix)):
            if unmatched_mask[i] and consensus_seq in seqs[i]:
                unmatched_mask[i] = False
                consensus_match_count += 1
        if consensus_match_count == 0 and len(consensus_seqs) > 0:
            break
        consensus_seqs.append(consensus_seq)
        consensus_matches.append(consensus_match_count)
        match_percent = consensus_match_count / len(seqs)
    return consensus_seqs, consensus_matches


def get_mutations(ref_seq, seq):
    # alignment_matrix = get_alignment_matrix([ref_seq, seq])
    alignment_matrix = get_alignment_matrix_mafft([ref_seq, seq])
    mutations = []
    aligned_ref_seq = alignment_matrix[0]
    aligned_seq = alignment_matrix[1]
    ins_count = 0
    for i in range(len(aligned_ref_seq)):
        if aligned_ref_seq[i] == '-':
            mutations.append(f'I{i - ins_count}{aligned_seq[i]}')
            ins_count += 1
            continue
        if aligned_seq[i] == '-':
            mutations.append(f'D{i}')
            continue
        if aligned_ref_seq[i] != aligned_seq[i]:
            mutations.append(f'M{i}{aligned_ref_seq[i]}{aligned_seq[i]}')
    return mutations

def translate_to_aa(seq, to_stop=True):
    seq = Seq(seq)
    aa_seq = seq.translate(to_stop=to_stop)
    return aa_seq

