
import pandas as pd

from utils.nt_seq import *

# Reader for ground-truth TCR sequences
def read_tcr_aa_seq(aa_ref_path):
    aa_df = pd.read_csv(aa_ref_path)
    tcra_aa_seqs = aa_df['TCR_alpha_full'].tolist()
    tcrb_aa_seqs = aa_df['TCR_beta_full'].tolist()
    return tcra_aa_seqs, tcrb_aa_seqs


def read_tcr_nt_seq(fastq_path, threshold=0.6):
    tcra_reads = []
    tcrb_reads = []

    total_reads, n_reads_with_kozak = 0, 0
    for record in SeqIO.parse(fastq_path, "fastq"):
        total_reads += 1
        seq = str(record.seq)
        rc_seq = str(Seq(seq).reverse_complement())

        for s in [seq, rc_seq]:
            trimmed_s = trim_plx(s)
            if trimmed_s is None:
                continue
            trimmed_s = trim_kozak(trimmed_s)
            if trimmed_s is None:
                continue
            n_reads_with_kozak += 1
            # trimmed_s = trim_translate(trimmed_s)
            # if trimmed_s is None:
            #     continue
            if is_tcra(trimmed_s):
                tcra_reads.append(Seq(trimmed_s))
            elif is_tcrb(trimmed_s):
                tcrb_reads.append(Seq(trimmed_s))

    tcra_cons_seqs, tcra_cons_matches = get_consensus_sequence(tcra_reads, threshold)
    tcrb_cons_seqs, tcrb_cons_matches = get_consensus_sequence(tcrb_reads, threshold)
    for i in range(len(tcra_cons_seqs)):
        print(f"TCR-alpha consensus sequence ({tcra_cons_matches[i] / len(tcra_reads) * 100:.2f}%): {tcra_cons_seqs[i]}")
    for i in range(len(tcrb_cons_seqs)):
        print(f"TCR-beta consensus sequence ({tcrb_cons_matches[i] / len(tcrb_reads) * 100:.2f}%): {tcrb_cons_seqs[i]}")

    return (total_reads, n_reads_with_kozak, len(tcra_reads), tcra_cons_seqs, tcra_cons_matches,
            len(tcrb_reads), tcrb_cons_seqs, tcrb_cons_matches)

def run(fastq_folder, aa_ref_path):
    tcra_aa_seqs, tcrb_aa_seqs = read_tcr_aa_seq(aa_ref_path)
    fastq_files = sorted(list(Path(fastq_folder).glob('*.fastq')))
    header = ['File_name', 'Total_reads', 'Kozak_reads', 'Type', 'Type_reads', 'Lib_seq', 'Cons_seq', 'Match', 'N_match', 'Percentage']
    rows = []
    for fastq_file in fastq_files:
        # if 'SMTC75_14_sample_14' not in fastq_file.stem:
        #     continue
        index = int(fastq_file.stem.split('_')[-1]) - 1
        tcra_aa_seq, tcrb_aa_seq = tcra_aa_seqs[index], tcrb_aa_seqs[index]

        (total_reads, n_reads_with_kozak, n_tcra, tcra_cons_seqs, tcra_cons_matches,
         n_tcrb, tcrb_cons_seqs, tcrb_cons_matches) = read_tcr_nt_seq(fastq_file, threshold=0.3)

        tcra_cons_aa_seqs = [translate_to_aa(seq, to_stop=True) for seq in tcra_cons_seqs]
        tcrb_cons_aa_seqs = [translate_to_aa(seq, to_stop=True) for seq in tcrb_cons_seqs]

        print(header)
        for tcra_cons_aa_seq, tcra_cons_match in zip(tcra_cons_aa_seqs, tcra_cons_matches):
            row = [fastq_file.stem, total_reads, n_reads_with_kozak, 'TCR-alpha', n_tcra, tcra_aa_seq,
                   str(tcra_cons_aa_seq), tcra_cons_aa_seq==tcra_aa_seq, tcra_cons_match, f'{tcra_cons_match / n_tcra * 100:.2f}']
            print(row)
            rows.append(row)
        for tcrb_cons_aa_seq, tcrb_cons_match in zip(tcrb_cons_aa_seqs, tcrb_cons_matches):
            row = [fastq_file.stem, total_reads, n_reads_with_kozak, 'TCR-beta', n_tcrb, tcrb_aa_seq,
                   str(tcrb_cons_aa_seq), tcrb_cons_aa_seq==tcrb_aa_seq, tcrb_cons_match, f'{tcrb_cons_match / n_tcrb * 100:.2f}']
            print(row)
            rows.append(row)
    result_df = pd.DataFrame(rows, columns=header)
    result_df.to_csv('./data/result.csv', index=False)


if __name__ == '__main__':
    run(fastq_folder='./data/fastq', aa_ref_path='./data/MS1089_out TCR primers_correct reference sequence.csv')