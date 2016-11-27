import argparse
import numpy as np
from itertools import groupby
import csv
from scipy.constants.codata import value


def fastaread(fasta_name):
    f = open(fasta_name)
    faiter = (x[1] for x in groupby(f, lambda line: line.startswith(">")))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


class ScoreMatrix:
    def __init__(self, matrix_path):
        with open(matrix_path, 'r') as tsvin:
            tsvin = csv.reader(tsvin, delimiter='\t')
            self.raw_matrix = []
            for row in tsvin:
                self.raw_matrix.append(row)

    def score(self, letter_a, letter_b):
        letters2idx = {'A': 1, 'C': 2, 'G': 3, 'T': 4, '-': 5}
        return int(self.raw_matrix[letters2idx[letter_a]][letters2idx[letter_b]])


def align(s, t, score_matrix, align_type):
    s = '-' + s
    t = '-' + t
    num_rows = len(s)
    num_cols = len(t)
    v = np.zeros([num_rows, num_cols])
    p = np.zeros([num_rows, num_cols])
    # init left column
    for i in range(1, num_rows):
        v[i, 0] = v[i-1, 0] + score_matrix.score(s[i], '-') if align_type == 'global' else 0
    # init top row
    for j in range(1, num_cols):
        v[0, j] = v[0, j-1] + score_matrix.score('-', t[j]) if align_type == 'global' else 0

    # dynamic program:
    restart_penalty = np.iinfo(np.int32).min if align_type == 'global' else 0
    for col in range(1, num_cols):
        for row in range(1, num_rows):
            vals = (v[row-1, col-1] + score_matrix.score(s[row], t[col]),
                    v[row-1, col] + score_matrix.score(s[row], '-'),
                    v[row, col-1] + score_matrix.score('-', t[col]),
                    restart_penalty)
            max_val = max(vals)
            p[row, col] = np.argmax(vals)
            v[row, col] = max_val

    # find optimal alignment:
    max_row = row = num_rows - 1 if align_type == 'global' else np.argmax(np.max(v, axis=1))
    max_col = col = num_cols - 1 if align_type == 'global' else np.argmax(np.max(v, axis=0))
    top = []
    bottom = []
    while row > 0 and col > 0:
        if p[row, col] == 0:
            top.insert(0, s[row])
            bottom.insert(0, t[col])
            row -= 1
            col -= 1
        elif p[row, col] == 1:
            top.insert(0, s[row])
            bottom.insert(0, '-')
            row -= 1
        elif p[row, col] == 2:
            top.insert(0, '-')
            bottom.insert(0, t[col])
            col -= 1
        elif p[row, col] == 3:  # should only happen with local alignment
            if align_type == 'global':
                raise Exception('Tried to restart during global alignment - very bad!')
            break

    # add missing letter alignments and get final alignments
    if align_type == 'global':
        top_str = '-'*col + s[1:row+1] + ''.join(top)
        bottom_str = '-'*row + t[1:col+1] + ''.join(bottom)
    if align_type == 'local':
        top_str = ''.join(top)
        bottom_str = ''.join(bottom)

    return top_str, bottom_str, v[max_row, max_col]


def create_blocks(top, bottom):
    block_length = 50
    blocks = []
    for n in range(0, len(top), block_length):
        blocks.append(top[n:n + block_length] + '\n' + bottom[n: n + block_length] + '\n')
    return '\n'.join(blocks)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    command_args = parser.parse_args()
    header_a, seq_a = list(n for n in fastaread(command_args.seq_a))[0]  # first sequence in seq_a
    header_b, seq_b = list(n for n in fastaread(command_args.seq_b))[0]  # first sequence in seq_b

    # create a utility ScoreMatrix for easily obtaining scores
    score_matrix = ScoreMatrix(command_args.score)
    # get an optimal alignment
    top, bottom, al_val = align(seq_a, seq_b, score_matrix, command_args.align_type)
    # convert al_val to int if possible (so we print "3" and not "3.0")
    al_val = int(al_val) if (np.round(al_val) == al_val) else al_val
    # blockify the alignment
    final_output = create_blocks(top, bottom)
    # output
    print(final_output)
    print('Score: ' + str(al_val))
if __name__ == '__main__':
    main()
