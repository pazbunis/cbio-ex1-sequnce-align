import argparse
import numpy as np
from itertools import groupby
import csv


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


def align_global(s, t, score_matrix):
    print(s)
    print(t)
    s = '-' + s
    t = '-' + t
    num_rows = len(s)
    num_cols = len(t)
    v = np.zeros([num_rows, num_cols])
    p = np.zeros([num_rows, num_cols])
    # init left column
    for i in range(1, num_rows):
        v[i, 0] = v[i-1, 0] + score_matrix.score(s[i], '-')
    # init top row
    for j in range(1, num_cols):
        v[0, j] = v[0, j-1] + score_matrix.score('-', t[j])

    # dynamic program:
    for col in range(1, num_cols):
        for row in range(1, num_rows):
            vals = (v[row-1, col-1] + score_matrix.score(s[row], t[col]),
                    v[row-1, col] + score_matrix.score(s[row], '-'),
                    v[row, col-1] + score_matrix.score('-', t[col]))
            max_val = max(vals)
            p[row, col] = np.argmax(vals)
            v[row, col] = max_val

    # print(v)
    # print(p)

    # find optimal alignment:
    row = num_rows - 1
    col = num_cols - 1
    top = []
    bottom = []
    while row > 0 and col > 0:
        if p[row, col] == 0:
            top.append(s[row])
            bottom.append(t[col])
            row -= 1
            col -= 1
        elif p[row, col] == 1:
            top.append(s[row])
            bottom.append('-')
            row -= 1
        elif p[row, col] == 2:
            top.append('-')
            bottom.append(t[col])
            col -= 1

    print(row)
    print(col)
    for n in range(1, row+1):
        top.append(s[n])
        bottom.append('-')
    for n in range(1, col+1):
        top.append('-')
        bottom.append(t[n])


    top_str = ''.join(top[::-1])
    bottom_str = ''.join(bottom[::-1])

    print(top_str)
    print(bottom_str)
    print(v[num_rows-1, num_cols-1])
    return top_str, bottom_str, v[num_rows-1, num_cols-1]


def align_local(seq_a, seq_b, score_matrix):
    print('local')
    print(seq_a)
    print(seq_b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('seq_a', help='Path to first FASTA file (e.g. fastas/HomoSapiens-SHH.fasta)')
    parser.add_argument('seq_b', help='Path to second FASTA file')
    parser.add_argument('--align_type', help='Alignment type (e.g. local)', required=True)
    parser.add_argument('--score', help='Score matrix in.tsv format (default is score_matrix.tsv) ', default='score_matrix.tsv')
    command_args = parser.parse_args()
    header_a, seq_a = list(n for n in fastaread(command_args.seq_a))[0]  # first sequence in seq_a
    header_b, seq_b = list(n for n in fastaread(command_args.seq_b))[0]  # first sequence in seq_b
    score_matrix = ScoreMatrix(command_args.score)
    if command_args.align_type == 'global':
        top, bottom, value = align_global(seq_a, seq_b, score_matrix)
    elif command_args.align_type == 'local':
        alignment = align_local(seq_a, seq_b, score_matrix)

if __name__ == '__main__':
    main()
