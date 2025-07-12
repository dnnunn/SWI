#!/usr/bin/env python3
"""
SWI Predictor Batch Wrapper
Standardizes output for benchmarking solubility predictors.

Usage:
  python swi_predict_wrapper.py --fasta <input_fasta> --out <output_csv>

Outputs CSV with columns:
  Accession, Sequence, Predictor, SolubilityScore, Probability_Soluble, Probability_Insoluble
"""
import argparse
import os
import numpy as np
import pandas as pd
from Bio import SeqIO

def fasta_reader(file):
    records = list(SeqIO.parse(file, "fasta"))
    if not records:
        raise ValueError("No sequences found in FASTA file.")
    data = {
        "Accession": [rec.id for rec in records],
        "Sequence": [str(rec.seq).upper().replace('U', 'C') for rec in records]
    }
    df = pd.DataFrame(data)
    return df

# SWI weights and logistic fit constants
weights = {'A': 0.8356471476582918,
           'C': 0.5208088354857734,
           'E': 0.9876987431418378,
           'D': 0.9079044671339564,
           'G': 0.7997168496420723,
           'F': 0.5849790194237692,
           'I': 0.6784124413866582,
           'H': 0.8947913996466419,
           'K': 0.9267104557513497,
           'M': 0.6296623675420369,
           'L': 0.6554221515081433,
           'N': 0.8597433107431216,
           'Q': 0.789434648348208,
           'P': 0.8235328714705341,
           'S': 0.7440908318492778,
           'R': 0.7712466317693457,
           'T': 0.8096922697856334,
           'W': 0.6374678690957594,
           'V': 0.7357837119163659,
           'Y': 0.6112801822947587}
A = 81.0581
B = -62.7775

def compute_swi(sequence):
    return np.mean([weights[aa] for aa in sequence if aa in weights])

def main():
    parser = argparse.ArgumentParser(description="Batch SWI predictor wrapper")
    parser.add_argument('--fasta', '-f', required=True, help='Input FASTA file')
    parser.add_argument('--out', '-o', required=True, help='Output CSV file')
    args = parser.parse_args()

    df = fasta_reader(args.fasta)
    df['SolubilityScore'] = df['Sequence'].apply(compute_swi)
    df['Probability_Soluble'] = 1 / (1 + np.exp(-(A * df['SolubilityScore'] + B)))
    df['Probability_Insoluble'] = 1 - df['Probability_Soluble']
    df['Predictor'] = 'SWI'
    df = df[['Accession', 'Sequence', 'Predictor', 'SolubilityScore', 'Probability_Soluble', 'Probability_Insoluble']]
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    df.to_csv(args.out, index=False)
    print(f"Results written to {args.out}")

if __name__ == '__main__':
    main()
