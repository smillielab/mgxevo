import os
import numpy as np
import pandas as pd
import pickle
import glob
import argparse

def get_bowtie_abundances(species_path, samples_path, output_path):
    # get species
    fl = [os.path.basename(file).replace('.aln.pickle', '') for file in glob.glob(f"{species_path}/*.pickle")]
    # get samples
    with open(samples_path, 'r') as f:
        samples = f.read().splitlines()
    # create counts (technically coverage matrix)
    counts = pd.DataFrame(0, index=samples, columns=fl)
    # unpickle and add to counts matrix
    for j, fj in enumerate(fl):
        with open(f"{species_path}/{fj}.aln.pickle", 'rb') as f:
            x = pickle.load(f)
        try:
            counts.iloc[x['i'],j] = x['x'].sum(axis=2).sum(axis=1)
        except:
            print(fj)
            break
    # save to file
    counts.to_csv(f"{output_path}/bowtie_abundances.csv", index=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Get bowtie abundances.')
    parser.add_argument('--species_path', required=True, help='Path to species files.')
    parser.add_argument('--samples_path', required=True, help='Path to samples file.')
    parser.add_argument('--output_path', required=True, help='Path to output file.')
    args = parser.parse_args()

    get_bowtie_abundances(args.species_path, args.samples_path, args.output_path)