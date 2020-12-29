"""
Identify matching haplotypes
"""
import argparse
from Bio import AlignIO

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Identify matching haplotypes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--alignment', type=str, metavar="FASTA", required=True, help="input FASTA alignment")
    args = parser.parse_args()

    aln = AlignIO.read(args.alignment, "fasta")
    ref = [a for a in aln if a.id=='Wuhan/Hu-1/2019'][0]

    ref_expected = ['G', 'C', 'C']
    positions = [28883, 913, 5388]
    wt_expected = ['C', 'T', 'A']

    for i in range(0, len(positions)):
        assert(ref[positions[i]-1].upper() == ref_expected[i])

    newwave_records = []
    for seq in aln:
        include = True
        for i in range(0, len(positions)):
            if seq[positions[i]-1].upper() != wt_expected[i]:
                include = False
        if include:
            newwave_records.append(seq.name)

    newwave_records.sort()
    for record in newwave_records:
        print(record)
