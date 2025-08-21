import os
import sys
import glob
from multiprocessing import Pool, cpu_count
from datetime import datetime

def process_fasta(file_path):
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            if not lines or not lines[0].startswith('>'):
                return None
            acc = lines[0].strip().split()[0][1:]  # Remove '>'
            seq = ''.join(line.strip() for line in lines[1:] if not line.startswith('>'))
            if acc and seq:
                return acc, '?', '1', seq
    except Exception:
        return None
def main(fasta_dir):
    fasta_files = glob.glob(os.path.join(fasta_dir, '*.fasta'))
    print(f"Found {len(fasta_files)} FASTA files.")

    with Pool(cpu_count()) as pool:
        results = pool.map(process_fasta, fasta_files)

    valid_results = sorted([r for r in results if r], key=lambda x: x[0])

    # save output in the 'data/' directory
    output_dir = 'data'
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"gb_mitogenomes_{datetime.today().strftime('%d-%m-%y')}.fst")

    with open(output_file, 'w') as out:
        out.write("#! ALL\n")
        for acc, qmark, one, seq in valid_results:
            out.write(f"{acc}\t{qmark}\t{one}\t{seq}\n")

    print(f"Wrote {len(valid_results)} records to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python fasta_dir_to_tsv.py /path/to/fasta_dir")
        sys.exit(1)
    main(sys.argv[1])
