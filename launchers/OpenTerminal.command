#!/bin/bash
# Manual launcher
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT" || exit 1
if [ ! -d "data/references" ] || [ ! -d "data/reads" ]; then
    mkdir -p data/references data/reads
    echo "Input folders created. Put reference FASTA files in data/references and reads in data/reads."
fi
echo "======================================================================================="
echo "Välkommen till 'short reads mapping program!'"
echo "*****enter 'python3 scripts/main.py -h' for help"
echo "*****Command example and explanation:*****"
echo "python3 scripts/main.py --k 15 refs0.fa reads0.fa Refs0_Reads0"
echo "python3 scripts/main.py                                ---interpreter and main function"
echo "                --k 15                                 ---length of kmer (default 10)"
echo "                       refs1.fa                        ---reference fasta filename"
echo "                                reads1.fa              ---read fasta filename"
echo "                                          Refs1_Reads1 ---output file prefix"
echo "======================================================================================="
echo "!!!!!Important NOTICES!!!!!"
echo "Filename extensions for the input files are needed"
echo "Do NOT use invalid character for the output file prefix! E.g. '/' or  '.'"
echo "Length of kmer should be no larger than 32 Preferably 10~30!"
echo "======================================================================================="
echo "*****About INPUT and OUTPUT file*****"
echo "Put reference FASTA files in 'data/references' and reads in 'data/reads'"
echo "All results will be stored in 'results/output'"
echo "======================================================================================="
exec $SHELL



