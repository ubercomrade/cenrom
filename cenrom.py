import os
import os.path
import sys
import argparse
from lib.pwm import PWM
from lib.functions import read_fasta, shuffle_fasta, get_number_of_sites


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta', action='store', help='path to FASTA file')
    parser.add_argument('matrix', action='store', help='path to PCM file')
    parser.add_argument('promoters', action='store', choices=['mm10', 'hg38', 'tair10'], metavar='N',
         help='promoters of organism (hg38, mm10, tair10)')
    parser.add_argument('-n', '--ntimes', action='store', default=100, type=int, dest='ntimes',
        help='N times FASTA for background')
    parser.add_argument('-f', '--format', action='store', choices=['homer', 'cisbp', 'hocomoco'], metavar='N',
        dest='format', default='homer', help='[homer, cisbp, hocomoco] format of input PFM (HOMER, CISBP) or PCM (HOCOMOCO)')
    parser.add_argument('-t', '--threshold', action='store', type=float, dest='threshold',
                        required=False, default=1.9*10**(-4), help='threshold based on FPR, def=1.9*10^(-4)')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    fasta_path = args.fasta
    matrix_path = args.matrix
    organism = args.promoters
    ntimes = args.ntimes
    matrix_format = args.format
    fpr_threshold = args.threshold
    
    this_dir, this_filename = os.path.split(__file__)
    if organism == 'mm10':
        path_to_promoters = os.path.join(this_dir, "promoters", "mm10.fasta")
    elif organism == 'hg38':
        path_to_promoters = os.path.join(this_dir, "promoters", "hg38.fasta")
    elif organism == 'tair10':
        path_to_promoters = os.path.join(this_dir, "promoters", "tair10.fasta")

    fasta = read_fasta(fasta_path)
    background = shuffle_fasta(fasta, ntimes)
    pwm = PWM(matrix_path, form=matrix_format)
    pwm.calculate_table(path_to_promoters)
    threshold = pwm.choose_threshold(fpr_threshold)
    scores = pwm.calculate_scores_upper_threshold(fasta, threshold)
    number_of_sites = get_number_of_sites(fasta, pwm.length)
    real_fraction = len(scores) / number_of_sites
    background_scores, background_number_of_sites = pwm.calculate_scores_upper_threshold(background, threshold)
    background_fraction = len(background_scores) / background_number_of_sites
    enrichment = real_fraction / background_fraction
    print(f'{real_fraction}\t{background_fraction}\t{enrichment}')

if __name__ == '__main__':
    main()
