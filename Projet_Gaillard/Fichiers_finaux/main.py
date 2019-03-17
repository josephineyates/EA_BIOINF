
import argparse
import methods as mt 

def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-align', default='NW_N', type=str, help="name of NW type to run, either 'NW_N' for multiple comparison or 'NW_struc' for structure information")
    parser.add_argument('-namefile', default='BBS11001.fasta', type=str, help="fasta file name to study alignments, select 'all' to get a global score compared to ClustalW")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    #TODO : ne pas oublier de rajouter les ID pour pouvoir enregistrer en FASTA
    if args.align == "NW_N":
        if args.namefile == "all":
            mt.eval_clustalw("score")
        else:
            sp_, tc_ = mt.align_score(args.namefile, algo="NW_blosum", silent=False)
    elif args.align == "NW_struct":
        if args.namefile == "all":
            mt.eval_clustalw("score_NW_struct")
        else:
            sp_, tc_ = mt.align_score(args.namefile, algo="NW_structural", silent=False)  


if __name__ == "__main__":
    main()
