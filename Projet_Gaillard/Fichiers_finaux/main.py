
import argparse
import methods as mt 

def parse_args():
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-align', default='NW_N', type=str, help="name of NW type to run, either 'NW_N' for multiple comparison or 'NW_struc' for structure information")
    parser.add_argument('-namefile', default='BBS11001.fasta', type=str, help="fasta file name to study alignments")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    BLOSUM = {}
    mt.blosum()
    sequences = mt.ouvertureFASTA("balibase/RV11.unaligned/"+args.namefile)
    #TODO : ne pas oublier de rajouter les ID pour pouvoir enregistrer en FASTA
    alignments=mt.alignNseq(sequences,args.align)
    #TODO écrire save_align_FASTA qui sauve les alignments en FASTA
    #save_align_FASTA(name_test_file,alignments)
    mt.bali_score("balibase/RV11.unaligned/"+args.namefile,name_test_file)     


if __name__ == "__main__":
    main()
