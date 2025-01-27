import argparse
import os
import sys
import random
import string
import pandas as pd
from Bio import SeqIO
import glob
import subprocess

length = 10


# cmd
# blastn -query $uhgv -db hg -num_threads $c -evalue 0.001 -perc_identity 70 -out ${sname}.tsv
# -outfmt "6 qseqid sseqid pident qlen slen length sstart send"
def break_contig(contigs, cname, start, end):
    """
    Returns:
        tuple: Two sequences representing the broken contig.
    """
    if cname not in contigs:
        raise ValueError(f"Contig '{cname}' not found.")

    contig = contigs[cname]

    if start < 1 or end > len(contig) or start > end:
        raise ValueError("Invalid start or end positions.")

    # Convert to 0-based indexing for slicing
    start -= 1
    end -= 1

    # Break the contig into two parts
    part1 = contig[:start]
    part2 = contig[end + 1:]

    return part1, part2


def merge_intervals(coord_df, dist=0):
    """
    Returns:
        data.frame: the extended coordinates intervals if they overlap.
    """
    if len(coord_df.index) == 0:
        return coord_df[['sstart', 'send']]

    coord_df = coord_df.sort_values(by="sstart").reset_index(drop=True)

    merged = []
    current_start, current_end = coord_df.loc[0, "sstart"], coord_df.loc[0, "send"]

    for i in range(1, len(coord_df)):
        start, end = coord_df.loc[i, "sstart"], coord_df.loc[i, "send"]

        # Check if intervals overlap
        if start <= current_end + dist:
            # Merge intervals
            current_end = max(current_end, end)
        else:
            # Append the previous interval and start a new one
            merged.append([current_start, current_end])
            current_start, current_end = start, end

    # Append the last interval
    merged.append([current_start, current_end])

    # Convert to DataFrame
    return pd.DataFrame(merged, columns=["sstart", "send"])


def write_cleaned(o, contigs, fwidth=80):
    print("Writing results")
    with open(o, 'w') as out_file:
        for cname, sequence in contigs.items():
            out_file.write('>' + cname + '\n')
            for i in range(0, len(sequence), fwidth):
                out_file.write(sequence[i:i + fwidth] + "\n")


def single(args):
    if not args.cname:
        print("You need to specify a contig name for cleaning")
        sys.exit()

    if "fasta" in args.ofile or "fna" in args.ofile or ".scaffolds" in args.ofile:
        o = args.ofile
    else:
        o = args.ofile + ".fna"

    if args.startp < args.endp:
        s = args.startp
        e = args.endp
    else:
        s = args.endp
        e = args.startp

    print("Loading fasta")
    fasta_sequences = SeqIO.parse(open(args.ifasta), 'fasta')
    contigs = dict()

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        contigs[name] = sequence

    c1, c2 = break_contig(contigs, args.cname, s, e)
    contigs.pop(args.cname)
    if len(c1) > args.clen:
        contigs[args.cname + args.suffix + "1"] = c1
    if len(c2) > args.clen:
        contigs[args.cname + args.suffix + "2"] = c2

    write_cleaned(args.ofile, contigs, args.fwidth)


def multi(args, clean_auto=False):
    print("Loading fasta")
    fasta_sequences = SeqIO.parse(open(args.ifasta), 'fasta')
    contigs = dict()

    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        contigs[name] = sequence

    df = pd.read_csv(args.tfile, sep='\t', header=None)
    df.columns = ['query', 'ref', 'perc_id', 'qlen', 'slen', 'alen', 'sstart', 'send']
    df['cov'] = (df['alen'] / df['qlen']) * 100
    df = df[(df["perc_id"] >= args.idn) & (df["cov"] >= args.cov)]
    if len(df.index) == 0:
        print("No clean effect with actual cov and idn values")
        sys.exit()

    df = df.sort_values(by=['query', 'ref', 'sstart'], ascending=[True, True, True])
    plusdf = df[df['sstart'] < df['send']]
    minusdf = df[df['sstart'] > df['send']]
    tmp = minusdf['sstart']
    minusdf.loc[:, "sstart"] = minusdf["send"]
    minusdf.loc[:, "send"] = tmp
    df = pd.concat([plusdf, minusdf])

    for ref in df['ref'].unique():
        ref_coords = merge_intervals(df[df['ref'] == ref], dist=args.clen)
        for i, row in enumerate(ref_coords.itertuples(index=False)):
            c1, c2 = break_contig(contigs, ref, row.sstart, row.send)
            if len(c1) > args.clen:
                contigs[ref + args.suffix + str(i + 1)] = c1
            if len(c2) > args.clen:
                contigs[ref + args.suffix + str(i + 2)] = c2

        contigs.pop(ref)

    write_cleaned(args.ofile, contigs, args.fwidth)
    if clean_auto:
        os.remove(args.tfile)
        os.remove('log.txt')
        os.remove('err.txt')


def auto(args):
    db = args.blastdb
    clean_auto = False
    procoutput = open('log.txt', 'w', buffering=1)
    erroutput = open('err.txt', 'w', buffering=1)
    if not db:
        print("Indexing file for blast")
        db = 'bloutdb'
        mkbdb = ['makeblastdb', '-in', args.cfasta, 'dbtype', 'nucl',
                 '-out', db, '-max_file_sz', '4GB']

        subprocess.run(mkbdb, check=True, stdout=procoutput, stderr=erroutput)
        clean_auto = True

    # add output logs CHECK VIROSCOPIOOOO
    rtsv = ''.join(random.choices(string.ascii_letters + string.digits, k=length)) + '.tsv'
    bln = ['blastn', '-query', args.cfasta, '-out', rtsv, '-db', db,
           '-evalue', '0.001', '-perc_identity', str(args.idn),
           '-outfmt', '"6 qseqid sseqid pident qlen slen length sstart send"',
           '-num_threads', str(args.nproc)
           ]
    print("Performing blast")
    subprocess.run(bln, check=True, stdout=procoutput, stderr=erroutput)

    procoutput.close()
    erroutput.close()

    args.tfile = rtsv

    if clean_auto:
        for bidx in glob.glob("bloutdb.n*"):
            os.remove(bidx)

    multi(args, clean_auto=True)



if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='seqcleaner',
        description='seqcleaner will remove sub sequences from specified contig, breaking it into two (or more) subcontigs. in case you have a allvsall: tsv file: -outfmt "6 qseqid sseqid pident qlen slen length sstart send"',
        epilog='In case of error, please open an issue at github @sanrrone',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #parser.set_defaults(func=lambda _: parser.print_help())

    subparser = parser.add_subparsers(help='modules available to clean your contigs', required=True)

    parser_single = subparser.add_parser('single',
                                         help='single module allows you to clear one coordinate on a specific contig')
    parser_single.add_argument('-i', dest='ifasta',
                               help='input fasta with the "contaminated" sequences you want to clean or decontaminate',
                               required=True)  # positional argument
    parser_single.add_argument('-c', dest='cname',
                               help='contig name to clean or decontaminate. Ignored if -t is specified')
    parser_single.add_argument('-s', dest='startp', help='start position of the contig specified with -c (1-idx)', type=int,
                               required=True)
    parser_single.add_argument('-e', dest='endp', help='end position of the contig specified with -c (1-idx)', type=int,
                               required=True)
    parser_single.add_argument('-l', dest='clen',
                               help='length of the contig to be kept after the cleaning, 0 = keep all',
                               type=int, default=1000)
    parser_single.add_argument('-o', dest='ofile', help='name of the outputfile', type=str, default="out.fna")
    parser_single.add_argument('-sx', dest='suffix', help='suffix for cleaned contigs', type=str, default="__c")
    parser_single.add_argument('-w', dest='fwidth', help='width of sequence for the output fasta after cleaning',
                               type=int,
                               default=80)

    parser_multi = subparser.add_parser('multi',
                                        help='multi module allows you to clear multiple coordinate on multiple contigs, it needs a specific blast tsv: '
                                             '-outfmt "6 qseqid sseqid pident qlen slen length sstart send"')
    parser_multi.add_argument('-i', dest='ifasta',
                              help='input fasta with the "contaminated" sequences you want to clean or decontaminate',
                              required=True)  # positional argument
    parser_multi.add_argument('-t', dest='tfile',
                              help='tsv with specific formated blast output: -outfmt "6 qseqid sseqid pident qlen length sstart send"',
                              type=str, required=True)
    parser_multi.add_argument('-cov', dest='cov', help='coverage of the alignment in case of tsv file', type=float,
                              default=80)
    parser_multi.add_argument('-idn', dest='idn', help='identity of the alignment in case of tsv file', type=float,
                              default=80)
    parser_multi.add_argument('-l', dest='clen',
                              help='length of the contig to be kept after the cleaning, 0 = keep all',
                              type=int, default=1000)
    parser_multi.add_argument('-o', dest='ofile', help='name of the outputfile', type=str, default="out.fna")
    parser_multi.add_argument('-sx', dest='suffix', help='suffix for cleaned contigs', type=str, default="__c")
    parser_multi.add_argument('-w', dest='fwidth', help='width of sequence for the output fasta after cleaning',
                              type=int,
                              default=80)

    parser_auto = subparser.add_parser('auto',
                                       help='auto module performs the blast for you and clean the sequences based on specified identity/coverage')
    parser_auto.add_argument('-i', dest='ifasta',
                             help='fasta with the "contaminated" sequences you want to clean or decontaminate',
                             required=True)  # positional argument
    parser_auto.add_argument('-c', dest='cfasta', help='clean sequences file to map against the contaminated sequences',
                             type=str)
    parser_auto.add_argument('-db', dest='blastdb',
                             help='blast index in case you have them already. If not the software will do it',
                             type=str)
    parser_auto.add_argument('-cov', dest='cov', help='coverage of the alignment in case of tsv file', type=float,
                             default=80)
    parser_auto.add_argument('-idn', dest='idn', help='identity of the alignment in case of tsv file', type=float,
                             default=80)
    parser_auto.add_argument('-l', dest='clen', help='length of the contig to be kept after the cleaning, 0 = keep all',
                             type=int, default=1000)
    parser_auto.add_argument('-p', dest='nproc', help='number of cores to use for blast',
                             type=int, default=4)
    parser_auto.add_argument('-o', dest='ofile', help='name of the outputfile', type=str, default="out.fna")
    parser_auto.add_argument('-sx', dest='suffix', help='suffix for cleaned contigs', type=str, default="__c")
    parser_auto.add_argument('-w', dest='fwidth', help='width of sequence for the output fasta after cleaning',
                             type=int,
                             default=80)

    args = parser.parse_args()
    args.func(args)
