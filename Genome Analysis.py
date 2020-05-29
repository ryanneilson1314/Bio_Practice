from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio import SeqIO

# Importing Gene using FASTA
for seq_record in SeqIO.parse("Mouse_Hemoglobin_Beta.fasta", "fasta"):
    # print(seq_record.id)
    Mouse_Hem_Gene = (str(seq_record.seq))
    # print(len(seq_record))
    # print(Mouse_Hem_Gene)

for seq_record in SeqIO.parse("Mouse_Myoglobin.fasta", "fasta"):
    # print(seq_record.id)
    Mouse_Myo_Gene = (str(seq_record.seq))
    # print(len(seq_record))
    # print(Mouse_Myo_Gene)

Sequence = ("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")  # For TroubleShooting


# print(len(Sequence))


# for DNA transcripton and translation
def DNA_to_mRNA(DNA_Seq):
    coding_dna = Seq(DNA_Seq, IUPAC.unambiguous_dna)
    template_dna = coding_dna.reverse_complement()
    mRNA = coding_dna.transcribe()
    template_dna.reverse_complement().transcribe()
    return mRNA


mRNA_Mouse_Hem_Gene = DNA_to_mRNA(Mouse_Hem_Gene)
mRNA_Mouse_Myo_Gene = DNA_to_mRNA(Mouse_Myo_Gene)
mRNA_Sequence = DNA_to_mRNA(Sequence)


# print(mRNA_Mouse_Hem_Gene)
# print(mRNA_Mouse_Myo_Gene)
# print(mRNA_Sequence)


def mRNA_Translation(mRNA):  # TODO: figure out how to only return AA of mRNA string
    if len(mRNA) % 3 == 0:
        return "The amino acid sequence is %s, * represents a stop codon." % mRNA.translate()
    else:
        return "Partial codons in mRNA"


AA_Mouse_Hem = mRNA_Translation(mRNA_Mouse_Hem_Gene)
AA_Mouse_Myo = mRNA_Translation(mRNA_Mouse_Myo_Gene)
AA_Sequence = mRNA_Translation(mRNA_Sequence)


print(AA_Mouse_Hem)
print(AA_Mouse_Myo)
print(AA_Sequence)


def DNA_Nuc_Content(DNA_Seq):
    count_A = []
    count_T = []
    count_C = []
    count_G = []
    for n in DNA_Seq:
        if n == "A":
            count_A.append(n)
        elif n == "T":
            count_T.append(n)
        elif n == "C":
            count_C.append(n)
        elif n == "G":
            count_G.append(n)
    Percent_of_A_in_DNA = (len(count_A) / len(DNA_Seq)) * 100
    Percent_of_T_in_DNA = (len(count_T) / len(DNA_Seq)) * 100
    Percent_of_C_in_DNA = (len(count_C) / len(DNA_Seq)) * 100
    Percent_of_G_in_DNA = (len(count_G) / len(DNA_Seq)) * 100
    return "There are %s adenines, %s thymines, %s cytosines, and %s guanines. " \
           "The percent of adenines is %s percent, thymines is %s percent, cytosines is %s percent, " \
           "and guanines is %s percent." \
           % (len(count_A), len(count_T), len(count_C), len(count_G), round(Percent_of_A_in_DNA, 2),
              round(Percent_of_T_in_DNA, 2), round(Percent_of_C_in_DNA, 2), round(Percent_of_G_in_DNA, 2))

# print(DNA_Nuc_Content(Sequence))
# print(DNA_Nuc_Content(Mouse_Hem_Gene))
# print(DNA_Nuc_Content(Mouse_Myo_Gene))
