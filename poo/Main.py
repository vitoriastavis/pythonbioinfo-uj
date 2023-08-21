import Seq

print(Seq.lines('ATG'+'GAC'*90+'TAA'))

seqs = Seq.from_file('seqs.fasta')