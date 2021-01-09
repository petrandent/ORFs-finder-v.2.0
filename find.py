import os
os.chdir("C:\\Users\\USER\\Desktop\\coursera")
print(os.getcwd())

start_codon=['ATG']
stop_codons=['TAA','TAG','TGA']


def read_fasta(filename='test.fasta'):
    """make a fasta file a dictionary of all fastas in the file"""
    file=open(filename)
    fast_dict={}
    for line in file:
        line=line.rstrip()
        if line[0]=='>':
            words=line.split()
            name=words[0][1:]
            fast_dict[name]=''
        else:
            fast_dict[name]=fast_dict[name]+line
    return fast_dict

def longest_shortest(filename='test.fasta'):
    """computes sequences with longest and shortest lengths and returns their length as a dictionary"""
    fast_dict=read_fasta(filename)
    length_checker={}
    lengths=[]
    final_result={}
    for i in fast_dict.items():
        length_checker[i[0]]=len(i[1])
    for i in length_checker.values():
        lengths.append(i)
    max_length,min_length=max(lengths),min(lengths)
    for i in length_checker.keys():
        if length_checker[i]==max_length:
            final_result[i]=length_checker[i]
        elif length_checker[i]==min_length:
            final_result[i]=length_checker[i]
        else: continue
    return(final_result)


def compl_dna(seq):
    """returns complementary DNA"""
    comp=""
    for i in seq:
        if i=="A":
            comp="T"+comp
        elif i=="T":
            comp="A"+comp
        elif i=="G":
            comp="C"+comp
        elif i=="C":
            comp="G"+comp
    return comp


def three_frames(seq):
    """three frames for a given seq"""
    seq1=seq[:]
    seq2=seq[1:]
    seq3=seq[2:]
    return [seq1,seq2,seq3]

def all_seqs(seq):
    """returns 2 variables one for template and one for complementary dna, each variable has a list of three sequences"""
    a1=three_frames(seq)
    a2=three_frames(compl_dna(seq))
    final=a1+a2
    return a1,a2


def dna_sparse(seq):
    """make a sequence a list of triplets"""
    triplets=[]
    for i in range(0,len(seq),3):
        triplets.append(seq[i:(i+3)])
    if len(triplets[-1])<3:
        del triplets[-1]
    return triplets

def first_orf(seq):
    """finds the first ORF in a sequence"""
    open=[]
    seq_sparse=dna_sparse(seq)
    clock=0
    while clock<len(seq_sparse):
        if seq_sparse[clock] in start_codon:
            open.append(clock)
            break
        else:
            clock+=1
    clock+=1
    while clock<len(seq_sparse):
        if seq_sparse[clock] in stop_codons:
            open.append(clock)
            break
        else:
            clock+=1
    if len(open)<2:
        del open
        return None
    else:
        return open



def orfs(seq):
    """returns reading frames in triplets position"""
    reading_frames=[]
    first=first_orf(seq)
    if first is None:
        return None
    else:
        reading_frames.append(first)
        seq=seq[(first[1]*3):]
        while True:
            first=first_orf(seq)
            if first is None:
                break
            else:
                reading_frames.append([i+reading_frames[-1][-1] for i in first])
                seq=seq[(first[1]*3):]
        return reading_frames

def lengths(seq):
    """computes the lenghts of reading frames in bases (NOT triplets)"""
    length_list=[]
    for i in orfs(seq):
        length_list.append(i[1]-i[0])
    return [i*3 for i in length_list]

def positions(seq):
    """returns the exact positions of ORFs, first position=0, last"""
    frames=orfs(seq)
    if frames is None:
        return None
    else:
        each=[]
        positions=[]
        for i in frames:
            start=int(i[0])*3
            end=int(i[1])*3+2
            each.append(start)
            each.append(end)
            positions.append(each)
            each=[]
        return positions



def orf_positions(seq):
    """computes all orfs in a sequence"""
    final=[]
    b1=three_frames(seq)
    for i in b1:
        final.append(positions(i))
    final1=[i for i in final if i!=None]
    final2=[j for i in final1 for j in i]
    return final2
