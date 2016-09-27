from itertools import tee, izip, izip_longest


def read_file(f):
    with open(f) as f_in:
        for _ in f_in:
            yield _.strip()


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)


def assess_terminal(s):
    pos = s.rfind('1')
    return s[:pos+1] + 'A' # + s[pos+2:]


def update(bit_arr, ASTI, stat):
    if len(bit_arr) > 1:
        e = bit_arr.pop()
        ASTI.append(e)
        stat = False if e == '0' else True
    else:
        ASTI.append('0')
    return ASTI, stat


def chk_asti(sub1, sub2):
    a = ''.join(sub1)
    b = ''.join(sub2)
    con1 = len(a) >  2 and len(b) > 2
    con2 = a.count('1') > 1 and b.count('1') > 1
    if con1 and con2:
        a, b = (b, a) if b < a else (a, b)
        return "%s,%s" % (a, b)
    else:
        return None


# ASTI = Alternative Splicing and Transcriptional Initiation
def extract_asti(seq1, seq2):
    seq1 = list(seq1)
    seq2 = list(seq2)
    seq1.reverse()
    seq2.reverse()
    stat1 = stat2 = True
    sub1 = []
    sub2 = []
    ASTIs = {}
    # After length of s1 and s2 became 1, the remaining is onyl 'A'
    # That means s1 and s2 were read overall
    while len(seq1) != 1 or len(seq2) != 1:
        sub1, stat1 = update(seq1, sub1, stat1)
        sub2, stat2 = update(seq2, sub2, stat2)
        if stat1 and stat2:
            if sub1 != sub2:
                pattern = chk_asti(sub1, sub2)
                if pattern is not None:
                    if pattern not in ASTIs:
                        ASTIs[pattern] = 1
                    else:
                        ASTIs[pattern] += 1
            sub1 = list(sub1[-1])
            sub2 = list(sub2[-1])
    return ASTIs


def generate_bit_array(V):
    #
    # variable exon_* explains whether there is exon in each region or not (0, 1)
    # variable intron generates the introns as 0 and '-' when two adjacency bounds was found.
    #
    exon_ref = map(lambda x: '1' if x[2] != 0 and x[3] != 0 else '0', V)
    exon_var = map(lambda x: '1' if x[4] != 0 and x[5] != 0 else '0', V)
    intron =  map(lambda rec: '' if abs(rec[0][1] - rec[1][0]) == 1 else '0', pairwise(V))
    bit_arr_ref = ''.join([a + b for (a, b) in izip_longest(exon_ref, intron, fillvalue='0')])
    bit_arr_var = ''.join([a + b for (a, b) in izip_longest(exon_var, intron, fillvalue='0')])
    return assess_terminal(bit_arr_ref), assess_terminal(bit_arr_var)


def main():
    V = {}
    trans_id = None
    for line in read_file('NM_002654.4.as'):
        if line.startswith('VAR'):
            trans_id = line.split()[1]
            V[trans_id] = []
        elif line.startswith('COMM'):
            _ = line.split()
            _ = map(lambda x : int(x), (_[2:8]))
            V[trans_id].append(_)

    for trans_id in V:
        print trans_id
        print "---"
        bt_r, bt_v = generate_bit_array(V[trans_id])
        print bt_r
        print bt_v
        print "---"
        for p in extract_asti(bt_r, bt_v):
            print p
        print ""
        # break

if __name__ == '__main__':
    main()