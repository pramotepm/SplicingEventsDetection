from itertools import tee, izip, izip_longest
from report import Stat
from os import makedirs
from os.path import isfile, join, exists
import util


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
    con1 = len(a) > 2 and len(b) > 2
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
    status1 = status2 = True
    sub1 = []
    sub2 = []
    ASTIs = {}
    # After length of s1 and s2 became 1, the remaining is onyl 'A'
    # That means s1 and s2 were read overall
    while len(seq1) != 1 or len(seq2) != 1:
        sub1, status1 = update(seq1, sub1, status1)
        sub2, status2 = update(seq2, sub2, status2)
        if status1 and status2:
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


def extract_event(stat, input_file_path, csv_file_path, stat_file_path):
    print input_file_path
    variant_comm_lines = {}
    ref_id = None
    trans_id = None
    for line in util.find_missing_exon(input_file_path):
        if line.startswith('VAR'):
            trans_id = line.split()[1]
            variant_comm_lines[trans_id] = []
        elif line.startswith('REF'):
            ref_id = line.split()[1].split('.')[0]
        elif line.startswith('COMM'):
            _ = line.split()
            _ = map(lambda x: int(x), (_[2:8]))
            variant_comm_lines[trans_id].append(_)

    with open(csv_file_path, 'a') as f_csv_out, open(stat_file_path + '.stat', 'w') as f_stat_out:
        for trans_id in variant_comm_lines:
            var_id = trans_id.split('.')[0]
            f_stat_out.write("%s <-> %s\n" % (ref_id, var_id))
            bt_r, bt_v = generate_bit_array(variant_comm_lines[trans_id])
            for bit_pattern, n_bit_pattern in extract_asti(bt_r, bt_v).items():
                stat.add(bit_pattern, n_bit_pattern)
                dec_pattern = util.reformat_to_dec(bit_pattern)
                f_csv_out.write("%s,%s,\"%s\",\"%s\",%d\n" % (ref_id, var_id, bit_pattern, dec_pattern, n_bit_pattern))
                f_stat_out.write("%s %s %d\n" % (bit_pattern, dec_pattern, n_bit_pattern))
            f_stat_out.write("//\n")
    return stat


def main(input_dir, output_dir):
    stat = Stat()

    if not exists(output_dir):
        makedirs(output_dir)

    result_dir = join(output_dir, 'out')
    if not exists(result_dir):
        makedirs(result_dir)

    stats_dir = join(output_dir, 'stats')
    if not exists(stats_dir):
        makedirs(stats_dir)

    # Create header/column name for an output file
    splicing_event_out_filename = 'splicing_events.csv'
    f_csv_out = join(output_dir, splicing_event_out_filename)
    f_out = open(f_csv_out, 'w')
    f_out.write('ref_id,var_id,binary,decimal,no_pattern\n')
    f_out.close()

    # test = 1010
    for f_path in util.read_dir(input_dir):
        f_path_in = join(input_dir, f_path)
        f_stat_out = join(stats_dir, f_path)

        stat = extract_event(stat, f_path_in, f_csv_out, f_stat_out)
        # if test == 0:
        #     break
        # test -= 1

    with open(join(result_dir, 'stat_report.log'), 'w') as f_out:
        f_out.write(stat.report())
    with open(join(result_dir, 'raw_data.txt'), 'w') as f_out:
        f_out.write(stat.report_raw_data())

if __name__ == '__main__':
    input_dir = './test'
    output_dir = './stat'

    # input_dir = "/Users/pramotepm/Desktop/mRNA/analyze/splicing_event/ascds/"
    # output_dir = "/Users/pramotepm/Desktop/mRNA/analyze/splicing_event/stat"
    main(input_dir, output_dir)