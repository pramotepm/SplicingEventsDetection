from os import listdir
from os.path import isfile, join
from itertools import chain
import pprint
import operator

def reformat_to_dec(a):
    return '(' + ','.join(map(lambda x: str(int(x, 2)), a.split(','))) + ')'


def read_file(file_path):
    with open(file_path) as f_in:
        for _ in f_in:
            yield _.strip()


def read_dir(dir_path):
    return [f for f in listdir(dir_path) if isfile(join(dir_path, f))]


def chk_not_overlap(s1, s2):
    c = set()
    l = list(chain(*s2))
    for _ in s1:
        a, b = _
        if a not in l and b not in l:
            c.add(_)
    return s1 - c

def adjust_boundary():
    var_id = None
    ref_bound = set()
    ref_swap = set()
    var_bound = set()
    comm_bound = set()
    is_ref = True
    chk = False
    for line in read_file('test/NM_002654.4.as'):
        if line.startswith("VAR"):
            var_id = line.split()[1]
            print var_id
            is_ref = False
            chk = True
        elif line.startswith('EXON'):
            _, __ = line.split()[4:6]
            if is_ref:
                ref_bound.add((_, __))
            else:
                var_bound.add((_, __))
            ref_swap = ref_bound.copy()
        elif line.startswith('COMM'):
            _, __ = line.split()[2:4]
            comm_bound.add(_)
            comm_bound.add(__)
        elif line.startswith('//') and chk:
            comm_bound -= set(chain(*ref_bound))
            comm_bound -= set(chain(*var_bound))

            ref_bound, var_bound = ref_bound - var_bound, var_bound - ref_bound

            ref_bound = chk_not_overlap(ref_bound, var_bound)
            var_bound = chk_not_overlap(var_bound, ref_bound)

            remaining = sorted([ (a,b,int(a)-int(b)+1) for a, b in ref_bound | var_bound ], key=operator.itemgetter(2))
            print remaining
            print comm_bound
            while len(remaining) > 0:
                st, ed, _ = remaining.pop()
                tmp_st = filter(lambda x: st in x[0], remaining)
                if tmp_st != []:
                    remaining.remove(tmp_st[0])
                    missing_bound = max(int(ed), int(tmp_st[0][1])) - 1
                    comm_bound.remove(str(missing_bound))
                tmp_ed = filter(lambda x: ed in x[1], remaining)
                if tmp_ed != []:
                    remaining.remove(tmp_ed[0])
                    missing_bound = min(int(st), int(tmp_ed[0][0])) + 1
                    comm_bound.remove(str(missing_bound))
            print comm_bound
            if len(comm_bound) != 0:
                print ">>>>" + var_id
            ref_bound = ref_swap.copy()
            var_bound.clear()
            comm_bound.clear()
            # break


if __name__ == '__main__':
    adjust_boundary()