from os import listdir
from os.path import isfile, join
from itertools import chain
import pprint


class Direction:
    # Ascending
    asc = 0
    # Descending
    des = 1


def reformat_to_dec(a):
    return '(' + ','.join(map(lambda x: str(int(x, 2)), a.split(','))) + ')'


def read_file(file_path):
    with open(file_path) as f_in:
        for _ in f_in:
            yield _.strip()


def read_dir(dir_path):
    return [f for f in listdir(dir_path) if isfile(join(dir_path, f))]


def __chk_stat_in_exon(cur_pos, bounds, direction):
    # start  ;return 1 or -1 if ascending
    # middle ;return 0
    # end    ;return -1 or 1 if ascending
    # extragenic region or intronic region ;return 2
    # current position that matched start and end, means its length is 1 (len = 1) ;return 3
    if direction == Direction.asc:
        compare_f = lambda x: x[0] <= cur_pos <= x[1]
        shift_before_start = -1
        shift_after_end = 1
    elif direction == Direction.des:
        compare_f = lambda x: x[0] >= cur_pos >= x[1]
        shift_before_start = 1
        shift_after_end = -1

    tmp = filter(compare_f, bounds)
    if tmp != []:
        start, end = tmp[0]
        if start == cur_pos == end:
            return 3
        #
        # *****|----| <- compare with this, and the current position of nucleotide matched with start position of this.
        # |---------|
        #
        if start == cur_pos:
            return shift_before_start
        #
        # |----|***** <- compare with this, and so on.
        # |---------|
        #
        elif end == cur_pos:
            return shift_after_end
        else:
            return 0
    return 2


def find_missing_exon(file_path):
    complete_lines = []
    org_lines = []
    ref_bound = []
    var_bound = []
    comm_bound = set()
    is_ref = True
    chk = False
    # var_id = None
    direction = None
    for line in read_file(file_path):
        if line.startswith('REF'):
            complete_lines.append(line)
        if line.startswith('VAR'):
            complete_lines.append(line)
            org_lines = []
            var_bound = []
            comm_bound.clear()
            # var_id = line.split()[1]
            # print var_id
            is_ref = False
            chk = True
            for i in range(0, len(ref_bound)):
                tmp = int(ref_bound[i][0]) - int(ref_bound[i][1])
                if tmp != 0:
                    direction = Direction.asc if tmp < 0 else Direction.des
                    break
        elif line.startswith('EXON'):
            _, __ = line.split()[4:6]
            if is_ref:
                ref_bound.append((int(_), int(__)))
            else:
                var_bound.append((int(_), int(__)))
        elif line.startswith('COMM'):
            org_lines.append(line)
            _, __ = line.split()[2:4]
            comm_bound.add((int(_), int(__)))
        elif line.startswith('//') and chk:
            c = set()
            all_bound = sorted(set(chain(*ref_bound)) | set(chain(*var_bound)), reverse=True if direction == Direction.des else False)
            # print all_bound
            start = None
            for bound in all_bound:
                # print bound
                shift_1 = __chk_stat_in_exon(bound, ref_bound, direction)
                shift_2 = __chk_stat_in_exon(bound, var_bound, direction)
                shift = shift_1 + shift_2
                if shift_1 == 0 or shift_2 == 0:
                    # shift_1 + shift_2 = 1 or -1 only
                    if (direction == Direction.asc and shift > 0) or (direction == Direction.des and shift < 0):
                        c.add((start, bound))
                        # print start, bound
                        start = bound + shift
                        # print 'start=', start
                    else:
                        c.add((start, bound + shift))
                        # print start, bound + shift
                        start = bound
                elif shift_1 == 3 or shift_2 == 3:
                        c.add((bound, bound))
                        # print bound, bound
                        # print start
                else:
                    if start is not None:
                        if (shift_1 == 1 and shift_2 == -1) or (shift_1 == -1 and shift_2 == 1):
                            if direction == Direction.asc:
                                c.add((start, bound - 1))
                                # print start, bound - 1
                                start = bound + 1
                                # print 'start=', start
                            elif direction == Direction.des:
                                c.add((start, bound + 1))
                                # print start, bound + 1
                                start = bound - 1
                                # print 'start=', start
                        else:
                            c.add((start, bound))
                            # print start, bound
                            start = None
                            # print 'start=', start
                    else:
                        start = bound
                        # print 'start=', start
            # print "---"
            if len(c) != len(comm_bound):
                # print ">>", var_id
                # print c - comm_bound
                # pprint.pprint(comm_bound)
                for missing_region in c - comm_bound:
                    if direction == Direction.asc:
                        in_ref = filter(lambda x: x[0] <= missing_region[0] and missing_region[1] <= x[1], ref_bound) != []
                        in_var = filter(lambda x: x[0] <= missing_region[0] and missing_region[1] <= x[1], var_bound) != []
                    elif direction == Direction.des:
                        in_ref = filter(lambda x: x[0] >= missing_region[0] and missing_region[1] >= x[1], ref_bound) != []
                        in_var = filter(lambda x: x[0] >= missing_region[0] and missing_region[1] >= x[1], var_bound) != []
                    fmt = "COMM  xx  {:d}  {:d}  {:>5d}  {:>5d}   {:>5d}  {:>5d}  {:>2d}  {:>2d}"
                    ref_p = -999 if in_ref else 0
                    var_p = -999 if in_var else 0
                    add_up = fmt.format(missing_region[0], missing_region[1], ref_p, ref_p, var_p, var_p, 1 if in_ref else -1, 1 if in_var else -1)
                    org_lines.append(add_up)
                complete_lines.extend(sorted(org_lines, key=lambda x: int(x.split()[2]), reverse=True if direction == Direction.des else False))
            else:
                complete_lines.extend(org_lines)
            complete_lines.append('//')
    return complete_lines

if __name__ == '__main__':
    # adjust_boundary()
    # print '\n'.join(find_missing_exon('/Users/pramotepm/PycharmProjects/SplicingEventsDetection/test/NM_021170.3.as'))
    print '\n'.join(find_missing_exon('test/NM_032329.4.as'))