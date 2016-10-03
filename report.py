import operator
import util
import pprint

class Stat:
    def __init__(self):
        self.m = {}

    def add(self, pattern, n):
        if pattern not in self.m:
            self.m[pattern] = n
        else:
            self.m[pattern] += n

    def __create_format_str(self, l_as):
        if len(l_as) == 0:
            # default numbers / dummy numbers
            bt_arr_len = 3
            de_len = 3
            n_len = 3
        else:
            bt_arr_len = max(map(lambda x: len(x[0]), l_as))
            de_len = len(min(map(lambda x: util.reformat_to_dec(x[0]), l_as)))
            n_len = len(str(max([ i[1] for i in l_as ])))
        fmt = "{bt:%ds}\t{de:%ds}\t{n:>%dd}({p:.1%%})" % (bt_arr_len, de_len, n_len)
        return fmt


    def __total_cal(self, as_l):
        return float(reduce(lambda x, y: (x + y[1]) if type(x) is int else (x[1] + y[1]), as_l)) if len(as_l) >= 2 else \
        as_l[0][1]

    def __print_content(self, as_l, topN):
        total = 0
        remaining = 0
        s = []
        total = self.__total_cal(as_l)
        sorted_m = sorted(as_l, key=operator.itemgetter(1), reverse=True)
        if len(as_l) < topN:
            topN = len(as_l)

        fmt = self.__create_format_str(sorted_m[:topN])
        for _ in sorted_m[:topN]:
            __ = _[0].split(',')
            s.append(fmt.format(bt=__[0], de=util.reformat_to_dec(_[0]), n=_[1], p=_[1]/total))
            s.append(__[1])
        if len(sorted_m[topN:]) > 0:
            remaining = reduce(lambda x,y: (x if type(x) is int else x[1]) + y[1], sorted_m[topN:])
        s.append("others = %d(%.1f%%)" % (remaining, remaining/total*100))
        return s


    # Split into 2 types:
    #   - AS types proposed by Breitbart el at. (1987)
    #   - ASTI
    def _split_splicing_type(self):
        a = []
        b = []
        for item in self.m.items():
            _1, _2 = item[0].split(',')
            if _1.startswith('0') or _2.startswith('0'):
                a.append(item)
            else:
                b.append(item)
        return a, b

    def report(self, topN=10):
        ati_type, as_type = self._split_splicing_type()

        ati_str = ""
        as_str = ""
        if ati_type != []:
            ati_str = self.__print_content(ati_type, topN)
        if as_type != []:
            as_str = self.__print_content(as_type, topN)

        s = []
        s.append("# Overview Statistical Report for Each Representation")
        s.append("# Bit_Array\tDecimal_Format\tno.(percentage)")
        s.append("#")
        s.append("# ATI Type")
        s.extend(ati_str)
        s.append("#")
        s.append("# AS Type")
        s.extend(as_str)
        return '\n'.join(s) + '\n'

    def report_raw_data(self):
        ati_type, as_type = self._split_splicing_type()
        s = []
        total_as = self.__total_cal(as_type)
        total_ati = self.__total_cal(ati_type)
        s.append('# ATI splicing type, total=%d' % total_ati)
        s.extend(map(lambda x: "%s\t%s\t%d" % (x[0], util.reformat_to_dec(x[0]), x[1]), sorted(ati_type, key=lambda x: x[0])))
        s.append('#')
        s.append('# AS splicing type, total=%d' % total_as)
        s.extend(map(lambda x: "%s\t%s\t%d" % (x[0], util.reformat_to_dec(x[0]), x[1]), sorted(as_type, key=lambda x: x[0])))
        return '\n'.join(s) + '\n'
