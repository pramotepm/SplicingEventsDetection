import operator
import util


class Stat:
    def __init__(self):
        self.m = {}
        self.fmt = None

    def add(self, pattern, n):
        if pattern not in self.m:
            self.m[pattern] = n
        else:
            self.m[pattern] += n

    def __create_format_str(self):
        bt_arr_len = max(map(lambda x: len(x.split(',')[0]), self.m.keys()))
        de_len = len(min(map(lambda x: util.reformat_to_dec(x), self.m.keys())))
        n_len = len(str(max(self.m.values())))
        fmt = "{bt:%ds}\t{de:%ds}\t{n:>%dd}({p:.1%%})" % (bt_arr_len, de_len, n_len)
        return fmt

    def __print_content(self, l_type, topN):
        remaining = 0
        s = []

        total = float(reduce(lambda x,y: x + y[1] if type(x) is int else x[1]+y[1], l_type))
        sorted_m = sorted(l_type, key=operator.itemgetter(1), reverse=True)
        if len(l_type) < topN:
            topN = len(l_type)

        for _ in sorted_m[:topN]:
            __ = _[0].split(',')
            s.append(self.fmt.format(bt=__[0], de=util.reformat_to_dec(_[0]), n=_[1], p=_[1]/total))
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
        self.fmt = self.__create_format_str()

        print len(as_type), len(ati_type)
        ati_str = self.__print_content(ati_type, topN)
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
