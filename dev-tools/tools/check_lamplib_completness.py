from ccdproc import ImageFileCollection


class CheckCompletness(object):
    keywords = ['wavmode', 'object', 'filter2']

    def __init__(self, path):
        self.path = path
        self.ic = None
        self.pd_ic = None

    def __call__(self, pattern, *args, **kwargs):
        self.ic = ImageFileCollection(self.path, keywords=self.keywords,
                                      glob_include=pattern)
        self.pd_ic = self.ic.summary.to_pandas()

        grouped = self.pd_ic.groupby(
            by=['wavmode']).size().reset_index().rename(columns={0: 'count'})

        for i in grouped.index:
            wavmode = grouped.iloc[i]['wavmode']
            lamps = self.pd_ic.object[self.pd_ic.wavmode == wavmode].tolist()
            print("Wavmode: {:s}, N: {:d}".format(wavmode, len(lamps)))
            for obj in sorted(lamps):
                print("\t{:s}".format(obj))
                # _object = grouped.iloc[i]['object']
                # self.pd_ic.object[self.pd_ic.wavmode == wavmode]
                #
                # print(i + 1, grouped.iloc[i]['wavmode'], grouped.iloc[i]['object'], grouped.iloc[i]['filter2'])

if __name__ == '__main__':
    check = CheckCompletness(path='/user/simon/data/soar/comp_lamp_lib/work/goodman_comp')
    print('ext_*fits')
    check(pattern='ext_*fits')
    # print('comb_*fits')
    # check(pattern='comb_*fits')
