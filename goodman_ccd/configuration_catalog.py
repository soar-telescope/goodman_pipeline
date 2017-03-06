from __future__ import print_function
from ccdproc import ImageFileCollection
import glob

path = '/user/simon/data/soar/raw/spectroscopy_engineering_night/'

eng_folders = glob.glob(path + '*')
all_gratings = []
all_configurations = []
for dfolder in sorted(eng_folders):
    print(dfolder)
    keywords_red = ['date',
                    'slit',
                    'date-obs',
                    'obstype',
                    'object',
                    'exptime',
                    'obsra',
                    'obsdec',
                    'grating',
                    'wavmode',
                    'cam_targ',
                    'grt_targ',
                    'filter',
                    'filter2']

    keywords_blue = ['date',
                     'slit',
                     'date-obs',
                     'obstype',
                     'object',
                     'exptime',
                     'obsra',
                     'obsdec',
                     'grating',
                     'cam_targ',
                     'grt_targ',
                     'filter',
                     'filter2']
    try:
        ic = ImageFileCollection(dfolder, keywords_red)
        pandas_df = ic.summary.to_pandas()
        f_pandas_df = pandas_df[pandas_df.obstype == 'OBJECT']
        # print(f_pandas_df.date)
        for i in f_pandas_df.index:
            # print(i)
            if str(f_pandas_df.wavmode[i]) != 'Imaging':
                try:
                    wavmode = str(f_pandas_df.wavmode[i]) + ' ' + str(f_pandas_df['cam_targ'][i]) + ' ' \
                                  + str(f_pandas_df['grt_targ'][i])
                    if 'nan' in wavmode:
                        wavmode = str(f_pandas_df.grating[i]).split('_')[1] + ' ' \
                                  + str(f_pandas_df['cam_targ'][i]) + ' ' \
                                  + str(f_pandas_df['grt_targ'][i])
                    elif 'Custom' in wavmode:
                        wavmode = str(f_pandas_df.grating[i]).split('_')[1] + '_' \
                                  +'Custom' + ' ' \
                                  + str(f_pandas_df['cam_targ'][i]) + ' ' \
                                  + str(f_pandas_df['grt_targ'][i])
                    filter2 = str(f_pandas_df['filter2'][i])
                    if filter2 == '<NO FILTER>':
                        filter2 = ''
                    class_string = wavmode + ' ' + filter2
                    if class_string not in all_configurations:
                        all_configurations.append(class_string)

                except IndexError:
                    print('IndexError')

    except KeyError:
        ic = ImageFileCollection(dfolder, keywords_blue)
        pandas_df = ic.summary.to_pandas()
        f_pandas_df = pandas_df[pandas_df.obstype == 'OBJECT']
        # print(f_pandas_df.date)
        for i in f_pandas_df.index:
            if str(f_pandas_df['filter2'][i]) != '<NO FILTER>':
                filter2 = str(f_pandas_df['filter2'][i])
            else:
                filter2 = ''

            class_string = + str(f_pandas_df.grating[i]) + ' ' \
                           + filter2 + ' ' \
                           + str(f_pandas_df.cam_targ[i]) + ' ' \
                           + str(f_pandas_df.grt_targ[i]) + ' ' \

                # + str(f_pandas_df['filter'][i])
            if class_string not in all_configurations:
                all_configurations.append(class_string)




    for grating in f_pandas_df.grating.unique():
        if grating != '<NO GRATING>' and grating not in all_gratings:
            all_gratings.append(grating)

# if len(all_gratings) > 0:
#     for i in sorted(all_gratings):
#         print(i)

print(" ")
for conf in sorted(all_configurations):
    print(conf)


    # night_object = Night()
    # bias_list = pandas_df[pandas_df.obstype == 'BIAS']
    # night_object.add_bias(bias_list)
    # if 'imaging' in path:
    #     all_flats_list = pandas_df[pandas_df.obstype == 'FLAT']
    #     # print(all_flats_list['filter'][pandas_df.obstype == 'FLAT'].unique())
    #     for i_filter in pandas_df['filter'][pandas_df.obstype == 'FLAT'].unique():
    #         flat_group = all_flats_list[all_flats_list['filter'] == i_filter]
    #         night_object.add_day_flats(flat_group)
    #         print(flat_group)
    #         # print(all_flats_list[all_flats_list.filter == i_filter])
