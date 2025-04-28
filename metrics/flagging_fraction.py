# Flagging fractions metric

__author__ = "E. A. K. Adams"

"""
Examine flagging fraction as a possible metric

Start with calibrators
"""

# import modules
import json
import statistics
import os
import glob
import numpy as np
import matplotlib.pyplot as plt

# define some global paths
# do this in a hacky way
this_dir, this_filename = os.path.split(__file__)
hbanldir = this_dir[:-7]
datadir = os.path.join(hbanldir, "data")
cal_json_file_path = os.path.join(datadir, 'cal_json')  # hope this relative path works


def get_cal_from_json(json_path):
    """Get calibrator information from json
    Don't just return the dictionary but drill down to a semi-useful level

    Returns field_name (string) and station_data (dictionary)
    """
    with open(json_path) as json_file:
        json_data = json.load(json_file)

    # drill down to a semi-useful level
    cal_data = json_data['metrics']['LINC']
    # close_sources = cal_data['close_sources']  # will maybe want this later
    field_name = cal_data['field_name']
    station_data = cal_data['stations']

    return field_name, station_data


def get_flagging_single_cal(json_path):
    """
    Load data from json, using helper function
    Get flagging fractions, this is for single calibrator
    :param json_path: Path to calibrator json file, string
    :return:
    """
    field_name, station_data = get_cal_from_json(json_path)
    # iterate through stations
    # want to store info by station_type
    flagged_data_core = []
    flagged_data_remote = []
    flagged_data_intl = []
    for s in station_data:
        # check station type
        station = s['station']
        if station[0:2] == 'CS':
            station_type = 'core'
        elif station[0:2] =='RS':
            station_type = 'remote'
        else:
            station_type = 'intl'
        # get flagged data
        flagged_frac = s['percentage_flagged']['final']
        # add to correct array
        if station_type == 'core':
            flagged_data_core.append(flagged_frac)
        elif station_type == 'remote':
            flagged_data_remote.append(flagged_frac)
        else:
            flagged_data_intl.append(flagged_frac)
    # return arrays of flagged data fractions
    # then can work with them elsewhere to explore metrics
    return flagged_data_core, flagged_data_remote, flagged_data_intl


def get_flagging_frac_metric(json_path):
    """For a single calibrator observation
    with output json given by json_path, get flagging metrics
    Put this in separate function so can easily update
    """
    # first get arrays of flagged fractions
    (flagged_data_core, flagged_data_remote,
     flagged_data_intl) = get_flagging_single_cal(json_path)

    # now get some metrics
    # place in a dictionary
    flagged_fracs = {}
    try:
        flagged_fracs['core'] = {'median': statistics.median(flagged_data_core),
                                 'mean': statistics.fmean(flagged_data_core)}
    except:
        flagged_fracs['core'] = {'median': np.nan, 'mean': np.nan}
    try:
        flagged_fracs['remote'] = {'median': statistics.median(flagged_data_remote),
                                   'mean': statistics.fmean(flagged_data_remote)}
    except:
        flagged_fracs['remote'] = {'median': np.nan, 'mean': np.nan}
    try:
        flagged_fracs['international'] = {'median': statistics.median(flagged_data_intl),
                                          'mean': statistics.fmean(flagged_data_intl)}
    except:
        flagged_fracs['international'] = {'median': np.nan, 'mean': np.nan}

    # and return the dictionary
    return flagged_fracs


def find_all_cals(cal_path=cal_json_file_path):
    """
    Search for the cal directory structure and create a list of all file paths.
    This will be used later for pulling statistics of everything to investigate metrics
    Also print an overall summary of how many cal files are available
    and how many for each calibrator
    :param cal_path:
    :return:
    """
    json_file_list = glob.glob(os.path.join(cal_path, 'L*/*calibrator_summary.json'))
    file_name_list = [f.split('/')[-1] for f in json_file_list]
    cal_name_list = [n.split('_')[0] for n in file_name_list]
    names, count = np.unique(cal_name_list, return_counts=True)
    print(f"The calibrators and their count are:")
    print([f"{c}: {n}" for c, n in zip(names, count)])
    print(f"The total number of calibrator json files is {len(json_file_list)}")
    # return the json_file_list with full paths
    return json_file_list


def examine_global_flagging_metrics(cal_path=cal_json_file_path):
    """
    Take a path to all the calibrators, get the file list (other function),
    iterate through to get flagging stats (other function),
    and then plot up the results
    Also return in arrays/lists so can work with further
    :param cal_path:
    :return:
    """
    json_file_list = find_all_cals(cal_path)
    # create a dictionary to put things in
    # want to keep structure and also ability to sort by type
    # can't actually figure out how to do this nicely, so just create massive arrays/lists
    cal_list = []
    obs_list = []
    fd_c_mean = []
    fd_c_median = []
    fd_r_mean = []
    fd_r_median = []
    fd_i_mean = []
    fd_i_median = []
    for f in json_file_list:
        ff = get_flagging_frac_metric(f)
        obs_list.append(f.split('/')[-2])
        cal_list.append(f.split('/')[-1].split('_')[0])
        # cal_dict = {'Observation ID': obs, 'Calibrator': cal,
        #             'Flagging fractions': ff}
        fd_c_mean.append(ff['core']['mean'])
        fd_c_median.append(ff['core']['median'])
        fd_r_mean.append(ff['remote']['mean'])
        fd_r_median.append(ff['remote']['median'])
        fd_i_mean.append(ff['international']['mean'])
        fd_i_median.append(ff['international']['median'])

    # now visualize things
    # plot histograms for everything
    # worry about parsing later
    fig, (ax1, ax2, ax3) = plt.subplots(1,3, figsize=(12, 4))
    ax1.hist([fd_c_mean, fd_c_median], 50,  histtype='step', fill=False, label=['Mean', 'Median'])
    ax1.set_title('Core stations flagged data fractions')
    ax2.hist([fd_r_mean, fd_r_median], 50, histtype='step', fill=False, label=['Mean', 'Median'])
    ax2.set_title('Remote stations flagged data fractions')
    ax3.hist([fd_i_mean, fd_i_median], 50, histtype='step', fill=False, label=['Mean', 'Median'])
    ax3.set_title('Intl stations flagged data fractions')
    ax3.legend()
    plt.savefig('fd_hist.pdf')











