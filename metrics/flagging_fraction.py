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
from astropy.table import Table
import sys

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


def get_flagging_frac_metric(json_path, threshold=70):
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
                                 'mean': statistics.fmean(flagged_data_core),
                                 'n_thresh': sum(i >= threshold for i in flagged_data_core),
                                 'n_pass': sum(i < threshold for i in flagged_data_core),
                                 'n_tot': len(flagged_data_core)}
    except:
        flagged_fracs['core'] = {'median': np.nan, 'mean': np.nan,
                                 'n_thresh': np.nan, 'n_tot': np.nan, 'n_pass': np.nan}
    try:
        flagged_fracs['remote'] = {'median': statistics.median(flagged_data_remote),
                                   'mean': statistics.fmean(flagged_data_remote),
                                   'n_thresh': sum(i >= threshold for i in flagged_data_remote),
                                   'n_pass': sum(i < threshold for i in flagged_data_remote),
                                   'n_tot': len(flagged_data_remote)
                                   }
    except:
        flagged_fracs['remote'] = {'median': np.nan, 'mean': np.nan,
                                   'n_thresh': np.nan, 'n_tot': np.nan, 'n_pass': np.nan}
    try:
        flagged_fracs['international'] = {'median': statistics.median(flagged_data_intl),
                                          'mean': statistics.fmean(flagged_data_intl),
                                          'n_thresh': sum(i >= threshold for i in flagged_data_intl),
                                          'n_pass': sum(i < threshold for i in flagged_data_intl),
                                          'n_tot': len(flagged_data_intl)
                                          }
    except:
        flagged_fracs['international'] = {'median': np.nan, 'mean': np.nan,
                                          'n_thresh': np.nan, 'n_tot': np.nan, 'n_pass': np.nan}

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


def examine_global_flagging_metrics(cal_path=cal_json_file_path, threshold=70,
                                    output='flagged_data_cals',
                                    n_pass_core=40, n_pass_remote=10,
                                    n_pass_intl=10):
    """
    Take a path to all the calibrators, get the file list (other function),
    iterate through to get flagging stats (other function),
    and then plot up the results
    Also return in arrays/lists so can work with further
    :param n_pass_intl: Number of international stations required to pass
    :param n_pass_remote: Number of remote stations required to pass
    :param n_pass_core: Number of core stations required to pass
    :param threshold: number, percentage data flagged threshold, default=70
    :param output: filename (no extension) for output image and table
    :param cal_path: path of all calibrator json files
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
    fd_c_thresh = []
    fd_c_pass = []
    fd_c_tot = []
    fd_r_mean = []
    fd_r_median = []
    fd_r_thresh = []
    fd_r_tot = []
    fd_r_pass = []
    fd_i_mean = []
    fd_i_median = []
    fd_i_thresh = []
    fd_i_tot = []
    fd_i_pass = []
    core_pass = []
    remote_pass = []
    intl_pass = []
    for f in json_file_list:
        ff = get_flagging_frac_metric(f, threshold=threshold)
        obs_list.append(f.split('/')[-2])
        cal_list.append(f.split('/')[-1].split('_')[0])
        # cal_dict = {'Observation ID': obs, 'Calibrator': cal,
        #             'Flagging fractions': ff}
        fd_c_mean.append(ff['core']['mean'])
        fd_c_median.append(ff['core']['median'])
        fd_c_thresh.append(ff['core']['n_thresh'])
        fd_c_tot.append(ff['core']['n_tot'])
        fd_c_pass.append(ff['core']['n_pass'])
        fd_r_mean.append(ff['remote']['mean'])
        fd_r_median.append(ff['remote']['median'])
        fd_r_thresh.append(ff['remote']['n_thresh'])
        fd_r_tot.append(ff['remote']['n_tot'])
        fd_r_pass.append(ff['remote']['n_pass'])
        fd_i_mean.append(ff['international']['mean'])
        fd_i_median.append(ff['international']['median'])
        fd_i_thresh.append(ff['international']['n_thresh'])
        fd_i_tot.append(ff['international']['n_tot'])
        fd_i_pass.append(ff['international']['n_pass'])
        # check and store if things pass
        core_pass.append(True) if ff['core']['n_pass'] >= n_pass_core else core_pass.append(False)
        remote_pass.append(True) if ff['remote']['n_pass'] >= n_pass_remote else remote_pass.append(False)
        intl_pass.append(True) if ff['international']['n_pass'] >= n_pass_intl else intl_pass.append(False)
    # calculate total numbers
    total_core_pass = np.sum(np.where(core_pass, 1, 0))
    total_remote_pass = np.sum(np.where(remote_pass, 1, 0))
    total_intl_pass = np.sum(np.where(intl_pass, 1, 0))
    total_obs = len(core_pass)
    # get an overall table to output
    metric_table = Table([obs_list, core_pass, remote_pass, intl_pass,
                          fd_c_pass, fd_r_pass, fd_i_pass],
                         names=('ObsID','Core Pass', 'Remote Pass', 'Intl Pass',
                                'Nstation good core', 'Nstation good remote',
                                'Nstation good intl'))
    metric_table.meta['comments'] = [f'Threshold of flagging is {threshold}%',
                                     f"{n_pass_core} core stations to pass",
                                     f"{n_pass_remote} remote stations to pass",
                                     f"{n_pass_intl} international stations to pass",
                                     f'For core stations, {total_core_pass}/{total_obs} observations pass',
                                     f'For remote stations, {total_remote_pass}/{total_obs} observations pass',
                                     f'For international stations, {total_intl_pass}/{total_obs} observations pass']
    metric_table.write(f'{output}.csv', format='csv', overwrite=True, comment='# ')

    # now visualize things
    # plot histograms for everything
    # worry about parsing later
    # okay do a quick check to ignore that are "POINT" - not sure they're good
    ind_cal = np.where(np.array(cal_list) != 'POINT')[0]
    fig, ((ax1, ax2, ax3),
          (ax4, ax5, ax6),
          (ax7, ax8, ax9)) = plt.subplots(3,3, figsize=(12, 12))
    ax1.hist([fd_c_mean, fd_c_median], 30,  histtype='step', fill=False, label=['Mean', 'Median'])
    ax1.hist([np.array(fd_c_mean)[ind_cal], np.array(fd_c_median)[ind_cal]], 30, histtype='step', fill=True,
             alpha=0.3, color=['#1f77b4', '#ff7f0e'])
    ax1.set_title('Core stations flagged data fractions')
    ax2.hist([fd_r_mean, fd_r_median], 30, histtype='step', fill=False, label=['Mean', 'Median'])
    ax2.hist([np.array(fd_r_mean)[ind_cal], np.array(fd_r_median)[ind_cal]], 30, histtype='step', fill=True,
             alpha=0.3, color=['#1f77b4', '#ff7f0e'])
    ax2.set_title('Remote stations flagged data fractions')
    ax3.hist([fd_i_mean, fd_i_median], 30, histtype='step', fill=False, label=['Mean', 'Median'])
    ax3.hist([np.array(fd_i_mean)[ind_cal], np.array(fd_i_median)[ind_cal]], 30, histtype='step', fill=True,
             alpha=0.3, color=['#1f77b4', '#ff7f0e'])
    ax3.set_title('Intl stations flagged data fractions')
    ax3.legend()
    ax4.scatter(fd_c_median, np.array(fd_c_thresh)/np.array(fd_c_tot), s=3)
    ax4.plot([threshold, threshold], [0,1], color='orange')
    ax4.plot([0, 100], [0.5, 0.5], color='orange', linestyle='--')
    ax4.set_ylabel('Fraction of stations above threshold')
    ax4.set_xlabel('Median flagged fraction')
    ax4.set_title(f'Core Threshold of {threshold}')
    ax5.scatter(fd_r_median, np.array(fd_r_thresh)/np.array(fd_r_tot), s=3)
    ax5.plot([threshold, threshold], [0, 1], color='orange')
    ax5.plot([0, 100], [0.5, 0.5], color='orange', linestyle='--')
    ax5.set_ylabel('Fraction of stations above threshold')
    ax5.set_xlabel('Median flagged fraction')
    ax5.set_title(f'Remote Threshold of {threshold}')
    ax6.scatter(fd_i_median,np.array(fd_i_thresh)/np.array(fd_i_tot), s=3)
    ax6.plot([threshold, threshold], [0, 1], color='orange')
    ax6.plot([0, 100], [0.5, 0.5], color='orange', linestyle='--')
    ax6.set_ylabel('Fraction of stations above threshold')
    ax6.set_xlabel('Median flagged fraction')
    ax6.set_title(f'Intl Threshold of {threshold}')
    # Look at numbers that pass
    ax7.scatter(fd_c_median, fd_c_pass, s=3)
    ax7.plot([threshold, threshold], [0, 50], color='orange')
    ax7.plot([0, 100], [n_pass_core, n_pass_core], color='orange', linestyle='--')
    ax7.text(10, 40, f'{total_core_pass}/{total_obs} observations pass')
    ax7.set_ylabel('Number of stations passing')
    ax7.set_xlabel('Median flagged fraction')
    ax8.scatter(fd_r_median, fd_r_pass, s=3)
    ax8.plot([threshold, threshold], [0, 15], color='orange')
    ax8.plot([0, 100], [n_pass_remote, n_pass_remote], color='orange', linestyle='--')
    ax8.text(10, 10, f'{total_remote_pass}/{total_obs} observations pass')
    ax8.set_ylabel('Number of stations passing')
    ax8.set_xlabel('Median flagged fraction')
    ax9.scatter(fd_i_median, fd_i_pass, s=3)
    ax9.plot([threshold, threshold], [0, 15], color='orange')
    ax9.plot([0, 100], [n_pass_intl, n_pass_intl], color='orange', linestyle='--')
    ax9.text(10, 10, f'{total_intl_pass}/{total_obs} observations pass')
    ax9.set_ylabel('Number of stations passing')
    ax9.set_xlabel('Median flagged fraction')
    plt.savefig(f'{output}.pdf')











