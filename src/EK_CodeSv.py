"""
EK_CodeSv.py

select data in regions or between lines and apply different stats to the data

jech
"""

import echopype as ep
from echopype import open_raw
import echoregions as er
from pathlib import Path
from sys import exit
#from dask.distributed import Client
import xarray as xr
from matplotlib.pyplot import figure, show, subplots_adjust, get_cmap, cm
#import matplotlib.ticker as mtick
import matplotlib.dates as mdates
import hvplot.xarray
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
import re
from datetime import datetime, time, timedelta
import itertools

# the Client gives lots of error messages, so I don't use it for now
#client = Client()

def parse_fnames(args):
    '''
    parse the arguments and create list files with the full path and file names
    for the data and line files

    input: command line arguments
    output: 
            list of data file names
            list of line file names
            dictionary with the json parameter settings
    '''
    EVL_dict = {}
    EKfn_dict = {}

    ### data files
    # if a data directory is specified and files are specified, create a list of
    # data file names
    if (args.datadirectory and args.filenames):
        ct = 0
        for f in args.filenames:
            EKfn_dict['f'+str(ct)] = {
                    'EK_data_path' : Path(args.datadirectory), 
                    'EK_data_filenames' : Path(f),
                    'save_path' : Path(args.datadirectory+'generated')
                    }
            #datafilelist.append(args.datadirectory / f)
            print('Data File: ', EKfn_dict['f'+str(ct)]['EK_data_path'] / 
                    EKfn_dict['f'+str(ct)]['EK_data_filenames'])
            ct += 1
        tmpdict.setdefault(linefn, {}).update({np.datetime64(tmp):round(float(dd), 2)})
    # if data directory is specified but no files, glob the entire directory
    if (args.datadirectory and not args.filenames):
        ct = 0
        for f in sorted(args.datadirectory.glob('**/*.nc')):
            EKfn_dict['f'+str(ct)] = {
                    'EK_data_path' : Path(args.datadirectory), 
                    'EK_data_filenames' : Path(f),
                    'save_path' : Path(args.datadirectory+'generated')
                    }
            print('Glob Data File: ', EKfn_dict['f'+str(ct)]['EK_data_path'] / 
                    EKfn_dict['f'+str(ct)]['EK_data_filenames'])
            ct += 1
            #print('Glob .nc file: ', f)
            #datafilelist.append(f)
    
    ### line files
    # if a line directory is specified and files are specified, create a list of
    # line file names
    if (args.linedirectory and args.filenames):
        ct = 0
        for f in args.filenames:
            EVL_dict['l'+str(ct)] = {
                    'path' : Path(args.datadirectory), 
                    'filenames' : Path(f),
                    'linecolor' : 'black',
                    'linewidth' : 1
                    }
            print('EVL File: ', EVL_dict['f'+str(ct)]['path'] / 
                    EVL_dict['f'+str(ct)]['filenames'])
            ct += 1
            #linefilelist.append(args.linedirectory / f)
            #print('Line File: ', args.linedirectory / f)
    # if line directory is specified but no files, glob the entire directory
    if (args.linedirectory and not args.filenames):
        for f in sorted(args.linedirectory.glob('**/*.evl')):
            EVL_dict['l'+str(ct)] = {
                    'path' : Path(args.datadirectory), 
                    'filenames' : Path(f),
                    'linecolor' : 'black',
                    'linewidth' : 1
                    }
            print('Glob EVL File: ', EVL_dict['f'+str(ct)]['path'] / 
                    EVL_dict['f'+str(ct)]['filenames'])
            ct += 1
            #print('Glob .evl file: ', f)
            #linefilelist.append(f)
    
    ### json
    # if a json file is specified, get the file names and paths
    # the json file overwrites any other specified files
    if (args.jsonfile):
        print('json provided')
        with open(args.jsonfile) as jsonfile:
            pars_dict = json.load(jsonfile)

    # data files
    # if data path and data filenames are specified, create a nested dictionary of
    # data file names
    if (pars_dict['path_config']):
            ct = 0
            for f in pars_dict['path_config']['EK_data_filenames']:
                EKfn_dict['f'+str(ct)] = {
                    'EK_data_path' : Path(pars_dict['path_config']['EK_data_path']),
                    'EK_data_filenames' : Path(f),
                    'save_path' : Path(pars_dict['path_config']['save_path'])
                    }
                print('Data File: ', EKfn_dict['f'+str(ct)]['EK_data_path'] / 
                    EKfn_dict['f'+str(ct)]['EK_data_filenames'])
                ct += 1
    # line files
    # if lines path and lines filenames are specified, create a dictionary of line files
    # the json EV_lines is a nested dictionary with line name as the primary
    # key, then other parameters as keys
    if (pars_dict['sub_selection']['EV_lines']): 
        for k in pars_dict['sub_selection']['EV_lines'].keys():
            if pars_dict['sub_selection']['EV_lines'][k]['linecolor']:
                lcolor = pars_dict['sub_selection']['EV_lines'][k]['linecolor'] 
            else:
                lcolor = 'black'
            if pars_dict['sub_selection']['EV_lines'][k]['linewidth']:
                lwidth = pars_dict['sub_selection']['EV_lines'][k]['linewidth'] 
            else:
                lwidth = 1 
            filelist = []
            for f in pars_dict['sub_selection']['EV_lines'][k]['filenames']:
                filelist.append(Path(f))
            EVL_dict[k] = {
                'path' : Path(pars_dict['sub_selection']['EV_lines'][k]['path']), 
                'filenames' : filelist,
                'linecolor' : lcolor,
                'linewidth' : lwidth
                }
            for f in EVL_dict[k]['filenames']:
                print('EVL File: ', EVL_dict[k]['path'] / f)
            ct += 1


    return EKfn_dict, EVL_dict, pars_dict


def convert_to_Sv(ed, wf='CW', encode='power'):
    '''
    apply calibrations and compute Sv
    it seems that for data collected with EK60 GPTs, but with EK80 software, you
    must specify CW and power.
    input: Echopype object
           wf - waveform mode, string, CW=continuous wave
           encode - encode mode, string, power=power
    Output: xarray with Sv data variable
    '''

    print('Converting to Sv with waveform {0} and encoding {1}'.format(wf, encode))
    #Sv = ep.calibrate.compute_Sv(ed, waveform_mode='CW', encode_mode='power')
    Sv = ep.calibrate.compute_Sv(ed, waveform_mode=wf, encode_mode=encode)

    return Sv


def BNR_Sv(Sv, s=10, p=5, snr='10.0dB', maxnoise='-66dB'):
    '''
    Remove background noise from Sv data

    input:  Sv xarray data variable
    optional:
            snr - signal-to-noise ratio (SNR) for background noise reduction,
                  string of the form '3.0dB'
            s - number of sample bins (samples), integer
            p - number of transmissions (pings), integer
            maxnoise - maximum noise threshold, string of the form '-160.0dB'
    output: Sv xarray data variable
    '''

    print('Sv background noise reduction with SNR {0}, #samples {1}, #pings {2}\
            max noise {3}'.format(snr, str(s), str(p), maxnoise))

    '''
    csel = 1
    psel = 10
    p = Sv['Sv'].sel(channel=Sv.channel[csel])[psel]
    erng = Sv['echo_range'].sel(channel=Sv['echo_range'].channel[csel])[psel]
    spreading_loss = 20 * np.log10(Sv["echo_range"].where(Sv["echo_range"] >= 1, other=1))
    sl = spreading_loss.sel(channel=spreading_loss.channel[csel])[psel]
    absorption_loss = 2 * Sv["sound_absorption"] * Sv["echo_range"]
    al = absorption_loss.sel(channel=absorption_loss.channel[csel])[psel]
    tl = sl+al
    #plt.plot(erng, tl)
    TVGnoise = -100
    Svnoise = TVGnoise+tl
    Svcorr = 10*np.log10(10**(p/10)-10**(Svnoise/10))
    #TVGnoise = tl+p[len(p)-1]-tl[len(tl)-1]
    plt.plot(erng, p)
    plt.plot(erng, Svcorr)
    plt.show(block=False)

    '''


    Sv = ep.clean.remove_background_noise(Sv, range_sample_num=s, ping_num=p, 
                                    SNR_threshold=snr,
                                    background_noise_max=maxnoise)

    return Sv

def MVBS_Sv(Sv, rv='echo_range', r='1m', t='1S'):
    '''
    compute mean volume backscattering strength (MVBS) from
    Sv data
    input:  Sv xarray data variable
    optional:
            rv - range variable, must be 'echo_range' or 'depth'
            r - range resolution (m), string, e.g., '1m' for one meter 
            t - time resolution (s), string, e.g., '5S' for five seconds
    output: MVBS
    '''

    print('MVBS with range var {0}, range bin {1} m, and time bin {2} s'.format(rv, 
          str(r), str(t)))
    mvbs = ep.commongrid.compute_MVBS(Sv, range_var=rv, range_bin=r, 
                                      ping_time_bin=t)

    return mvbs


def get_frequency_list(Sv):
    '''
    get a list of the available frequencies in the Sv data
    this returns the frequency list in original order, i.e., the order of the
    frequencies in the data file, 
    input: Sv xarray data set or data array with "frequency_nominal" 
           as a dimension or coordinate
    outut: frequency list in original order
           frequency list sorted by low to high frequency
           index list that maps the original order to a sorted order
    '''

    f_Hz = list(Sv.frequency_nominal.values)
    f_Hz_sorted = sorted(f_Hz)
    print('f_Hz: ', f_Hz)
    fidx = [f_Hz.index(i) for i in f_Hz_sorted]

    return(f_Hz, f_Hz_sorted, fidx)


def display_echograms(xrds, display_var, nrow=2, ncol=1, rv='echo_range', 
                      evl=None, pars=None):
    '''
    Display coded echograms on the monitor using imshow
    input:  xrds, xarray data set wtih the code variable
            display_variable, variable to display
                              must be 'Sv' or 'PWcode'
    optional:
            nrow - number of rows per page
            ncol - number of columns per page
            rv - range variable. Must be 'echo_range' or 'depth'
            evl - dictionary of lines
            pars - dictionary of the parameters, use for plotting
    output: None
    ### this is hard coded for now to 2 plots per page
    '''
    
    # check for lines
    if (evl != None):
        # lines exist
        # number of lines
        nl = len(evl.keys())

    # get a list of frequencies
    f_Hz, f_Hz_sorted, fidx = get_frequency_list(xrds)
    nf = len(f_Hz)
    f_kHz = [f/1000 for f in f_Hz]

    # range/depth limits
    # the range values are the upper limit of the range/depth bin
    ymin = xrds[rv].values[0]
    zoomymin = ymin
    nrange = xrds.dims[rv]
    rint = float(re.findall(r"[-+]?(?:\d*\.*\d+)", 
                 xrds[display_var].attrs['range_meter_interval'])[0])
    # add the range interval to the last range
    ymax = xrds[rv].values[nrange-1]+rint
    zoomymax = ymax
    ylims = [ymin, ymax]

    # zoom the y-dimension if desired
    if (pars['plotting'][display_var]['zoom']):
        if (pars['plotting'][display_var]['depth_min'] != None):
            zoomymin = float(pars['plotting'][display_var]['depth_min'])
        if (pars['plotting'][display_var]['depth_max'] != None):
            zoomymax = float(pars['plotting'][display_var]['depth_max'])

    # time limits
    # xmit time interval
    xmin = xrds.ping_time.values[0]
    zoomxmin = xmin
    dtstr = xrds[display_var].attrs['ping_time_interval']
    dt = int(re.findall(r'\d+', dtstr)[0])
    tunit = re.findall('[a-zA-Z]+', dtstr)[0]
    # number of pings
    nping = xrds.dims['ping_time']
    xmax = xrds.ping_time.values[nping-1]
    zoomxmax = xmax
    xlims = [xmin, xmax]

    # zoom the x-dimension if desired
    if (pars['plotting'][display_var]['zoom']):
        if (pars['plotting'][display_var]['time_min'] != None):
            zoomxmin = pd.to_datetime(pars['plotting'][display_var]['time_min'])
        if (pars['plotting'][display_var]['time_max'] != None):
            zoomxmax = pd.to_datetime(pars['plotting'][display_var]['time_max'])

    # convert xlimits from time to number for the x-axis
    xlims = mdates.date2num(xlims)

    # display extent
    dextent = [xlims[0], xlims[1], ylims[1], ylims[0]] 

    # display limits
    # default values
    dBmin = -90
    dBmax = -30
    if (pars['plotting'][display_var]['dB_min'] != None):
        dBmin = pars['plotting'][display_var]['dB_min']
    if (pars['plotting'][display_var]['dB_max'] != None):
        dBmax = pars['plotting'][display_var]['dB_max']

    # other parameters
    colormap = 'jet'
    ptitle = ''
    if (pars['plotting'][display_var]['color_map']):
        colormap = pars['plotting'][display_var]['color_map']
    if (pars['plotting'][display_var]['title']):
        ptitle = pars['plotting'][display_var]['title']

    if (display_var == 'PWcode'):
        colormap = cm.get_cmap(colormap, abs(dBmax-dBmin))

    nrow = 2
    ncol = 1
    n2pg, n1pg = divmod(nf, nrow)
    ct = 0
    for i in range(n2pg):
        fig, (ax0, ax1), = plt.subplots(nrows=nrow, ncols=ncol, figsize=(9,9))
        #fig.subplots_adjust(hspace = 0.3)
        if (display_var == 'Sv'):
            cset0 = ax0.imshow(xrds[display_var].sel(frequency_nominal=f_Hz[fidx[ct]]).transpose(), 
                               cmap=colormap, vmin=dBmin, vmax=dBmax, aspect='auto', 
                               interpolation='none', extent=dextent)
        else:
            cset0 = ax0.imshow(xrds[display_var].astype(int).transpose(), 
                               cmap=colormap, vmin=dBmin, vmax=dBmax, aspect='auto', 
                               interpolation='none', extent=dextent)
        if (pars['plotting'][display_var]['zoom']):
            ax0.set_ylim(zoomymax, zoomymin)
            ax0.set_xlim(zoomxmin, zoomxmax)
        ax0.set_xlabel('Date-Time')
        ax0.set_ylabel('Depth (m)')
        ax0.set_title(str(f_kHz[fidx[ct]])+' kHz')
        ax0.xaxis_date()
        date_format = mdates.DateFormatter('%Y/%m/%d\n%H:%M:%S')
        ax0.xaxis.set_major_formatter(date_format)
        fig.autofmt_xdate()
        cbar = fig.colorbar(cset0, ax=ax0, shrink=0.75)
        cbar.set_label('dB re $m^2 m^{-3}$')
        if (nl > 0):
            for l in evl.keys():
                xline = list(evl[l]['ping_data'].time)
                yline = list(evl[l]['ping_data'].depth)
                ax0.plot(xline, yline, color=evl[l]['linecolor'], 
                        linewidth=evl[l]['linewidth'])
        ct += 1
        if (display_var == 'Sv'):
            cset1 = ax1.imshow(xrds[display_var].sel(frequency_nominal=f_Hz[fidx[ct]]).transpose(), 
                               cmap=colormap, vmin=dBmin, vmax=dBmax, aspect='auto', 
                               interpolation='none', extent=dextent)
        else:
            cset1 = ax1.imshow(xrds[display_var].astype(int).transpose(), 
                               cmap=colormap, vmin=dBmin, vmax=dBmax, aspect='auto', 
                               interpolation='none', extent=dextent)
        if (pars['plotting'][display_var]['zoom']):
            ax1.set_ylim(zoomymax, zoomymin)
            ax1.set_xlim(zoomxmin, zoomxmax)
        ax1.set_xlabel('Date-Time')
        ax1.set_ylabel('Depth (m)')
        ax1.set_title(str(f_kHz[fidx[ct]])+' kHz')
        ax1.xaxis_date()
        date_format = mdates.DateFormatter('%Y/%m/%d\n%H:%M:%S')
        ax1.xaxis.set_major_formatter(date_format)
        fig.autofmt_xdate()
        cbar = fig.colorbar(cset1, ax=ax1, shrink=0.75)
        cbar.set_label('dB re $m^2 m^{-3}$')
        if (nl > 0):
            for l in evl.keys():
                xline = list(evl[l]['ping_data'].time)
                yline = list(evl[l]['ping_data'].depth)
                ax1.plot(xline, yline, color=evl[l]['linecolor'], 
                        linewidth=evl[l]['linewidth'])

        plt.suptitle(ptitle)
        plt.show(block=False)
        if (pars['plotting'][display_var]['save_fig']):
            # path and filename
            figpath = Path(pars['path_config']['fig_path'])
            createOutDir(figpath)
            fn = Path(pars['plotting'][display_var]['fig_fn'])
            if (pars['plotting'][display_var]['zoom']):
                figfn = figpath / Path(fn.stem+'_zoom_'+str(ct)+fn.suffix)
            else:
                figfn = figpath / Path(fn.stem+'_'+str(ct)+fn.suffix)
            plt.savefig(str(figfn))

        ct += 1
    if (n1pg > 0):
        fig, (ax0, ax1), = plt.subplots(nrows=nrow, ncols=ncol, figsize=(9,9))
        #fig.subplots_adjust(hspace = 0.3)
        if (display_var == 'Sv'):
            cset0 = ax0.imshow(xrds[display_var].sel(frequency_nominal=f_Hz[fidx[ct]]).transpose(), 
                               cmap=colormap, vmin=dBmin, vmax=dBmax, aspect='auto', 
                               interpolation='none', extent=dextent)
        else:
            cset0 = ax0.imshow(xrds[display_var].astype(int).transpose(), 
                               cmap=colormap, vmin=dBmin, vmax=dBmax, aspect='auto', 
                               interpolation='none', extent=dextent)
        if (pars['plotting'][display_var]['zoom']):
            ax0.set_ylim(zoomymax, zoomymin)
            ax0.set_xlim(zoomxmin, zoomxmax)
        ax0.set_xlabel('Date-Time')
        ax0.set_ylabel('Depth (m)')
        ax0.set_title(str(f_kHz[fidx[ct]])+' kHz')
        ax0.xaxis_date()
        date_format = mdates.DateFormatter('%Y/%m/%d\n%H:%M:%S')
        ax0.xaxis.set_major_formatter(date_format)
        fig.autofmt_xdate()
        cbar = fig.colorbar(cset0, ax=ax0, shrink=0.75)
        cbar.set_label('dB re $m^2 m^{-3}$')
        if (nl > 0):
            for l in evl.keys():
                xline = list(evl[l]['ping_data'].time)
                yline = list(evl[l]['ping_data'].depth)
                ax0.plot(xline, yline, color=evl[l]['linecolor'], 
                        linewidth=evl[l]['linewidth'])

        plt.suptitle(ptitle)
        plt.show(block=False)
        if (pars['plotting'][display_var]['save_fig']):
            # path and filename
            figpath = Path(pars['path_config']['fig_path'])
            createOutDir(figpath)
            fn = Path(pars['plotting'][display_var]['fig_fn'])
            if (pars['plotting'][display_var]['zoom']):
                figfn = figpath / Path(fn.stem+'_zoom_'+str(ct)+fn.suffix)
            else:
                figfn = figpath / Path(fn.stem+'_'+str(ct)+fn.suffix)
            plt.savefig(str(figfn))

def display_echograms_hvplot(Sv):
    '''
    Display echograms on the monitor using hvplot
    input: Sv xarray data variable
    optional:
            nrow - number of rows per page
            ncol - number of columns per page
    output: None
    ### this is hard coded for now to 2 plots per page
    '''
    # get a list of frequencies
    f_Hz = list(Sv.frequency_nominal.values)
    #print('f_Hz: ', f_Hz)
    nf = len(f_Hz)
    f_kHz = [f/1000 for f in f_Hz]

    #Sv.Sv.sel(frequency_nominal=f_Hz[1]).hvplot.image(Dynamic=False)
    eg=Sv.Sv.sel(frequency_nominal=f_Hz[1]).hvplot.image(x='ping_time', y='echo_range', color='Sv',
                 rasterize=True, cmap='jet', clim=(-90, -30), xlabel='Time (UTC)', 
                 ylabel='Depth (m)').options(width=800, invert_yaxis=True)

    #hvplot.show


def read_EVL(linedict):
    '''
    read Echoview line files
    input: dictinary with Echoview line file names
    output: dictionary with date-time as keys and depth as the values
            concatenate multiple files for the same line
    '''

    for k in linedict.keys():
        evldatalist = []
        for fn in linedict[k]['filenames']:
            tmp = er.read_evl(linedict[k]['path'] / fn)
            evldatalist.append(tmp.data)
        evlpd = pd.concat(evldatalist).astype({'time':'datetime64[ns]'})
        # sort the line data to fix time reversals
        evlpd = evlpd.sort_values(by='time')
        linedict[k].update({'data':evlpd[['time', 'depth']]})

    return(linedict)


def EVL_to_Sv(evldata, Sv, lname):
    '''
    the EVL time will be aligned to the pds (primary data set)
    input: 
        Echoview line data
        Sv data set from which we'll get the coordinates
    output: Echoview line Panda dataframe
    '''

    ###
    # merge and align the EVL times and depths with the Sv times
    # using Pandas
    Svtime_df = pd.DataFrame({'time':Sv['ping_time']})
    # limit to the EVL times
    minevltime = min(evldata['data'].time)
    maxevltime = max(evldata['data'].time)
    Svtime_df = Svtime_df.query('@minevltime <= time <= @maxevltime')
    LLsub = pd.merge_asof(Svtime_df, evldata['data'], on='time', direction='nearest')

    '''
    # this creates an xarray data set
    evl_ds = xr.Dataset(
            data_vars = {
                lname : (['ping_time'], LLsub.depth)
            },
            coords = {
                'ping_time': LLsub.time,
            },
            attrs = {
                'description' : 'Echoview line'
            },
    )
    '''

    return LLsub

def get_list_index(lst, val):
    '''
    get the index of the value closest to a value in the list
    input:
        lst = list of values
        val = value to find the index in the list
    output: index of the value in the list
    '''

    indx = min(range(len(lst)), key=lambda i: abs(lst[i]-val))

    return indx

def createOutDir(outputdir):
    '''
    create the output directory
    input: output directory path
    '''
    success = True
    if (outputdir.exists()):
        print('Directory %s exists' % outputdir)
        success = True
    else:
        try:
            outputdir.mkdir()
            success = True
        except OSError:
            print('Unable to create output directory %s' % outputdir)
            success = False
            #exit()
        else:
            print('Output directory created %s' % outputdir)
            success = True
    return(success)


def fixed_line(Sv, z=10.0, ab='above', rv='echo_range'):
    '''
    set values above or below a fixed line to nan
    This is hard coded for Sv!
    input:  Echopype Sv data set
            z - range/depth in meters of the line
            ab - above/below. must be 'above' or 'below'
            rv - range variable. must be 'echo_range' or 'depth'
    output: Sv data array with values above/below the line set to nan 
    '''        
    
    if (ab == 'above'):
        # set values above a depth to nan
        Sv['Sv'] = Sv['Sv'].where(Sv[rv] >= z, drop=False)
    elif (ab == 'below'):
        # set values below a depth to nan
        Sv['Sv'] = Sv['Sv'].where(Sv[rv] <= z, drop=False)
    else:
        print('{0} is an incorrect range variable. Must be echo_range or depth'.format(rv))
     
    return(Sv)


def between_lines(Sv, line_dict, ulname, llname, rv='echo_range'):
    '''
    Slice an Sv data set based on line segments
    This is extremely hard coded for specific variables!
    This is hard coded for "between" lines
    input:  Echopype Sv data set
            Line dictionary
            line name for the upper (shallower) line
            line name for the lower (deeper) line
            rv - range variable. must be 'echo_range' or 'depth'
    output: sliced subset of the Sv data
    '''        
    ###
    # select the data between lines
    # this is klugy
    UL_da = xr.DataArray(
            dims=('ping_time'),
            coords={'ping_time':line_dict[ulname]['ping_data'].time},
            data=line_dict[ulname]['ping_data'].depth,)
    LL_da = xr.DataArray(
            dims=('ping_time'),
            coords={'ping_time':line_dict[llname]['ping_data'].time},
            data=line_dict[llname]['ping_data'].depth,)
    # subset the times
    Svsub_t = Sv.where((Sv.ping_time >= min(UL_da.ping_time)) &
                       (Sv.ping_time <= max(UL_da.ping_time)), drop=True)
    # subset the depths
    Svsub = Svsub_t.where((Svsub_t[rv] >= UL_da) &
                          (Svsub_t[rv] <= LL_da), drop=True)

    return(Svsub)


def seabed_line(Sv, line_dict, sblname, rv='echo_range'):
    '''
    apply a seabed line to the Sv echograms
    this assumes the values from the seabed line to the bottom of the Sv
    array will be set to nan
    input:  Sv - Sv xarray data set
            line_dict - line dictinary
            sblname - seabed line name
            rv - range variable, must be echo_range or depth
    output: Sv data set
    '''
    SBL_da = xr.DataArray(
             dims=('ping_time'),
             coords={'ping_time':line_dict[sblname]['ping_data'].time},
             data=line_dict[sblname]['ping_data'].depth,)
    # set values below the line to nan, this is done by setting drop=False
    Svsub = Sv.where(Sv[rv] <= SBL_da, drop=False)

    return(Svsub)


def Sv_to_dict(Svsub, rv='echo_range'):
    '''
    make a dictionary of the Sv data between the lines
    input:  Sv data set
            rv - range variable. must be 'echo_range' or 'depth'
    output: dictionary
    '''

    Svsub_dict = {}
    #for pidx in range(2):
    for pidx in range(Svsub.dims['ping_time']):
        tmpt = Svsub.isel(ping_time=pidx)
        tmptf = tmpt.isel(frequency_nominal=0)
        Svmask = ~np.isnan(tmptf.Sv)
        tmptfm = tmptf.where(Svmask, drop=True)
        for d in list(tmptfm[rv].values):
            tmp = tmpt.where(tmpt[rv] == d, drop=True)
            tmpSv = [round(x[0],2) for x in tmp.Sv.values]
            tmpc = tmp.PWcode.values[0] 
            # dictionary with time as primary key, depth/range as the
            # secondary key, and a list of Sv by frequency
            #Svsub_dict.setdefault(tmpt.ping_time.values, {}).update({er: tmp})
            Svsub_dict.setdefault(tmpt.ping_time.values, {}).update(
                    {d: {'Svf':tmpSv, 'PWcode':tmpc}})

    return(Svsub_dict)


def pair_wise(ilist, eqval):
    '''
    generate a code based on pair-wise relations
    input:  list
    output: coded list
    '''

    # TODO: check for nan!!!!
    outlist = ['2' if abs(x) <= eqval else '1' if x > eqval \
               else '3' if x < -1*eqval else '0' for x in ilist]
               #else '3' if x < -1*eqval else np.nan for x in ilist]

    return(outlist)


def PW_code(Svarr, eqval):
    '''
    label each Sv bin by the frequency response
    TODO: generalize for n-sequence combinations 
    the code is 1, 2, or 3
    nan indicates missing or below-threshold data, i.e., no data
    '1' indicates the second value is less than the first (descending) 
    '2' indicates the values are "equal" (within the eqval) (flat) 
    '3' indicates the second value is greater than the first (ascending) 
    the input dictionary is modified
    input:  Sv xarray data set
            eqval
    output: array the same size as the input Sv arrays
    '''
    # get a list of the available frequencies
    f_Hz, f_Hz_sorted, fidx = get_frequency_list(Svarr)
    nf = len(f_Hz)

    ###
    # generate codes based on frequency relationships
    # generate the list of combinations of frequency pairs
    # the "2" is used for pairs. TODO: generalise to n sequences.
    nseq = 2
    fcomb = list(itertools.combinations(range(nf),nseq))
    #fcomb[0], fcomb[0][0], fcomb[0][1]
    # sequential combinations
    fi = list(range(nf))
    fcomb_seq = list(zip(*(fi[i:] for i in range(nseq))))

    Svshape = np.shape(Svarr.Sv)
    ccode = []
    for i in range(len(fcomb_seq)):
        dSv = Svarr.Sv.isel(frequency_nominal=fcomb_seq[i][0]) - \
              Svarr.Sv.isel(frequency_nominal=fcomb_seq[i][1])
        tmp = pair_wise(dSv.data.ravel(), eqval)
        if (len(ccode) > 0):
            ccode = [''.join(x) for x in zip(ccode, tmp)]
        else:
            ccode = tmp

    if (len(Svshape) > 2):
        codearray = np.reshape(ccode, (Svshape[1], Svshape[2]))
    else:
        codearray = np.reshape(ccode, (Svshape[0], Svshape[1]))

    return(codearray)



def  write_csv(Svsub_dict, pars):
    '''
    write the codes to a file
    input:  Sv data dictionary
            pars dictionary
    output: CSV file
    '''

    # create the output directory if needed
    if (create_directory(Path(pars['path_config']['save_path']))):
        # create the output file name
        fnout = Path(pars['plotting']['Sv']['title']+'.csv')
        outfile = pars['path_config']['save_path'] / fnout 

        # write a header line
        hdr = 'DateTime,depth_m' 
        for f in pars['analysis']['frequency_list']:
            hdr = ','.join([hdr, f+'_Sv'])
        hdr = hdr+',code\n'
        # write the data
        with outfile.open('w') as outfl:
            outfl.write(hdr)
            for k in Svsub_dict.keys():
                nk = pd.to_datetime(k, unit='ns').strftime('%Y-%m-%d_%H:%M:%S.%f')[:-4]
                #print('time: {0}'.format(k))
                for j in Svsub_dict[k].keys():
                    oline = ','.join([str(nk), str(j)])
                    for f in Svsub_dict[k][j]['Svf']:
                        oline = ','.join([oline, str(f)])
                    oline = ','.join([oline, Svsub_dict[k][j]['PWcode']])
                    outfl.write(oline+'\n')
    else:
        print('Not able to write: {0}'.format(str(outfile)))


def create_directory(outdir):
    '''
    create a directory
    input:  directory name with path
    output: success/failure
    '''

    yn = True
    if (outdir.exists()):
        print('  Output Directory %s exists' % str(outdir))
    else:
        try:
            outdir.mkdir()
        except OSError:
            print('  Unable to create directory %s' % str(outdir))
            yn = False
            #exit()
        else:
            print('  Created output directory %s' % str(outdir))
            yn = True

    return(yn)


def Sv_subset_f(Svf, pars):
    '''
    Subset the Sv data set by selected frequencies
    input:  Sv xarray data set
            pars dictionary
    output: Sv data set with selected frequencies
    '''

    f_Hz, f_Hz_sorted, fidx = get_frequency_list(Svf)
    nf = len(f_Hz)
    falist = pars['analysis']['frequency_list'] 
    f_Hz_anlz = []
    for f in falist:
        blah, fq, units = re.split('(\d+)', f)
        if (units == 'kHz'):
            f_Hz_anlz.append(float(fq)*1000)
        elif (units == 'Hz'):
            f_Hz_anlz.append(float(fq))
        else:
            print('frequency units for analysis are not Hz or kHz!')

    if (len(f_Hz_anlz) > 0):
        # get the frequency_nominal indices for the select frequencies
        fidx_anlz = [f_Hz.index(i) for i in f_Hz_anlz]
        Svsub = Svf.isel(frequency_nominal = fidx_anlz) 

    return(Svsub)
        

def main(args):
    # make a list of the file names with full path
    EKfn_dict, EVL_dict, pars_dict = parse_fnames(args)
    # check that there are data files
    if (len(EKfn_dict) == 0):
        print('No files provided!')
        exit()

    # create the echo data object
    edlist = []
    for k in EKfn_dict.keys():
        fname = EKfn_dict[k]['EK_data_path'] / EKfn_dict[k]['EK_data_filenames'] 
        edlist.append(ep.open_converted(str(fname)))
    #ed = ep.combine_echodata(edlist, client=client)
    ed = ep.combine_echodata(edlist)
    #sonar_model = ed['Sonar'].attrs['sonar_model']

    # generate Sv
    Sv = convert_to_Sv(ed, wf=pars_dict['waveform'],
                           encode=pars_dict['encode'])

    # add depth (i.e., relative to the water surface)
    if (pars_dict['data_reduction']['range_var'] == 'depth'):
        Sv = ep.consolidate.add_depth(Sv, depth_offset=Sv.water_level.values)

    # new/old? version of Echopype has beam as a dimension, select beam=1
    # assume there is only one beam for now
    if [dim for dim in Sv.Sv.dims if dim == 'beam']:
        Sv = Sv.sel(beam='1').drop_vars('beam')

    # background noise reduction
    if pars_dict['noise_removal']['BNR']['clean']:
        Sv = BNR_Sv(Sv, s=pars_dict['noise_removal']['BNR']['samples'],
                        p=pars_dict['noise_removal']['BNR']['pings'],
                        snr=pars_dict['noise_removal']['BNR']['SNR'],
                        maxnoise=pars_dict['noise_removal']['BNR']['max_noise'])
        # remove the Sv array
        Sv = Sv.drop_vars('Sv')
        # rename the noise correct array to Sv
        Sv = Sv.rename({'Sv_corrected' : 'Sv'})

    # calculate MVBS
    Sv = MVBS_Sv(Sv, rv=pars_dict['data_reduction']['range_var'], 
                     r=pars_dict['data_reduction']['range_meter_bin'], 
                     t=pars_dict['data_reduction']['ping_time_bin'])

    Sv = ep.consolidate.swap_dims_channel_frequency(Sv)

    # apply a bubble layer line where values above it are nan
    if (pars_dict['data_reduction']['bubble_line_apply']):
        Sv = fixed_line(Sv, z=pars_dict['data_reduction']['bubble_line_m'],
                            ab=pars_dict['data_reduction']['bubble_line_ab'],
                            rv=pars_dict['data_reduction']['range_var'])
    # apply maximum depth lines where values below it are set to nan
    # these are frequency dependent
    #if (pars_dict['data_reduction']['maxdepth_apply']):
    #    Sv = fixed_line(Sv, z=pars_dict['data_reduction']['bubble_line_m'],
    #                        ab=pars_dict['data_reduction']['bubble_line_ab'],
    #                        rv=pars_dict['data_reduction']['range_var'])


    #HB1603_SpermWhaleDive_Span0.2_07032016_1853_UTC_depth.evl read Echoview line files
    # EVLdict has top level key of the file name and a dictionary under
    # that with time as keys and depths as values
    # use echoregions to read and manipulate the EVL file and data
    if (len(EVL_dict) > 0):
        EVL_dict = read_EVL(EVL_dict)

    # combine the lines into an Panda data frame
    # we will merge this with the Sv data, so set the EVL coordinates to the Sv
    # coordinates
    # add the dataset to the EVL_dict
    if (pars_dict['sub_selection']['align_line_time']):
        for k in EVL_dict.keys():
            evl_pd = EVL_to_Sv(EVL_dict[k], Sv, k)
            EVL_dict[k].update({'ping_data':evl_pd})

    # apply a seabed line. values below the line are nan
    if (pars_dict['data_reduction']['seabed_line_apply']):
        sblname = 'ev_bottom'
        Sv = seabed_line(Sv, EVL_dict, sblname,
                         rv=pars_dict['data_reduction']['range_var']) 

    # display echograms
    if (pars_dict['plotting']['Sv']['plot'] == True):
        display_echograms(Sv, 'Sv', nrow=2, ncol=1,
                          rv=pars_dict['data_reduction']['range_var'], 
                          evl=EVL_dict, pars=pars_dict)

    # select Sv data that match the frequencies for analysis
    if (len(pars_dict['analysis']['frequency_list']) > 0):
        Sv_sub = Sv_subset_f(Sv, pars_dict)

    # code the Sv_sub arrays
    codearr = PW_code(Sv_sub, pars_dict['analysis']['eqval'])
    # add the code array to the Sv data set
    Sv_sub = Sv_sub.assign(PWcode=(('ping_time',
                           pars_dict['data_reduction']['range_var']), codearr,
                           Sv_sub.Sv.attrs))

    # echogram
    if (pars_dict['plotting']['PWcode']['plot'] == True):
        display_echograms(Sv_sub, 'PWcode', nrow=2, ncol=1,
                          rv=pars_dict['data_reduction']['range_var'], 
                          evl=EVL_dict, pars=pars_dict)

    # get the data between the U99CI and L99CI
    ulname = 'U99CI'
    llname = 'L99CI'
    Svsub_bl = between_lines(Sv_sub, EVL_dict, ulname, llname,
                              rv=pars_dict['data_reduction']['range_var'])

    # echogram
    #if (pars_dict['plotting']['Sv']['plot'] == True):
    #    display_echograms(Sv_sub, 'Sv', nrow=2, ncol=1,
    #                      rv=pars_dict['data_reduction']['range_var'], 
    #                      evl=EVL_dict, pars=pars_dict)

    # add the subset Sv & code data to a dictionary
    Svsub_dict = Sv_to_dict(Svsub_bl, rv=pars_dict['data_reduction']['range_var'])

    # permutations
    #fperm = list(itertools.permutations(range(nf),2))

    # write the codes to a .csv file
    if pars_dict['output_to_csv']:
        write_csv(Svsub_dict, pars_dict)


if __name__=='__main__':
    import sys
    import argparse
    parser = argparse.ArgumentParser(description='Display echograms & lines')
    parser.add_argument('-dd', '--datadirectory', type=Path, help='data file directory')
    parser.add_argument('-f', '--filenames', type=Path, default=[], 
                        nargs='+', help='data file name(s): returns a list')
    parser.add_argument('-ld', '--linedirectory', type=Path, help='line file directory')
    parser.add_argument('-l', '--linenames', type=Path, default=[], 
                        nargs='+', help='line file name(s): returns a list')
    parser.add_argument('-j', '--jsonfile', type=Path, help='json file full\
                        path and file name')
    args = parser.parse_args()

    if any(vars(args).values()):
        main(args)
    else:
        print('No Arguments Provided! Try the -h option')


