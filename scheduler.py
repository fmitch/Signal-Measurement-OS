import sys
import pyfftw
import threading, time, signal
import yaml, pickle
import numpy as np
import uhd
import logging
import scipy.signal
import math
import matplotlib
import os.path
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.dates import HourLocator
import constants as C
import definitions as D

from datetime import timedelta
from datetime import datetime

import pdb

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

class ProgramKilled(Exception):
  pass

 
def signal_handler(signum, frame):
  raise ProgramKilled

    
def handler ():
  return

class Radio:
  def __init__ (self, sample_rate, gain, usrp):
    self.sample_rate = sample_rate
    self.gain = gain
    self.usrp = usrp
 
  def acquire_samples(self, frequency, n, nskip=0, retries=5):  # -> np.ndarray:
    """Aquire nskip+n samples and return the last n"""

    # Compute the linear gain
    db_gain = -40
    linear_gain = 10 ** (db_gain / 20.0)

    # Try to acquire the samples
    max_retries = retries
    while True:
      nsamps = n + nskip

      samples = self.usrp.recv_num_samps(
          nsamps,  # number of samples
          frequency,  # center frequency in Hz
          self.sample_rate,  # sample rate in samples per second
          [0],  # channel list
          self.gain,  # gain in dB
          )

      assert samples.dtype == np.complex64
      assert len(samples.shape) == 2 and samples.shape[0] == 1
      data = samples[0]  # isolate data for channel 0
      data_len = len(data)

      data = data[nskip:]

      if not len(data) == n:
        if retries > 0:
          msg = "USRP error: requested {} samples, but got {}."
          logger.warning(msg.format(n + nskip, data_len))
          logger.warning("Retrying {} more times.".format(retries))
          retries = retries - 1
        else:
          err = "Failed to acquire correct number of samples "
          err += "{} times in a row.".format(max_retries)
          raise RuntimeError(err)
      else:
        logger.info("Successfully acquired {} samples.".format(n))

        # Scale the data back to RF power and return it
        data /= linear_gain
        return data

  def psd_welch (self, frequency, nfft, samples, Fs=30.72e6):
    """ use welch's method to estimate PSD
    """
    center_freq = frequency

    labels0, psd0lin = scipy.signal.welch(samples, Fs, nperseg=nfft, return_onesided=False)
    #pdb.set_trace()

    psd0 = np.nan_to_num(10.0 * np.log10(psd0lin))
    labels0 = labels0 + center_freq

    # put DC in the center
    psd0 = np.fft.fftshift(psd0)
    labels0 = np.fft.fftshift(labels0)

    logger.info ("Each bin : {}".format (labels0[1]-labels0[0]))
    return psd0, labels0


  def plot_psd (self, psd, freq):
       plt.plot (freq, psd)
       plt.ylim(-90, -30)
       plt.xlabel('frequency [Hz]')
       plt.ylabel('PSD [dB]')
       plt.grid(True)
       plt.show()

  def plot_live_psd (self, fig, ax, psd, freq):
      freq[:] = [x / 1e6 for x in freq]
      ax.clear()
      ax.set_xlabel('frequency [MHz]')
      ax.set_ylabel('PSD [dB]')
      ax.plot(freq, psd)
      fig.canvas.draw()


def get_threshold_from_file (thres_file, loc_name):
  if (os.path.exists (thres_file) == False):
    print ("No such threshold file")
    return

  thres_dict = {}
  with open (thres_file) as tfp:
      line = tfp.readline ()

      while (line):
        if (line.isspace()):
          line = tfp.readline ()
          continue

        line = line.rstrip ("\n")
        if (line[0] != '#'):
          line_lst = line.split (",")
          if (line_lst[0] == loc_name):
            key = line_lst[1] + "-" + line_lst[2]
            
            thres_dict[key] = []
            thres_dict[key].append (line_lst[0])
            for i in range (3, 6):
              thres_dict[key].append (line_lst[i]) 
             
        line = tfp.readline ()

  return thres_dict
  

def threshold_lookup_from_dict (thres_dict, fcs_s_mz):
  for key, val in thres_dict.items ():
    key_s_e = key.split ("-")

    if (fcs_s_mz >= float (key_s_e[0]) and fcs_s_mz <= float (key_s_e[1])):
      return float (val[3])
  
 
def save_frequency_dict (lfreq_dict):
  plot_lst = [] 
  for key, val in lfreq_dict.items():
    plot_lst = plot_lst + val.freq_energy_lst
    
  plot_np_lst = np.array(plot_lst)

  start_time = plot_np_lst[0][0].strftime('%y%m%d%H%M%S')
  end_time = plot_np_lst[-1][0].strftime('%y%m%d%H%M%S')
  with open('data/data_%s_%s.pkl' % (start_time, end_time),'ab+') as f:
    pickle.dump (plot_np_lst, f)



def plot_usage_scatter_plot (lfreq_dict):
  plot_lst = [] 
  for key, val in lfreq_dict.items():
    plot_lst.append (val.freq_energy_lst[0])
    
  plot_np_lst = np.array(plot_lst)

  #Uncomment it for debugging
  #pdb.set_trace()
  plot_np_lst = plot_np_lst [plot_np_lst[:,5] == 1]
  x = plot_np_lst[:, 1]
  y = plot_np_lst[:, 0]

  ax = plt.subplot()
  ax.plot(x, y, 'o', color='red')
  plt.grid(True)
  ax.yaxis.set_major_locator(HourLocator())
  ax.yaxis.set_major_formatter(DateFormatter('%H:%M'))
  plt.show()

  return       


def sweep_frequencies (measurement_list, map_list, top, plot_psd_flag):
  sample_rate = float (measurement_list ['sample_rates'])
  measurement_duration = measurement_list['duration_ms'] / 1000
  num_sub_intervals = measurement_list['duration_ms'] // measurement_list['sub_duration_ms']
  if measurement_list['duration_ms'] < measurement_list['sub_duration_ms']:
      num_sub_intervals = 1
  nsamps = int(sample_rate * measurement_list['duration_ms'] * 1e-3)
  usrp_args = measurement_list['usrp_args']
  gain = measurement_list['gains']
  fcs_start = float (measurement_list ['fcs_start'])
  fcs_end = float (measurement_list ['fcs_end'])

  loc_type = measurement_list ['loc_type']
  #pdb.set_trace()

  if (loc_type == 'static'):
    loc = measurement_list['loc']
    top.loc_dict[loc] = D.location ()

    loc_obj = top.loc_dict[loc]
    loc_obj.loc_name = loc
    loc_obj.gps_coordinate = map_list[loc]


  usrp = top.usrp

  # Drop ~10 ms of samples
  nskip = int(0.01 * sample_rate)

  # Arrange the frequencies in the range
  freq_lst =  np.arange(fcs_start, fcs_end, C.ANALOG_BW)

  radio = Radio (sample_rate, gain, usrp)
  loc_obj.radio = radio

  if (plot_psd_flag):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.ion()

    fig.show()
    fig.canvas.draw()

  fcs_bw = float(measurement_list['fcs_diff'])
  for idx, center_fcs_val in enumerate (freq_lst):
    center_fcs = center_fcs_val + C.ANALOG_BW/2
    if (center_fcs <= fcs_end):
      no_bins = int (int(C.ANALOG_BW)/fcs_bw)
      fcs_bw = int(fcs_bw)

      for s_idx in range (no_bins):

        fcs_s = center_fcs_val + fcs_bw*s_idx + fcs_bw/2
        fcs_s_mz = float(fcs_s)/1e6

        lfreq_dict = loc_obj.freq_dict

        if (fcs_s_mz not in lfreq_dict):
          lfreq_dict[fcs_s_mz] = D.frequency()

        lfreq_dict[fcs_s_mz].center_fcs = fcs_s_mz
        lfreq_dict[fcs_s_mz].bandwidth = fcs_bw

  pdb.set_trace()
  for idx, center_fcs_val in enumerate (freq_lst):
    center_fcs = center_fcs_val + C.ANALOG_BW/2

    if (center_fcs > fcs_end):
      save_frequency_dict (lfreq_dict)  
      return  

    time_start = datetime.now()
    #Time domain data
    t0 = time.time()
    acq = radio.acquire_samples(center_fcs, nsamps, nskip=nskip).astype(np.complex64)
    time_end = datetime.now()
    time_delta = (time_end - time_start) / num_sub_intervals
    t1 = time.time()

    #Frequency domain data
    ####f_acq = (np.fft.fft(acq))/nsamps
    int_size = int(nsamps / num_sub_intervals)
    print('Intervals', num_sub_intervals)
    for sub_int in range(num_sub_intervals):
      f_acq = (pyfftw.interfaces.numpy_fft.fft(acq[sub_int*int_size:(sub_int+1)*int_size], threads=2))/(nsamps/num_sub_intervals)

      #T = 1 * (nsamps/num_sub_intervals)/sample_rate
      #dw = 2*math.pi/T

      # put DC in the center
      f_acq = np.fft.fftshift(f_acq)

      for s_idx in range (no_bins):
        start_i = int(s_idx * fcs_bw * measurement_duration)
        end_i = int((fcs_bw + s_idx*fcs_bw) * measurement_duration - 1)

        fcs_s = center_fcs_val + fcs_bw*s_idx + fcs_bw/2
        fcs_s_mz = float(fcs_s)/1e6

        temp_lst = []
        temp_lst.append (time_start+time_delta*sub_int)
        temp_lst.append (fcs_s_mz)
        temp_lst.append (fcs_bw)
        P_fft = np.sum(np.abs(f_acq[start_i:end_i+1])**2)/fcs_bw

        P_fft_log = np.nan_to_num (10.0 * np.log10(P_fft))

        temp_lst.append (P_fft_log)

        #Dummy threshold and signal presence values, updated in plot_graph
        temp_lst.append (0)
        temp_lst.append (0)


        lfreq_dict[fcs_s_mz].freq_energy_lst.append(temp_lst)

    logger.info("DFT Power of freq {} MHz = {}, {} dbm".format(fcs_s_mz, P_fft, P_fft_log))
    t4 = time.time()
    print('Measurement Shape', acq.shape)
    print('Measurement Time', t1-t0)
    print('Processing Time', t4-t1)

    #plot live radio psd
    if (plot_psd_flag):
      S_xx_welch, freq = radio.psd_welch (center_fcs, 2**10, acq, sample_rate)
      radio.plot_live_psd (fig, ax, S_xx_welch, freq)


  pdb.set_trace()
  save_frequency_dict (lfreq_dict)  
  #plot_usage_scatter_plot (lfreq_dict)
  return 


class Job(threading.Thread):
  def __init__(self, interval, sleep_time, execute, timer, *args, **kwargs):
    threading.Thread.__init__(self)
    self.daemon = False
    self.stopped = threading.Event()
    self.interval = interval
    self.execute = execute
    self.args = args
    self.kwargs = kwargs
    self.timer = timer
    self.start_time = None
    self.sleep_time = sleep_time 
        
  def stop(self):
    self.stopped.set()
    self.join()

  def run(self):
    while (time.time()- self.start_time <= self.interval):
      self.execute(*self.args, **self.kwargs)
      time.sleep (self.sleep_time)



if __name__ == "__main__":

    len_arg = len (sys.argv)
    plot_psd_flag = False

    #Enter command python3.6 scheduler.py -p
    if (len_arg > 1):
      plot_psd_flag = True

    with open(r'measurement.yaml') as file:
      # The FullLoader parameter handles the conversion from YAML
      # scalar values to Python the dictionary format
       measurement_list = yaml.load(file, Loader=yaml.FullLoader)

    interval = measurement_list['interval']
    sleep_time = measurement_list['sleep_time'] 
    usrp_args = measurement_list['usrp_args']
    timer = threading.Timer(interval, handler)
    top = D.top ()

    with open(r'map.yaml') as file:
      # The FullLoader parameter handles the conversion from YAML
      # scalar values to Python the dictionary format
      map_list = yaml.load(file, Loader=yaml.FullLoader)


    try:
      usrp = uhd.usrp.MultiUSRP(usrp_args)
    except RuntimeError:
      err = "No device found matching search parameters {!r}\n"
      err = err.format(usrp_args)
      raise RuntimeError(err)

    usrp.set_rx_gain (40, 0)
    logger.info ("Using the following USRP:")                             
    logger.info (usrp.get_pp_string())
    logger.info ("RX-Bandwidth: {}".format (usrp.get_rx_bandwidth()))
    logger.info ("RX-Bandwidth Range: {}".format (usrp.get_rx_bandwidth_range()))
    logger.info ("RX-Gain: {}".format (usrp.get_rx_gain()))
    logger.info ("RX-Rate: {}".format (usrp.get_rx_rate()))
    logger.info ("RX-Info: {}".format (usrp.get_usrp_rx_info(1)))
    logger.info ("RX-Antenna: {}".format (usrp.get_rx_antenna(0)))

    top.usrp = usrp

    job = Job(interval, sleep_time, sweep_frequencies, timer, measurement_list, map_list, top, plot_psd_flag)

    job.start_time = time.time()
    timer.start()
    job.start()

"""
    signal.signal(signal.SIGTERM, signal_handler)
    signal.signal(signal.SIGINT, signal_handler)

    while True:
          try:
              time.sleep(1)
          except ProgramKilled:
              print ("Program killed: running cleanup code")
              job.stop()
              break
"""
