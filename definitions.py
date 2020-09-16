class top: 
  def __init__ (self):
    self.loc_dict = {}
    self.usrp = None


class location:
  def __init__ (self):
    self.loc_name = None
    self.gps_coordinate = None
    self.radio = None
    self.freq_dict = {}
    self.threshold_dict = {}



class frequency:
  def __init__ (self):
    self.center_fcs = None
    self.bandwidth = None
    self.freq_energy_lst = []
