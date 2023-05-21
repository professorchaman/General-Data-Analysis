import spe2py as spe
from io import StringIO
import spe_loader as sl
import matplotlib.pyplot as plt
import numpy as np

class spe_reader():
    
    def __init__(self, file_name) -> None:
        self._file_name = file_name
        self.spe_files = sl.load_from_files([self._file_name])
        
        self.raw_filepath = None
        self.bkg_data = None
        self.exposure_time = None
        self.accumulation_value = None
        self.accumulation_method = None
        self.grating_groove_value = None
        self.center_wavelength = None

    
    def get_background_info(self, metadata):
        
        filepath = self.spe_files.filepath
        self.raw_filepath = filepath[:-4] + "-raw.spe"
        # print(raw_filepath)
        
        bkg_class = metadata.SpeFormat.DataHistories.DataHistory.Origin.Experiment.Devices.Cameras.Camera.Experiment.OnlineCorrections.BackgroundCorrection
        bkg_clt = bkg_class.Enabled.cdata
        
        if bkg_clt == 'True':
            raw_file = sl.load_from_files([raw_filepath])

            raw_data = self.raw_filepath.data[0][0][0]
            self.bkg_data =  self.spe_files.data[0][0][0] - raw_data

            # bkg_data = bkg_class.ReferenceFile.cdata
        else:
            self.bkg_data = None
        
        return self.bkg_data

    def get_exposure_time(self, metadata):
        self.exposure_time = metadata.SpeFormat.DataHistories.DataHistory.Origin.Experiment.Devices.Cameras.Camera.ShutterTiming.ExposureTime.cdata
        return self.exposure_time

    def get_accumulation_info(self, metadata):
        acc_class = metadata.SpeFormat.DataHistories.DataHistory.Origin.Experiment.Devices.Cameras.Camera.Experiment.OnlineProcessing.FrameCombination
        self.accumulation_value = acc_class.FramesCombined.cdata
        self.accumulation_method = acc_class.Method.cdata
        
        return self.accumulation_value, self.accumulation_method

    def get_grating_info(self, metadata):
        
        grating_class = metadata.SpeFormat.DataHistories.DataHistory.Origin.Experiment.Devices.Spectrometers.Spectrometer.Grating
        self.center_wavelength = float(grating_class.CenterWavelength.cdata)
        self.grating_groove_value = int(self.grating_grove(grating_class))
        #print(grating_grove(grating_class) + ' g/mm')
        #print(cw_val + ' nm')
        
        return self.grating_groove_value, self.center_wavelength

    def grating_grove(self, grating_class):
        grating_str = grating_class.Selected.cdata
        
        string_start = '['
        string_end = ']'
        string_to_search = ','

        grating_start = grating_str.find(string_to_search,grating_str.find(string_start),grating_str.find(string_end))+1
        grating_end = grating_str.find(string_end)
        grating_val = grating_str[grating_start:grating_end]
        #print(grating_val)
        
        return grating_val

if __name__ == "main":
    pass