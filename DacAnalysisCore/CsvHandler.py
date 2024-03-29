import os, datetime, time, csv, re
from itertools import zip_longest
from .ddn import DictDotNotation as ddn
from .exceptions import *
from .utils import get_range


class InsightCsvHandler:
    """
    InsightCsvHandlar, for parse usage.

    - self.elapsed : time information
    - self.probes : probe list (see InsightCsvHandler.get_probe_info.__doc__)
    """
    
    def __init__(self, path):
        """
        args:
        - path : csv path of insight exported csv.
        
        returns:
        - None
        
        raises:
        - AnalyticalException : When parse was failed.
        """
        self.datapath = os.path.abspath(path)
        self.rawdata = self.open_csv(self.datapath)
        self.check_csv(self.rawdata)
        data_starts, amount_probe = self.get_datastarts(self.rawdata)
        t = list(zip_longest(*self.rawdata[data_starts-1:], fillvalue=""))
        self.timestamp = [datetime.datetime.strptime(i, "%Y-%m-%d %H:%M:%S") for i in list(t[0])[1:]]
        time_start = self.timestamp[0]
        self.elapsed = ddn({
            "unit": "s",
            "elapsed": [(i-time_start).total_seconds() for i in self.timestamp]
        })
        self.probes = self.get_probe_info(self.rawdata, t, data_starts, amount_probe)
        self.metadata = ddn({
            "insight_version": self.rawdata[0][1],
            "generated": datetime.datetime.strptime(self.rawdata[1][1], "%Y-%m-%d %H:%M:%S"),
            "started": datetime.datetime.strptime(self.rawdata[2][1], "%Y-%m-%d %H:%M:%S"),
            "ended": datetime.datetime.strptime(self.rawdata[3][1], "%Y-%m-%d %H:%M:%S"),
            "sampling_point": len(self.elapsed.elapsed)
        })
    
    def check_csv(self, rawdata):
        """
        Check provided csv was generated by Vaisala insight or not.

        args:
        - rawdata : 2 dimentional list, return value of InsightCsvHandler.open_csv(path)
        
        returns:
        - None
        """
        if "Insight" not in rawdata[0][0]:
            raise AnalyticalException("Only Insight exported csv file is acceptable.")
        return
        
    @staticmethod
    def open_csv(path):
        """
        args:
        - path : path of the insight exported csv file
        
        returns:
        - all_data : 2 dimentional list of csv content
        """
        with open(path, "r", encoding='utf_8_sig') as f:
            csvreader = csv.reader(f)
            all_data = [row for row in csvreader]
        return all_data
    
    @staticmethod
    def get_datastarts(rawdata):
        """
        arg:
        - rawdata : 2 dimentional list, return value of InsightCsvHandler.open_csv(path)
        
        returns:
        - datastarts : data start index of inputted raw data
        - amount_probe : detected probe quantity
        """
        for count, i in enumerate(rawdata):
            try:
                datetime.datetime.strptime(i[0], "%Y-%m-%d %H:%M:%S")
            except (ValueError, IndexError):
                pass
            else:
                return count, (count-6)//2
    
    @staticmethod
    def get_probe_info(rawdata, zipped_data, data_starts, probe_amount):
        """
        arg:
        - rawdata : 2 dimentional list, return value of InsightCsvHandler.open_csv(path)
        - zipped_data : data that contains measurement data
        - data_starts : data start index of inputted raw data, return value of InsightCsvHandler.get_datastarts(rawdata)
        - probe_amount : detected probe quantity, return value of InsightCsvHandler.get_datastarts(rawdata)
        
        returns:
        - probes : list of probe infomation, see below for what infomation is include
            -probe : "GMP251" or "GMP252"
            -serial : serial code of probe
            -probe_version : I don't know how this number works
            -unit: unit of the values, "%" or "ppm"
            -value: list of measurement values
        """
        probes = []
        for i in range(probe_amount):
            probe_info = ddn()
            probe_info_index = 5+i*2
            probe_info.probe = rawdata[probe_info_index][0]
            probe_info.serial = rawdata[probe_info_index][1]
            probe_info.probe_version = rawdata[probe_info_index][2]
            zipped = list(zipped_data[1+i])
            probe_info.unit = re.findall("\(.+?\)", zipped[0])[0].replace("(", "").replace(")", "")
            probe_info.label = zipped[0]
            probe_info.value = [InsightCsvHandler.num_convert(j) for j in zipped[1:]]
            probe_info.range = get_range(probe_info.value)
            probes.append(probe_info)
        return probes
    
    @staticmethod
    def num_convert(value):
        """
        arg
        - value : string raw data of measurement value
        
        returns:
        - r : converted float data. If "null" was given, returns None.
        """
        try:
            r = float(value)
        except ValueError:
            return None
        else:
            return r