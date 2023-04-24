from .ddn import DictDotNotation as ddn
from .utils import *
from copy import copy, deepcopy
import numpy as np



class CO2DspExp:
    """
    CO2AbsExp like data trees, apis, processings
    """

    """ utilities """
    def __init__(self):
        """
        args:
        - None
        
        returns:
        - None
        """
        self.exp_type = "Desorption"
        self.original = ddn()
        self.processed = ddn()
        self.result = ddn()
        self.format = False
        self.proc = False
        self.data_len = None
        self.data_imported = ddn({"time": False, "gas_conc": False, "flow_rate": False})
        self.calculated = ddn({"diff": False, "dsp": False})
        self.form_setting = []
    
    def delete_analysed(self):
        """
        delate processed data.

        args:
        - None
        
        returns:
        - None
        """
        self.processed = ddn()
        self.calculated = ddn({"diff": False, "dsp": False})
        self.format = False
        self.proc = True
        return self
    
    def get_total_dsp(self):
        if not self.calculated.dsp:
            self.full_analyse()
        return self.result.co2_desorbed.integral_co2_dsp[-1]*1000
    
    
    """ data import section """
    def set_time_axis(self, data):
        """
        Import time data. provided data will converted to "m" (minutes) unit.

        args:
        - data : specify InsightCsvHandler.elapsed
        
        returns:
        - self
        """
        if self.proc:
            raise AnalyticalException("Processed data already exist, please make new instance for other analysis")
        unit = data.unit
        elapsed = data.elapsed
        if self.data_len is not None:
            if self.data_len != len(elapsed):
                raise AnalyticalException("Unmatch data length. Data length should be " + str(self.data_len))
        if unit not in ["s", "m", "h"]:
            raise AnalyticalException("Unknow time unit")
        t_axis = ddn({
            "unit": "m",
            "data": change_time_unit(elapsed, unit, "m")
        })
        self.original.elapsed = t_axis
        #self.processed.elapsed = t_axis
        if self.data_len is None:
            self.data_len = len(elapsed)
        self.data_imported.time = True
        return self
    
    def set_gas_conc(self, upstream, upstream_unit, downstream, downstream_unit):
        """
        Import CO2 concentration data. provided data will converted to "%" unit.

        args:
        - upstream : list of float, upstream CO2 concentration
        - upstream_unit : CO2 concentration unit of upstream flow
        - downstream : list of float, downstream CO2 concentration
        - downstream_unit : CO2 concentration unit of downstream flow
        
        returns:
        - self
        """
        if self.proc:
            raise AnalyticalException("processed data already exist, please make new instance for other analysis")
        if isinstance(upstream, (int, float)):
            upstream = [upstream]*len(downstream)
        if self.data_len is not None:
            if (self.data_len != len(upstream)) or (self.data_len != len(downstream)):
                raise AnalyticalException("Unmatch data length. Data length should be " + str(self.data_len))
        if len(upstream) != len(downstream):
            raise AnalyticalException("Unmatch data length between given two datas.")
        if (upstream_unit not in ["ppm", "%"]) or (downstream_unit not in ["ppm", "%"]):
            raise AnalyticalException("Unknow time unit")
        gas_dict = ddn()
        upstream_percent = change_conc_unit(upstream, upstream_unit, "%")
        downstream_percent = change_conc_unit(downstream, downstream_unit, "%")
        gas_dict.upstream = ddn({"unit": "%", "range": get_range(upstream_percent), "data": upstream_percent})
        gas_dict.downstream = ddn({"unit": "%", "range": get_range(downstream_percent), "data": downstream_percent})
        self.original.gas_conc = gas_dict
        if self.data_len is None:
            self.data_len = len(upstream)
        self.data_imported.gas_conc = True
        return self
    
    def set_flow_rate(self, flow_rate):
        """
        Import flow rate data

        args:
        - flow_rate : int or float or list of numbers, in cc/m unit
        
        returns:
        - self
        """
        if self.proc:
            raise AnalyticalException("processed data already exist, please make new instance for other analysis")
        if isinstance(flow_rate, (int, float)):
            if self.data_len is None:
                raise AnalyticalException("Data length is unknow")
            else:
                flow_rate = [flow_rate]*self.data_len
        if self.data_len is not None:
            if self.data_len != len(flow_rate):
                raise AnalyticalException("Unmatch data length. Data length should be " + str(self.data_len))
        self.original.flow_rate = flow_rate
        #self.processed.flow_rate = flow_rate
        if self.data_len is None:
            self.data_len = len(flow_rate)
        self.data_imported.flow_rate = True
        return self
    
    
    """ data format section """
    def set_normalize_setting(self, min_value, max_value, unit):
        """
        Add normalize setting to self.form_setting
        use min-max normalization

        args:
        - min_value : min
        - max_value : max
        - unit : unit of above 2 arguments

        returns:
        - self
        """
        self.form_setting.append(ddn({
            "type": "normalize",
            "min_value": min_value,
            "max_value": max_value,
            "unit": unit
        }))
        return self

    def set_normalize_setting(self, stream_type ,min_value, max_value, unit):
        """
        Add normalize setting to self.form_setting
        use min-max normalization

        args:
        - stream_type : "up" or "down"
        - min_value : min
        - max_value : max
        - unit : unit of above 2 arguments

        returns:
        - self
        """
        self.form_setting.append(ddn({
            "type": "normalize",
            "stream": stream_type,
            "min_value": min_value,
            "max_value": max_value,
            "unit": unit
        }))
        return self
    
    def exec_format(self):
        if all([i for i in self.data_imported.values()]) is False:
            raise AnalyticalException("some data was lacking. "+str(self.data_imported))
        
        self.processed = self.original.copy()
        for i in self.form_setting:
            if i.type == "normalize":
                self._set_normalize(i.stream, i.min_value, i.max_value, i.unit)
            elif i.type == "interval_average":
                self._set_interval_average(i.interval, i.unit)
        self.result.format_setting = deepcopy(self.form_setting)
        self.form_setting = []
        self.format = True
        return self


    def _set_normalize(self, stream, min_value, max_value, unit):
        setting_value_percent = change_conc_unit([min_value, max_value], unit, "%")
        min_value_percent = setting_value_percent[0]
        max_value_percent = setting_value_percent[1]

        bef_conc = self.processed.gas_conc
        if stream == "up":
            bef_min = bef_conc.upstream.range.min
            bef_diff = bef_conc.upstream.range.diff
            stream_data = bef_conc.upstream.data
        elif stream == "down":
            bef_min = bef_conc.downstream.range.min
            bef_diff = bef_conc.downstream.range.diff
            stream_data = bef_conc.downstream.data
        norm = [(i-bef_min)/bef_diff*(max_value_percent-min_value_percent)+min_value_percent if i is not None else None for i in stream_data]
        update_value = ddn({
                "unit": "%",
                "range": ddn({"max": max_value_percent, "min": min_value_percent, "diff": max_value_percent-min_value_percent}),
                "data": norm
        })
        
        if stream == "up":
            self.processed.gas_conc.upstream = update_value
        elif stream == "down":
            self.processed.gas_conc.downstream = update_value
        return self

    def _set_interval_average(self, interval, unit):
        proc = self.processed
        
        interval_min = change_time_unit([interval], unit, "m")[0]
        
        # 時間でくぎる
        data_separation = {}
        for t, up, down, fr in zip(proc.elapsed.data, proc.gas_conc.upstream.data, proc.gas_conc.downstream.data, proc.flow_rate):
            interval_t = (t//interval_min)*interval_min
            try:
                data_separation[interval_t]
            except KeyError:
                data_separation[interval_t] = ddn({"up": [], "down": [], "fr": []})
            data_separation[interval_t].up.append(up)
            data_separation[interval_t].down.append(down)
            data_separation[interval_t].fr.append(fr)
        
        # 区切った後にそれぞれ平均する
        conv = ddn({"t": [], "up": [], "down": [], "fr":[]})
        for interval, prop in data_separation.items():
            conv.t.append(interval)
            conv.up.append(np.mean([i for i in prop.up if i is not None]))
            conv.down.append(np.mean([i for i in prop.down if i is not None]))
            conv.fr.append(np.mean(prop.fr))
        
        proc.elapsed.data = conv.t
        proc.gas_conc.upstream.data = conv.up
        proc.gas_conc.downstream.data = conv.down
        proc.flow_rate = conv.fr
        
        return self
    
    
    """ calculate section """
    def calc_co2_diff(self):
        if self.format is False:
            raise AnalyticalException("Data format is not executed. please execute self.exec_format() even if there are no data formatting.")
        def calc_diff(up, down):
            try:
                return down-up
            except TypeError:
                # none value handle
                return None
        diff_list = [calc_diff(up, down) for up, down in zip(self.processed.gas_conc.upstream.data, self.processed.gas_conc.downstream.data)]

        self.result.co2_diff = ddn({
            "unit": "%",
            "data": diff_list
        })
        self.calculated.diff = True
        return self
        

    def calc_co2_dsp(self):
        if self.format is False:
            raise AnalyticalException("Data format is not executed. please execute self.exec_format() even if there are no data formatting.")
        if self.calculated.diff is False:
            self.calc_co2_diff()
        proc = self.processed
        
        partial_list = [0]
        integral_list = [0]
        formrate_list = [0]
        integral = 0
        
        length = len(proc.elapsed.data)
        
        for i in range(0, length-1):
            flow = np.mean(proc.flow_rate[i:i+2])
            elp = (proc.elapsed.data[i+1]-proc.elapsed.data[i]) #min
            co2_conc = np.mean([0 if k is None else k for k in self.result.co2_diff.data[i:i+2]])
            form_rate = co2_conc*flow/22.4/(10**5) #mol/min
            formrate_list.append(form_rate)
            partial = elp*form_rate #mol
            partial_list.append(partial)
            integral += partial
            integral_list.append(integral)
        
        self.result.co2_desorbed = ddn({
            "partial_co2_dsp": partial_list,
            "integral_co2_dsp": integral_list,
            "form_rate": formrate_list
        })
        self.calculated.dsp = True
        return self

    
    def full_analyse(self):
        self.calc_co2_dsp()
        return self
    
    def __add__(self, other):
        return self.get_total_dsp()+other.get_total_dsp()
    
    def __sub__(self, other):
        return self.get_total_dsp()-other.get_total_dsp()

    """ export setting """
    """
    def set_time_unit(self, unit):
        if not self.data_imported.time:
            raise AnalyticalException("time information was lacking.")
        
        if unit not in ["s", "m", "h"]:
            raise AnalyticalException("Unknow time unit")
        
        processed = self.change_time_unit(self.processed.elapsed.data, self.processed.elapsed.unit, unit)
        
        self.processed.elapsed = ddn({
            "unit": unit,
            "data": processed
        })
        return self
    
    def set_conc_unit(self, unit):
        if not self.data_imported.gas_conc:
            raise AnalyticalException("CO2 concentration data was lacking.")
        
        bef_conc = self.processed.gas_conc
        bef_conc.upstream.data = change_conc_unit(bef_conc.upstream.data, bef_conc.upstream.unit, unit)
        bef_conc.upstream.unit = unit
        bef_conc.downstream.data = change_conc_unit(bef_conc.downstream.data, bef_conc.downstream.unit, unit)
        bef_conc.downstream.unit = unit
        return self
    """