""" applied dict definition to use dot notation """

class DictDotNotation():
    def __init__(self, *args, **kwargs):
        for i in args:
            if not isinstance(i, dict):
                raise ValueError(f"To instantiate 'DictDotNotation' object, argument must be 'dict' object")
            for key, val in i.items():
                if isinstance(val, dict):
                    val = DictDotNotation(val)
                self.__setatrr__(key, val)
        for key, val in kwargs.items():
            if isinstance(val, dict):
                val = DictDotNotation(val)
            self.__setatrr__(key, val)
    
    def __setatrr__(self, key, value):
        setattr(self, key, value)
    
    def keys(self):
        return self.__dict__.keys()

    def values(self):
        return self.__dict__.values()

    def items(self, *args, **kwargs):
        return self.__dict__.items(*args, **kwargs)
    
    def __repr__(self):
        return f"DictDotNotation({self.__dict__})"
    
    def __str__(self):
        return str(self.to_flatdict())
    
    def to_flatdict(self):
        return_dict = {}
        for key, value in self.__dict__.items():
            if isinstance(value, DictDotNotation):
                value = value.to_flatdict()
            return_dict[key] = value
        return return_dict
    
    def copy(self):
        return DictDotNotation(self.to_flatdict())