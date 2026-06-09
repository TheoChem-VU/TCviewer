class Atom:
    @register_setting('Material', 'Color', 'color', 
        limits={'R': (0, 255), 'G': (0, 255), 'B': (0, 255)}, 
        defaults={'radius': AtomSpecificColor()})
    def set_color(self, R: int = None, G: int = None, B: int = None):
        if R is not None:
            self.settings['color'][0] = R/255
        if G is not None:
            self.settings['color'][1] = G/255
        if B is not None:
            self.settings['color'][2] = B/255
        self._reset_properties()

def register_setting(category: str, variable_name: str, internal_name: str, limits: dict = None, defaults=None):
    # print(category, variable_name, internal_name, limits, defaults)
    def inner_dec(func):
        cls = func.__qualname__.split('.')[-2]
        if defaults is None:
            signature = inspect.signature(func)
            new_defaults = {k: v.default for k, v in signature.parameters.items() if v.default is not inspect._empty}
        else:
            new_defaults = defaults
        _registered_settings.append((category, variable_name, internal_name, limits, new_defaults, func, cls))
        return func
    return inner_dec


