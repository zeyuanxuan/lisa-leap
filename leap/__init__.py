# filepath leap/__init__.py

# import from core.py
try:
    from .core import LISAeccentric as _CoreEngine
    from .core import CompactBinary

    from .core import set_output_control, set_verbose

except ImportError as e:
    raise ImportError(f"LISAeccentric package initialization failed. Could not import 'core.py'.\nDetails: {e}")


_default_instance = _CoreEngine()

GN = _default_instance.GN
GC = _default_instance.GC
Field = _default_instance.Field
Waveform = _default_instance.Waveform
Noise = _default_instance.Noise
getMWcatalog = _default_instance.getMWcatalog

CompactBinary = CompactBinary

# public interface
__all__ = [
    'GN',
    'GC',
    'Field',
    'Waveform',
    'Noise',
    'CompactBinary',
    'set_output_control',
    'set_verbose',
    'getMWcatalog'
]

# print("leap package initialized successfully.")