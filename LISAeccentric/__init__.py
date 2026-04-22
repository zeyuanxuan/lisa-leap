import sys
import warnings

try:
    import leap
except ImportError:
    raise ImportError("Failed to find the new 'leap' package. Make sure you are in the project root.")

sys.modules['LISAeccentric'] = leap

warnings.warn(
    "LISAeccentric is deprecated. Using 'leap' instead. "
    "Please update your imports to 'import leap' to avoid this warning.",
    DeprecationWarning,
    stacklevel=2
)