"""Calls into spglib that behave the same on spglib 2.4 to 3.x.

spglib >= 2.7 warns on every call ("Set OLD_ERROR_HANDLING to false ...") because it is
moving from returning ``None`` on failure to raising ``SpglibError``. The callers here
handle a ``None`` result, so failures are mapped to ``None`` on every version and that
one warning is silenced locally (spglib's global setting is left alone, so other code in
the same process is not affected).
"""
import warnings


def call(name, *args, **kwargs):
    """Return ``spglib.<name>(*args, **kwargs)``, or ``None`` if spglib fails."""
    import spglib

    error = getattr(getattr(spglib, "error", None), "SpglibError", None) or getattr(spglib, "SpglibError", None)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=".*OLD_ERROR_HANDLING.*", category=DeprecationWarning)
        try:
            return getattr(spglib, name)(*args, **kwargs)
        except Exception as exc:
            if error is not None and isinstance(exc, error):
                return None
            raise
