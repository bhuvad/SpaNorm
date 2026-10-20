# Make read-only `stream` properties on logging handlers assignable.
#
# reticulate's output capture (rpytools/output.py) assigns `handler.stream = ...` (via
# StreamHandler.setStream) on every logger handler whose stream is the captured stdout/stderr.
# Several handler classes expose `stream` as a read-only property that tracks sys.stderr: the
# stdlib's logging._StderrHandler (logging.lastResort) and, crucially, PyTorch's own
# torch._logging._internal._StderrHandler, attached to ~40 torch.* loggers as soon as torch
# is imported (scanpy/spanormpy import it). Assigning to them raises
#   AttributeError: property 'stream' of '_StderrHandler' object has no setter
# and every Python chunk fails. This gives each such class a setter: an assigned stream is
# remembered as an override, and assigning the value the getter would return anyway (the live
# sys.stderr) clears it so the handler tracks sys.stderr again.
#
# Classes are patched both now (all existing Handler subclasses) and when later defined
# (an __init_subclass__ hook on logging.Handler), so it works even though torch is only
# imported after this file runs.
import logging

if not getattr(logging.Handler, "_rpy_stream_shim", False):

    def _patch(cls):
        attr = cls.__dict__.get("stream")
        if not isinstance(attr, property) or attr.fset is not None or attr.fget is None:
            return
        fget = attr.fget

        def _get(self, _fget=fget):
            override = self.__dict__.get("_stream_override")
            return _fget(self) if override is None else override

        def _set(self, value, _fget=fget):
            self.__dict__["_stream_override"] = None if value is _fget(self) else value

        cls.stream = property(_get, _set, attr.fdel, attr.__doc__)

    def _walk(cls):
        for sub in cls.__subclasses__():
            _patch(sub)
            _walk(sub)

    _orig = logging.Handler.__dict__.get("__init_subclass__")

    def __init_subclass__(cls, **kwargs):
        super(logging.Handler, cls).__init_subclass__(**kwargs)
        _patch(cls)

    logging.Handler.__init_subclass__ = classmethod(__init_subclass__)
    logging.Handler._rpy_stream_shim = True
    _walk(logging.Handler)
