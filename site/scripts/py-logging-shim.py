# Make logging._StderrHandler.stream assignable.
#
# reticulate's output capture (rpytools/output.py) assigns `handler.stream = ...` on every
# logging handler whose stream is the captured stdout/stderr. The stdlib's
# logging.lastResort is a _StderrHandler whose `stream` is a read-only property, so as soon
# as any package attaches it to a logger (and sys.stderr is aliased to the captured stream)
# every Python chunk fails with "property 'stream' of '_StderrHandler' object has no setter".
# This gives the property a setter: an assigned stream is remembered as an override, and
# assigning the live sys.stderr (or None) clears it so the handler tracks sys.stderr again.
import logging
import sys

_cls = logging._StderrHandler

if getattr(_cls.stream, "fset", None) is None:
    def _get(self):
        override = self.__dict__.get("_stream_override")
        return sys.stderr if override is None else override

    def _set(self, value):
        self.__dict__["_stream_override"] = None if value is sys.stderr else value

    _cls.stream = property(_get, _set)
