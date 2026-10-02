"""Runtime reporting adapters for Reflower.

The bioinformatics command modules should report status, progress, and
subprocess output through this boundary instead of writing directly to curses
panels.  The default UI factory prefers Textual, but the required recovery path
is always plain stdout/stderr + logfile output.
"""

from __future__ import annotations

import os
import sys
import traceback
from typing import List, Optional, Protocol, Tuple

from logfile_manager import logfile_manager


class Reporter(Protocol):
    def set_logfile(self, path: Optional[str]) -> None: ...
    def info(self, message: str, indent: int = 0, color: str = "c") -> None: ...
    def warning(self, message: str, indent: int = 0) -> None: ...
    def error(self, message: str, indent: int = 0, exc_info: bool = False, color: str = "r") -> None: ...
    def command(self, argv) -> None: ...
    def progress(self, message: str) -> None: ...
    def status(self, message: str) -> None: ...
    def stdout(self, line: str) -> None: ...
    def stderr(self, line: str) -> None: ...
    def close(self) -> None: ...


def safe_report(reporter: Reporter, method: str, *args, **kwargs) -> None:
    """Best-effort reporting for exception and cleanup paths.

    Reporter failures must not replace the command exception that is already
    being handled.  Emit a plain stderr warning if the reporter itself fails.
    """

    try:
        getattr(reporter, method)(*args, **kwargs)
    except Exception as exc:
        print(f"WARN: Reporter.{method} failed: {exc}", file=sys.stderr)


class BaseReporter:
    stream = sys.stdout
    error_stream = sys.stderr

    def __init__(self, logfile: Optional[str] = None, log_level: str = "INFO"):
        self.log_level = log_level.upper()
        self._closed = False
        self._logfile: Optional[str] = None
        if logfile:
            self.set_logfile(logfile)

    def set_logfile(self, path: Optional[str]) -> None:
        self._logfile = path
        # Preserve the legacy global path for older log helpers, while reporter
        # users read instance state so repeated main(argv) calls do not leak
        # paths across reporters.
        logfile_manager.logfile = path

    @property
    def logfile(self) -> Optional[str]:
        return self._logfile

    @logfile.setter
    def logfile(self, path: Optional[str]) -> None:
        self.set_logfile(path)

    def _format(self, level: str, message: str, indent: int = 0) -> str:
        prefix = " " * max(indent, 0)
        return f"{level}: {prefix}{message}" if level else f"{prefix}{message}"

    def _write_log(self, line: str) -> None:
        if not self.logfile:
            return
        os.makedirs(os.path.dirname(self.logfile) or ".", exist_ok=True)
        with open(self.logfile, "a", encoding="utf-8") as handle:
            handle.write(f"{line}\n")

    def _emit(self, line: str, *, err: bool = False) -> None:
        print(line, file=self.error_stream if err else self.stream)
        self._write_log(line)

    def info(self, message: str, indent: int = 0, color: str = "c") -> None:
        self._emit(self._format("INFO", message, indent))

    def warning(self, message: str, indent: int = 0) -> None:
        self._emit(self._format("WARN", message, indent))

    def error(self, message: str, indent: int = 0, exc_info: bool = False, color: str = "r") -> None:
        self._emit(self._format("ERROR", message, indent), err=True)
        if exc_info:
            for line in traceback.format_exc().rstrip().splitlines():
                self._emit(self._format("TRACE", line), err=True)

    def command(self, argv) -> None:
        if isinstance(argv, (list, tuple)):
            message = " ".join(str(part) for part in argv)
        else:
            message = str(argv)
        self._emit(self._format("COMMAND", message))

    def progress(self, message: str) -> None:
        self._emit(self._format("PROGRESS", message))

    def status(self, message: str) -> None:
        self._emit(self._format("STATUS", message))

    def stdout(self, line: str) -> None:
        self._emit(self._format("STDOUT", line))

    def stderr(self, line: str) -> None:
        self._emit(self._format("STDERR", line), err=True)

    def close(self) -> None:
        self._closed = True


class PlainReporter(BaseReporter):
    """Plain stdout/stderr reporter and required fallback path."""


class MemoryReporter(BaseReporter):
    """Test reporter that captures events without terminal dependencies."""

    def __init__(self, logfile: Optional[str] = None, log_level: str = "INFO"):
        self.events: List[Tuple[str, str]] = []
        super().__init__(logfile=logfile, log_level=log_level)

    def _emit(self, line: str, *, err: bool = False) -> None:
        self.events.append(("stderr" if err else "stdout", line))
        self._write_log(line)


class TextualReporter(PlainReporter):
    """Textual-backed reporter boundary.

    This first-pass adapter lazily validates Textual availability and keeps
    command modules behind the Reporter API.  It intentionally does not expose
    Textual widgets to command code; a later pass can route these methods
    through Textual workers/messages without changing command modules again.
    """

    def __init__(self, logfile: Optional[str] = None, log_level: str = "INFO"):
        # Lazy import so --help/plain mode do not require Textual startup.
        from textual.app import App  # noqa: F401

        super().__init__(logfile=logfile, log_level=log_level)
        self.info("Textual UI reporter initialized", color="g")


class LegacyCursesReporter(PlainReporter):
    """Legacy explicit curses selection placeholder.

    Curses is no longer the required recovery path.  Until a true legacy
    adapter is needed, explicit curses selection degrades to plain reporting
    with a warning instead of risking fragile terminal state.
    """

    def __init__(self, logfile: Optional[str] = None, log_level: str = "INFO"):
        super().__init__(logfile=logfile, log_level=log_level)
        self.warning("Legacy curses UI is disabled in this stability pass; using plain reporter.")


def build_reporter(ui: str = "textual", logfile: Optional[str] = None, log_level: str = "INFO") -> Reporter:
    ui = (ui or "textual").lower()
    if ui == "plain":
        return PlainReporter(logfile=logfile, log_level=log_level)
    if ui == "curses":
        return LegacyCursesReporter(logfile=logfile, log_level=log_level)
    if ui == "textual":
        try:
            return TextualReporter(logfile=logfile, log_level=log_level)
        except Exception as exc:
            reporter = PlainReporter(logfile=logfile, log_level=log_level)
            reporter.warning(f"Textual UI unavailable; falling back to plain reporter: {exc}")
            return reporter
    raise ValueError(f"Unsupported UI mode: {ui}")
