"""Textual dashboard runtime for Reflower.

The command modules still talk only to ``Reporter``. This module owns the
Textual app, routes reporter events to panes, and runs command dispatch in a
background worker thread so the UI can keep repainting while long subprocesses
run.
"""

from __future__ import annotations

import sys
import threading
import time
import traceback
from typing import Callable, Optional

from reporting import BaseReporter


class TextualDashboardUnavailable(RuntimeError):
    """Raised when the Textual dashboard cannot be started safely."""


class TextualDashboardReporter(BaseReporter):
    """Reporter implementation that posts thread-safe events to a Textual app."""

    def __init__(self, app, logfile: Optional[str] = None, log_level: str = "INFO"):
        super().__init__(logfile=logfile, log_level=log_level)
        self._app = app
        self._last_status_update = 0.0
        self.status_throttle_seconds = 0.25

    def _emit_to_pane(self, pane: str, level: str, message: str, *, err: bool = False, indent: int = 0) -> None:
        line = self._format(level, message, indent)
        self._write_log(line)
        self._app.call_from_thread(self._app.write_reporter_event, pane, line)

    def info(self, message: str, indent: int = 0, color: str = "c") -> None:
        self._emit_to_pane("status", "INFO", message, indent=indent)

    def warning(self, message: str, indent: int = 0) -> None:
        self._emit_to_pane("status", "WARN", message, indent=indent)

    def error(self, message: str, indent: int = 0, exc_info: bool = False, color: str = "r") -> None:
        self._emit_to_pane("errors", "ERROR", message, err=True, indent=indent)
        if exc_info:
            for line in traceback.format_exc().rstrip().splitlines():
                self._emit_to_pane("errors", "TRACE", line, err=True)

    def command(self, argv) -> None:
        if isinstance(argv, (list, tuple)):
            message = " ".join(str(part) for part in argv)
        else:
            message = str(argv)
        self._emit_to_pane("output", "COMMAND", message)

    def progress(self, message: str) -> None:
        self._emit_to_pane("status", "PROGRESS", message)

    def status(self, message: str) -> None:
        line = self._format("STATUS", message)
        now = time.monotonic()
        if now - self._last_status_update < self.status_throttle_seconds:
            return
        self._last_status_update = now
        # Status is ephemeral spinner/current-state UI. Do not append it to the
        # status log or logfile, otherwise 0.1s spinner ticks flood both.
        self._app.call_from_thread(self._app.update_current_status, line)

    def stdout(self, line: str) -> None:
        self._emit_to_pane("output", "STDOUT", line)

    def stderr(self, line: str) -> None:
        self._emit_to_pane("errors", "STDERR", line, err=True)

    def close(self) -> None:
        super().close()
        self._app.call_from_thread(self._app.update_current_status, "STATUS: reporter closed")


def run_textual_dashboard(args, dispatch: Callable, *, hold: str = "never", log_level: str = "INFO"):
    """Run ``dispatch(args, reporter)`` inside a Textual dashboard.

    Textual imports and app startup errors become ``TextualDashboardUnavailable``
    so callers can safely fall back to plain reporting. Command failures are
    re-raised after the dashboard exits, preserving the original traceback.
    """

    try:
        from textual.app import App, ComposeResult
        from textual.containers import Horizontal, Vertical
        from textual.widgets import Footer, Header, RichLog, Static
    except Exception as exc:  # pragma: no cover - exercised in environments without Textual
        raise TextualDashboardUnavailable(str(exc)) from exc

    class ReflowerTextualApp(App):
        CSS = """
        Screen {
            layout: vertical;
        }
        #summary {
            height: 3;
            padding: 0 1;
            border: heavy $accent;
            content-align: left middle;
        }
        #current-status {
            height: 3;
            padding: 0 1;
            border: round $primary;
            content-align: left middle;
        }
        #panes {
            height: 1fr;
        }
        .pane {
            width: 1fr;
            height: 1fr;
            border: round $surface-lighten-1;
        }
        .pane-title {
            height: 1;
            padding: 0 1;
            background: $surface-lighten-1;
            text-style: bold;
        }
        RichLog {
            height: 1fr;
            padding: 0 1;
        }
        """
        BINDINGS = [("q", "request_quit", "Quit")]

        def __init__(self):
            super().__init__()
            self.command_started = False
            self.command_finished = False
            self.command_result = None
            self.command_exc_info = None
            self.reporter = TextualDashboardReporter(self, log_level=log_level)
            self.worker_thread: Optional[threading.Thread] = None

        def compose(self) -> ComposeResult:
            yield Header(show_clock=True)
            yield Static("Reflower running — q quits after command completion", id="summary")
            yield Static("STATUS: waiting to start", id="current-status")
            with Horizontal(id="panes"):
                with Vertical(classes="pane"):
                    yield Static("Status / Progress", classes="pane-title")
                    yield RichLog(id="status-log", wrap=True, markup=False, auto_scroll=True)
                with Vertical(classes="pane"):
                    yield Static("Commands / stdout", classes="pane-title")
                    yield RichLog(id="output-log", wrap=True, markup=False, auto_scroll=True)
                with Vertical(classes="pane"):
                    yield Static("Errors / stderr", classes="pane-title")
                    yield RichLog(id="error-log", wrap=True, markup=False, auto_scroll=True)
            yield Footer()

        def on_mount(self) -> None:
            self.command_started = True
            self.update_current_status("STATUS: dashboard initialized")
            self.worker_thread = threading.Thread(target=self._run_command, name="reflower-command", daemon=True)
            self.worker_thread.start()

        def _run_command(self) -> None:
            result = None
            exc_info = None
            try:
                result = dispatch(args, self.reporter)
            except BaseException:
                exc_info = sys.exc_info()

            try:
                self.reporter.close()
            except BaseException as close_exc:
                if exc_info is None:
                    exc_info = sys.exc_info()
                else:
                    self.call_from_thread(
                        self.write_reporter_event,
                        "errors",
                        f"WARN: Reporter cleanup failed: {close_exc}",
                    )

            self.call_from_thread(self.command_finished_callback, result, exc_info)

        def write_reporter_event(self, pane: str, line: str) -> None:
            self._write(pane, line)

        def update_current_status(self, line: str) -> None:
            try:
                self.query_one("#current-status", Static).update(line)
            except Exception:
                pass

        def command_finished_callback(self, result, exc_info) -> None:
            self.command_finished = True
            self.command_result = result
            self.command_exc_info = exc_info
            failed = exc_info is not None
            if failed:
                exc = exc_info[1]
                self._write("errors", f"ERROR: command failed: {exc}")
                self.query_one("#summary", Static).update("Reflower failed — press q to quit")
            else:
                self.update_current_status("STATUS: command completed successfully")
                self.query_one("#summary", Static).update("Reflower completed successfully — press q to quit")

            if hold == "always" or (hold == "on-error" and failed):
                return
            self.exit(result)

        def action_request_quit(self) -> None:
            if not self.command_finished:
                self.notify("Command is still running; wait for completion before quitting.", severity="warning")
                return
            self.exit(self.command_result)

        def _write(self, pane: str, line: str) -> None:
            target = {
                "status": "#status-log",
                "output": "#output-log",
                "errors": "#error-log",
            }.get(pane, "#status-log")
            try:
                self.query_one(target, RichLog).write(line)
            except Exception:
                # During early shutdown there may be no mounted widgets. Keep the
                # command result authoritative rather than failing in UI rendering.
                pass

    app = ReflowerTextualApp()
    try:
        result = app.run()
    except Exception as exc:
        if not app.command_started:
            raise TextualDashboardUnavailable(str(exc)) from exc
        raise

    if app.command_exc_info is not None:
        exc_type, exc, tb = app.command_exc_info
        raise exc.with_traceback(tb)
    return result
