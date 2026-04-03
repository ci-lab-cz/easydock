import atexit
import json
import logging
import queue
import shlex
import subprocess
import threading
import time
from itertools import count
from typing import Any, Dict, Optional, Sequence, Union


class JsonLineProcessClient:
    """
    Lightweight client for long-lived subprocesses using JSON Lines over STDIN/STDOUT.

    Protocol:
      request: {"id": <int>, "command": <str>, "payload": <dict>}
      response: {"id": <int>, ...}
    """

    def __init__(
        self,
        command: Union[str, Sequence[str]],
        startup_timeout: float = 120.0,
        request_timeout: Optional[float] = None,
        cwd: Optional[str] = None,
        env: Optional[Dict[str, str]] = None,
        stderr_tail_lines: int = 100,
    ):
        if isinstance(command, str):
            command = shlex.split(command)
        if not isinstance(command, (list, tuple)) or not command:
            raise ValueError("command must be a non-empty string or sequence of strings")

        self._command = list(command)
        self._request_timeout = request_timeout
        self._stderr_tail_lines = max(1, int(stderr_tail_lines))
        self._stderr_lines = []  # type: List[str]
        self._stderr_lock = threading.Lock()
        self._response_queue = queue.Queue()  # type: queue.Queue
        self._request_lock = threading.Lock()
        self._request_id = count(1)
        self._closed = False

        self._proc = subprocess.Popen(
            self._command,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            cwd=cwd,
            env=env,
        )
        if self._proc.stdin is None or self._proc.stdout is None or self._proc.stderr is None:
            raise RuntimeError("Failed to create subprocess pipes for JSON client")

        self._stdout_thread = threading.Thread(target=self._read_stdout, daemon=True)
        self._stderr_thread = threading.Thread(target=self._read_stderr, daemon=True)
        self._stdout_thread.start()
        self._stderr_thread.start()

        self._wait_started(startup_timeout)
        atexit.register(self.close)

    def _wait_started(self, startup_timeout: float) -> None:
        deadline = time.time() + float(startup_timeout)
        while time.time() < deadline:
            if self._proc.poll() is not None:
                raise RuntimeError(
                    "Docking server exited during startup with code {}. STDERR tail:\n{}".format(
                        self._proc.returncode, self._stderr_tail()
                    )
                )
            # If process is alive after a short grace period, assume startup success.
            time.sleep(0.05)
            return
        raise TimeoutError("Docking server startup timed out")

    def _read_stdout(self) -> None:
        assert self._proc.stdout is not None
        for line in self._proc.stdout:
            self._response_queue.put(line.rstrip("\n"))
        self._response_queue.put(None)

    def _read_stderr(self) -> None:
        assert self._proc.stderr is not None
        for line in self._proc.stderr:
            with self._stderr_lock:
                self._stderr_lines.append(line.rstrip("\n"))
                if len(self._stderr_lines) > self._stderr_tail_lines:
                    self._stderr_lines = self._stderr_lines[-self._stderr_tail_lines:]

    def _stderr_tail(self) -> str:
        with self._stderr_lock:
            if not self._stderr_lines:
                return "<empty>"
            return "\n".join(self._stderr_lines[-self._stderr_tail_lines:])

    def is_alive(self) -> bool:
        return (not self._closed) and self._proc.poll() is None

    def _raise_if_dead(self) -> None:
        if self._proc.poll() is not None:
            raise RuntimeError(
                "Docking server is not running (exit code {}). STDERR tail:\n{}".format(
                    self._proc.returncode, self._stderr_tail()
                )
            )

    def request(
        self,
        command: str,
        payload: Optional[Dict[str, Any]] = None,
        timeout: Optional[float] = None,
    ) -> Dict[str, Any]:
        if self._closed:
            raise RuntimeError("Client is closed")

        req_id = next(self._request_id)
        body = {
            "id": req_id,
            "command": command,
            "payload": payload if payload is not None else {},
        }
        request_line = json.dumps(body, ensure_ascii=True)
        wait_timeout = self._request_timeout if timeout is None else timeout

        with self._request_lock:
            self._raise_if_dead()
            assert self._proc.stdin is not None
            self._proc.stdin.write(request_line + "\n")
            self._proc.stdin.flush()
            return self._wait_response(req_id, wait_timeout)

    def _wait_response(self, req_id: int, timeout: Optional[float]) -> Dict[str, Any]:
        deadline = None if timeout is None else (time.time() + float(timeout))
        while True:
            remaining = None if deadline is None else max(0.0, deadline - time.time())
            if remaining == 0:
                raise TimeoutError("Timeout waiting for response from docking server")

            try:
                raw = self._response_queue.get(timeout=remaining)
            except queue.Empty:
                raise TimeoutError("Timeout waiting for response from docking server")

            if raw is None:
                self._raise_if_dead()
                raise RuntimeError("Docking server STDOUT stream closed unexpectedly")

            try:
                response = json.loads(raw)
            except ValueError:
                logging.warning("Ignoring non-JSON line from docking server stdout: %s", raw)
                continue

            if response.get("id") != req_id:
                logging.warning(
                    "Ignoring out-of-order response from docking server (expected id=%s, got id=%s)",
                    req_id,
                    response.get("id"),
                )
                continue
            return response

    def close(self, terminate_timeout: float = 5.0) -> None:
        if self._closed:
            return
        self._closed = True

        if self._proc.poll() is None:
            try:
                self._proc.terminate()
                self._proc.wait(timeout=terminate_timeout)
            except Exception:
                try:
                    self._proc.kill()
                except Exception:
                    pass

