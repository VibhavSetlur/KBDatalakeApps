from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Callable, Optional, Tuple, Dict
import time
import traceback
from concurrent.futures import ThreadPoolExecutor, Future


class TaskStatus(str, Enum):
    PENDING = "PENDING"
    RUNNING = "RUNNING"
    SUCCESS = "SUCCESS"
    FAILED = "FAILED"


@dataclass
class TaskHandle:
    fn: Callable[..., Any]
    args: Tuple[Any, ...]
    kwargs: Dict[str, Any]

    status: TaskStatus = field(default=TaskStatus.PENDING, init=False)
    result: Any = field(default=None, init=False)
    error: Optional[BaseException] = field(default=None, init=False)
    traceback: Optional[str] = field(default=None, init=False)
    started_at: Optional[float] = field(default=None, init=False)
    finished_at: Optional[float] = field(default=None, init=False)

    _future: Optional[Future] = field(default=None, init=False, repr=False)

    def done(self) -> bool:
        return self.status in (TaskStatus.SUCCESS, TaskStatus.FAILED)

    def wait(self, timeout: Optional[float] = None) -> "TaskHandle":
        """
        If running async, wait for completion. If already sync, no-op.
        """
        if self._future is not None:
            self._future.result(timeout=timeout)  # re-raises only inside future; we store error too
        return self


class TaskExecutor:

    def __init__(self, max_workers: int = 4):
        self._pool = ThreadPoolExecutor(max_workers=max_workers)

    def run_task(self, fn: Callable[..., Any], *args: Any, **kwargs: Any) -> TaskHandle:
        h = TaskHandle(fn=fn, args=args, kwargs=kwargs)
        h._future = self._pool.submit(self._run_into_handle, h)
        return h

    def shutdown(self, wait: bool = True) -> None:
        self._pool.shutdown(wait=wait)

    @staticmethod
    def _run_into_handle(handle: TaskHandle) -> None:
        handle.status = TaskStatus.RUNNING
        handle.started_at = time.time()
        try:
            handle.result = handle.fn(*handle.args, **handle.kwargs)
            handle.status = TaskStatus.SUCCESS
        except BaseException as e:
            handle.error = e
            handle.traceback = traceback.format_exc()
            handle.status = TaskStatus.FAILED
        finally:
            handle.finished_at = time.time()
