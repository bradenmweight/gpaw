from __future__ import annotations

import contextlib
import os
import sys
from pathlib import Path
from typing import IO, Any

from gpaw.mpi import MPIComm, world


def obj2str(obj: Any, indentation: str = '') -> str:
    """Convert Python object to string.

    >>> print(obj2str({'a': {'b': 42}}))
    a:
      b: 42
    """
    if isinstance(obj, dict):
        i = indentation
        txt = f'\n{i}'.join(f'{k}: {obj2str(v, i + "  ")}'
                            for k, v in obj.items())
        if i:
            return '\n' + i + txt
        return txt.replace(': \n', ':\n')
    return repr(obj)


def indent(text: Any, indentation='  ') -> str:
    r"""Indent text blob.

    >>> indent('line 1\nline 2', '..')
    '..line 1\n..line 2'
    """
    if not isinstance(text, str):
        text = str(text)
    return indentation + text.replace('\n', '\n' + indentation)


class Logger:
    def __init__(self,
                 filename: str | Path | IO[str] | None = '-',
                 comm: MPIComm | None = None):
        self.comm = comm or world

        self.fd: IO[str]

        if self.comm.rank > 0 or filename is None:
            self.fd = open(os.devnull, 'w', encoding='utf-8')
            self.close_fd = True
        elif filename == '-':
            self.fd = sys.stdout
            self.close_fd = False
        elif isinstance(filename, (str, Path)):
            self.fd = open(filename, 'w', encoding='utf-8')
            self.close_fd = True
        else:
            self.fd = filename
            self.close_fd = False

        self.indentation = ''

    def __del__(self) -> None:
        if self.close_fd:
            self.fd.close()

    @contextlib.contextmanager
    def indent(self, text):
        self(text)
        self.indentation += '  '
        yield
        self.indentation = self.indentation[2:]

    @contextlib.contextmanager
    def comment(self):
        self.indentation += '# '
        yield
        self.indentation = self.indentation[2:]

    def __call__(self, *args, **kwargs) -> None:
        if not self.fd.closed:
            i = self.indentation
            if kwargs:
                for kw, arg in kwargs.items():
                    assert kw not in ['end', 'sep', 'flush', 'file'], kw
                    print(f'{i}{kw}: {obj2str(arg, i + "  ")}',
                          file=self.fd)
            else:
                text = ' '.join(str(arg) for arg in args)
                if i:
                    text = i + text.replace('\n', '\n' + i)
                print(text, file=self.fd)
