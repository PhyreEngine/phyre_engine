"""
This module contains components that interface with Pymol via the built-in
Remote Process Call (RPC) interface.

Pymol is quite slow to start, so this module provides components to start a Pymol
server, then a component to run commands using that interpreter, and a component
to cause the interpreter to exit.

:py:class:`.Init`
    Starts a new Pymol instance in the background. This instance of Pymol will
    listen for remote commands. Initialising Pymol is fairly expensive, so in
    most cases it makes sense to spawn a single instance and then run multiple
    commands later.

:py:class:`.Run`
    Run a Pymol command (or commands), connecting to an existing Pymol server
    spawned by the :py:class:`.Init` component.

:py:class:`.Quit`
    Cause an existing Pymol instance to exit.

.. note::

    All commands executed by components in this module are formatted with the
    fields in the pipeline state. This allows you to say, for example, ``load
    {model}``, with ``model`` being taken from the corresponding field of the
    pipeline state. To use a literal brace (``{`` or ``}``) in the Pymol
    command, you *must* escape it with another brace. That is, ``{{`` will
    become a ``{`` and ``}}`` a ``}``.
"""
import atexit
import os
import socket
import subprocess
import time
import xmlrpc.client

from phyre_engine.component import Component

class Init(Component):
    """
    Starts a new Pymol instance listening for remote commands in the
    background.

    Initialising Pymol is fairly expensive, so in most cases it makes sense to
    spawn a single instance and then run multiple commands later.

    This component adds the ``pymol`` field to the pipeline state. This field
    is a dictionary, containing the following keys:

    ``pid``
        Process ID of the launched process.

    ``port``
        Network port on which the RPC server is hosted.

    .. note::

        This component installs an :py:mod:`atexit` handler to terminate Pymol
        if the python interpreter running Phyre Engine quits. This can be
        disabled by passing the parameter `auto_exit = False`, but beware that
        this can potentially cause a Pymol instance to hang around, which will
        not make you popular with your system administrator.

    :param str command: Pymol command (or multiple commands, each on a separate
        line) to run when the server is initialised.

    :param str pymol: Path to the pymol executable.

    :param int port: Port on which to host the RPC server. A random port will
        be chosen if this is not specified.

    :param int max_retries: When using a random port, the maximum number of
        retries before raising a :py:exc:`.PymolCreationError`.

    :param bool quiet: If `True`, pass the ``-Q`` flag to pymol.

    :ivar server: RPC proxy to the Pymol instance, instantiated when
        :py:meth:`.run` is called.
    :vartype server: :py:class:`xmlrpc.client.ServerProxy`
    """
    REQUIRED = []
    ADDS = ["pymol"]
    REMOVES = []

    class PymolCreationError(Exception):
        """Raised when a pymol instance cannot be created."""
        pass

    # Nasty hack so we can start the internal RPC server on a different port.
    # We set nToTry to 1 to stop Pymol from incrementing the port number, and
    # instead exit with a status of 1 if the port is taken.
    PYMOL_RPC_INIT = (
        "import pymol.rpc; "
        "pymol.rpc.launch_XMLRPC(port={port:d}, nToTry=1); "
        "if not pymol.rpc.serv: cmd.quit(1);")

    PYMOL_ARGS = [
        "-c", # Command-line only
        "-K", # Keep alive, required because of "-c"
    ]

    def __init__(self, command, pymol="pymol", port=None, max_retries=5,
                 quiet=True):
        self.command = command
        self.pymol = pymol
        self.port = port
        self.max_retries = max_retries
        self.server = None
        if quiet:
            self.PYMOL_ARGS.append("-Q")

    @staticmethod
    def _free_port():
        """
        Bind to an arbitrary free port, close the socket and return the port
        number. There is a race condition here: it is entirely possible for
        this port to be used by the time Pymol actually gets around to binding
        it. In that case, we will try again with another port.
        """
        sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        sock.bind(('', 0))
        _, port = sock.getsockname()
        sock.close()
        return port

    def _healthy_instance(self, process, port):
        """
        Wait until Pymol has either started sucessfully or failed. The server
        has started successfully when we can connect to it via RPC, and has
        failed if an exit status is returned by poll().
        """
        exit_code = process.poll()
        self.logger.debug("Pymol exit status: %s", exit_code)
        while exit_code is None:
            try:
                conn = xmlrpc.client.ServerProxy(pymol_rpc_addr(port))
                conn.do("1 + 1")
                self.logger.debug("Connected to port %d", port)
                return True
            except ConnectionRefusedError as err:
                self.logger.debug(
                    "Pymol connection error on port %d", port)
                time.sleep(1)
                exit_code = process.poll()
                self.logger.debug("Pymol exit status: %s", exit_code)
        return False

    def _start_pymol(self, port):
        """Start pymol instance. Returns result of Popen"""
        init_cmd = self.PYMOL_RPC_INIT.format(port=port)
        pymol_args = [self.pymol] + self.PYMOL_ARGS + ["-d", init_cmd]
        self.logger.info(
            "Starting pymol server on port %d: (Command: %s)",
            port, pymol_args)
        proc = subprocess.Popen(pymol_args, close_fds=True)
        return proc

    def start_server(self):
        """
        Start a Pymol server.

        If `self.port` is `None`, then random ports are tried until one
        succeeds or `self.max_retries` is exceeded.

        :raises PymolCreationError: When a defined port is taken or the
            maximum number of retries is exceeded.
        """

        if self.port is None:
            # Repeat until we find a free port
            for _ in range(self.max_retries):
                port = self._free_port()
                proc = self._start_pymol(port)
                if self._healthy_instance(proc, port):
                    return port, proc
                else:
                    self.logger.warning(
                        "Could not create Pymol server on port %d", port)
            raise self.PymolCreationError("Could not find a port for pymol.")
        else:
            proc = self._start_pymol(self.port)
            if not self._healthy_instance(proc, self.port):
                raise self.PymolCreationError(
                    "Could not create pymol server on port {}".format(
                        self.port))
            return self.port, proc

    def run(self, data, config=None, pipeline=None):
        """Start Pymol server."""
        port, process = self.start_server()
        self.server = xmlrpc.client.ServerProxy(pymol_rpc_addr(port))
        atexit.register(exit_handler, self.server)

        rpc_proxy = xmlrpc.client.ServerProxy(pymol_rpc_addr(port))
        rpc_proxy.do(self.command.format(**data))
        data["pymol"] = {"port": port, "pid": process.pid}
        return data

class Run(Component):
    """
    Run commands using an existing Pymol server.

    Each field in the pipeline state is made available using Python's string
    formatting. See :py:mod:`phyre_engine.component.pymol` for details.

    :param list[str] command: Pymol commands to run.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = ["pymol"]

    def __init__(self, commands):
        self.commands = commands

    def run(self, data, config=None, pipeline=None):
        """Run Pymol command."""
        pymol_details = self.get_vals(data)

        rpc_addr = pymol_rpc_addr(pymol_details["port"])
        rpc_proxy = xmlrpc.client.ServerProxy(rpc_addr)
        for cmd in self.commands:
            rpc_proxy.do(cmd.format(**data))
        return data

class Quit(Component):
    """Close an existing Pymol instance."""
    ADDS = []
    REMOVES = ["pymol"]
    REQUIRED = ["pymol"]

    def run(self, data, config=None, pipeline=None):
        """Exit Pymol server."""
        pymol_details = self.get_vals(data)

        rpc_addr = pymol_rpc_addr(pymol_details["port"])
        rpc_proxy = xmlrpc.client.ServerProxy(rpc_addr)
        try:
            rpc_proxy.do("quit")
        except ConnectionRefusedError:
            # Expected, since we are closing the connection
            pass
        os.waitpid(pymol_details["pid"], 0)
        del data["pymol"]
        return data

def pymol_rpc_addr(port):
    """Returns the address of the Pymol RPC server."""
    return "http://localhost:{port:d}".format(port=port)

def exit_handler(server):
    """Attempt to force Pymol to quit

    Pymol spawns off a python interpreter in a different process group, so it's
    not possible to kill it with core python functions. If this ends up being an
    issue, we'll depend on psutil to do it. As it is, we'll just call the "quit"
    command explicitly.
    """
    try:
        server.do("quit")
    except ConnectionRefusedError:
        # We will always get an error, because "quit" causes pymol to close the
        # connection
        pass
