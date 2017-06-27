"""
This module contains tools for automatically updating tools that can be
downloaded via HTTP. It expects to be given a download page, which it then
scrapes for a file to download via a custom regular expression. If the version
number is not equal to a previously-downloaded version number, then a
customisable series of build commands is run.
"""
from phyre_engine.component.component import Component
from pathlib import Path
import urllib.request
import re
from html.parser import HTMLParser
from . errors import UpdateError

import yaml
import subprocess
try:
    from yaml import CSafeLoader as SafeLoader, CSafeDumper as SafeDumper
except ImportError:
    from yaml import SafeLoader, SafeDumper

class HTTPUpdater(Component):
    """
    :param str name: Name of the tool to be downloaded.
    :param str version_db: Path of a YAML file containing tool versions.
    :param str download_dir: Directory in which to save the archive of the tool.
    :param str install_dir: Directory into which the tool will be installed.
    :param str download_page: URL of the page containing links to download the
        tool.
    :param str link_regex: Regular expression matching the URL of the download
        link itself. The regex must contain the ``version`` named capture.
    :param list[str] build_commands: List of shell commands that will be applied
        to build the tool.

    .. warning::
        Commands are executed with a shell. You should *never* supply build
        commands from untrusted sources.
    """
    ADDS = []
    REMOVES = []
    REQUIRED = []


    class _LinkParser(HTMLParser):
        """Used to find the download link. Calls "callback" with the match
        object when a link matches "regex"."""
        def __init__(self, regex, callback):
            super().__init__()
            self.regex = regex
            self.callback = callback

        def handle_starttag(self, tag, attrs):
            if tag == "a":
                attr_dict = {a[0]: a[1] for a in attrs}
                if "href" in attr_dict:
                    match = self.regex.search(attr_dict["href"])
                    if match:
                        self.callback(match)

        def error(self, message):
            print(message)

    class _MetaCharsetParser(HTMLParser):
        """Finds any meta tag giving the charset. Calls "callback" with the
        charset."""
        def __init__(self, callback):
            super().__init__()
            self.callback = callback

        def handle_starttag(self, tag, attrs):
            if tag == "meta":
                if ("http-equiv" in attrs and "content" in attrs
                        and attrs["http-equiv"] == "content-type"):

                    charset = HTTPUpdater.parse_charset(attrs["content"])
                    if charset is not None:
                        self.callback(charset)
                elif "charset" in attrs:
                    self.callback(attrs["charset"])

        def error(self, message):
            print(message)

    def __init__(
            self, name, version_db, download_dir, install_dir,
            download_page, link_regex, build_commands):
        self.name = name
        self.version_db = Path(version_db)
        self.download_dir = Path(download_dir)
        self.install_dir = Path(install_dir)
        self.download_page = download_page
        self.link_regex = re.compile(link_regex)
        self.build_commands = build_commands


    def read_current_version(self):
        """
        Read current version from the version database.
        :return: Current version.
        :rtype: str
        """
        with self.version_db.open("r") as version_db_in:
            version_db = yaml.load(version_db_in, SafeLoader)
            if version_db is None:
                version_db = {}
        current_ver = version_db[self.name] if self.name in version_db else None
        return current_ver

    def _guess_charset(self, page_bytes, headers):
        """Guess the charset, either from headers or by parsing meta tags."""
        # Pick a content type
        charset = "UTF-8"
        if "Content-Type" in headers:
            charset = self.parse_charset(headers["Content-Type"])
        else:
            def setter(ch):
                nonlocal charset
                charset = ch
            ascii_page = page_bytes.decode("ASCII", ignore_errors=True)
            self._MetaCharsetParser(setter).feed(ascii_page)
        return charset

    def read_download_page(self):
        """Retrieve the download page as a string."""
        # Now, let's get the download page and try to find the link
        with urllib.request.urlopen(self.download_page) as download_in:
            page_bytes = download_in.read()
            headers = download_in.info()

        # Do our best to guess the charset and decode the page
        charset = self._guess_charset(page_bytes, headers)
        page = page_bytes.decode(charset)
        return page

    def find_download_link(self, page):
        """
        Finds the download link from the download page. The ``href`` of the
        link must match the regular expression ``link_regex`` supplied when
        this object was constructed.

        :return: Tuple containing the version string, archive name (extracted
            from the link href) and absolute URL pointing to the archive.
        """
        # Look for the download link
        download_matches = []
        def appender(match):
            download_matches.append(match)
        self._LinkParser(self.link_regex, appender).feed(page)
        if len(download_matches) != 1:
            raise UpdateError("One matching download link required: {}".format(
                download_matches))
        download_url = download_matches[0].group(0)
        new_version = download_matches[0].group("version")

        # Get the name of the archive file to save
        parsed_url = urllib.parse.urlparse(download_url)
        # Archive name given by final component of path
        archive_name = parsed_url.path.split("/")[-1]

        # Download URL may need to be converted to absolute URL
        if not parsed_url.netloc:
            download_url = urllib.parse.urljoin(
                self.download_page, download_url)

        return new_version, archive_name, download_url


    def get_archive(self, archive_name, download_url):
        """
        Retrieve the archive containing the software.

        :param archive_name: :py:class:`pathlib.Path` pointing to the location
            in which the archive will be saved.
        :param download_url: String containing the full URL of the archive.
        """
        with (self.download_dir / archive_name).open("wb") as archive_out:
            with urllib.request.urlopen(download_url) as url_in:
                archive_out.write(url_in.read())

    def build(self, new_version, archive_name):
        """
        Execute all build commands. The following keys may be used in each
        command with standard python string formatting:

        archive:
            Name (without directories) of the saved archive.

        download_dir:
            Path of the download directory (passed as ``download_dir`` to the
            constructor).

        install_dir:
            Directory in which to install the software (the ``install_dir``
            argument to the constructor).

        .. note::
            Python string formatting is applied to each command, so if you wish
            to use a brace (``{`` or ``}``), you must double them, so ``{``
            would become ``{{``.

        .. warning::
            Commands are executed with a shell. You should *never* supply build
            commands from untrusted sources.

        """
        # Start executing build commands.
        for cmd in self.build_commands:
            formatted_cmd = cmd.format(version=new_version,
                archive=archive_name,
                download_dir=str(self.download_dir),
                install_dir=str(self.install_dir))
            print(formatted_cmd)
            subprocess.run([formatted_cmd], shell=True, check=True)


    def update_db(self, new_version):
        """Update the version database with the new version."""
        with self.version_db.open("r+") as db_fh:
            version_db = yaml.load(db_fh, SafeLoader)
            if version_db is None:
                version_db = {}
            version_db[self.name] = new_version
            db_fh.truncate()
            yaml.dump(version_db, db_fh, SafeDumper, default_flow_style=False)

    def run(self, data, config=None, pipeline=None):
        current_ver = self.read_current_version()
        page = self.read_download_page()
        new_version, archive_name, download_url = self.find_download_link(page)

        # Don't do anything if the version is the same
        if current_ver is not None and new_version == current_ver:
            return data

        self.get_archive(archive_name, download_url)
        self.build(new_version, archive_name)
        self.update_db(new_version)

        return data


    @staticmethod
    def parse_charset(content_type_header):
        match = re.search(r"charset=(\S+)", content_type_header)
        if match is not None:
            return match.group(1)
        return None
