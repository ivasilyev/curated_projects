
import os


def is_file_valid(file: str, report: bool = False):
    if not os.path.exists(file):
        if report:
            print("Not found: '{}'".format(file))
        return False
    if not os.path.isfile(file):
        if report:
            print("Not a file: '{}'".format(file))
        return False
    if os.path.getsize(file) == 0:
        if report:
            print("Empty file: '{}'".format(file))
        return False
    return True


def scan_whole_dir(dir_name: str):
    out = []
    for root, dirs, files in os.walk(dir_name):
        for file in files:
            out.append(os.path.join(root, file))
    return sorted(out)


def find_file_by_tail(dir_name: str, tail: str, multiple: bool = False):
    files = [i for i in scan_whole_dir(dir_name) if i.endswith(tail)]
    if len(files) == 0:
        return ""
    if multiple:
        return files
    return files[0]


def decompress_file(file: str, directory: str = "", remove: bool = True):
    from subprocess import getoutput
    _ARCHIVE_EXTENSIONS = ("tar", "gz", "bz2", "zip")
    if not any(file.endswith(".{}".format(i)) for i in _ARCHIVE_EXTENSIONS):
        print("Nothing to extract: '{}'".format(file))
        return

    if len(directory) == 0:
        directory = os.path.splitext(file)[0]
    directory = os.path.normpath(directory)
    os.makedirs(directory, exist_ok=True)

    cmd = 'cd "{o}" && '
    if file.endswith(".tar.gz"):
        cmd += 'tar -xzf "{i}" -C {o}/'
    elif file.endswith(".tar.bz2"):
        cmd += 'tar -xjf "{i}" -C {o}/'
    elif file.endswith(".tar"):
        cmd += 'tar -xf "{i}" -C {o}/'
    elif file.endswith(".gz"):
        cmd += 'gzip -dkq "{i}"'
    elif file.endswith(".zip"):
        cmd += 'unzip "{i}" -dq {o}/'
    elif file.endswith(".rar"):
        cmd += 'unrar e "{i}" {o}/'

    print(getoutput(cmd.format(i=file, o=directory)))

    print("Extracting completed: '{}'".format(file))
    if os.path.isfile(file) and remove:
        print("Removing file: '{}'".format(file))
        os.remove(file)

