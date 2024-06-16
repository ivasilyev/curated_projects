
import os


def filename_only(s: str):
    return os.path.splitext(os.path.basename(s))[0]


def get_file_extension(file: str, deep: int = 1):
    split = os.path.basename(file.strip()).split(".")[::-1]
    out = []
    for sub in split:
        if 5 >= len(sub) > 1:
            out.append(str(sub))
        else:
            break
    return ".{}".format(".".join(out[:deep][::-1]))


def get_script_dir():
    import sys
    return os.path.dirname(os.path.realpath(sys.argv[0]))


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


def scan_top_level_directories(dir_name: str):
    return sorted(
        [j for j in [os.path.join(dir_name, i) for i in os.listdir(dir_name)] if os.path.isdir(j)]
    )


def find_file_by_tail(dir_name: str, tail: str, multiple: bool = False):
    files = [i for i in scan_whole_dir(dir_name) if i.endswith(tail)]
    if len(files) == 0:
        return ""
    if multiple:
        return files
    return files[0]


def find_by_regex(regex: str, dir_name: str):
    from meta.utils.primitive import safe_findall, remove_empty_values
    nodes = scan_whole_dir(dir_name)
    matches = remove_empty_values([safe_findall(regex, i) for i in nodes])
    return matches


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


def backup_file(file: str):
    from shutil import copy2
    backup_file = f"{file}.bak"
    while is_file_valid(backup_file):
        backup_file = f"{backup_file}.bak"
    copy2(file, backup_file)
    return backup_file


def symlink(source: str, destination: str):
    if not os.path.exists(destination):
        os.symlink(source, destination)