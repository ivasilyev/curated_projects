
import os


def get_page(url: str, header: str = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.19 Safari/537.36"):
    from requests import get
    response = get(url, headers={"User-Agent": header})
    content = response.content
    if response.status_code != 200:
        print(f"Got response with status {response.status_code} for '{url}'")
    return content


def get_file(url: str, file: str = "", force: bool = True, header: str = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.19 Safari/537.36"):
    from requests import get
    from meta.utils.file_system import is_file_valid
    if is_file_valid(file) and not force:
        print(f"Skip already downloaded file: '{file}'")
        return file
    response = get(url, headers={"User-Agent": header}, stream=True)
    if response.status_code != 200:
        print(f"Got response with status {response.status_code} for '{url}'")
    if len(file) == 0:
        file = os.path.basename(url)
    else:
        os.makedirs(os.path.dirname(file), exist_ok=True)
    with open(file, "wb") as f:
        for data in response.iter_content():
            f.write(data)
        f.close()
        print(f"Downloaded: '{file}'")
    return file


def download_file_to_dir(url: str, directory: str, **kwargs):
    file = os.path.join(directory, os.path.basename(url))
    return get_file(url=url, file=file, **kwargs)


def get_soup(*args, **kwargs):
    from bs4 import BeautifulSoup
    return BeautifulSoup(get_page(*args, **kwargs), features="lxml")
