
import os

DEFAULT_HEADER = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/90.0.4430.19 Safari/537.36"


def is_url(s: str):
    from urllib.parse import urlparse
    if not isinstance(s, str):
        return False
    try:
        result = urlparse(s)
        return all([result.scheme, result.netloc])
    except AttributeError:
        return False


def get_page(url: str, header: str = DEFAULT_HEADER, empty_content_retries: int = 5):
    from requests import get
    for retry in range(empty_content_retries):
        response = get(url, headers={"User-Agent": header})
        if response.status_code != 200:
            print(f"Got response with status {response.status_code} for '{url}'")
        content = response.content
        if len(content) > 0:
            return content
    print("Exceeded empty content retries count for the URL: '{}'".format(url))
    return b""


def get_file(url: str, file: str = "", force: bool = True, header: str = DEFAULT_HEADER):
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


def parse_links_from_soup(soup, prefix: str = ""):
    from urllib.parse import urljoin
    return [urljoin(prefix, j) for j in
                [i.get("href") for i in soup.find_all("a")]
            if j is not None and len(j) > 0]


def parse_table(soup):
    out = dict()
    for row_soup in soup.find_all("tr"):
        column_soups = row_soup.find_all("td")
        key = column_soups[0].text.strip()
        values = [i.text.strip() for i in column_soups[1:]]
        if len(key) > 0 and sum([len(i) for i in values]) > 0:
            out[key] = values
    return out

