# Path manipulations

import os


def get_path(folder_name):
    """Get the path of the directory 'folder_name', assuming it is abowe cwd."""
    path = os.getcwd()
    while (os.path.split(path)[1]) != folder_name:
        path = os.path.split(path)[0]

    return path

